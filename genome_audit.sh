#!/bin/bash

# ==============================================================================
# Genome Audit  
# ==============================================================================

OS_TYPE="$(uname -s)"
# Cleanup Trap: Ensures temp folders are deleted on exit, crash, or Ctrl+C
trap '[[ -n "$TEMP_EXTRACT_DIR" ]] && rm -rf "$TEMP_EXTRACT_DIR"' EXIT

echo "➜ [Audit] OS: $OS_TYPE" >&2

# --- TOOLS CHECK & SAFE INSTALL ---
can_install_tools() {
    # 1. If we are not in an interactive shell (e.g., cron, batch job), 
    #    we cannot prompt for a password. Fail early to avoid hanging.
    if [[ ! -t 0 ]]; then
        # Check if we have passwordless sudo
        if ! sudo -n true 2>/dev/null; then
            return 1 # Cannot install non-interactively
        fi
    fi
    return 0 # Can proceed to prompt or use passwordless sudo
}

install_tools() {
    echo "⚠️  [Setup] Required tools (samtools/bcftools) missing." >&2
    
    if ! can_install_tools; then
        echo "❌ [Error] Missing tools and cannot install (non-interactive shell)." >&2
        echo "   Please install 'samtools' and 'bcftools' manually." >&2
        exit 1
    fi

    echo "➜ [Setup] Attempting automatic installation..." >&2
    if [[ "$OS_TYPE" == "Darwin" ]]; then
        command -v brew >/dev/null 2>&1 || { /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"; }
        brew install bcftools samtools htslib
    elif [[ "$OS_TYPE" == "Linux" ]]; then
        # Check sudo explicitly
        if command -v sudo >/dev/null; then
            sudo apt-get update && sudo apt-get install -y bcftools samtools htslib-test unzip
        else
             echo "❌ [Error] 'sudo' not found. Cannot auto-install tools." >&2
             exit 1
        fi
    fi
}

command -v samtools >/dev/null 2>&1 || install_tools
command -v bcftools >/dev/null 2>&1 || install_tools

INPUT_FILE="$1"
if [[ -z "$INPUT_FILE" ]]; then echo "Usage: ./genome_audit.sh <file>" >&2; exit 1; fi

# Normalize path
ABS_PATH=$(cd "$(dirname "$INPUT_FILE")" && pwd)/$(basename "$INPUT_FILE")
WORKING_FILE="$ABS_PATH"

detect_build() {
    local LEN=$1
    if [[ -n "$LEN" && "$LEN" -gt 249000000 ]]; then echo "GRCh37"; else echo "GRCh38"; fi
}

# ==============================================================================
# PHASE 0: PRE-FLIGHT (Handling Indices & Containers)
# ==============================================================================
echo "➜ [Audit] Pre-flight Check..." >&2

# 1. HANDLE INDICES (.tbi, .bai, .crai)
if [[ "$WORKING_FILE" == *.tbi || "$WORKING_FILE" == *.bai || "$WORKING_FILE" == *.crai ]]; then
    BASE_FILE="${WORKING_FILE%.*}"
    if [[ ! -f "$BASE_FILE" ]]; then BASE_FILE="${BASE_FILE%.*}"; fi

    if [[ -f "$BASE_FILE" ]]; then
        echo "⚠️  [Smart-Fix] You provided an Index file ($INPUT_FILE)." >&2
        echo "➜ [Smart-Fix] Switching to Data file: $(basename "$BASE_FILE")" >&2
        WORKING_FILE="$BASE_FILE"
    else
        echo "❌ [Error] You provided an Index file, but the Data file is missing." >&2
        exit 1
    fi
fi

# 2. HANDLE ZIP CONTAINERS (.zip) - NOW SAFE
MAGIC_HEX_PRE=$(od -t x1 -N 4 "$WORKING_FILE" 2>/dev/null | head -n 1 | awk '{$1=""; print $0}' | xargs)
if [[ "$MAGIC_HEX_PRE" == "50 4b 03 04" ]]; then
    echo "⚠️  [Container] Detected ZIP archive." >&2
    echo "➜ [Container] Extracting to temp storage..." >&2
    
    # Create temp dir and register it for the TRAP cleanup
    TEMP_EXTRACT_DIR=$(mktemp -d)
    
    unzip -q -o "$WORKING_FILE" -d "$TEMP_EXTRACT_DIR"
    
    # Find largest file
    LARGEST_FILE=$(find "$TEMP_EXTRACT_DIR" -type f -not -path '*/.*' -exec ls -S {} + | head -n 1)
    
    if [[ -f "$LARGEST_FILE" ]]; then
        echo "➜ [Container] Analyzing extracted file: $(basename "$LARGEST_FILE")" >&2
        WORKING_FILE="$LARGEST_FILE"
    else
        echo "❌ [Error] ZIP file appeared empty or invalid." >&2; exit 1
    fi
fi

# ==============================================================================
# PHASE 1: MAGIC BYTE DETECTION
# ==============================================================================
echo "➜ [Audit] Inspecting File Signature..." >&2

MAGIC_HEX=$(od -t x1 -N 4 "$WORKING_FILE" 2>/dev/null | head -n 1 | awk '{$1=""; print $0}' | xargs)
MAGIC_ASCII=$(head -c 4 "$WORKING_FILE")

TYPE="UNKNOWN"

if [[ "$MAGIC_ASCII" == "CRAM" ]]; then
    TYPE="CRAM"
elif [[ "$MAGIC_HEX" == "1f 8b"* ]]; then
    if samtools view -H "$WORKING_FILE" 2>/dev/null | head -n 1 | grep -q "@HD"; then
        TYPE="BAM"
    else
        TYPE="VCF"
    fi
elif [[ "$MAGIC_ASCII" == "@HD" || "$MAGIC_ASCII" == "@SQ" ]]; then
    TYPE="SAM"
elif [[ "$MAGIC_ASCII" == "##fi" || "$MAGIC_ASCII" == "##FI" ]]; then
    TYPE="VCF_TEXT"
elif [[ "$MAGIC_ASCII" == "IDAT" ]]; then
    TYPE="IDAT"
fi

echo "   [Signature] Detected Real Type: $TYPE" >&2

# ==============================================================================
# PHASE 2: NORMALIZATION & HEALING
# ==============================================================================

if [[ "$TYPE" == "SAM" || "$TYPE" == "BAM" || "$TYPE" == "CRAM" ]]; then
    if [[ "$TYPE" == "SAM" ]]; then
        echo "⚠️  [Heal] Converting Text SAM to Binary BAM..." >&2
        NEW_FILE="${WORKING_FILE%.*}.fixed.bam"
        samtools view -Sb "$WORKING_FILE" > "$NEW_FILE"
        WORKING_FILE="$NEW_FILE"
        TYPE="BAM"
    fi

    if [[ "$TYPE" == "CRAM" && "$WORKING_FILE" != *.cram ]]; then
        mv "$WORKING_FILE" "${WORKING_FILE%.*}.cram"
        WORKING_FILE="${WORKING_FILE%.*}.cram"
    fi

    HEADER=$(samtools view -H "$WORKING_FILE" 2>/dev/null)
    if ! echo "$HEADER" | grep -q "SO:coordinate"; then
        echo "⚠️  [Heal] File is Unsorted. Sorting now (this may take time)..." >&2
        SORTED_FILE="${WORKING_FILE%.*}.sorted.bam"
        samtools sort -@ 4 -o "$SORTED_FILE" "$WORKING_FILE"
        WORKING_FILE="$SORTED_FILE"
    fi

    # Indexing
    if ! samtools index "$WORKING_FILE" 2>/dev/null; then
        echo "⚠️  [Heal] Index missing or failed. Re-sorting and Indexing..." >&2
        REPAIRED="${WORKING_FILE%.*}.repaired.bam"
        samtools sort -@ 4 -o "$REPAIRED" "$WORKING_FILE" && samtools index "$REPAIRED"
        WORKING_FILE="$REPAIRED"
    fi

    # Detect Build
    LEN=$(echo "$HEADER" | grep -E "SN:(chr)?1[[:space:]]" | head -n 1 | sed 's/.*LN:\([0-9]*\).*/\1/')
    BUILD=$(detect_build "$LEN")

elif [[ "$TYPE" == "VCF" || "$TYPE" == "VCF_TEXT" ]]; then
    if [[ "$TYPE" == "VCF_TEXT" ]]; then
        bgzip -c "$WORKING_FILE" > "${WORKING_FILE}.gz"
        WORKING_FILE="${WORKING_FILE}.gz"
    fi
    
    if ! bcftools index -t "$WORKING_FILE" 2>/dev/null; then
         echo "⚠️  [Heal] Index failed. Re-compressing to BGZF..." >&2
         RECOMP="${WORKING_FILE%.*}.bgz.vcf.gz"
         gunzip -c "$WORKING_FILE" | bgzip > "$RECOMP"
         WORKING_FILE="$RECOMP"
         bcftools index -t "$WORKING_FILE"
    fi
    
    LEN=$(bcftools view -h "$WORKING_FILE" | grep -E "##contig=<ID=(chr)?1," | head -n 1 | sed 's/.*length=\([0-9]*\).*/\1/')
    BUILD=$(detect_build "$LEN")

elif [[ "$TYPE" == "IDAT" ]]; then
    if [[ "$WORKING_FILE" == *"_Grn.idat" && ! -f "${WORKING_FILE%_Grn.idat}_Red.idat" ]]; then
        echo "❌ [Error] Orphan Green IDAT found." >&2; exit 1
    fi
    BUILD="unknown"
else
    echo "❌ [Error] Unknown File Format." >&2; exit 1
fi

echo "➜ [Audit] Final Status: READY" >&2

# OUTPUT
echo "DETECTED_TYPE=$TYPE"
echo "FINAL_PATH=$WORKING_FILE"
echo "BUILD=$BUILD"