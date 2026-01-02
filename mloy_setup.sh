#!/bin/bash
# ==============================================================================
# mLOY Unified Setup  
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Ensure executables
chmod +x "$SCRIPT_DIR/mloy" "$SCRIPT_DIR/mloy_caller.R" "$SCRIPT_DIR/install_deps.R" "$SCRIPT_DIR/genome_audit.sh" 2>/dev/null

INPUT="$1"
shift
EXTRA_ARGS="$@"

if [[ -z "$INPUT" ]]; then
    echo "Usage: ./mloy_setup.sh <file_or_directory>" >&2
    exit 1
fi

# === OPTIMIZATION: BATCH PROCESSING FOR IDATS ===
# If input is a directory containing IDATs, we SKIP the per-file loop.
# This prevents restarting R 1000 times and lets R handle the chunking natively.
if [[ -d "$INPUT" ]]; then
    # Peek for IDAT files
    IDAT_COUNT=$(find "$INPUT" -maxdepth 2 -name "*.[iI][dD][aA][tT]" | head -n 1 | wc -l)
    
    if [[ "$IDAT_COUNT" -gt 0 ]]; then
        echo "➜ [Batch Mode] Detected IDAT directory. Delegating entirely to R engine..." >&2
        "$SCRIPT_DIR/mloy" "$INPUT" --type meth $EXTRA_ARGS
        exit 0
    fi
fi

# === FALLBACK: PER-FILE PROCESSING (Safe for BAMs/VCFs) ===
# We keep this for BAMs/VCFs because 'genome_audit.sh' performs valuable 
# repairs (fixing headers, creating indices) that R cannot do easily.

process_file() {
    local FILE=$1
    echo "------------------------------------------------" >&2
    echo "➜ Analyzing: $(basename "$FILE")" >&2

    # 1. CALL THE UNIVERSAL AUDITOR
    AUDIT_OUT=$("$SCRIPT_DIR/genome_audit.sh" "$FILE")
    local RET=$?
    
    if [[ $RET -ne 0 ]]; then
        echo "❌ Skipping: Audit failed." >&2
        return
    fi

    # 2. PARSE AUDIT RESULTS
    local FINAL_PATH=""
    local DETECTED_TYPE=""
    local BUILD=""

    for line in $AUDIT_OUT; do
        if [[ "$line" == FINAL_PATH=* ]]; then FINAL_PATH="${line#*=}"; fi
        if [[ "$line" == DETECTED_TYPE=* ]]; then DETECTED_TYPE="${line#*=}"; fi
        if [[ "$line" == BUILD=* ]]; then BUILD="${line#*=}"; fi
    done

    # 3. MAP TYPE TO ENGINE FLAG
    local ENGINE_TYPE="auto"
    if [[ "$DETECTED_TYPE" == "BAM" || "$DETECTED_TYPE" == "SAM" || "$DETECTED_TYPE" == "CRAM" ]]; then
        ENGINE_TYPE="bam"
    elif [[ "$DETECTED_TYPE" == "VCF" || "$DETECTED_TYPE" == "VCF_TEXT" ]]; then
        ENGINE_TYPE="vcf"
    elif [[ "$DETECTED_TYPE" == "IDAT" ]]; then
        ENGINE_TYPE="meth" 
    fi

    # 4. EXECUTE R ENGINE
    "$SCRIPT_DIR/mloy" "$FINAL_PATH" --type "$ENGINE_TYPE" --build "$BUILD" $EXTRA_ARGS
}

if [[ -d "$INPUT" ]]; then
    # Loop for non-IDAT directories (e.g. folder of BAMs)
    for f in "$INPUT"/*; do [[ -e "$f" ]] && process_file "$f"; done
else
    process_file "$INPUT"
fi