#!/usr/bin/env Rscript

# ==============================================================================
# mLOY Dependency Installer (v2.0)
# ==============================================================================

# 1. Ensure BiocManager is present
if (!require("BiocManager", quietly = TRUE)) {
    message("➜ Installing BiocManager...")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# 2. Define Critical Packages
# Added 'sesameData' which is explicitly required by the new mloy_caller.R
pkgs <- c(
  "BiocGenerics",         # Core
  "S4Vectors",            # Core
  "IRanges",              # Core
  "GenomicRanges",        # Critical for range overlaps
  "SummarizedExperiment", # Required by sesame
  "sesame",               # Methylation analysis
  "sesameData",           # <--- NEW: Annotation data for sesame
  "illuminaio",           # IDAT parsing
  "oligo",                # CEL/Array parsing
  "Rsamtools",            # BAM/CRAM interfacing
  "VariantAnnotation",    # VCF parsing
  "argparse",             # CLI argument parsing
  "data.table",           # Fast data manipulation
  "dplyr"                 # Tidy data manipulation
)

message("➜ Checking and installing Bioconductor packages...")
message("   (This may take some time if compiling from source)")

# 3. Install & Update
# update=TRUE is crucial to fix version mismatches in GenomicRanges/S4Vectors
BiocManager::install(pkgs, update=TRUE, ask=FALSE)

# 4. Verification
missing <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(missing) > 0) {
    message(paste("❌ Error: The following packages failed to install:", paste(missing, collapse=", ")))
    quit(status=1)
} else {
    message("➜ Dependency check complete. All packages ready.")
}