#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs <- c(
  "sesame",           # Methylation
  "illuminaio",       # Genotyping
  "oligo",            # CEL
  "Rsamtools",        # BAM/CRAM
  "VariantAnnotation",# VCF/BCF
  "GenomicRanges", 
  "argparse", 
  "data.table", 
  "dplyr"
)

message("Installing required packages...")
BiocManager::install(pkgs, update=FALSE, ask=FALSE)
message("Done.")
