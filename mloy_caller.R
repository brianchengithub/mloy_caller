#!/usr/bin/env Rscript

# ==============================================================================
# mLOY Caller  
# ==============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(parallel)
  library(data.table)
  library(dplyr)
  library(sesame)
  library(sesameData)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(Rsamtools)
})

# ==============================================================================
# CONFIGURATION & CONSTANTS
# ==============================================================================

GENOME_RULES <- list(
  GRCh37 = list(
    Y_NONPAR_START = 2649520, Y_NONPAR_END = 59033286, Y_LABELS = c("Y", "chrY", "24", "chr24"),
    X_LABELS = c("X", "chrX", "23", "chr23"), X_CHECK_START = 10000000, X_CHECK_END = 20000000
  ),
  GRCh38 = list(
    Y_NONPAR_START = 2781479, Y_NONPAR_END = 56887139, Y_LABELS = c("Y", "chrY", "24", "chr24"),
    X_LABELS = c("X", "chrX", "23", "chr23"), X_CHECK_START = 10000000, X_CHECK_END = 20000000
  ),
  hg19 = list(
    Y_NONPAR_START = 2649520, Y_NONPAR_END = 59033286, Y_LABELS = c("Y", "chrY", "24", "chr24"),
    X_LABELS = c("X", "chrX", "23", "chr23"), X_CHECK_START = 10000000, X_CHECK_END = 20000000
  ),
  hg38 = list(
    Y_NONPAR_START = 2781479, Y_NONPAR_END = 56887139, Y_LABELS = c("Y", "chrY", "24", "chr24"),
    X_LABELS = c("X", "chrX", "23", "chr23"), X_CHECK_START = 10000000, X_CHECK_END = 20000000
  )
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

format_mocha_row <- function(id, gender, chrom, beg, end, n_sites, rel_cov, type, cf, x_lrr) {
  paste(id, gender, chrom, beg, end, end - beg, n_sites, 
        sprintf("%.4f", x_lrr),   
        sprintf("%.4f", rel_cov), 
        type, 
        sprintf("%.4f", cf),      
        sep = "\t")
}

get_manifest_dt <- function(platform, build_key) {
    cache_dir <- file.path(Sys.getenv("HOME"), ".mloy_cache")
    if (!dir.exists(cache_dir)) dir.create(cache_dir)
    filename <- paste0(platform, ".", build_key, ".manifest.tsv.gz")
    local_file <- file.path(cache_dir, filename)
    if (!file.exists(local_file) || file.size(local_file) < 1000) {
        base_url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno"
        url <- paste(base_url, platform, filename, sep="/")
        tryCatch({ options(timeout=600); download.file(url, local_file, quiet=TRUE) }, 
                 error = function(e) { unlink(local_file); stop("Failed to download manifest.") })
    }
    header <- names(fread(local_file, nrows=0))
    col_map <- list()
    if ("Probe_ID" %in% header) col_map[["Probe_ID"]] <- "Probe_ID"
    if ("CpG_chrm" %in% header) col_map[["CpG_chrm"]] <- "Chrom" else if ("chrm" %in% header) col_map[["chrm"]] <- "Chrom" else if ("Chr" %in% header) col_map[["Chr"]] <- "Chrom"
    if ("CpG_beg" %in% header) col_map[["CpG_beg"]] <- "Start" else if ("start" %in% header) col_map[["start"]] <- "Start" else if ("MapInfo" %in% header) col_map[["MapInfo"]] <- "Start"
    if (length(col_map) < 3) return(NULL)
    dt <- fread(local_file, select=names(col_map))
    setnames(dt, names(col_map), unlist(col_map))
    dt <- dt[!is.na(Start) & !is.na(Chrom)]
    setkey(dt, Probe_ID)
    return(dt)
}

# ==============================================================================
# FILE PROCESSORS
# ==============================================================================

# --- VCF Processor ---
get_vcf_depth <- function(vcf_obj) {
    if ("DP" %in% names(geno(vcf_obj))) {
        dp <- geno(vcf_obj)$DP
        if (!is.null(dp) && length(dp) > 0) return(as.numeric(dp))
    }
    if ("AD" %in% names(geno(vcf_obj))) {
        ad <- geno(vcf_obj)$AD
        if (!is.null(ad)) return(sapply(ad, function(x) sum(x, na.rm=TRUE)))
    }
    return(NULL)
}

process_vcf <- function(file_path, build) {
  tryCatch({
    hdr <- scanVcfHeader(file_path)
    seq_names <- seqlevels(hdr)
    r <- GENOME_RULES[[build]]
    y_chrom <- intersect(r$Y_LABELS, seq_names)[1]
    x_chrom <- intersect(r$X_LABELS, seq_names)[1]
    auto_chrom <- intersect(c("chr20", "20", "chr22", "22"), seq_names)[1]
    
    if (is.na(y_chrom) || is.na(auto_chrom) || is.na(x_chrom)) return(NULL)
    
    y_med <- median(get_vcf_depth(readVcf(file_path, genome=build, param=ScanVcfParam(which=GRanges(y_chrom, IRanges(r$Y_NONPAR_START, r$Y_NONPAR_END)), geno=c("DP", "AD")))), na.rm=TRUE)
    x_med <- median(get_vcf_depth(readVcf(file_path, genome=build, param=ScanVcfParam(which=GRanges(x_chrom, IRanges(r$X_CHECK_START, r$X_CHECK_END)), geno=c("DP", "AD")))), na.rm=TRUE)
    a_med <- median(get_vcf_depth(readVcf(file_path, genome=build, param=ScanVcfParam(which=GRanges(auto_chrom, IRanges(1, 60000000)), geno=c("DP", "AD")))), na.rm=TRUE)
    
    return(c(y_med, a_med, x_med))
  }, error = function(e) return(NULL))
}

# --- Universal Alignment Processor (BAM/CRAM/SAM) ---
process_alignment <- function(file_path, build) {
  tryCatch({
    # Capture 'samtools idxstats' output
    cmd <- paste("samtools idxstats", shQuote(file_path))
    stats <- tryCatch({
        fread(cmd = cmd, header = FALSE, col.names = c("seqnames", "seqlength", "mapped", "unmapped"))
    }, error = function(e) return(NULL))
    
    if (is.null(stats) || nrow(stats) == 0) return(NULL)

    r <- GENOME_RULES[[build]]

    # 1. Total Mapped Reads per category
    y_reads <- sum(stats[seqnames %in% r$Y_LABELS, mapped])
    x_reads <- sum(stats[seqnames %in% r$X_LABELS, mapped])
    sex_mito_regex <- "^(chr)?([XYM]|MT)$" # Regex to exclude Sex/Mito
    auto_stats <- stats[!grepl(sex_mito_regex, seqnames) & !grepl("_", seqnames) & seqlength > 1000000]
    auto_reads <- sum(auto_stats$mapped)

    if (auto_reads == 0) return(NULL)

    # 2. Dynamic Genomic Length Calculation (reads per BP)
    # We sum the actual lengths of the autosomes found in THIS bam header
    total_auto_bp <- sum(as.numeric(auto_stats$seqlength))
    total_y_bp    <- sum(as.numeric(stats[seqnames %in% r$Y_LABELS, seqlength]))
    total_x_bp    <- sum(as.numeric(stats[seqnames %in% r$X_LABELS, seqlength]))
    
    # Avoid division by zero
    if(total_y_bp == 0) total_y_bp <- 57000000 
    if(total_x_bp == 0) total_x_bp <- 156000000

    # 3. Calculate Normalized Coverage Density (reads per bp)
    norm_y <- y_reads / total_y_bp
    norm_x <- x_reads / total_x_bp
    norm_auto <- auto_reads / total_auto_bp

    # Return density values
    return(c(norm_y, norm_auto, norm_x)) 
  }, error = function(e) return(NULL))
}

# --- Methylation (IDAT) Processor ---
process_meth <- function(file_path, build, manifest_dt) {
  prefix <- sub("_(Grn|Red).idat$", "", file_path, ignore.case = TRUE)
  tryCatch({
    intensities <- sesame::totalIntensities(sesame::noob(sesame::readIDATpair(prefix)))
    common_ids <- intersect(names(intensities), manifest_dt$Probe_ID)
    if(length(common_ids) < 1000) return(NULL)
    sub_man <- manifest_dt[common_ids, on="Probe_ID"]
    sub_ints <- intensities[common_ids]
    r <- GENOME_RULES[[build]]
    
    auto_mask <- grepl("^(chr)?(1[0-9]?|2[0-2]?|[3-9])$", sub_man$Chrom)
    auto_sig  <- median(sub_ints[auto_mask], na.rm=TRUE)
    y_mask <- (sub_man$Chrom %in% r$Y_LABELS & sub_man$Start >= r$Y_NONPAR_START & sub_man$Start <= r$Y_NONPAR_END)
    y_sig  <- median(sub_ints[y_mask], na.rm=TRUE)
    x_mask <- (sub_man$Chrom %in% r$X_LABELS)
    x_sig <- median(sub_ints[x_mask], na.rm=TRUE)
    
    return(c(y_sig, auto_sig, x_sig))
  }, error = function(e) return(NULL))
}

process_cel <- function(file_path, build) {
  tryCatch({
    requireNamespace("oligo", quietly=TRUE)
    data <- oligo::read.celfiles(file_path)
    ints <- oligo::exprs(data)
    y_med <- median(ints[grep("Y", rownames(ints)), ], na.rm=TRUE)
    a_med <- median(ints[grep("chr20", rownames(ints)), ], na.rm=TRUE)
    x_med <- median(ints[grep("X", rownames(ints)), ], na.rm=TRUE)
    return(c(y_med, a_med, x_med))
  }, error = function(e) return(NULL))
}

process_txt <- function(file_path, build) return(NULL)

# ==============================================================================
# METRIC CALCULATION
# ==============================================================================

calculate_metrics <- function(stats_vec, id, build) {
  # Expects vector: [Y_Signal, Auto_Signal, X_Signal]
  # For BAMs, these are now DENSITIES (reads/bp). For others, they are INTENSITIES.
  # The ratio logic remains valid (Density ratio vs Intensity ratio).
  
  y <- stats_vec[1]; auto <- stats_vec[2]; x <- stats_vec[3]
  
  if (is.na(auto) || auto == 0) return(NULL)
  if (is.na(y)) y <- 0
  if (is.na(x)) x <- 0
  
  # --- STEP 1: GENDER ANCHORING ---
  # X/Auto Ratio
  # BAM (Density): Male (1X vs 2A) -> 0.5. Female (2X vs 2A) -> 1.0
  # IDAT (Intensity): Similar scaling usually applies roughly.
  x_ratio <- x / auto
  x_lrr <- log2(x_ratio / 0.5) 
  
  if (x_ratio > 0.75) gender <- "F" else gender <- "M"
  
  # --- STEP 2: mLOY DETERMINATION ---
  y_ratio <- y / auto
  lrr <- log2(y_ratio / 0.5) 
  
  cf <- max(0, min(1, 1 - (2^lrr)))
  type <- "Undetermined"
  
  if (gender == "F") {
      type <- "Female"; lrr <- NA; cf <- 0
  } else {
      if (lrr < -0.2) type <- "Loss" else if (lrr <= 0.2) type <- "Normal" else type <- "Gain"
  }
  
  r <- GENOME_RULES[[build]]
  format_mocha_row(id, gender, "Y", r$Y_NONPAR_START, r$Y_NONPAR_END, 0, lrr, type, cf, x_lrr)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

main <- function() {
  parser <- ArgumentParser(description="mLOY Detector")
  parser$add_argument("input", help="Input Directory or File")
  parser$add_argument("--build", default="GRCh38")
  parser$add_argument("--cores", default=4, type="integer")
  parser$add_argument("--manifest", default=NULL)
  parser$add_argument("--chunk", default=50, type="integer")
  parser$add_argument("--quiet", action="store_true")
  parser$add_argument("--type", default="auto") 
  args <- parser$parse_args()
  
  if (dir.exists(args$input)) {
    files <- list.files(args$input, full.names=TRUE, recursive=TRUE)
  } else { files <- args$input }
  
  files <- files[grepl("\\.(idat|bam|cram|vcf|bcf|vcf\\.gz|cel|txt)$", files, ignore.case=TRUE)]
  files <- files[!grepl("_Red.idat$", files, ignore.case=TRUE)]
  if (length(files) == 0) stop("No supported files found.")
  
  meth_manifest_dt <- NULL
  if (any(grepl("idat", files, ignore.case=TRUE))) {
     first_idat <- files[grep("idat", files, ignore.case=TRUE)][1]
     tryCatch({
          tmp_dat <- sesame::readIDATpair(sub("_(Grn|Red).idat$", "", first_idat, ignore.case=TRUE))
          tmp_int <- sesame::totalIntensities(sesame::noob(tmp_dat))
          platform <- attr(tmp_int, "platform")
          if (is.null(platform)) platform <- "EPIC"
          build_key <- args$build
          if (grepl("37|19", build_key)) build_key <- "hg19"
          if (grepl("38", build_key)) build_key <- "hg38"
          meth_manifest_dt <- get_manifest_dt(platform, build_key)
      }, error=function(e) { if(!args$quiet) message(paste("Note: Manifest load failed:", e$message)) })
  }

  cat(paste("sample_id", "computed_gender", "chrom", "beg_pos", "end_pos", "length", "n_sites", "x_lrr", "y_lrr", "type", "cf", sep="\t"), "\n")
  
  chunks <- split(files, ceiling(seq_along(files) / args$chunk))
  
  for (i in seq_along(chunks)) {
    chunk_files <- chunks[[i]]
    results <- mclapply(chunk_files, function(f) {
      ext <- tolower(tools::file_ext(f))
      b_key <- args$build
      if (grepl("37|19", b_key)) b_key <- "GRCh37"
      if (grepl("38", b_key)) b_key <- "GRCh38"
      
      res <- NULL
      if (grepl("idat", ext)) {
        if (!is.null(meth_manifest_dt)) res <- process_meth(f, b_key, meth_manifest_dt)
      } else if (grepl("vcf|bcf|gz", ext)) {
         res <- process_vcf(f, b_key)
      } else if (grepl("bam|cram", ext)) {
         res <- process_alignment(f, b_key)
      } else if (grepl("cel", ext)) {
         res <- process_cel(f, b_key)
      }
      
      if (!is.null(res)) return(calculate_metrics(res, basename(f), b_key))
      return(NULL)
    }, mc.cores = args$cores)
    
    for (r in results) { if (!is.null(r)) cat(r, "\n") }
  }
}

if (!interactive()) main()