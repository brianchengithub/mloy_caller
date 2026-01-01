#!/usr/bin/env Rscript

# ==============================================================================
# mLOY Caller (v2.0 - Persistent Caching)
# ==============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(parallel)
  library(data.table)
  library(dplyr)
  library(sesame)
  library(sesameData)
})

# --- Genome Constants ---
GENOME_RULES <- list(
  GRCh37 = list(Y_NONPAR_START = 2649520, Y_NONPAR_END = 59033286, Y_LABELS = c("Y", "chrY", "24")),
  GRCh38 = list(Y_NONPAR_START = 2781479, Y_NONPAR_END = 56887139, Y_LABELS = c("Y", "chrY", "24"))
)

# --- Output Helper ---
format_mocha_row <- function(id, gender, chrom, beg, end, n_sites, rel_cov, type, cf) {
  paste(id, gender, chrom, beg, end, end - beg, n_sites, "NA", 
        sprintf("%.4f", rel_cov), type, sprintf("%.4f", cf), sep = "\t")
}

# --- Progress Helper ---
update_progress <- function(current, total) {
  pct <- round(current / total * 100, 1)
  width <- 30
  n_fill <- round((current / total) * width)
  bar <- paste0("[", strrep("=", n_fill), strrep(" ", width - n_fill), "]")
  cat(sprintf("\r%s %d/%d files (%s%%)", bar, current, total, pct), file = stderr())
  if (current == total) cat("\n", file = stderr())
}

# --- Persistent Manifest Handler ---
get_manifest_dt <- function(platform, build_key) {
    # 1. Define Cache Directory (~/.mloy_cache)
    cache_dir <- file.path(Sys.getenv("HOME"), ".mloy_cache")
    if (!dir.exists(cache_dir)) dir.create(cache_dir)
    
    filename <- paste0(platform, ".", build_key, ".manifest.tsv.gz")
    local_file <- file.path(cache_dir, filename)
    
    # 2. Check Local Cache First
    if (file.exists(local_file) && file.size(local_file) > 1000) {
        if (!interactive()) message(paste0("➜ Using cached manifest: ", local_file))
    } else {
        # 3. Download if missing
        base_url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno"
        url <- paste(base_url, platform, filename, sep="/")
        
        message(paste0("➜ Downloading manifest to cache: ", url))
        tryCatch({
            options(timeout=600)
            download.file(url, local_file, quiet=TRUE)
        }, error = function(e) {
            unlink(local_file) # Clean up partial file
            stop("Failed to download manifest. Check internet connection.")
        })
    }
    
    # 4. Read & Standardize
    header <- names(fread(local_file, nrows=0))
    col_map <- list()
    
    if ("Probe_ID" %in% header) col_map[["Probe_ID"]] <- "Probe_ID"
    else if ("probeID" %in% header) col_map[["probeID"]] <- "Probe_ID"
    
    if ("CpG_chrm" %in% header) col_map[["CpG_chrm"]] <- "Chrom"
    else if ("chrm" %in% header) col_map[["chrm"]] <- "Chrom"
    
    if ("CpG_beg" %in% header) col_map[["CpG_beg"]] <- "Start"
    else if ("start" %in% header) col_map[["start"]] <- "Start"
    
    if (length(col_map) < 3) stop(paste("Manifest missing required columns. Found:", paste(header, collapse=", ")))
    
    # Read optimized
    dt <- fread(local_file, select=names(col_map))
    setnames(dt, names(col_map), unlist(col_map))
    
    # Filter
    dt <- dt[!is.na(Start) & !is.na(Chrom)]
    setkey(dt, Probe_ID)
    
    if (!interactive()) message(paste0("➜ Manifest loaded: ", nrow(dt), " probes."))
    return(dt)
}

# ==============================================================================
# Worker Functions
# ==============================================================================

process_meth <- function(file_path, build, manifest_dt, ...) {
  prefix <- sub("_(Grn|Red).idat$", "", file_path, ignore.case = TRUE)
  
  tryCatch({
    # 1. Read Intensity
    intensities <- totalIntensities(noob(readIDATpair(prefix)))
    idat_names <- names(intensities)
    
    # 2. Intersect
    common_ids <- intersect(idat_names, manifest_dt$Probe_ID)
    
    # --- DEBUGGING SAFETY NET ---
    if(length(common_ids) < 1000) return(NULL)
    # ----------------------------

    # 3. Filter Data
    sub_man <- manifest_dt[common_ids, on="Probe_ID"]
    sub_ints <- intensities[common_ids]

    # 4. Calculate Signals
    auto_mask <- grepl("^(chr)?(1[0-9]?|2[0-2]?|[3-9])$", sub_man$Chrom)
    auto_sig  <- median(sub_ints[auto_mask], na.rm=TRUE)
    
    r <- GENOME_RULES[[build]]
    y_mask <- (sub_man$Chrom %in% r$Y_LABELS & 
               sub_man$Start >= r$Y_NONPAR_START & 
               sub_man$Start <= r$Y_NONPAR_END)
               
    y_sig  <- median(sub_ints[y_mask], na.rm=TRUE)
    
    return(c(y_sig, auto_sig))
    
  }, error = function(e) { 
      message(paste0("Error processing ", basename(file_path), ": ", e$message))
      return(NULL) 
  })
}

# Placeholders
process_geno <- function(f, b, m) { return(NULL) } 
process_bam <- function(f, b) { return(NULL) }
process_vcf <- function(f, b) { return(NULL) }
process_cel <- function(f, b) { return(NULL) }

calculate_metrics <- function(y, auto, id, build) {
  if (is.null(y) || is.na(auto) || auto == 0) return(NULL)
  ratio <- y / auto
  lrr   <- log2(ratio / 0.5)
  
  # --- CORRECTED FORMULA ---
  cf    <- max(0, min(1, 1 - (2^lrr)))
  # -------------------------
  
  type <- "Undetermined"; gender <- "U"
  if (lrr < -0.2) { type <- "Loss"; gender <- "M" }
  else if (lrr > -0.15 && lrr < 0.15) { type <- "Normal"; gender <- "M" }
  else { gender <- "F"; if (lrr < -2.0) type <- "Normal_F" }
  r <- GENOME_RULES[[build]]
  format_mocha_row(id, gender, "Y", r$Y_NONPAR_START, r$Y_NONPAR_END, 0, lrr, type, cf)
}

# ==============================================================================
# Main
# ==============================================================================

main <- function() {
  parser <- ArgumentParser(description="mLOY Detector")
  parser$add_argument("input", help="Input Directory or File")
  parser$add_argument("--build", default="GRCh38", help="Genome Build (default: GRCh38)")
  parser$add_argument("--cores", default=4, type="integer")
  parser$add_argument("--manifest", default=NULL)
  parser$add_argument("--chunk", default=50, type="integer")
  parser$add_argument("--quiet", action="store_true")
  
  args <- parser$parse_args()
  
  # 1. File Discovery
  if (dir.exists(args$input)) {
    files <- list.files(args$input, full.names=TRUE, recursive=TRUE)
  } else { files <- args$input }
  
  files <- files[grepl("\\.(idat|bam|cram|vcf|bcf|vcf\\.gz|cel)$", files, ignore.case=TRUE)]
  files <- files[!grepl("_Red.idat$", files, ignore.case=TRUE)]
  
  if (length(files) == 0) stop("No supported files found.")
  
  # 2. Pre-Flight
  first_idat <- files[grep("idat", files)][1]
  meth_manifest_dt <- NULL
  
  if (!is.na(first_idat)) {
      if(!args$quiet) message("Detecting platform from first sample...")
      suppressMessages({
          tmp_dat <- sesame::readIDATpair(sub("_(Grn|Red).idat$", "", first_idat, ignore.case=TRUE))
          tmp_int <- sesame::totalIntensities(sesame::noob(tmp_dat))
          platform <- attr(tmp_int, "platform")
          if (is.null(platform)) platform <- "EPIC"
      })
      if(!args$quiet) message(paste0("Platform detected: ", platform))
      
      build_key <- ifelse(args$build=="GRCh37", "hg19", "hg38")
      meth_manifest_dt <- get_manifest_dt(platform, build_key)
  }

  # 3. Execution
  chunks <- split(files, ceiling(seq_along(files) / args$chunk))
  total_files <- length(files)
  processed_count <- 0
  
  cat(paste("sample_id", "computed_gender", "chrom", "beg_pos", "end_pos", 
            "length", "n_sites", "bdev", "rel_cov", "type", "cf", sep="\t"), "\n")
  
  if(!args$quiet) message(paste0("Starting analysis of ", total_files, " files..."))
  
  for (i in seq_along(chunks)) {
    chunk_files <- chunks[[i]]
    
    results <- mclapply(chunk_files, function(f) {
      ext <- tolower(tools::file_ext(f))
      res <- NULL
      
      if (grepl("idat", ext)) {
        if (!is.null(meth_manifest_dt)) res <- process_meth(f, args$build, meth_manifest_dt)
      } 
      
      if (!is.null(res)) return(calculate_metrics(res[1], res[2], basename(f), args$build))
      return(NULL)
    }, mc.cores = args$cores)
    
    for (r in results) { if (!is.null(r)) cat(r, "\n") }
    
    processed_count <- processed_count + length(chunk_files)
    if(!args$quiet) update_progress(processed_count, total_files)
    
    rm(results); gc(verbose=FALSE)
  }
}

if (!interactive()) main()