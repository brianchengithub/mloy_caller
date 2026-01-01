#!/usr/bin/env Rscript

# ==============================================================================
# mLOY Caller  
# ==============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(parallel)
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
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

# --- Manifest Handler ---
get_manifest <- function(platform, build_key) {
    # 1. Construct URL (GitHub raw link)
    base_url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno"
    filename <- paste0(platform, ".", build_key, ".manifest.tsv.gz")
    url <- paste(base_url, platform, filename, sep="/")
    
    tmp_file <- file.path(tempdir(), filename)
    
    # 2. Download (Blocking - Single Thread)
    if (!file.exists(tmp_file) || file.size(tmp_file) < 1000) {
        message(paste0("➜ Downloading manifest: ", url))
        tryCatch({
            options(timeout=600)
            download.file(url, tmp_file, quiet=TRUE)
        }, error = function(e) stop("Failed to download manifest. Check internet connection."))
    }
    
    # 3. Read & Standardize Columns
    # Check header first
    header <- names(fread(tmp_file, nrows=0))
    
    col_map <- list()
    # Map Probe ID
    if ("Probe_ID" %in% header) col_map[["Probe_ID"]] <- "probeID"
    else if ("probeID" %in% header) col_map[["probeID"]] <- "probeID"
    
    # Map Chromosome
    if ("CpG_chrm" %in% header) col_map[["CpG_chrm"]] <- "seqnames"
    else if ("chrm" %in% header) col_map[["chrm"]] <- "seqnames"
    
    # Map Start Position
    if ("CpG_beg" %in% header) col_map[["CpG_beg"]] <- "start"
    else if ("start" %in% header) col_map[["start"]] <- "start"
    
    # Validation
    if (length(col_map) < 3) stop(paste("Manifest missing required columns. Found:", paste(header, collapse=", ")))
    
    # Read Data
    dt <- fread(tmp_file, select=names(col_map))
    setnames(dt, names(col_map), unlist(col_map))
    
    # Filter & Convert
    dt <- dt[!is.na(start) & !is.na(seqnames)]
    gr <- GRanges(seqnames = dt$seqnames, ranges = IRanges(start = dt$start, width = 1), names = dt$probeID)
    
    message(paste0("➜ Manifest loaded: ", length(gr), " probes."))
    return(gr)
}

# ==============================================================================
# Worker Functions
# ==============================================================================

process_meth <- function(file_path, build, probes_gr, ...) {
  prefix <- sub("_(Grn|Red).idat$", "", file_path, ignore.case = TRUE)
  
  tryCatch({
    # 1. Read Intensity (Using sesame)
    intensities <- totalIntensities(noob(readIDATpair(prefix)))
    
    # 2. Intersect (Fast set operation)
    common <- intersect(names(intensities), names(probes_gr))
    
    # Safety Check
    if(length(common) < 1000) {
       return(NULL) 
    }
    
    # 3. Extract Signals
    seqs <- as.character(seqnames(probes_gr[common]))
    ints <- intensities[common]
    pos  <- start(probes_gr[common])
    
    # 4. Compute Median Signals
    auto_mask <- grep("^(chr)?(1[0-9]?|2[0-2]?|[3-9])$", seqs)
    auto_sig  <- median(ints[auto_mask], na.rm=TRUE)
    
    r <- GENOME_RULES[[build]]
    y_mask <- which(seqs %in% r$Y_LABELS & pos >= r$Y_NONPAR_START & pos <= r$Y_NONPAR_END)
    y_sig  <- median(ints[y_mask], na.rm=TRUE)
    
    return(c(y_sig, auto_sig))
    
  }, error = function(e) { 
      message(paste0("Error processing ", basename(file_path), ": ", e$message))
      return(NULL) 
  })
}

# Placeholders for other types to maintain script integrity
process_geno <- function(f, b, m) { return(NULL) } 
process_bam <- function(f, b) { return(NULL) }
process_vcf <- function(f, b) { return(NULL) }
process_cel <- function(f, b) { return(NULL) }

calculate_metrics <- function(y, auto, id, build) {
  if (is.null(y) || is.na(auto) || auto == 0) return(NULL)
  ratio <- y / auto
  lrr   <- log2(ratio / 0.5)
  cf    <- max(0, min(1, 2 * (1 - (2^lrr))))
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
  # Exclude Red files
  files <- files[!grepl("_Red.idat$", files, ignore.case=TRUE)]
  
  if (length(files) == 0) stop("No supported files found.")
  
  # 2. Pre-Flight: Detect Platform & Download Manifest (ONCE)
  # We peek at the first IDAT file to determine the platform
  first_idat <- files[grep("idat", files)][1]
  meth_manifest <- NULL
  
  if (!is.na(first_idat)) {
      if(!args$quiet) message("Detecting platform from first sample...")
      suppressMessages({
          # Read just header/intensities to get attr
          tmp_dat <- sesame::readIDATpair(sub("_(Grn|Red).idat$", "", first_idat, ignore.case=TRUE))
          tmp_int <- sesame::totalIntensities(sesame::noob(tmp_dat))
          platform <- attr(tmp_int, "platform")
          if (is.null(platform)) platform <- "EPIC"
      })
      if(!args$quiet) message(paste0("Platform detected: ", platform))
      
      # Download/Load Manifest BEFORE parallel loop
      build_key <- ifelse(args$build=="GRCh37", "hg19", "hg38")
      meth_manifest <- get_manifest(platform, build_key)
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
        if (!is.null(meth_manifest)) res <- process_meth(f, args$build, meth_manifest)
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