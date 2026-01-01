#!/usr/bin/env Rscript

# ==============================================================================
# mLOY Caller (v1.2 - Robust Manifest Support)
# ==============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(parallel)
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
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

# ==============================================================================
# Worker Functions
# ==============================================================================

process_meth <- function(file_path, build, ...) {
  suppressPackageStartupMessages(library(sesame))
  prefix <- sub("_(Grn|Red).idat$", "", file_path, ignore.case = TRUE)
  tryCatch({
    # 1. Read Intensity Data
    intensities <- totalIntensities(noob(readIDATpair(prefix)))
    platform <- attr(intensities, "platform") 
    if (is.null(platform)) platform <- "EPIC"
    
    # 2. Resolve Manifest (Robust Fallback Logic)
    build_key <- ifelse(build=="GRCh37", "hg19", "hg38")
    man_key <- paste0(platform, ".", build_key, ".manifest")
    
    # Check if standard manifest exists; if not, look for KYCG alternative
    if (!man_key %in% sesameDataList()) {
        # Fallback for newer sesameData versions
        alt_key <- paste0("KYCG.", platform, ".chromosome.", build_key, ".20210630")
        if (alt_key %in% sesameDataList()) {
            man_key <- alt_key
        } else {
            stop(paste("No manifest found for", platform, build_key))
        }
    }
    
    probes_gr <- sesameDataGet(man_key)
    
    # Handle GRangesList (e.g. if KYCG returns a list of chromosomes)
    if (is(probes_gr, "GRangesList")) {
        probes_gr <- unlist(probes_gr)
    }
    
    # 3. Intersect and Extract Signals
    common <- intersect(names(intensities), names(probes_gr))
    if(length(common) < 100) return(NULL)
    
    seqs <- as.character(seqnames(probes_gr[common]))
    ints <- intensities[common]
    pos  <- start(probes_gr[common])
    
    auto_mask <- grep("^(chr)?(1[0-9]?|2[0-2]?|[3-9])$", seqs)
    auto_sig  <- median(ints[auto_mask], na.rm=TRUE)
    
    r <- GENOME_RULES[[build]]
    y_mask <- which(seqs %in% r$Y_LABELS & pos >= r$Y_NONPAR_START & pos <= r$Y_NONPAR_END)
    y_sig  <- median(ints[y_mask], na.rm=TRUE)
    
    return(c(y_sig, auto_sig))
  }, error = function(e) { 
      message(paste0("\n❌ FAIL: ", basename(file_path), " -> ", e$message))
      return(NULL) 
  })
}

process_geno <- function(file_path, build, manifest_dt, ...) {
  suppressPackageStartupMessages(library(illuminaio))
  prefix <- sub("_(Grn|Red).idat$", "", file_path, ignore.case = TRUE)
  grn_f <- paste0(prefix, "_Grn.idat"); red_f <- paste0(prefix, "_Red.idat")
  if (!file.exists(grn_f) || !file.exists(red_f)) return(NULL)
  
  g_dat <- readIDAT(grn_f)$Quants; r_dat <- readIDAT(red_f)$Quants
  ids <- rownames(g_dat)
  total_sig <- g_dat[ids, "Mean"] + r_dat[ids, "Mean"]
  
  dt_sig <- data.table(Address = as.integer(ids), Signal = total_sig)
  merged <- manifest_dt[dt_sig, on = "Address", nomatch = 0]
  
  auto_sig <- median(merged[grepl("^(chr)?(1[0-9]?|2[0-2]?|[3-9])$", Chromosome), Signal], na.rm=TRUE)
  r <- GENOME_RULES[[build]]
  y_sig <- median(merged[Chromosome %in% r$Y_LABELS & Position >= r$Y_NONPAR_START & Position <= r$Y_NONPAR_END, Signal], na.rm=TRUE)
  
  return(c(y_sig, auto_sig))
}

process_bam <- function(file_path, build, ...) {
  suppressPackageStartupMessages(library(Rsamtools))
  if (!file.exists(paste0(file_path, ".bai")) && !file.exists(sub(".bam$", ".bai", file_path))) return(NULL)
  r <- GENOME_RULES[[build]]
  
  auto_gr <- GRanges(paste0("chr", 1:5), IRanges(50000000, 51000000))
  y_gr    <- GRanges("chrY", IRanges(r$Y_NONPAR_START, r$Y_NONPAR_END))
  
  get_d <- function(gr) {
    p <- ScanBamParam(which=gr, flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE))
    cnt <- countBam(file_path, param=p)
    sum(cnt$records) / sum(width(gr))
  }
  return(c(get_d(y_gr), get_d(auto_gr)))
}

process_vcf <- function(file_path, build, ...) {
  suppressPackageStartupMessages(library(VariantAnnotation))
  r <- GENOME_RULES[[build]]
  auto_gr <- GRanges(paste0("chr", 1:5), IRanges(50000000, 51000000))
  y_gr    <- GRanges("chrY", IRanges(r$Y_NONPAR_START, r$Y_NONPAR_END))
  
  get_dp <- function(gr) {
    p <- ScanVcfParam(which=gr, info="DP", geno=c("DP","AD"))
    vcf <- readVcf(file_path, genome=build, param=p)
    if ("DP" %in% names(geno(vcf))) return(median(as.numeric(geno(vcf)$DP), na.rm=TRUE))
    if ("AD" %in% names(geno(vcf))) return(median(sapply(geno(vcf)$AD, sum, na.rm=TRUE), na.rm=TRUE))
    return(NA)
  }
  return(c(get_dp(y_gr), get_dp(auto_gr)))
}

process_cel <- function(file_path, build, ...) {
  suppressPackageStartupMessages(library(oligo))
  cel <- read.celfiles(file_path)
  pm_mat <- pm(cel)
  tryCatch({
    annot <- pData(getNetAffx(cel, "probeset"))
    auto_idx <- which(annot$seqname %in% as.character(1:22))
    y_idx <- which(annot$seqname == "Y")
    auto_sig <- median(pm_mat[auto_idx, 1], na.rm=TRUE)
    r <- GENOME_RULES[[build]]
    valid_y <- annot$start[y_idx] >= r$Y_NONPAR_START & annot$start[y_idx] <= r$Y_NONPAR_END
    y_sig <- median(pm_mat[y_idx[valid_y], 1], na.rm=TRUE)
    return(c(y_sig, auto_sig))
  }, error = function(e) return(NULL))
}

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
  
  # --- FILTER: Exclude Red files to prevent double-processing pairs ---
  files <- files[!grepl("_Red.idat$", files, ignore.case=TRUE)]
  
  if (length(files) == 0) stop("No supported files found.")
  
  # 2. Manifest Loading
  manifest_dt <- NULL
  if (!is.null(args$manifest)) {
    if(!args$quiet) message("Loading manifest...")
    manifest_dt <- fread(args$manifest, select=c("Address", "Chromosome", "Position"))
    setkey(manifest_dt, Address)
  }
  
  # 3. Execution
  chunks <- split(files, ceiling(seq_along(files) / args$chunk))
  total_files <- length(files)
  processed_count <- 0
  
  cat(paste("sample_id", "computed_gender", "chrom", "beg_pos", "end_pos", 
            "length", "n_sites", "bdev", "rel_cov", "type", "cf", sep="\t"), "\n")
  
  if(!args$quiet) message(paste0("Starting analysis of ", total_files, " files (Build: ", args$build, ")..."))
  
  for (i in seq_along(chunks)) {
    chunk_files <- chunks[[i]]
    
    results <- mclapply(chunk_files, function(f) {
      ext <- tolower(tools::file_ext(f))
      res <- NULL
      if (grepl("idat", ext)) {
        if (!is.null(manifest_dt)) res <- process_geno(f, args$build, manifest_dt)
        else res <- process_meth(f, args$build)
      } 
      else if (grepl("bam|cram", ext)) res <- process_bam(f, args$build)
      else if (grepl("vcf|bcf", ext)) res <- process_vcf(f, args$build)
      else if (grepl("cel", ext)) res <- process_cel(f, args$build)
      
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