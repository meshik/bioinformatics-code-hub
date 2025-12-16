# download_geo_gse49667.R
# Programmatic dataset downloader for: GEO GSE49667 (Illumina HumanMethylation450)
# Creates:
#   - <out_dir>/idats/        (decompressed .idat files)
#   - <out_dir>/targets.csv   (minfi-compatible targets with Basename + phenotype)
#
# Two public functions:
#   - download(out_dir, gse="GSE49667")
#   - download_if_missing(out_dir, gse="GSE49667")

.extract_characteristic <- function(pheno, key) {
  # Extract values like "donor: M28" from any characteristics_ch1* column.
  ch_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
  if (length(ch_cols) == 0) return(rep(NA_character_, nrow(pheno)))

  out <- rep(NA_character_, nrow(pheno))
  for (col in ch_cols) {
    vals <- as.character(pheno[[col]])
    hit <- grepl(paste0("^", key, "\\s*:"), vals, ignore.case = TRUE)
    out[is.na(out) & hit] <- sub(paste0("^", key, "\\s*:\\s*"), "", vals[hit], ignore.case = TRUE)
  }
  out
}

.sanitize_title <- function(x) {
  # GEO titles may contain non-ASCII characters (e.g., "naïve"); normalize to a safe token.
  x2 <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x2 <- gsub("[^A-Za-z0-9_]+", "_", x2)
  x2 <- gsub("_+", "_", x2)
  x2 <- gsub("^_|_$", "", x2)
  x2
}

download <- function(out_dir, gse = "GSE49667") {
  # Download and prepare GSE49667 methylation IDATs + targets.csv
  #
  # Args:
  #   out_dir: dataset folder (will be created)
  #   gse: GEO series accession
  #
  # Returns:
  #   normalizePath(out_dir)

  # Dependencies
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Missing package GEOquery. Install with BiocManager::install('GEOquery').")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Missing package Biobase. Install with BiocManager::install('Biobase').")
  }
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Missing package R.utils. Install with install.packages('R.utils').")
  }

  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  idat_dir <- file.path(out_dir, "idats")
  dir.create(idat_dir, recursive = TRUE, showWarnings = FALSE)

  # 1) Fetch supplementary files (includes GSE49667_RAW.tar)
  message("[download] Downloading supplementary files from GEO for ", gse)
  GEOquery::getGEOSuppFiles(GEO = gse, baseDir = out_dir, makeDirectory = FALSE, fetch_files = TRUE)

  tar_files <- list.files(out_dir, pattern = paste0("^", gse, "_RAW\\.tar$"), full.names = TRUE)
  if (length(tar_files) != 1) {
    # Fallback: search recursively (sometimes GEOquery writes into nested folders)
    tar_files <- list.files(out_dir, pattern = paste0("^", gse, "_RAW\\.tar$"), full.names = TRUE, recursive = TRUE)
  }
  if (length(tar_files) < 1) stop("Could not find ", gse, "_RAW.tar after download.")
  tar_path <- tar_files[[1]]
  message("[download] Found raw tar: ", tar_path)

  # 2) Extract TAR -> (likely) gzipped idats
  message("[download] Extracting TAR to: ", idat_dir)
  utils::untar(tar_path, exdir = idat_dir)

  # 3) Decompress any .gz (typical for GEO raw archives)
  gz_files <- list.files(idat_dir, pattern = "\\.gz$", full.names = TRUE, recursive = TRUE)
  if (length(gz_files) > 0) {
    message("[download] Decompressing ", length(gz_files), " .gz files")
    for (f in gz_files) {
      R.utils::gunzip(f, overwrite = TRUE, remove = TRUE)
    }
  }

  # 4) Build Basename table from IDATs
  grn <- list.files(idat_dir, pattern = "_Grn\\.idat$", full.names = TRUE, recursive = TRUE)
  red <- list.files(idat_dir, pattern = "_Red\\.idat$", full.names = TRUE, recursive = TRUE)
  if (length(grn) == 0 || length(red) == 0) {
    stop("IDAT files not found after extraction. Expected *_Grn.idat and *_Red.idat in ", idat_dir)
  }

  basenames <- sort(sub("_Grn\\.idat$", "", grn))
  gsm <- sub("^.*(GSM[0-9]+).*$", "\\1", basename(basenames))
  basename_tbl <- data.frame(
    geo_accession = gsm,
    Basename = basenames,
    stringsAsFactors = FALSE
  )

  # 5) Pull phenotype data from GEO Series Matrix
  message("[download] Fetching GEO metadata (Series Matrix) for ", gse)
  gse_obj <- GEOquery::getGEO(gse, GSEMatrix = TRUE)
  eset <- if (is.list(gse_obj)) gse_obj[[1]] else gse_obj
  pheno <- Biobase::pData(eset)

  if (!"geo_accession" %in% colnames(pheno)) {
    stop("Expected 'geo_accession' in GEO phenotype table; not found.")
  }
  pheno_ids <- as.character(pheno$geo_accession)

  title <- if ("title" %in% colnames(pheno)) as.character(pheno$title) else pheno_ids
  title_clean <- .sanitize_title(title)
  cell_type <- sub("_rep[0-9]+$", "", title_clean)

  donor <- .extract_characteristic(pheno, "donor")
  if (all(is.na(donor))) donor <- .extract_characteristic(pheno, "individual")

  idx <- match(basename_tbl$geo_accession, pheno_ids)

  # Combine basenames + full phenotype table (keeps IDAT order)
  pheno_aligned <- pheno[idx, , drop = FALSE]

  targets <- cbind(basename_tbl, pheno_aligned)

  # Standard columns used throughout the notebook
  targets$Sample_Name <- title_clean[idx]
  targets$Sample_Group <- cell_type[idx]
  targets$Donor <- donor[idx]

  # Bring key columns to the front (keep extra GEO columns at the back)
  key <- c("Sample_Name", "Sample_Group", "Donor", "geo_accession", "Basename")
  targets <- targets[, c(key, setdiff(colnames(targets), key))]
  targets_path <- file.path(out_dir, "targets.csv")
  utils::write.csv(targets, targets_path, row.names = FALSE)
  message("[download] Wrote targets: ", targets_path)

  invisible(out_dir)
}

download_if_missing <- function(out_dir, gse = "GSE49667") {
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  targets_path <- file.path(out_dir, "targets.csv")
  idat_dir <- file.path(out_dir, "idats")

  has_idats <- dir.exists(idat_dir) &&
    length(list.files(idat_dir, pattern = "_Grn\\.idat$|_Red\\.idat$", recursive = TRUE)) > 0
  has_targets <- file.exists(targets_path)

  if (has_idats && has_targets) {
    message("[download_if_missing] Dataset already present: ", out_dir)
    return(invisible(out_dir))
  }
  download(out_dir, gse = gse)
}
