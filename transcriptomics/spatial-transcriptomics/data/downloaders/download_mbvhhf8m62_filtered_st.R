# - Put two files in ds_dir with stable names:
#     filtered_matrix.tsv.gz
#     meta_data.tsv.gz

download <- function(
  ds_dir,
  zip_url = "https://data.mendeley.com/public-files/datasets/mbvhhf8m62/files/f76ec6ad-addd-41c3-9eec-56e31ddbac71/file_downloaded"
) {
  # --- Ensure output directory exists ---
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Download ZIP to a temp file ---
  tmp_zip <- tempfile(fileext = ".zip")
  tmp_dir <- tempfile(pattern = "mbvhhf8m62_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(c(tmp_zip, tmp_dir), recursive = TRUE, force = TRUE), add = TRUE)

  message("Downloading ZIP...")
  utils::download.file(zip_url, destfile = tmp_zip, mode = "wb", quiet = FALSE)

  # --- Inspect ZIP contents and pick the first match for each required file ---
  zip_names <- utils::unzip(tmp_zip, list = TRUE)$Name

  m_in_zip <- zip_names[grepl("filtered_matrix\\.tsv(\\.gz)?$", zip_names)]
  d_in_zip <- zip_names[grepl("meta_data\\.tsv(\\.gz)?$", zip_names)]
  stopifnot(length(m_in_zip) >= 1, length(d_in_zip) >= 1)

  # --- Extract only those two files (faster + cleaner than extracting the entire archive) ---
  message("Extracting required files from ZIP...")
  utils::unzip(tmp_zip, files = c(m_in_zip[1], d_in_zip[1]), exdir = tmp_dir)

  src_matrix <- file.path(tmp_dir, m_in_zip[1])
  src_meta   <- file.path(tmp_dir, d_in_zip[1])

  # --- Write standardized outputs (always .tsv.gz) ---
  # Why: the analysis notebook should not care whether the original archive stored .tsv or .tsv.gz
  write_gz <- function(src, dst_gz) {
    if (grepl("\\.gz$", src, ignore.case = TRUE)) {
      file.copy(src, dst_gz, overwrite = TRUE)
      return(invisible(TRUE))
    }
    in_con  <- file(src, "rb")
    out_con <- gzfile(dst_gz, "wb")
    on.exit({ try(close(in_con), silent = TRUE); try(close(out_con), silent = TRUE) }, add = TRUE)

    repeat {
      buf <- readBin(in_con, what = "raw", n = 1024 * 1024)
      if (!length(buf)) break
      writeBin(buf, out_con)
    }
    invisible(TRUE)
  }

  message("Writing standardized outputs...")
  write_gz(src_matrix, file.path(ds_dir, "filtered_matrix.tsv.gz"))
  write_gz(src_meta,   file.path(ds_dir, "meta_data.tsv.gz"))

  invisible(TRUE)
}

download_if_missing <- function(ds_dir) {
  # Idempotency: do nothing if both standardized outputs exist
  fm <- file.path(ds_dir, "filtered_matrix.tsv.gz")
  md <- file.path(ds_dir, "meta_data.tsv.gz")
  if (file.exists(fm) && file.exists(md)) return(invisible(TRUE))
  download(ds_dir)
  invisible(TRUE)
}
