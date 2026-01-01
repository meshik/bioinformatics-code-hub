# spatial/data/downloaders/download_mbvhhf8m62_filtered_st.R

#' Download filtered spatial transcriptomics matrices for developing human heart (Mendeley mbvhhf8m62).
#'
#' This script is intentionally idempotent:
#' - download_if_missing() will do nothing if expected files already exist.
#' - download() performs download + extraction into the target dataset directory.
#'
#' Notes:
#' - Mendeley Data landing pages are sometimes JS-rendered. For reproducible automation, this script
#'   supports providing a direct ZIP URL via an environment variable:
#'     MENDELEY_MBVHHF8M62_V2_ZIP_URL
#' - That URL should point to the "Download All" ZIP for version 2 (mbvhhf8m62/2),
#'   which contains Filtered/filtered_ST_matrix_and_meta_data/{filtered_matrix,meta_data}.tsv(.gz)
#'
#' @param ds_dir Character. Output directory to populate with filtered_matrix.tsv.gz and meta_data.tsv.gz.
#' @param zip_url Character.
#' @return Invisibly TRUE on success, error otherwise.
download <- function(
  ds_dir,
  zip_url = "https://data.mendeley.com/public-files/datasets/mbvhhf8m62/files/f76ec6ad-addd-41c3-9eec-56e31ddbac71/file_downloaded"
) {
  if (!dir.exists(ds_dir)) dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  message("Downloading dataset ZIP...")
  tmp_zip <- tempfile(fileext = ".zip")
  utils::download.file(zip_url, destfile = tmp_zip, mode = "wb", quiet = FALSE)

  message("Extracting ZIP...")
  tmp_dir <- tempfile(pattern = "mbvhhf8m62_v2_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(tmp_zip, exdir = tmp_dir)

  # Locate the two files we care about (some zips store as .tsv, some as .tsv.gz)
  fm_candidates <- list.files(
    tmp_dir,
    pattern = "filtered_matrix\\.tsv(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
  md_candidates <- list.files(
    tmp_dir,
    pattern = "meta_data\\.tsv(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(fm_candidates) < 1 || length(md_candidates) < 1) {
    stop(
      paste(
        "Could not find filtered_matrix.tsv(.gz) and/or meta_data.tsv(.gz) in the extracted ZIP.",
        "\nFound filtered_matrix candidates: ", paste(fm_candidates, collapse = ", "),
        "\nFound meta_data candidates: ", paste(md_candidates, collapse = ", "),
        sep = ""
      )
    )
  }

  filtered_matrix_src <- fm_candidates[1]
  meta_data_src <- md_candidates[1]

  # Standardize names in ds_dir
  filtered_matrix_dst <- file.path(ds_dir, "filtered_matrix.tsv.gz")
  meta_data_dst <- file.path(ds_dir, "meta_data.tsv.gz")

  # Copy and gzip if needed
  copy_to_gz <- function(src, dst_gz) {
    if (grepl("\\.gz$", src)) {
      file.copy(src, dst_gz, overwrite = TRUE)
      return(invisible(TRUE))
    }
    # gzip via connection (no extra deps)
    con_in <- file(src, open = "rb")
    con_out <- gzfile(dst_gz, open = "wb")
    on.exit({
      try(close(con_in), silent = TRUE)
      try(close(con_out), silent = TRUE)
    }, add = TRUE)
    while (length(buf <- readBin(con_in, what = "raw", n = 1024 * 1024)) > 0) {
      writeBin(buf, con_out)
    }
    invisible(TRUE)
  }

  message("Writing standardized outputs...")
  copy_to_gz(filtered_matrix_src, filtered_matrix_dst)
  copy_to_gz(meta_data_src, meta_data_dst)

  message("Done.")
  invisible(TRUE)
}

#' Idempotent wrapper: download only if expected files are missing.
#' @param ds_dir Character. Dataset directory.
#' @return Invisibly TRUE if present or downloaded.
download_if_missing <- function(ds_dir) {
  fm <- file.path(ds_dir, "filtered_matrix.tsv.gz")
  md <- file.path(ds_dir, "meta_data.tsv.gz")

  if (file.exists(fm) && file.exists(md)) {
    message("Dataset already present: ", ds_dir)
    return(invisible(TRUE))
  }

  message("Dataset missing; downloading to: ", ds_dir)
  download(ds_dir = ds_dir)
  invisible(TRUE)
}
