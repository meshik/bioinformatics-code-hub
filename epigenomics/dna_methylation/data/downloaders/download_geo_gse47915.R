download <- function(out_dir, gse = "GSE47915") {
  ds_dir   <- file.path(out_dir, "geo_gse47915_prostate_450k")
  idat_dir <- file.path(ds_dir, "idats")

  dir.create(idat_dir, recursive = TRUE, showWarnings = FALSE)

  # Download supplementary files (GSE47915_RAW.tar)
  GEOquery::getGEOSuppFiles(GEO = gse, baseDir = ds_dir, makeDirectory = TRUE)

  # Find the RAW tar wherever GEOquery put it
  raw_tar <- list.files(ds_dir, pattern = paste0("^", gse, "_RAW\\.tar$"),
                        recursive = TRUE, full.names = TRUE)[1]

  # Extract everything into idats/ (works whether tar contains idat.gz or nested tars)
  utils::untar(raw_tar, exdir = idat_dir)

  # If there are nested archives (some GEO series do this), extract them too
  nested_archives <- list.files(
    idat_dir,
    pattern = "\\.(tar|tgz|tar\\.gz)$",
    recursive = TRUE,
    full.names = TRUE
  )
  for (f in nested_archives) {
    utils::untar(f, exdir = idat_dir)
  }

  # Decompress *.idat.gz to *.idat
  idat_gz <- list.files(
    idat_dir,
    pattern = "\\.idat\\.gz$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(idat_gz) > 0) {
    if (!requireNamespace("R.utils", quietly = TRUE)) {
      install.packages("R.utils")
    }
    vapply(idat_gz, function(f) R.utils::gunzip(f, overwrite = TRUE, remove = FALSE), character(1))
  }

    # 5) Write targets.csv with Basename RELATIVE to ds_dir (no machine paths)
    # List the actual IDAT files (these exist, so normalizePath works)
    grn_files <- list.files(
    idat_dir,
    pattern = "_Grn\\.idat$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
    )
    red_files <- list.files(
    idat_dir,
    pattern = "_Red\\.idat$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
    )

    # Normalize only existing file paths
    ds_norm  <- normalizePath(ds_dir, winslash = "/", mustWork = TRUE)
    grn_norm <- normalizePath(grn_files, winslash = "/", mustWork = TRUE)
    red_norm <- normalizePath(red_files, winslash = "/", mustWork = TRUE)

    # Make them relative to ds_dir
    grn_rel <- sub(paste0("^", ds_norm, "/"), "", grn_norm)
    red_rel <- sub(paste0("^", ds_norm, "/"), "", red_norm)

    # Strip channel suffix to get Basename prefixes (these prefixes don't need to "exist")
    grn_base_rel <- sub("_Grn\\.idat$", "", grn_rel, ignore.case = TRUE)
    red_base_rel <- sub("_Red\\.idat$", "", red_rel, ignore.case = TRUE)

    bases_rel <- intersect(grn_base_rel, red_base_rel)

    targets <- data.frame(
    Sample_Name = basename(bases_rel),
    Basename    = bases_rel,
    stringsAsFactors = FALSE
    )

    utils::write.csv(targets, file.path(ds_dir, "targets.csv"), row.names = FALSE)



  invisible(ds_dir)
}

download_if_missing <- function(out_dir, gse = "GSE47915") {
  ds_dir   <- file.path(out_dir, "geo_gse47915_prostate_450k")
  idat_dir <- file.path(ds_dir, "idats")
  targets_path <- file.path(ds_dir, "targets.csv")

  have_idats <- length(list.files(idat_dir, pattern = "_Grn\\.idat(\\.gz)?$",
                                  recursive = TRUE, ignore.case = TRUE)) > 0
  if (file.exists(targets_path) && have_idats) {
    return(invisible(ds_dir))
  }
  download(out_dir, gse = gse)
}
