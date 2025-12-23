download <- function(out_dir, gse = "GSE47915") {
  ds_dir   <- file.path(out_dir, "geo_gse47915_prostate_450k")
  idat_dir <- file.path(ds_dir, "idats")
  dir.create(idat_dir, recursive = TRUE, showWarnings = FALSE)

  # 1) Download + extract (handles common "RAW.tar" + nested tars)
  GEOquery::getGEOSuppFiles(GEO = gse, baseDir = ds_dir, makeDirectory = TRUE)
  raw_tar <- list.files(ds_dir, pattern = paste0("^", gse, "_RAW\\.tar$"),
                        recursive = TRUE, full.names = TRUE)[1]
  utils::untar(raw_tar, exdir = idat_dir)

  for (i in 1:2) { # two passes covers typical nested tar layouts
    nested <- list.files(idat_dir, pattern = "\\.(tar|tgz|tar\\.gz)$",
                         recursive = TRUE, full.names = TRUE)
    for (f in nested) utils::untar(f, exdir = idat_dir)
  }

  # 2) Gunzip any *.idat.gz
  if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
  gz <- list.files(idat_dir, pattern = "\\.idat\\.gz$", recursive = TRUE,
                   full.names = TRUE, ignore.case = TRUE)
  for (f in gz) R.utils::gunzip(f, overwrite = TRUE, remove = FALSE)

  # 3) Build Basename RELATIVE to ds_dir (supports nested folders; no copying)
  ds_norm <- normalizePath(ds_dir, winslash = "/", mustWork = TRUE)

  grn_files <- list.files(idat_dir, pattern = "_Grn\\.idat$", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)
  grn_rel <- sub(paste0("^", ds_norm, "/"), "",
                 normalizePath(grn_files, winslash = "/", mustWork = TRUE))

  basenames_rel <- sub("_Grn\\.idat$", "", grn_rel, ignore.case = TRUE)
  has_red <- file.exists(file.path(ds_dir, paste0(basenames_rel, "_Red.idat")))
  basenames_rel <- basenames_rel[has_red]

  sample_name <- basename(basenames_rel)
  geo_accession <- sub("^(GSM\\d+).*", "\\1", sample_name)

  # 4) Pull GEO titles + label benign/tumor
  gsm_pd <- Biobase::pData(GEOquery::getGEO(gse, GSEMatrix = TRUE)[[1]])
  meta <- data.frame(
    geo_accession = gsm_pd$geo_accession,
    title = gsm_pd$title,
    stringsAsFactors = FALSE
  )
  title <- meta$title[match(geo_accession, meta$geo_accession)]

  Sample_Group <- ifelse(grepl("benign", title, ignore.case = TRUE), "benign",
                  ifelse(grepl("tumou?r", title, ignore.case = TRUE), "tumor", NA))

  targets <- data.frame(
    Sample_Name   = sample_name,
    geo_accession = geo_accession,
    title         = title,
    Sample_Group  = Sample_Group,
    Basename      = basenames_rel,
    stringsAsFactors = FALSE
  )

  utils::write.csv(targets, file.path(ds_dir, "targets.csv"), row.names = FALSE)
  invisible(ds_dir)
}

download_if_missing <- function(out_dir, gse = "GSE47915") {
  ds_dir <- file.path(out_dir, "geo_gse47915_prostate_450k")
  ok <- file.exists(file.path(ds_dir, "targets.csv")) &&
    length(list.files(file.path(ds_dir, "idats"),
                      pattern = "_Grn\\.idat(\\.gz)?$",
                      recursive = TRUE, ignore.case = TRUE)) > 0
  if (!ok) download(out_dir, gse = gse)
  invisible(ds_dir)
}
