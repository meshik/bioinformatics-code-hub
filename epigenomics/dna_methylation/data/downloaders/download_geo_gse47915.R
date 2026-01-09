download <- function(out_dir, gse = "GSE47915") {
  ds_dir <- file.path(out_dir, "geo_gse47915_prostate_450k")
  idat_dir <- file.path(ds_dir, "idats")
  dir.create(idat_dir, recursive = TRUE, showWarnings = FALSE)

  GEOquery::getGEOSuppFiles(GEO = gse, baseDir = ds_dir)  # download GEO supplementary files
  raw_tar <- file.path(ds_dir, gse, paste0(gse, "_RAW.tar"))
  utils::untar(raw_tar, exdir = idat_dir)  # extract into idats/

  gz <- list.files(idat_dir, pattern = "\\.idat\\.gz$", full.names = TRUE)  # find compressed IDATs
  for (f in gz) R.utils::gunzip(f, overwrite = TRUE, remove = FALSE)  # decompress to *.idat

  grn <- list.files(idat_dir, pattern = "_Grn\\.idat$", full.names = TRUE)  # Green files define sample basenames
  sample_name <- sub("_Grn\\.idat$", "", basename(grn))  # strip channel suffix
  geo_accession <- sub("_.*$", "", sample_name)  # assume GSM is before first underscore

  pd <- Biobase::pData(GEOquery::getGEO(gse, GSEMatrix = TRUE)[[1]])  # GEO sample metadata table
  title <- pd$title[match(geo_accession, pd$geo_accession)]  # map GSM -> title

  Sample_Group <- ifelse(grepl("benign", title, ignore.case = TRUE),
                         "benign", "tumor")

  # minfi sample sheet
  targets <- data.frame( 
    Sample_Name   = sample_name,
    geo_accession = geo_accession,
    Sample_Group  = Sample_Group,
    Basename      = file.path("idats", sample_name)  # relative to ds_dir
  )
  utils::write.csv(targets, file.path(ds_dir, "targets.csv"), row.names = FALSE)
  invisible(ds_dir)  # return dataset folder
}

download_if_missing <- function(out_dir, gse = "GSE47915") {
  ds_dir <- file.path(out_dir, "geo_gse47915_prostate_450k")  # dataset root folder
  if (!dir.exists(ds_dir)) download(out_dir, gse = gse)  # assume folder presence means data exists
  invisible(ds_dir)  # return dataset folder
}
