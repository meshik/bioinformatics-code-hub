# download_ubilength.R
#
# Idempotent downloader for the UbiLength LFQ proteomics example dataset.
# Writes plain-text "contract" files into:
#   <out_dir>/ubilength_ubiquitin_interactors/
#     - proteinGroups.tsv
#     - experimental_design.csv
#     - SOURCE.txt
#
# Why a contract here?
# - proteinGroups.tsv is exactly what MaxQuant produces (so this mirrors real life).
# - experimental_design.csv is the minimal sample annotation most tools need.
# - Keeping these as plain text makes notebooks language-agnostic and easy to adapt.

#' Download the UbiLength label-free proteomics dataset
#'
#' @param out_dir Character scalar. Parent directory for the dataset folder.
#' @param dep_version Character scalar. Version of DEP source tarball to use as a data container.
#' @param tar_url Optional. Override URL for the DEP source tarball.
#' @return Invisibly returns the dataset directory path.
download <- function(out_dir,
                     dep_version = "1.32.0",
                     tar_url = NULL) {
  stopifnot(is.character(out_dir), length(out_dir) == 1)
  stopifnot(is.character(dep_version), length(dep_version) == 1)

  ds_dir  <- file.path(out_dir, "ubilength_ubiquitin_interactors")
  raw_dir <- file.path(ds_dir, "raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  # DEP is deprecated, but it’s a convenient container for this tiny example dataset.
  # We only use it to materialize UbiLength into plain-text files; analysis uses proDA.
  url_candidates <- c(
    sprintf("https://bioconductor.org/packages/release/bioc/src/contrib/DEP_%s.tar.gz", dep_version),
    # Common fallback if the release URL changes or is archived:
    sprintf("https://bioconductor.org/packages/3.22/bioc/src/contrib/Archive/DEP/DEP_%s.tar.gz", dep_version)
  )
  if (!is.null(tar_url)) url_candidates <- c(tar_url, url_candidates)

  tar_path <- file.path(raw_dir, sprintf("DEP_%s.tar.gz", dep_version))

  if (!file.exists(tar_path)) {
    ok <- FALSE
    last_err <- NULL
    for (u in url_candidates) {
      tryCatch({
        utils::download.file(u, tar_path, mode = "wb", quiet = TRUE)
        ok <- TRUE
      }, error = function(e) {
        last_err <<- e
      })
      if (ok && file.exists(tar_path) && file.info(tar_path)$size > 0) break
    }
    if (!ok) {
      stop(
        "Failed to download the DEP source tarball from all candidate URLs.\n",
        "Tried:\n- ", paste(url_candidates, collapse = "\n- "), "\n\n",
        "Last error:\n", as.character(last_err)
      )
    }
  }

  tmp_dir <- tempfile("DEP_src_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  utils::untar(tar_path, exdir = tmp_dir)

  dep_data_dir <- file.path(tmp_dir, "DEP", "data")
  rda_files <- c("UbiLength.rda", "UbiLength_ExpDesign.rda")
  for (f in rda_files) {
    src <- file.path(dep_data_dir, f)
    if (!file.exists(src)) {
      stop("Expected file not found inside the tarball: ", src)
    }
    file.copy(src, file.path(raw_dir, f), overwrite = TRUE)
  }

  e <- new.env(parent = emptyenv())
  base::load(file.path(raw_dir, "UbiLength.rda"), envir = e)
  base::load(file.path(raw_dir, "UbiLength_ExpDesign.rda"), envir = e)

  if (!exists("UbiLength", envir = e, inherits = FALSE)) stop("Missing object: UbiLength")
  if (!exists("UbiLength_ExpDesign", envir = e, inherits = FALSE)) stop("Missing object: UbiLength_ExpDesign")

  protein_groups <- get("UbiLength", envir = e, inherits = FALSE)
  exp_design     <- get("UbiLength_ExpDesign", envir = e, inherits = FALSE)

  # Validate and reduce the experimental design to a minimal schema.
  needed <- c("label", "condition", "replicate")
  if (!all(needed %in% colnames(exp_design))) {
    stop(
      "UbiLength_ExpDesign is missing expected columns. Found: ",
      paste(colnames(exp_design), collapse = ", ")
    )
  }
  exp_design_min <- exp_design[, needed, drop = FALSE]

  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.table(
    protein_groups,
    file.path(ds_dir, "proteinGroups.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.csv(
    exp_design_min,
    file.path(ds_dir, "experimental_design.csv"),
    row.names = FALSE
  )

  provenance <- c(
    "UbiLength example dataset (ubiquitin interactors; MaxQuant-style proteinGroups).",
    "Analysis notebook uses proDA (missing-value-aware differential abundance).",
    sprintf("Data extracted from DEP_%s source tarball (used only as a data container).", dep_version),
    sprintf("First successful download URL: %s", url_candidates[1])
  )
  writeLines(provenance, con = file.path(ds_dir, "SOURCE.txt"))

  invisible(ds_dir)
}

#' Download UbiLength if missing (idempotent)
#'
#' @param out_dir Character scalar.
#' @param dep_version Character scalar.
#' @param tar_url Optional URL override.
#' @return Invisibly returns the dataset directory path.
download_if_missing <- function(out_dir,
                                dep_version = "1.32.0",
                                tar_url = NULL) {
  stopifnot(is.character(out_dir), length(out_dir) == 1)

  ds_dir <- file.path(out_dir, "ubilength_ubiquitin_interactors")
  contract_ok <-
    file.exists(file.path(ds_dir, "proteinGroups.tsv")) &&
    file.exists(file.path(ds_dir, "experimental_design.csv"))

  if (!contract_ok) {
    download(out_dir = out_dir, dep_version = dep_version, tar_url = tar_url)
  }

  invisible(ds_dir)
}
