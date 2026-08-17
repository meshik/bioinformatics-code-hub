#' Download the UbiLength LFQ proteomics dataset
#'
#' Downloads the UbiLength example data from the pinned DEP source package
#' and writes the proteinGroups table and experimental design as plain-text files.
#'
#' @param out_dir Character scalar. Parent data directory.
#' @return Invisibly returns the dataset directory.
download <- function(out_dir) {

  # Pin the exact DEP package version that contains the UbiLength dataset.
  dep_version <- "1.32.0"
  # Pin the Bioconductor release that distributed DEP 1.32.0.
  bioc_version <- "3.22"
  # Final location for the materialized dataset files.
  ds_dir <- file.path(
    out_dir,
    "ubilength_ubiquitin_interactors"
  )

  # Try the normal version-specific Bioconductor URL first.
  # The archive URL is a fallback in case the retired package is moved.
  urls <- c(
    sprintf(
      "https://bioconductor.org/packages/%s/bioc/src/contrib/DEP_%s.tar.gz",
      bioc_version,
      dep_version
    ),
    sprintf(
      "https://bioconductor.org/packages/%s/bioc/src/contrib/Archive/DEP/DEP_%s.tar.gz",
      bioc_version,
      dep_version
    )
  )

  # Create a temporary working directory.
  # The DEP tarball and .rda files are only needed while constructing the
  # analysis-ready TSV/CSV files
  tmp_dir <- tempfile("ubilength_")
  dir.create(tmp_dir)

  # Remove the temporary directory automatically when download() finishes,
  # including when the function exits because of an error.
  on.exit(
    unlink(tmp_dir, recursive = TRUE),
    add = TRUE
  )

  # Local temporary path for the downloaded DEP source tarball.
  tar_path <- file.path(
    tmp_dir,
    sprintf("DEP_%s.tar.gz", dep_version)
  )

  # Try each candidate URL until one download succeeds.
  for (url in urls) {

    # download.file() normally throws an error when a URL fails.
    # try(..., silent = TRUE) lets us attempt the fallback URL instead.
    status <- try(
      utils::download.file(
        url,
        tar_path,
        mode = "wb",  # binary mode avoids corrupting compressed files
        quiet = TRUE
      ),
      silent = TRUE
    )

    # A successful download returns status code 0.
    # Stop trying URLs as soon as one succeeds.
    if (!inherits(status, "try-error") && status == 0) {
      break
    }
  }

  # Fail clearly if neither Bioconductor location worked.
  if (inherits(status, "try-error") || status != 0) {
    stop(
      "Failed to download DEP ",
      dep_version
    )
  }

  # Extract the source package into the temporary directory.
  # This creates a DEP/ directory containing the package source and bundled data.
  utils::untar(
    tar_path,
    exdir = tmp_dir
  )

  # DEP stores its bundled .rda dataset objects under DEP/data/.
  dep_data_dir <- file.path(
    tmp_dir,
    "DEP",
    "data"
  )

  # Load the .rda objects into an isolated environment rather than adding
  # UbiLength and UbiLength_ExpDesign directly to the caller's workspace.
  data_env <- new.env()

  # Load the MaxQuant-style proteinGroups example table.
  base::load(
    file.path(
      dep_data_dir,
      "UbiLength.rda"
    ),
    envir = data_env
  )

  # Load the corresponding sample-level experimental design.
  base::load(
    file.path(
      dep_data_dir,
      "UbiLength_ExpDesign.rda"
    ),
    envir = data_env
  )

  # Extract the proteinGroups-style table from the temporary environment.
  protein_groups <- data_env$UbiLength

  # Keep only the sample metadata needed by the downstream notebook:
  # - label: matches the LFQ intensity column names
  # - condition: biological group
  # - replicate: replicate number within condition
  experimental_design <- data_env$UbiLength_ExpDesign[
    ,
    c(
      "label",
      "condition",
      "replicate"
    )
  ]

  # Create the final dataset directory.
  # recursive = TRUE also creates out_dir if necessary.
  dir.create(
    ds_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Write the MaxQuant-style protein table as tab-separated text.
  # TSV preserves the wide table structure and is easy to read from R,
  # Python, command-line tools, and spreadsheet software.
  utils::write.table(
    protein_groups,
    file.path(
      ds_dir,
      "proteinGroups.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Write the compact sample metadata table as CSV.
  utils::write.csv(
    experimental_design,
    file.path(
      ds_dir,
      "experimental_design.csv"
    ),
    row.names = FALSE
  )

  # Record basic provenance so the generated files remain traceable to their
  # original study and exact package source.
  writeLines(
    c(
      "Dataset: UbiLength",
      "Description: LFQ data for interactors of linear ubiquitin baits of different lengths.",
      "Original study: Zhang et al., Molecular Cell (2017).",
      "DOI: 10.1016/j.molcel.2017.01.004",
      sprintf(
        "Source: DEP %s, Bioconductor %s",
        dep_version,
        bioc_version
      ),
      sprintf(
        "Downloaded from: %s",
        url
      )
    ),
    file.path(
      ds_dir,
      "SOURCE.txt"
    )
  )

  # Return the dataset directory without printing it automatically.
  invisible(ds_dir)
}


#' Download UbiLength if missing
#'
#' @param out_dir Character scalar. Parent data directory.
#' @return Invisibly returns the dataset directory.
download_if_missing <- function(out_dir) {

  # Expected location of the materialized dataset.
  ds_dir <- file.path(
    out_dir,
    "ubilength_ubiquitin_interactors"
  )

  # These files define a complete local copy of the dataset for this workflow.
  required_files <- c(
    "proteinGroups.tsv",
    "experimental_design.csv",
    "SOURCE.txt"
  )

  # Download again if any required file is absent.
  # Checking the files themselves is safer than checking only whether the
  # dataset directory exists, because a previous download may have stopped early.
  if (!all(
    file.exists(
      file.path(
        ds_dir,
        required_files
      )
    )
  )) {
    download(out_dir)
  }

  # Return the dataset directory without printing it automatically.
  invisible(ds_dir)
}