#' Load and preprocess developing human heart spatial transcriptomics into Seurat.
#'
#' Robust to:
#' - metadata files with an extra blank/unnamed first column (common in TSV exports)
#' - spot/barcode IDs getting changed when used as column names (e.g., "-" -> ".")
#'
#' @param expression_matrix Path to filtered_matrix.tsv.gz (or .tsv)
#' @param metadata_path Path to meta_data.tsv.gz (or .tsv)
#' @return Seurat object with UMAP + clusters + a custom "spatial" reduction from new_x/new_y (if present)
fetal_human_heart <- function(expression_matrix, metadata_path) {
  # ---- helpers ----
  canon_ids <- function(x) {
    # Normalize IDs for matching without permanently altering originals
    x <- trimws(as.character(x))
    # Undo common make.names-style changes
    x <- sub("^X", "", x)               # leading X
    x <- gsub("\\.", "-", x)            # dots back to dashes (common for barcodes)
    x
  }

  drop_all_na_cols <- function(dt) {
    na_cols <- names(dt)[vapply(dt, function(col) all(is.na(col)), logical(1))]
    if (length(na_cols) > 0) dt[, (na_cols) := NULL]
    dt
  }

  # ---- read counts ----
  counts_dt <- data.table::fread(
    expression_matrix,
    fill = TRUE,
    check.names = FALSE
  )

  # First column is typically gene ID
  counts_id_col <- names(counts_dt)[1]
  counts_mat <- as.matrix(counts_dt[, -1, with = FALSE])
  rownames(counts_mat) <- as.character(counts_dt[[counts_id_col]])
  suppressWarnings(storage.mode(counts_mat) <- "numeric")

  # Candidate spot IDs might be rownames OR colnames depending on orientation
  counts_rows <- canon_ids(rownames(counts_mat))
  counts_cols <- canon_ids(colnames(counts_mat))

  # ---- read metadata ----
  meta_dt <- data.table::fread(
    metadata_path,
    fill = TRUE,
    check.names = FALSE
  )
  meta_dt <- drop_all_na_cols(meta_dt)

  # If fread warned about header mismatch, it often creates V1 for an extra leading column.
  # We'll just auto-detect the ID column instead of assuming "first column".
  best <- list(overlap = -1L, id_col = NA_character_, orient = NA_character_)
  for (mc in names(meta_dt)) {
    vals <- canon_ids(meta_dt[[mc]])
    ov_cols <- length(intersect(vals, counts_cols))
    ov_rows <- length(intersect(vals, counts_rows))

    if (ov_cols > best$overlap) best <- list(overlap = ov_cols, id_col = mc, orient = "cols")
    if (ov_rows > best$overlap) best <- list(overlap = ov_rows, id_col = mc, orient = "rows")
  }

  message("Metadata ID column chosen: ", best$id_col)
  message("Counts orientation chosen: spots in ", ifelse(best$orient == "cols", "columns", "rows"))
  message("Matched IDs (overlap): ", best$overlap)

  if (best$overlap < 10) {
    # Print a little debug context to make this self-service for students
    message("Example counts col IDs: ", paste(head(counts_cols, 5), collapse = ", "))
    message("Example counts row IDs: ", paste(head(counts_rows, 5), collapse = ", "))
    message("Example metadata IDs:   ", paste(head(canon_ids(meta_dt[[best$id_col]]), 5), collapse = ", "))
    stop("Too few overlapping IDs between counts and metadata after canonicalization.")
  }

  spot_ids <- canon_ids(meta_dt[[best$id_col]])

  meta_df <- as.data.frame(meta_dt[, setdiff(names(meta_dt), best$id_col), with = FALSE])
  rownames(meta_df) <- spot_ids

  # Make counts genes x spots
  if (best$orient == "rows") {
    message("Transposing counts to genes x spots.")
    counts_mat <- t(counts_mat)
    colnames(counts_mat) <- canon_ids(colnames(counts_mat))
  } else {
    colnames(counts_mat) <- counts_cols
  }

  # Align
  common_spots <- intersect(colnames(counts_mat), rownames(meta_df))
  if (length(common_spots) < 10) stop("Still too few overlapping spots after alignment; check files.")

  counts_mat <- counts_mat[, common_spots, drop = FALSE]
  meta_df <- meta_df[common_spots, , drop = FALSE]

  counts_sparse <- Matrix::Matrix(counts_mat, sparse = TRUE)

  seu_obj <- Seurat::CreateSeuratObject(
    counts = counts_sparse,
    meta.data = meta_df,
    min.cells = 3,
    min.features = 200
  )

  # Spatial coordinates in metadata
  if (all(c("new_x", "new_y") %in% colnames(seu_obj@meta.data))) {
    spatial_coords <- seu_obj@meta.data[, c("new_x", "new_y")]
    spatial_coords <- apply(spatial_coords, 2, as.numeric)
    rownames(spatial_coords) <- rownames(seu_obj@meta.data)
    colnames(spatial_coords) <- c("spatial_1", "spatial_2")

    seu_obj[["spatial"]] <- Seurat::CreateDimReducObject(
      embeddings = as.matrix(spatial_coords),
      key = "spatial_"
    )
  } else {
    message("No new_x/new_y found; skipping spatial reduction.")
  }

  # Standard pipeline
  seu_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu_obj, pattern = "^MT-")
  seu_obj <- Seurat::NormalizeData(seu_obj)
  seu_obj <- Seurat::FindVariableFeatures(seu_obj)
  seu_obj <- Seurat::ScaleData(seu_obj)
  seu_obj <- Seurat::RunPCA(seu_obj)
  seu_obj <- Seurat::FindNeighbors(seu_obj, dims = 1:10)
  seu_obj <- Seurat::FindClusters(seu_obj, resolution = 0.5)
  seu_obj <- Seurat::RunUMAP(seu_obj, dims = 1:10)

  seu_obj
}
