library(data.table, quietly = TRUE)
library(Seurat, quietly = TRUE) # general single-cell processing
library(Matrix, quietly = TRUE) # sparse matrix
library(tibble, quietly = TRUE) # data manipulation


fetal_human_heart <- function(expression_matrix, metadata_path) {
  #### Load data ####
  data <- fread(expression_matrix) %>%
    tibble::column_to_rownames("V1") %>%
    as.matrix() %>%
    Matrix(sparse = TRUE)
  print(head(data, 1))

  metadata <- fread(metadata_path) %>%
    tibble::column_to_rownames("V1")
  print(head(metadata, 1))

  seu_obj <- CreateSeuratObject(
    counts = data,
    meta.data = metadata,
    min.cells = 3,
    min.features = 200
  )

  # spatial data is in 'new_x' and 'new_y' in the metadata:
  spatial_coords <- metadata[, c("new_x", "new_y")]
  # Rename columns to match the key
  colnames(spatial_coords) <- paste0("spatial_", seq_len(ncol(spatial_coords)))
  print(head(spatial_coords))

  seu_obj[["spatial"]] <- CreateDimReducObject(
    embeddings = as.matrix(spatial_coords),
    key = "spatial_"
  )



  #### Standard Seurat pipeline ####
  seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
  seu_obj <- NormalizeData(
    seu_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  seu_obj <- FindVariableFeatures(
    seu_obj,
    selection.method = "vst",
    nfeatures = 2000
  )
  seu_obj <- ScaleData(seu_obj)
  seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
  seu_obj <- FindNeighbors(seu_obj, dims = 1:10)
  seu_obj <- FindClusters(seu_obj, resolution = 0.5)

  seu_obj <- RunUMAP(seu_obj, dims = 1:10)
}

matrix_data <- file.path(
  "spatial-transcriptomics",
  "data",
  "fetal_human_heart",
  "filtered_matrix.tsv.gz"
)
meta_data <- file.path(
  "spatial-transcriptomics",
  "data",
  "fetal_human_heart",
  "meta_data.tsv.gz"
)
x_test <- fetal_human_heart(matrix_data, meta_data)
