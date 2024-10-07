Sys.setenv(VROOM_CONNECTION_SIZE = 131072 * 2)

library(readr, quietly = TRUE)
library(Seurat, quietly = TRUE) # general single-cell processing
library(Matrix, quietly = TRUE) # sparse matrix
library(tibble, quietly = TRUE) # data manipulation

adult_heart <- function(expression_matrix, metadata_path) {

  #### Load data ####
  data <- readr::read_csv(expression_matrix) %>%
    tibble::column_to_rownames("...1") %>%
    as.matrix() %>%
    Matrix(sparse = TRUE)

  metadata <- readr::read_tsv(metadata_path) %>%
    tibble::column_to_rownames("ID")

  seu_obj <- CreateSeuratObject(
    counts = data,
    meta.data = metadata,
    min.cells = 3,
    min.features = 200
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
