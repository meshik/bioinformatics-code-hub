#### Install packages ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Seurat", "monocle3", "Matrix"))

#### Load packages ####
library(Seurat)  # general single-cell processing
library(monocle3)  # pseudotime analysis
library(Matrix)  # sparse matrix
library(readr)

#### Load data ####
file.path("single-cell", "data", "GSE109816_normal_heart_umi_matrix.csv.gz")
data <- readr::read_csv(
  "single-cell\\data\\GSE109816_normal_heart_umi_matrix.csv.gz") %>%
  tibble::column_to_rownames("...1") %>%
  as.matrix() %>%
  Matrix(sparse = TRUE)

metadata <- readr::read_tsv(paste0(BASE_DIR, SC_METADATA_PATH)) %>%
  tibble::column_to_rownames("ID")

seu_obj <- CreateSeuratObject(counts = data, meta.data = metadata,
                              min.cells = 3, min.features = 200)

seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj,pattern = "^MT-")

return(seu_obj)