BiocManager::install("ComplexHeatmap")  # required for CellChat
remotes::install_github("sqjin/CellChat")
#### Load packages ####
library(Seurat) # general single-cell processing
library(Matrix) # sparse matrix
library(readr)
library(tibble)
library(dplyr)
library(CellChat)

Sys.setenv(VROOM_CONNECTION_SIZE = 131072 * 2)

#### Load data ####
expression_matrix <- file.path(
  "single-cell",
  "data",
  "GSE109816_normal_heart_umi_matrix.csv.gz"
)
metadata_path <- file.path(
  "single-cell",
  "data",
  "GSE109816_normal_heart_cell_cluster_info.txt"
)

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

## Increment cluster labels by 1 because cellchat hates 0's
seu_obj$seurat_clusters <- as.factor(as.numeric(seu_obj$seurat_clusters) + 1)

seu_obj <- RunUMAP(seu_obj, dims = 1:10)
DimPlot(seu_obj, reduction = "umap", label = TRUE)


#### Create a CellChat object ####
cellchat <- createCellChat(object = seu_obj, group.by = "seurat_clusters")

# Set the default database (Secreted Signaling or ECM-Receptor datasets)
# there is also mouse and zebrafish databases available
# Assign the database to the CellChat object
cellchat@DB <- CellChatDB.human  # If working with human data


#### Identify Overexpressed Genes and Ligand-Receptor Pairs ####
# Subset expression data to overexpressed genes
cellchat <- subsetData(cellchat) # Filter out genes based on their expression
# Identify overexpressed ligands/receptors
cellchat <- identifyOverExpressedGenes(cellchat)
# Identify overexpressed interactions between cell types
cellchat <- identifyOverExpressedInteractions(cellchat)


#### Calculate Communication Probabilities ####
# Calculate communication probabilities
cellchat <- computeCommunProb(cellchat)

# Filter out communications with low probabilities
cellchat <- filterCommunication(cellchat, min.cells = 10)


#### Infer the communication network at the pathway level ####
# Infer the communication network at the pathway level
cellchat <- computeCommunProbPathway(cellchat)
# Aggregate the inferred cell-cell communications
cellchat <- aggregateNet(cellchat)


#### Visualize the communication network ####
# Plot the aggregated network
netVisual_circle(
  cellchat@net$count,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale = TRUE,
  label.edge = TRUE
)

# Chord diagram for a specific signaling pathway, e.g., "WNT"
netVisual_chord_gene(cellchat, signaling = "VEGF")
# Visualize contribution of ligands/receptors in a specific signaling pathway
netAnalysis_contribution(cellchat, signaling = "VEGF")


# Heatmap for all signaling pathways
netVisual_heatmap(cellchat)
# Hierarchical plot for a specific signaling pathway, e.g., "WNT"
cellchat <- netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)
netAnalysis_signalingRole_network(
  cellchat,
  signaling = "VEGF",
  slot.name = "netP"
)


# Heatmap showing the role of different cell types as senders or receivers
# View both outgoing and incoming roles for VEGF pathway
netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing"
) # For outgoing (sending)
netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming"
) # For incoming (receiving)
