BiocManager::install("ComplexHeatmap")  # required for CellChat
install.packages("remotes")
remotes::install_github("sqjin/CellChat")

#### Load packages ####
library(dplyr)
library(CellChat)

#### Load data ####
source(file.path(
  "single-cell",
  "preprocessing",
  "load_and_prepare_seurat.R"
))
seu_obj <- adult_heart()

## Increment cluster labels by 1 because cellchat hates 0's
# Ideally you should 
seu_obj$seurat_clusters <- as.factor(as.numeric(seu_obj$seurat_clusters) + 1)
DimPlot(seu_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

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
