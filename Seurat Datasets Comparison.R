library(Seurat)
library(patchwork)
library(metap)
library(ggplot2)

# Load data
# RA data
countData <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Right Atrium/output/RA_human_from_mouse_sc_data.csv", header = TRUE, row.names =  1)
cell_type <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Right Atrium/output/RA_human_from_mouse_sc_celltype.csv", header = TRUE)
rownames(cell_type) <- cell_type[,1]
ra_cells <- CreateSeuratObject(counts = countData, project = "SC_project", min.cells = 3, min.features = 200, meta.data = cell_type)
ra_cells <- AddMetaData(ra_cells, metadata = "Right_Atrium", col.name = "Tissue")


# SN data
countData <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Sinus Node/output/SN_human_from_mouse_sc_data.csv", header = TRUE, row.names =  1)
cell_type <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Sinus Node/output/SN_human_from_mouse_sc_celltype.csv", header = TRUE)
rownames(cell_type) <- cell_type[,1]
sn_cells <- CreateSeuratObject(counts = countData, project = "SC_project", min.cells = 3, min.features = 200, meta.data = cell_type)
sn_cells <- AddMetaData(sn_cells, metadata = "Sinus_Node", col.name = "Tissue")


# Merging the two datasets
datasets.list <- list(ra_cells, sn_cells)

# Normalize and identify variable features for each dataset independently
datasets.list <- lapply(X = datasets.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Integration feature selection
features <- SelectIntegrationFeatures(object.list = datasets.list)

# Anchors identification and dataset integration
cardiac.anchors <- FindIntegrationAnchors(object.list = datasets.list, anchor.features = features)
cardiac.combined <- IntegrateData(anchorset = cardiac.anchors)

# Working of integrated corrected data
DefaultAssay(cardiac.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
cardiac.combined <- ScaleData(cardiac.combined, verbose = FALSE)
cardiac.combined <- RunPCA(cardiac.combined)

#Determine the 'dimensionality' of the dataset
cardiac.combined <- JackStraw(cardiac.combined, num.replicate = 100)
cardiac.combined <- ScoreJackStraw(cardiac.combined, dims = 1:20)

#Visualise with Straw plot and Elbow plot
JackStrawPlot(cardiac.combined, dims = 1:20)
ElbowPlot(cardiac.combined, ndims = 20)

cardiac.combined <- FindNeighbors(cardiac.combined, reduction = "pca", dims = 1:12)
cardiac.combined <- FindClusters(cardiac.combined, resolution = 0.2)


# Visualization using UMAP
cardiac.combined <- RunUMAP(cardiac.combined, reduction = "pca", dims = 1:12)
p1 <- DimPlot(cardiac.combined, reduction = "umap", group.by = "Tissue", label.size = 5)
p2 <- DimPlot(cardiac.combined, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5)
p1 + p2
umap <- p1+p2
ggsave(umap, file = "UMAP 2 plots.png", width = 13, height = 5, units = 'in', dpi = 600 )

# Visualization using TSNE
cardiac.combined <- RunTSNE(cardiac.combined, reduction = "pca", dims = 1:12)
p3 <- DimPlot(cardiac.combined, reduction = "tsne", group.by = "Tissue", order = c('Right Atrium','Sinus Node'), label.size = 5)
p4 <- DimPlot(cardiac.combined, reduction = "tsne", label = TRUE, repel = FALSE, label.size = 5)
tsne <- p3 + p4
tsne
ggsave(tsne, file = "tSNE 2 plots.png", width = 13, height = 5, units = 'in', dpi = 600 )