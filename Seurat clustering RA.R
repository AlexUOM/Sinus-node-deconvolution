library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

#Loading data
countData <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Right Atrium/output/RA_human_from_mouse_sc_data.csv", header = TRUE, row.names =  1)
cell_type <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Right Atrium/output/RA_human_from_mouse_sc_celltype.csv", header = TRUE)
cell_type <- cell_type[,2:3]
rownames(cell_type) <- cell_type[,1]
ra_cells <- CreateSeuratObject(counts = countData, project = "ra_project", min.cells = 3, min.features = 200, meta.data = cell_type)

#Calculate percentage of mitochondrial genes in each cell
ra_cells[["percent.mt"]] <- PercentageFeatureSet(ra_cells, pattern = "MT-")
VlnPlot(ra_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualise feature-feaature relationships
plot1 <- FeatureScatter(ra_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ra_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Remove cells that have a given number of genes OR
#more than a certain % of mitochondrial genes
ra_cells <- subset(ra_cells, subset = nFeature_RNA > 200 & percent.mt < 25)

#Normalize data with LogNormalize
ra_cells <- NormalizeData(ra_cells, normalization.method = "LogNormalize", scale.factor = 10000)

#Find highly variable genes 
ra_cells <- FindVariableFeatures(ra_cells, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ra_cells), 10)
top10
# plot variable features with and without labels
plot <- VariableFeaturePlot(ra_cells)
plot1 <- LabelPoints(plot = plot, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1

#Scale data
all.genes <- rownames(ra_cells)
ra_cells <- ScaleData(ra_cells, features = all.genes)

#Perform linear dimension reduction...
ra_cells <- RunPCA(ra_cells, features = VariableFeatures(object = ra_cells))

#...And visualise PCA in some ways
VizDimLoadings(ra_cells, dims = 1:2, reduction = "pca")
DimPlot(ra_cells, reduction = "pca")
DimHeatmap(ra_cells, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(ra_cells, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the 'dimensionality' of the dataset
ra_cells <- JackStraw(ra_cells, num.replicate = 100)
ra_cells <- ScoreJackStraw(ra_cells, dims = 1:20)

#Visualise with Straw plot
JackStrawPlot(ra_cells, dims = 1:20)

#Or with Elbow plot
ElbowPlot(ra_cells)

#Cluster cells
ra_cells <- FindNeighbors(ra_cells, dims = 1:12)
ra_cells <- FindClusters(ra_cells, resolution = 0.1)

#Run UMAP/tSNEs
ra_cells <- RunUMAP(ra_cells, umap.method = 'umap-learn', metric = 'correlation', dims = 1:12)
ra_cells <- RunTSNE(ra_cells, dims = 1:12)
p1 <- DimPlot(ra_cells, reduction = "umap", label = T, label.box = F, repel =T)
p2 <- DimPlot(ra_cells, reduction = "tsne", label = F, label.box = F, repel= T)
p1 +p2
p2
# Save tSNE or UMAP
ggsave(p2, file = "tSNE Plot.png", width =12 , height = 5, units = 'in', dpi = 600, bg= "white" )

# Find markers for every cluster compared to all remaining cells, report only the positive
# ones
ra_cells.markers <- FindAllMarkers(ra_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ra_cells.markers %>% #Load dplyer for pipe operator
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC)%>% 
  write.csv(file='cell_markers_res0.1_dims12.csv')
cluster0.markers <- FindMarkers(ra_cells, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Heathmap of marker genes of each cluster
ra_cells.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ra_cells, features = top10$gene) + NoLegend()

# Final tSNE plot with labelled names
new.cluster.ids <- c('Pericytes 1', "Fibroblasts 1", "Endothelial cells 1", 
                     "Atrial myocytes 1", "Smooth muscle cells", "Myeloid cells 1", "Endocardial cells",
                     "Atrial myocytes 2", "Pericytes 2", "Myeloid cells 2",
                     "Myeloid cells 3", "Atrial myocytes 3", "Myeloid cells 4",
                     "Atrial myocytes 4", 'Endothelial cells 2', "Myeloid cells 5",
                     "Lymphoid cells", "Fibroblasts 2")

names(new.cluster.ids) <- levels(ra_cells)
new.cluster.ids <- sort(new.cluster.ids)
ra_cells <- RenameIdents(ra_cells, new.cluster.ids)
labeled_tsne <- DimPlot(ra_cells, reduction = "tsne", label = T, pt.size = 0.5, repel = F)
labeled_tsne <- labeled_tsne + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) +
  theme(axis.text = element_text(size = 19))+
  theme(legend.text = element_text(size = 16))
labeled_tsne
ggsave(labeled_tsne, file = "labelled tSNE.png", width = 7.5, height = 5.5, units = 'in', dpi = 600, bg= 'white' )

# Dotplot of marker genes
markers.to.plot <- c("TNNT2", "ACTC1", 'NPPA',
                     "PLVAP", "EMCN", 'PECAM1',
                     "FABP4", "MGLL",
                     "LGALS3", "PLAC8",
                     "DCN","GSN",
                     "CD83", "LTB",
                     "KCNJ8", "RGS5",
                     "CD74", "CD14",
                     "RAC2", "ROGDI",
                     "IL1B", "HDC",
                     "HBEGF", "GSTO1",
                     "C1QB", "C1QA",
                     "ACTA2", "MYL9")
                                     
new.order <- c("Atrial myocytes 1", "Atrial myocytes 2", "Atrial myocytes 3",
               "Atrial myocytes 4", "Endocardial cells", "Endothelial cells 1",
               "Endothelial cells 2", "Fibroblasts 1", "Fibroblasts 2", 
               "Lymphoid cells", "Pericytes 1", "Pericytes 2", "Myeloid cells 1",
               "Myeloid cells 2", "Myeloid cells 3", "Myeloid cells 4", "Myeloid cells 5", 
               "Smooth muscle cells")
levels(ra_cells) <- new.order
dotplot <- DotPlot(ra_cells, features = markers.to.plot, idents = cell.order, col.min = 0, dot.scale = 9) +
  theme(axis.text = element_text(size = 16)) +
  RotatedAxis()
dotplot

# Saving the plot with good resolution
ggsave(dotplot, file = "Markers Dotplot.png", width = 14, height = 8, units = 'in', dpi = 600, bg= 'white' )
