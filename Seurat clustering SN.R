library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

#Loading data
countData <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Sinus Node/output/SN_human_from_mouse_sc_data.csv", header = TRUE, row.names =  1)
cell_type <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Sinus Node/output/SN_human_from_mouse_sc_celltype.csv", header = TRUE)
cell_type <- cell_type[,2:3]
rownames(cell_type) <- cell_type[,1]
sn_cells <- CreateSeuratObject(counts = countData, project = "SN_project", min.cells = 3, min.features = 200, meta.data = cell_type)

#Calculate percentage of mitochondrial genes in each cell
sn_cells[["percent.mt"]] <- PercentageFeatureSet(sn_cells, pattern = "MT-")
VlnPlot(sn_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualise feature-feaature relationships
plot1 <- FeatureScatter(sn_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sn_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Remove cells that have a given number of genes OR
#more than a certain % of mitochondrial genes
sn_cells <- subset(sn_cells, subset = nFeature_RNA > 200 & percent.mt < 25)

#Normalize data with LogNormalize
sn_cells <- NormalizeData(sn_cells, normalization.method = "LogNormalize", scale.factor = 1000)

#Find highly variable genes 
sn_cells <- FindVariableFeatures(sn_cells, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sn_cells), 10)
top10

# plot variable features with and without labels
plot <- VariableFeaturePlot(sn_cells)
plot1 <- LabelPoints(plot = plot, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1

#Scale data
all.genes <- rownames(sn_cells)
sn_cells <- ScaleData(sn_cells, features = all.genes)

#Perform linear dimension reduction...
sn_cells <- RunPCA(sn_cells, features = VariableFeatures(object = sn_cells))

#...And visualise PCA in some ways
VizDimLoadings(sn_cells, dims = 1:2, reduction = "pca")
DimPlot(sn_cells, reduction = "pca")
DimHeatmap(sn_cells, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sn_cells, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the 'dimensionality' of the dataset
sn_cells <- JackStraw(sn_cells, num.replicate = 100)
sn_cells <- ScoreJackStraw(sn_cells, dims = 1:20)

#Visualise with Straw plot
JackStrawPlot(sn_cells, dims = 1:20)

#Or with Elbow plot
ElbowPlot(sn_cells)

#Cluster cells
sn_cells <- FindNeighbors(sn_cells, dims = 1:12)
sn_cells <- FindClusters(sn_cells, resolution = 0.2)


#Run UMAP/tSNEs
sn_cells <- RunUMAP(sn_cells, umap.method = 'umap-learn', metric = 'correlation', dims = 1:12)
sn_cells <- RunTSNE(sn_cells, dims = 1:12)
p1 <- DimPlot(sn_cells, reduction = "umap", label = T, label.box = F, repel =T)
p2 <- DimPlot(sn_cells, reduction = "tsne", label = F, label.box = F, repel= T)
p1 +p2
p2
# Save tSNE or UMAP
ggsave(p2, file = "tSNE Plot labelled.png", width =6 , height = 5, units = 'in', dpi = 600, bg= "white" )

# Find markers for every cluster compared to all remaining cells, report only the positive
# ones
sn_cells.markers <- FindAllMarkers(sn_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sn_cells.markers %>% #Load dplyer for pipe operator
  group_by(cluster) %>%
  slice_max(n = 1000, order_by = avg_log2FC)%>% 
  write.csv(file='cell_markers_res0.2_dims12_labelled.csv')

sn_cells.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(sn_cells, features = top10$gene) + NoLegend()

# Annotate clusters with cell type names
new.cluster.ids <- c(`0` = "Adipocytes 1", # 0
                     `1` ="Fibroblasts", # 1
                     `2` ="Sinus node myocytes 1", # 2
                     `3` ="Sinus node myocytes 2", # 3
                     `4` ="Adipocytes 2", # 4
                     `5` ="Macrophages", # 5
                     `6` ="Endothelial cells", # 6
                     `7` ="Pacemaker myocytes", # 7
                     `8` ="Lymphoid cells 1", # 8
                     `9` ="Lymphoid cells 3", # 9
                     `10` ="Lymphoid cells 2") # 10



sn_cells <- RenameIdents(sn_cells, new.cluster.ids)

labelled <- DimPlot(sn_cells, reduction = "tsne", label = T, pt.size = 0.2, repel = F)
labelled
ggsave(labelled, file = "tSNE Plot labelled.png", width =7 , height = 5, units = 'in', dpi = 600, bg= "white" )


# Dotplot
markers.to.plot <- c('HCN1', 'HCN4', 'CACNA1D', 'SHOX2', 'TBX3', #P CMs
                     'MYH6', 'KCNJ3', 'CTNNA3', 'RYR2', 'TBX5', # Remote CMs
                     'ICAM1','PECAM1', 'VWF', 'CDH5', 'ENG', # EC
                     'GPAM', 'PPARG', 'FASN', 'LEP', 'PCK1','PLIN1', # Adipo
                     'CD163', 'C1QA', 'LYVE1', 'IGF1', # Macro
                     'ALCAM', 'CD84', 'ATP8B4', 'ITGA4', 'CD44', #Lymphoid
                     'IKZF3', 'BANK1', 'CARD11', 'PRKCB', #Lymphoid 1
                     'COL4A1', 'COL4A2', 'COL6A6', 'COL6A3','FBN1') #Fibroblasts
                      
sn_cells@active.ident <- factor(sn_cells@active.ident,
                                levels = c("Pacemaker myocytes",
                                           "Sinus node myocytes 1",
                                           "Sinus node myocytes 2",
                                           "Endothelial cells",
                                           "Adipocytes 1",
                                           "Adipocytes 2",
                                           "Macrophages",
                                           "Lymphoid cells 1",
                                           "Lymphoid cells 2",
                                           "Lymphoid cells 3",
                                           "Fibroblasts"))

dotplot <- DotPlot(sn_cells, features = markers.to.plot, col.min = 0, dot.scale = 8) +
  theme(axis.text = element_text(size = 16)) +
  theme(legend.text = element_text(size = 13.5)) +
  RotatedAxis()
dotplot

# Saving the plot with good resolution
ggsave(dotplot, file = "Markers Dotplot SN.png", width = 15, height = 6, units = 'in', dpi = 600, bg= 'white' )

# DGE between the three SN cardiomyocyte populations
# Step 1: calculate DGE
pacemaking_vs_sn_cms.markers <- FindMarkers(sn_cells, ident.1 = 'Pacemaker myocytes', ident.2 = c('Sinus node myocytes 1','Sinus node myocytes 2'), logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
pacemaking_vs_all.markers <- FindMarkers(sn_cells, ident.1 = 'Pacemaker myocytes', logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)

sn_cms1.markers <- FindMarkers(sn_cells, ident.1 = 'Sinus node myocytes 1', ident.2 = c('Pacemaker myocytes','Sinus node myocytes 2'), logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
sn_cms2.markers <- FindMarkers(sn_cells, ident.1 = 'Sinus node myocytes 2', ident.2 = c('Sinus node myocytes 1','Sinus node myocytes 1'), logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)

# Step 2: Export results
write.csv(pacemaking_vs_all.markers, file = 'Pacemaker CMs DGE.csv')
write.csv(sn_cms1.markers, file = 'SN CMs 1 DGE.csv')
write.csv(sn_cms2.markers, file = 'SN CMs 2 DGE.csv')

