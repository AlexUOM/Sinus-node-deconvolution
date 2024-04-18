library(SingleCellExperiment)
library(Seurat)
library(anndata)
library(rhdf5)
library(data.table)
library(dplyr)

#Load data and annotations (SN Linscheid et al. 2020)
setwd('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Sinus Node')
sinus.file <- read.csv('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Sinus Node//Mouse scREF/GSE130710_normdata.csv')
annotations <- read.csv('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Sinus Node//Mouse scREF/GSE130710_metadata.csv')
sinus.file

# Remove the doublets
sinus.file <- sinus.file[,annotations$pANNPredictions == 'Singlet']
annotations <- annotations[,10] # Keep cell type column

#sc data and meta file names
sc_data_filename <- 'SN_sc_data.csv'
sc_meta_filename <- "SN_sc_meta.csv"

#-------------------------------------------------------------------------------
#Load data and annotations (Tabula Muris Heart)
setwd('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Right Atrium')
heart.file <- readRDS('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Right Atrium/Tabula Muris Heart/facs.normalized.Heart.rds')
annotations <- read.csv('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Right Atrium/Tabula Muris Heart/tabula-muris-senis-facs-official-raw-obj__cell-metadata__cleaned_ids.csv')
heart.file

#Filter for cells present in the RA
RA <- colData(heart.file)$subtissue == 'RA'
heart.file <- heart.file[,RA]

# Quick look at the number of cells in the RA
table(heart.file$free_annotation)

#Exclude ventricular cardiomyocytes
exclude_Ventric.CMs <- colData(heart.file)$free_annotation != "ventricular cardiomyocyte"
heart.file <- heart.file[,exclude_Ventric.CMs]

#Plot the number of cells in the RA
cells <- aggregate(data.frame(count = colData(heart.file)$free_annotation), list(value = colData(heart.file)$free_annotation), length)
barplot(height = cells$count, names = cells$value, las=2)

#sc data and meta file names
sc_data_filename <- 'RA_sc_data.csv'
sc_meta_filename <- "RA_sc_meta.csv"
#-------------------------------------------------------------------------------
## Tabula Sapiens Heart (Tabula Sapiens Consortium et al. 2022) ##
# Load the Human Heart scRNA-Seq data from h5ad (using anndata library)
setwd('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Validation deconvolution/TS_Heart/')
heart_file <- read_h5ad('TS_heart.h5ad') # anndata package needed

# Subset the AnnData object to only include cells from the atria
heart_file <- heart_file[(heart_file$obs$anatomical_information == "atria") | 
                            (heart_file$obs$anatomical_information == "Atria"), ]

table(heart_file$obs$free_annotation) #Check cell populations
heart_file <- heart_file[(heart_file$obs$free_annotation != "Hepatocyte"), ] #Discard Hepatocytes
annotations <- heart_file$obs$free_annotation # Extract cell type information

#sc data and meta file names
sc_data_filename <- 'RA_TS_sc_data.csv'
sc_meta_filename <- "RA_TS_sc_meta.csv"
#-------------------------------------------------------------------------------
## SAN cells from Kanemaru et al 2023 ##
# Load the Human Heart scRNA-Seq data from h5ad (using anndata library)
setwd('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Validation deconvolution/Litvinukova global heart/')
heart_file <- read_h5ad('Global_lognormalised.h5ad') # anndata package needed

# Filter the global dataset
heart_file <-  heart_file[(heart_file$obs$region == "SAN"), ] # Keep SAN cells
heart_file <-  heart_file[(heart_file$obs$cell_type != "Ventricular Cardiomyocyte"), ] # Discard the Ventricular CMs
heart_file <-  heart_file[(heart_file$obs$cell_state != "unclassified"), ] # Discard unclassified CMs


# Keep Central node cells as shown in Kanemaru et al. (From Extended Data 4a)
heart_file <- heart_file[(heart_file$obs$cell_state == "SAN_P_cell") | 
                           (heart_file$obs$cell_state == "NC2_glial_NGF+") |
                           (heart_file$obs$cell_state == "FB6") |
                           (heart_file$obs$cell_state == "PC3_str") |
                           (heart_file$obs$cell_state == "FB4_activated") |
                           (heart_file$obs$cell_state == "aCM3") |
                           (heart_file$obs$cell_state == "PC1_vent") , ]
annotations <- heart_file$obs$cell_state

#sc data and meta file names
sc_data_filename <- 'SAN_complete_sc_data.csv'
sc_meta_filename <- "SAN_complete_sc_meta.csv"
#-------------------------------------------------------------------------------
# Generate the sc_data
heart_file <- as.data.frame(as.matrix(heart_file$X))

sc_data <- transpose(heart_file) # data.table package needed 
rownames(sc_data) <- colnames(heart_file)
colnames(sc_data) <- rownames(heart_file)

counter <- 0
cell_vector <- character(0)
for (j in 0:length(colnames(sc_data))){
  cell_vector[counter] <- paste("C", j, sep = "_")
  counter <- counter+1
}

colnames(sc_data) <-  cell_vector

# Generate sc_meta
sc_meta <- as.data.frame(colnames(sc_data))
rownames(sc_meta) <- sc_meta[,1]
colnames(sc_meta)[1] <- "Cell"
sc_meta["Cell_type"] <- annotations

#-------------------------------------------------------------------------------
# USE FOR SN FROM LINSCHEID ET AL. ONLY!!
## Generate sc_data
sc_data <- as.data.frame(as.matrix(sinus.file)) #Remove "assay" function for SN
rownames(sc_data) <- sc_data$gene # Run for SN only
sc_data <- sc_data[,-1] # Run for SN only
rownames(annotations) <- annotations[,1] # Set the rownames as the cell_IDs

annotations <- annotations[colnames(sc_data),]

counter <- 0
cell_vector <- character(0)
for (j in 0:length(colnames(sc_data))){
  cell_vector[counter] <- paste("C", j, sep = "_")
  counter <- counter+1
}

colnames(sc_data) <-  cell_vector

#Generate sc_meta
sc_meta <- as.data.frame(colnames(sc_data))
rownames(sc_meta) <- sc_meta[,1]
colnames(sc_meta)[1] <- "Cell"
cell_type <- annotations[,8] # Extract the cell type from the appropriate column in annotations
sc_meta["Cell_type"] <- cell_type

#--------------------------------------------------------------------------------
#Write sc_data and sc_meta
write.csv(sc_data, file=sc_data_filename)
write.csv(sc_meta, file= sc_meta_filename)
