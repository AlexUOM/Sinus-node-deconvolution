library(anndata)

## Global heart scRNA-Seq (Litvinuková et al. 2020) ##
# Load the Human RA scRNA-Seq Ref from h5ad (using anndata library)
setwd('C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Validation deconvolution/Kanemaru human SN scRNA-Seq/')
file <- 'C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Halina Dobrzynski/Validation deconvolution/Litvinukova global heart/Global_lognormalised.h5ad'
heart_file <- read_h5ad(file) # anndata package needed

# Filer the global dataset
heart_file <-  heart_file[(heart_file$obs$region == "SAN"), ] # Keep SAN cells
heart_file <-  heart_file[(heart_file$obs$cell_type != "Ventricular Cardiomyocyte"), ] # Discard the 2 Ventricular CMs

# Save data and meta files
write.csv(as.matrix(heart_file$X), file='SN_only_data.csv')
write.csv(heart_file$obs$cell_type, file='SN_only_meta.csv')
