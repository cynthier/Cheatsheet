############# read Loom file
library(loomR)
library(SeuratDisk)

data <- Connect("/home/rawdata/Public/scRNA//GSM4504449_E18_10X_18_065.loom", mode = "r") 
data.ser <- as.Seurat(data)


########### convert h5ad to seurat object
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# library(rhdf5) ###### have a look on the data
# h5f <- H5Fopen("/home//rawdata/Public/scRNA/human_fetal_gut/ENS_log_counts.h5seurat")
# # h5ls(h5f, recursive = TRUE)
# H5Fclose(h5f)library(SeuratDisk)

Convert("/home/rawdata/Public/scRNA/final_fetal_object_cellxgene.h5ad", dest = "h5seurat", overwrite = F)
obj <- LoadH5Seurat("final_fetal_object_cellxgene.h5seurat")
