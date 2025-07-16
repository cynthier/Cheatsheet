############# read Loom file
library(loomR)
library(SeuratDisk)

data <- Connect("/home/rawdata/Public/scRNA//GSM4504449_E18_10X_18_065.loom", mode = "r") 
data.ser <- as.Seurat(data)


########### convert h5ad to seurat object
library(Seurat)
library(SeuratData)
library(SeuratDisk)

library(SeuratDisk)
Convert("/home/rawdata/Public/scRNA/final_fetal_object_cellxgene.h5ad", dest = "h5seurat", overwrite = F)
obj <- LoadH5Seurat("final_fetal_object_cellxgene.h5seurat")
