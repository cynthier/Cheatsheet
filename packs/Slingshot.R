library(slingshot)
library(SingleCellExperiment)
library(Seurat)
obj <- readRDS("./outputdata/obj_normalized_con_ko_final.rds")

#### Run slingshot
sce <- as.SingleCellExperiment(obj, assay = "RNA") 
sce_ss <- slingshot(sce,      
                     reducedDim = obj@reductions$pca@cell.embeddings[,c(1,3)], 
                     clusterLabels = sce$final.annotation, 
                     start.clus = 'BP', 
                     end.clus = "BranchA",
                     approx_points = 150)
obj$pseudotime <- sce_ss@colData[colnames(obj), "slingPseudotime_1"]

# saveRDS(obj@meta.data, file = "./outputdata/pseudotime_normalized.rds")

#### Get metadata
meta <- reducedDim(sce_ss, type = "PCA", withDimnames = TRUE) %>% as.data.frame()
meta$pseudotime <- sce_ss@colData[rownames(meta), "slingPseudotime_1"]
meta$stage <- sce_ss@colData[rownames(meta), "stage"]
meta$group <- sce_ss@colData[rownames(meta), "group"]
meta$anno <- sce_ss@colData[rownames(meta), "final.annotation"]

#### plot pseudotime in dimensions
options(repr.plot.width = 7, repr.plot.height = 5)
p <- ggplot(meta, aes(x = PC_1, y = PC_3, color = pseudotime)) + 
    geom_point(size = 0.8) + scale_color_gradientn(colours = c(colorRampPalette(c("#5b51a3","#79c9a4","#f2faac","#fdb465","#a4104d"))(90))) + 
    theme_bw(base_size = 16) + 
    theme(axis.text = element_text(size = 14), 
          axis.title =  element_text(size = 14), 
          strip.text = element_text(size = 14)) & 
    xlab("PC_1") &  ylab("PC_3"); print(p)

#### plot marker expression
options(repr.plot.width = 10, repr.plot.height = 8)
dims <- c(1,3)
FeaturePlot(obj, features = c("Sox10", "Ube2c", "Tubb3", "Cartpt"),reduction = "pca", dims = dims, pt.size = 0.001, order = T, ncol = 2)


#### density plot for the pseudotime
options(repr.plot.width = 6, repr.plot.height = 6)
ggplot(meta, aes(x = pseudotime, fill = group)) + 
  geom_density(alpha = 0.3) + 
  scale_fill_manual(values = colors) + 
  facet_grid(stage~.) & theme_bw() & 
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14))
