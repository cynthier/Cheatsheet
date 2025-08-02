pdf(paste(unique(genes$term), "pathway_score.pdf"), height = 6, width = 6)
for(i in 1:nrow(datalist)){ 
    name <- paste0(datalist[i, "stage"], "_", datalist[i, "celltype"])
    obj.temp <- subset(obj, final.annotation == datalist[i, "celltype"] & 
                       stage == datalist[i, "stage"])
    count.mat <- GetAssayData(object = obj.temp, slot = "counts", assay = "RNA")
    cells.rankings <- AUCell_buildRankings(count.mat, plotStats = FALSE)

    cells_AUC <- AUCell_calcAUC(list("geneset" = genes$SYMBOL), cells.rankings)
    auc.mat <- getAUC(cells_AUC)
    obj.temp$pathway_score <- auc.mat[1,colnames(obj.temp), drop = T]
    
    comparison <- list(c("Control", "Vcl"))
    obj.temp@meta.data$group <- factor(obj.temp@meta.data$group, 
                                   levels = c("Control", "Vcl cKO"))

    p <- ggboxplot(obj.temp@meta.data, 
              x = "group", 
              y = "pathway_score",
              color = "group", 
              palette = "jco",
              add = "jitter", 
              add.params = list(alpha = 0.2),  
              title = paste0(name, "\n", unique(genes$term), "\n", unique(genes$GOALL))) + 
        stat_compare_means(method = "wilcox.test") +
        theme_bw(base_size = 14) + 
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none") + 
        ylab("Pathway score") + xlab(NULL)
    print(p)
        
    #### determine thresholds to allocate the cells
    # cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign=TRUE)
    # thr = cells_assignment$geneSet$aucThr$selected;thr
    # new_cells <- names(which(getAUC(cells_AUC)["geneset",]> thr))
     
}
dev.off()


###############################################

library(Seurat)
library(AUCell)
library(ggplot2)
library(ggpubr)
colors <- readRDS("/home/lfliuhku/cheatsheet/colors.rds")

count.mat <- GetAssayData(object = obj, slot = "counts", assay = "RNA")
meta <- obj@meta.data
i <- names(genes)

pdf( "pathway_score.pdf", height = 4, width = 12)
options(repr.plot.width = 12, repr.plot.height = 4)
for(i in names(genes)){ 

    cells.rankings <- AUCell_buildRankings(count.mat, plotStats = FALSE)
    cells_AUC <- AUCell_calcAUC(list("geneset" = genes[[i]]), cells.rankings)
    auc.mat <- getAUC(cells_AUC)
    meta[, i] <- auc.mat[1,rownames(meta), drop = T]
    
    comparison <- list(c("Control", "Mutant"))
    meta$group <- factor(meta$group, levels = c("Control", "Mutant"))
    p <- ggplot(meta, aes(x = group, y = .data[[i]], color = group)) + 
        geom_violin() + 
        stat_summary(fun = "mean", geom = "crossbar", width = 0.3, size = 0.5, aes(color = group))+ 
        facet_grid(.~subcluster) + 
        scale_color_manual(values = colors) +
        stat_compare_means(method = "wilcox.test", size = 3) +
        theme_bw(base_size = 14) + 
        theme(plot.title = element_text(hjust = 0.5), 
              legend.position = "right", 
              axis.text.x = element_blank()) + 
        ggtitle(pathinfo[i, "name"]) + labs(x = NULL, y = "Pathway score")
    print(p)  
}
dev.off()


