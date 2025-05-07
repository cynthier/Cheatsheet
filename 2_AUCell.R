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
