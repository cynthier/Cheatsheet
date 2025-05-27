###### 
library(ggpubr)

summ <- data.frame()

pdf("Expression_Ptn_Receptors.pdf", height = 16, width = 10)
for(st in c("E13.5", "E15.5")){
    genes <- c("Ptn",'Sdc2','Sdc3','Ncam1','Fgfr1','Itga5',"Itgb1", "Cdh2", "Cadm1", "Angpt2") 
    
    obj <- subset(objs, stage == st)
    obj <- ScaleData(obj, features = genes, assay = "RNA")
    exp.mat <- GetAssayData(obj, layer = "data", assay = "RNA")
    
    ################## Prepare data 
    cellanno <- data.frame(group = obj$group, 
                      row.names = colnames(obj),
                      cluster = obj$final.annotation)

    cellannoExp <- cbind(cellanno, t(as.matrix(exp.mat)[genes, rownames(cellanno)]))
    cellannoExpMat <- reshape2::melt(cellannoExp, id.vars = c("group", "cluster"))
    colnames(cellannoExpMat) <- c('group', 'cluster', 'gene', 'expression')
    

    #### expression matrix scale to 0-1 (gene level)
    cellannoExpMat2 <- data.frame()
    cellannoExpMat$scale <- 0
    for (i in unique(cellannoExpMat$gene)) {
        tmp.df <- subset(cellannoExpMat, gene == i)
        # tmp.df <- tmp.df %>% filter(expression < quantile(tmp.df$expression, 1))
        tmp.df$scale <- (tmp.df$expression - min(tmp.df$expression)) / (max(tmp.df$expression) - min(tmp.df$expression))
        tmp.df[is.na(tmp.df)] <- 0 
        cellannoExpMat2 <- rbind(cellannoExpMat2, tmp.df)
    }

    
    ##### text label for geom_text
    data_text <- data.frame()
    for (i in unique(cellannoExpMat2$cluster)) {
        for (j in unique(cellannoExpMat2$gene)) {
            direction <- "DOWN"
            mean.expression.ko <- mean(subset(cellannoExpMat2, 
                                  cluster == i & 
                                  group == "Vcl cKO" & 
                                  gene == j)$scale)
            mean.expression.control <- mean(subset(cellannoExpMat2, 
                                  cluster == i & 
                                  group == "Control" & 
                                  gene == j)$scale)
            
            log2fc <- log2(mean.expression.ko/mean.expression.control) 
            if(log2fc > 0){direction <- "UP"}
            # print(paste0("log2fc",log2fc, " ", direction))
            data_text <- rbind(data_text, c(i, j, mean.expression.control, 
                                            mean.expression.ko, log2fc,  direction))
        }
    }
    colnames(data_text) <- c("cluster", "gene", "mean.expression.control", "mean.expression.ko", "log2FC","label")
    data_text$gene <- factor(data_text$gene, levels = genes)
    data_text$cluster <- factor(data_text$cluster, levels = c('GP','BP','Neuroblast','BranchA','BranchB','ENMFB'))
    # data_text[data_text$gene %in% "Ptn",]
    
    data_text$group <- "test"
    data_text$test <- "Wilcoxon.test"
    data_text$alternative <- "two.sided"
    data_text$stage <- st

    summ <- rbind(summ, data_text)

    ################## plot data
    options(repr.plot.width = 10, repr.plot.height = 16)
    cellannoExpMat2$cluster <- factor(cellannoExpMat2$cluster, 
                                      levels = c('GP','BP','Neuroblast','BranchA','BranchB','ENMFB'))
    
    p1 <- ggplot(cellannoExpMat2, aes(x = group, y = scale, fill = group)) + 
        geom_violin(trim = FALSE, alpha = 0.8) + 
        stat_summary(fun = "mean", geom = "crossbar", width = 0.3, size = 0.5, color = "black") +
        geom_text(data = data_text, x = 1.5, y = 0.8, size = 3, 
                  mapping = aes(label = label),
                  color = ifelse(as.numeric(data_text$log2FC) < 0, 'blue', 'red')) +
        facet_grid(gene ~ cluster, switch = "y", scales = "free") +
        scale_y_continuous(position = "right", limits = c(-0.2, 1.2)) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) +
        labs(x = "", y = "", title = "") +
        stat_compare_means(label = "p.signif", 
                           method = "wilcox.test", 
                           label.x = 1.35, 
                           label.y = 1.0) +
        theme(strip.background = element_rect(fill = "gray97", color = NA), 
              strip.placement = "outside", 
              strip.text.x = element_text(face = "bold", 
                                          size = 10),
              strip.text.y = element_text(face = "italic", 
                                          size = 11)) +
        theme(axis.ticks.x = element_blank(), 
              axis.ticks = element_line(size = 0.1), 
              axis.text.x = element_text(face = "plain", 
                                         size = 10, 
                                         color = "black", 
                                         vjust = 0.5),
              axis.text.y = element_text(face = "plain", 
                                         size = 8, 
                                         color = "black"),
              axis.title = element_text(size = 12)) +
        scale_fill_manual(values = c("Vcl cKO" = "#f82401", "Control" = "#0f7dbd")) &
        ggtitle(st)
    
    print(p1)  
}
dev.off()
