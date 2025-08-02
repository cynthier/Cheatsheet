library(CytoTRACE2)
result <- cytotrace2(counts(sce_ss), 
                     species = "human", 
                    is_seurat = FALSE,
                    slot_type = "counts")

plot.data <- meta %>% group_by(subcluster) %>% filter(CytoTRACE2_Score <= quantile(CytoTRACE2_Score, 1, na.rm = TRUE)) %>% ungroup()
options(repr.plot.width = 8, repr.plot.height = 4)
comparison <- list(c("Control", "Mutant"))

my_viol_compare.plt(meta = plot.data,
                    x = "sample",
                    y = "CytoTRACE2_Score", 
                    color = "sample", 
                    title = "Potency score",
                    ylab = "Potency score", 
                    facet = FALSE, 
                    comparison = comparison) & ylim(0, 0.8)

my_viol_compare.plt <- function(
                               meta = meta,
                               x = x, 
                               y = y, 
                               color = color, 
                               title = "cytoTRACE", 
                               ylab = "Potency score", 
                               facet = TRUE, 
                               facet.var = NULL,
                               comparison = NULL){ 
    library(ggplot2)
    library(ggpubr)
    
    colors <- readRDS("/home/lfliuhku/cheatsheet/colors.rds")
    if(facet){ 
        p <- ggplot(meta, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) + 
            geom_violin() + 
            stat_summary(fun = "mean", geom = "crossbar", width = 0.3, size = 0.5, aes(color = .data[[color]]))+ 
            facet_grid(.~.data[[facet.var]]) + 
            scale_color_manual(values = colors) +
            stat_compare_means(method = "wilcox.test", size = 4, comparison = comparison) +
            theme_bw(base_size = 16) + 
            theme(plot.title = element_text(hjust = 0.5), 
                  legend.position = "right", 
                  axis.text.x = element_blank()) + 
            ggtitle(title) + labs(x = NULL, y = ylab)
            return(p)  
    }else{ 
        p <- ggplot(meta, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) + 
        geom_violin() + 
        stat_summary(fun = "mean", geom = "crossbar", width = 0.3, size = 0.5, aes(color = .data[[color]]))+ 
        # facet_grid(.~.data[[facet.var]]) + 
        scale_color_manual(values = colors) +
        stat_compare_means(method = "wilcox.test", size = 4, comparison = comparison) +
        theme_bw(base_size = 16) + 
        theme(plot.title = element_text(hjust = 0.5), 
              legend.position = "right", 
              axis.text.x = element_blank()) + 
        ggtitle(title) + labs(x = NULL, y = ylab)
        return(p)    
    }
}
