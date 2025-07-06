
####################### Dotplot marker genes
options(repr.plot.width = 9, repr.plot.height = 4)
DotPlot(subset(obj, doublet == "Singlet"), 
        features = unique(genes), 
        group.by = "celltype") & 
    theme(axis.text = element_text(size = 14, 
                                   color = "black"), 
          axis.text.x = element_text(face = "italic", 
                                     angle = 30, 
                                     hjust = 1),
          axis.text.y = element_text(), 
          axis.title = element_text(size = 16), 
          plot.title = element_text(size = 18, 
                                    face = "bold", 
                                    hjust = 0.5)) & 
    ylab(NULL) & xlab("Marker genes") & xlab(NULL) &
    scale_size(range = c(0.1,8)) & 
    ggtitle("Marker genes expression")
