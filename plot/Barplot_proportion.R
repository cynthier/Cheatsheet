colors <- readRDS("/home/lfliuhku/cheatsheet/colors.rds")

plot.var1 <- "group"
plot.var2 <- "subcluster"
levels.var2 <- c('new C3','new C2','new C5','new C4','new C1','new C0')
titlename <- "HDAC newcluster"

plotdata <- prop.table(table(temp@meta.data[,c(plot.var1, plot.var2)]), margin = 1) %>% as.data.frame()

plotdata[,plot.var2]<- factor(plotdata[,plot.var2], levels = levels.var2)
plotdata$label <- paste0(round(plotdata$Freq, 3) * 100, "%")

p <- ggplot(plotdata, aes(fill = .data[[plot.var2]], y = Freq, x = .data[[plot.var1]])) + 
    geom_bar(position = "fill", stat="identity", alpha = 0.8) + 
    theme_classic() + 
    geom_text(aes(label=label), position = position_stack(vjust = 0.5)) + 
    scale_fill_manual(values = colors) + 
    theme(axis.text = element_text(size = 14, face = "bold", color = "black"), 
          axis.title = element_text(size = 14, face = "bold", color = "black") ) & 
    ggtitle(titlename) & labs(x = "", y = "Percentage")
print(p)

