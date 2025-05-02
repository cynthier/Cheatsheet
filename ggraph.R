library(igraph)
library(ggraph)

colors <- c("unique" = "#f4fa58", "gene" = "#1dbcb6", "Vcl" = "red")
for(i in 1:nrow(datalist)){
    name <- paste0(datalist[i, "stage"], "_", datalist[i, "celltype"])
    temp <- read.csv(paste0("./tradseq_output/string/string_", name, ".csv"))
    temp <- temp[temp$'X.node1' %in% genes[[name]] & temp$'node2' %in% genes[[name]],]

    node <- unique(c(temp$'X.node1', temp$'node2'))
    nodes <- data.frame(
        name = node,
        group = ifelse(node %in% unique_genes_set[[name]], "unique", ifelse(node == "Vcl", "Vcl", "gene")))
    
    edges <- data.frame(
        from = c(temp$'X.node1'),
        to = c(temp$'node2'),
        type = c("PPI")
    )
    graph <- graph_from_data_frame(edges, vertices = nodes)


    pdf(paste0(name, "_PPI_network.pdf"),height = datalist[i, "height"], width = datalist[i, "width"])
    
    options(repr.plot.width = 14, repr.plot.height = 12)
    p <- ggraph(graph, layout = "kk") +  
      geom_edge_link(aes(color = type), arrow = arrow(length = unit(1, "mm")), 
                     width = 0.05, color = "#cecece", alpha = 0.5)+
      geom_node_point(aes(color = group), size = 8, alpha = 0.5) + 
      geom_node_text(aes(label = name), 
                     size = 2, 
                     color = "#333333", 
                     check_overlap = F, 
                     nudge_x = 0, 
                     nudge_y = 0) + 
      scale_color_manual(values = colors)  +
      theme_void()  +
      theme(legend.position = "right")
    print(p)
    dev.off()
}

