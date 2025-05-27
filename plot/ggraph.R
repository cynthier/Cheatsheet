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
    V(graph)$label <- as.character(1:vcount(graph)) 


    pdf(paste0(name, "_PPI_network.pdf"),height = datalist[i, "height"], width = datalist[i, "width"])
    options(repr.plot.width = 14, repr.plot.height = 12)

    ########### plot
    p <- ggraph(graph, layout = "stress") +  
      geom_edge_link(aes(color = type), 
                     arrow = arrow(length = unit(1, "mm")), 
                     # color = "#cecece", 
                     # width = ifelse(E(graph)$type == "pp", 1.2, 0.6),
                     edge_alpha = 0.3
                    )+
        scale_edge_colour_manual(values = colors)  +
      geom_node_point(aes(color = group), 
                      size = 6.5, 
                      alpha = ifelse(V(graph)$name %in% c(names(temp), gene_large), 1, 0.5)) + 
      geom_node_text(aes(label = name),
                     size = ifelse(V(graph)$name %in% c(names(temp), "Vcl"),2, 
                                   ifelse(V(graph)$name %in% gene_large, 2, 1)),
                     check_overlap = F, 
                     face =  ifelse(V(graph)$name %in% c(names(temp), gene_large), "bold", "plain"),
                     nudge_x = 0, 
                     nudge_y = 0,
                     color = ifelse(V(graph)$name %in% gene_large, "#4374fe", "black")
                    ) + 
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
    
      theme_void()  +
      theme(legend.position = "right")
    print(p) 
}

