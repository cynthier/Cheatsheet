library(ape)
library(ggtree)

all_objs.integrated <- subset(all_objs.integrated, group == "HDACdata")
all_objs.integrated <- FindVariableFeatures(all_objs.integrated,nfeatures = 3000) all_objs.integrated <- subset(all_objs.integrated, !seurat_clusters %in% c(11,10)) 

all_objs.integrated$cluster <- paste0("C", all_objs.integrated$seurat_clusters)


ave_exp <- AverageExpression(object = all_objs.integrated, 
                  assays = "RNA",
                 layer = "data",
                  group.by = "cluster")$RNA[VariableFeatures(all_objs.integrated),]
dist_matrix <- dist(t(ave_exp), method = "euclidean")
nj_tree <- nj(dist_matrix)

options(repr.plot.width = 5, repr.plot.height = 7)
ggtree(nj_tree, layout="rectangular"
       ,branch.length = "none"
      ) +  # rectangular
      geom_tiplab(angle = 0, hjust = 0.5, size = 6, vjust = 1) + theme_tree()
