

################################
################################ GO profile cluster 
library(org.Mm.eg.db)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(GO.db)
library(AnnotationDbi) 

source("https://raw.githubusercontent.com/cynthier/Functions/main/my_functions.R")
result <- compareCluster(
    geneCluster = genes,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    ont = "BP",  
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE) 

go1 <- result@compareClusterResult %>% filter(Cluster == "E13.5_BranchA") %>% arrange(desc(Count))
go2 <- result@compareClusterResult %>% filter(Cluster == "E13.5_BranchA") %>% arrange(desc(Count)) 
terms <- c("mitotic cell cycle phase transition", "chromosome segregation", "regulation of cell cycle phase transition", "nuclear division", "DNA replication")

go <- result@compareClusterResult %>% filter(Description %in% terms)


options(repr.plot.width = 6, repr.plot.height = 4.5) 
p1 <- dotplt.enriched(data = go, 
                      title = "target genes of Ezh2 \n (Top pathways)",  
                      size = "Count", 
                      x = "Cluster") & ylab(NULL)
p1


wb <- createWorkbook()
for(name in unique(all.result$Cluster)){ 
    go_result <- subset(all.result, Cluster == name)
    addWorksheet(wb, sheetName = name)
    temp <- go_result[go_result$Cluster == name,] %>% arrange(desc(Count))
    writeData(wb, sheet =  name, x = temp)
}

saveRDS(all.result, file = "GO_results.rds")
saveWorkbook(wb = wb, file = "GO_results.xlsx", overwrite = T)

################################# get pathway genes accroding to the GO ID
#################################

library(org.Mm.eg.db)
library(openxlsx)
library(dplyr)
library(GO.db)
library(AnnotationDbi)
library(foreach)
library(doParallel)
registerDoParallel(cores = 5)
set.seed(2025)

datalist <- data.frame("stage" = c("E13.5","E15.5", "E13.5", "E15.5"), 
                      "celltype" = c("BP","Neuroblast", "BranchA", "BranchA"))

packages <- c("AnnotationDbi", "clusterProfiler", "openxlsx", "AnnotationDbi", "dplyr", "org.Mm.eg.db") 
gene.tem <- list()
path_genes <- list()

for(name in paste0(datalist$stage, "_", datalist$celltype)){ 
    temp <- all.result %>% filter(Cluster == name)  %>% filter(p.adjust < 0.05) %>% arrange(desc(Count)) %>% head(100) 
    path.genes.temp <- foreach(i = temp$ID, 
                               .packages = packages) %dopar% {
        genes <- AnnotationDbi::select(org.Mm.eg.db,
                                       keys = i,
                                       columns = c("ENTREZID", "SYMBOL"), 
                                       keytype = "GOALL") %>% filter(ONTOLOGYALL == "BP")
        genes$SYMBOL
    } 
    names(path.genes.temp) <- temp$ID 
    path_genes[[name]] <- path.genes.temp
}

stopImplicitCluster()
saveRDS(path_genes, file = "all_Top100_pathway_genes.rds")
