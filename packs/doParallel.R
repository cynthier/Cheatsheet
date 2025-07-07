library(foreach)
library(doParallel)
registerDoParallel(cores = 5)
set.seed(2025)
packages <- c("AUCell", "Seurat", "dplyr")
path.score <- foreach(i = names(path_genes[[name]]),
                     .combine = "cbind",
                     .packages = packages) %dopar% { 
            
   }

stopImplicitCluster()
saveRDS(path.scores, file = "Top50_pathway_score.rds")
saveRDS(test.scores, file = "Top50_pathway_score_test.rds")
