library(tradeSeq)
library(dplyr)
library(monocle)
library(Seurat)
library(SingleCellExperiment)

extract_monocle_info <- function(cds) {
  if (cds@dim_reduce_type != "DDRTree") {
    stop(paste0("For now tradeSeq only support Monocle with DDRTree",
                "reduction. If you want to use another type",
                "please use another format for tradeSeq inputs."))
  }
  # Get the reduced dimension of DDRT
  rd <- t(monocle::reducedDimS(cds)) %>% as.data.frame()
  
  # Get the various lineages info for weights and pseudotime
  y_to_cells <- cds@auxOrderingData[["DDRTree"]]
  y_to_cells <- y_to_cells$pr_graph_cell_proj_closest_vertex %>%
    as.data.frame()
  y_to_cells$cells <- rownames(y_to_cells)
  y_to_cells$Y <- y_to_cells$V1
  root <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root <- y_to_cells$Y[y_to_cells$cells == root]
  mst <- monocle::minSpanningTree(cds)
  endpoints <- names(which(igraph::degree(mst) == 1))
  endpoints <- endpoints[endpoints != paste0("Y_", root)]
  cellWeights <- lapply(endpoints, function(endpoint) {
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    df <- y_to_cells[y_to_cells$Y %in% path, ]
    df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
    colnames(df) <- endpoint
    return(df)
  }) %>% do.call(what = 'cbind', args = .)
  pseudotime <- sapply(cellWeights, function(w) cds$Pseudotime)
  rownames(cellWeights) <- rownames(pseudotime) <- colnames(cds)
  # Get the lineages representation
  edges_rd <- t(monocle::reducedDimK(cds)) %>% as.data.frame()
  rd_lineages <- lapply(endpoints, function(endpoint){
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    path <- paste("Y", path, sep = "_")
    return(edges_rd[path, ])
  })
  return(list("pseudotime" = pseudotime,
              "cellWeights" = as.matrix(cellWeights)))
}
#### 1. select informative genes for analysis 
set.seed(3)

cds <- readRDS("/home/projects/HDAC1/2_human/neuronal_lineage/monocle2/cds.neuronal.rds")
# expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
# diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~newcluster")     
# saveRDS(diff_test_res, file = "diff_test_res.rds")

diff_test_res <- readRDS("/home/projects/HDAC1/2_human/neuronal_lineage/monocle2//diff_test_res.rds")
ordering_genes <- rownames(diff_test_res)[order(diff_test_res$qval)][1:10000]


#### 2. evaluate the K parameters (select the elbow value)
info <- extract_monocle_info(cds)
# icmat <- evaluateK(counts = Biobase::exprs(cds)[ordering_genes, ],
#                    cellWeights = info$cellWeights,
#                    pseudotime = info$pseudotime,
#                    conditions = factor(pData(cds)$group)
#                   )

#### 3. fit the model
sce <- fitGAM(counts = Biobase::exprs(cds)[ordering_genes, ],
              cellWeights = info$cellWeights,
              pseudotime = info$pseudotime, 
              conditions = factor(pData(cds)$group), 
              nknots = 5, parallel = T)
table((rowData(sce)$tradeSeq$converged))

#### 4. identify the DEG genes along pseudotime across group 
assocRes <- associationTest(
    sce, 
    lineages = TRUE, #test for all lineages
    l2fc = log2(2)) # log2 fold change threshold to test against)

rowData(sce)$assocRes <- assocRes

saveRDS(sce, file = "sce_10000.rds")

