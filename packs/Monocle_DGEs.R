library(monocle)
library(Seurat)
library(dplyr)
library(foreach)
library(doParallel)
library(monocle)

# Set up parallel backend
num_cores <- 5  
registerDoParallel(cores = num_cores)

# Parallel loop
foreach(i = samples), .packages = c("monocle", "Seurat", "dplyr")) %dopar% {
  
  metadata <- subset(obj@meta.data, sample == i)
    
  cells <- rownames(metadata)
  matrix <- GetAssayData(obj, layer = "counts")[, cells]
  
  # Gene annotation
  gene_anno <- data.frame("gene_short_name" = rownames(matrix))
  rownames(gene_anno) <- gene_anno$gene_short_name
  
  # Subset cells by group
  cells_ko <- rownames(subset(metadata, group == "cKO"))
  cells_con <- rownames(subset(metadata, group == "Control"))
  
  # Create AnnotatedDataFrame objects
  pd <- new("AnnotatedDataFrame", data = metadata)
  fd <- new("AnnotatedDataFrame", data = gene_anno)
  
  # Create CellDataSet
  cds <- newCellDataSet(
    matrix,
    phenoData = pd,
    featureData = fd,
    expressionFamily = negbinomial.size()
  )
  
  # Estimate size factors and dispersions
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  # Differential gene expression test
  res <- differentialGeneTest(cds, fullModelFormulaStr = "~group")
  res <- res %>% mutate(qval = p.adjust(pval, method = "fdr"))
  
  # Add log2 fold change
  a <- 1 + apply(matrix[rownames(res), cells_ko], 1, mean)
  b <- 1 + apply(matrix[rownames(res), cells_con], 1, mean)
  res$log2FC_mean <- log2(a / b)[res$gene]
  
  # Add percentage expressed
  res$bcg.pct <- (rowSums(matrix[, cells_con] > 0) / length(cells_con))[res$gene]
  res$obj.pct <- (rowSums(matrix[, cells_ko] > 0) / length(cells_ko))[res$gene]
  
  # Save results
    name <- "BranchA"
  saveRDS(res, file = paste0("./tradseq_output/3. monocles_DEGs_", name, ".rds"))
}
stopImplicitCluster()
