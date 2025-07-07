  metadata.sub <- subset(meta.data, group %in% i)
    matrix <- counts(sce.hdac)[, rownames(metadata.sub)]
    
    gene_anno <- data.frame("gene_short_name" = rownames(matrix))
    rownames(gene_anno) <- gene_anno$gene_short_name
    
    pd <- new("AnnotatedDataFrame", data = metadata.sub)
    fd <- new("AnnotatedDataFrame", data = gene_anno)
     
    cds <- newCellDataSet(
                        matrix, 
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
    
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    
    res <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(pseudotime)")
    res <- res %>% dplyr::mutate(qval = p.adjust(pval, method = "fdr")) 
    res <- dplyr::filter(res, qval < 0.05)
    
    cells <- rownames(metadata.sub)
    
    res$pct <- (rowSums(matrix[,cells] > 0) / length(cells))[res$gene]
    
    saveRDS(res, file = paste0("deg_along_trajectory_", i, ".rds"))
