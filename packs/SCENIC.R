################# Prepare LoomR data
library(dplyr)
library(Seurat)
library(ggplot2)
library(SCopeLoomR)

objs.sub <- readRDS(paste0("../outputdata/obj_normalized_con_ko_final.rds"))
cells <- subset(objs.sub@meta.data,  anno2 %in% type) %>% rownames()
raw.mat <- GetAssayData(objs.sub[,cells], layer = "counts", assay = "RNA")[genes,]

exprMat_filtered <- raw.mat
loom <- SCopeLoomR::build_loom(
    file.name = paste0("./SCENIC/",name ,".loom"),
    dgem = exprMat_filtered,
    default.embedding = NULL  
)
############################## read the results 
split_genes <- function(target){ 
    genes <- regmatches(target, gregexpr("(?<=')[A-Z0-9\\-]+(?:_[A-Z0-9]+)?(?=')", target, perl = TRUE))[[1]]
    values  <- as.numeric(regmatches(target, gregexpr("np\\.float64\\((\\d+\\.\\d+)\\)", target, perl = TRUE))[[1]])
    return(data.frame("genes" = genes, "weight" = values))
}

names.col <- c('TF',
'MotifID',
'AUC',
'NES', ### enrichment score
'MotifSimilarityQvalue',# q < 0.05
'OrthologousIdentity',
'Annotation', # direct or indirect annotation
'Context', # genomic context, promoter or enhancers
'TargetGenes',
'RankAtMax')

reg <- readr::read_csv(paste0("./SCENIC/", name, "_reg.csv"), col_names = FALSE, skip = 1) %>% 
    setNames(names.col) %>%
    select(TF, TargetGenes,NES,MotifID, MotifSimilarityQvalue) %>% na.omit() 
    # filter(MotifSimilarityQvalue < 0.05)

##### read the regulon activity
activity <- read_csv("./SCENIC/neuronal_auc_mtx.csv", ) 



################# pyScenic predict the regulatory relationship

qsub -I -q medium_ext -l nodes=1:ppn=12,walltime=10:00:00,mem=50gb
ml miniconda3
source activate
conda activate /home/lfliuhku/conda/pkgs/pyscenic

path="/home/lfliuhku/projects/Vcl_mouse/Merged/5_neuronBranch/Alignment/"
cd ${path}

ref_tf="/home/lfliuhku/projects/Vcl_mouse/Merged/0_processing/branchA/1_Alignment/SCENIC/allTFs_mm.txt"

ref_feather="/home/lfliuhku/projects//Vcl_mouse/Merged/0_processing/branchA/1_Alignment/cisTarget_databases/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather" 
#ref_feather="/home/lfliuhku/projects/Vcl_mouse/Merged/0_processing/branchA/1_Alignment/cisTarget_databases/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather" 

ref_tbl="/home/lfliuhku/projects/Vcl_mouse/Merged/0_processing/branchA/1_Alignment/SCENIC/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

name = "all"
pyscenic grn \
--seed 2015 \
--num_workers 3 \
--output ./SCENIC/${name}_adj.tsv \
--method grnboost2 \
./SCENIC/${name}-Copy1.loom \
${ref_tf}


pyscenic ctx \
./SCENIC/${name}_adj.tsv \
${ref_feather} \
--annotations_fname ${ref_tbl} \
--expression_mtx_fname ./SCENIC/${name}-Copy1.loom \
--output ./SCENIC/${name}_reg.csv \
--num_workers 3 

pyscenic aucell \
./SCENIC/${name}-Copy1.loom \
./SCENIC/${name}_reg.csv \
-o ./${name}_auc_mtx.csv \
--num_workers 3  
  

pyscenic ctx \
./SCENIC/${name}_adj.tsv \
${ref_feather} \
--annotations_fname ${ref_tbl} \
--expression_mtx_fname ./SCENIC/${name}-Copy1.loom \
--output ./SCENIC/${name}_reg.gmt \
--num_workers 3

######################### heatmap to visualize the activity of the regulon
auc <- read.csv(paste0("./SCENIC/", name, "_auc_mtx.csv"),
                row.names = "Cell", 
                check.names = FALSE) %>% t()
tfs <- c("Klf7","Egr1","E2f1") %>% rev()
auc <- auc[paste0(tfs, "(+)"),] 

cells_con <- pse[cells_con,] %>% arrange(pseudotime) %>% rownames()
cells_ko <- pse[cells_ko,] %>% arrange(pseudotime) %>% rownames()

genes <- rownames(auc)
anno.color <- list(celltype = c("Neuroblast" = "#fdcee6", "BP" = "#6b853e", "BranchA" = "#bca9f5"))
anno.col <- data.frame(celltype = pse[colnames(auc),"anno"], 
                       row.names = colnames(auc))
rownames(anno.col) <- colnames(auc)

auc <- auc[,c(rev(cells_con), cells_ko)]
auc <- apply(auc,MARGIN = 1,function(x){ 
    (x-min(x))/(max(x)-min(x))
}) %>% t() 

auc[is.na(auc)] <- 0
auc[auc > 0.9] <- 0.9

p1 <- pheatmap(
    auc[genes,],
    scale = "none",
    cluster_rows = F, 
    cluster_cols = F, 
    show_colnames = F,
    show_rownames = T, 
    fontsize_row = 8, 
    annotation_colors = anno.color,
    annotation_col = anno.col,
    gaps_col = length(cells_con),
    main = paste0("<-- Control | Vcl cKO -->"), 
);

options(repr.plot.width = 5, repr.plot.height = 2)
print(p1)



######################### binary the regulon
auc <- t(auc)
cells_assig <- AUCell_exploreThresholds(auc,
                                        plotHist = F,
                                        thrP = 0.1, # expected 10%
                                        assignCells = T)

thresholds <- getThresholdSelected(cells_assig)

regulonsCells <- setNames(lapply(names(thresholds), function(x) {
    trh <- thresholds[x]
    names(which(auc[x,] > trh))}),names(thresholds))

regulonActivity <- melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2])) %>% as.matrix()
binaryRegulonActivity[1:4,1:4]

ptions(repr.plot.width = 7, repr.plot.height = 6)
library(pheatmap)
p <- pheatmap(
    binaryRegulonActivity[genes,c(rev(cells_con), cells_ko)],
    scale = "none",
    cluster_rows = T, 
    cluster_cols = F, 
    show_colnames = F,
    show_rownames = T, 
    fontsize_row = 8, 
    # annotation_row = row.anno,
    annotation_colors = anno.color,
    annotation_col = anno.col,
    gaps_col = length(cells_con),
    color = colorRampPalette(rev(c("#d73027", "#a3bfdc")))(10),
    main = paste0("<-- Control | Vcl cKO -->"), 
);


############################## get target genes
gmt_file <- paste0("/SCENIC/", time, "_reg.gmt")
    gene_sets <- GSA.read.gmt(gmt_file)
    tf <- "Ezh2"
    temp <- gene_sets$genesets
    names(temp) <- gene_sets$geneset.names
