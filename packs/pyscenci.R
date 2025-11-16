library(SCopeLoomR)
library(dplyr)


meta <- readRDS("./meta.neuronal.rds")
mat <- readRDS("./mat.nor.neuronal.rds") %>% as.matrix()

meta <- meta[-which(is.na(meta$anno)),]
mat <- mat[,rownames(meta)]

dir.create("./SCENIC")
name = "neuronal"
loom <- SCopeLoomR::build_loom(
  file.name = paste0("./SCENIC/",name ,".loom"),
  dgem = mat,
  default.embedding = NULL  
)


############################################# in the terminal
qsub -I -q medium_ext -l nodes=1:ppn=12,walltime=2:00:00,mem=50gb
ml miniconda3
source activate
conda activate /home/lfliuhku/conda/pkgs/pyscenic

path="/home/lfliuhku/projects/HDAC1/2_human/2_neuronal_lineage/HDAC_target"
cd ${path}

ref_tf="/home/lfliuhku/reference/TF/allTFs_human.txt"
# ref_tf="/home/lfliuhku/reference/TF/allTFs_mouse.txt"


ref_feather="/home/lfliuhku/reference/CisTarget/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather" 
#ref_feather="/home/lfliuhku/projects/Vcl_mouse/Merged/0_processing/branchA/1_Alignment/cisTarget_databases/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather" 

ref_tbl="/home/lfliuhku/reference/CisTarget/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

name="neuronal"

pyscenic grn \
--seed 2015 \
--num_workers 3 \
--output ./SCENIC/${name}_adj.tsv \
--method grnboost2 \
./SCENIC/${name}_copy1.loom \
${ref_tf}


pyscenic ctx \
./SCENIC/${name}_adj.tsv \
${ref_feather} \
--annotations_fname ${ref_tbl} \
--expression_mtx_fname ./SCENIC/${name}_copy1.loom \
--output ./SCENIC/${name}_reg.csv \
--num_workers 3 

pyscenic aucell \
./SCENIC/${name}_copy1.loom \
./SCENIC/${name}_reg.csv \
-o ./${name}_auc_mtx.csv \
--num_workers 3  
  
pyscenic ctx \
./SCENIC/${name}_adj.tsv \
${ref_feather} \
--annotations_fname ${ref_tbl} \
--expression_mtx_fname ./SCENIC/${name}_copy1.loom \
--output ./SCENIC/${name}_reg.gmt \
--num_workers 3
