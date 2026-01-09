#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=120g
#PBS -l walltime=30:00:00
#PBS -m abe
#PBS -N ENCC_control_velocyto
#PBS -q large
#PBS -o /home/lfliuhku/job.out/ENCC_control_velocyto.out
#PBS -e /home/lfliuhku/job.out/ENCC_control_velocyto.err
#PBS -r n

# Load the miniconda3 module
ml miniconda3

# Activate the conda environment that has velocyto installed
source activate velocyto_env

# Move to the working directory
cd /home/lfliuhku/projects/HDAC1/2_human/1_integration/0_Velocity/day9

# ----------------------------------------------------------------------------
# To run velocyto run10x on scRNA-seq data from multiome experiments (joint scATAC+scRNA)
# you need the following:
#   1. Cell Ranger output directory from the RNA data (outs/ from 10x)
#   2. Reference genome annotation GTF file for genes
#   3. Optionally, a repeat mask GTF file (for masking repetitive regions)

# Input requirement:
#   - The RNA (gene expression) data directory must contain the original 10x Genomics output
#     with the following structure:
#       <sample_path>/outs/filtered_feature_bc_matrix/
#       <sample_path>/outs/possorted_genome_bam.bam
#       <sample_path>/outs/possorted_genome_bam.bam.bai
#     If you used cellranger-arc, this structure is available in the "outs" folder.
#   - This is necessary for velocyto to extract spliced/unspliced information.
# ----------------------------------------------------------------------------

# Running velocyto on the RNA count data from your multiome (scATAC+scRNA) experiment:
velocyto run10x \
  -m /home/lfliuhku/reference/rmsk/hg38_rmsk.gtf \
  /home/lfliuhku/rawdata/HDAC1/hNP_day9/ \
  /home/lfliuhku/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf

# NOTE:
# - Velocyto will only analyze the RNA (gene expression) data. It does not use ATAC data.
# - Please make sure the directory "/home/lfliuhku/rawdata/HDAC1/hNP_day9/" contains the full
#   CellRanger (or CellRanger-arc) output, especially "outs/possorted_genome_bam.bam".
# - If you used cellranger-arc, point to the "outs" directory inside your sample folder
#   if needed. For example:
#     velocyto run10x -m ... /path/to/sample/outs /path/to/genes.gtf
# - ATAC data cannot be directly used by velocyto; only the RNA bam + annotations are used.

# Output:
# - Velocyto will create a loom file with spliced/unspliced matrices in your working directory.

