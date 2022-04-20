#!/bin/sh
###
# preprocess gene expression based on count data
###
#BSUB -J preprocess_expression
#BSUB -oo preprocess_expression.o
#BSUB -eo preprocess_expression.e
#BSUB -B
#BSUB -N
#BSUB -M 30000
#BSUB -u Shuo.Zhang@Pennmedicine.upenn.edu
#BSUB -W 100:00

module add python/3.6.3 

cd /project/eheller_itmat_lab/shuo/Mousiplier/data
outfile="preprocessed_NAc_PFC_VTA_counts.txt"

python /project/eheller_itmat_lab/shuo/Mousiplier/mousiplier/src/3_preprocess_expression.py NAc_PFC_VTA_counts.txt gene_lengths.tsv plier_pathways.tsv $outfile
