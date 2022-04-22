#!/bin/sh
###
# preprocess gene expression based on count data
###
#BSUB -J preprocess_expression
#BSUB -oo preprocess_expression.o
#BSUB -eo preprocess_expression.e
#BSUB -B
#BSUB -N



outfile="../data/preprocessed_NAc_PFC_VTA_counts.txt"

python 3_preprocess_expression.py ../data/NAc_PFC_VTA_counts.txt ../data/gene_lengths.tsv ../data/plier_pathways.tsv $outfile
