#!/bin/sh
###
# preprocess gene expression based on count data
###


outfile="../data/preprocessed_NAc_PFC_VTA_counts.txt"

python 3_preprocess_expression.py ../data/NAc_PFC_VTA_counts.txt ../data/gene_lengths.tsv ../data/plier_pathways.tsv $outfile
