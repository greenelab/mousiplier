#!/bin/sh
####
# To transform the gene expression into latent space
####




outfile="../data/NAc_PFC_VTA_Lvs.txt"
python 10a_NAc_PFC_VTA_transform.py ../output/Z.tsv  ../output/lambda.txt ../data/preprocessed_NAc_PFC_VTA_counts.txt $outfile


