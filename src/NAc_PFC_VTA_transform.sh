#!/bin/sh
####
# To transform the gene expression into latent space
####
#BSUB -J NAc_PFC_VTA_transform
#BSUB -oo NAc_PFC_VTA_transform.o
#BSUB -eo NAc_PFC_VTA_transform.e
#BSUB -B
#BSUB -N
#BSUB -M 10000
#BSUB -u Shuo.Zhang@Pennmedicine.upenn.edu
#BSUB -W 80:00
module add python/3.6.3

date

cd /project/eheller_itmat_lab/shuo/Mousiplier/mousiplier/src
outfile="/project/eheller_itmat_lab/shuo/Mousiplier/data/NAc_PFC_VTA_Lvs.txt"
python NAc_PFC_VTA_transform.py ../output/Z.tsv  ../output/lambda.txt /project/eheller_itmat_lab/shuo/Mousiplier/data/preprocessed_NAc_PFC_VTA_counts.txt $outfile

date
