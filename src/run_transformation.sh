#!/bin/sh
#BSUB -J run_transformation
#BSUB -oo run_transformation.o
#BSUB -eo run_transformation.e
#BSUB -B
#BSUB -N
#BSUB -M 10000
#BSUB -u Shuo.Zhang@Pennmedicine.upenn.edu
#BSUB -W 80:00
module add python/3.6.3 
module add R/4.0.2

### set the working directory
workDir="/project/eheller_itmat_lab/shuo/mousiplier"
cd $workDir

### prepare count table from feature count
## remove header and only retain geneid and counts
grep -v '^#' data/NA_PFC_VTA_1d_featureCounts.txt | cut -f 1,7- | awk '{if($1=="Geneid") $1=""; print $0}' > data/day1_counts.txt
sed -i'' -e 's/  */\t/g' data/day1_counts.txt
grep -v '^#' data/NAc_PFC_VTA_28d_featureCounts.txt | cut -f 1,7- | awk '{if($1=="Geneid") $1=""; print $0}'> data/day28_counts.txt
sed -i'' -e 's/  */\t/g' data/day28_counts.txt

python src/8a_reformat_counts.py data/day1_counts.txt data/day28_counts.txt data/name_map.txt data/NAc_PFC_VTA_counts.txt
rm data/day1_counts.txt data/day28_counts.txt

### preprocess gene expression based on count data
outfile="data/preprocessed_NAc_PFC_VTA_counts.txt"
python src/3_preprocess_expression.py data/NAc_PFC_VTA_counts.txt data/gene_lengths.tsv data/plier_pathways.tsv $outfile

### transform the gene expression into latent space
outLV="output/NAc_PFC_VTA_LVs.txt"
python src/10a_NAc_PFC_VTA_transform.py output/Z.tsv  output/lambda.txt data/preprocessed_NAc_PFC_VTA_counts.txt $outLV

### reformat LVs
python src/11a_reformat_LVs.py $outLV output/reformated_NAc_PFC_VTA_LVs.txt


### one-way anova
Rscript src/12_differential_LVs.R
