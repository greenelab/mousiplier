#!/bin/sh
#BSUB -J Aging_AD_transformation
#BSUB -oo Aging_AD_transformation.o
#BSUB -eo Aging_AD_transformation.e
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
grep -v '^#' data/Aging_AD_featureCounts.txt | cut -f 1,7- | awk '{if($1=="Geneid") $1=""; print $0}' > data/Aging_AD_featureCounts_no_header.txt
sed -i'' -e 's/  */\t/g' data/Aging_AD_featureCounts_no_header.txt
sed -i'' -e 's/Aligned.sortedByCoord.out.bam//g' data/Aging_AD_featureCounts_no_header.txt

### tranpose counts
transposedCounts="data/transposed_Aging_AD_featureCounts.txt"
python src/3a_count_transpose.py data/Aging_AD_featureCounts_no_header.txt $transposedCounts
rm data/Aging_AD_featureCounts_no_header.txt

### preprocess gene expression based on count data
outfile="data/preprocessed_Aging_AD_counts.txt"
python src/3_preprocess_expression.py $transposedCounts data/gene_lengths.tsv data/plier_pathways.tsv $outfile

### transform the gene expression into latent space
outLV="output/Aging_AD_microglia_LVs.txt"
python src/10_NAc_PFC_VTA_transform.py output/filtered_Z.tsv output/filtered_lambda.txt $outfile $outLV

### reformat LVs
python src/11_reformat_LVs.py $outLV output/reformatted_Aging_AD_LVs.txt


