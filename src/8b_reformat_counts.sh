#!/bin/sh
###
# reformat count results from featureCount
# The main code is 8a_reformat_counts.py
###


## remove header and only retain geneid and counts
grep -v '^#' ../data/NA_PFC_VTA_1d_featureCounts.txt | cut -f 1,7- | awk '{if($1=="Geneid") $1=""; print $0}' > ../data/day1_counts.txt
sed -I '' -e 's/  */\t/g' ../data/day1_counts.txt
grep -v '^#' ../data/NAc_PFC_VTA_28d_featureCounts.txt | cut -f 1,7- | awk '{if($1=="Geneid") $1=""; print $0}'> ../data/day28_counts.txt
sed -I '' -e 's/  */\t/g' ../data/day28_counts.txt

python reformat_counts.py ../data/day1_counts.txt ../data/day28_counts.txt ../data/name_map.txt ../data/NAc_PFC_VTA_counts.txt

