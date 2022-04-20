#!/bin/sh
###
# reformat count results from featureCount
# The main code is reformat_counts.py
###

### working directory might need be changed
cd /Users/szhang32/Desktop/Mousiplier/data

## remove header and only retain geneid and counts
grep -v '^#' NA_PFC_VTA_1d_featureCounts.txt | cut -f 1,7- | awk '{if($1=="Geneid") $1=""; print $0}' > day1_counts.txt
sed -I '' -e 's/  */\t/g' day1_counts.txt
grep -v '^#' NAc_PFC_VTA_28d_featureCounts.txt | cut -f 1,7- | awk '{if($1=="Geneid") $1=""; print $0}'> day28_counts.txt
sed -I '' -e 's/  */\t/g' day28_counts.txt

python ../scripts/reformat_counts.py day1_counts.txt day28_counts.txt name_map.txt NAc_PFC_VTA_counts.txt

