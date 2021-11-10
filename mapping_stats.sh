#!/bin/bash
for i in $(ls mapped | grep '.bam')
do
IND=$(echo $i | sed 's/.bam//g')
echo $IND >> mapping_stats.txt
samtools flagstat $i >> mapping.tmp
MAPPED=$( grep "mapped (" mapping.tmp | cut -f 2 -d '(' | cut -f 1 -d '%' )
echo "$IND $MAPPED" >> mapping_summary.txt
cat mapping.tmp >> mapping_stats.txt
rm mapping.tmp
done
