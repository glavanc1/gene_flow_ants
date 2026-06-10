#!/bin/bash

module load gcc
module load vcftools

for DP in {18..20}
do
for MISS in $(seq -f "%f" 0.05 0.05 1)
do
vcftools --vcf populations.snps.vcf --minDP $DP --max-missing $MISS --mac 3 --stdout > vcftools_all.out
done
done

grep 'minDP' vcftools_all.out | cut -f 2 | cut -f 2 -d ' ' > DP
grep 'max-missing' test_filtering.out | cut -f 2 | cut -f 2 -d ' ' > miss
grep 'out of a possible' test_filtering.out | cut -f 4 -d ' ' > kept

paste DP miss kept > summary_filtering.txt

rm DP
rm miss
rm kept
