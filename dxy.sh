#!/bin/bash

## This script takes 3 arguments:

## $1 is the path to the working directory

## $2 is the ABSOLUTE path to a gvcf

## $3 is the path to a population map
## which is a tab-delimited file where the first column is individual ID (the same as in the provided VCF)
## and the second column is a population identifier.

## Exemple:
## do sbatch ~/scripts/pi_estimation.sh /path/to/wd /path/to/populations.all.vcf path/to/popmap.txt
## for some reason using $pwd instead of giving the actual path to the wd as $1 didn't work.

module load gcc
module load samtools
module load vcftools
module load htslib

cd $1

## 1. Format the data

# make a list of individuals to retain
cut -f 1 $3 > keep.tmp

# subset the vcf to retain only the individuals included in the popmap
vcftools --gzvcf $2 --keep keep.tmp --recode --out unsorted

# Sort the gvcf, because tabix needs sorted vcf by position,
# but populations (stacks) outputs variants on the - strand in decreasing order.

# Header lines must not be sorted
grep '^#' unsorted.recode.vcf > clean.vcf

# Sort non-header lines by first, then second columns
grep -v '^#' unsorted.recode.vcf | awk '{ if ($1 !~ /^#/) print $0 }' | sort -k1,1 -k2,2n >> clean.vcf

rm unsorted.recode.vcf

# bgzip (tabix requires bgzipped file)
bgzip clean.vcf

# index the vcf
tabix clean.vcf.gz

## 2. Run pixy
conda activate pixy

# change end of line character (because the file was generated in R on windows)
dos2unix $3

# set window size to 100k because that's what is done in the temp files anyway.
pixy --stats dxy pi fst --vcf clean.vcf.gz --populations $3 --window_size 100000 --n_cores 10 --output_prefix dxy_all
conda deactivate

# Erase temporary files (the vcfs are very heavy!)
rm clean.vcf.gz
rm clean.vcf.gz.tbi
