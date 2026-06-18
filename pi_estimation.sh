#!/bin/bash --login

## This script takes 3 arguments:

## $1 is the path to the working directory

## $2 is the ABSOLUTE path to a gvcf

## $3 is the name of a file (without path) with lists of
## individuals of one species for one introgression threshold
## expected name format: Genus_species-maxQlist-threshold.txt

# There must be a directory named maxQ_lists in the working directory
# containing all maxQ lists

## Exemple:
## for i in $(ls path/to/wd/maxQ_lists); do sbatch ~/scripts/pi_estimation.sh /path/to/wd /path/to/populations.all.vcf $i ; done
## for some reason using $pwd instead of giving the actual path to the wd as $1 didn't work.

module load gcc
module load samtools
module load vcftools
module load htslib

cd $1

## 1. Format the data

# get the name of the list of individuals.
FILENAME=$3

# define species name
SP=$(echo $FILENAME | cut -f 1 -d '-')

# define threshold
THRESHOLD=$(echo $FILENAME | cut -f 1,2 -d '.' | cut -f 3 -d '-')

# subset the vcf to retain only the individuals corresponding to the threshold
vcftools --gzvcf $2 --keep $1/maxQ_lists/$FILENAME --recode --out unsorted_${SP}_${THRESHOLD}

# Sort the gvcf, because tabix needs sorted vcf by position,
# but populations (stacks) outputs variants on the - strand in decreasing order.

# Header lines must not be sorted
grep '^#' unsorted_${SP}_${THRESHOLD}.recode.vcf > ${SP}_${THRESHOLD}.vcf

# Sort non-header lines by first, then second columns
grep -v '^#' unsorted_${SP}_${THRESHOLD}.recode.vcf | awk '{ if ($1 !~ /^#/) print $0 }' | sort -k1,1 -k2,2n >> ${SP}_${THRESHOLD}.vcf

rm unsorted_${SP}_${THRESHOLD}.recode.vcf

# bgzip (tabix requires bgzipped file)
bgzip ${SP}_${THRESHOLD}.vcf

# index the vcf
tabix ${SP}_${THRESHOLD}.vcf.gz

## 2 Prepare the populations file.
# We want all individuals to belong to the same population.
# We already activate the conda environment here because dos2unix is also in there
conda activate pixy

# change end of line character (because the files were generated in R on windows)
dos2unix $1/maxQ_lists/$FILENAME

LENGTH=$(cat $1/maxQ_lists/$FILENAME | wc -l)
yes 1 | head -n $LENGTH > ones.tmp
paste $1/maxQ_lists/$FILENAME ones.tmp > popmap_${SP}_${THRESHOLD}.tmp

## 3 Run pixy
# set window size to 100k because that's what is done in the temp files anyway.
pixy --stats pi --vcf ${SP}_${THRESHOLD}.vcf.gz --populations popmap_${SP}_${THRESHOLD}.tmp --window_size 100000 --n_cores 48 --output_prefix ${SP}_${THRESHOLD}
conda deactivate

# Erase temporary files (the vcfs are very heavy!)
rm ${SP}_${THRESHOLD}.vcf.gz
rm ${SP}_${THRESHOLD}.vcf.gz.tbi
rm popmap_${SP}_${THRESHOLD}.tmp
