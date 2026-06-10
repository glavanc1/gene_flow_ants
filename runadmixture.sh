#!/bin/bash

### This script takes 4 arguments:
### $1 is the input file in vcf format, NOT vcf.gz!
### $2 is the minimum number of K that we want to test.
### $3 is the maximum number of K that we want to test.
### $4 is the number of replicates that we want to run.
### Usage example: sbatch admixture.sh cleaned.vcf 1 6 10

### First: "clean" the vcf to change chromosome names to "1"
### Because plink is dumb and doesn't link non-numeric or high chromosome values.

# Isolate the headers
grep '#' $1 > headers.tmp

# remove the first column (chromosome names)
grep -v '#' $1 | cut -f 2- > data.tmp

# define a variable with the length of the data
# that's the number of lines with "1" we want
LENGTH=$(wc -l data.tmp | cut -f 1 -d ' ')

# write one "1" per column for $LENGTH columns
yes 1 | head -n $LENGTH > ones.tmp

paste ones.tmp data.tmp > core.tmp

cat headers.tmp core.tmp > cleaned.vcf


### Then: run plink to convert vcf into bed format

module load gcc
module load plink2

plink2 --vcf cleaned.vcf --make-pgen --out cleaned --allow-extra-chr --sort-vars
plink2 --pfile cleaned --make-bed --out cleaned

### make one folder per replicate
for i in $(seq 1 $4)
do
mkdir rep${i}
cd rep${i}

### For each replicate, run admixture with each k
for k in $( seq $2 $3 )
do
echo "K=${k}"
~/software/admixture -s time --cv=10 ../cleaned.bed $k &> log_K${k}.out
done
