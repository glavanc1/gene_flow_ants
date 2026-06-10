#!/bin/bash

### This script takes 4 arguments:
### $1 is the input file in vcf format.
### $2 is the number of K that we want.
### $3 is the number of replicates that we want to run.
### $4 is a population map, where the first column has individual names and second column population info. Unknown is denoted as "-".

### Usage example: sbatch 2nd_step_admixture.sh cleaned.vcf 6 10 popmap.pop

### First: "clean" the vcf to change chromosome names to "1"
### Because plink is dumb and doesn't like non-numeric or high chromosome values.

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

### Filter the vcf to retain only individuals that are in the popmap ($4)
### and convert the vcf into bed format

module load gcc/10.4.0
module load plink-ng/2.00a3.6

cut -d ' ' -f 1 $4 > keep.tmp

plink2 --vcf cleaned.vcf --keep keep.tmp  --make-pgen --sort-vars --out cleaned
plink2 --pfile cleaned --make-bed --out cleaned

### create a popmap for admixture
cut -f 2 -d ' ' $4 > cleaned.pop


### Run admixture for $3 replicates

for i in $(seq $3)
do

mkdir rep${i}
cd rep${i}

~/software/admixture -s time --cv=10 ../cleaned.bed $2 --supervised &> log_rep${i}.out

cd ../

grep -h CV rep*/log*.out > CV_results.txt

done

### Clean stuff
rm *.tmp
