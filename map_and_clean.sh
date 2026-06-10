#!/bin/bash

## $1 is the sample name
## $2 is the path where the trimmed reads are stored
## $3 is the reference genome assembly, which must already have been indexed with bwa index
## $4 is the output path. Defaults to .
## Usage example:
## for i in $(cat individuals.txt); do sbatch ~/scripts/map_and_clean.sh $i /path/to/trimmed/reads /path/to/genome/database /path/to/outputs; done


module load gcc
module load bwa
module load samtools


NAME=$1
OUT_PATH="${4:-$.}"

# Map
bwa mem -t 4 $3 ${2}/${NAME}_1P.fq.gz ${2}/${NAME}_2P.fq.gz | \

# fills in mate coordinates and insert size fields
samtools fixmate -m - - | \

# Sort by coordinates
samtools sort - | \

# Mark and remove duplicates
# 2500 is the distance recommended for NovaSeq by the samtools manual
samtools markdup -r -S -d 2500 - - | \

# edit read group
samtools addreplacerg -r "ID:${1}" -r "SM:${1}" - - | \

# compress to bam and write to file
samtools view -b - > ${OUT_PATH}/${NAME}_clean.bam
