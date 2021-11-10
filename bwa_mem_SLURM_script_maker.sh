#!/bin/bash

## This script makes a separate slurm script for each sample.
## They can be submitted separately through a for loop.
## The script requires three arguments:
dir=$1 ## the directory containing the demultiplexed fastq reads (I usually call it 'clean')
genome=$2 ## the path and name of the reference genome
THREADS=$3 ## the number of threads requested per job

## Usage example:
## bash bwa_mem_SLURM_script_maker.sh $(pwd)/clean/ $(pwd)/genome/1_Tdi_b3v08.fasta 1

cd $dir/..

mkdir bwa_mem_SLURM_scripts
mkdir mapped

for i in $(ls $dir/*.fq.gz)
do

        j=$(echo $i | tr '/' '\n' | grep '.fq.gz' | sed 's/.fq.gz//g')
#       echo "Sample = $j
        echo "#!/bin/bash" > ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH -p cpu" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --account=tschwand_default" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --job-name=bwa_mem_${j}" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --mem=8G" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --output=%x_%j.out" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --error=%x_%j.err" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --time=1:00:00" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --mail-type=NONE" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --nodes=1" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --ntasks=1" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "#SBATCH --cpus-per-task=${THREADS}" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "module load gcc/9.3.0" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "module load bwa/0.7.17" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "module load samtools/1.4" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh
        echo "bwa mem -t ${THREADS} ${genome} $i | samtools sort --output-fmt=BAM --threads ${THREADS} -o $(pwd)/mapped/${j}.bam -" >> ./bwa_mem_SLURM_scripts/${j}_bwa_mem_slurm.sh

done
