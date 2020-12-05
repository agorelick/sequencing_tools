#!/bin/bash

#SBATCH -J bowtie2
#SBATCH -c 32
#SBATCH --error=bowtie2.err
#SBATCH --output=bowtie2.out
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2G

## Summary:
## This script expects paired-end reads which have been run through trimmomatic to remove adapter sequences. 
## Running trimmommatic on paired-end data results in 4 files per fastq-pair: 
## *_1_paired.*, *_1_unpaired.*, *_2_paired.*, *_2_unpaired.*;
## As such, this script should be run the directory containing (at least) the *_1_paired* and *_2_paired* files.
## It obtains the unique file prefixes for each set of trimmed FASTQs (i.e. the prefix corresponding to a single experimental sample)
## It then runs bowtie2 to align the two paired FASTQs, and pipes the result to samtools in order to compress the output to a BAM file.

## load required modules
module load biology
module load bowtie2
module load samtools

## specify the reference genome (needs to be indexed using bowtie2)
ref=$GROUP_HOME/assemblies/GRCh38/GRCh38

for f in $(ls *_1_paired.fq.gz | sed 's/_1_paired.fq.gz//' | sort -u)
do
    bam_file=${f}.bam
    if [ ! -f $bam_file ]; then
        bowtie2 -p 4 -x $ref -1 ${f}_1_paired.fq.gz -2 ${f}_2_paired.fq.gz | samtools view -bS - > $bam_file
    fi
done




