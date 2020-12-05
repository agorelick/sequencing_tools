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
## This script executes the following workflow on all files for which an indexed BAM does not already exist:
## 1. Obtain the unique file prefixes for each set of trimmed FASTQs (i.e. the prefix corresponding to a single experimental sample)
## 2. Run bowtie2 to align the two paired FASTQs, and pipes the result to samtools in order to compress the output to a BAM file, while filtering out all unmapped reads
## 3. Sort the resulting BAM file with samtools sort
## 4. Index the BAM
## 5. Remove temporary unsorted BAM 


## load required modules
module load biology
module load bowtie2
module load samtools

## specify the reference genome (needs to be indexed using bowtie2)
ref=$GROUP_HOME/assemblies/GRCh38/GRCh38

for f in $(ls *_1_paired.fq.gz | sed 's/_1_paired.fq.gz//' | sort -u)
do
    unsorted_bam=${f}.unsorted.bam
    sorted_bam=${f}.bam
    index=${f}.bam.bai

    if [ ! -f $index ]; then
        bowtie2 -p 4 -x $ref -1 ${f}_1_paired.fq.gz -2 ${f}_2_paired.fq.gz | samtools view -hbS -F 4- > $unsorted_bam
        samtools sort -o $sorted_bam $unsorted_bam
        samtools index $sorted_bam
        rm $unsorted_bam
    fi
done




