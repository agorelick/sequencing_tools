#!/bin/bash

#SBATCH -J bowtie2
#SBATCH -c 20
#SBATCH --error=bowtie2.err
#SBATCH --output=bowtie2.out
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G

## Summary:
## This script executes the following workflow on all files for which an indexed BAM does not already exist:
## 1. Run bowtie2 to align the two paired FASTQs, and pipes the result to samtools in order to compress the output to a BAM file, while filtering out all unmapped reads
## 2. Sort the resulting BAM file with samtools sort
## 3. Index the BAM
## 4. Remove temporary unsorted BAM 

## To run:
## 1. copy this script to the directory containing your FASTQ files
## 2. update the "suffix1/2" variables to match the end of the FASTQ filenames.
## 3. run this script by executing the command "sbatch < run_bowtie2.sh"

## load required modules
module load biology
module load bowtie2
module load samtools

## specify the reference genome (needs to have been indexed using bowtie2)
ref=$GROUP_HOME/assemblies/GRCh38/GRCh38

## modify these to suit the filenames for paired-end FASTQ files
suffix1='_R1.fastq.gz'
suffix2='_R2.fastq.gz'

## for each pair of FASTQs for a sample, use the filename prefix based on the first read file, then run bowtie2/samtools for the pair of files.
for f in $(ls *$suffix1 | sed "s/$suffix1//" | sort -u)
do
    unsorted_bam=${f}.unsorted.bam
    sorted_bam=${f}.bam
    index=${f}.bam.bai

    if [ ! -f $index ]; then
        ## here we run bowtie2 with 4 parallel processes, and directly pipe the results to samtools, which using the -F 4- argument filters out reads that are unmapped.
        bowtie2 -p 4 -x $ref -1 ${f}${suffix1} -2 ${f}${suffix2} | samtools view -hbS -F 4 > $unsorted_bam
        samtools sort -o $sorted_bam $unsorted_bam
        samtools index $sorted_bam
        rm $unsorted_bam
    fi
done
