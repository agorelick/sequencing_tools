#!/bin/bash

#SBATCH -J bowtie2
#SBATCH -c 5
#SBATCH --error=bowtie2.err
#SBATCH --output=bowtie2.out
#SBATCH --time=8:00:00
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


## NB, samtools flag:
## not primary alignment (0x100)
## read fails platform/vendor quality checks (0x200)
## read is PCR or optical duplicate (0x400)
## supplementary alignment (0x800)


## load required modules
module load biology
module load bowtie2
module load samtools
touch read_counts.txt


## specify the reference genome (needs to be indexed using bowtie2)
ref=$GROUP_HOME/assemblies/GRCh38/GRCh38
for f in $(ls *_1_paired.fq.gz | sed 's/_1_paired.fq.gz//' | sort -u)
do
    fq1=${f}_1_paired.fq.gz
    fq2=${f}_2_paired.fq.gz
    index=bams/${f}.bam.bai
    bam_original=bams/${f}.original.bam
    bam_original_sorted=bams/${f}.original.sorted.bam
    bam_original_sorted_rmdup=bams/${f}.original.sorted.rmdup.bam

    if [ ! -f $index ]; then
        bowtie2 -p 8 -x $ref -1 $fq1 -2 $fq2 | samtools view -hbS -F 3840- > $bam_original
        reads_orig=$(samtools view -c $bam_original)
        samtools sort -o $bam_original_sorted $bam_original
        samtools rmdup $bam_original_sorted $bam_original_sorted_rmdup
        reads_rmdup=$(samtools view -c $bam_original_sorted_rmdup)
        samtools view -hbS -F 4- $bam_original_sorted_rmdup > ${f}.bam
        reads_mapped=$(samtools view -c ${f}.bam) 
        printf '%s\t%s\t%s\t%s\n' $reads_orig $reads_rmdup $reads_mapped $f >> read_counts.txt
        samtools index ${f}.bam

        rm $bam_original_sorted; rm $bam_original_sorted_rmdup
    fi
done




