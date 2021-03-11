#!/bin/bash

#SBATCH -J LIC4L
#SBATCH -p jgreiter
#SBATCH --error=LIC4L.err
#SBATCH --output=LIC4L.log
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G

module load biology
module load samtools
module load bowtie2
module load fastqc
module load py-cutadapt

## run fastqc and cutadapt to clean FASTQ files before Bismark
trim_galore --paired --trim1 LIC4L/LIC4L_R1.fastq.gz LIC4L/LIC4L_R2.fastq.gz --cores 4 --basename LIC4L --fastqc

## run Bismark on the result
bismark --parallel 4 --genome $GROUP_HOME/assemblies/hg19/ -1 LIC4L/LIC4L_val_1.fq.gz -2 LIC4L/LIC4L_val_2.fq.gz

## extract methylation from the result
bismark_methylation_extractor LIC4L/LIC4L_val_1_bismark_bt2_pe.bam --paired-end --comprehensive --merge_non_CpG --parallel 4 --no_overlap --gzip --bedGraph -o LIC4L/

