#!/bin/bash

#SBATCH -J trim
#SBATCH -c 4
#SBATCH --error=trim.err
#SBATCH --output=trim.out
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=2G


## Summary:
## This script expects raw sequencing data with paired-end reads in the form of two FASTQs per sample, with filename format e.g. <PREFIX>_1.fq.gz, <PREFIX>_2.fq.gz.
## It will run the trimmomatic software to remove adapter sequences from the reads, and then for each pair of FASTQs generate 4 files:
## <PREFIX>_1_paired.fq.gz <PREFIX>_1_unpaired.fq.gz <PREFIX>_2_paired.fq.gz <PREFIX>_2_unpaired.fq.gz
## The *_paired.fq.gz files will then be used for the alignment step to create a single BAM file for the sample.
##
## This script should be run the directory containing the <PREFIX>_1.fq.gz, <PREFIX>_2.fq.gz files. It will run in a loop for each pair of samples.
## It will skip previously-processed samples (so delete them if you want to re-run and replace them)
##
## NB: The input argument 'ILLUMINACLIP' specifies that an illumina sequencer was used, and the adapter sequences to remove are specified in the *TruSeq3-PE.fa. This should be determined for each sample by running FASTQC (for ducmentation, see: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
##
## To run:
## 1. Copy this script to the directory containing FASTQs for paired-end sequencing data.
## 2. Modify the suffix1/2 variables to match the end of the FASTQ filenames 
## 3. (Alex to do) Based on the FASTQC report (or if you otherwise know the sequencing platform used), replace the ILLUMINACLIP argument for the adapter that was detected by FASTQC 
## 4. Run command 'sbatch < run_trimmomatic.sh' to execute this script as a slurm job


## specify path to where trimmomatic was downloaded and unpacked
trimmomatic_dir=/home/users/gorelica/install/Trimmomatic-0.39

## modify these to suit the filenames for paired-end FASTQ files
suffix1='_R1.fastq.gz'
suffix2='_R2.fastq.gz'

## for each sample with paired end reads, run trimmommatic. When it completes, continue to next sample, until finished.
for f in $(ls *$suffix1 | sed "s/$suffix1//" | sort -u)
do
    trim1_paired=${f}_1_paired.fq.gz
    trim2_paired=${f}_2_paired.fq.gz
    trim1_unpaired=${f}_1_unpaired.fq.gz
    trim2_unpaired=${f}_2_unpaired.fq.gz

    if [ ! -f $trim1_paired ]; then
        java -jar $trimmomatic_dir/trimmomatic-0.39.jar PE -phred33 ${f}${suffix1} ${f}${suffix2} $trim1_paired $trim1_unpaired $trim2_paired $trim2_unpaired ILLUMINACLIP:$trimmomatic_dir/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
    fi
done


