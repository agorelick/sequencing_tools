#!/bin/bash

#SBATCH -J bin
#SBATCH -c 2
#SBATCH --error=bin.err
#SBATCH --output=bin.out
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH -p jgreiter

# This script invokes a software called mosdepth to calculate depth of coverage for each bam file with the given filename suffix. The options provided are for: 
# (1) running with 4 threads
# (2) only consider reads mapped with MAPQ = 10+ (meaning <= 10% likelihood read is incorrectly mapped)
# (3) Not sure, run faster?
# (4) Don't print coverage for EVERY basepair
# (5) Get coverage in 500kb tiles
#
# Use bioconda to install mosdepth: conda install -c bioconda mosdepth 
# PMID: 29096012

cd /scratch/groups/jgreiter/bams/misc

## Make sure bam files are indexed!
for f in $(ls *.original.sorted.rmdup.bam | sed 's/.original.sorted.rmdup.bam//' | sort -u)
do
    mosdepth -t 4 --mapq 10 --fast-mode -n --by 500000 ${f}.original.sorted.rmdup ${f}.original.sorted.rmdup.bam
done



