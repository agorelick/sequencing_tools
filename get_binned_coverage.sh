#!/bin/bash

#SBATCH -J depth
#SBATCH --error=depth.err
#SBATCH --output=depth.out
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH -p jgreiter

# This script invokes a software called featureCounts to quickly count the number of PE reads in a BAM falling in predefined regions (from a GTF file)
# In the GTF file, I split hg19 into 100KB bins.
#
# install Subread to get the 'featureCounts' software, and make sure it is in your $PATH (or ammend the line below to its absolute path)
# http://bioinf.wehi.edu.au/featureCounts/
# PMID: 24227677 

module load biology
module load samtools

for f in $(ls *_PS.bam | sed 's/_PS.bam//' | sort -u)
do
    if [ ! -f ${f}_counts_100kb.txt ]; then
        ## NB: using 'samtools view -q 20' to filter out poorly mapping reads (MAPQ < 20) and then pipe the result to featureCounts
        samtools view -h -q 20 ${f}_PS.bam | featureCounts -p -a hg19_100kb_regions.gtf -F SAF -o ${f}_counts_100kb.txt
    fi
done



