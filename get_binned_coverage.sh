#!/bin/bash

#SBATCH -J bin
#SBATCH -c 4
#SBATCH --error=bin.err
#SBATCH --output=bin.out
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH -p jgreiter

module load biology
module load py-deeptools
cd /home/users/gorelica/box/FanDemirci_ctDNA_exosomes/Inititial_six-samples/bams/orig


for f in $(ls *.original.sorted.rmdup.bam | sed 's/.original.sorted.rmdup.bam//' | sort -u)
do
    bamCoverage --bam ${f}.original.sorted.rmdup.bam -o ${f}.original.sorted.rmdup.bindepth.bigwig --binSize 500000 --outFileFormat bedgraph
done


