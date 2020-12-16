#!/bin/bash

#SBATCH -J maprate
#SBATCH -c 4
#SBATCH --error=maprate.err
#SBATCH --output=maprate.out
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH -p jgreiter

module load biology
module load samtools

touch mapping_rate.txt

for f in $(ls *.original.bam | sed 's/.original.bam//' | sort -u) 
do

    if [ ! -f ${f}.original.sorted.bam ]; then
        samtools sort -o ${f}.original.sorted.bam ${f}.original.bam
    fi

    if [ ! -f ${f}.original.sorted.rmdup.bam ]; then
        samtools rmdup ${f}.original.sorted.bam ${f}.original.sorted.rmdup.bam
    fi

    total_reads=$(samtools view -c ${f}.original.sorted.rmdup.bam)
    primary_mapped_reads=$(samtools view -c -F 260 ${f}.original.sorted.rmdup.bam)
    printf '%s\t%s\t%s\n' $total_reads $primary_mapped_reads $f >> mapping_rate.txt
done



