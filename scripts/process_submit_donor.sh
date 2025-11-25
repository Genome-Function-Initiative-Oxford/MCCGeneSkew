#!/bin/bash

#SBATCH --job-name=bam_counts
#SBATCH --ntasks=15
#SBATCH --mem=200G
#SBATCH --partition=long
#SBATCH --mail-user=<edward.sanders@imm.ox.ac.uk>
#SBATCH --time=5-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

DONOR=$1
echo $DONOR
./process_bams_chromosome.sh $DONOR chrX &
for j in {1..22}; do
    CHROM="chr${j}"
    ./process_bams_chromosome.sh $DONOR $CHROM &
done

wait

./process_concat.sh $DONOR
