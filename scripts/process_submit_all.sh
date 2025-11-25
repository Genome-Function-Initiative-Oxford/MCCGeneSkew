#!/bin/bash

for i in {1..3}; do
    DONOR="HUVEC_D${i}"
    echo $DONOR
    for j in {22..1}; do
        CHROM="chr${j}"
        sbatch process_bams_chromosome.sh $DONOR $CHROM
    done
    sbatch process_bams_chromosome.sh $DONOR chrX
done

