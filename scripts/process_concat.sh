#!/bin/bash

#SBATCH --job-name=concat
#SBATCH --ntasks=1
#SBATCH --mem=200M
#SBATCH --partition=short
#SBATCH --mail-user=<edward.sanders@imm.ox.ac.uk>
#SBATCH --time=00-01:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

module load bcftools
module load htslib

DONOR=$1
bcftools concat work_dir/${DONOR}.chr{1..22}.BAM_counts.vcf.gz work_dir/${DONOR}.chrX.BAM_counts.vcf.gz -Oz -o ${DONOR}.BAM_counts.vcf.gz
tabix -p vcf ${DONOR}.BAM_counts.vcf.gz
