#!/bin/bash

#SBATCH --job-name=bam_counts
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --partition=short,long
#SBATCH --mail-user=<edward.sanders@imm.ox.ac.uk>
#SBATCH --time=00-23:59:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err


# module load python-cbrg


source /project/hugheslab/esanders/mambaforge/bin/activate skew_structure 

DONOR=$1
CHROM=$2

VCF="/project/Wellcome_Discovery/shared/01_in_house_data/02_HUVEC/00_WGS/00_PackBio/01_phased_VCFs_no_RefCalls/${DONOR}_phased_variants.vcf.gz"

OUT_DIR=$(pwd)

cd /project/Wellcome_Discovery/esanders/code/SkewStructure/src/

echo $VCF
OUT="${OUT_DIR}/work_dir/${DONOR}.${CHROM}.BAM_counts.vcf.gz"
CONFIG="${OUT_DIR}/${DONOR}.json"

if [ ! -f $OUT.tbi ] ; then
    echo "output file" $OUT
    python new_skew.py --input-vcf $VCF --chrom $CHROM --config $CONFIG --out $OUT

    module load htslib
    tabix -p vcf "${OUT_DIR}/work_dir/${DONOR}.${CHROM}.BAM_counts.vcf.gz"
else
    echo "output file" $OUT.tbi "already exists. Delete file first to rewrite"
fi
