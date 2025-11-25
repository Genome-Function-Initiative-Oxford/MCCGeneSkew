#!/bin/bash

#SBATCH --job-name=run_wasp
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --partition=short,long
#SBATCH --mail-user=<edward.sanders@imm.ox.ac.uk>
#SBATCH --time=00-23:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err


# module load python-cbrg


source /project/hugheslab/esanders/mambaforge/bin/activate WASP 


ls HUVEC*.tsv > cht_input_files.txt

WASP_PYTHON_LOCATION="/project/Wellcome_Discovery/shared/03_publicly_available_tools/WASP/WASP"

# Don't run with min_as_counts as can filter at end.
# We are combining from multiple experiment approaches so need all 
# sites evaluated
python ${WASP_PYTHON_LOCATION}/CHT/combined_test.py \
	--as_disp /project/Wellcome_Discovery/esanders/01_huvec_skew/06_wasp_files/POINTSeq_forward/cht_as_coef.txt \
	--as_only \
	cht_input_files.txt \
	cht_results.txt
	
# --bnb_disp ${INPUT_DIR_ATAC}/cht_bnb_coef.txt \
