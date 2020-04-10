#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1

processed_expression_dir="$1"
gene_annotation_file="$2"
genotype_data_dir="$3"
pseudobulk_eqtl_dir="$4"





###########################
# Prepare input data for pseudobulk eqtl analysis
###########################
if false; then
python prepare_pseudobulk_eqtl_analysis_input_data.py $processed_expression_dir $gene_annotation_file $genotype_data_dir $pseudobulk_eqtl_dir
fi