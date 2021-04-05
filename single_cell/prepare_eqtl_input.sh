#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=80GB

#SBATCH --nodes=1

gene_annotation_file="$1"
processed_expression_dir="$2"
eqtl_input_dir="$3"
genotype_data_dir="$4"
pseudobulk_eqtl_dir="$5"
single_cell_eqtl_dir="$6"

module load python/3.7.4-anaconda

python prepare_eqtl_input.py $gene_annotation_file $processed_expression_dir $eqtl_input_dir $genotype_data_dir $pseudobulk_eqtl_dir $single_cell_eqtl_dir