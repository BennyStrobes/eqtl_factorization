#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1

gene_annotation_file="$1"
processed_expression_dir="$2"
eqtl_input_dir="$3"
genotype_data_dir="$4"


python prepare_eqtl_input.py $gene_annotation_file $processed_expression_dir $eqtl_input_dir $genotype_data_dir