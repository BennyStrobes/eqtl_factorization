#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1

input_h5py_file="$1"
processed_expression_dir="$2"
visualize_processed_expression_dir="$3"
gene_annotation_file="$4"





if false; then
python preprocess_single_cell_expression.py $input_h5py_file $processed_expression_dir $gene_annotation_file
fi
if false; then
Rscript visualize_processed_single_cell_expression.R $processed_expression_dir $visualize_processed_expression_dir
fi

















