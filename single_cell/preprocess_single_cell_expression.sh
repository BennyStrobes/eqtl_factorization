#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1

raw_umi_count_dir="$1"
meta_data_dir="$2"
processed_expression_dir="$3"
visualize_processed_expression_dir="$4"


Rscript preprocess_single_cell_expression.R $raw_umi_count_dir $meta_data_dir $processed_expression_dir $visualize_processed_expression_dir