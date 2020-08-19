#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=100GB

#SBATCH --nodes=1


pre_processed_data_dir="$1"
gene_annotation_file="$2"
per_time_step_eqtl_input_data_dir="$3"

python prepare_single_cell_eqtl_input_data.py $pre_processed_data_dir $gene_annotation_file $per_time_step_eqtl_input_data_dir