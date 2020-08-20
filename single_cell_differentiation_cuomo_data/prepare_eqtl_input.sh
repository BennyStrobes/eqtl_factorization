#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=100GB

#SBATCH --nodes=1


gene_annotation_file="$1"
pre_processed_data_dir="$2"
eqtl_factorization_input_dir="$3"
per_time_step_eqtl_dir="$4"


python prepare_eqtl_input.py $gene_annotation_file $pre_processed_data_dir $eqtl_factorization_input_dir $per_time_step_eqtl_dir