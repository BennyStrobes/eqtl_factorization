#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=10GB


ase_input_data_dir="$1"
tissues_file="$2"
gtex_individual_information_file="$3"
cell_type_decomposition_hlv_file="$4"
output_root="$5"


python preprocess_ase_data.py $ase_input_data_dir $tissues_file $gtex_individual_information_file $cell_type_decomposition_hlv_file $output_root