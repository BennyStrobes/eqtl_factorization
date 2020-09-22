#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=10GB


ase_input_data_dir="$1"
het_prob_file="$2"
annotated_samples_file="$3"
pca_file="$4"
output_root="$5"


python preprocess_ase_data.py $ase_input_data_dir $het_prob_file $annotated_samples_file $pca_file $output_root