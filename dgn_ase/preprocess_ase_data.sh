#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=10GB



ase_input_data_dir="$1"
dgn_technical_covariates="$2"
dgn_biological_covariates="$3"
gencode_gene_annotation_file="$4"
output_root="$5"



python preprocess_ase_data.py $ase_input_data_dir $dgn_technical_covariates $dgn_biological_covariates $gencode_gene_annotation_file $output_root