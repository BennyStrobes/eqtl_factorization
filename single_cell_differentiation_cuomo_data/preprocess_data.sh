#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=100GB

#SBATCH --nodes=1


normalized_expression_file="$1"
meta_data_file="$2"
genotype_dir="$3"
gene_annotation_file="$4"
pre_processed_data_dir="$5"
visualize_pre_processed_data_dir="$6"

if false; then
python preprocess_genotype_data.py $genotype_dir $pre_processed_data_dir
fi
if false; then
python preprocess_data.py $normalized_expression_file $meta_data_file $gene_annotation_file $pre_processed_data_dir
fi
if false; then
Rscript visualize_processed_data.R $pre_processed_data_dir $visualize_pre_processed_data_dir
fi