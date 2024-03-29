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
genotype_pc_file="$5"
pre_processed_data_dir="$6"
visualize_pre_processed_data_dir="$7"
g2m_cell_cycle_genes_file="$8"
s_cell_cycle_genes_file="$9"
cell_modules_file="${10}"


if false; then
python preprocess_genotype_data.py $genotype_dir $pre_processed_data_dir
fi

genotype_file=$pre_processed_data_dir"genotype_mean_inputed.txt"
python preprocess_data.py $normalized_expression_file $meta_data_file $genotype_file $gene_annotation_file $genotype_pc_file $pre_processed_data_dir $g2m_cell_cycle_genes_file $s_cell_cycle_genes_file $cell_modules_file


if false; then
Rscript visualize_processed_data.R $pre_processed_data_dir $visualize_pre_processed_data_dir
fi