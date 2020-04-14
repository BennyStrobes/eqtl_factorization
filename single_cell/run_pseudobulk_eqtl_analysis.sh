#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1

processed_expression_dir="$1"
gene_annotation_file="$2"
genotype_data_dir="$3"
pseudobulk_eqtl_dir="$4"



###########################
# Prepare input data for pseudobulk eqtl analysis
###########################
if false; then
python prepare_pseudobulk_eqtl_analysis_input_data.py $processed_expression_dir $gene_annotation_file $genotype_data_dir $pseudobulk_eqtl_dir
fi


###########################
# Run pseudobulk eqtl analysis in each cell type
###########################
cell_type_file=$pseudobulk_eqtl_dir"cell_types.txt"
while read cell_type; do
	echo "eQTL analysis in "$cell_type
	# Input files
	expression_file=$pseudobulk_eqtl_dir$cell_type"_eqtl_input_expression.txt"
	genotype_file=$pseudobulk_eqtl_dir$cell_type"_eqtl_input_genotype.txt"
	test_names_file=$pseudobulk_eqtl_dir$cell_type"_eqtl_input_variant_gene_pairs.txt"
	covariate_file=$pseudobulk_eqtl_dir$cell_type"_pca_scores.txt"
	# Output root
	output_root=$pseudobulk_eqtl_dir$cell_type"_pseudobulk_eqtl_analysis_"
	# Run pseudobulk eqtl analysis in this cell type
	python run_eqtl_analysis.py $expression_file $genotype_file $test_names_file $covariate_file $output_root
done <$cell_type_file


