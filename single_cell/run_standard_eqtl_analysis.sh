#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1

processed_expression_dir="$1"
gene_annotation_file="$2"
genotype_data_dir="$3"
pseudobulk_eqtl_dir="$4"
single_cell_eqtl_dir="$5"
visualize_pseudobulk_eqtl_dir="$6"


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
if false; then
while read cell_type; do
	echo "pseudo-bulk eQTL analysis in "$cell_type
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
fi

###########################
# Prepare input data for single-cell eqtl analysis
###########################
if false; then
python prepare_single_cell_eqtl_analysis_input_data.py $processed_expression_dir $gene_annotation_file $genotype_data_dir $single_cell_eqtl_dir
fi



###########################
# Run single-cell eqtl analysis in each cell type
###########################
cell_type="B_cells"
num_pcs="15"
# How many nodes to run in parallel
total_jobs="50"
echo "single-cell eQTL analysis in "$cell_type" with "$num_pcs" PCs"
# Input files
expression_file=$single_cell_eqtl_dir$cell_type"_eqtl_input_expression.txt"
genotype_file=$single_cell_eqtl_dir$cell_type"_eqtl_input_genotype.txt"
test_names_file=$single_cell_eqtl_dir$cell_type"_eqtl_input_variant_gene_pairs.txt"
covariate_file=$processed_expression_dir$cell_type"_pca_scores_sle_individuals.txt"
sample_overlap_file=$single_cell_eqtl_dir$cell_type"_eqtl_input_sample_overlap.txt"
# Output root
output_root=$single_cell_eqtl_dir$cell_type"_sc_eqtl_analysis_"$num_pcs"_pcs_"

if false; then
for job_number in $(seq 0 `expr $total_jobs - "1"`); do
	sbatch run_eqtl_analysis_with_random_effects_in_parallel.sh $expression_file $genotype_file $test_names_file $covariate_file $sample_overlap_file $num_pcs $output_root $job_number $total_jobs
done
fi

###########################
# Visualize pseudobulk eqtl results
###########################
if false; then
Rscript visualize_pseudobulk_eqtls.R $processed_expression_dir $pseudobulk_eqtl_dir $visualize_pseudobulk_eqtl_dir
fi