#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=100GB

#SBATCH --nodes=1


pre_processed_data_dir="$1"
gene_annotation_file="$2"
per_time_step_eqtl_input_data_dir="$3"
per_time_step_eqtl_dir="$4"

if false; then
python prepare_single_cell_eqtl_input_data.py $pre_processed_data_dir $gene_annotation_file $per_time_step_eqtl_input_data_dir
fi



###########################
# Run single-cell eqtl analysis in each cell type
###########################

num_pcs="10"
total_jobs="20"


day_num="3"
echo "single-cell eQTL analysis in day"$day_num" cells with "$num_pcs" PCs"

expression_file=$per_time_step_eqtl_input_data_dir"day_"$day_num"_eqtl_input_expression.txt"
genotype_file=$per_time_step_eqtl_input_data_dir"day_"$day_num"_eqtl_input_genotype.txt"
test_names_file=$per_time_step_eqtl_input_data_dir"day_"$day_num"_eqtl_input_variant_gene_pairs.txt"
sample_overlap_file=$per_time_step_eqtl_input_data_dir"day_"$day_num"_eqtl_input_sample_overlap.txt"
covariate_file=$pre_processed_data_dir"standardized_normalized_per_day_"$day_num"_expression_pca_loadings.txt"
# Output root
output_root=$per_time_step_eqtl_dir"sc_per_time_step_eqtl_analysis_"$day_num"_day_"$num_pcs"_pcs_"
if false; then
for job_number in $(seq 0 `expr $total_jobs - "1"`); do
	sbatch run_eqtl_analysis_with_random_effects_in_parallel.sh $expression_file $genotype_file $test_names_file $covariate_file $sample_overlap_file $num_pcs $output_root $job_number $total_jobs
done
fi


if false; then
for day in $(seq 0 3); do 
	echo $day
	output_root=$per_time_step_eqtl_dir"sc_per_time_step_eqtl_analysis_"$day"_day_"$num_pcs"_pcs_"
	python merge_parallelized_eqtl_calls.py $output_root"all_variant_gene_pairs_" $total_jobs
done
fi

