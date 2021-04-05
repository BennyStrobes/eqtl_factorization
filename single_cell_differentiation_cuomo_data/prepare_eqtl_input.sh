#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=60GB


#SBATCH --nodes=1


gene_annotation_file="$1"
pre_processed_data_dir="$2"
eqtl_factorization_input_dir="$3"
per_time_step_eqtl_dir="$4"
bootstrapped_eqtl_stability_dir="$5"

#####################################
# Dynamic eQTL analysis
######################################
num_lines="466373"
total_jobs="100"
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
	sbatch run_dynamic_eqtl_analysis.sh $bootstrapped_eqtl_stability_dir $num_lines $job_number $total_jobs
done
fi
if false; then
python merge_parallelized_dynamic_eqtls.py $bootstrapped_eqtl_stability_dir"dynamic_eqtl_results_" $total_jobs
fi


module load python/2.7-anaconda
python prepare_eqtl_input.py $gene_annotation_file $pre_processed_data_dir $eqtl_factorization_input_dir $bootstrapped_eqtl_stability_dir


if false; then
module load R/3.5.1
input_root="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/eqtl_factorization_input/"
genotype_training_file=$input_root"single_cell_r_squared_thresh_0.01_standardized_genotype_training_data_corrected_r_squared_pruned.txt"
expression_training_file=$input_root"single_cell_r_squared_thresh_0.01_expression_training_data_corrected_r_squared_pruned.txt"
sample_overlap_file=$input_root"single_cell_r_squared_thresh_0.01_individual_id.txt"
k="4"
output_root=$eqtl_factorization_input_dir"mixture_no_intercept_"$k"_results_"
cell_covariates_file=$pre_processed_data_dir"cell_covariates_appeneded.txt"
Rscript eqtl_mixture_model.R $genotype_training_file $expression_training_file $sample_overlap_file $k $cell_covariates_file $output_root 

fi






























#####################
# OLD
######################
if false; then
python prepare_bootstrapped_eqtl_stability_data.py $gene_annotation_file $pre_processed_data_dir $bootstrapped_eqtl_stability_dir
fi

######################################
# Covariate Modulated eQTL analysis
######################################
num_lines="466373"
total_jobs="100"
job_number="0"
if false; then
sh run_covariate_modulated_eqtl_analysis.sh $bootstrapped_eqtl_stability_dir $num_lines $job_number $total_jobs
fi
if false; then
for job_number in $(seq 0 $(($total_jobs-1))); do 
	sbatch run_covariate_modulated_eqtl_analysis.sh $bootstrapped_eqtl_stability_dir $num_lines $job_number $total_jobs
done
fi

if false; then
python merge_parallelized_covariate_modulated_eqtls.py $bootstrapped_eqtl_stability_dir"covariate_modulated_eqtl_results2_" $total_jobs
fi

#######################