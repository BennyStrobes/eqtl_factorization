#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem

#SBATCH --nodes=1


gene_annotation_file="$1"
pre_processed_data_dir="$2"
eqtl_factorization_input_dir="$3"
per_time_step_eqtl_dir="$4"
bootstrapped_eqtl_stability_dir="$5"

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

######################################
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

python prepare_eqtl_input.py $gene_annotation_file $pre_processed_data_dir $eqtl_factorization_input_dir $bootstrapped_eqtl_stability_dir
