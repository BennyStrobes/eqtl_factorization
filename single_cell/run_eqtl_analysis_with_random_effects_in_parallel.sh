#!/bin/bash -l

#SBATCH
#SBATCH --time=15:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


expression_file="$1"
genotype_file="$2"
test_names_file="$3"
covariate_file="$4"
sample_overlap_file="$5"
num_pcs="$6"
output_root="$7"
job_number="$8"
total_jobs="$9"


source activate sklmer

python run_eqtl_analysis_with_random_effects_in_parallel.py $expression_file $genotype_file $test_names_file $covariate_file $sample_overlap_file $num_pcs $output_root $job_number $total_jobs
