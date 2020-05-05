#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


expression_file="$1"
genotype_file="$2"
test_names_file="$3"
covariate_file="$4"
sample_overlap_file="$5"
library_size_file="$6"
num_pcs="$7"
output_root="$8"
job_number="$9"
total_jobs="${10}"


source activate sklmer

python run_eqtl_analysis_with_negative_binomial_in_parallel.py $expression_file $genotype_file $test_names_file $covariate_file $sample_overlap_file $library_size_file $num_pcs $output_root $job_number $total_jobs
