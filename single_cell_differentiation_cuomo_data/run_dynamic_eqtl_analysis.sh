#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


bootstrapped_eqtl_stability_dir="$1"
num_lines="$2"
job_number="$3"
total_jobs="$4"


Rscript run_dynamic_eqtl_analysis.R $bootstrapped_eqtl_stability_dir $num_lines $job_number $total_jobs
