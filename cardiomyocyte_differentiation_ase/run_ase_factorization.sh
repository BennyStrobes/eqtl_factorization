#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=10GB



ase_file="$1"
covariate_file="$2"
k="$3"
eqtl_results_dir="$4"

module load python/2.7-anaconda53

python run_ase_factorization.py $ase_file $covariate_file $k $eqtl_results_dir