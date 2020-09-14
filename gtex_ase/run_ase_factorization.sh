#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=10GB



ase_file="$1"
k="$2"
eqtl_results_dir="$3"

module load python/2.7-anaconda53

python run_ase_factorization.py $ase_file $k $eqtl_results_dir