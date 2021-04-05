#!/bin/bash -l

#SBATCH
#SBATCH --time=24:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --cpus-per-task=24


ase_file="$1"
covariate_file="$2"
sample_overlap_file="$3"
batch_overlap_file="$4"
k="$5"
model_name="$6"
eqtl_results_dir="$7"

# For pystan
module load python/2.7-anaconda
module load gcc/6.4.0



python run_ase_factorization.py $ase_file $covariate_file $sample_overlap_file $batch_overlap_file $k $model_name $eqtl_results_dir