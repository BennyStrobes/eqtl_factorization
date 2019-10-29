#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --mem=30GB


sample_overlap_file="$1"
expression_training_file="$2"
genotype_training_file="$3"
expression_testing_file="$4"
genotype_testing_file="$5"
num_latent_factors="$6"
file_stem="$7"
eqtl_results_dir="$8"
lasso_param_u="$9"
lasso_param_v="${10}"
initialization="${11}"
seed="${12}"
model_name="${13}"

python eqtl_factorization.py $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $lasso_param_u $lasso_param_v $initialization $seed $model_name

