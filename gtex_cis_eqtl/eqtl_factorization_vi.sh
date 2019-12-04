#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --mem=10GB


sample_overlap_file="$1"
expression_training_file="$2"
genotype_training_file="$3"
expression_testing_file="$4"
genotype_testing_file="$5"
num_latent_factors="$6"
file_stem="$7"
eqtl_results_dir="$8"
seed="$9"
model_name="${10}"
bernoulli_prob="${11}"

python run_eqtl_factorization_vi.py $sample_overlap_file $expression_testing_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name $bernoulli_prob
