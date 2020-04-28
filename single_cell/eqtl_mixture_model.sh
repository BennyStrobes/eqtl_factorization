#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=50GB
#SBATCH --nodes=1



genotype_training_file="$1"
expression_training_file="$2"
sample_overlap_file="$3"
k="$4"
output_root="$5"

date
Rscript eqtl_mixture_model.R $genotype_training_file $expression_training_file $sample_overlap_file $k $output_root
date