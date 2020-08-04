#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=50GB
#SBATCH --nodes=1

output_dir="$1"

module load python/3.7-anaconda

python debug_sc_eqtl_factorization_data_processor.py $output_dir


if false; then
Rscript debug_sc_eqtl_factorization_data_visualizer.R $output_dir
fi