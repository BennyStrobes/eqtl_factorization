#!/bin/bash -l

#SBATCH
#SBATCH --time=14:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=10GB
#SBATCH --nodes=1


#########################
# Output directories
#########################
# Root directory for simulation analysis
simulation_root_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/simulation/"
# Directory containing simulated data
simulated_data_dir=$simulation_root_dir"simulated_data/"
# Directory containing eqtl results on simulated data
eqtl_results_dir=$simulation_root_dir"eqtl_results/"
# Directory containing visualizations
visualization_dir=$simulation_root_dir"visualize_simulation/"



#########################
# Model paramaters
#########################
# Number of variant gene pairs in simulation analysis
num_tests="5000"
# Number of individuals in simulation analysis
num_individuals="100"
# Numbers of cells per individual in simulation analysis
num_cells_per_individual="100"
# Number of latent factors in simulation analysis
num_latent_factors="3"
# Minor allele frequency
maf=".3"
# Output file stem (based on model parameters)
file_stem="eqtl_factorization_"$num_tests"_tests_"$maf"_maf_"$num_individuals"_individuals_"$num_cells_per_individual"_cells_per_individual_"$num_latent_factors"_factors_"


#########################
# Generate simulated data
#########################
if false; then
python simulate_eqtl_factorization_data.py $num_tests $maf $num_individuals $num_cells_per_individual $num_latent_factors $simulated_data_dir $file_stem
fi



#########################
# Run eqtl factorization model
#########################
# eqtl factorization input files (generated in 'simulate_eqtl_factorization_data.py')
sample_overlap_file=$simulated_data_dir$file_stem"simulated_individual_id.txt"
# TRAINING
expression_training_file=$simulated_data_dir$file_stem"simulated_expression.txt"
genotype_training_file=$simulated_data_dir$file_stem"simulated_genotype.txt"
# TESTING
expression_testing_file=$simulated_data_dir$file_stem"simulated_expression.txt"
genotype_testing_file=$simulated_data_dir$file_stem"simulated_genotype.txt"
if false; then
python eqtl_factorization.py $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir
fi


Rscript visualize_eqtl_factorization.R $simulated_data_dir $eqtl_results_dir $visualization_dir $file_stem