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
root_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl/"
# Directory containing input data
input_data_dir=$root_dir"input_data/"
# Directory containing simulated data
processed_data_dir=$root_dir"processed_data/"
# Directory containing simulated data
processed_20_tissue_data_dir=$root_dir"processed_data_20_tissues/"
# Directory containing eqtl results on simulated data
eqtl_results_dir=$root_dir"eqtl_results/"
# Directory containing visualizations
visualization_dir=$root_dir"visualize/"


#########################
# Input Data
#########################
gtex_expression_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/"
gtex_tpm_dir="/work-zfs/abattle4/lab_data/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm/"
gtex_covariate_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/"
gtex_genotype_dir="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/"
gtex_egene_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/"

tissues_file=$input_data_dir"tissues.txt"
if false; then
python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $processed_20_tissue_data_dir

fi
gtex_expression_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/"
gtex_tpm_dir="/work-zfs/abattle4/lab_data/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm/"
gtex_covariate_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/"
gtex_genotype_dir="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/"
gtex_egene_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/"

tissues_file=$input_data_dir"tissues_small_2.txt"
if false; then
python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $processed_data_dir
fi


#########################
# Run eqtl factorization model
#########################
# eqtl factorization input files (generated in 'simulate_eqtl_factorization_data.py')
sample_overlap_file=$processed_data_dir"individual_id.txt"
# TRAINING
expression_training_file=$processed_data_dir"expr.txt"
genotype_training_file=$processed_data_dir"genotype.txt"
# TESTING
expression_testing_file=$processed_data_dir"expr.txt"
genotype_testing_file=$processed_data_dir"genotype.txt"

num_latent_factors="3"
file_stem="eqtl_factorization_on_4_tissue_gtex_data_"$num_latent_factors"_factors"



################################
# Run eqtl factorization over a number of parameters
#lasso_param_us=( "0.0001" "0.001" "0.01" "0.1" "1")
#asso_param_vs=( "0.0" "0.0001" "0.001" "0.01" "0.1" "1")
initializations=("random")
lasso_param_us=( "0.0001" "0.001" "0.01" "0.1" "1" "10" "100")
lasso_param_vs=( "0.0001" "0.001" "0.01" "0.1" "1")
################################
# Loop through covariate methods
for lasso_param_u in "${lasso_param_us[@]}"; do
	# for lasso_param_v in "${lasso_param_vs[@]}"; do
		for initialization in "${initializations[@]}"; do
			sh eqtl_factorization.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $lasso_param_u $lasso_param_u $initialization 
		done
	# done
done



if false; then
python initialization_analysis.py $expression_training_file $genotype_training_file $num_latent_factors $eqtl_results_dir $processed_data_dir"sample_tissue_names.txt"
fi
if false; then
Rscript visualize_eqtl_factorization.R $processed_data_dir $eqtl_results_dir $visualization_dir
fi
