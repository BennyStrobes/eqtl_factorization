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
gtex_tissue_colors_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl/input_data/gtex_colors.txt"


#########################
# Preprocess data
#########################
## 4 tissues case
tissues_file=$input_data_dir"tissues_subset_4.txt"
output_dir=$processed_data_dir"tissues_subset_4_"
if false; then

python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $output_dir

## 10 tissues case
tissues_file=$input_data_dir"tissues_subset_10.txt"
output_dir=$processed_data_dir"tissues_subset_10_"
python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $output_dir

## 20 tissues case
tissues_file=$input_data_dir"tissues_subset_20.txt"
output_dir=$processed_data_dir"tissues_subset_20_"
python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $output_dir
fi

#########################
# Run eqtl factorization model
#########################
tissue_subset_name="tissues_subset_4_"

# eqtl factorization input files (generated in 'simulate_eqtl_factorization_data.py')
sample_overlap_file=$processed_data_dir$tissue_subset_name"individual_id.txt"
# TRAINING
expression_training_file=$processed_data_dir$tissue_subset_name"expr.txt"
genotype_training_file=$processed_data_dir$tissue_subset_name"genotype.txt"
# TESTING
expression_testing_file=$processed_data_dir$tissue_subset_name"expr.txt"
genotype_testing_file=$processed_data_dir$tissue_subset_name"genotype.txt"


################################
# Paramaters
initialization="random"
seed="0"
model_name="alm"  # can either be alm or almm (alternating least model or alternating linear mixed model)
################################
# Run eqtl factorization over a number of parameters
lasso_params=( "0.001" )
num_latent_factor_arr=("4")
################################
# Loop through covariate methods
if false; then
for lasso_param in "${lasso_params[@]}"; do
	for num_latent_factors in "${num_latent_factor_arr[@]}"; do
			lasso_param_v=$lasso_param
			lasso_param_u=$lasso_param
			file_stem="eqtl_factorization_"$tissue_subset_name"gtex_data_"$num_latent_factors"_factors_"$model_name"_lasso_U_"$lasso_param_u"_lasso_V_"$lasso_param_v"_initialization_"$initialization"_"$seed
			sh eqtl_factorization.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $lasso_param_u $lasso_param_v $initialization $seed $model_name
	done
done
fi
model_name="vi_shared_effect"
model_name="vi"
model_name="vi_ard_loadings_only"
model_name="vi_no_ard"
num_latent_factors="4"
file_stem="eqtl_factorization_"$tissue_subset_name"gtex_data_"$num_latent_factors"_factors_"$model_name"_model_"$seed"_seed"
sh eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name





if false; then
python initialization_analysis.py $expression_training_file $genotype_training_file $num_latent_factors $eqtl_results_dir $processed_data_dir"sample_tissue_names.txt"

Rscript visualize_eqtl_factorization.R $processed_data_dir $eqtl_results_dir $visualization_dir $gtex_tissue_colors_file
fi