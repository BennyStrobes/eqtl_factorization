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
gtex_eqtl_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/"
gtex_tissue_colors_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl/input_data/gtex_colors.txt"
gtex_individual_information_file="/work-zfs/abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
cell_type_decomposition_hlv_file="/work-zfs/abattle4/marios/GTEx_v8/coloc/scRNA_decon/CIBERSORT_GTEx_LV.csv"
african_ancestry_gtex_eqtl_dir="/work-zfs/abattle4/surya/worksets/for_ben/GTEx_Analysis_v8_AFR_eQTL/"
european_ancestry_gtex_eqtl_dir="/work-zfs/abattle4/surya/worksets/for_ben/GTEx_Analysis_v8_EUR_eQTL/"

#########################
# Preprocess data
#########################
if false; then

## 1 tissues case: HLV
tissues_file=$input_data_dir"tissues_subset_1_Heart_Left_Ventricle.txt"
output_dir=$processed_data_dir"tissues_subset_1_Heart_Left_Ventricle_"
python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_eqtl_dir $cell_type_decomposition_hlv_file $output_dir

## 4 tissues case
tissues_file=$input_data_dir"tissues_subset_4.txt"
output_dir=$processed_data_dir"tissues_subset_4_"
python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_eqtl_dir $cell_type_decomposition_hlv_file $output_dir

## 10 tissues case
tissues_file=$input_data_dir"tissues_subset_10.txt"
output_dir=$processed_data_dir"tissues_subset_10_"
python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_eqtl_dir $cell_type_decomposition_hlv_file $output_dir

## 20 tissues case
tissues_file=$input_data_dir"tissues_subset_20.txt"
output_dir=$processed_data_dir"tissues_subset_20_"
python preprocess_gtex_data_for_eqtl_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_eqtl_dir $cell_type_decomposition_hlv_file $output_dir
fi

#########################
# Get ancestry specific eQTL calls for 4 tissue case
#########################
test_names_file=$processed_data_dir"tissues_subset_4_test_names.txt"
tissues_file=$input_data_dir"tissues_subset_4.txt"
output_root=$processed_data_dir"tissues_subset_4_ancestry_specific_"
if false; then
module load python/3.7.4-anaconda
python extract_ancestry_specific_eqtl_calls.py $test_names_file $tissues_file $european_ancestry_gtex_eqtl_dir $african_ancestry_gtex_eqtl_dir $output_root
fi
#########################
# Run eqtl factorization model for 4 gtex tissues
#########################
tissue_subset_name="tissues_subset_4_"
# eqtl factorization input files (generated in 'simulate_eqtl_factorization_data.py')
sample_overlap_file=$processed_data_dir$tissue_subset_name"individual_id.txt"
# TRAINING
expression_training_file=$processed_data_dir$tissue_subset_name"expr.txt"
genotype_training_file=$processed_data_dir$tissue_subset_name"genotype_standardized.txt"
# TESTING
expression_testing_file=$processed_data_dir$tissue_subset_name"expr.txt"
genotype_testing_file=$processed_data_dir$tissue_subset_name"genotype_standardized.txt"


# Paramaters
model_name="eqtl_factorization_vi_spike_and_slab"
num_latent_factors="30"
random_effects="False"
svi="True"
parrallel="False"

seeds=("0" "1" "2" "3" )
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	file_stem="eqtl_factorization_standardized_genotype_"$tissue_subset_name"gtex_data_"$num_latent_factors"_factors_"$model_name"_model_"$random_effects"_re_"$svi"_svi_"$seed"_seed"
	sbatch eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name $random_effects $svi $parrallel
done

svi="False"
seeds=("0" "1" "2" "3")
if false; then
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	file_stem="eqtl_factorization_standardized_genotype_"$tissue_subset_name"gtex_data_"$num_latent_factors"_factors_"$model_name"_model_"$random_effects"_re_"$svi"_svi_"$seed"_seed"
	sbatch eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name $random_effects $svi $parrallel
done
fi




#########################
# Run eqtl factorization model for 10 gtex tissues
#########################
tissue_subset_name="tissues_subset_10_"
# eqtl factorization input files (generated in 'simulate_eqtl_factorization_data.py')
sample_overlap_file=$processed_data_dir$tissue_subset_name"individual_id.txt"
# TRAINING
expression_training_file=$processed_data_dir$tissue_subset_name"expr.txt"
genotype_training_file=$processed_data_dir$tissue_subset_name"genotype_standardized.txt"
# TESTING
expression_testing_file=$processed_data_dir$tissue_subset_name"expr.txt"
genotype_testing_file=$processed_data_dir$tissue_subset_name"genotype_standardized.txt"


# Paramaters
model_name="eqtl_factorization_vi_spike_and_slab"
num_latent_factors="30"
random_effects="False"
svi="True"
seeds=("0" "1")
if false; then
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	file_stem="eqtl_factorizations_"$tissue_subset_name"gtex_data_"$num_latent_factors"_factors_"$model_name"_model_"$random_effects"_re_"$svi"_svi_"$seed"_seed"
	sbatch eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name $random_effects $svi
done
fi

svi="False"
seeds=("0" "1" )
if false; then
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	file_stem="eqtl_factorizations_"$tissue_subset_name"gtex_data_"$num_latent_factors"_factors_"$model_name"_model_"$random_effects"_re_"$svi"_svi_"$seed"_seed"
	sbatch eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name $random_effects $svi
done
fi


#########################
# Run eqtl factorization model for 20 gtex tissues
#########################
tissue_subset_name="tissues_subset_20_"
# eqtl factorization input files (generated in 'simulate_eqtl_factorization_data.py')
sample_overlap_file=$processed_data_dir$tissue_subset_name"individual_id.txt"
# TRAINING
expression_training_file=$processed_data_dir$tissue_subset_name"expr.txt"
genotype_training_file=$processed_data_dir$tissue_subset_name"genotype_standardized.txt"
# TESTING
expression_testing_file=$processed_data_dir$tissue_subset_name"expr.txt"
genotype_testing_file=$processed_data_dir$tissue_subset_name"genotype_standardized.txt"


# Paramaters
model_name="eqtl_factorization_vi_spike_and_slab"
num_latent_factors="30"
random_effects="False"
svi="True"
parrallel="True"
seeds=("0")
if false; then
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	file_stem="eqtl_factorization_"$tissue_subset_name"gtex_data_"$num_latent_factors"_factors_"$model_name"_model_"$random_effects"_re_"$svi"_svi_"$seed"_seed"
	sh eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name $random_effects $svi $parrallel
done
fi




#########################
# Visualize results
#########################
if false; then
Rscript visualize_eqtl_factorization.R $processed_data_dir $eqtl_results_dir $visualization_dir $gtex_tissue_colors_file
fi














































###################################
# OLD ANALYSIS: NO LONGER USED
###################################
