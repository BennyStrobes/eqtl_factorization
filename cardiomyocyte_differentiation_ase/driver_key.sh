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
root_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/cardiomyocyte_differentiation_ase/"
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
ase_input_data_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/cardiomyocyte_differentiation_ase/input_data/"
het_prob_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/cardiomyocyte_differentiation_ase/input_data/YRI_het_prob_genotype.vcf"
annotated_samples_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/cardiomyocyte_differentiation_ase/input_data/processed_covariates_categorical.txt"
pca_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/cardiomyocyte_differentiation_ase/input_data/principal_components_10.txt"
if false; then
sh preprocess_ase_data.sh $ase_input_data_dir $het_prob_file $annotated_samples_file $pca_file $processed_data_dir
fi

ase_file=$processed_data_dir"filtered_ase_counts_0.35_0.5_min_allele.txt"
covariate_file="NA"
covariate_file=$processed_data_dir"covariates.txt"
sample_overlap_file=$processed_data_dir"sample_overlap.txt"
k="3"
model_name="ase_factorization_via_stan_vb"
model_name="ase_factorization_via_pymc3_vb"
model_name="ase_factorization_via_pymc3_mvn_vb"
model_name="ase_factorization_via_pymc3_lmm_vb"
model_name="ase_factorization"
model_name="ase_factorization_old"

if false; then
random_seed="2"
model_name="ase_factorization_smart_init"
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $k $model_name $eqtl_results_dir$model_name"_cardiomyocyte_differentiation_"$k"_seed_"$random_seed"_" $random_seed

random_seed="3"
model_name="ase_factorization_smart_init"
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $k $model_name $eqtl_results_dir$model_name"_cardiomyocyte_differentiation_"$k"_seed_"$random_seed"_" $random_seed
fi

if false; then
random_seed="2"
model_name="ase_factorization"
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $k $model_name $eqtl_results_dir$model_name"_cardiomyocyte_differentiation_"$k"_seed_"$random_seed"_" $random_seed

random_seed="2"
model_name="ase_factorization_old"
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $k $model_name $eqtl_results_dir$model_name"_cardiomyocyte_differentiation_"$k"_seed_"$random_seed"_" $random_seed


random_seed="3"
model_name="ase_factorization"
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $k $model_name $eqtl_results_dir$model_name"_cardiomyocyte_differentiation_"$k"_seed_"$random_seed"_" $random_seed

random_seed="3"
model_name="ase_factorization_old"
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $k $model_name $eqtl_results_dir$model_name"_cardiomyocyte_differentiation_"$k"_seed_"$random_seed"_" $random_seed

random_seed="4"
model_name="ase_factorization"
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $k $model_name $eqtl_results_dir$model_name"_cardiomyocyte_differentiation_"$k"_seed_"$random_seed"_" $random_seed

random_seed="4"
model_name="ase_factorization_old"
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $k $model_name $eqtl_results_dir$model_name"_cardiomyocyte_differentiation_"$k"_seed_"$random_seed"_" $random_seed
fi


if false; then
Rscript visualize_ase_factorization.R $eqtl_results_dir $processed_data_dir $processed_data_dir"annotated_samples.txt" $visualization_dir
fi


































































###################################
# OLD ANALYSIS: NO LONGER USED
###################################




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
if false; then
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	file_stem="eqtl_factorization_standardized_genotype_"$tissue_subset_name"gtex_data_"$num_latent_factors"_factors_"$model_name"_model_"$random_effects"_re_"$svi"_svi_"$seed"_seed"
	sbatch eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name $random_effects $svi $parrallel
done
fi
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


