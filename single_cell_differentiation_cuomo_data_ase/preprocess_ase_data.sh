#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


ase_input_data_dir="$1"
annotated_samples_file="$2"
processed_data_dir="$3"
gencode_gene_annotation_file="$4"
go_terms_file="$5"
standardized_total_expression_file="$6"
cell_cycle_file="$7"
genotype_file="$8"
cell_state_file="$9"

module load python/2.7-anaconda

if false; then
python generate_go_terms_cell_loadings.py $standardized_total_expression_file $go_terms_file $gencode_gene_annotation_file $processed_data_dir
fi

if false; then
python preprocess_ase_data.py $ase_input_data_dir $annotated_samples_file $processed_data_dir $gencode_gene_annotation_file $cell_cycle_file $cell_state_file
fi
if false; then

num_bins="30"
ase_binned_by_pseudotime_root=$processed_data_dir"ase_"$num_bins"_binned_by_pseudotime_"
python generate_ase_binned_by_pseudotime.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5.txt" $num_bins $ase_binned_by_pseudotime_root



num_bins="30"
distribution="folded_binomial"
ase_binned_by_U_root=$processed_data_dir"ase_"$distribution"_mean_"$num_bins"_binned_by_mixture_U_"
loading_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/eqtl_results/ase_factorization_via_pymc3_lmm_mixture_vb_subsampled_high_biallelic_fraction_only_endoderm_differentiation_3_2_ase_factorization_temper_U.txt"
# For pystan
module load python/2.7-anaconda
module load gcc/6.4.0
python generate_ase_mean_binned_by_U.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only_subsampled.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled.txt" $num_bins $loading_file $ase_binned_by_U_root $distribution


num_bins="30"
distribution="binomial"
ase_binned_by_U_root=$processed_data_dir"ase_"$distribution"_mean_"$num_bins"_binned_by_mixture_U_"
loading_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/eqtl_results/ase_factorization_via_pymc3_lmm_mixture_vb_subsampled_high_biallelic_fraction_only_endoderm_differentiation_3_2_ase_factorization_temper_U.txt"
module load python/2.7-anaconda
module load gcc/6.4.0
python generate_ase_mean_binned_by_U.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only_subsampled.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled.txt" $num_bins $loading_file $ase_binned_by_U_root $distribution

num_bins="30"
distribution="folded_beta_binomial"
ase_binned_by_U_root=$processed_data_dir"ase_"$distribution"_mean_"$num_bins"_binned_by_mixture_U_"
loading_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/eqtl_results/ase_factorization_via_pymc3_lmm_mixture_vb_subsampled_high_biallelic_fraction_only_endoderm_differentiation_3_2_ase_factorization_temper_U.txt"
module load python/2.7-anaconda
module load gcc/6.4.0
python generate_ase_mean_binned_by_U.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only_subsampled.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled.txt" $num_bins $loading_file $ase_binned_by_U_root $distribution

num_bins="30"
distribution="beta_binomial"
ase_binned_by_U_root=$processed_data_dir"ase_"$distribution"_mean_"$num_bins"_binned_by_mixture_U_"
loading_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/eqtl_results/ase_factorization_via_pymc3_lmm_mixture_vb_subsampled_high_biallelic_fraction_only_endoderm_differentiation_3_2_ase_factorization_temper_U.txt"
module load python/2.7-anaconda
module load gcc/6.4.0
python generate_ase_mean_binned_by_U.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only_subsampled.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled.txt" $num_bins $loading_file $ase_binned_by_U_root $distribution
fi
module load R/3.5.1
Rscript visualize_processed_allelic_counts.R $processed_data_dir







if false; then
python check_for_sample_contamination.py $processed_data_dir"filtered_ase_counts_0.3_0.5.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5.txt" $genotype_file $processed_data_dir
fi


num_bins="40"
ase_binned_by_pseudotime_root=$processed_data_dir"ase_"$num_bins"_binned_by_pseudotime_"
if false; then
python generate_ase_binned_by_pseudotime.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5.txt" $num_bins $ase_binned_by_pseudotime_root
fi

num_bins="30"
ase_binned_by_pseudotime_root=$processed_data_dir"ase_"$num_bins"_binned_by_pseudotime_"
if false; then
python generate_ase_binned_by_pseudotime.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5.txt" $num_bins $ase_binned_by_pseudotime_root
fi

num_bins="30"
ase_binned_by_U_root=$processed_data_dir"ase_"$num_bins"_binned_by_mixture_U_"
loading_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/eqtl_results/ase_factorization_via_pymc3_lmm_mixture_vb_subsampled_high_biallelic_fraction_only_endoderm_differentiation_3_2_ase_factorization_temper_U.txt"
if false; then
python generate_ase_binned_by_U.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only_subsampled.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled.txt" $num_bins $loading_file $ase_binned_by_U_root
fi


ase_binned_by_category_root=$processed_data_dir"ase_binned_by_donor_long_id_"
categorical_variable_name="donor_long_id"
if false; then
python generate_non_min_ase_binned_by_categorical_variable.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5.txt" $categorical_variable_name $ase_binned_by_category_root
fi

ase_binned_by_category_root=$processed_data_dir"ase_binned_by_experiment_"
categorical_variable_name="experiment"
if false; then
python generate_ase_binned_by_categorical_variable.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5.txt" $categorical_variable_name $ase_binned_by_category_root
fi

ase_binned_by_category_root=$processed_data_dir"ase_binned_by_plate_id_"
categorical_variable_name="plate_id"
if false; then
python generate_ase_binned_by_categorical_variable.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5.txt" $categorical_variable_name $ase_binned_by_category_root
fi

if false; then
python summarize_ase_counts.py $processed_data_dir
fi

if false; then
Rscript visualize_processed_allelic_counts.R $processed_data_dir
fi
