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

module load python/2.7-anaconda

if false; then
python generate_go_terms_cell_loadings.py $standardized_total_expression_file $go_terms_file $gencode_gene_annotation_file $processed_data_dir
fi

python preprocess_ase_data.py $ase_input_data_dir $annotated_samples_file $processed_data_dir $gencode_gene_annotation_file $cell_cycle_file


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
ase_binned_by_U_root=$processed_data_dir"ase_"$num_bins"_binned_by_U_"
loading_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/eqtl_results/ase_factorization_via_pymc3_lmm_vb_high_biallelic_fraction_only_subsampled_endoderm_differentiation_3_ase_factorization1_temper_U.txt"
if false; then
python generate_ase_binned_by_U.py $processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only_subsampled.txt" $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled.txt" $num_bins $loading_file $ase_binned_by_U_root
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
