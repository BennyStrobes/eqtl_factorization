

################################
# Input data
#################################
# File containing raw, log2(CPM+1) data
normalized_expression_file="/work-zfs/abattle4/lab_data/sc_endo_diff/counts.tsv"

# File containing meta-data for each cell
meta_data_file="/work-zfs/abattle4/lab_data/sc_endo_diff/cell_metadata_cols.tsv"

# File containing vcf files for each individual
genotype_dir="/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/"

# Gencode hg19 gene annotation file
gene_annotation_file="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"

# Genotype PC file
genotype_pc_file="/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/genotype_PCs.txt"


################################
# Output Directories
################################
# Root directory
root_directory="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/"
# Directory containing pre-processed results
pre_processed_data_dir=$root_directory"preprocessed_data/"
# Directory containing visualizations of pre-processed results
visualize_pre_processed_data_dir=$root_directory"visualize_preprocessed_data/"
# Directory containing input data for per-time step eQTL analysi
per_time_step_eqtl_input_data_dir=$root_directory"per_time_step_eqtl_input_data/"
# Directory containing output for per-time step eqtl analysis
per_time_step_eqtl_dir=$root_directory"per_time_step_eqtl_results/"
# Directory containing eQTL Factorization input data
eqtl_factorization_input_dir=$root_directory"eqtl_factorization_input/"
# Directory containing eqtl factorization results
eqtl_factorization_results_dir=$root_directory"eqtl_factorization_results/"
# Directory to visualize eqtl factorization results
visualize_eqtl_factorization_results_dir=$root_directory"visualize_eqtl_factorization_results/"




################################
# Pre-process data
################################
sh preprocess_data.sh $normalized_expression_file $meta_data_file $genotype_dir $gene_annotation_file $genotype_pc_file $pre_processed_data_dir $visualize_pre_processed_data_dir



################################
# Run eQTL analysis in each day, seperately
################################
if false; then
sh per_time_step_eqtl_analysis.sh $pre_processed_data_dir $gene_annotation_file $per_time_step_eqtl_input_data_dir $per_time_step_eqtl_dir
fi


################################
# Prepare input for eqtl factorization
################################
if false; then
sh prepare_eqtl_input.sh $gene_annotation_file $pre_processed_data_dir $eqtl_factorization_input_dir $per_time_step_eqtl_dir
fi


######################
# Run eqtl-factorization on single cell data
######################
# eqtl factorization input files (generated in 'prepare_eqtl_input.sh')
sample_overlap_file=$eqtl_factorization_input_dir"single_cell_sig_tests_in_each_day_individual_id.txt"
# TRAINING
expression_training_file=$eqtl_factorization_input_dir"single_cell_sig_tests_in_each_day_expression_training_data_uncorrected_zero_centered_r_squared_pruned.h5"
genotype_training_file=$eqtl_factorization_input_dir"single_cell_sig_tests_in_each_day_standardized_genotype_training_data_uncorrected_r_squared_pruned.h5"
# TESTING
expression_testing_file=$eqtl_factorization_input_dir"single_cell_sig_tests_in_each_day_expression_training_data_uncorrected_zero_centered_r_squared_pruned.h5"
genotype_testing_file=$eqtl_factorization_input_dir"single_cell_sig_tests_in_each_day_standardized_genotype_training_data_uncorrected_r_squared_pruned.h5"

covariate_file=$eqtl_factorization_input_dir"single_cell_sig_tests_in_each_day_covariate_subset_10.txt"


# Paramaters
model_name="eqtl_factorization_vi_factor_loading_spike_and_slab_with_re"
num_latent_factors="5"
random_effects="False"
svi="False"
parrallel="False"
lasso_param_v="1"

seeds=("0")
if false; then
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	file_stem="eqtl_factorization_single_cell_sig_tests_50_pc_min_expressed_cells_$transformation_type_transform_data_corrected_genotype_"$num_latent_factors"_factors_"$model_name"_model_"$random_effects"_re_"$svi"_svi_"$seed"_seed_"$lasso_param_v"_lasso_param"
	sbatch eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_factorization_results_dir $seed $model_name $random_effects $svi $parrallel $lasso_param_v $covariate_file
done
fi

if false; then
Rscript visualize_single_cell_eqtl_factorization.R $pre_processed_data_dir $eqtl_factorization_input_dir $eqtl_factorization_results_dir $visualize_eqtl_factorization_results_dir 
fi
