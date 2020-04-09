


######################
# Input data
######################
# Input single cell expression data (emailed by Meena Subramaniam on Nov. 18, 2019)
input_h5py_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/input_data/CLUESImmVar_nonorm.V6.h5ad"

# Gene annotation file (hg19)
gene_annotation_file="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"

# Directory containing genotype data
genotype_data_dir="/work-zfs/abattle4/lab_data/Ye_Lab_Data/interactionQTL_data/genotypes/"

######################
# Output directories
######################
# Root of all output directoreis
output_root="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/"

# Directory containing processed single cell expression
processed_expression_dir=$output_root"processed_expression/"

# Directory containing visualizations of processed single cell expression
visualize_processed_expression_dir=$output_root"visualize_processed_expression/"

# Directory containing pre-processed eqtl files
eqtl_input_dir=$output_root"eqtl_input/"

# Directory containing pre-processed eqtl files
eqtl_results_dir=$output_root"eqtl_factorization_results/"

# Directory containing pre-processed eqtl files
eqtl_visualization_dir=$output_root"visualize_eqtl_factorization/"



######################
# Preprocess single cell expression
######################
if false; then
sh preprocess_single_cell_expression.sh $input_h5py_file $processed_expression_dir $visualize_processed_expression_dir $gene_annotation_file
fi

######################
# Get single cell expression and genotype data into a format to run eqtl-factorization
######################
if false; then
sh prepare_eqtl_input.sh $gene_annotation_file $processed_expression_dir $eqtl_input_dir $genotype_data_dir
fi

######################
# Run eqtl-factorization
######################
# eqtl factorization input files (generated in 'prepare_eqtl_input.sh')
sample_overlap_file=$eqtl_input_dir"sc_individual_id.txt"
# TRAINING
expression_training_file=$eqtl_input_dir"sc_expression_training_data_corrected_10000_bp_0.5_r_squared_pruned.h5"
genotype_training_file=$eqtl_input_dir"sc_genotype_training_data_corrected_10000_bp_0.5_r_squared_pruned.h5"
# TESTING
expression_testing_file=$eqtl_input_dir"sc_expression_training_data_corrected_10000_bp_0.5_r_squared_pruned.h5"
genotype_testing_file=$eqtl_input_dir"sc_genotype_training_data_corrected_10000_bp_0.5_r_squared_pruned.h5"


# Paramaters
model_name="eqtl_factorization_vi_spike_and_slab_loadings_ard_factors"
num_latent_factors="20"
random_effects="True"
svi="True"
parrallel="True"

seeds=("0" )
for seed in "${seeds[@]}"; do
	echo "Seed: "$seed
	file_stem="eqtl_factorization_sc_datas_"$num_latent_factors"_factors_"$model_name"_model_"$random_effects"_re_"$svi"_svi_"$seed"_seed"
	sh eqtl_factorization_vi.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $seed $model_name $random_effects $svi $parrallel
done


if false; then
Rscript visualize_single_cell_eqtl_factorization.R $processed_expression_dir $eqtl_input_dir $eqtl_results_dir $eqtl_visualization_dir
fi



































#############################
# Old data (no longer used)
##############################