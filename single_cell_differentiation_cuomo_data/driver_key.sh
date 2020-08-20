

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






################################
# Pre-process data
################################
if false; then
sh preprocess_data.sh $normalized_expression_file $meta_data_file $genotype_dir $gene_annotation_file $pre_processed_data_dir $visualize_pre_processed_data_dir
fi


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
