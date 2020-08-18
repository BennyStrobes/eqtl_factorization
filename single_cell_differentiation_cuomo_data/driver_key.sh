

################################
# Input data
#################################
# File containing raw, log2(CPM+1) data
normalized_expression_file="/work-zfs/abattle4/lab_data/sc_endo_diff/counts.tsv"

# File containing meta-data for each cell
meta_data_file="/work-zfs/abattle4/lab_data/sc_endo_diff/cell_metadata_cols.tsv"

# File containing genotype data for each individual
genotype_file="/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/genotypes_all.txt"

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






################################
# Run Analysis
################################
sh preprocess_data.sh $normalized_expression_file $meta_data_file $genotype_file $gene_annotation_file $pre_processed_data_dir $visualize_pre_processed_data_dir
