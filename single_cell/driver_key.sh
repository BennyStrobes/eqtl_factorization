


######################
# Input data
######################
# Data from Ye lab
# Directory containing Raw UMI counts for each cell type
raw_umi_count_dir="/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data_for_ben/raw_umi/"

# Data from Ye lab
# Directory containing meta data for each cell
meta_data_dir="/work-zfs/abattle4/prashanthi/Single_cell_eQTL/data_for_ben/meta_data/per_cell/"





######################
# Output directories
######################
# Root of all output directoreis
output_root="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/"

# Directory containing processed single cell expression
processed_expression_dir=$output_root"processed_expression/"

# Directory containing visualizations of processed single cell expression
visualize_processed_expression_dir=$output_root"visualize_processed_expression/"




sh preprocess_single_cell_expression.sh $raw_umi_count_dir $meta_data_dir $processed_expression_dir $visualize_processed_expression_dir