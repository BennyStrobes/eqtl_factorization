


######################
# Input data
######################
# Input single cell expression data (emailed by Meena Subramaniam on Nov. 18, 2019)
input_h5py_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/input_data/CLUESImmVar_processed.V6.1.h5ad"

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


if false; then
sh preprocess_single_cell_expression.sh $input_h5py_file $processed_expression_dir $visualize_processed_expression_dir $gene_annotation_file
fi

sh prepare_eqtl_input.sh $gene_annotation_file $processed_expression_dir $eqtl_input_dir $genotype_data_dir














































#############################
# Old data (no longer used)
##############################