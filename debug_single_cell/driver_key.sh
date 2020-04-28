


##########################
# Input Directories
##########################
# Root of all input directoreis
input_root="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/"

# Directory containing processed single cell expression
processed_expression_dir=$input_root"processed_expression/"

# Directory containing pseudobulk eqtl data
pseudobulk_eqtl_dir=$input_root"pseudobulk_eqtl/"

# Directory containing single-cell eqtl data
single_cell_eqtl_dir=$input_root"single_cell_eqtl/"


# Raw genoptype data dir
genotype_data_dir="/work-zfs/abattle4/lab_data/Ye_Lab_Data/interactionQTL_data/genotypes/"

# Raw expression file
raw_expression_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/processed_expression/scanpy_processed_single_cell_data.h5ad"

##########################
# Output Directories
##########################
# Root of all ouptut directoreis
output_root="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_eqtl_debug/"

# Wrangle data a bit
data_dir=$output_root"organized_data/"

# Wrangle data a bit
visualization_dir=$output_root"visualize/"


if false; then
python organize_data_for_debugging.py $processed_expression_dir $pseudobulk_eqtl_dir $single_cell_eqtl_dir $data_dir
fi

if false; then
python error_check_to_make_sure_genotypes_are_correct.py $genotype_data_dir $data_dir
fi

python error_check_to_make_sure_expression_is_correct.py $data_dir $raw_expression_file

################
# 