


raw_gtex_data="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl_ancestry/processed_data/Adipose_Subcutaneous.v8.normalized_expression.bed"

cov_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl_ancestry/processed_data/Adipose_Subcutaneous_sample_covariates.txt"

pc_file="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/Adipose_Subcutaneous.v8.covariates.txt"

gene_subset="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl_ancestry/processed_data/Adipose_Subcutaneous_test_names.txt"




output_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/confounder_debug/"

python process_data.py $raw_gtex_data $cov_file $pc_file $gene_subset $output_dir

Rscript visualize_confounder.R $output_dir
