
#########################
# Output directories
#########################
# Root directory for simulation analysis
root_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl_ancestry/"
# Directory containing input data
input_data_dir=$root_dir"input_data/"
# Directory containing simulated data
processed_data_dir=$root_dir"processed_data/"
# Directory containing eqtl results on simulated data
eqtl_results_dir=$root_dir"eqtl_results/"
# Directory containing visualizations
visualization_dir=$root_dir"visualize/"


#########################
# Input Data
#########################
gtex_expression_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/"
gtex_tpm_dir="/work-zfs/abattle4/lab_data/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm/"
gtex_covariate_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/"
gtex_genotype_dir="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/"
gtex_egene_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/"
gtex_tissue_colors_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl/input_data/gtex_colors.txt"
gtex_individual_information_file="/work-zfs/abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"



#########################
# Preprocess data
#########################
tissues_file=$input_data_dir"gtex_v8_tissues.txt"
if false; then
python preprocess_gtex_data_for_eqtl_ancestry_factorization.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $processed_data_dir
fi


################################
# Run eqtl factorization over a number of parameters
tissue_name_arr=("Adipose_Subcutaneous")
initializations=("residual_clustering" "random1" )
lasso_param_us=( "0.001"  )
num_latent_factor_arr=("2" "3" "4" "5" "6")
################################
# Loop through covariate methods
if false; then
for tissue_name in "${tissue_name_arr[@]}"; do
	for lasso_param_u in "${lasso_param_us[@]}"; do
		for initialization in "${initializations[@]}"; do
			for num_latent_factors in "${num_latent_factor_arr[@]}"; do
				
				file_stem="eqtl_ancestry_factorization_"$tissue_name_"gtex_data_"$num_latent_factors"_factors"
				lasso_param_v=$lasso_param_u

				sample_overlap_file=$processed_data_dir$tissue_name"_sample_overlap.txt"
				expression_training_file=$processed_data_dir$tissue_name"_expr.txt"
				genotype_training_file=$processed_data_dir$tissue_name"_genotype.txt"
				
				expression_testing_file=$processed_data_dir$tissue_name"_expr.txt"
				genotype_testing_file=$processed_data_dir$tissue_name"_genotype.txt"


				sh eqtl_factorization.sh $sample_overlap_file $expression_training_file $genotype_training_file $expression_testing_file $genotype_testing_file $num_latent_factors $file_stem $eqtl_results_dir $lasso_param_u $lasso_param_v $initialization 
			done
		done
	done
done
fi
