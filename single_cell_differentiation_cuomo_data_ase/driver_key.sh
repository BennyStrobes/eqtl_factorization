
#########################
# Input Data
#########################
ase_input_data_dir="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/input_data/ase_aggregated_by_donor_open_access_lines/"
annotated_samples_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/input_data/cell_covariates.txt"
gencode_gene_annotation_file="/work-zfs/abattle4/lab_data/hg19/gencode_gene_annotations/gencode.v19.annotation.gtf.gz"
go_terms_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/input_data/c5.bp.v5.1.symbols.gmt.txt"
standardized_total_expression_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/preprocessed_data/standardized_normalized_expression_all_cells.txt"
cell_cycle_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/preprocessed_data/cell_cycle_scores.txt"
genotype_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/preprocessed_data/genotype_mean_inputed.txt"
#########################
# Output Data
#########################
output_root="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/"
processed_data_dir=$output_root"processed_data/"
eqtl_results_dir=$output_root"eqtl_results/"
visualization_dir=$output_root"visualize/"
if false; then
sh preprocess_ase_data.sh $ase_input_data_dir $annotated_samples_file $processed_data_dir $gencode_gene_annotation_file $go_terms_file $standardized_total_expression_file $cell_cycle_file $genotype_file
fi

ase_file=$processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only.txt"
ase_file=$processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only_subsampled.txt"

covariate_file="NA"
sample_overlap_file=$processed_data_dir"cell_overlap_0.3_0.5_high_biallelic_fraction_only.txt"
sample_overlap_file=$processed_data_dir"cell_overlap_0.3_0.5_high_biallelic_fraction_only_subsampled.txt"

plate_overlap_file=$processed_data_dir"cell_plate_overlap_0.3_0.5_high_biallelic_fraction_only_subsampled.txt"
experiment_overlap_file=$processed_data_dir"cell_experiment_overlap_0.3_0.5_high_biallelic_fraction_only_subsampled.txt"


k="3"

model_name="ase_factorization_via_pymc3_lmm_mb_vb"
model_name="ase_factorization_via_pymc3_lmm_dirichlet_vb"
model_name="ase_factorization_via_pymc3_lmm_exponential_vb"
model_name="ase_factorization_via_pymc3_lmm_vb"
model_name="ase_factorization_via_pymc3_lmm_exponential_vb"
model_name="ase_factorization_via_pymc3_lmm_horseshoe_vb"
model_name="ase_factorization_via_pca"
model_name="ase_factorization_via_pymc3_lmm_horseshoe_vb"
model_name="ase_factorization_via_pca_regress_out_cell_line"
model_name="ase_factorization_via_pymc3_lmm_vb"
model_name="ase_factorization_via_als"
model_name="ase_factorization_via_pymc3_double_lmm_vb"
model_name="ase_factorization_via_pymc3_binomial_double_lmm_vb"
model_name="ase_factorization_via_pymc3_binomial_double_lmm_vb"
if false; then
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $plate_overlap_file $k $model_name $eqtl_results_dir$model_name"_subsampled_high_biallelic_fraction_only_donor_and_plate_endoderm_differentiation_"$k
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $experiment_overlap_file $k $model_name $eqtl_results_dir$model_name"_subsampled_high_biallelic_fraction_only_donor_and_experiment_endoderm_differentiation_"$k

model_name="ase_factorization_via_pymc3_binomial_lmm_vb"

sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $experiment_overlap_file $k $model_name $eqtl_results_dir$model_name"_subsampled_high_biallelic_fraction_only_donor_endoderm_differentiation_"$k
fi


ase_file=$processed_data_dir"filtered_ase_counts_0.3_0.5_final_subsampled.txt"
covariate_file="NA"
sample_overlap_file=$processed_data_dir"cell_overlap_0.3_0.5_subsampled.txt"
k="3"

plate_overlap_file=$processed_data_dir"cell_plate_overlap_0.3_0.5_subsampled.txt"
experiment_overlap_file=$processed_data_dir"cell_experiment_overlap_0.3_0.5_subsampled.txt"
model_name="ase_factorization_via_pymc3_double_lmm_vb"
if false; then
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $plate_overlap_file $k $model_name $eqtl_results_dir$model_name"_subsampled_donor_and_plate_endoderm_differentiation_"$k

sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $experiment_overlap_file $k $model_name $eqtl_results_dir$model_name"_subsampled_donor_and_experiment_endoderm_differentiation_"$k
fi


ase_file=$processed_data_dir"filtered_ase_counts_0.3_0.5_final_high_biallelic_fraction_only_subsampled_small.txt"
covariate_file="NA"
sample_overlap_file=$processed_data_dir"cell_overlap_0.3_0.5_high_biallelic_fraction_only_subsampled_small.txt"
k="3"

plate_overlap_file=$processed_data_dir"cell_plate_overlap_0.3_0.5_high_biallelic_fraction_only_subsampled_small.txt"
experiment_overlap_file=$processed_data_dir"cell_experiment_overlap_0.3_0.5_high_biallelic_fraction_only_subsampled_small.txt"
model_name="ase_factorization_via_pymc3_double_lmm_vb"
if false; then
sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $plate_overlap_file $k $model_name $eqtl_results_dir$model_name"_subsampled_small_donor_and_plate_endoderm_differentiation_"$k

sbatch run_ase_factorization.sh $ase_file $covariate_file $sample_overlap_file $experiment_overlap_file $k $model_name $eqtl_results_dir$model_name"_subsampled_small_donor_and_experiment_endoderm_differentiation_"$k
fi



Rscript visualize_ase_factorization.R $eqtl_results_dir $processed_data_dir $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled_small.txt" $processed_data_dir"go_terms_cell_loadings_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled_small.txt" $visualization_dir"high_biallelic_fraction_subsampled_small_"

if false; then
Rscript visualize_ase_factorization.R $eqtl_results_dir $processed_data_dir $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled.txt" $processed_data_dir"go_terms_cell_loadings_after_filtering_0.3_0.5_high_biallelic_fraction_only_subsampled.txt" $visualization_dir"high_biallelic_fraction_subsampled_"
fi
if false; then
Rscript visualize_ase_factorization.R $eqtl_results_dir $processed_data_dir $processed_data_dir"cell_info_after_filtering_0.3_0.5_high_biallelic_fraction_only.txt" $processed_data_dir"go_terms_cell_loadings_after_filtering_0.3_0.5_high_biallelic_fraction_only.txt" $visualization_dir"high_biallelic_fraction_"
fi