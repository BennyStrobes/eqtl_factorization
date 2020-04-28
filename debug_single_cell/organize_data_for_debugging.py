import numpy as np 
import os
import sys
import pdb

def check_make_sure_variant_gene_pairs_are_same(pseudobulk_eqtl_result_file, sc_eqtl_result_file):
	f = open(pseudobulk_eqtl_result_file)
	g = open(sc_eqtl_result_file)

	head_count = 0
	for line in f:
		line = line.rstrip()
		data_f = line.split()
		data_g = g.next().rstrip().split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		test_f = data_f[0] + '_' + data_f[1]
		test_g = data_g[0] + '_' + data_g[1]
		if test_f != test_g:
			print('assumption error!')
	f.close()
	g.close()

def get_test_names_for_debugging_analysis(pseudobulk_eqtl_result_file, sig_pseudobulk_eqtl_results_file):
	f = open(sig_pseudobulk_eqtl_results_file)
	test_names = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		test_names[data[0] + '_' + data[1]] = 1
	f.close()
	test_indices = []
	test_indices_dicti = {}
	ordered_test_names = []
	f = open(pseudobulk_eqtl_result_file)
	head_count = 0
	counter = -1
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		counter = counter + 1
		test_name = data[0] + '_' + data[1]
		if test_name in test_names:
			test_indices.append(counter)
			ordered_test_names.append(test_name)
			test_indices_dicti[counter] = 1
	f.close()
	return np.asarray(ordered_test_names), test_indices_dicti

def filter_to_test_indices(input_file, filtered_output_file, test_indices, header_boolean):
	f = open(input_file)
	t = open(filtered_output_file, 'w')
	if header_boolean == True:
		counter = -1
	else:
		counter = 0
	for line in f:
		line = line.rstrip()
		if counter in test_indices:
			t.write(line + '\n')
		counter = counter + 1
	f.close()
	t.close()

#####################
# Command Line Args
#####################
# Input directories
processed_expression_dir = sys.argv[1]
pseudobulk_eqtl_dir = sys.argv[2]
single_cell_eqtl_dir = sys.argv[3]
# output directory
data_dir = sys.argv[4]

cell_type = 'B_cells'

#######################
# Input files
#######################
pseudobulk_eqtl_result_file = pseudobulk_eqtl_dir + cell_type + '_pseudobulk_eqtl_analysis_all_variant_gene_pairs.txt'
sc_eqtl_result_file = single_cell_eqtl_dir + cell_type + '_sc_eqtl_analysis_15_pcs_all_variant_gene_pairs_merged.txt'

sig_pseudobulk_eqtl_results_file = pseudobulk_eqtl_dir + cell_type + '_pseudobulk_eqtl_analysis_multiple_testing_bf_bh_0.1_fdr_.txt'

pseudobulk_expression_file = pseudobulk_eqtl_dir + cell_type + '_eqtl_input_expression.txt'
sc_expression_file = single_cell_eqtl_dir + cell_type + '_eqtl_input_expression.txt'

pseudobulk_genotype_file = pseudobulk_eqtl_dir + cell_type + '_eqtl_input_genotype.txt'
sc_genotype_file = single_cell_eqtl_dir + cell_type + '_eqtl_input_genotype.txt'

pseudobulk_raw_expression_file = pseudobulk_eqtl_dir + cell_type + '_eqtl_input_raw_expression.txt'
sc_raw_expression_file = single_cell_eqtl_dir + cell_type + '_eqtl_input_raw_expression.txt'

pseudobulk_sample_covariate_file = pseudobulk_eqtl_dir + cell_type + '_sample_covariates.txt'
cell_covariates_file = processed_expression_dir + cell_type + '_cell_covariates_sle_individuals.txt'

cell_sample_overlap_file = single_cell_eqtl_dir + cell_type + '_eqtl_input_sample_overlap.txt'

pseudobulk_covariate_file = pseudobulk_eqtl_dir + cell_type + '_pca_scores.txt'
sc_covariate_file = processed_expression_dir + cell_type + '_pca_scores_sle_individuals.txt'

# First make sure variant-gene pairs from single cell analysis and pseudobulk analysis are the same
check_make_sure_variant_gene_pairs_are_same(pseudobulk_eqtl_result_file, sc_eqtl_result_file)


# Get list of test names to run analysis with
test_names, test_indices = get_test_names_for_debugging_analysis(pseudobulk_eqtl_result_file, sig_pseudobulk_eqtl_results_file)

# Save pseudo-bulk covariates
pseudobulk_sample_info = np.loadtxt(pseudobulk_sample_covariate_file, dtype=str, delimiter='\t')
output_file = data_dir + cell_type + '_pseudobulk_sample_covariates.txt'
np.savetxt(output_file, pseudobulk_sample_info, fmt="%s", delimiter='\t')

# Save cell covariates
cell_covariates = np.loadtxt(cell_covariates_file, dtype=str, delimiter='\t')
output_file = data_dir + cell_type + '_sc_cell_covariates.txt'
np.savetxt(output_file, cell_covariates, fmt="%s", delimiter='\t')

# Save cell sample overlap 
output_file = data_dir + cell_type + '_sc_sample_overlap.txt'
sample_overlap = np.loadtxt(cell_sample_overlap_file)
np.savetxt(output_file, sample_overlap.astype(int), fmt="%s", delimiter='\n')

# Save pseudobulk_eqtl_covariates
output_file = data_dir + cell_type + '_pseudobulk_eqtl_covariates.txt'
covs = np.loadtxt(pseudobulk_covariate_file)
np.savetxt(output_file, covs, fmt="%s", delimiter='\t')

# Save sc_eqtl_covariates
output_file = data_dir + cell_type + '_sc_eqtl_covariates.txt'
covs = np.loadtxt(sc_covariate_file)
np.savetxt(output_file, covs[:,:15], fmt="%s", delimiter='\t')


######
# Filter full files to just test_indices

#1: Filter pseudo-bulk test names file to just these indices
pseudobulk_eqtl_result_debug_file = data_dir + cell_type + '_pseudobulk_eqtl_variant_gene_pairs.txt'
filter_to_test_indices(pseudobulk_eqtl_result_file, pseudobulk_eqtl_result_debug_file, test_indices, True)

#2: Filter sc test names file to just these indices
sc_eqtl_result_debug_file = data_dir + cell_type + '_sc_eqtl_variant_gene_pairs.txt'
filter_to_test_indices(sc_eqtl_result_file, sc_eqtl_result_debug_file, test_indices, True)

#3: Filter pseudo-bulk expression file to just these indices
pseudobulk_expression_debug_file = data_dir + cell_type + '_pseudobulk_expression.txt'
filter_to_test_indices(pseudobulk_expression_file, pseudobulk_expression_debug_file, test_indices, False)

#4: Filter sc expression file to just these indices
sc_expression_debug_file = data_dir + cell_type + '_sc_expression.txt'
filter_to_test_indices(sc_expression_file, sc_expression_debug_file, test_indices, False)

#5: Filter pseudo-bulk genotype file to just these indices
pseudobulk_genotype_debug_file = data_dir + cell_type + '_pseudobulk_genotype.txt'
filter_to_test_indices(pseudobulk_genotype_file, pseudobulk_genotype_debug_file, test_indices, False)

#6: Filter sc genotype file to just these indices
sc_genotype_debug_file = data_dir + cell_type + '_sc_genotype.txt'
filter_to_test_indices(sc_genotype_file, sc_genotype_debug_file, test_indices, False)

#7: Filter raw pseudo-bulk expression file to just these indices
pseudobulk_raw_expression_debug_file = data_dir + cell_type + '_pseudobulk_raw_expression.txt'
filter_to_test_indices(pseudobulk_raw_expression_file, pseudobulk_raw_expression_debug_file, test_indices, False)

#8: Filter raw sc expression file to just these indices
sc_raw_expression_debug_file = data_dir + cell_type + '_sc_raw_expression.txt'
filter_to_test_indices(sc_raw_expression_file, sc_raw_expression_debug_file, test_indices, False)


