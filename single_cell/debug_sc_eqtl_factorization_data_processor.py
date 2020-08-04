import numpy as np 
import pickle
import h5py
import eqtl_factorization_als
import sys
import pdb
from sklearn.linear_model import LinearRegression

def compute_rss_of_each_gene(Y, G, U, V, intercept):
	num_genes = Y.shape[1]
	rss = []
	F = V[0,:]
	V_comp = V[1:,:]

	UV = np.dot(U, V_comp)
	# Iterate across genes
	for gene_num in range(num_genes):
		pred_expr = intercept[gene_num] + UV[:, gene_num]*G[:, gene_num] + F[gene_num]*G[:, gene_num]
		gene_rss = sum(np.square(Y[:, gene_num] - pred_expr))
		rss.append(gene_rss)
	return rss


def compute_r_squared_of_each_gene(Y, G, U, V, intercept, null_F,  null_intercept):
	num_genes = Y.shape[1]
	r_squared = []
	F = V[0,:]
	V_comp = V[1:,:]

	UV = np.dot(U, V_comp)
	# Iterate across genes
	for gene_num in range(num_genes):
		pred_expr = intercept[gene_num] + UV[:, gene_num]*G[:, gene_num] + F[gene_num]*G[:, gene_num]
		gene_rss = sum(np.square(Y[:, gene_num] - pred_expr))
		null_pred_expr = null_intercept[gene_num] + null_F[gene_num]*G[:, gene_num]
		gene_tss = sum(np.square(Y[:,gene_num] - null_pred_expr))
		gene_r_squared = 1.0 - (gene_rss/gene_tss)
		r_squared.append(gene_r_squared)
	return r_squared

def extract_vector_of_known_cell_types_for_cells(cell_covariates_file):
	ct = []
	f = open(cell_covariates_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ct.append(data[1])
	return np.asarray(ct)

def fit_eqtl_model(Y, G):
	###
	# First make loading matrix
	num_cells = Y.shape[0]
	num_tests = Y.shape[1]
	intercept = []
	F = []
	# Iterate through tests
	for test_num in range(num_tests):
		print(test_num)
		Y_test = Y[:, test_num]
		G_test = G[:, test_num]
		cov = np.zeros((num_cells, 1))
		cov[:,0] = G_test
		reg = LinearRegression().fit(cov, Y_test)
		intercept.append(reg.intercept_)
		F.append(reg.coef_[0])
	return np.asarray(F), np.asarray(intercept)

def fit_model_where_loadings_are_one_hot_encodings_of_known_cell_types(Y, G, ct):
	###
	# First make loading matrix
	unique_ct = np.unique(ct)
	num_cells = Y.shape[0]
	num_tests = Y.shape[1]
	U = np.zeros((num_cells, len(unique_ct)))
	dicti = {}
	for i, ct_name in enumerate(unique_ct):
		dicti[ct_name] = i
	for cell_num in range(num_cells):
		cells_cell_type = ct[cell_num]
		col_index = dicti[cells_cell_type]
		U[cell_num, col_index] = 1.0
	###
	# Initialize factor matrix, F, and intercept
	intercept = []
	V_comp = np.zeros((len(unique_ct) + 1, num_tests))
	# Iterate through tests
	for test_num in range(num_tests):
		print(test_num)
		Y_test = Y[:, test_num]
		G_test = G[:, test_num]
		cov = np.zeros((num_cells, len(unique_ct)))
		for ct_num in range(len(unique_ct)):
			cov[:, ct_num] = U[:, ct_num]*G_test
		reg = LinearRegression().fit(np.hstack((np.transpose(np.asmatrix(G_test)),cov)), Y_test)
		intercept.append(reg.intercept_)
		V_comp[:, test_num] = reg.coef_
	return U, V_comp, np.asarray(intercept)

def get_gene_stats(Y_raw):
	percent_cells_expresssed = []
	total_reads = []
	num_tests = Y_raw.shape[1]
	num_cells = float(Y_raw.shape[0])

	for test_num in range(num_tests):
		lib_size = sum(Y_raw[:,test_num])
		total_reads.append(lib_size)
		frac_expressed = sum(Y_raw[:,test_num] > 0)/num_cells
		percent_cells_expresssed.append(frac_expressed)
	return np.asarray(percent_cells_expresssed), np.asarray(total_reads)

output_dir = sys.argv[1]

########################################
# Input data
#########################################
# Raw expression file
raw_expression_file = '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_input/single_cell_random_subset_sig_tests_50_pc_raw_expression_training_data_uncorrected_r_squared_pruned.h5' 
Y_raw = np.transpose(np.asarray(h5py.File(raw_expression_file,'r')['data']))
print('loaded')

# Cell covariates file
cell_covariates_file = '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/processed_expression/cell_covariates_sle_individuals_random_subset.txt' 
ct = extract_vector_of_known_cell_types_for_cells(cell_covariates_file)

print('got ct')
# Pre-trained models
file_seed_1 = '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_factorization_results/eqtl_factorization_single_cell_sig_tests_50_pc_data_8_factors_eqtl_factorization_als_model_False_re_False_svi_1_seed_model'
file_seed_2 = '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_factorization_results/eqtl_factorization_single_cell_sig_tests_50_pc_data_8_factors_eqtl_factorization_als_model_False_re_False_svi_2_seed_model'
eqtl_vi_1 = pickle.load(open(file_seed_1, 'rb'))
eqtl_vi_2 = pickle.load(open(file_seed_2, 'rb'))
print('loaded models')


pdb.set_trace()


######################################
# Save correlation between tests
######################################
geno_corr = np.corrcoef(np.transpose(eqtl_vi_1.G))
np.savetxt(output_dir + 'geno_corr.txt', geno_corr, fmt="%s", delimiter='\t')

########################################
# Generate test info (stats on each gene)
#########################################
'''
percent_cells_expresssed, total_reads = get_gene_stats(Y_raw) 
output_file = output_dir + 'percent_cells_expresssed_per_test.txt'
np.savetxt(output_file, percent_cells_expresssed, fmt="%s", delimiter='\t')
output_file = output_dir + 'total_reads_per_test.txt'
np.savetxt(output_file, total_reads, fmt="%s", delimiter='\t')
'''

null_F, null_intercept = fit_eqtl_model(eqtl_vi_1.Y, eqtl_vi_1.G)

########################################
# Fit model where loadings are one hot encodings of known cell types
######################################### 
ct_loading_U, ct_loading_V, ct_loading_intercept = fit_model_where_loadings_are_one_hot_encodings_of_known_cell_types(eqtl_vi_1.Y, eqtl_vi_1.G, ct)



########################################
# Compute correlation matrix amongst U betweeen 2 runs
######################################### 
U_correlation_mat = np.corrcoef(np.transpose(eqtl_vi_1.U), np.transpose(eqtl_vi_2.U))
output_file = output_dir + 'U_correlation_between_runs.txt'
np.savetxt(output_file, U_correlation_mat, fmt="%s", delimiter='\t')

########################################
# Compute correlation matrix amongst V between 2 runs
#########################################
V_correlation_mat = np.corrcoef(eqtl_vi_1.V, eqtl_vi_2.V)
output_file = output_dir + 'V_correlation_between_runs.txt'
np.savetxt(output_file, V_correlation_mat, fmt="%s", delimiter='\t')


########################################
# Compute RSS of each test according to eqtl factorization model 1
#########################################
rss_genes_1 = compute_rss_of_each_gene(eqtl_vi_1.Y, eqtl_vi_1.G, eqtl_vi_1.U, eqtl_vi_1.V, eqtl_vi_1.intercept)
output_file = output_dir + 'rss_run_1.txt'
np.savetxt(output_file, rss_genes_1, fmt="%s", delimiter='\t')


########################################
# Compute RSS of each test according to eqtl factorization model 2
#########################################
rss_genes_2 = compute_rss_of_each_gene(eqtl_vi_2.Y, eqtl_vi_2.G, eqtl_vi_2.U, eqtl_vi_2.V, eqtl_vi_2.intercept)
output_file = output_dir + 'rss_run_2.txt'
np.savetxt(output_file, rss_genes_2, fmt="%s", delimiter='\t')


########################################
# Compute RSS of each test according to cell-type one hot model
#########################################
rss_genes_ct_loading = compute_rss_of_each_gene(eqtl_vi_2.Y, eqtl_vi_2.G, ct_loading_U, ct_loading_V, ct_loading_intercept)
output_file = output_dir + 'rss_run_cell_type_loading.txt'
np.savetxt(output_file, rss_genes_ct_loading, fmt="%s", delimiter='\t')

########################################
# Compute R^2 of each test according to eqtl factorization model 1
#########################################
r_squared_genes_1 = compute_r_squared_of_each_gene(eqtl_vi_1.Y, eqtl_vi_1.G, eqtl_vi_1.U, eqtl_vi_1.V, eqtl_vi_1.intercept, null_F, null_intercept)
output_file = output_dir + 'r_squared_run_1.txt'
np.savetxt(output_file, r_squared_genes_1, fmt="%s", delimiter='\t')


########################################
# Compute R^2 of each test according to eqtl factorization model 2
#########################################
r_squared_genes_2 = compute_r_squared_of_each_gene(eqtl_vi_2.Y, eqtl_vi_2.G, eqtl_vi_2.U, eqtl_vi_2.V, eqtl_vi_2.intercept, null_F, null_intercept)
output_file = output_dir + 'r_squared_run_2.txt'
np.savetxt(output_file, r_squared_genes_2, fmt="%s", delimiter='\t')

########################################
# Compute R^2 of each test according to cell-type one hot model
#########################################
r_squared_genes_ct_loading = compute_r_squared_of_each_gene(eqtl_vi_2.Y, eqtl_vi_2.G, ct_loading_U, ct_loading_V, ct_loading_intercept, null_F, null_intercept)
output_file = output_dir + 'r_squared_run_cell_type_loading.txt'
np.savetxt(output_file, r_squared_genes_ct_loading, fmt="%s", delimiter='\t')



