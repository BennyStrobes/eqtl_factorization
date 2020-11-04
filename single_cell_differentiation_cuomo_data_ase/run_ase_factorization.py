import numpy as np 
import os
import sys
import pdb
#import ase_factorization
#import ase_factorization_via_stan_vb
#import ase_factorization_via_pymc3_lmm_vb
import ase_factorization_via_pymc3_lmm_vb
import ase_factorization_via_pymc3_double_lmm_vb
import ase_factorization_via_pymc3_binomial_double_lmm_vb
import ase_factorization_via_pymc3_binomial_lmm_vb

#import ase_factorization_via_pymc3_lmm_mb_vb
#import ase_factorization_via_pymc3_lmm_vb
#import ase_factorization_via_pymc3_lmm_dirichlet_vb
#import ase_factorization_via_pymc3_lmm_exponential_vb
#import ase_factorization_via_pymc3_lmm_horseshoe_vb
#import ase_factorization_via_pca
#import ase_factorization_via_pca_regress_out_cell_line
#import ase_factorization_via_als



def load_in_ase_data(ase_file):
	full_data = np.loadtxt(ase_file, dtype=str, delimiter='\t', comments='*')
	count_data = full_data[1:,1:]
	row_num, col_num = count_data.shape
	allelic_counts = np.zeros((row_num, col_num))
	total_counts = np.zeros((row_num, col_num))
	for ii in range(row_num):
		for jj in range(col_num):
			if count_data[ii,jj] == 'NA':
				allelic_counts[ii,jj] = np.nan
				total_counts[ii,jj] = np.nan
			else:
				ref = int(count_data[ii,jj].split('/')[0])
				tot = int(count_data[ii,jj].split('/')[1])
				ref_min = np.min((ref, tot-ref))
				allelic_counts[ii,jj] = ref_min
				total_counts[ii,jj] = tot
	return np.transpose(allelic_counts), np.transpose(total_counts)

def add_intercept_column_to_matrix(X):
	n,m = X.shape # for generality
	X0 = np.ones((n,1))
	Xnew = np.hstack((X0, X))
	return Xnew

def make_cell_line_vector_into_matrix(Z):
	num_cell_lines = len(np.unique(Z))
	num_cells = len(Z)
	Z_mat = np.zeros((num_cells, num_cell_lines-1))
	for n in range(num_cells):
		line_index = Z[n]
		if line_index != (num_cell_lines-1):
			Z_mat[n, int(line_index)] = 1.0
	return Z_mat

def train_ase_factorization_model(ase_file, covariate_file, sample_overlap_file, batch_overlap_file, k, model_name, output_dir):
	if model_name == 'ase_factorization_via_pymc3_lmm_mb_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz = np.loadtxt(sample_overlap_file)
		ase_factorization_obj = ase_factorization_via_pymc3_lmm_mb_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z=zz)
	elif model_name == 'ase_factorization_via_pymc3_lmm_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz = np.loadtxt(sample_overlap_file)
		ase_factorization_obj = ase_factorization_via_pymc3_lmm_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z=zz)
	elif model_name == 'ase_factorization_via_pymc3_double_lmm_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz_sample = np.loadtxt(sample_overlap_file)
		zz_batch = np.loadtxt(batch_overlap_file)
		ase_factorization_obj = ase_factorization_via_pymc3_double_lmm_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z_sample=zz_sample, z_batch=zz_batch)
	elif model_name == 'ase_factorization_via_pymc3_binomial_double_lmm_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz_sample = np.loadtxt(sample_overlap_file)
		zz_batch = np.loadtxt(batch_overlap_file)
		ase_factorization_obj = ase_factorization_via_pymc3_binomial_double_lmm_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z_sample=zz_sample, z_batch=zz_batch)
	elif model_name == 'ase_factorization_via_pymc3_binomial_lmm_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz_sample = np.loadtxt(sample_overlap_file)
		#zz_batch = np.loadtxt(batch_overlap_file)
		ase_factorization_obj = ase_factorization_via_pymc3_binomial_lmm_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z_sample=zz_sample)
	elif model_name == 'ase_factorization_via_pymc3_lmm_dirichlet_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz = np.loadtxt(sample_overlap_file)
		ase_factorization_obj = ase_factorization_via_pymc3_lmm_dirichlet_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z=zz)
	elif model_name == 'ase_factorization_via_pymc3_lmm_exponential_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz = np.loadtxt(sample_overlap_file)
		ase_factorization_obj = ase_factorization_via_pymc3_lmm_exponential_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z=zz)
	elif model_name == 'ase_factorization_via_pymc3_lmm_horseshoe_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz = np.loadtxt(sample_overlap_file)
		ase_factorization_obj = ase_factorization_via_pymc3_lmm_horseshoe_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z=zz)
	elif model_name == 'ase_factorization_via_pca':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz = np.loadtxt(sample_overlap_file)
		ase_factorization_obj = ase_factorization_via_pca.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z=zz)
	elif model_name == 'ase_factorization_via_pca_regress_out_cell_line':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz = np.loadtxt(sample_overlap_file)
		ase_factorization_obj = ase_factorization_via_pca_regress_out_cell_line.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept, z=zz)
	elif model_name == 'ase_factorization_via_als':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		zz = np.loadtxt(sample_overlap_file)
		zz_mat = make_cell_line_vector_into_matrix(zz)
		full_cov = np.hstack((cov_plus_intercept, zz_mat))
		ase_factorization_obj = ase_factorization_via_als.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=full_cov)

ase_file = sys.argv[1]
covariate_file = sys.argv[2]
sample_overlap_file = sys.argv[3]
batch_overlap_file = sys.argv[4]
k = int(sys.argv[5])
model_name = sys.argv[6]
output_dir = sys.argv[7]


train_ase_factorization_model(ase_file, covariate_file, sample_overlap_file, batch_overlap_file, k, model_name, output_dir)