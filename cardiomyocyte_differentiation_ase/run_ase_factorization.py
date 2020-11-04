import numpy as np 
import os
import sys
import pdb
#import ase_factorization
#import ase_factorization_via_stan_vb
import ase_factorization_via_pymc3_vb
import ase_factorization_via_pymc3_mvn_vb
import ase_factorization_via_pymc3_lmm_vb


def load_in_ase_data(ase_file):
	full_data = np.loadtxt(ase_file, dtype=str, delimiter='\t')
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
				allelic_counts[ii,jj] = int(count_data[ii,jj].split('/')[0])
				total_counts[ii,jj] = int(count_data[ii,jj].split('/')[1])
	return np.transpose(allelic_counts), np.transpose(total_counts)

def add_intercept_column_to_matrix(X):
	n,m = X.shape # for generality
	X0 = np.ones((n,1))
	Xnew = np.hstack((X0, X))
	return Xnew

def train_ase_factorization_model(ase_file, covariate_file, sample_overlap_file, k, model_name, output_dir):
	if model_name == 'ase_factorization':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		ase_factorization_obj = ase_factorization.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept)
	elif model_name == 'ase_factorization_via_stan_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		ase_factorization_obj = ase_factorization_via_stan_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept)
	elif model_name == 'ase_factorization_via_pymc3_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		ase_factorization_obj = ase_factorization_via_pymc3_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept)
	elif model_name == 'ase_factorization_via_pymc3_mvn_vb':
		allelic_counts, total_counts = load_in_ase_data(ase_file)
		if covariate_file != 'NA':
			cov = np.loadtxt(covariate_file)
			cov_plus_intercept = add_intercept_column_to_matrix(cov)
		else:
			cov_plus_intercept = np.ones((allelic_counts.shape[0], 1))
		ase_factorization_obj = ase_factorization_via_pymc3_mvn_vb.ASE_FACTORIZATION(K=k, output_root=output_dir + '_ase_factorization')
		ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts, cov=cov_plus_intercept)
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

ase_file = sys.argv[1]
covariate_file = sys.argv[2]
sample_overlap_file = sys.argv[3]
k = int(sys.argv[4])
model_name = sys.argv[5]
output_dir = sys.argv[6]


train_ase_factorization_model(ase_file, covariate_file, sample_overlap_file, k, model_name, output_dir)