import numpy as np 
import os
import sys
import pdb
import ase_factorization

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

def train_ase_factorization_model(ase_file, k, output_dir):
	allelic_counts, total_counts = load_in_ase_data(ase_file)
	ase_factorization_obj = ase_factorization.ASE_FACTORIZATION(K=k, output_root=output_dir + 'ase_factorization')
	ase_factorization_obj.fit(allelic_counts=allelic_counts, total_counts=total_counts)




ase_file = sys.argv[1]
k = int(sys.argv[2])
output_dir = sys.argv[3]


train_ase_factorization_model(ase_file, k,  output_dir)