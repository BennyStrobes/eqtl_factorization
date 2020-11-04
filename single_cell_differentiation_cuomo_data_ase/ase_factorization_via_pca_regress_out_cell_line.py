import numpy as np 
import os
import sys
import pdb
import pickle
from ppca import PPCA
from sklearn.linear_model import LinearRegression

# scale each column
def scale_allelic_ratios(rat):
	scaled_rat = np.copy(rat)
	N, P = scaled_rat.shape
	# loop through columns
	for column_index in range(P):
		column_mean = np.nanmean(rat[:,column_index])
		column_std = np.nanstd(rat[:,column_index])
		for sample_index in range(N):
			if np.isnan(rat[sample_index, column_index]) == False:
				standardized_rat = (rat[sample_index, column_index] - column_mean)/column_std
				scaled_rat[sample_index, column_index] = standardized_rat
	return scaled_rat

def regress_out_cell_line(scaled_rat, Z):
	# First make 1 hot encoding of Z
	num_cell_lines = len(np.unique(Z))
	num_cells = len(Z)
	Z_mat = np.zeros((num_cells, num_cell_lines-1))
	for n in range(num_cells):
		line_index = Z[n]
		if line_index != num_cell_lines-1:
			Z_mat[n, line_index] = 1.0
	final_rat = np.copy(scaled_rat)
	N, P = scaled_rat.shape
	for column_index in range(P):
		valid_row_indices = np.isnan(scaled_rat[:,column_index]) == False
		reg = LinearRegression().fit(Z_mat[valid_row_indices,:], scaled_rat[valid_row_indices,column_index])
		predicted = reg.predict(Z_mat[valid_row_indices,:])
		final_rat[valid_row_indices,column_index] = scaled_rat[valid_row_indices,column_index] - predicted
	return final_rat
def filter_lowly_expressed_sites(k, n, count_thresh):
		N, S = k.shape
		new_rat = np.zeros((N,S))
		for row_num in range(N):
			for col_num in range(S):
				if np.isnan(n[row_num, col_num]):
					new_rat[row_num, col_num] = np.nan
				elif n[row_num,col_num] < count_thresh:
					new_rat[row_num, col_num] = np.nan
				else:
					new_rat[row_num, col_num] = n[row_num, col_num]
		return new_rat

# BB_GLM = pystan.StanModel(file = "betabinomial_glm.stan")
#BB_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm.pkl', 'rb'))
# BB_GLM_FIXED_CONC = pystan.StanModel(file = "betabinomial_glm_fixed_conc.stan")
#BB_GLM_FIXED_CONC = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm_fixed_conc.pkl', 'rb'))
#BB_FACTORIZATION = pystan.StanModel(file='/work-zfs/abattle4/bstrober/temp/ase_factorization.stan')

class ASE_FACTORIZATION(object):
	def __init__(self, K=5, concShape=1.0001, concRate=1e-4, max_iter=1000, output_root='temp'):
		self.max_iter = max_iter
		self.K = K
		self.concShape = concShape
		self.concRate = concRate
		self.iter = 0
		self.output_root = output_root + '3'
	def fit(self, allelic_counts, total_counts, cov, z):
		""" Fit the model.
			Args:
		"""
		self.allelic_counts = allelic_counts
		self.total_counts = total_counts
		self.cov = cov
		self.Z = z.astype(int)
		N = self.allelic_counts.shape[0]
		T = self.allelic_counts.shape[1]
		I = len(np.unique(self.Z))

		self.run_factorization(N, T, self.cov, self.Z, I, self.K, self.allelic_counts, self.total_counts)
		
		#data = dict(N=N, T=T, K=self.K, num_cov=self.num_cov, cov=self.cov, ys=np.transpose(self.allelic_counts.astype(int)), ns=np.transpose(self.total_counts.astype(int)), concShape=self.concShape, concRate=self.concRate)
		#print("START")
		#aa = BB_FACTORIZATION.vb(data=data)
		#pickle.dump(aa, open(self.output_root + '_model', 'wb'))

	def run_factorization(self, N, S, X, Z, I, K, k, n):
		# Smart initialization
		#rat2 = k/n
		rat = filter_lowly_expressed_sites(k,n, 7)
		nans = np.isnan(rat)

		scaled_rat = scale_allelic_ratios(rat)
		scaled_residual_rat = regress_out_cell_line(scaled_rat, Z)
		rescaled_residual_rat = scale_allelic_ratios(scaled_residual_rat)
		ppca = PPCA()
		ppca.fit(data=np.transpose(rescaled_residual_rat), d=K, verbose=True, tol=1e-6)
		U = ppca.C
		V = ppca.transform()
		pickle.dump(ppca, open(self.output_root + '_model', 'wb'))
		np.savetxt(self.output_root + '_temper_U.txt', U, fmt="%s", delimiter='\t')
		np.savetxt(self.output_root + '_temper_V.txt', V.T, fmt="%s", delimiter='\t')




