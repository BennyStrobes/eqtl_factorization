import numpy as np 
import os
import sys
import pdb
import pickle
import pymc3 as pm
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

class ASE_FACTORIZATION(object):
	def __init__(self, K=5, concShape=1.0001, concRate=1e-4, max_iter=1000, output_root='temp'):
		self.max_iter = max_iter
		self.K = K
		self.concShape = concShape
		self.concRate = concRate
		self.iter = 0
		self.output_root = output_root + '10'
	def fit(self, allelic_counts, total_counts, cov, z_sample, z_batch):
		""" Fit the model.
			Args:
		"""
		self.allelic_counts = allelic_counts
		self.total_counts = total_counts
		self.cov = cov
		self.Z_sample = z_sample.astype(int)
		self.Z_batch = z_batch.astype(int)
		self.initialize_variables()

		self.run_factorization()
		
		#data = dict(N=N, T=T, K=self.K, num_cov=self.num_cov, cov=self.cov, ys=np.transpose(self.allelic_counts.astype(int)), ns=np.transpose(self.total_counts.astype(int)), concShape=self.concShape, concRate=self.concRate)
		#print("START")
		#aa = BB_FACTORIZATION.vb(data=data)
		#pickle.dump(aa, open(self.output_root + '_model', 'wb'))

	def run_factorization(self):
		rat = self.allelic_counts/self.total_counts
		nans = np.isnan(rat)
		# Run bb-mf
		with pm.Model() as bb_glm:
			CONC = pm.HalfCauchy('CONC', beta=5, shape=(1,self.S), testval=self.conc_init)
			BETA = pm.Normal('BETA', mu=0, tau=(1/1000000.0), shape=(self.S, self.num_cov), testval=self.beta_init)
			U = pm.Normal('U', mu=0, tau=(1.0/1.0), shape=(self.N, self.K), testval=self.U_init)
			V = pm.Normal('V', mu=0, tau=(1.0/1.0), shape=(self.S, self.K), testval=self.V_init)

			MU_A_SAMPLE = pm.Normal("MU_A_SAMPLE", mu=0., sd=100**2, shape=(1,self.S), testval=self.mu_a_sample_init)
			SIGMA_A_SAMPLE = pm.HalfCauchy("SIGMA_A_SAMPLE", beta=5.0, shape=(1,self.S), testval=self.sigma_a_sample_init)
			mu_a_sample_mat = pm.math.dot(np.ones((self.I_sample,1)), MU_A_SAMPLE)
			sigma_a_sample_mat = pm.math.dot(np.ones((self.I_sample,1)), SIGMA_A_SAMPLE)
			A_SAMPLE = pm.Normal('A_SAMPLE', mu=mu_a_sample_mat, sigma=sigma_a_sample_mat, shape=(self.I_sample,self.S), testval=self.A_sample_init)

			MU_A_BATCH = pm.Normal("MU_A_BATCH", mu=0., sd=100**2, shape=(1,self.S), testval=self.mu_a_batch_init)
			SIGMA_A_BATCH = pm.HalfCauchy("SIGMA_A_BATCH", beta=5.0, shape=(1,self.S), testval=self.sigma_a_batch_init)
			mu_a_batch_mat = pm.math.dot(np.ones((self.I_batch,1)), MU_A_BATCH)
			sigma_a_batch_mat = pm.math.dot(np.ones((self.I_batch,1)), SIGMA_A_BATCH)
			A_BATCH = pm.Normal('A_BATCH', mu=mu_a_batch_mat, sigma=sigma_a_batch_mat, shape=(self.I_batch,self.S), testval=self.A_batch_init)

			p = pm.math.invlogit(pm.math.dot(self.cov, BETA.T) + pm.math.dot(U,V.T) + A_SAMPLE[self.Z_sample,:] + A_BATCH[self.Z_batch,:])
			conc_mat = pm.math.dot(np.ones((self.N,1)), CONC)
			R = pm.BetaBinomial('like',alpha=(p*conc_mat)[~nans], beta=((1.0-p)*conc_mat)[~nans], n=self.total_counts[~nans], observed=self.allelic_counts[~nans])
			approx = pm.fit(method='advi', n=100000)
		pickle.dump(approx, open(self.output_root + '_model', 'wb'))
		#approx = pickle.load( open(self.output_root + '_model', "rb" ) )
		means_dict = approx.bij.rmap(approx.params[0].eval())
		np.savetxt(self.output_root + '_temper_U.txt', (means_dict['U']), fmt="%s", delimiter='\t')
		np.savetxt(self.output_root + '_temper_V.txt', (means_dict['V'].T), fmt="%s", delimiter='\t')
		np.savetxt(self.output_root + '_temper_BETA.txt', (means_dict['BETA'].T), fmt="%s", delimiter='\t')

	def initialize_variables(self):
		self.I_sample = len(np.unique(self.Z_sample))
		self.I_batch = len(np.unique(self.Z_batch))
		self.N = self.allelic_counts.shape[0]
		self.S = self.allelic_counts.shape[1]
		self.num_cov = self.cov.shape[1]
		# Randomly initialize (only U matters)
		self.U_init = np.random.randn(self.N, self.K)
		self.V_init = np.random.randn(self.S, self.K)  # plus 1 for intercept
		rat = self.allelic_counts/self.total_counts
		nans = np.isnan(rat)
		self.conc_init = np.zeros((1, self.S))
		self.beta_init = np.zeros((self.num_cov, self.S))
		for index_s in range(self.S):
			column_rat = rat[:, index_s]
			column_nans = np.isnan(column_rat)
			valid_rat = column_rat[~column_nans]
			conc_init_s = min(1.0/np.var(valid_rat), 1000.0)
			m_init = min(max(np.mean(valid_rat), 1.0/1000 ), 1.0-(1.0/1000))
			self.conc_init[0, index_s] = conc_init_s
			self.beta_init[0, index_s] = np.log(m_init/(1.0-m_init))
		self.beta_init = np.transpose(self.beta_init)
		self.A_sample_init = np.zeros((self.I_sample, self.S))
		self.sigma_a_sample_init = np.ones((1,self.S))
		self.mu_a_sample_init = np.zeros((1, self.S))
		self.A_batch_init = np.zeros((self.I_batch, self.S))
		self.sigma_a_batch_init = np.ones((1,self.S))
		self.mu_a_batch_init = np.zeros((1, self.S))

		'''
		for n in range(self.N):
			for t in range(self.T):
				if np.isnan(self.allelic_counts[n,t]) or np.isnan(self.total_counts[n,t]):
					self.allelic_counts[n,t] = 0
					self.total_counts[n,t] = 0
		'''



