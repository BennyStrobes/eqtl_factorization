import numpy as np 
import os
import sys
import pdb
import pystan
import pickle


# BB_GLM = pystan.StanModel(file = "betabinomial_glm.stan")
BB_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm.pkl', 'rb'))
# BB_GLM_FIXED_CONC = pystan.StanModel(file = "betabinomial_glm_fixed_conc.stan")
BB_GLM_FIXED_CONC = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm_fixed_conc.pkl', 'rb'))


class ASE_FACTORIZATION(object):
	def __init__(self, K=5, concShape=1.0001, concRate=1e-4, max_iter=1000, output_root='temp'):
		self.max_iter = max_iter
		self.K = K
		self.concShape = concShape
		self.concRate = concRate
		self.iter = 0
		self.output_root = output_root
	def fit(self, allelic_counts, total_counts, cov):
		""" Fit the model.
			Args:
		"""
		self.allelic_counts = allelic_counts
		self.total_counts = total_counts
		self.cov = cov
		self.initialize_variables()
		for vi_iter in range(self.max_iter):
			self.update_V_and_conc()
			self.update_U()
			# Save to output every five iters
			if np.mod(vi_iter, 5) == 0: 
				np.savetxt(self.output_root + '_iter.txt', np.asmatrix(vi_iter), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_U.txt', (self.U), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_V.txt', (self.V), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_conc.txt', (self.conc), fmt="%s", delimiter='\t')

	def update_U(self):
		covariate_predicted = np.dot(self.cov, self.C)
		for sample_iter in range(self.N):
			observed_indices = np.isnan(self.allelic_counts[sample_iter,:]) == False
			merged_intercept = self.V[0,:] + covariate_predicted[sample_iter, :]
			data = dict(N=sum(observed_indices), P=self.U.shape[1], x=np.transpose(self.V[1:,observed_indices]), intercept=merged_intercept[observed_indices], ys=self.allelic_counts[sample_iter, observed_indices].astype(int), ns=self.total_counts[sample_iter, observed_indices].astype(int), conc=self.conc[observed_indices])
			########################################
			# Initialize
			beta_init = np.zeros(data['P'])
			init = dict(beta=beta_init)
			########################################
			# Run optimization
			pdb.set_trace()
			op = BB_GLM_FIXED_CONC.optimizing(data = data, verbose=False,iter=5000,seed=1, init=init)
			self.U[sample_iter, :] = op['beta']

	def update_V_and_conc(self):
		# Add column of ones (intercept to U)
		X0 = np.ones((self.N,1))
		U_new = np.hstack((X0, self.U, self.cov))
		for test_iter in range(self.T):
			observed_indices = np.isnan(self.allelic_counts[:,test_iter]) == False
			data = dict(N=sum(observed_indices), P=U_new.shape[1], x=U_new[observed_indices,:], ys=self.allelic_counts[observed_indices, test_iter].astype(int), ns=self.total_counts[observed_indices, test_iter].astype(int), concShape=self.concShape, concRate=self.concRate)
			########################################
			# Get MOM estimates for initialization
			rat = data['ys']/data['ns'].astype(float)
			if np.sum(np.isnan(rat)) > 0:
				pdb.set_trace()
			# moment estimator of concentration parameter
			conc_init = min(1.0/np.var(rat), 1000.0)
			# thresholded moment estimator of the mean 
			m_init = min(max(np.mean(rat), 1.0/1000 ), 1.0-(1.0/1000))
			beta_init = np.zeros(data['P'])
			beta_init[0] = np.log(m_init/(1.0-m_init))
			init = dict(conc=conc_init, beta=beta_init)
			########################################
			# Run optimization
			op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1, init=init)
			self.V[:, test_iter] = op['beta'][:(self.K+1)]
			self.C[:, test_iter] = op['beta'][(self.K+1):]
			self.conc[test_iter] = np.asmatrix(op['conc'])[0,0]
			pdb.set_trace()

	def initialize_variables(self):
		self.N = self.allelic_counts.shape[0]
		self.T = self.allelic_counts.shape[1]
		self.num_cov = self.cov.shape[1]
		# Randomly initialize (only U matters)
		self.U = np.random.randn(self.N, self.K)
		self.V = np.random.randn((self.K + 1), self.T)  # plus 1 for intercept
		self.conc = np.random.randn(self.T)
		self.C = np.random.randn(self.num_cov, self.T)
