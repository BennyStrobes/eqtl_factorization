import numpy as np 
import os
import sys
import pdb
import scipy.special as special
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import time
import sklearn.decomposition
from joblib import Parallel, delayed
import multiprocessing
import time

def initialize_U_with_known_cell_types(U, known_cell_types):
	f = open(known_cell_types)
	dicti = {}
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[1])
	f.close()
	uni_arr = np.unique(arr)
	for i, ele in enumerate(uni_arr):
		dicti[ele] = i
	for sample_num in range(U.shape[0]):
		U[sample_num,:] = U[sample_num,:]*0.0
		ct = arr[sample_num]
		U[sample_num,dicti[ct]] = 1.0
	return U

def sigmoid_function(x):
	return 1.0/(1.0 + np.exp(-x))

def outside_update_V_t(g_test, y_test, K, U, lasso_param):
	# Get U scaled by genotype for this test
	U_scaled = U*g_test[:,None]
		
	covariates = np.hstack((np.ones((U_scaled.shape[0],1)), np.transpose(np.asmatrix(g_test)), U_scaled))

	model = sm.OLS(y_test, covariates)
	alpha_param = np.zeros(covariates.shape[1]) + lasso_param
	alpha_param[0] = 0
	alpha_param[1] = 0
	fit = model.fit_regularized(method='elastic_net', alpha=alpha_param,L1_wt=0.0)

	return fit.params


def outside_update_U_n(g_sample, y_sample, K, V, lasso_param):
	# Get V scaled by genotype for this sample
	V_scaled = np.transpose(V)*g_sample[:,None]

	clf = linear_model.ElasticNet(alpha=lasso_param, l1_ratio=1.0, positive=True, fit_intercept=False)
	clf.fit(V_scaled, y_sample)

	return clf.coef_


class EQTL_FACTORIZATION_ALS(object):
	def __init__(self, K=5, max_iter=1000, lasso_param_u=.0001, lasso_param_v=.0001, parrallel_boolean=False, num_test_cores=24, num_sample_cores=24):
		self.max_iter = max_iter
		self.K = K
		self.lasso_param_u = lasso_param_u
		self.lasso_param_v = lasso_param_v
		self.iter = 0
		self.parrallel_boolean = parrallel_boolean
		self.num_sample_cores = num_sample_cores
		self.num_test_cores = num_test_cores
	def fit(self, G, Y, z):
		""" Fit the model.
			Args:
			G: A genotype matrix of floats with shape [num_samples, num_tests].
			Y: An expression matrix of floats with shape [num_samples, num_tests].
			z: groupings of length num_samples
		"""
		self.G_full = G
		self.Y_full = Y
		self.z_full = np.asarray(z)
		self.initialize_variables()
		#self.update_elbo()
		# Loop through VI iterations
		for iter_num in range(self.max_iter):
			start_time = time.time()
			# Update parameter estimaters via ALS
			self.update_V()
			self.update_U()
			self.iter = self.iter + 1
			# Remove irrelevent factors
			if np.mod(iter_num, 50) == 0 and iter_num > 0:
				# UPDATE remove irrelevent_factors TO BE IN TERMS OF *_FULL (ie re-learn theta_U on all data)
				np.savetxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_U_als.txt', (self.U), fmt="%s", delimiter='\t')
				np.savetxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_V_als.txt', (self.V), fmt="%s", delimiter='\t')
				np.savetxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_intercept_als.txt', (self.intercept), fmt="%s", delimiter='\t')
				#np.savetxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_psi.txt', (self.psi_alpha/self.psi_beta), fmt="%s", delimiter='\t')

			print('Variational Inference iteration: ' + str(iter_num))
			####################
			end_time = time.time()
			print(end_time-start_time)
			print(self.U[0:8, :])
			print('##############')
			print('##############')
	def compute_variance_explained_of_factors(self):
		# Based on bottom of P21 of https://arxiv.org/pdf/1802.06931.pdf

		variance_effect = self.N*np.sum(self.tau_beta/self.tau_alpha)
	

		F_terms = self.G*np.dot(np.ones((self.N,1)),[self.F_mu])
		shared_effect = np.sum(np.square(F_terms))

		# Initailize array to keep track of variance explained from each factor
		U_S = self.U_mu*self.S_U
		V_S = self.V_mu
		factor_effects = []
		for k in range(self.K):
			componenent_effects = np.sum(np.square(self.G*(np.dot(np.transpose([U_S[:,k]]), [V_S[k,:]]))))
			factor_effects.append(componenent_effects)
		denominator = np.sum(factor_effects) + shared_effect + variance_effect
		shared_pve = shared_effect/denominator
		factor_pve = factor_effects/denominator
		return shared_pve, factor_pve
	def update_V(self):
		###################
		# UPDATE V
		###################
		# Keep track of variables
		V_update_data = []

		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				V_update_data.append(outside_update_V_t(self.G[:, test_index], self.Y[:, test_index], self.K, self.U, self.lasso_param_v))
		elif self.parrallel_boolean == True:
			V_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_V_t)(self.G[:, test_index], self.Y[:, test_index], self.K, self.U, self.lasso_param_v) for test_index in range(self.T))

		# Convert to array
		V_update_data = np.asarray(V_update_data).T

		self.intercept = V_update_data[0,:]
		self.V = V_update_data[1:,:]
	def update_U(self):
		Y_scaled = self.Y - self.intercept - (self.V[0,:]*self.G)

		U_update_data = []
		###################
		# UPDATE U
		###################
		# Don't parrallelize
		if self.parrallel_boolean == False:
			for sample_index in range(self.N):
				U_update_data.append(outside_update_U_n(self.G[sample_index, :], Y_scaled[sample_index, :], self.K, self.V[1:,:], self.lasso_param_u))
		# Parrallelize
		elif self.parrallel_boolean == True:
			U_update_data = Parallel(n_jobs=self.num_sample_cores)(delayed(outside_update_U_n)(self.G[sample_index, :], Y_scaled[sample_index, :], self.K, self.V[1:,:], self.lasso_param_u) for sample_index in range(self.N))

		# Convert to array
		U_update_data = np.asarray(U_update_data)
		#self.sample_intercept = U_update_data[:,0]
		self.U = U_update_data


	def initialize_variables(self):
		# Add model dimensions to object
		self.N_full = self.Y_full.shape[0]
		self.T_full = self.Y_full.shape[1]

		# Do standard variational inference
		self.N = self.Y_full.shape[0]
		self.T = self.Y_full.shape[1]
		self.G = np.copy(self.G_full)
		self.Y = np.copy(self.Y_full)
		self.z = np.copy(self.z_full)
		self.U = np.random.random(size=(self.N, self.K))
		# self.U = initialize_U_with_known_cell_types(self.U, '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/processed_expression/cell_covariates_sle_individuals_random_subset.txt')
		self.V = np.zeros((self.K, self.T))
		self.intercept = np.zeros(self.T)
		#self.sample_intercept = np.zeros(self.N)
