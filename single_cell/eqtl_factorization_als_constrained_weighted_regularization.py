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
import pystan
import pickle


#CONSTRAINED_LM = pystan.StanModel(file='update_V_constrained_als.stan')
CONSTRAINED_LM = pickle.load(open('update_V_constrained_als.pkl', 'rb'))
SIMPLEX_CONSTRAINED_LM = pickle.load(open('update_U_constrained_als.pkl', 'rb'))


class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        for fd in self.null_fds + self.save_fds:
            os.close(fd)

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

def outside_update_V_t(g_test, y_test, K, U, lasso_param, genotype_intercept):
	# Get U scaled by genotype for this test
	U_scaled = U*g_test[:,None]
		
	if genotype_intercept == True:
		covariates = np.hstack((np.ones((U_scaled.shape[0],1)), np.transpose(np.asmatrix(g_test)), U_scaled))
		model = sm.OLS(y_test, covariates)
		alpha_param = np.zeros(covariates.shape[1]) + lasso_param
		alpha_param[0] = 0
		alpha_param[1] = 0
		alpha_param[2:] = alpha_param[2:]*(np.sum(U,axis=0)/U.shape[0])
		fit = model.fit_regularized(method='elastic_net', alpha=alpha_param,L1_wt=0.0)
		params = fit.params
	elif genotype_intercept == False:
		covariates = np.hstack((np.ones((U_scaled.shape[0],1)), U_scaled))
		model = sm.OLS(y_test, covariates)
		alpha_param = np.zeros(covariates.shape[1]) + lasso_param
		alpha_param[0] = 0
		alpha_param[1:] = alpha_param[1:]*(np.sum(U,axis=0)/U.shape[0])
		fit = model.fit_regularized(method='elastic_net', alpha=alpha_param,L1_wt=0.0)
		params = fit.params
	return params


def outside_update_U_n(g_sample, y_sample, K, V):
	# Get V scaled by genotype for this sample
	V_scaled = np.transpose(V)*g_sample[:,None]

	data = dict(N=V_scaled.shape[0], K=V_scaled.shape[1], x=V_scaled, y=y_sample)


	try:
		with suppress_stdout_stderr():
			op = SIMPLEX_CONSTRAINED_LM.optimizing(data = data,verbose=False, iter=5000,seed=1)
	except RuntimeError:
		try: 
			op = SIMPLEX_CONSTRAINED_LM.optimizing(data = data,verbose=False, iter=5000, algorithm="Newton",seed=1)
		except RuntimeError:
			try: 
				op = SIMPLEX_CONSTRAINED_LM.optimizing(data = data,verbose=False, iter=5000,seed=2)
			except RuntimeError:
				try: 
					op = SIMPLEX_CONSTRAINED_LM.optimizing(data = data,verbose=False, iter=5000, algorithm="Newton",seed=2)
				except RuntimeError:
					print('Errorr')
	#clf = linear_model.ElasticNet(alpha=lasso_param, l1_ratio=1.0, positive=True, fit_intercept=False)
	#clf.fit(V_scaled, y_sample)

	return op['beta']


class EQTL_FACTORIZATION_ALS_CONSTRAINED(object):
	def __init__(self, K=5,genotype_intercept=True,max_iter=1000, lasso_param_v=.1, parrallel_boolean=False, num_test_cores=24, num_sample_cores=24):
		self.max_iter = max_iter
		self.K = K
		self.lasso_param_v = lasso_param_v
		self.iter = 0
		self.parrallel_boolean = parrallel_boolean
		self.num_sample_cores = num_sample_cores
		self.num_test_cores = num_test_cores
		self.genotype_intercept = genotype_intercept
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
		print(self.lasso_param_v)
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
				np.savetxt('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/temp_output/temper_U_als_constrained_weighted_regularization_' + str(self.lasso_param_v) + '_' + str(self.genotype_intercept) + '_' + str(iter_num) + '.txt', (self.U), fmt="%s", delimiter='\t')
				np.savetxt('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/temp_output/temper_V_als_constrained_weighted_regularization_' + str(self.lasso_param_v) + '_' + str(self.genotype_intercept) + '_' + str(iter_num) + '.txt', (self.V), fmt="%s", delimiter='\t')
				np.savetxt('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/temp_output/temper_intercept_als_constrained_weighted_regularization_' + str(self.lasso_param_v) + '_' + str(self.genotype_intercept) + '_' + str(iter_num) + '.txt', (self.intercept), fmt="%s", delimiter='\t')
				#np.savetxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_psi.txt', (self.psi_alpha/self.psi_beta), fmt="%s", delimiter='\t')

			print('Variational Inference iteration: ' + str(iter_num))
			####################
			end_time = time.time()
			print(end_time-start_time)
			print(self.U[0:10, :])
			print(self.V[:,0:10])
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
				V_update_data.append(outside_update_V_t(self.G[:, test_index], self.Y[:, test_index], self.K, self.U, self.lasso_param_v, self.genotype_intercept))
		elif self.parrallel_boolean == True:
			V_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_V_t)(self.G[:, test_index], self.Y[:, test_index], self.K, self.U, self.lasso_param_v, self.genotype_intercept) for test_index in range(self.T))

		# Convert to array
		V_update_data = np.asarray(V_update_data).T

		self.intercept = V_update_data[0,:]
		self.V = V_update_data[1:,:]
	def update_U(self):
		if self.genotype_intercept == True:
			Y_scaled = self.Y - self.intercept - (self.V[0,:]*self.G)
			temp_V = self.V[1:,:]
		elif self.genotype_intercept == False:
			Y_scaled = self.Y - self.intercept
			temp_V = self.V[:,:]

		U_update_data = []
		###################
		# UPDATE U
		###################
		# Don't parrallelize
		if self.parrallel_boolean == False:
			for sample_index in range(self.N):
				U_update_data.append(outside_update_U_n(self.G[sample_index, :], Y_scaled[sample_index, :], self.K, temp_V))
		# Parrallelize
		elif self.parrallel_boolean == True:
			U_update_data = Parallel(n_jobs=self.num_sample_cores)(delayed(outside_update_U_n)(self.G[sample_index, :], Y_scaled[sample_index, :], self.K, temp_V) for sample_index in range(self.N))

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
		# Initialize randomly onto simplex
		self.U = np.random.random(size=(self.N, self.K))
		for n in range(self.N):
			self.U[n,:] = self.U[n,:]/np.sum(self.U[n,:])
		# self.U = initialize_U_with_known_cell_types(self.U, '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/processed_expression/cell_covariates_sle_individuals_random_subset.txt')
		self.V = np.zeros((self.K, self.T))
		self.intercept = np.zeros(self.T)
		#self.sample_intercept = np.zeros(self.N)
