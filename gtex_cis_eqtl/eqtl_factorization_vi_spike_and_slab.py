import numpy as np 
import os
import sys
import pdb
import scipy.special as special
from sklearn.linear_model import LinearRegression
import time
import sklearn.decomposition
from joblib import Parallel, delayed
import multiprocessing
import time

def sigmoid_function(x):
	return 1.0/(1.0 + np.exp(-x))

def weighted_SVI_updated(old_parameter, new_parameter, step_size):
	updated_parameter = ((1.0 - step_size)*old_parameter) + (step_size*new_parameter)
	return updated_parameter

def run_linear_model_for_initialization(Y, G):
	num_tests = Y.shape[1]
	betas = []
	for test_number in range(num_tests):
		y_vec = Y[:,test_number]
		g_vec = G[:,test_number]
		reg = LinearRegression().fit(np.transpose(np.asmatrix(g_vec)), np.transpose(np.asmatrix(y_vec)))
		betas.append(reg.coef_[0][0])
	return np.asarray(betas)


def compute_kl_divergence_of_gaussian_bernoulli(S, W_mu, W_var, W_var_s_0, gamma_expected, theta_a, theta_b, K):
	num_feat = W_mu.shape[1]
	# Relevent expectations
	#log_gamma_expected = special.digamma(gamma_alpha) - np.log(gamma_beta)
	#gamma_expected = gamma_alpha/gamma_beta
	log_theta_expected_val = special.digamma(theta_a) - special.digamma(theta_a + theta_b)
	log_1_minus_theta_expected_val = special.digamma(theta_b) - special.digamma(theta_a + theta_b)
	#W_var_s_0_temp = np.dot(np.transpose([(gamma_beta/gamma_alpha)]),np.ones((1,num_feat)))
	W_squared_expected_val = (S*(np.square(W_mu) + W_var)) + ((1.0-S)*W_var_s_0)

	# Initialize variables
	likelihood_term_a = 0
	likelihood_term_b = 0
	likelihood_term_c = 0
	likelihood_term_d = 0
	entropy_term_a = 0
	entropy_term_b = 0

	for k in range(K):
		#likelihood_term_a = likelihood_term_a + (num_feat/2.0)*log_gamma_expected[k]
		likelihood_term_b = likelihood_term_b - (np.sum(W_squared_expected_val[k,:])*gamma_expected/2.0)
		likelihood_term_c = likelihood_term_c + (np.sum(S[k,:])*log_theta_expected_val[k])
		likelihood_term_d = likelihood_term_d + (np.sum(1.0-S[k,:])*log_1_minus_theta_expected_val[k])
		#entropy_term_a = entropy_term_a + (.5)*np.sum(np.log((S[k,:]*W_var[k,:]) + ((1-S[k,:])*W_var_s_0[k,:])))
		entropy_term_a = entropy_term_a - (.5)*np.sum(S[k,:]*np.log(W_var[k,:]) + (1.0-S[k,:])*np.log(W_var_s_0[k,:]))

		temp_term_b = ((1.0-S[k,:])*np.log(1.0-S[k,:])) + (S[k,:]*np.log(S[k,:]))
		temp_term_b[np.isnan(temp_term_b)] = 0.
		entropy_term_b = entropy_term_b + np.sum(temp_term_b)
		#entropy_term_b = entropy_term_b + np.sum(((1.0-S[k,:])*np.log(log_const+1.0-S[k,:])) - (S[k,:]*np.log(log_const+S[k,:])))
	
	kl_divergence = entropy_term_a + entropy_term_b - likelihood_term_a - likelihood_term_b - likelihood_term_c - likelihood_term_d

	return kl_divergence

def compute_kl_divergence_of_gaussian(W_mu, W_var, gamma_alpha, gamma_beta, K):
	num_feat = W_mu.shape[1]
	# Relevent expectations
	log_gamma_expected = special.digamma(gamma_alpha) - np.log(gamma_beta)
	gamma_expected = gamma_alpha/gamma_beta
	#log_theta_expected_val = special.digamma(theta_a) - special.digamma(theta_a + theta_b)
	#log_1_minus_theta_expected_val = special.digamma(theta_b) - special.digamma(theta_a + theta_b)
	#W_var_s_0_temp = np.dot(np.transpose([(gamma_beta/gamma_alpha)]),np.ones((1,num_feat)))
	W_squared_expected_val = ((np.square(W_mu) + W_var))

	# Initialize variables
	likelihood_term_a = 0
	likelihood_term_b = 0
	likelihood_term_c = 0
	likelihood_term_d = 0
	entropy_term_a = 0
	entropy_term_b = 0

	for k in range(K):
		likelihood_term_a = likelihood_term_a + (num_feat/2.0)*log_gamma_expected[k]
		likelihood_term_b = likelihood_term_b - (np.sum(W_squared_expected_val[k,:])*gamma_expected[k]/2.0)
		#likelihood_term_c = likelihood_term_c + (np.sum(S[k,:])*log_theta_expected_val[k])
		#likelihood_term_d = likelihood_term_d + (np.sum(1.0-S[k,:])*log_1_minus_theta_expected_val[k])
		#entropy_term_a = entropy_term_a + (.5)*np.sum(np.log((S[k,:]*W_var[k,:]) + ((1-S[k,:])*W_var_s_0[k,:])))
		entropy_term_a = entropy_term_a - (.5)*np.sum(np.log(W_var[k,:]))
		#temp_term_b = ((1.0-S[k,:])*np.log(1.0-S[k,:])) + (S[k,:]*np.log(S[k,:]))
		#temp_term_b[np.isnan(temp_term_b)] = 0.
		#entropy_term_b = entropy_term_b + np.sum(temp_term_b)
		#entropy_term_b = entropy_term_b + np.sum(((1.0-S[k,:])*np.log(log_const+1.0-S[k,:])) - (S[k,:]*np.log(log_const+S[k,:])))
	
	kl_divergence = entropy_term_a  - likelihood_term_a - likelihood_term_b - likelihood_term_c - likelihood_term_d

	return kl_divergence

def compute_kl_divergence_of_gaussian_fixed_variance(W_mu, W_var, gamma_expected, K):
	num_feat = W_mu.shape[1]
	# Relevent expectations
	#log_gamma_expected = special.digamma(gamma_alpha) - np.log(gamma_beta)
	#gamma_expected = gamma_alpha/gamma_beta
	#log_theta_expected_val = special.digamma(theta_a) - special.digamma(theta_a + theta_b)
	#log_1_minus_theta_expected_val = special.digamma(theta_b) - special.digamma(theta_a + theta_b)
	#W_var_s_0_temp = np.dot(np.transpose([(gamma_beta/gamma_alpha)]),np.ones((1,num_feat)))
	W_squared_expected_val = ((np.square(W_mu) + W_var))

	# Initialize variables
	likelihood_term_a = 0
	likelihood_term_b = 0
	likelihood_term_c = 0
	likelihood_term_d = 0
	entropy_term_a = 0
	entropy_term_b = 0

	for k in range(K):
		#likelihood_term_a = likelihood_term_a + (num_feat/2.0)*log_gamma_expected[k]
		likelihood_term_b = likelihood_term_b - (np.sum(W_squared_expected_val[k,:])*gamma_expected/2.0)
		#likelihood_term_c = likelihood_term_c + (np.sum(S[k,:])*log_theta_expected_val[k])
		#likelihood_term_d = likelihood_term_d + (np.sum(1.0-S[k,:])*log_1_minus_theta_expected_val[k])
		#entropy_term_a = entropy_term_a + (.5)*np.sum(np.log((S[k,:]*W_var[k,:]) + ((1-S[k,:])*W_var_s_0[k,:])))
		entropy_term_a = entropy_term_a - (.5)*np.sum(np.log(W_var[k,:]))
		#temp_term_b = ((1.0-S[k,:])*np.log(1.0-S[k,:])) + (S[k,:]*np.log(S[k,:]))
		#temp_term_b[np.isnan(temp_term_b)] = 0.
		#entropy_term_b = entropy_term_b + np.sum(temp_term_b)
		#entropy_term_b = entropy_term_b + np.sum(((1.0-S[k,:])*np.log(log_const+1.0-S[k,:])) - (S[k,:]*np.log(log_const+S[k,:])))
	
	kl_divergence = entropy_term_a  - likelihood_term_b - likelihood_term_c - likelihood_term_d

	return kl_divergence



def compute_kl_divergence_of_gamma(alpha_prior, beta_prior, gamma_alpha, gamma_beta):
	# Relevent expectations
	log_gamma_expected = special.digamma(gamma_alpha) - np.log(gamma_beta)
	gamma_expected = gamma_alpha/gamma_beta
	# Compute kl divergence
	likelihood_term = np.sum(alpha_prior*np.log(beta_prior) + (alpha_prior-1.0)*log_gamma_expected - beta_prior*gamma_expected - special.gammaln(alpha_prior))
	entropy_term = np.sum(gamma_alpha*np.log(gamma_beta) + (gamma_alpha-1.0)*log_gamma_expected - gamma_beta*gamma_expected - special.gammaln(gamma_alpha))
	kl_divergence = entropy_term - likelihood_term
	return kl_divergence

def compute_kl_divergence_of_beta(a_prior, b_prior, theta_a, theta_b):
	# Relevent expectations
	ln_theta_expected_val = special.digamma(theta_a) - special.digamma(theta_a + theta_b)
	ln_1_minus_theta_expected_val = special.digamma(theta_b) - special.digamma(theta_a + theta_b)
	# Compuate kl divergence
	likelihood_term = np.sum((a_prior-1.0)*ln_theta_expected_val + (b_prior-1.0)*ln_1_minus_theta_expected_val - special.betaln(a_prior, b_prior))
	entropy_term = np.sum((theta_a-1.0)*ln_theta_expected_val + (theta_b-1.0)*ln_1_minus_theta_expected_val - special.betaln(theta_a, theta_b))
	kl_divergence = entropy_term - likelihood_term
	return kl_divergence


def outside_update_U_n(U_mu, S_U, U_var, U_var_s_0, G_slice, Y_slice, K, V_S_expected_val, V_S_squared_expected_val, F_S_expected_val, intercept_mu, gamma_v, tau_expected_val, theta_U_a, theta_U_b):
	for k in range(K):
		# Compute relevent expectations
		U_S_expected_val = U_mu*S_U
		V_k_S_k_expected_val = V_S_expected_val[k,:]
		theta_U_expected_val = theta_U_a[k]/(theta_U_a[k] + theta_U_b[k])
		ln_theta_U_expected_val = special.digamma(theta_U_a[k]) - special.digamma(theta_U_a[k]+theta_U_b[k])  # expectation of ln(1-X)
		ln_1_minus_theta_U_expected_val = special.digamma(theta_U_b[k]) - special.digamma(theta_U_a[k]+theta_U_b[k])
		# Compute expectations on other components
		other_components_expected = (U_S_expected_val@V_S_expected_val) - U_S_expected_val[k]*V_S_expected_val[k,:]
		# Update variance of q(U|s=1)
		a_term = np.sum(tau_expected_val*np.square(G_slice)*V_S_squared_expected_val[k,:]) + gamma_v
		U_var[k] = 1.0/a_term
		# Update variance of q(U|s=0)
		U_var_s_0[k] = 1.0/gamma_v
		# Update mean of q(U|s=1)
		resid = Y_slice - intercept_mu - G_slice*(F_S_expected_val + other_components_expected)
		b_term = np.sum(tau_expected_val*G_slice*V_k_S_k_expected_val*resid)
		U_mu[k] = U_var[k]*b_term
		# Now update q(S_U=1)
		z_term = ln_theta_U_expected_val - ln_1_minus_theta_U_expected_val + .5*np.log(gamma_v) - .5*np.log(a_term) + (np.square(b_term)/(2.0*a_term))
		S_U[k] = sigmoid_function(z_term)
	return np.hstack((U_mu, S_U, U_var, U_var_s_0))

def outside_update_V_t(V_mu, V_var, G_slice, Y_slice, K, U_S_expected_val, U_S_squared_expected_val, F_S_t_expected_val, intercept_mu, gamma_v, tau_t_expected_val, sample_batch_fraction, step_size, SVI):
	for k in range(K):
		# Compute expectations on other components
		other_components_expected = (U_S_expected_val@V_mu) - U_S_expected_val[:, k]*V_mu[k]
		# Update variance of q(V|s=1)
		a_term = gamma_v + (1.0/sample_batch_fraction)*(tau_t_expected_val*np.sum(np.square(G_slice)*U_S_squared_expected_val[:,k]))
		# Update mean of q(U|s=1)
		resid = Y_slice - intercept_mu - G_slice*(other_components_expected + F_S_t_expected_val)
		b_term = (1.0/sample_batch_fraction)*np.sum(tau_t_expected_val*G_slice*U_S_expected_val[:,k]*resid)

		new_var = 1.0/a_term
		new_mu = new_var*b_term
		if SVI == False:
			V_var[k] = new_var
			V_mu[k] = new_mu
		elif SVI == True:
			V_var[k] = weighted_SVI_updated(V_var[k], new_var, step_size)
			V_mu[k] = weighted_SVI_updated(V_mu[k], new_mu, step_size)
	return np.hstack((V_mu, V_var))

def outside_update_intercept_t(intercept_mu, intercept_var, G_slice, Y_slice, N, U_S_expected_val, V_S_t_expected_val, F_S_t_expected_val, tau_t_expected_val, sample_batch_fraction, step_size, SVI):
	# Compute relevent expectations
	# Compute expectations on other components
	other_components_expected = U_S_expected_val@V_S_t_expected_val
	resid = Y_slice - G_slice*(F_S_t_expected_val + other_components_expected)

	new_var = 1.0/((1.0/sample_batch_fraction)*N*tau_t_expected_val)
	new_mu = new_var*tau_t_expected_val*np.sum(resid)*(1.0/sample_batch_fraction)

	if SVI == False:
		intercept_var = new_var
		intercept_mu = new_mu
	elif SVI == True:
		intercept_var = weighted_SVI_updated(intercept_var, new_var, step_size)
		intercept_mu = weighted_SVI_updated(intercept_mu, new_mu, step_size)
	return np.hstack((intercept_mu, intercept_var))

def outside_update_F_t(F_mu, F_var, G_slice, Y_slice, U_S_expected_val, V_S_t_expected_val, intercept_mu_t, gamma_f_expected_val, tau_t_expected_val, sample_batch_fraction, step_size, SVI):
	# Compute expectations on other components
	other_components_expected = U_S_expected_val@V_S_t_expected_val

	# Update variance of q(F|s=1)
	a_term = gamma_f_expected_val + (1.0/sample_batch_fraction)*tau_t_expected_val*np.sum(np.square(G_slice))
	# Update mean of q(F|s=1)
	resid = Y_slice - intercept_mu_t - G_slice*(other_components_expected)
	b_term = np.sum(tau_t_expected_val*G_slice*resid*(1.0/sample_batch_fraction))
	new_var = 1.0/a_term
	new_mu = new_var*b_term
	if SVI == False:
		F_var = new_var
		F_mu = new_mu
	elif SVI == True:
		F_var = weighted_SVI_updated(F_var, new_var, step_size)
		F_mu = weighted_SVI_updated(F_mu, new_mu, step_size)
	return np.hstack((F_mu, F_var))

def outside_update_tau_t(tau_alpha, tau_beta, G_slice, Y_slice, N, U_S, V_S_t, F_S_t, intercept_t, V_S_t_squared, F_S_t_squared, U_S_squared, intercept_t_squared, alpha_prior, beta_prior, sample_batch_fraction, step_size, SVI):
	# Compute Relevent expectations
	squared_factor_terms = U_S_squared@V_S_t_squared
	factor_terms = U_S@V_S_t

	# First add together square terms
	resid = np.square(Y_slice) + intercept_t_squared + np.square(G_slice)*(F_S_t_squared + squared_factor_terms)
	# Now add terms with Y
	resid = resid - (2.0*Y_slice*(intercept_t + G_slice*factor_terms + G_slice*F_S_t))

	resid = resid + 2.0*intercept_t*(G_slice*(factor_terms + F_S_t))
	# Now add terms with factors
	resid = resid + 2.0*G_slice*factor_terms*G_slice*F_S_t
	# Now add terms with interactions between factors
	resid = resid + (np.square(G_slice)*(factor_terms*factor_terms - np.sum(np.square(U_S*V_S_t),axis=1)))

	# Make Updates
	new_alpha = alpha_prior + ((N/2.0)*(1.0/sample_batch_fraction))
	new_beta = beta_prior + ((np.sum(resid)/2.0)*(1.0/sample_batch_fraction))

	if SVI == False:
		tau_alpha = new_alpha
		tau_beta = new_beta
	elif SVI == True:
		tau_alpha = weighted_SVI_updated(tau_alpha, new_alpha, step_size)
		tau_beta = weighted_SVI_updated(tau_beta, new_beta, step_size)
	return np.hstack((tau_alpha, tau_beta))


class EQTL_FACTORIZATION_VI(object):
	def __init__(self, K=25, alpha=1e-3, beta=1e-3, a=1, b=1, gamma_v=1.0, max_iter=1000, delta_elbo_threshold=1e-8, SVI=False, parrallel_boolean=False, sample_batch_fraction=.3, learning_rate=.4, forgetting_rate=.01, num_test_cores=24, num_sample_cores=24):
		self.alpha_prior = alpha
		self.beta_prior = beta
		self.a_prior = a 
		self.b_prior = b
		self.max_iter = max_iter
		self.K = K
		self.gamma_v = gamma_v
		self.iter = 0
		self.delta_elbo_threshold = delta_elbo_threshold
		self.SVI = SVI
		self.parrallel_boolean = parrallel_boolean
		self.num_sample_cores = num_sample_cores
		self.num_test_cores = num_test_cores
		if SVI == True:
			self.sample_batch_fraction = sample_batch_fraction
			self.learning_rate = learning_rate
			self.forgetting_rate = forgetting_rate
		elif SVI == False:
			self.sample_batch_fraction = 1.0
			self.step_size = 1.0
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
		self.update_elbo()
		# Loop through VI iterations
		for vi_iter in range(self.max_iter):
			start_time = time.time()
			# Update step size (If using SVI)
			self.update_step_size()
			# Update parameter estimaters via coordinate ascent
			print('Step size: ' + str(self.step_size))
			self.update_U()
			self.update_V()
			self.update_intercept()
			self.update_F()
			self.update_theta_U()
			self.update_tau()
			self.iter = self.iter + 1
			# Remove irrelevent factors
			if np.mod(vi_iter, 50) == 0 and vi_iter > 0:
				# UPDATE remove irrelevent_factors TO BE IN TERMS OF *_FULL (ie re-learn theta_U on all data)
				self.remove_irrelevent_factors()

			
			# Compute ELBO after update
			print('Variational Inference iteration: ' + str(vi_iter))
			'''
			if np.mod(vi_iter, 20) == 0:
				self.update_elbo()
				current_elbo = self.elbo[len(self.elbo)-1]
				delta_elbo = (current_elbo - self.elbo[len(self.elbo)-2])
				print('delta ELBO: ' + str(delta_elbo))
			'''
			####################
			print(self.theta_U_a/(self.theta_U_a + self.theta_U_b))
			end_time = time.time()
			print(end_time-start_time)
			print('##############')
			print('##############')
	def update_step_size(self):
		# Only needs to be done for SVI
		if self.SVI == True:
			self.step_size = self.learning_rate/(np.power((1.0 + (self.forgetting_rate*self.iter)), .75))
	def remove_irrelevent_factors(self):
		#shared_pve, factor_pve = self.compute_variance_explained_of_factors()
		factor_sparsity = self.theta_U_a/(self.theta_U_a + self.theta_U_b)
		#num_factors = len(np.where(factor_sparsity > 0.01)[0])
		#factor_ordering = np.flip(np.argsort(factor_sparsity))[:num_factors]
		factor_ordering = np.where(factor_sparsity > 0.01)[0]
		print(factor_ordering)
		self.U_mu = self.U_mu[:, factor_ordering]
		self.U_var = self.U_var[:, factor_ordering]
		self.U_var_s_0 = self.U_var_s_0[:, factor_ordering]
		self.S_U = self.S_U[:, factor_ordering]
		if self.SVI == True:
			self.U_mu_full = self.U_mu_full[:, factor_ordering]
			self.U_var_full = self.U_var_full[:, factor_ordering]
			self.U_var_s_0_full = self.U_var_s_0_full[:, factor_ordering]
			self.S_U_full = self.S_U_full[:, factor_ordering]

		self.V_mu = self.V_mu[factor_ordering, :]
		self.V_var = self.V_var[factor_ordering, :]
		self.theta_U_a = self.theta_U_a[factor_ordering]
		self.theta_U_b = self.theta_U_b[factor_ordering]
		self.K = len(factor_ordering)
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
		# Precompute quantities
		U_S_expected_val = self.U_mu*self.S_U
		U_S_squared_expected_val = (np.square(self.U_mu) + self.U_var)*self.S_U
		tau_expected_val = self.tau_alpha/self.tau_beta
		V_mu_copy = np.copy(self.V_mu)
		V_var_copy = np.copy(self.V_var)

		# Keep track of variables
		V_update_data = []

		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				V_update_data.append(outside_update_V_t(V_mu_copy[:, test_index], V_var_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.K, U_S_expected_val, U_S_squared_expected_val, self.F_mu[test_index], self.intercept_mu[test_index], self.gamma_v, tau_expected_val[test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		elif self.parrallel_boolean == True:
			V_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_V_t)(V_mu_copy[:, test_index], V_var_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.K, U_S_expected_val, U_S_squared_expected_val, self.F_mu[test_index], self.intercept_mu[test_index], self.gamma_v, tau_expected_val[test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))

		# Convert to array
		V_update_data = np.asarray(V_update_data).T
		# Fill in data structures
		self.V_mu = V_update_data[(self.K*0):(1*self.K), :]
		self.V_var = V_update_data[(self.K*1):(2*self.K), :]
		'''
		for test_index in range(self.T):
			self.update_V_t(test_index, U_S_expected_val, U_S_squared_expected_val)
		pdb.set_trace()
		'''
	'''
	def update_V_t(self, test_index, U_S_expected_val, U_S_squared_expected_val):
		for k in range(self.K):
			self.update_V_kt(k, test_index, U_S_expected_val, U_S_squared_expected_val[:,k])
	def update_V_kt(self, k, test_index, U_S_expected_val, U_k_S_k_squared_expected_val):
		# Compute relevent expectations
		gamma_k_v_expected_val = self.gamma_v
		tau_t_expected_val = self.tau_alpha[test_index]/self.tau_beta[test_index]
		#U_S_expected_val = self.U_mu*self.S_U
		# U_k_S_k_expected_val = U_S_expected_val[:,k]
		# U_k_S_k_squared_expected_val = (np.square(self.U_mu[:,k]) + self.U_var[:,k])*self.S_U[:,k]
		V_S_t_expected_val = self.V_mu[:,test_index]
		F_S_t_expected_val = self.F_mu[test_index]
		# Compute expectations on other components
		other_components_expected = (U_S_expected_val@V_S_t_expected_val) - U_S_expected_val[:, k]*V_S_t_expected_val[k]
		# Update variance of q(V|s=1)
		a_term = gamma_k_v_expected_val + (1.0/self.sample_batch_fraction)*(tau_t_expected_val*np.sum(np.square(self.G[:,test_index])*U_k_S_k_squared_expected_val))
		#self.V_var[k, test_index] = 1.0/a_term
		# Update variance of q(V|s=1)
		#self.V_var_s_0[k,test_index] = self.gamma_V_beta[k]/self.gamma_V_alpha[k]
		# Update mean of q(U|s=1)
		resid = self.Y[:,test_index] - self.intercept_mu[test_index] - self.G[:,test_index]*(other_components_expected + F_S_t_expected_val)

		b_term = (1.0/self.sample_batch_fraction)*np.sum(tau_t_expected_val*self.G[:,test_index]*U_S_expected_val[:,k]*resid)

		new_var = 1.0/a_term
		new_mu = new_var*b_term
		if self.SVI == False:
			self.V_var[k, test_index] = new_var
			self.V_mu[k, test_index] = new_mu
		elif self.SVI == True:
			self.V_var[k,test_index] = weighted_SVI_updated(self.V_var[k,test_index], new_var, self.step_size)
			self.V_mu[k, test_index] = weighted_SVI_updated(self.V_mu[k,test_index], new_mu, self.step_size)
	'''
	def update_U(self):
		if self.SVI == True:
			# Randomly generate indices for SVI
			svi_sample_indices = self.get_svi_sample_indices()
			# Subset matrices
			self.G = np.copy(self.G_full[svi_sample_indices, :])
			self.Y = np.copy(self.Y_full[svi_sample_indices, :])
			self.z = np.copy(self.z_full[svi_sample_indices])
			self.U_mu = np.copy(self.U_mu_full[svi_sample_indices, :])
			self.S_U = np.copy(self.S_U_full[svi_sample_indices, :])
			self.U_var = np.copy(self.U_var_full[svi_sample_indices, :])
			self.U_var_s_0 = np.copy(self.U_var_s_0_full[svi_sample_indices, :])
		#G_squared = np.square(self.G)
		###################
		# UPDATE U
		###################
		V_S_squared_expected_val = (np.square(self.V_mu) + self.V_var)
		U_mu_copy = np.copy(self.U_mu)
		S_U_copy = np.copy(self.S_U)
		U_var_copy = np.copy(self.U_var)
		U_var_s_0_copy = np.copy(self.U_var_s_0)
		U_update_data = []

		# Don't parrallelize
		if self.parrallel_boolean == False:
			for sample_index in range(self.N):
				U_update_data.append(outside_update_U_n(U_mu_copy[sample_index,:], S_U_copy[sample_index,:], U_var_copy[sample_index,:], U_var_s_0_copy[sample_index,:], self.G[sample_index, :], self.Y[sample_index, :], self.K, self.V_mu, V_S_squared_expected_val, self.F_mu, self.intercept_mu, self.gamma_v, self.tau_alpha/self.tau_beta, self.theta_U_a, self.theta_U_b))
		# Parrallelize
		elif self.parrallel_boolean == True:
			U_update_data = Parallel(n_jobs=self.num_sample_cores)(delayed(outside_update_U_n)(U_mu_copy[sample_index,:], S_U_copy[sample_index,:], U_var_copy[sample_index,:], U_var_s_0_copy[sample_index,:], self.G[sample_index, :], self.Y[sample_index, :], self.K, self.V_mu, V_S_squared_expected_val, self.F_mu, self.intercept_mu, self.gamma_v, self.tau_alpha/self.tau_beta, self.theta_U_a, self.theta_U_b) for sample_index in range(self.N))

		# Convert to array
		U_update_data = np.asarray(U_update_data)
		# Fill in data structures
		self.U_mu = U_update_data[:,(self.K*0):(1*self.K)]
		self.S_U = U_update_data[:,(self.K*1):(2*self.K)]
		self.U_var = U_update_data[:,(self.K*2):(3*self.K)]
		self.U_var_s_0 = U_update_data[:,(self.K*3):(4*self.K)]
		'''
		for sample_index in range(self.N):
			self.update_U_n(sample_index, V_S_squared_expected_val)
		'''
		if self.SVI == True:
			self.U_mu_full[svi_sample_indices, :] = np.copy(self.U_mu)
			self.S_U_full[svi_sample_indices, :] = np.copy(self.S_U)
			self.U_var_full[svi_sample_indices, :] = np.copy(self.U_var)
			self.U_var_s_0_full[svi_sample_indices, :] = np.copy(self.U_var_s_0)
	'''
	def update_U_n(self, sample_index, V_S_squared_expected_val):
		for k in range(self.K):
			self.update_U_nk(sample_index, k, V_S_squared_expected_val[k,:])
	def update_U_nk(self, sample_index, k, V_k_S_k_squared_expected_val):
		# Compute relevent expectations
		gamma_k_u_expected_val = self.gamma_v
		tau_expected_val = self.tau_alpha/self.tau_beta
		V_S_expected_val = self.V_mu
		U_S_expected_val = self.U_mu[sample_index,:]*self.S_U[sample_index,:]
		V_k_S_k_expected_val = V_S_expected_val[k,:]
		#V_k_S_k_squared_expected_val = (np.square(self.V_mu[k,:]) + self.V_var[k,:])
		F_S_expected_val = self.F_mu
		theta_U_expected_val = self.theta_U_a[k]/(self.theta_U_a[k] + self.theta_U_b[k])
		ln_theta_U_expected_val = special.digamma(self.theta_U_a[k]) - special.digamma(self.theta_U_a[k]+self.theta_U_b[k])  # expectation of ln(1-X)
		ln_1_minus_theta_U_expected_val = special.digamma(self.theta_U_b[k]) - special.digamma(self.theta_U_a[k]+self.theta_U_b[k])
		# Compute expectations on other components
		other_components_expected = (U_S_expected_val@V_S_expected_val) - U_S_expected_val[k]*V_S_expected_val[k,:]
		# Update variance of q(U|s=1)
		a_term = np.sum(tau_expected_val*np.square(self.G[sample_index,:])*V_k_S_k_squared_expected_val) + gamma_k_u_expected_val
		self.U_var[sample_index, k] = 1.0/a_term
		# Update variance of q(U|s=0)
		self.U_var_s_0[sample_index, k] = 1.0/self.gamma_v
		# Update mean of q(U|s=1)
		resid = self.Y[sample_index,:] - self.intercept_mu - self.G[sample_index,:]*(F_S_expected_val + other_components_expected)
		b_term = np.sum(tau_expected_val*self.G[sample_index,:]*V_k_S_k_expected_val*resid)
		self.U_mu[sample_index, k] = self.U_var[sample_index, k]*b_term
		# Now update q(S_U=1)
		z_term = ln_theta_U_expected_val - ln_1_minus_theta_U_expected_val + .5*np.log(gamma_k_u_expected_val) - .5*np.log(a_term) + (np.square(b_term)/(2.0*a_term))
		self.S_U[sample_index, k] = sigmoid_function(z_term)
	'''
	'''
	def update_U_k(self, k, V_S_expected_val, V_k_S_k_squared_expected_val, G_squared):
		gamma_u_expected_val = self.gamma_v
		tau_expected_val = self.tau_alpha/self.tau_beta
		U_S_expected_val = self.U_mu*self.S_U
		theta_U_expected_val = self.theta_U_a[k]/(self.theta_U_a[k] + self.theta_U_b[k])
		ln_theta_U_expected_val = special.digamma(self.theta_U_a[k]) - special.digamma(self.theta_U_a[k]+self.theta_U_b[k])  # expectation of ln(1-X)
		ln_1_minus_theta_U_expected_val = special.digamma(self.theta_U_b[k]) - special.digamma(self.theta_U_a[k]+self.theta_U_b[k])

		a_term = (tau_expected_val*G_squared)@V_k_S_k_squared_expected_val + gamma_u_expected_val
		self.U_var[:, k] = 1.0/a_term
		self.U_var_s_0[:, k] = 1.0/self.gamma_v

		other_components_expected = np.dot(U_S_expected_val, V_S_expected_val) - np.dot(np.transpose([U_S_expected_val[:,k]]), [V_S_expected_val[k,:]])
		resid = self.Y - self.intercept_mu - ((self.G)*(other_components_expected + self.F_mu))
		b_term = np.sum(((resid*self.G)*tau_expected_val)*V_S_expected_val[k,:],axis=1)
		self.U_mu[:, k] = self.U_var[:, k]*b_term

		z_term = ln_theta_U_expected_val - ln_1_minus_theta_U_expected_val + .5*np.log(gamma_u_expected_val) - .5*np.log(a_term) + (np.square(b_term)/(2.0*a_term))
		self.S_U[:, k] = sigmoid_function(z_term)
	'''
	def update_F(self):
		U_S_expected_val = self.U_mu*self.S_U

		tau_expected_val = self.tau_alpha/self.tau_beta
		F_mu_copy = np.copy(self.F_mu)
		F_var_copy = np.copy(self.F_var)

		F_update_data = []
		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				F_update_data.append(outside_update_F_t(F_mu_copy[test_index], F_var_copy[test_index], self.G[:, test_index], self.Y[:, test_index], U_S_expected_val, self.V_mu[:,test_index], self.intercept_mu[test_index], self.gamma_v, tau_expected_val[test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		if self.parrallel_boolean == True:
			F_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_F_t)(F_mu_copy[test_index], F_var_copy[test_index], self.G[:, test_index], self.Y[:, test_index], U_S_expected_val, self.V_mu[:,test_index], self.intercept_mu[test_index], self.gamma_v, tau_expected_val[test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))
		F_update_data = np.asarray(F_update_data)
		self.F_mu = F_update_data[:,0]
		self.F_var = F_update_data[:,1]
		'''
		for test_index in range(self.T):
			self.update_F_t(test_index, U_S_expected_val)
		'''
	'''
	def update_F_t(self, test_index, U_S_expected_val):
		# Compute relevent expectations
		gamma_f_expected_val = 1.0
		tau_t_expected_val = self.tau_alpha[test_index]/self.tau_beta[test_index]
		#theta_F_expected_val = self.theta_F_a/(self.theta_F_a + self.theta_F_b)
		#ln_theta_F_expected_val = special.digamma(self.theta_F_a) - special.digamma(self.theta_F_a+self.theta_F_b)  # expectation of ln(1-X)
		#ln_1_minus_theta_F_expected_val = special.digamma(self.theta_F_b) - special.digamma(self.theta_F_a+self.theta_F_b)
		#U_S_expected_val = self.U_mu*self.S_U
		V_S_t_expected_val = self.V_mu[:,test_index]

		# Compute expectations on other components

		other_components_expected = U_S_expected_val@V_S_t_expected_val

		# Update variance of q(F|s=1)
		a_term = gamma_f_expected_val + (1.0/self.sample_batch_fraction)*tau_t_expected_val*np.sum(np.square(self.G[:,test_index]))
		#self.F_var[test_index] = 1.0/a_term
		# Update variance of q(F|s=1)
		#self.F_var_s_0[test_index] = self.gamma_F_beta/self.gamma_F_alpha
		# Update mean of q(F|s=1)
		resid = self.Y[:,test_index] - self.intercept_mu[test_index] - self.G[:,test_index]*(other_components_expected)
		b_term = np.sum(tau_t_expected_val*self.G[:,test_index]*resid*(1.0/self.sample_batch_fraction))
		#self.F_mu[test_index] = self.F_var[test_index]*b_term
		new_var = 1.0/a_term
		new_mu = new_var*b_term
		if self.SVI == False:
			self.F_var[test_index] = new_var
			self.F_mu[test_index] = new_mu
		elif self.SVI == True:
			self.F_var[test_index] = weighted_SVI_updated(self.F_var[test_index], new_var, self.step_size)
			self.F_mu[test_index] = weighted_SVI_updated(self.F_mu[test_index], new_mu, self.step_size)
	'''
	def update_intercept(self):
		U_S_expected_val = self.U_mu*self.S_U
		tau_expected_val = self.tau_alpha/self.tau_beta
		intercept_mu_copy = np.copy(self.intercept_mu)
		intercept_var_copy = np.copy(self.intercept_var)

		intercept_update_data = []
		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				intercept_update_data.append(outside_update_intercept_t(intercept_mu_copy[test_index], intercept_var_copy[test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S_expected_val, self.V_mu[:,test_index], self.F_mu[test_index], tau_expected_val[test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		if self.parrallel_boolean == True:
				intercept_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_intercept_t)(intercept_mu_copy[test_index], intercept_var_copy[test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S_expected_val, self.V_mu[:,test_index], self.F_mu[test_index], tau_expected_val[test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))

		intercept_update_data = np.asarray(intercept_update_data)
		self.intercept_mu = intercept_update_data[:,0]
		self.intercept_var = intercept_update_data[:,1]
		'''
		for test_index in range(self.T):
			self.update_intercept_t(test_index, U_S_expected_val)
		'''
	'''
	def update_intercept_t(self, test_index, U_S_expected_val):
		# Compute relevent expectations
		tau_t_expected_val = self.tau_alpha[test_index]/self.tau_beta[test_index]
		F_S_t_expected_val = self.F_mu[test_index]
		#U_S_expected_val = self.U_mu*self.S_U
		V_S_t_expected_val = self.V_mu[:,test_index]

		# Compute expectations on other components
		other_components_expected = U_S_expected_val@V_S_t_expected_val
		resid = self.Y[:, test_index] - self.G[:, test_index]*(F_S_t_expected_val + other_components_expected)

		new_var = 1.0/((1.0/self.sample_batch_fraction)*self.N*tau_t_expected_val)
		new_mu = new_var*tau_t_expected_val*np.sum(resid)*(1.0/self.sample_batch_fraction)

		if self.SVI == False:
			self.intercept_var[test_index] = new_var
			self.intercept_mu[test_index] = new_mu
		elif self.SVI == True:
			self.intercept_var[test_index] = weighted_SVI_updated(self.intercept_var[test_index], new_var, self.step_size)
			self.intercept_mu[test_index] = weighted_SVI_updated(self.intercept_mu[test_index], new_mu, self.step_size)
	'''
	def update_gamma_U(self):
		# Loop through factors
		for k in range(self.K):
			#U_squared_k_expected_val = ((np.square(self.U_mu[:,k]) + self.U_var[:,k])*self.S_U[:,k]) + (1.0-self.S_U[:,k])*(self.gamma_U_beta[k]/self.gamma_U_alpha[k])
			U_squared_k_expected_val = ((np.square(self.U_mu[:,k]) + self.U_var[:,k])*self.S_U[:,k]) + (1.0-self.S_U[:,k])*(self.U_var_s_0[:,k])
			self.gamma_U_alpha[k] = self.alpha_prior + (self.N/2.0)
			self.gamma_U_beta[k] = self.beta_prior + np.sum(U_squared_k_expected_val)/2.0
	def update_theta_U(self):
		# Loop through factors
		for k in range(self.K):
			#self.theta_U_a[k] = self.a_prior + np.sum(self.S_U[:,k])
			#self.theta_U_b[k] = self.b_prior + self.N - np.sum(self.S_U[:,k])
			new_theta_U_a = self.a_prior + (1.0/self.sample_batch_fraction)*np.sum(self.S_U[:,k])
			new_theta_U_b = self.b_prior + (1.0/self.sample_batch_fraction)*(self.N - np.sum(self.S_U[:,k]))
			if self.SVI == False:
				self.theta_U_a[k] = new_theta_U_a
				self.theta_U_b[k] = new_theta_U_b
			elif self.SVI == True:
				self.theta_U_a[k] = weighted_SVI_updated(self.theta_U_a[k], new_theta_U_a, self.step_size)
				self.theta_U_b[k] = weighted_SVI_updated(self.theta_U_b[k], new_theta_U_b, self.step_size)

	def update_tau(self):
		tau_alpha_copy = np.copy(self.tau_alpha)
		tau_beta_copy = np.copy(self.tau_beta)

		# Precompute quantities
		F_S_squared = np.square(self.F_mu) + self.F_var
		intercept_squared = np.square(self.intercept_mu) + self.intercept_var
		V_S_squared = np.square(self.V_mu) + self.V_var
		U_S_squared = ((np.square(self.U_mu) + self.U_var)*self.S_U)
		U_S = (self.U_mu)*(self.S_U)
		# Loop through tests
		tau_update_data = []
		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				tau_update_data.append(outside_update_tau_t(tau_alpha_copy[test_index], tau_beta_copy[test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S, self.V_mu[:,test_index], self.F_mu[test_index], self.intercept_mu[test_index], V_S_squared[:, test_index], F_S_squared[test_index], U_S_squared, intercept_squared[test_index], self.alpha_prior, self.beta_prior, self.sample_batch_fraction, self.step_size, self.SVI))
		if self.parrallel_boolean == True:
			tau_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_tau_t)(tau_alpha_copy[test_index], tau_beta_copy[test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S, self.V_mu[:,test_index], self.F_mu[test_index], self.intercept_mu[test_index], V_S_squared[:, test_index], F_S_squared[test_index], U_S_squared, intercept_squared[test_index], self.alpha_prior, self.beta_prior, self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))
		tau_update_data = np.asarray(tau_update_data)
		self.tau_alpha = tau_update_data[:,0]
		self.tau_beta = tau_update_data[:,1]
		'''
		for test_index in range(self.T):
			self.update_tau_t(test_index, F_S_squared[test_index], intercept_squared[test_index], V_S_squared[:, test_index], U_S_squared, U_S)
		'''
	'''
	def update_tau_t(self, test_index, F_S_t_squared, intercept_t_squared, V_S_t_squared, U_S_squared, U_S):
		# Compute Relevent expectations
		#F_S_t_squared = ((np.square(self.F_mu[test_index]) + self.F_var[test_index]))
		F_S_t = self.F_mu[test_index]
		#intercept_t_squared = ((np.square(self.intercept_mu[test_index]) + self.intercept_var[test_index]))
		intercept_t = self.intercept_mu[test_index]
		#V_S_t_squared = ((np.square(self.V_mu[:,test_index]) + self.V_var[:,test_index]))
		V_S_t = self.V_mu[:,test_index]
		#U_S_squared = ((np.square(self.U_mu) + self.U_var)*self.S_U)
		#U_S = (self.U_mu)*(self.S_U)

		squared_factor_terms = U_S_squared@V_S_t_squared
		factor_terms = U_S@V_S_t

		# First add together square terms
		resid = np.square(self.Y[:,test_index]) + intercept_t_squared + np.square(self.G[:,test_index])*(F_S_t_squared + squared_factor_terms)
		# Now add terms with Y
		resid = resid - (2.0*self.Y[:,test_index]*(intercept_t + self.G[:,test_index]*factor_terms + self.G[:,test_index]*F_S_t))

		resid = resid + 2.0*intercept_t*(self.G[:,test_index]*(factor_terms + F_S_t))
		# Now add terms with factors
		resid = resid + 2.0*self.G[:,test_index]*factor_terms*self.G[:,test_index]*F_S_t
		# Now add terms with interactions between factors
		resid = resid + (np.square(self.G[:,test_index])*(factor_terms*factor_terms - np.sum(np.square(U_S*V_S_t),axis=1)))

		# Make Updates
		new_alpha = self.alpha_prior + ((self.N/2.0)*(1.0/self.sample_batch_fraction))
		new_beta = self.beta_prior + ((np.sum(resid)/2.0)*(1.0/self.sample_batch_fraction))

		if self.SVI == False:
			self.tau_alpha[test_index] = new_alpha
			self.tau_beta[test_index] = new_beta
		elif self.SVI == True:
			self.tau_alpha[test_index] = weighted_SVI_updated(self.tau_alpha[test_index], new_alpha, self.step_size)
			self.tau_beta[test_index] = weighted_SVI_updated(self.tau_beta[test_index], new_beta, self.step_size)
	'''
	def update_elbo(self):
		data_likelihood_term = self.compute_elbo_log_likelihood_term()
		kl_V_S = self.compute_kl_divergence_of_V_S()
		kl_U_S = self.compute_kl_divergence_of_U_S()
		kl_F_S = self.compute_kl_divergence_of_F_S()
		# kl_gamma_u = self.compute_kl_divergence_of_gamma_u()
		kl_tau = self.compute_kl_divergence_of_tau()
		#kl_theta_v = self.compute_kl_divergence_of_theta_v()
		kl_theta_u = self.compute_kl_divergence_of_theta_u()
		#kl_theta_f = self.compute_kl_divergence_of_theta_f()

		kl_divergence = kl_V_S + kl_U_S + kl_F_S + kl_tau + kl_theta_u 

		elbo = data_likelihood_term - kl_divergence
		self.elbo.append(elbo)

	def compute_kl_divergence_of_theta_u(self):
		a_prior = self.a_prior
		b_prior = self.b_prior
		theta_a = self.theta_U_a 
		theta_b = self.theta_U_b
		kl_divergence = compute_kl_divergence_of_beta(a_prior, b_prior, theta_a, theta_b)
		return kl_divergence
	def compute_kl_divergence_of_tau(self):
		alpha_prior = self.alpha_prior
		beta_prior = self.beta_prior
		gamma_alpha = self.tau_alpha
		gamma_beta = self.tau_beta
		kl_divergence = compute_kl_divergence_of_gamma(alpha_prior, beta_prior, gamma_alpha, gamma_beta)
		return kl_divergence
	def compute_kl_divergence_of_gamma_u(self):
		alpha_prior = self.alpha_prior
		beta_prior = self.beta_prior
		gamma_alpha = self.gamma_U_alpha
		gamma_beta = self.gamma_U_beta
		kl_divergence = compute_kl_divergence_of_gamma(alpha_prior, beta_prior, gamma_alpha, gamma_beta)
		return kl_divergence
	def compute_kl_divergence_of_F_S(self):
		W_mu = np.asarray([self.F_mu])
		W_var = np.asarray([self.F_var])
		expected_gamma = 1.0
		kl_divergence = compute_kl_divergence_of_gaussian_fixed_variance(W_mu, W_var, expected_gamma, 1)
		return kl_divergence
	def compute_kl_divergence_of_V_S(self):
		W_mu = self.V_mu
		W_var = self.V_var
		expected_gamma_v = self.gamma_v
		kl_divergence = compute_kl_divergence_of_gaussian_fixed_variance(W_mu, W_var, expected_gamma_v, self.K)
		return kl_divergence
	def compute_kl_divergence_of_U_S(self):
		if self.SVI == True:
			S = np.transpose(self.S_U_full)
			W_mu = np.transpose(self.U_mu_full)
			W_var = np.transpose(self.U_var_full)
			W_var_s_0 = np.transpose(self.U_var_s_0_full)
		elif self.SVI == False:
			S = np.transpose(self.S_U)
			W_mu = np.transpose(self.U_mu)
			W_var = np.transpose(self.U_var)
			W_var_s_0 = np.transpose(self.U_var_s_0)
		gamma_expected = self.gamma_v
		theta_a = self.theta_U_a
		theta_b = self.theta_U_b
		kl_divergence = compute_kl_divergence_of_gaussian_bernoulli(S, W_mu, W_var, W_var_s_0, gamma_expected, theta_a, theta_b, self.K)
		return kl_divergence
	def compute_elbo_log_likelihood_term(self):
		# Compute expectation of log of gamma variables
		log_tau_expected = special.digamma(self.tau_alpha) - np.log(self.tau_beta)
		# Compute expectation of gamma variable
		tau_expected = self.tau_alpha/self.tau_beta
		# Other relevent expectations
		if self.SVI == True:
			U_S = (self.U_mu_full)*(self.S_U_full)
		elif self.SVI == False:
			U_S = (self.U_mu)*(self.S_U)
		V_S = (self.V_mu)
		F_S = (self.F_mu)
		# alpha_squared = np.square(self.alpha_big_mu) + self.alpha_big_var
		# alpha = self.alpha_big_mu
		F_S_squared = ((np.square(self.F_mu) + self.F_var))
		V_S_squared = ((np.square(self.V_mu) + self.V_var))
		if self.SVI == True:
			U_S_squared = ((np.square(self.U_mu_full) + self.U_var_full)*self.S_U_full)
		else:
			U_S_squared = ((np.square(self.U_mu) + self.U_var)*self.S_U)
		intercept_squared = np.square(self.intercept_mu) + self.intercept_var
		intercept = self.intercept_mu

		componenent_squared_terms = np.dot(U_S_squared, V_S_squared) + np.dot(np.ones((self.N_full,1)),[F_S_squared])
		componenent_terms = np.dot(U_S, V_S)
		F_terms = np.dot(np.ones((self.N_full,1)),[F_S])
		intercept_terms = np.dot(np.ones((self.N_full,1)),[intercept])
		intercept_squared_terms = np.dot(np.ones((self.N_full,1)),[intercept_squared])


		# Terms of interest in likelihood
		term_a = -np.log(2.0*np.pi)*(self.N_full*self.T_full/2.0)
		term_b = (self.N_full/2.0)*np.sum(log_tau_expected)
		# Compute residual matrix
		#residual_mat = self.Y - self.G*(np.dot(U_S, V_S) + np.dot(np.ones((self.N,1)),[F_S])) - self.alpha_big_mu
		squared_residual_mat = np.square(self.Y_full) + intercept_squared_terms + np.square(self.G_full)*componenent_squared_terms
		squared_residual_mat = squared_residual_mat - 2.0*self.Y_full*(intercept_terms + self.G_full*(componenent_terms+ F_terms))
		squared_residual_mat = squared_residual_mat + 2.0*intercept_terms*(self.G_full*(componenent_terms + F_terms))
		squared_residual_mat = squared_residual_mat + 2.0*np.square(self.G_full)*componenent_terms*F_terms
		#squared_residual_mat = squared_residual_mat + 2.0*np.square(self.G)*componenent_terms*componenent_terms
		squared_residual_mat = squared_residual_mat + np.square(self.G_full)*(componenent_terms*componenent_terms)
		for k in range(self.K):
			squared_residual_mat = squared_residual_mat - np.square(self.G_full)*np.square(np.dot(np.transpose([U_S[:,k]]), [V_S[k,:]]))

		term_c = np.sum(squared_residual_mat*tau_expected)/2.0
		data_likelihood_term = term_a + term_b - term_c
		return data_likelihood_term
	# Randomly generate indices
	def get_svi_sample_indices(self):
		return np.random.choice(self.N_full, size=self.N, replace=False)
	def initialize_variables(self):
		# Initialize array to keep track of ELBO
		self.elbo = []

		# Add model dimensions to object
		self.N_full = self.Y_full.shape[0]
		self.T_full = self.Y_full.shape[1]

		# Do standard variational inference
		if self.SVI == False:
			self.N = self.Y_full.shape[0]
			self.T = self.Y_full.shape[1]
			self.G = np.copy(self.G_full)
			self.Y = np.copy(self.Y_full)
			self.z = np.copy(self.z_full)
			pca = sklearn.decomposition.PCA(n_components=self.K, whiten=True)
			pca.fit(np.random.randn(self.N, 9999).T)
			self.U_mu = pca.components_.T*10.0
			#self.U_mu = np.random.randn(self.N, self.K)
			self.U_var = np.ones((self.N, self.K))*(1.0/self.gamma_v) 
			self.U_var_s_0 = np.ones((self.N, self.K))*(1.0/self.gamma_v)
			self.S_U = np.ones((self.N,self.K))
		# Do stochastic variational inference
		elif self.SVI == True:
			self.T = self.Y_full.shape[1]
			#self.N = round((self.Y_full.shape[0])*self.sample_batch_fraction)
			self.N = self.Y_full.shape[0]
			# Randomly generate indices for SVI
			#svi_sample_indices = self.get_svi_sample_indices()
			# For initialization we don't randomly subset
			svi_sample_indices = np.arange(self.N_full)
			actual_sample_batch_fraction = self.sample_batch_fraction
			self.sample_batch_fraction = 1.0
			# Subset matrices
			self.G = np.copy(self.G_full)[svi_sample_indices, :]
			self.Y = np.copy(self.Y_full)[svi_sample_indices, :]
			self.z = np.copy(self.z_full)[svi_sample_indices]
			pca = sklearn.decomposition.PCA(n_components=self.K, whiten=True)
			pca.fit(np.random.randn(self.N_full, 9999).T)
			self.U_mu_full = pca.components_.T*10.0
			#self.U_mu_full = np.random.randn(self.N_full, self.K)
			self.U_var_full = np.ones((self.N_full, self.K))*(1.0/self.gamma_v)
			self.U_var_s_0_full = np.ones((self.N_full, self.K))*(1.0/self.gamma_v)
			self.S_U_full = np.ones((self.N_full,self.K))
			self.U_mu = np.copy(self.U_mu_full[svi_sample_indices, :])
			self.U_var = np.copy(self.U_var_full[svi_sample_indices, :])
			self.U_var_s_0 = np.copy(self.U_var_s_0_full[svi_sample_indices, :])
			self.S_U = np.copy(self.S_U_full[svi_sample_indices, :])

		# Fit model with limited priors for initialization
		#eqtl_vi_init = eqtl_factorization_vi_shared_effect_no_priors_for_initialization.EQTL_FACTORIZATION_VI(K=self.K, alpha=self.alpha_prior, beta=self.beta_prior, gamma_v=self.gamma_v, max_iter=3, delta_elbo_threshold=.01)
		#eqtl_vi_init.fit(G=self.G, Y=self.Y, z=self.z)

		#self.U_mu = eqtl_vi_init.U_mu
		#self.U_mu = np.random.randn(self.N, self.K)
		# Initialize with orthogonal factors

		
		#self.V_mu = eqtl_vi_init.V_mu
		#self.V_mu = np.random.randn(self.K, self.T)
		pca = sklearn.decomposition.PCA(n_components=self.K, whiten=True)
		pca.fit(np.random.randn(self.T, 9999).T)
		self.V_mu = pca.components_*10.0
		self.V_var = np.ones((self.K, self.T))*(1.0/self.gamma_v)

		betas = run_linear_model_for_initialization(self.Y, self.G)
		#self.F_mu = eqtl_vi_init.F_mu
		self.F_mu = betas
		self.F_var = np.ones(self.T)
		# Masks
		#self.S_U = np.random.uniform(0, 1, size=(self.N, self.K))
		# Intercepts
		# self.intercept_mu = eqtl_vi_init.intercept_mu
		self.intercept_mu = np.zeros(self.T)
		self.intercept_var = np.ones(self.T)

		# Variances
		self.tau_alpha = np.ones(self.T)*self.alpha_prior
		self.tau_beta = np.ones(self.T)*self.beta_prior
		# bernoulli probs
		self.theta_U_a = np.ones(self.K)*self.a_prior + 9
		self.theta_U_b = np.ones(self.K)*self.b_prior

		# Initialize other variables based around U
		self.step_size = 1.0
		self.update_V()
		self.update_intercept()
		self.update_F()
		self.update_tau()

		if self.SVI == True:
			self.sample_batch_fraction = actual_sample_batch_fraction
			self.N = round((self.Y_full.shape[0])*self.sample_batch_fraction)

		#self.gamma_U_alpha = np.ones(self.K)*self.alpha_prior
		#self.gamma_U_beta = np.ones(self.K)*self.beta_prior

