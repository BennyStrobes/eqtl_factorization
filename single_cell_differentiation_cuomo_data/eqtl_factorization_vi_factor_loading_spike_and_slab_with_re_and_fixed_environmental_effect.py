import numpy as np 
import os
import sys
import pdb
import scipy.special as special
from sklearn.linear_model import LinearRegression
import time
import sklearn.decomposition
from pymer4.models import Lmer
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing
#import time

def sigmoid_function(x):
	return 1.0/(1.0 + np.exp(-x))

def weighted_SVI_updated(old_parameter, new_parameter, step_size):
	updated_parameter = ((1.0 - step_size)*old_parameter) + (step_size*new_parameter)
	return updated_parameter

def run_linear_model_for_initialization(Y, G, cov, U):
	num_tests = Y.shape[1]
	betas = []
	c_vec = []
	v_vec = []
	delta_vec = []
	intercepts = []
	num_cov = cov.shape[1]
	for test_number in range(num_tests):
		y_vec = Y[:,test_number]
		g_vec = G[:,test_number]
		K = U.shape[1]
		U_scaled = np.copy(U)
		for k in range(K):
			U_scaled[:, k] = U[:,k]*g_vec
		model_covariates = np.hstack((np.transpose(np.asmatrix(g_vec)), cov, U_scaled, U))
		reg = LinearRegression().fit(model_covariates, np.transpose(np.asmatrix(y_vec)))
		betas.append(reg.coef_[0][0])
		c_vec.append(reg.coef_[0][1:(1+num_cov)])
		v_vec.append(reg.coef_[0][(1+num_cov):(1+num_cov+U.shape[1])])
		delta_vec.append(reg.coef_[0][(1+num_cov+U.shape[1]):])
		intercepts.append(reg.intercept_[0])
	return np.asarray(v_vec), np.asarray(betas), np.asarray(c_vec), np.asarray(intercepts), np.asarray(delta_vec)

def run_linear_mixed_model_one_test(y_vec, g_vec, U, cov, z):
	num_cov = cov.shape[1]
	num_latent_factors = U.shape[1]
	U_scaled = np.copy(U)
	for k in range(num_latent_factors):
		U_scaled[:, k] = U[:,k]*g_vec
	# Covariate matrix
	X = np.vstack((y_vec, z, g_vec, cov.T, U_scaled.T)).T
	# Create column names
	cov_names = ['cov'+str(i) for i in range(num_cov)]
	latent_factor_names = ['lf'+str(i) for i in range(num_latent_factors)]
	col_names = ['y', 'group', 'g'] + cov_names + latent_factor_names

	# Make df
	df = pd.DataFrame(X, columns=col_names)
	# Make formula for LMM
	formula = 'y ~ g + ' + ' + '.join(cov_names) + ' + ' + ' + '.join(latent_factor_names) + ' + (1 | group)'

	model = Lmer(formula, data=df)
	model.fit()
	return np.asarray(model.coefs['Estimate'])


def run_linear_mixed_model_for_initialization(Y, G, cov, U, z):
	num_tests = Y.shape[1]
	betas = []
	c_vec = []
	v_vec = []
	intercepts = []
	num_cov = cov.shape[1]
	num_latent_factors = U.shape[1]
	'''
	results = []
	for test_number in range(num_tests):
		results.append(run_linear_mixed_model_one_test(Y[:,test_number], G[:,test_number], U, cov, z))
	results = np.asarray(results)
	'''

	#U_update_data = Parallel(n_jobs=self.num_sample_cores)(delayed(outside_update_U_n)(U_mu_copy[sample_index,:],S_U_copy[sample_index,:], U_var_copy[sample_index,:], U_var_s_0_copy[sample_index,:], self.G[sample_index, :], self.Y[sample_index, :], self.K, V_S_expected_val, V_S_squared_expected_val, self.F_mu, self.S[sample_index,:], self.intercept_mu, self.gamma_u, self.tau_alpha/self.tau_beta, self.cov[sample_index,:], self.C_mu, self.alpha_big_mu[sample_index,:], self.theta_U_a, self.theta_U_b) for sample_index in range(self.N))
	
	results = Parallel(n_jobs=10)(delayed(run_linear_mixed_model_one_test)(Y[:,test_number], G[:,test_number], U, cov, z) for test_number in range(num_tests))
	results = np.asarray(results)
	intercepts = results[:, 0]
	betas = results[:, 1]
	c_vec = results[:,2:(2+num_cov)]
	v_vec = results[:, (2+num_cov):]

			#intercepts.append(model.coefs['Estimate'][0])
		#betas.append(model.coefs['Estimate'][1])
		#c_vec.append(model.coefs['Estimate'][2:(2+num_cov)])
		#v_vec.append(model.coefs['Estimate'][(2+num_cov):])
		#print(test_number)
		#y_vec = Y[:,test_number]
		#g_vec = G[:,test_number]
		#3K = U.shape[1]
		#U_scaled = np.copy(U)
		#for k in range(K):
		#	U_scaled[:, k] = U[:,k]*g_vec
		# Covariate matrix
		#X = np.vstack((y_vec, z, g_vec, cov.T, U_scaled.T)).T
		# Create column names
		#cov_names = ['cov'+str(i) for i in range(num_cov)]
		#latent_factor_names = ['lf'+str(i) for i in range(num_latent_factors)]
		#col_names = ['y', 'group', 'g'] + cov_names + latent_factor_names

		# Make df
		#df = pd.DataFrame(X, columns=col_names)
		# Make formula for LMM
		#formula = 'y ~ g + ' + ' + '.join(cov_names) + ' + ' + ' + '.join(latent_factor_names) + ' + (1 | group)'


		#model = Lmer(formula, data=df)
		#model.fit()


		#intercepts.append(model.coefs['Estimate'][0])
		#betas.append(model.coefs['Estimate'][1])
		#c_vec.append(model.coefs['Estimate'][2:(2+num_cov)])
		#v_vec.append(model.coefs['Estimate'][(2+num_cov):])
	
		#beta = model.coefs['Estimate'][1]
		#standard_error = model.coefs['SE'][1]
		#pvalue = model.coefs['P-val'][1]
		#reg = LinearRegression().fit(model_covariates, np.transpose(np.asmatrix(y_vec)))
		#betas.append(reg.coef_[0][0])
		#c_vec.append(reg.coef_[0][1:(1+num_cov)])
		#v_vec.append(reg.coef_[0][(1+num_cov):])
		#intercepts.append(reg.intercept_[0])
	return np.asarray(v_vec), np.asarray(betas), np.asarray(c_vec), np.asarray(intercepts)


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



def outside_update_U_n(U_mu, S_U, U_var,U_var_s_0, G_slice, Y_slice, K, V_S_expected_val, V_S_squared_expected_val, F_S_expected_val, S_slice, intercept_mu, gamma_u, tau_expected_val, cov_slice, C_mu, alpha_i_expected_val, theta_U_a, theta_U_b, Delta_expected_val, Delta_squared_expected_val):
	covariate_predictions = cov_slice@C_mu
	for k in range(K):
		# Compute relevent expectations
		U_S_expected_val = U_mu
		V_k_S_k_expected_val = V_S_expected_val[k,:]
		theta_U_expected_val = theta_U_a[k]/(theta_U_a[k] + theta_U_b[k])
		ln_theta_U_expected_val = special.digamma(theta_U_a[k]) - special.digamma(theta_U_a[k]+theta_U_b[k])  # expectation of ln(1-X)
		ln_1_minus_theta_U_expected_val = special.digamma(theta_U_b[k]) - special.digamma(theta_U_a[k]+theta_U_b[k])

		# Compute expectations on other components
		other_components_expected = (U_S_expected_val@V_S_expected_val) - U_S_expected_val[k]*V_S_expected_val[k,:]
		other_non_genetic_components_expected = (U_S_expected_val@Delta_expected_val) - U_S_expected_val[k]*Delta_expected_val[k,:]
		# Update variance of q(U|s=1)
		a_term = np.sum(tau_expected_val*S_slice*(np.square(G_slice)*V_S_squared_expected_val[k,:] + Delta_squared_expected_val[k,:] + G_slice*V_S_expected_val[k,:]*Delta_expected_val[k,:])) + gamma_u
		U_var[k] = 1.0/a_term

		U_var_s_0[k] = 1.0/gamma_u
		# Update mean of q(U|s=1)
		#resid = Y_slice - intercept_mu - alpha_i_expected_val - (covariate_predictions) - G_slice*(F_S_expected_val + other_components_expected)
		#b_term = np.sum(tau_expected_val*G_slice*V_k_S_k_expected_val*S_slice*resid)
		resid = Y_slice - intercept_mu - alpha_i_expected_val - (covariate_predictions) - G_slice*(F_S_expected_val)
		temp1 = ((G_slice*V_k_S_k_expected_val + Delta_expected_val[k,:])*resid) - np.square(G_slice)*other_components_expected*V_k_S_k_expected_val - other_non_genetic_components_expected*Delta_expected_val[k,:] - .5*other_non_genetic_components_expected*G_slice*V_k_S_k_expected_val - .5*other_components_expected*G_slice*Delta_expected_val[k,:]
		b_term = np.sum(tau_expected_val*temp1)

		U_mu[k] = U_var[k]*b_term
		# Now update q(S_U=1)
		z_term = ln_theta_U_expected_val - ln_1_minus_theta_U_expected_val + .5*np.log(gamma_u) - .5*np.log(a_term) + (np.square(b_term)/(2.0*a_term))
		S_U[k] = sigmoid_function(z_term)

	return np.hstack((U_mu, S_U, U_var, U_var_s_0))

def outside_update_V_t(V_mu, V_var, S_V, V_var_s_0, G_slice, Y_slice, K, U_S_expected_val, U_S_squared_expected_val, F_S_t_expected_val, intercept_mu, S_slice, gamma_v, tau_t_expected_val, cov, C_mu_t, alpha_t_mu, Delta_t_mu, theta_V_a, theta_V_b, sample_batch_fraction, step_size, SVI):
	covariates_expected = (cov@C_mu_t) + (U_S_expected_val@Delta_t_mu)
	for k in range(K):
		ln_theta_V_expected_val = special.digamma(theta_V_a[k]) - special.digamma(theta_V_a[k]+theta_V_b[k])  # expectation of ln(1-X)
		ln_1_minus_theta_V_expected_val = special.digamma(theta_V_b[k]) - special.digamma(theta_V_a[k]+theta_V_b[k])
		# Compute expectations on other components
		other_components_expected = (U_S_expected_val@V_mu) - U_S_expected_val[:, k]*V_mu[k]
		# Update variance of q(V|s=1)
		a_term = gamma_v + (1.0/sample_batch_fraction)*(tau_t_expected_val*np.sum(S_slice*np.square(G_slice)*U_S_squared_expected_val[:,k]))
		# Update mean of q(U|s=1)
		resid = Y_slice - alpha_t_mu - intercept_mu - covariates_expected - G_slice*(other_components_expected + F_S_t_expected_val)
		b_term = (1.0/sample_batch_fraction)*np.sum(tau_t_expected_val*G_slice*U_S_expected_val[:,k]*S_slice*resid)

		z_term = ln_theta_V_expected_val - ln_1_minus_theta_V_expected_val + .5*np.log(gamma_v) - .5*np.log(a_term) + (np.square(b_term)/(2.0*a_term))
		new_var = 1.0/a_term
		new_mu = new_var*b_term
		if SVI == False:
			V_var[k] = new_var
			V_mu[k] = new_mu
			V_var_s_0[k] = 1.0/gamma_v
			S_V[k] = sigmoid_function(z_term)
		elif SVI == True:
			V_var[k] = weighted_SVI_updated(V_var[k], new_var, step_size)
			V_mu[k] = weighted_SVI_updated(V_mu[k], new_mu, step_size)
			pdb.set_trace()
	return np.hstack((V_mu, V_var, S_V, V_var_s_0))

def outside_update_Delta_t(Delta_mu, Delta_var, G_slice, Y_slice, K, U_S_expected_val, U_S_squared_expected_val, F_S_t_expected_val, intercept_mu, S_slice, gamma_v, tau_t_expected_val, cov, C_mu_t, alpha_t_mu, V_mu_t, sample_batch_fraction, step_size, SVI):
	interaction_components_expected = U_S_expected_val@V_mu_t
	covariates_expected = (cov@C_mu_t)
	for k in range(K):
		other_components_expected = (U_S_expected_val@Delta_mu) - U_S_expected_val[:, k]*Delta_mu[k]
		a_term = gamma_v + (1.0/sample_batch_fraction)*(tau_t_expected_val*np.sum(S_slice*U_S_squared_expected_val[:,k]))
		resid = Y_slice - intercept_mu - alpha_t_mu - other_components_expected - covariates_expected - G_slice*(interaction_components_expected + F_S_t_expected_val)
		b_term = (1.0/sample_batch_fraction)*np.sum(tau_t_expected_val*S_slice*U_S_expected_val[:,k]*resid)
		new_var = 1.0/a_term
		new_mu = new_var*b_term
		if SVI == False:
			Delta_var[k] = new_var
			Delta_mu[k] = new_mu
		elif SVI == True:
			Delta_var[k] = weighted_SVI_updated(Delta_var[k], new_var, step_size)
			Delta_mu[k] = weighted_SVI_updated(Delta_mu[k], new_mu, step_size)
	return np.hstack((Delta_mu, Delta_var))

def outside_update_C_t(C_mu, C_var, G_slice, Y_slice, K, U_S_expected_val, V_S_t_expected_val, F_S_t_expected_val, intercept_mu, S_slice, gamma_v, tau_t_expected_val, cov, alpha_t_mu, delta_t_mu, sample_batch_fraction, step_size, SVI):
	num_cov = cov.shape[1]
	interaction_components_expected = U_S_expected_val@V_S_t_expected_val
	components_expected = U_S_expected_val@delta_t_mu
	for j in range(num_cov):
		# Compute expectations on other components
		other_covariates_expected = (cov@C_mu) - cov[:, j]*C_mu[j]
		# Update variance of q(V|s=1)
		a_term = (1.0/sample_batch_fraction)*(tau_t_expected_val*np.sum(S_slice*np.square(cov[:, j])))

		resid = Y_slice - intercept_mu - alpha_t_mu - other_covariates_expected - components_expected - G_slice*(interaction_components_expected + F_S_t_expected_val)
		b_term = (1.0/sample_batch_fraction)*np.sum(tau_t_expected_val*S_slice*cov[:,j]*resid)
		
		new_var = 1.0/a_term
		new_mu = new_var*b_term
		if SVI == False:
			C_var[j] = new_var
			C_mu[j] = new_mu
		elif SVI == True:
			C_var[j] = weighted_SVI_updated(C_var[j], new_var, step_size)
			C_mu[j] = weighted_SVI_updated(C_mu[j], new_mu, step_size)
	return np.hstack((C_mu, C_var))

def outside_update_intercept_t(intercept_mu, intercept_var, G_slice, Y_slice, N, U_S_expected_val, V_S_t_expected_val, F_S_t_expected_val, tau_t_expected_val, S_slice, cov, C_mu_t, alpha_t_mu, delta_t_mu, sample_batch_fraction, step_size, SVI):
	# Compute relevent expectations
	# Compute expectations on other components
	interaction_components_expected = U_S_expected_val@V_S_t_expected_val
	components_expected = U_S_expected_val@delta_t_mu
	covariates_expected = cov@C_mu_t
	resid = (S_slice*Y_slice) - (S_slice*alpha_t_mu) - (S_slice*G_slice*(F_S_t_expected_val + interaction_components_expected)) - (S_slice*covariates_expected) - (S_slice*components_expected)

	new_var = 1.0/((1.0/sample_batch_fraction)*sum(S_slice)*tau_t_expected_val)
	new_mu = new_var*tau_t_expected_val*np.sum(resid)*(1.0/sample_batch_fraction)

	if SVI == False:
		intercept_var = new_var
		intercept_mu = new_mu
	elif SVI == True:
		intercept_var = weighted_SVI_updated(intercept_var, new_var, step_size)
		intercept_mu = weighted_SVI_updated(intercept_mu, new_mu, step_size)
	return np.hstack((intercept_mu, intercept_var))

def outside_update_alpha_t(alpha_mu_copy, alpha_var_copy, G_slice, Y_slice, I, U_S_expected_val, V_S_t_expected_val, F_S_t_expected_val, intercept_mu_t, cov, C_mu_t, tau_t_expected_val, psi_t_expected_val, Delta_t_expected_val, individual_to_sample_indices, individual_to_number_full_indices, step_size, SVI):
	other_components_expected = U_S_expected_val@V_S_t_expected_val

	covariates_expected = (cov@C_mu_t) + (U_S_expected_val@Delta_t_expected_val)
	resid = Y_slice - intercept_mu_t - covariates_expected - G_slice*(F_S_t_expected_val + other_components_expected)
	# Loop through individuals
	for individual_index in range(I):
		# Indices of samples corresponding to this label
		sample_indices = individual_to_sample_indices[individual_index]
		# Number of indices corresponding to this sample (w/o subsetting)
		num_full_indices = individual_to_number_full_indices[individual_index]
		# Number of indices corresponding to this individaul
		n_i = len(sample_indices)
		individual_batch_fraction = n_i/num_full_indices
		# Update variance of q(alpha_it)
		new_var = 1.0/((1.0/individual_batch_fraction)*n_i*tau_t_expected_val + psi_t_expected_val)
		new_mu = new_var*tau_t_expected_val*np.sum(resid[sample_indices])*(1.0/individual_batch_fraction)
		if SVI == False:
			alpha_var_copy[individual_index] = new_var
			alpha_mu_copy[individual_index] = new_mu
		elif SVI == True:
			alpha_var_copy[individual_index] = weighted_SVI_updated(alpha_var_copy[individual_index], new_var, step_size)
			alpha_mu_copy[individual_index] = weighted_SVI_updated(alpha_mu_copy[individual_index], new_mu, step_size)
	return np.hstack((alpha_mu_copy, alpha_var_copy))

def outside_update_F_t(F_mu, F_var, G_slice, Y_slice, U_S_expected_val, V_S_t_expected_val, intercept_mu_t, gamma_f_expected_val, tau_t_expected_val, S_slice, cov, C_mu_t, alpha_t_mu, delta_t_mu, sample_batch_fraction, step_size, SVI):
	# Compute expectations on other components
	components_expected = U_S_expected_val@V_S_t_expected_val
	covariates_expected = (cov@C_mu_t) + (U_S_expected_val@delta_t_mu)

	# Update variance of q(F|s=1)
	a_term = gamma_f_expected_val + (1.0/sample_batch_fraction)*tau_t_expected_val*np.sum(S_slice*np.square(G_slice))
	# Update mean of q(F|s=1)
	resid = Y_slice - alpha_t_mu - intercept_mu_t - covariates_expected - G_slice*(components_expected)
	b_term = np.sum(tau_t_expected_val*G_slice*S_slice*resid*(1.0/sample_batch_fraction))
	new_var = 1.0/a_term
	new_mu = new_var*b_term
	if SVI == False:
		F_var = new_var
		F_mu = new_mu
	elif SVI == True:
		F_var = weighted_SVI_updated(F_var, new_var, step_size)
		F_mu = weighted_SVI_updated(F_mu, new_mu, step_size)
	return np.hstack((F_mu, F_var))

def outside_update_tau_t(tau_alpha, tau_beta, G_slice, Y_slice, N, U_S, V_S_t, F_S_t, intercept_t, V_S_t_squared, F_S_t_squared, U_S_squared, intercept_t_squared, S_slice, alpha_prior, beta_prior, cov, cov_squared, C_t, C_squared_t, alpha_mu_t, alpha_var_t, delta_t, delta_squared_t, sample_batch_fraction, step_size, SVI):
	# Compute Relevent expectations
	squared_factor_terms = U_S_squared@V_S_t_squared
	factor_terms = U_S@V_S_t

	squared_covariate_terms = cov_squared@C_squared_t
	covariate_terms = cov@C_t

	squared_non_genetic_factor_terms = U_S_squared@delta_squared_t
	non_genetic_factor_terms = U_S@delta_t

	alpha_t_squared = np.square(alpha_mu_t) + alpha_var_t

	# First add together square terms
	resid = S_slice*np.square(Y_slice) + S_slice*intercept_t_squared + S_slice*alpha_t_squared + S_slice*squared_covariate_terms + S_slice*squared_non_genetic_factor_terms + S_slice*np.square(G_slice)*(F_S_t_squared + squared_factor_terms)
	# Now add terms with Y
	resid = resid - (2.0*Y_slice*S_slice*(alpha_mu_t + intercept_t + covariate_terms + non_genetic_factor_terms + G_slice*factor_terms + G_slice*F_S_t))

	resid = resid + 2.0*alpha_mu_t*S_slice*(intercept_t + covariate_terms + non_genetic_factor_terms + G_slice*(factor_terms + F_S_t))

	resid = resid + 2.0*intercept_t*S_slice*(covariate_terms + non_genetic_factor_terms + G_slice*(factor_terms + F_S_t))

	resid = resid + 2.0*covariate_terms*S_slice*(non_genetic_factor_terms + G_slice*(factor_terms + F_S_t))

	resid = resid + 2.0*non_genetic_factor_terms*S_slice*(G_slice*(factor_terms + F_S_t))
	# Now add terms with factors
	resid = resid + 2.0*G_slice*factor_terms*G_slice*F_S_t*S_slice

	# WE ARE HERE
	# Now add terms with interactions between factors
	resid = resid + (np.square(G_slice)*S_slice*(factor_terms*factor_terms - np.sum(np.square(U_S*V_S_t),axis=1)))
	resid = resid + (S_slice*(covariate_terms*covariate_terms - np.sum(np.square(cov*C_t),axis=1)))
	resid = resid + (S_slice*(non_genetic_factor_terms*non_genetic_factor_terms - np.sum(np.square(U_S*delta_t),axis=1)))

	# Make Updates
	new_alpha = alpha_prior + ((np.sum(S_slice)/2.0)*(1.0/sample_batch_fraction))
	new_beta = beta_prior + ((np.sum(resid)/2.0)*(1.0/sample_batch_fraction))

	if SVI == False:
		tau_alpha = new_alpha
		tau_beta = new_beta
	elif SVI == True:
		tau_alpha = weighted_SVI_updated(tau_alpha, new_alpha, step_size)
		tau_beta = weighted_SVI_updated(tau_beta, new_beta, step_size)
	return np.hstack((tau_alpha, tau_beta))


def outside_update_S_t(S_copy_slice, G_slice, Y_slice, N, U_S, V_S_t, F_S_t, intercept_t, V_S_t_squared, F_S_t_squared, U_S_squared, intercept_t_squared, tau_t, theta_S_a_t, theta_S_b_t, cov, cov_squared, C_t, C_squared_t, sample_batch_fraction, step_size, SVI):
	# Compute Relevent expectations
	squared_factor_terms = U_S_squared@V_S_t_squared
	factor_terms = U_S@V_S_t

	squared_covariate_terms = cov_squared@C_squared_t
	covariate_terms = cov@C_t

	# Contribution of prior	
	ln_theta_S_expected_val = special.digamma(theta_S_a_t) - special.digamma(theta_S_a_t+theta_S_b_t)  # expectation of ln(1-X)
	ln_1_minus_theta_S_expected_val = special.digamma(theta_S_b_t) - special.digamma(theta_S_a_t+theta_S_b_t)
	# Contribution of data likelihood
	resid = np.square(Y_slice) + intercept_t_squared + squared_covariate_terms + np.square(G_slice)*(F_S_t_squared + squared_factor_terms)

	resid = resid - (2.0*Y_slice*(intercept_t + covariate_terms + G_slice*factor_terms + G_slice*F_S_t))

	resid = resid + 2.0*intercept_t*(covariate_terms + G_slice*(factor_terms + F_S_t))

	resid = resid + 2.0*covariate_terms*(G_slice*(factor_terms + F_S_t))

	resid = resid + 2.0*G_slice*factor_terms*G_slice*F_S_t

	resid = resid + (np.square(G_slice)*(factor_terms*factor_terms - np.sum(np.square(U_S*V_S_t),axis=1)))

	resid = resid + ((covariate_terms*covariate_terms - np.sum(np.square(cov*C_t),axis=1)))

	z_term = ln_theta_S_expected_val - ln_1_minus_theta_S_expected_val - (resid*tau_t/2.0) - .5*np.log((2.0*np.pi)/tau_t) - np.log((np.abs(Y_slice) < 1e-14))

	new_s = sigmoid_function(z_term)
	return new_s


class EQTL_FACTORIZATION_VI(object):
	def __init__(self, K=25, alpha=1e-3, beta=1e-3, a=1, b=1, gamma_v=1.0, gamma_u=1.0, max_iter=1000, delta_elbo_threshold=1e-8, SVI=False, parrallel_boolean=False, sample_batch_fraction=.3, learning_rate=.4, forgetting_rate=.01, num_test_cores=24, num_sample_cores=24, warmup_iters=50, output_root=''):
		self.alpha_prior = alpha
		self.beta_prior = beta
		self.a_prior = a 
		self.b_prior = b
		self.max_iter = max_iter
		self.K = K
		self.gamma_v = gamma_v
		self.gamma_u = gamma_u
		self.iter = 0
		self.delta_elbo_threshold = delta_elbo_threshold
		self.SVI = SVI
		self.parrallel_boolean = parrallel_boolean
		self.num_sample_cores = num_sample_cores
		self.num_test_cores = num_test_cores
		self.output_root = output_root
		self.warmup_iters = warmup_iters
		if SVI == True:
			self.sample_batch_fraction = sample_batch_fraction
			self.learning_rate = learning_rate
			self.forgetting_rate = forgetting_rate
		elif SVI == False:
			self.sample_batch_fraction = 1.0
			self.step_size = 1.0
	def fit(self, G, Y, cov, z):
		""" Fit the model.
			Args:
			G: A genotype matrix of floats with shape [num_samples, num_tests].
			Y: An expression matrix of floats with shape [num_samples, num_tests].
			z: groupings of length num_samples
		"""
		self.G_full = G
		self.Y_full = Y
		self.z_full = np.asarray(z)
		self.cov_full = cov
		self.initialize_variables()
		first_time = True
		#self.update_elbo()
		# Loop through VI iterations
		for vi_iter in range(self.max_iter):
			start_time = time.time()
			# Update step size (If using SVI)
			self.update_step_size()
			# Update parameter estimaters via coordinate ascent
			print('Step size: ' + str(self.step_size))
			#self.update_S()
			#self.update_theta_S()
			self.update_intercept()
			self.update_alpha()
			self.update_psi()
			self.update_C()
			self.update_F()
			#if vi_iter > self.warmup_iters:
			#	if first_time == True:
			#		self.U_mu = np.copy(self.U_init)
			#		first_time = False
			self.update_U()
			self.update_V()
			self.update_Delta()
			if vi_iter > 4:
				self.update_theta_U()
				#self.update_theta_V()
			self.update_tau()
			self.iter = self.iter + 1
			print('day0')
			print((self.U_mu*self.S_U)[790:795,:])
			print('day1')
			print((self.U_mu*self.S_U)[0:5,:])
			print('day2')
			print((self.U_mu*self.S_U)[1710:1715,:])
			print('day3')
			print((self.U_mu*self.S_U)[1761:1765,:])
			print('##############')
			print(self.theta_U_a/(self.theta_U_a + self.theta_U_b))
			print(self.theta_V_a/(self.theta_V_a + self.theta_V_b))



			# Remove irrelevent factors
			if np.mod(vi_iter, 5) == 0: 
				print(vi_iter)
				#self.remove_irrelevent_factors()
				np.savetxt(self.output_root + '_temper_U_S.txt', (self.U_mu*self.S_U), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_S_U.txt', (self.S_U), fmt="%s", delimiter='\t')
				#np.savetxt(self.output_root + '_temper2_S.txt', (self.S), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_V.txt', (self.V_mu*self.S_V), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_F.txt', (self.F_mu), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_Delta.txt', (self.Delta_mu), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_intercept.txt', (self.intercept_mu), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_tau.txt', (self.tau_alpha/self.tau_beta), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_theta_U.txt', self.theta_U_a/(self.theta_U_a + self.theta_U_b), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_theta_V.txt', self.theta_V_a/(self.theta_V_a + self.theta_V_b), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_C.txt', (self.C_mu), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_alpha.txt', (self.alpha_mu), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_psi.txt', (self.psi_alpha/self.psi_beta), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_iter.txt', np.asmatrix(vi_iter), fmt="%s", delimiter='\t')
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
		loading_sparsity = self.theta_V_a/(self.theta_V_a + self.theta_V_b)
		#num_factors = len(np.where(factor_sparsity > 0.01)[0])
		#factor_ordering = np.flip(np.argsort(factor_sparsity))[:num_factors]
		#factor_ordering = np.where(factor_sparsity > 0.01 and loading_sparsity > 0.01)[0]
		factor_ordering = []
		for k in range(self.K):
			if factor_sparsity[k] > .01 and loading_sparsity[k] > .01:
				factor_ordering.append(k)
		factor_ordering = np.asarray(factor_ordering)
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

		self.Delta_mu = self.Delta_mu[factor_ordering, :]
		self.Delta_var = self.Delta_var[factor_ordering, :]
		self.V_mu = self.V_mu[factor_ordering, :]
		self.V_var = self.V_var[factor_ordering, :]
		self.V_var_s_0 = self.V_var_s_0[factor_ordering, :]
		self.S_V = self.S_V[factor_ordering, :]
		self.theta_U_a = self.theta_U_a[factor_ordering]
		self.theta_U_b = self.theta_U_b[factor_ordering]
		self.theta_V_a = self.theta_V_a[factor_ordering]
		self.theta_V_b = self.theta_V_b[factor_ordering]
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
	def update_C(self):
		###################
		# UPDATE V
		###################
		# Precompute quantities
		U_S_expected_val = self.U_mu*self.S_U
		tau_expected_val = self.tau_alpha/self.tau_beta
		C_mu_copy = np.copy(self.C_mu)
		C_var_copy = np.copy(self.C_var)
		gamma_v = self.gamma_v
		V_S_expected_val = self.V_mu*self.S_V

		# Keep track of variables
		C_update_data = []

		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				C_update_data.append(outside_update_C_t(C_mu_copy[:, test_index], C_var_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.K, U_S_expected_val, V_S_expected_val[:, test_index], self.F_mu[test_index], self.intercept_mu[test_index], self.S[:, test_index], gamma_v, tau_expected_val[test_index], self.cov, self.alpha_big_mu[:, test_index], self.Delta_mu[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		elif self.parrallel_boolean == True:
			C_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_C_t)(C_mu_copy[:, test_index], C_var_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.K, U_S_expected_val, V_S_expected_val[:, test_index], self.F_mu[test_index], self.intercept_mu[test_index], self.S[:, test_index], gamma_v, tau_expected_val[test_index], self.cov,self.alpha_big_mu[:, test_index], self.Delta_mu[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))

		J = self.cov.shape[1]
		# Convert to array
		C_update_data = np.asarray(C_update_data).T
		# Fill in data structures
		self.C_mu = C_update_data[(J*0):(1*J), :]
		self.C_var = C_update_data[(J*1):(2*J), :]
	def update_alpha(self):
		U_S_expected_val = self.U_mu*self.S_U
		tau_expected_val = self.tau_alpha/self.tau_beta
		psi_expected_val = self.psi_alpha/self.psi_beta
		alpha_mu_copy = np.copy(self.alpha_mu)
		alpha_var_copy = np.copy(self.alpha_var)
		V_S_expected_val = self.V_mu*self.S_V

		alpha_update_data = []
		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				alpha_update_data.append(outside_update_alpha_t(alpha_mu_copy[:, test_index], alpha_var_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.I, U_S_expected_val, V_S_expected_val[:, test_index], self.F_mu[test_index], self.intercept_mu[test_index], self.cov, self.C_mu[:, test_index], tau_expected_val[test_index], psi_expected_val[test_index], self.Delta_mu[:, test_index], self.individual_to_sample_indices, self.individual_to_number_full_indices, self.step_size, self.SVI))
		else:
			alpha_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_alpha_t)(alpha_mu_copy[:, test_index], alpha_var_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.I, U_S_expected_val, V_S_expected_val[:, test_index], self.F_mu[test_index], self.intercept_mu[test_index], self.cov, self.C_mu[:, test_index], tau_expected_val[test_index], psi_expected_val[test_index], self.Delta_mu[:, test_index], self.individual_to_sample_indices, self.individual_to_number_full_indices, self.step_size, self.SVI) for test_index in range(self.T))
		
		alpha_update_data = np.transpose(np.asarray(alpha_update_data))
		self.alpha_mu = alpha_update_data[:(self.I),:]
		self.alpha_var = alpha_update_data[(self.I):, :]
		# Now fill in big matrix
		self.alpha_big_mu = np.zeros((self.N, self.T))
		self.alpha_big_var = np.zeros((self.N, self.T))
		for sample_num, z_label in enumerate(self.z):
			self.alpha_big_mu[sample_num,:] = self.alpha_mu[self.z_mapping[z_label], :]
			self.alpha_big_var[sample_num,:] = self.alpha_var[self.z_mapping[z_label], :]

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
		S_V_copy = np.copy(self.S_V)
		V_var_s_0_copy = np.copy(self.V_var_s_0)
		gamma_v = self.gamma_v

		# Keep track of variables
		V_update_data = []

		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				V_update_data.append(outside_update_V_t(V_mu_copy[:, test_index], V_var_copy[:, test_index], S_V_copy[:, test_index], V_var_s_0_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.K, U_S_expected_val, U_S_squared_expected_val, self.F_mu[test_index], self.intercept_mu[test_index], self.S[:, test_index], gamma_v, tau_expected_val[test_index], self.cov, self.C_mu[:, test_index], self.alpha_big_mu[:, test_index], self.Delta_mu[:, test_index], self.theta_V_a, self.theta_V_b, self.sample_batch_fraction, self.step_size, self.SVI))
		elif self.parrallel_boolean == True:
			V_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_V_t)(V_mu_copy[:, test_index], V_var_copy[:, test_index], S_V_copy[:, test_index], V_var_s_0_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.K, U_S_expected_val, U_S_squared_expected_val, self.F_mu[test_index], self.intercept_mu[test_index], self.S[:, test_index], gamma_v, tau_expected_val[test_index], self.cov, self.C_mu[:, test_index], self.alpha_big_mu[:, test_index], self.Delta_mu[:, test_index], self.theta_V_a, self.theta_V_b, self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))

		# Convert to array
		V_update_data = np.asarray(V_update_data).T
		# Fill in data structures
		self.V_mu = V_update_data[(self.K*0):(1*self.K), :]
		self.V_var = V_update_data[(self.K*1):(2*self.K), :]
		#if self.iter > 4:
			#self.S_V = V_update_data[(self.K*2):(3*self.K), :]
		#self.V_var_s_0 = V_update_data[(self.K*3):(4*self.K), :]
	def update_Delta(self):
		###################
		# UPDATE V
		###################
		# Precompute quantities
		U_S_expected_val = self.U_mu*self.S_U
		U_S_squared_expected_val = (np.square(self.U_mu) + self.U_var)*self.S_U
		tau_expected_val = self.tau_alpha/self.tau_beta
		Delta_mu_copy = np.copy(self.Delta_mu)
		Delta_var_copy = np.copy(self.Delta_var)
		gamma_v = self.gamma_v

		# Keep track of variables
		Delta_update_data = []

		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				Delta_update_data.append(outside_update_Delta_t(Delta_mu_copy[:, test_index], Delta_var_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.K, U_S_expected_val, U_S_squared_expected_val, self.F_mu[test_index], self.intercept_mu[test_index], self.S[:, test_index], gamma_v, tau_expected_val[test_index], self.cov, self.C_mu[:, test_index], self.alpha_big_mu[:, test_index], self.V_mu[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		elif self.parrallel_boolean == True:
			Delta_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_Delta_t)(Delta_mu_copy[:, test_index], Delta_var_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.K, U_S_expected_val, U_S_squared_expected_val, self.F_mu[test_index], self.intercept_mu[test_index], self.S[:, test_index], gamma_v, tau_expected_val[test_index], self.cov, self.C_mu[:, test_index], self.alpha_big_mu[:, test_index], self.V_mu[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))

		# Convert to array
		Delta_update_data = np.asarray(Delta_update_data).T
		# Fill in data structures
		self.Delta_mu = Delta_update_data[(self.K*0):(1*self.K), :]
		self.Delta_var = Delta_update_data[(self.K*1):(2*self.K), :]
	def update_U(self):
		if self.SVI == True:
			print('error: SVI currently not implemented')
			pdb.set_trace()
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
			# Convert random effects matrix to samplesXtests instead of groupsXtest
			self.alpha_big_mu = np.zeros((self.N, self.T))
			self.alpha_big_var = np.zeros((self.N, self.T))
			for sample_num, z_label in enumerate(self.z):
				self.alpha_big_mu[sample_num,:] = self.alpha_mu[self.z_mapping[z_label], :]
				self.alpha_big_var[sample_num,:] = self.alpha_var[self.z_mapping[z_label], :]
			self.individual_to_sample_indices = []
			for ii in range(self.I):
				# z_label corresponding to this individual
				z_label = self.z_inverse_mapping[ii]
				sample_indices = np.where(np.asarray(self.z) == z_label)[0]
				if len(sample_indices) == 0:
					print('sample indices update U assumption error')
					pdb.set_trace()
				self.individual_to_sample_indices.append(sample_indices)
		#G_squared = np.square(self.G)
		###################
		# UPDATE U
		###################
		V_S_squared_expected_val = (np.square(self.V_mu) + self.V_var)*self.S_V
		V_S_expected_val = self.V_mu*self.S_V
		Delta_squared_expected_val = np.square(self.Delta_mu) + self.Delta_var
		U_mu_copy = np.copy(self.U_mu)
		S_U_copy = np.copy(self.S_U)
		U_var_copy = np.copy(self.U_var)
		U_var_s_0_copy = np.copy(self.U_var_s_0)
		U_update_data = []

		# Don't parrallelize
		if self.parrallel_boolean == False:
			for sample_index in range(self.N):
				U_update_data.append(outside_update_U_n(U_mu_copy[sample_index,:],S_U_copy[sample_index,:], U_var_copy[sample_index,:],U_var_s_0_copy[sample_index,:], self.G[sample_index, :], self.Y[sample_index, :], self.K, V_S_expected_val, V_S_squared_expected_val, self.F_mu, self.S[sample_index, :], self.intercept_mu, self.gamma_u, self.tau_alpha/self.tau_beta, self.cov[sample_index,:], self.C_mu, self.alpha_big_mu[sample_index,:], self.theta_U_a, self.theta_U_b, self.Delta_mu, Delta_squared_expected_val))
		# Parrallelize
		elif self.parrallel_boolean == True:
			U_update_data = Parallel(n_jobs=self.num_sample_cores)(delayed(outside_update_U_n)(U_mu_copy[sample_index,:],S_U_copy[sample_index,:], U_var_copy[sample_index,:], U_var_s_0_copy[sample_index,:], self.G[sample_index, :], self.Y[sample_index, :], self.K, V_S_expected_val, V_S_squared_expected_val, self.F_mu, self.S[sample_index,:], self.intercept_mu, self.gamma_u, self.tau_alpha/self.tau_beta, self.cov[sample_index,:], self.C_mu, self.alpha_big_mu[sample_index,:], self.theta_U_a, self.theta_U_b, self.Delta_mu, Delta_squared_expected_val) for sample_index in range(self.N))

		# Convert to array
		U_update_data = np.asarray(U_update_data)
		# Fill in data structures
		self.U_mu = U_update_data[:,(self.K*0):(1*self.K)]
		if self.iter > 4:
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
		V_S_expected_val = self.V_mu*self.S_V
		tau_expected_val = self.tau_alpha/self.tau_beta
		F_mu_copy = np.copy(self.F_mu)
		F_var_copy = np.copy(self.F_var)

		F_update_data = []
		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				F_update_data.append(outside_update_F_t(F_mu_copy[test_index], F_var_copy[test_index], self.G[:, test_index], self.Y[:, test_index], U_S_expected_val, V_S_expected_val[:,test_index], self.intercept_mu[test_index], self.gamma_v, tau_expected_val[test_index], self.S[:, test_index], self.cov, self.C_mu[:, test_index], self.alpha_big_mu[:, test_index], self.Delta_mu[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		if self.parrallel_boolean == True:
			F_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_F_t)(F_mu_copy[test_index], F_var_copy[test_index], self.G[:, test_index], self.Y[:, test_index], U_S_expected_val, V_S_expected_val[:,test_index], self.intercept_mu[test_index], self.gamma_v, tau_expected_val[test_index], self.S[:, test_index], self.cov, self.C_mu[:, test_index], self.alpha_big_mu[:, test_index], self.Delta_mu[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))
		F_update_data = np.asarray(F_update_data)
		self.F_mu = F_update_data[:,0]
		self.F_var = F_update_data[:,1]
	def update_intercept(self):
		U_S_expected_val = self.U_mu*self.S_U
		V_S_expected_val = self.V_mu*self.S_V
		tau_expected_val = self.tau_alpha/self.tau_beta
		intercept_mu_copy = np.copy(self.intercept_mu)
		intercept_var_copy = np.copy(self.intercept_var)

		intercept_update_data = []
		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				intercept_update_data.append(outside_update_intercept_t(intercept_mu_copy[test_index], intercept_var_copy[test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S_expected_val, V_S_expected_val[:,test_index], self.F_mu[test_index], tau_expected_val[test_index], self.S[:, test_index], self.cov, self.C_mu[:, test_index], self.alpha_big_mu[:, test_index], self.Delta_mu[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		if self.parrallel_boolean == True:
				intercept_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_intercept_t)(intercept_mu_copy[test_index], intercept_var_copy[test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S_expected_val, V_S_expected_val[:,test_index], self.F_mu[test_index], tau_expected_val[test_index], self.S[:, test_index], self.cov, self.C_mu[:, test_index], self.alpha_big_mu[:, test_index], self.Delta_mu[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))

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
	def update_gamma_V(self):
		V_S_squared_expected_val = np.square(self.V_mu) + self.V_var
		# Loop through tests
		'''
		for k in range(self.K):
			self.gamma_V_alpha[k] = self.alpha_prior + (self.T/2.0)
			self.gamma_V_beta[k] = self.beta_prior + np.sum(V_S_squared_expected_val[k, :])/2.0
		'''
		self.gamma_V_alpha = self.alpha_prior + (self.T*self.K/2.0)
		self.gamma_V_beta = self.beta_prior + np.sum(V_S_squared_expected_val)/2.0
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
	def update_theta_V(self):
		# Loop through factors
		for k in range(self.K):
			#self.theta_U_a[k] = self.a_prior + np.sum(self.S_U[:,k])
			#self.theta_U_b[k] = self.b_prior + self.N - np.sum(self.S_U[:,k])
			new_theta_V_a = self.a_prior + (1.0/self.sample_batch_fraction)*np.sum(self.S_V[k,:])
			new_theta_V_b = self.b_prior + (1.0/self.sample_batch_fraction)*(self.T - np.sum(self.S_V[k,:]))
			if self.SVI == False:
				self.theta_V_a[k] = new_theta_V_a
				self.theta_V_b[k] = new_theta_V_b
			elif self.SVI == True:
				self.theta_V_a[k] = weighted_SVI_updated(self.theta_V_a[k], new_theta_V_a, self.step_size)
				self.theta_V_b[k] = weighted_SVI_updated(self.theta_V_b[k], new_theta_V_b, self.step_size)
	def update_theta_S(self):
		for t in range(self.T):
			new_theta_S_t_a = self.a_prior + (1.0/self.sample_batch_fraction)*np.sum(self.S[:,t])
			new_theta_S_t_b = self.b_prior + (1.0/self.sample_batch_fraction)*(self.N - np.sum(self.S[:,t]))
			if self.SVI == False:
				self.theta_S_a[t] = new_theta_S_t_a
				self.theta_S_b[t] = new_theta_S_t_b
			elif self.SVI == True:
				self.theta_S_a[t] = weighted_SVI_updated(self.theta_S_a[t], new_theta_S_t_a, self.step_size)
				self.theta_S_b[t] = weighted_SVI_updated(self.theta_S_b[t], new_theta_S_t_b, self.step_size)

	def update_psi(self):
		alpha_squared_expected_value = np.square(self.alpha_mu) + self.alpha_var
		# Loop through tests
		for test_index in range(self.T):
			self.psi_alpha[test_index] = self.alpha_prior + (self.I/2.0)
			self.psi_beta[test_index] = self.beta_prior + (np.sum(alpha_squared_expected_value[:,test_index])/2.0)
	def update_tau(self):
		tau_alpha_copy = np.copy(self.tau_alpha)
		tau_beta_copy = np.copy(self.tau_beta)

		# Precompute quantities
		F_S_squared = np.square(self.F_mu) + self.F_var
		intercept_squared = np.square(self.intercept_mu) + self.intercept_var
		V_S_squared = (np.square(self.V_mu) + self.V_var)*self.S_V
		U_S_squared = ((np.square(self.U_mu) + self.U_var)*self.S_U)
		U_S = (self.U_mu)*self.S_U
		V_S = self.V_mu*self.S_V
		cov_squared = np.square(self.cov)
		C_squared = np.square(self.C_mu) + self.C_var
		Delta_squared = np.square(self.Delta_mu) + self.Delta_var
		# Loop through tests
		tau_update_data = []
		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				tau_update_data.append(outside_update_tau_t(tau_alpha_copy[test_index], tau_beta_copy[test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S, V_S[:,test_index], self.F_mu[test_index], self.intercept_mu[test_index], V_S_squared[:, test_index], F_S_squared[test_index], U_S_squared, intercept_squared[test_index], self.S[:, test_index], self.alpha_prior, self.beta_prior, self.cov, cov_squared, self.C_mu[:, test_index], C_squared[:, test_index], self.alpha_big_mu[:, test_index], self.alpha_big_var[:, test_index], self.Delta_mu[:, test_index], Delta_squared[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		if self.parrallel_boolean == True:
			tau_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_tau_t)(tau_alpha_copy[test_index], tau_beta_copy[test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S, V_S[:,test_index], self.F_mu[test_index], self.intercept_mu[test_index], V_S_squared[:, test_index], F_S_squared[test_index], U_S_squared, intercept_squared[test_index], self.S[:, test_index], self.alpha_prior, self.beta_prior, self.cov, cov_squared, self.C_mu[:, test_index], C_squared[:, test_index], self.alpha_big_mu[:, test_index], self.alpha_big_var[:, test_index], self.Delta_mu[:, test_index], Delta_squared[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))
		tau_update_data = np.asarray(tau_update_data)
		self.tau_alpha = tau_update_data[:,0]
		self.tau_beta = tau_update_data[:,1]
	def update_S(self):
		S_copy = np.copy(self.S)
		# Precompute quantities
		F_S_squared = np.square(self.F_mu) + self.F_var
		intercept_squared = np.square(self.intercept_mu) + self.intercept_var
		V_S_squared = np.square(self.V_mu) + self.V_var
		U_S_squared = ((np.square(self.U_mu) + self.U_var))
		C_squared = np.square(self.C_mu) + self.C_var
		U_S = (self.U_mu)
		cov_squared = np.square(self.cov)
		tau_expected_val = self.tau_alpha/self.tau_beta
		# Loop through tests
		S_update_data = []
		if self.parrallel_boolean == False:
			for test_index in range(self.T):
				S_update_data.append(outside_update_S_t(S_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S, self.V_mu[:,test_index], self.F_mu[test_index], self.intercept_mu[test_index], V_S_squared[:, test_index], F_S_squared[test_index], U_S_squared, intercept_squared[test_index], tau_expected_val[test_index], self.theta_S_a[test_index], self.theta_S_b[test_index], self.cov, cov_squared, self.C_mu[:, test_index], C_squared[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI))
		if self.parrallel_boolean == True:
			S_update_data = Parallel(n_jobs=self.num_test_cores)(delayed(outside_update_S_t)(S_copy[:, test_index], self.G[:, test_index], self.Y[:, test_index], self.N, U_S, self.V_mu[:,test_index], self.F_mu[test_index], self.intercept_mu[test_index], V_S_squared[:, test_index], F_S_squared[test_index], U_S_squared, intercept_squared[test_index], tau_expected_val[test_index], self.theta_S_a[test_index], self.theta_S_b[test_index], self.cov, cov_squared, self.C_mu[:, test_index], C_squared[:, test_index], self.sample_batch_fraction, self.step_size, self.SVI) for test_index in range(self.T))
		self.S = np.transpose(np.asarray(S_update_data))

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
		svi_sample_indices = []
		for individual_index in range(self.I):
			individual_indices = self.individual_to_sample_indices_full[individual_index]
			n_i = len(individual_indices)
			subset_n_i = int(np.ceil(n_i*self.sample_batch_fraction))
			randos = np.random.choice(individual_indices, size=int(subset_n_i), replace=False)
			for rando in randos:
				svi_sample_indices.append(rando)
		svi_sample_indices = np.asarray(svi_sample_indices)
		if len(svi_sample_indices) > self.N:
			svi_sample_indices = np.random.choice(svi_sample_indices,size=self.N, replace=False)
		if len(svi_sample_indices) != self.N:
			print('svi_sample_indices assumption error1')
			pdb.set_trace()
		if len(np.unique(svi_sample_indices)) != self.N:
			print('svi_sample_indices assumption error2')
			pdb.set_trace()
		return svi_sample_indices
		# return np.random.choice(self.N_full, size=self.N, replace=False)
	def initialize_variables(self):
		# Initialize array to keep track of ELBO
		self.elbo = []

		# Add model dimensions to object
		self.N_full = self.Y_full.shape[0]
		self.T_full = self.Y_full.shape[1]
		self.num_cov = self.cov_full.shape[1]

		# Do standard variational inference
		if self.SVI == False:
			self.N = self.Y_full.shape[0]
			self.T = self.Y_full.shape[1]
			self.G = np.copy(self.G_full)
			self.Y = np.copy(self.Y_full)
			self.z = np.copy(self.z_full)
			self.cov = np.copy(self.cov_full)
			pca = sklearn.decomposition.PCA(n_components=self.K, whiten=True)
			pca.fit(np.random.randn(self.N, 9999).T)
			self.U_mu = pca.components_.T
			for k in range(self.K):
				self.U_mu[:,k] = ((self.U_mu[:,k]-np.mean(self.U_mu[:,k]))/(np.std(self.U_mu[:,k])))
			self.U_var = np.ones((self.N, self.K))*(1.0/self.gamma_u)*100.0
			self.U_var_s_0 = np.ones((self.N, self.K))*(1.0/self.gamma_u)
			self.S = np.ones((self.N,self.T))
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
			#pca = sklearn.decomposition.PCA(n_components=self.K, whiten=True)
			#pca.fit(np.random.randn(self.N_full, 9999).T)
			#self.U_mu_full = pca.components_.T
			#for k in range(self.K):
			#	self.U_mu_full[:,k] = ((self.U_mu_full[:,k]-np.mean(self.U_mu_full[:,k]))/np.std(self.U_mu_full[:,k]))
			self.U_mu_full = np.random.randn(self.N_full, self.K)
			self.U_var_full = np.ones((self.N_full, self.K))*(1.0/self.gamma_v)
			self.S_full = np.ones((self.N_full,self.T_full))
			self.U_mu = np.copy(self.U_mu_full[svi_sample_indices, :])
			self.U_var = np.copy(self.U_var_full[svi_sample_indices, :])*100.0
			self.S = np.copy(self.S_full[svi_sample_indices, :])

		self.U_init = np.copy(self.U_mu)
		# Random effects
		self.z_mapping = {}
		self.z_inverse_mapping = {}
		# Create mapping from grouping to index
		_, idx = np.unique(self.z, return_index=True)
		unique_groups = np.asarray(self.z)[np.sort(idx)]
		for i, label in enumerate(unique_groups):
			self.z_mapping[label] = i
			self.z_inverse_mapping[i] = label
		self.I = len(np.unique(self.z))
		self.individual_to_sample_indices = []
		self.individual_to_sample_indices_full = []
		self.individual_to_number_full_indices = []
		for ii in range(self.I):
			# z_label corresponding to this individual
			z_label = self.z_inverse_mapping[ii]
			sample_indices = np.where(np.asarray(self.z) == z_label)[0]
			self.individual_to_sample_indices.append(sample_indices)
			self.individual_to_sample_indices_full.append(sample_indices)
			self.individual_to_number_full_indices.append(float(len(sample_indices)))
		# Variances
		self.psi_alpha = np.ones(self.T)*self.alpha_prior
		self.psi_beta = np.ones(self.T)*self.beta_prior
		# Random effects
		self.alpha_mu = np.zeros((self.I, self.T))
		self.alpha_var = (np.zeros((self.I, self.T)) + 1.0)
		# Convert random effects matrix to samplesXtests instead of groupsXtest
		self.alpha_big_mu = np.zeros((self.N, self.T))
		self.alpha_big_var = np.zeros((self.N, self.T))
		for sample_num, z_label in enumerate(self.z):
			self.alpha_big_mu[sample_num,:] = self.alpha_mu[self.z_mapping[z_label], :]
			self.alpha_big_var[sample_num,:] = self.alpha_var[self.z_mapping[z_label], :]

		#pca = sklearn.decomposition.PCA(n_components=self.K, whiten=True)
		#pca.fit(np.random.randn(self.T, 9999).T)
		#self.V_mu = pca.components_
		#for k in range(self.K):
		#	self.V_mu[k,:] = ((self.V_mu[k,:]-np.mean(self.V_mu[k,:]))/np.std(self.V_mu[k,:]))
		#self.V_var = np.ones((self.K, self.T))*(1.0/self.gamma_v)

		V_init, betas, C_init, intercepts, delta_init = run_linear_model_for_initialization(self.Y, self.G, self.cov, self.U_init)
		#V_init, betas, C_init, intercepts = run_linear_mixed_model_for_initialization(self.Y, self.G, self.cov, self.U_mu, self.z)
		self.V_mu = np.transpose(V_init)
		self.V_var = np.ones((self.K, self.T))*100.0
		self.V_var_s_0 = np.ones((self.K, self.T))*(1.0/self.gamma_u)

		self.Delta_mu = np.transpose(delta_init)
		self.Delta_var = np.ones((self.K, self.T))*100.0
		

		self.S = np.ones((self.N,self.T))
		#self.F_mu = eqtl_vi_init.F_mu
		self.F_mu = betas
		self.F_var = np.ones(self.T)*100.0

		self.C_mu = np.transpose(C_init)
		self.C_var = np.ones(self.C_mu.shape)*100.0
		# Masks
		self.S_U = np.ones((self.N, self.K))
		self.S_V = np.ones((self.K, self.T))

		# Intercepts
		# self.intercept_mu = eqtl_vi_init.intercept_mu
		self.intercept_mu = intercepts
		self.intercept_var = np.ones(self.T)*100.0

		#self.gamma_V_alpha = np.ones(self.K)
		#self.gamma_V_beta = np.ones(self.K)
		#self.gamma_V_alpha = 1.0
		#self.gamma_V_beta = 1.0

		# Variances
		self.tau_alpha = np.ones(self.T)*self.alpha_prior
		self.tau_beta = np.ones(self.T)*self.beta_prior
		# bernoulli probs
		self.theta_S_a = np.ones(self.T)*self.a_prior + 9
		self.theta_S_b = np.ones(self.T)*self.b_prior

		# bernoulli probs
		self.theta_U_a = np.ones(self.K)*self.a_prior + 99
		self.theta_U_b = np.ones(self.K)*self.b_prior
		self.theta_V_a = np.ones(self.K)*self.a_prior + 99
		self.theta_V_b = np.ones(self.K)*self.b_prior	

		# Initialize other variables based around U
		self.step_size = 1.0
		#self.update_S()
		#self.update_theta_S()		
		#self.update_V()
		#self.update_intercept()
		#self.update_F()
		#self.update_tau()

		if self.SVI == True:
			self.sample_batch_fraction = actual_sample_batch_fraction
			self.N = round((self.Y_full.shape[0])*self.sample_batch_fraction)

		#self.gamma_U_alpha = np.ones(self.K)*self.alpha_prior
		#self.gamma_U_beta = np.ones(self.K)*self.beta_prior

