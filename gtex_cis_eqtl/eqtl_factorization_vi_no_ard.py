import numpy as np 
import os
import sys
import pdb
import scipy.special as special

def sigmoid_function(x):
	return 1.0/(1.0 + np.exp(-x))

def compute_kl_divergence_of_gaussian_bernoulli(S, W_mu, W_var, W_var_s_0, precision, p, K):
	num_feat = W_mu.shape[1]
	# Relevent expectations
	log_gamma_expected = np.log(precision)
	gamma_expected = precision
	log_theta_expected_val = np.log(p)
	log_1_minus_theta_expected_val = np.log(1.0-p)
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
		likelihood_term_a = likelihood_term_a + (num_feat/2.0)*log_gamma_expected
		likelihood_term_b = likelihood_term_b - (np.sum(W_squared_expected_val[k,:])*gamma_expected/2.0)
		likelihood_term_c = likelihood_term_c + (np.sum(S[k,:])*log_theta_expected_val)
		likelihood_term_d = likelihood_term_d + (np.sum(1.0-S[k,:])*log_1_minus_theta_expected_val)
		#entropy_term_a = entropy_term_a + (.5)*np.sum(np.log((S[k,:]*W_var[k,:]) + ((1-S[k,:])*W_var_s_0[k,:])))
		entropy_term_a = entropy_term_a - (.5)*np.sum(S[k,:]*np.log(W_var[k,:]) + (1.0-S[k,:])*np.log(W_var_s_0[k,:]))

		temp_term_b = ((1.0-S[k,:])*np.log(1.0-S[k,:])) + (S[k,:]*np.log(S[k,:]))
		temp_term_b[np.isnan(temp_term_b)] = 0.
		entropy_term_b = entropy_term_b + np.sum(temp_term_b)
		#entropy_term_b = entropy_term_b + np.sum(((1.0-S[k,:])*np.log(log_const+1.0-S[k,:])) - (S[k,:]*np.log(log_const+S[k,:])))
	
	kl_divergence = entropy_term_a + entropy_term_b - likelihood_term_a - likelihood_term_b - likelihood_term_c - likelihood_term_d
	return kl_divergence

def compute_kl_divergence_of_gaussian(W_mu, W_var, expected_gamma_v, K):
	num_feat = W_mu.shape[1]
	# Relevent expectations
	W_squared_expected_val = ((np.square(W_mu) + W_var))

	# Initialize variables
	likelihood_term_a = 0
	entropy_term_a = 0

	for k in range(K):
		likelihood_term_a = likelihood_term_a - (np.sum(W_squared_expected_val[k,:])*expected_gamma_v/2.0)
		#entropy_term_a = entropy_term_a + (.5)*np.sum(np.log((S[k,:]*W_var[k,:]) + ((1-S[k,:])*W_var_s_0[k,:])))
		entropy_term_a = entropy_term_a - (.5)*np.sum(np.log(W_var[k,:]))

		#entropy_term_b = entropy_term_b + np.sum(((1.0-S[k,:])*np.log(log_const+1.0-S[k,:])) - (S[k,:]*np.log(log_const+S[k,:])))
	
	kl_divergence = entropy_term_a - likelihood_term_a

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

class EQTL_FACTORIZATION_VI(object):
	def __init__(self, K, alpha, beta, p_u, lambda_u, lambda_v, max_iter, delta_elbo_threshold):
		self.alpha_prior = alpha
		self.beta_prior = beta
		self.max_iter = max_iter
		self.K = K
		self.iter = 0
		self.lambda_v = lambda_v # Precision of prior on weight matrix (V)
		self.lambda_u = lambda_u
		self.p_u = p_u
		self.delta_elbo_threshold = delta_elbo_threshold
	def fit(self, G, Y, z):
		""" Fit the model.
			Args:
			G: A genotype matrix of floats with shape [num_samples, num_tests].
			Y: An expression matrix of floats with shape [num_samples, num_tests].
			z: groupings of length num_samples
		"""
		self.G = G
		self.Y = Y
		self.z = z
		self.initialize_variables()
		self.update_elbo()
		# Loop through VI iterations
		for vi_iter in range(self.max_iter):
			self.iter = self.iter + 1
			old_U = np.copy(self.U_mu) 
			old_V = np.copy(self.V_mu)
			old_alpha = np.copy(self.alpha_mu)
			# Update parameter estimaters via coordinate ascent
			self.update_intercept()
			self.update_V()
			self.update_U()
			#self.update_alpha()
			#self.update_psi()
			self.update_tau()
			# Compute ELBO after update
			


			#if self.iter % 10 == 0:
			self.update_elbo()
			elbo_iter = len(self.elbo) - 1
			current_elbo = self.elbo[elbo_iter]
			delta_elbo = (current_elbo - self.elbo[elbo_iter-1])
			if delta_elbo < 0:
				print('elbo decrease')
				np.savetxt('U_S.txt', self.U_mu*self.S_U, fmt="%s", delimiter='\t')
				pdb.set_trace()
			print('delta ELBO: ' + str(delta_elbo))	

			####################
			pred_U = self.U_mu*self.S_U
			print(pred_U[0:4,:])
			print(pred_U[641:645,:])
			print(pred_U[1090:1094,:])
			print(pred_U[1640:1644,:])
			print(self.tau_alpha/self.tau_beta)
			print('Variational Inference iteration: ' + str(vi_iter))
			print('##############')

	def update_U(self):
		for sample_index in range(self.N):
			# Could pararallelize here (not at lower level)
			for k in range(self.K):
				self.update_U_nk(sample_index, k)
	def update_U_nk(self, sample_index, k):
		# Compute relevent expectations
		gamma_k_u_expected_val = self.lambda_u
		tau_expected_val = self.tau_alpha/self.tau_beta
		V_S_expected_val = self.V_mu
		U_S_expected_val = self.U_mu[sample_index,:]*self.S_U[sample_index,:]
		V_k_S_k_expected_val = V_S_expected_val[k,:]
		V_k_S_k_squared_expected_val = (np.square(self.V_mu[k,:]) + self.V_var[k,:])
		individual_index = self.z_mapping[self.z[sample_index]]
		alpha_i_expected_val = self.alpha_mu[individual_index,:]
		theta_U_expected_val = self.p_u
		ln_theta_U_expected_val = np.log(self.p_u)  # expectation of ln(1-X)
		ln_1_minus_theta_U_expected_val = np.log(1.0-self.p_u)
		# Compute expectations on other components
		other_components_expected = np.zeros(self.T)
		for j in range(self.K):
			if j != k:
				other_components_expected = other_components_expected + V_S_expected_val[j,:]*U_S_expected_val[j]
		# Update variance of q(U|s=1)
		a_term = np.sum(tau_expected_val*np.square(self.G[sample_index,:])*V_k_S_k_squared_expected_val) + gamma_k_u_expected_val
		self.U_var[sample_index, k] = 1.0/a_term
		# Update variance of q(U|s=0)
		self.U_var_s_0[sample_index, k] = 1.0/self.lambda_u
		# Update mean of q(U|s=1)
		resid = self.Y[sample_index,:] - self.intercept_mu - alpha_i_expected_val - self.G[sample_index,:]*(other_components_expected)
		b_term = np.sum(tau_expected_val*self.G[sample_index,:]*V_k_S_k_expected_val*resid)
		self.U_mu[sample_index, k] = self.U_var[sample_index, k]*b_term
		# Now update q(S_U=1)
		#z_term = np.log(theta_U_expected_val/(1.0-theta_U_expected_val)) + .5*np.log(gamma_k_u_expected_val) - .5*np.log(a_term) + (np.square(b_term)/(2.0*a_term))
		z_term = ln_theta_U_expected_val - ln_1_minus_theta_U_expected_val + .5*np.log(gamma_k_u_expected_val) - .5*np.log(a_term) + (np.square(b_term)/(2.0*a_term))
		self.S_U[sample_index, k] = sigmoid_function(z_term)
	def update_V(self):
		for test_index in range(self.T):
			# Could pararallelize here (not at lower level)
			for k in range(self.K):
				self.update_V_kt(k, test_index)
	def update_V_kt(self, k, test_index):
		# Compute relevent expectations
		gamma_k_v_expected_val = self.lambda_v
		tau_t_expected_val = self.tau_alpha[test_index]/self.tau_beta[test_index]
		U_S_expected_val = self.U_mu*self.S_U
		U_k_S_k_expected_val = U_S_expected_val[:,k]
		U_k_S_k_squared_expected_val = (np.square(self.U_mu[:,k]) + self.U_var[:,k])*self.S_U[:,k]
		V_S_t_expected_val = self.V_mu[:,test_index]
		alpha_t_expected_val = self.alpha_big_mu[:,test_index]

		# Compute expectations on other components
		other_components_expected = np.zeros(self.N)
		for j in range(self.K):
			if j != k:
				other_components_expected = other_components_expected + U_S_expected_val[:,j]*V_S_t_expected_val[j]
		# Update variance of q(V|s=1)
		a_term = gamma_k_v_expected_val + tau_t_expected_val*np.sum(np.square(self.G[:,test_index])*U_k_S_k_squared_expected_val)
		self.V_var[k, test_index] = 1.0/a_term
		# Update variance of q(V|s=1)
		#self.V_var_s_0[k,test_index] = self.gamma_V_beta[k]/self.gamma_V_alpha[k]
		# Update mean of q(U|s=1)
		resid = self.Y[:,test_index] - self.intercept_mu[test_index] - alpha_t_expected_val - self.G[:,test_index]*(other_components_expected)
		b_term = np.sum(tau_t_expected_val*self.G[:,test_index]*U_k_S_k_expected_val*resid)
		self.V_mu[k, test_index] = self.V_var[k, test_index]*b_term
		# Now update q(S_V=1)
		# z_term = ln_theta_V_expected_val - ln_1_minus_theta_V_expected_val + .5*np.log(gamma_k_v_expected_val) - .5*np.log(a_term) + (np.square(b_term)/(2.0*a_term))
		# self.S_V[k, test_index] = sigmoid_function(z_term)
	def update_alpha(self):
		# Update alphas
		for test_index in range(self.T):
			self.update_alpha_t(test_index)
		# BOTTOM LOOP MAPPING BACK TO ALPHA_BIG
		for sample_num, z_label in enumerate(self.z):
			self.alpha_big_mu[sample_num,:] = self.alpha_mu[self.z_mapping[z_label], :]
			self.alpha_big_var[sample_num,:] = self.alpha_var[self.z_mapping[z_label], :]

	def update_alpha_t(self, test_index):
		# Compute relevent expectations
		tau_t_expected_val = self.tau_alpha[test_index]/self.tau_beta[test_index]
		psi_t_expected_val = self.psi_alpha[test_index]/self.psi_beta[test_index]
		U_S_expected_val = self.U_mu*self.S_U
		V_S_t_expected_val = self.V_mu[:,test_index]

		# Compute expectations on other components
		other_components_expected = np.zeros(self.N)
		for j in range(self.K):
			other_components_expected = other_components_expected + U_S_expected_val[:,j]*V_S_t_expected_val[j]
		resid = self.Y[:, test_index] - self.intercept_mu[test_index] - self.G[:, test_index]*(other_components_expected)

		# Loop through individuals
		for individual_index in range(self.I):
			# Indices of samples corresponding to this label
			sample_indices = self.individual_to_sample_indices[individual_index]
			# Number of indices corresponding to this individaul
			n_i = len(sample_indices)
			# Update variance of q(alpha_it)
			self.alpha_var[individual_index, test_index] = 1.0/(n_i*tau_t_expected_val + psi_t_expected_val)
			# Update mu of q(alpha_it)
			self.alpha_mu[individual_index, test_index] = self.alpha_var[individual_index, test_index]*tau_t_expected_val*np.sum(resid[sample_indices])
	def update_intercept(self):
		for test_index in range(self.T):
			self.update_intercept_t(test_index)
	def update_intercept_t(self, test_index):
		# Compute relevent expectations
		tau_t_expected_val = self.tau_alpha[test_index]/self.tau_beta[test_index]
		U_S_expected_val = self.U_mu*self.S_U
		V_S_t_expected_val = self.V_mu[:,test_index]

		# Compute expectations on other components
		other_components_expected = np.zeros(self.N)
		for j in range(self.K):
			other_components_expected = other_components_expected + U_S_expected_val[:,j]*V_S_t_expected_val[j]
		resid = self.Y[:, test_index] - self.alpha_big_mu[:,test_index] - self.G[:, test_index]*(other_components_expected)

		self.intercept_var[test_index] = 1.0/(self.N*tau_t_expected_val)
		self.intercept_mu[test_index] = self.intercept_var[test_index]*tau_t_expected_val*np.sum(resid)
	def update_psi(self):
		alpha_squared_expected_value = np.square(self.alpha_mu) + self.alpha_var
		# Loop through tests
		for test_index in range(self.T):
			self.psi_alpha[test_index] = self.alpha_prior + (self.I/2.0)
			self.psi_beta[test_index] = self.beta_prior + (np.sum(alpha_squared_expected_value[:,test_index])/2.0)
	def update_gamma_V(self):
		# Loop through factors
		for k in range(self.K):
			#V_squared_k_expected_val = ((np.square(self.V_mu[k,:]) + self.V_var[k,:])*self.S_V[k,:]) + (1.0-self.S_V[k,:])*(self.gamma_V_beta[k]/self.gamma_V_alpha[k])
			V_squared_k_expected_val = ((np.square(self.V_mu[k,:]) + self.V_var[k,:])*self.S_V[k,:]) + (1.0-self.S_V[k,:])*(self.V_var_s_0[k,:])
			self.gamma_V_alpha[k] = self.alpha_prior + (self.T/2.0)
			self.gamma_V_beta[k] = self.beta_prior + np.sum(V_squared_k_expected_val)/2.0
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
			self.theta_U_a[k] = self.a_prior + np.sum(self.S_U[:,k])
			self.theta_U_b[k] = self.b_prior + self.N - np.sum(self.S_U[:,k])
	def update_theta_V(self):
		# Loop through factors
		for k in range(self.K):
			self.theta_V_a[k] = self.a_prior + np.sum(self.S_V[k,:])
			self.theta_V_b[k] = self.b_prior + self.T - np.sum(self.S_V[k,:])
	def update_tau(self):
		# Other relevent expectations
		U_S = (self.U_mu)*(self.S_U)
		V_S = (self.V_mu)
		alpha_squared = np.square(self.alpha_big_mu) + self.alpha_big_var
		alpha = self.alpha_big_mu
		V_S_squared = ((np.square(self.V_mu) + self.V_var))
		U_S_squared = ((np.square(self.U_mu) + self.U_var)*self.S_U)
		intercept_squared = np.square(self.intercept_mu) + self.intercept_var
		intercept = self.intercept_mu

		componenent_squared_terms = np.dot(U_S_squared, V_S_squared)
		componenent_terms = np.dot(U_S, V_S)
		intercept_terms = np.dot(np.ones((self.N,1)),[intercept])
		intercept_squared_terms = np.dot(np.ones((self.N,1)),[intercept_squared])

		# Compute residual matrix
		#residual_mat = self.Y - self.G*(np.dot(U_S, V_S) + np.dot(np.ones((self.N,1)),[F_S])) - self.alpha_big_mu
		squared_residual_mat = np.square(self.Y) + intercept_squared_terms + alpha_squared + np.square(self.G)*componenent_squared_terms
		squared_residual_mat = squared_residual_mat - 2.0*self.Y*(intercept_terms + alpha + self.G*(componenent_terms))
		squared_residual_mat = squared_residual_mat + 2.0*intercept_terms*(alpha + self.G*(componenent_terms))
		squared_residual_mat = squared_residual_mat + 2.0*alpha*self.G*(componenent_terms)
		#squared_residual_mat = squared_residual_mat + 2.0*np.square(self.G)*componenent_terms*componenent_terms
		for k in range(self.K):
			for j in range(self.K):
				if k != j:
					squared_residual_mat = squared_residual_mat + np.square(self.G)*np.dot(np.transpose([U_S[:,k]]), [V_S[k,:]])*np.dot(np.transpose([U_S[:,j]]), [V_S[j,:]])
		self.tau_alpha =  self.tau_alpha*0.0 + self.alpha_prior + (self.N/2.0)
		self.tau_beta = self.beta_prior + (np.sum(squared_residual_mat,axis=0)/2.0)

	def update_elbo(self):
		data_likelihood_term = self.compute_elbo_log_likelihood_term()
		#kl_alpha = self.compute_kl_divergence_of_alpha()
		kl_tau = self.compute_kl_divergence_of_tau()
		#kl_psi = self.compute_kl_divergence_of_psi()
		kl_intercept = self.compute_kl_divergence_of_intercept()
		kl_V = self.compute_kl_divergence_of_V()
		kl_U_S = self.compute_kl_divergence_of_U_S()
		#kl_divergence = kl_V + kl_U_S + kl_alpha + kl_gamma_u + kl_tau + kl_psi + kl_theta_u + kl_intercept
		kl_divergence = kl_V + kl_U_S  + kl_tau + kl_intercept
		elbo = data_likelihood_term - kl_divergence
		self.elbo.append(elbo)
		if len(self.elbo) == 1:
			self.elbo_increase = True
		else:
			itera = len(self.elbo)-1
			if self.elbo[itera] < self.elbo[itera-1]:
				self.elbo_increase = False
				print(self.elbo[itera] - self.elbo[itera-1])
			else:
				self.elbo_increase = True

	def compute_kl_divergence_of_theta_f(self):
		a_prior = self.a_prior
		b_prior = self.b_prior
		theta_a = np.asarray([self.theta_F_a])
		theta_b = np.asarray([self.theta_F_b])
		kl_divergence = compute_kl_divergence_of_beta(a_prior, b_prior, theta_a, theta_b)
		return kl_divergence

	def compute_kl_divergence_of_theta_v(self):
		a_prior = self.a_prior
		b_prior = self.b_prior
		theta_a = self.theta_V_a 
		theta_b = self.theta_V_b
		kl_divergence = compute_kl_divergence_of_beta(a_prior, b_prior, theta_a, theta_b)
		return kl_divergence

	def compute_kl_divergence_of_theta_u(self):
		a_prior = self.a_prior
		b_prior = self.b_prior
		theta_a = self.theta_U_a 
		theta_b = self.theta_U_b
		kl_divergence = compute_kl_divergence_of_beta(a_prior, b_prior, theta_a, theta_b)
		return kl_divergence

	def compute_kl_divergence_of_psi(self):
		alpha_prior = self.alpha_prior
		beta_prior = self.beta_prior
		gamma_alpha = self.psi_alpha
		gamma_beta = self.psi_beta
		kl_divergence = compute_kl_divergence_of_gamma(alpha_prior, beta_prior, gamma_alpha, gamma_beta)
		return kl_divergence

	def compute_kl_divergence_of_tau(self):
		alpha_prior = self.alpha_prior
		beta_prior = self.beta_prior
		gamma_alpha = self.tau_alpha
		gamma_beta = self.tau_beta
		kl_divergence = compute_kl_divergence_of_gamma(alpha_prior, beta_prior, gamma_alpha, gamma_beta)
		return kl_divergence

	def compute_kl_divergence_of_gamma_f(self):
		alpha_prior = self.alpha_prior
		beta_prior = self.beta_prior
		gamma_alpha = np.asarray([self.gamma_F_alpha])
		gamma_beta = np.asarray([self.gamma_F_beta])
		kl_divergence = compute_kl_divergence_of_gamma(alpha_prior, beta_prior, gamma_alpha, gamma_beta)
		return kl_divergence

	def compute_kl_divergence_of_gamma_u(self):
		alpha_prior = self.alpha_prior
		beta_prior = self.beta_prior
		gamma_alpha = self.gamma_U_alpha
		gamma_beta = self.gamma_U_beta
		kl_divergence = compute_kl_divergence_of_gamma(alpha_prior, beta_prior, gamma_alpha, gamma_beta)
		return kl_divergence

	def compute_kl_divergence_of_gamma_v(self):
		alpha_prior = self.alpha_prior
		beta_prior = self.beta_prior
		gamma_alpha = self.gamma_V_alpha
		gamma_beta = self.gamma_V_beta
		kl_divergence = compute_kl_divergence_of_gamma(alpha_prior, beta_prior, gamma_alpha, gamma_beta)
		return kl_divergence

	def compute_kl_divergence_of_alpha(self):
		# Relevent expectations
		log_psi_expected = special.digamma(self.psi_alpha) - np.log(self.psi_beta)
		psi_expected = self.psi_alpha/self.psi_beta
		alpha_squared_expected_value = np.square(self.alpha_mu) + self.alpha_var

		likelihood_term_a = np.sum(log_psi_expected)*self.I/2.0
		likelihood_term_b = -np.sum(alpha_squared_expected_value*psi_expected)/2.0
		entropy_term_a = -np.sum(np.log(self.alpha_var))/2.0
		kl_divergence = entropy_term_a - likelihood_term_a - likelihood_term_b
		return kl_divergence

	def compute_kl_divergence_of_intercept(self):
		# Relevent expectations
		entropy_term_a = -np.sum(np.log(self.intercept_var))/2.0
		kl_divergence = entropy_term_a 
		return kl_divergence

	def compute_kl_divergence_of_F_S(self):
		S = np.asarray([self.S_F])
		W_mu = np.asarray([self.F_mu])
		W_var = np.asarray([self.F_var])
		W_var_s_0 = np.asarray([self.F_var_s_0])
		gamma_alpha = np.asarray([self.gamma_F_alpha])
		gamma_beta = np.asarray([self.gamma_F_beta])
		theta_a = np.asarray([self.theta_F_a])
		theta_b = np.asarray([self.theta_F_b])
		kl_divergence = compute_kl_divergence_of_gaussian_bernoulli(S, W_mu, W_var, W_var_s_0, gamma_alpha, gamma_beta, theta_a, theta_b, 1)
		return kl_divergence

	def compute_kl_divergence_of_V(self):
		W_mu = self.V_mu
		W_var = self.V_var
		expected_gamma_v = self.lambda_v
		kl_divergence = compute_kl_divergence_of_gaussian(W_mu, W_var, expected_gamma_v, self.K)
		return kl_divergence

	def compute_kl_divergence_of_U_S(self):
		S = np.transpose(self.S_U)
		W_mu = np.transpose(self.U_mu)
		W_var = np.transpose(self.U_var)
		W_var_s_0 = np.transpose(self.U_var_s_0)
		precision = self.lambda_u
		p = self.p_u
		kl_divergence = compute_kl_divergence_of_gaussian_bernoulli(S, W_mu, W_var, W_var_s_0,precision, p, self.K)
		return kl_divergence

	def compute_elbo_log_likelihood_term(self):
		# Compute expectation of log of gamma variables
		log_tau_expected = special.digamma(self.tau_alpha) - np.log(self.tau_beta)
		# Compute expectation of gamma variable
		tau_expected = self.tau_alpha/self.tau_beta
		# Other relevent expectations
		U_S = (self.U_mu)*(self.S_U)
		V_S = (self.V_mu)
		alpha_squared = np.square(self.alpha_big_mu) + self.alpha_big_var
		alpha = self.alpha_big_mu
		V_S_squared = ((np.square(self.V_mu) + self.V_var))
		U_S_squared = ((np.square(self.U_mu) + self.U_var)*self.S_U)
		intercept_squared = np.square(self.intercept_mu) + self.intercept_var
		intercept = self.intercept_mu

		componenent_squared_terms = np.dot(U_S_squared, V_S_squared)
		componenent_terms = np.dot(U_S, V_S)
		intercept_terms = np.dot(np.ones((self.N,1)),[intercept])
		intercept_squared_terms = np.dot(np.ones((self.N,1)),[intercept_squared])


		# Terms of interest in likelihood
		term_a = -np.log(2.0*np.pi)*(self.N*self.T/2.0)
		term_b = (self.N/2.0)*np.sum(log_tau_expected)
		# Compute residual matrix
		#residual_mat = self.Y - self.G*(np.dot(U_S, V_S) + np.dot(np.ones((self.N,1)),[F_S])) - self.alpha_big_mu
		squared_residual_mat = np.square(self.Y) + intercept_squared_terms + alpha_squared + np.square(self.G)*componenent_squared_terms
		squared_residual_mat = squared_residual_mat - 2.0*self.Y*(intercept_terms + alpha + self.G*(componenent_terms))
		squared_residual_mat = squared_residual_mat + 2.0*intercept_terms*(alpha + self.G*(componenent_terms))
		squared_residual_mat = squared_residual_mat + 2.0*alpha*self.G*(componenent_terms)
		#squared_residual_mat = squared_residual_mat + 2.0*np.square(self.G)*componenent_terms*componenent_terms
		for k in range(self.K):
			for j in range(self.K):
				if k != j:
					squared_residual_mat = squared_residual_mat + np.square(self.G)*np.dot(np.transpose([U_S[:,k]]), [V_S[k,:]])*np.dot(np.transpose([U_S[:,j]]), [V_S[j,:]])
		
		term_c = np.sum(squared_residual_mat*tau_expected)/2.0
		data_likelihood_term = term_a + term_b - term_c
		return data_likelihood_term
	def initialize_variables(self):
		# Initialize mapping from z-label to individual
		self.z_mapping = {}
		self.z_inverse_mapping = {}

		# Initialize array to keep track of ELBO
		self.elbo = []

		# Create mapping from grouping to index
		_, idx = np.unique(self.z, return_index=True)
		unique_groups = np.asarray(self.z)[np.sort(idx)]
		for i, label in enumerate(unique_groups):
			self.z_mapping[label] = i
			self.z_inverse_mapping[i] = label

		# Add model dimensions to object
		self.N = self.Y.shape[0]
		self.I = len(np.unique(self.z))
		self.T = self.Y.shape[1]

		self.individual_to_sample_indices = []
		# Create mapping from individual to sample indices
		for ii in range(self.I):
			# z_label corresponding to this individual
			z_label = self.z_inverse_mapping[ii]
			sample_indices = np.where(np.asarray(self.z) == z_label)[0]
			self.individual_to_sample_indices.append(sample_indices)

		# Initialize variational parameters
		# Us and Vs
		self.U_mu = np.random.randn(self.N, self.K)
		self.U_var = np.ones((self.N, self.K))
		self.U_var_s_0 = np.zeros((self.N, self.K)) + .001
		self.V_mu = np.random.randn(self.K, self.T)/10.0
		self.V_var = np.ones((self.K, self.T)) 
		# Masks
		self.S_U = np.ones((self.N,self.K))

		# Intercepts
		self.intercept_mu = np.zeros(self.T)
		self.intercept_var = np.ones(self.T) -.99
		# Random effects
		self.alpha_mu = np.zeros((self.I, self.T))
		self.alpha_var = np.ones((self.I, self.T))-.99
		# Convert random effects matrix to samplesXtests instead of groupsXtest
		self.alpha_big_mu = np.zeros((self.N, self.T))
		self.alpha_big_var = np.zeros((self.N, self.T))
		for sample_num, z_label in enumerate(self.z):
			self.alpha_big_mu[sample_num,:] = self.alpha_mu[self.z_mapping[z_label], :]
			self.alpha_big_var[sample_num,:] = self.alpha_var[self.z_mapping[z_label], :]
		# Variances
		self.psi_alpha = np.ones(self.T)
		self.psi_beta = np.ones(self.T)
		self.tau_alpha = np.ones(self.T)
		self.tau_beta = np.ones(self.T)
