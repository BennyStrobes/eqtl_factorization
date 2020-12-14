import numpy as np 
import os
import sys
import pdb
import pystan
import pickle
from joblib import Parallel, delayed
import multiprocessing
from ppca import PPCA
from sklearn.linear_model import LinearRegression
import scipy.special

# BB_GLM = pystan.StanModel(file = "betabinomial_glm.stan")
#BB_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm.pkl', 'rb'))
# BB_GLM_FIXED_CONC = pystan.StanModel(file = "betabinomial_glm_fixed_conc.stan")
#BB_GLM_FIXED_CONC = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm_fixed_conc.pkl', 'rb'))

BB_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/soft_weighted_folded_beta_binomial_glm.pkl', 'rb'))
BB_GLM_INTERCEPT = pickle.load(open('/work-zfs/abattle4/bstrober/temp/soft_weighted_folded_beta_binomial_w_intercept_glm.pkl', 'rb'))
BB_FIXED_CONC_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/soft_weighted_folded_beta_binomial_fixed_conc_glm.pkl', 'rb'))
SKEWED_BB = pickle.load(open('/work-zfs/abattle4/bstrober/temp/soft_weighted_skewed_beta_binomial.pkl', 'rb'))


def beta_binomial_log_p(n, k, alpha, beta):
	logp = binomln(n, k) + betaln(k + alpha, n - k + beta) - betaln(alpha, beta)
	return logp
def factln(n):
    return scipy.special.gammaln(n + 1)
def betaln(x, y):
    return scipy.special.gammaln(x) + scipy.special.gammaln(y) - scipy.special.gammaln(x + y)
def binomln(n, k):
    return factln(n) - factln(k) - factln(n - k)

def sum_over_rows(matty):
	nrow, ncol = matty.shape
	summation = np.squeeze(np.asarray(np.nansum(matty, axis=0)))
	true_nans = np.squeeze(np.asarray(np.sum(np.isnan(matty),axis=0) == nrow))
	summation[true_nans] = np.nan
	return summation

def outside_update_V_conc_C_t(U_new, ys, ns, phi, z, concShape, concRate, cauchy_scale):
	########################################
	# load in your data
	########################################
	observed_indices = np.isnan(ys) == False
	phi_full = phi[z,:]
	# quick error checking
	if np.sum(np.isnan(phi_full[:,0])) > np.sum(np.isnan(ys)):
		print('assumption error')
		pdb.set_trace()
	if np.sum(np.isnan(phi_full[observed_indices,0])) > 0:
		print('assumption error')
		pdb.set_trace()
	########################################
	# Part 1: Update V and Conc
	## Run optimization with pystan
	########################################
	random_subset = False
	if random_subset == False:
		data = dict(N=sum(observed_indices), P=U_new.shape[1], x=U_new[observed_indices,:], ys=ys[observed_indices].astype(int), ns=ns[observed_indices].astype(int), prob1=phi_full[observed_indices,0], prob2=phi_full[observed_indices,1], concShape=concShape, concRate=concRate)
	if random_subset == True:
		subset_size = 1000
		if sum(observed_indices) < subset_size:
			observed_subsetted_indices = np.copy(observed_indices)
		else:
			indexed_observed_indices = np.where(observed_indices ==True)[0]
			subset_indexes = np.random.choice(indexed_observed_indices, size=len(indexed_observed_indices) - subset_size, replace=False)
			observed_subsetted_indices = np.copy(observed_indices)
			observed_subsetted_indices[subset_indexes] = False
			data = dict(N=sum(observed_subsetted_indices), P=U_new.shape[1], x=U_new[observed_subsetted_indices,:], ys=ys[observed_subsetted_indices].astype(int), ns=ns[observed_subsetted_indices].astype(int), prob1=phi_full[observed_subsetted_indices,0], prob2=phi_full[observed_subsetted_indices,1], concShape=concShape, concRate=concRate)
	try:
		with suppress_stdout_stderr():
			op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1)
	except RuntimeError:
		try: 
			op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=2)
		except RuntimeError:
			try:
				op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1, algorithm="Newton")
			except RuntimeError:
				try:
					op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=2, algorithm="Newton")
				except RuntimeError:
					print('error')
					pdb.set_trace()
	conc_new = op['conc']
	beta_new = op['beta']
	########################################
	# Part 2: Update alpha_a1
	## Run optimization with pystan
	########################################
	data2 = dict(N=sum(observed_indices), ys=ys[observed_indices].astype(int), ns=ns[observed_indices].astype(int), prob1=phi_full[observed_indices,2], cauchy_scale=cauchy_scale)
	try:
		with suppress_stdout_stderr():
			op2 = SKEWED_BB.optimizing(data = data2, verbose=False,iter=5000,seed=1)
	except RuntimeError:
		try: 
			op2 = SKEWED_BB.optimizing(data = data2, verbose=False,iter=5000,seed=2)
		except RuntimeError:
			try:
				op2 = SKEWED_BB.optimizing(data = data2, verbose=False,iter=5000,seed=1, algorithm="Newton")
			except RuntimeError:
				try:
					op2 = SKEWED_BB.optimizing(data = data2, verbose=False,iter=5000,seed=2, algorithm="Newton")
				except RuntimeError:
					print('error')
					pdb.set_trace()
	a1_new = op2['a1']
	########################################
	# Part 3: Update alpha_a2
	## Run optimization with pystan
	########################################	
	data3 = dict(N=sum(observed_indices), ys=ns[observed_indices].astype(int) - ys[observed_indices].astype(int), ns=ns[observed_indices].astype(int), prob1=phi_full[observed_indices,3], cauchy_scale=cauchy_scale)
	try:
		with suppress_stdout_stderr():
			op3 = SKEWED_BB.optimizing(data = data3, verbose=False,iter=5000,seed=1)
	except RuntimeError:
		try: 
			op3 = SKEWED_BB.optimizing(data = data3, verbose=False,iter=5000,seed=2)
		except RuntimeError:
			try:
				op3 = SKEWED_BB.optimizing(data = data3, verbose=False,iter=5000,seed=1, algorithm="Newton")
			except RuntimeError:
				try:
					op3 = SKEWED_BB.optimizing(data = data3, verbose=False,iter=5000,seed=2, algorithm="Newton")
				except RuntimeError:
					print('error')
					pdb.set_trace()
	a2_new = op3['a1']
	return np.hstack((conc_new, a1_new, a2_new, beta_new))

def outside_update_U_n(V_new, ys, ns, conc, covariate_predicted, phi_full):
	observed_indices = np.isnan(ys) == False

	data = dict(N=sum(observed_indices), P=V_new.shape[0], x=np.transpose(V_new[:,observed_indices]), intercept=covariate_predicted[observed_indices], ys=ys[observed_indices].astype(int), ns=ns[observed_indices].astype(int), conc=conc[observed_indices], prob1=phi_full[observed_indices,0], prob2=phi_full[observed_indices,1])
	########################################
	# Run optimization
	try:
		with suppress_stdout_stderr():
			op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1)
	except RuntimeError:
		try:
			op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=2)
		except RuntimeError:
			try:
				op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=10, algorithm="BFGS")
			except RuntimeError:
				try:
					op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=11, algorithm="BFGS")
				except RuntimeError:
					try:
						op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=3)
					except RuntimeError:
						try:
							op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=4)
						except RuntimeError:
							print('eror')
							pdb.set_trace()
	return op['beta']

def outside_update_V_t(U_new, ys, ns, phi, z, covariate_predicted, concShape, concRate):
	########################################
	# load in your data
	########################################
	observed_indices = np.isnan(ys) == False
	phi_full = phi[z,:]
	# quick error checking
	if np.sum(np.isnan(phi_full[:,0])) > np.sum(np.isnan(ys)):
		print('assumption error')
		pdb.set_trace()
	if np.sum(np.isnan(phi_full[observed_indices,0])) > 0:
		print('assumption error')
		pdb.set_trace()
	########################################
	# Part 1: Update V and Conc
	## Run optimization with pystan
	########################################
	random_subset = False

	data = dict(N=sum(observed_indices), P=U_new.shape[1], x=U_new[observed_indices,:], ys=ys[observed_indices].astype(int), ns=ns[observed_indices].astype(int), prob1=phi_full[observed_indices,0], prob2=phi_full[observed_indices,1], intercept=covariate_predicted[observed_indices], concShape=concShape, concRate=concRate)

	try:
		with suppress_stdout_stderr():
			op = BB_GLM_INTERCEPT.optimizing(data = data, verbose=False,iter=5000,seed=1)
	except RuntimeError:
		try: 
			op = BB_GLM_INTERCEPT.optimizing(data = data, verbose=False,iter=5000,seed=2)
		except RuntimeError:
			try:
				op = BB_GLM_INTERCEPT.optimizing(data = data, verbose=False,iter=5000,seed=1, algorithm="Newton")
			except RuntimeError:
				try:
					op = BB_GLM_INTERCEPT.optimizing(data = data, verbose=False,iter=5000,seed=2, algorithm="Newton")
				except RuntimeError:
					print('error')
					pdb.set_trace()
	beta_new = op['beta']
	conc_new = op['conc']
	return np.hstack((conc_new, beta_new))
	

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


class ASE_FACTORIZATION(object):
	def __init__(self, K=5, concShape=1.0001, concRate=1e-4, cauchy_scale=2, max_iter=1000, output_root='temp', random_seed=1):
		self.max_iter = max_iter
		self.K = K
		self.concShape = concShape
		self.concRate = concRate
		self.cauchy_scale = cauchy_scale
		self.iter = 0
		self.random_seed = random_seed
		self.output_root = output_root + '_no_ppca_seed_' + str(self.random_seed) + '_'
	def filter_out_outlier_genes(self):
		num_genes = self.allelic_counts.shape[1]
		valid_genes = []
		for gene_num in range(num_genes):
			rat = self.allelic_counts[:, gene_num]/self.total_counts[:,gene_num]
			nans = np.isnan(rat)
			rat = rat[~nans]
			if 2.0*(sum(rat < .1) + sum(rat>.9)) > len(rat):
				continue
			valid_genes.append(gene_num)
		valid_genes = np.asarray(valid_genes)
		self.allelic_counts = self.allelic_counts[:, valid_genes]
		self.total_counts = self.total_counts[:,valid_genes]
		print(self.total_counts.shape)
	def fit(self, allelic_counts, total_counts, cov, z):
		""" Fit the model.
			Args:
		"""
		np.random.seed(self.random_seed)
		self.allelic_counts = allelic_counts
		self.total_counts = total_counts
		self.filter_out_outlier_genes()
		self.cov = cov
		self.z = z.astype(int)
		self.observed_data_log_likelihoods = []
		self.initialize_variables()
		for vi_iter in range(self.max_iter):
			print('ITERATION ' + str(self.iter))
			
			observed_data_log_likelihood = self.compute_observed_data_log_likelihood()
			print('Observed data log like: ' + str(observed_data_log_likelihood))
			self.observed_data_log_likelihoods.append(observed_data_log_likelihood)
			if vi_iter > 0:
				diff = observed_data_log_likelihood - self.observed_data_log_likelihoods[-2]
				print('Delta log like: ' + str(diff))
			# E-STEP: Loop through individuals
			print('Starting E-STEP')
			self.e_step()
			print('Starting M-STEP')
			self.update_w()
			self.update_V_and_conc_and_C()
			self.update_U()
			for inner_iter in range(15):
				self.update_V()
				self.update_U()
			# Save to output every five iters
			if np.mod(vi_iter, 1) == 0: 
				np.savetxt(self.output_root + '_iter.txt', np.asmatrix(vi_iter), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root +'_temper_U.txt', (self.U), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_V.txt', (self.V), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_conc.txt', (self.conc), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_w.txt', (self.w), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_C.txt', (self.C), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_alpha_a1.txt', (self.alpha_a1), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_alpha_a2.txt', (self.alpha_a2), fmt="%s", delimiter='\t')
				pickle.dump(self.phi, open(self.output_root + "temper_phi.pkl", "wb" ) )
			self.iter = self.iter + 1
	def update_w(self):
		# Initialize matrix to keep track of w for this iteration
		posterior_sums = np.zeros((self.T, 4)) + 5.0
		# Loop through individuals
		for individual_index in range(self.I):
			for test_index in range(self.T):
				if np.isnan(self.phi[individual_index][test_index, 0]):
					continue
				posterior_sums[test_index,:] = posterior_sums[test_index, :] + self.phi[individual_index][test_index,:]
		for test_index in range(self.T):
			self.w[test_index,:] = posterior_sums[test_index,:]/np.sum(posterior_sums[test_index,:])
	def compute_observed_data_log_likelihood(self):
		# Initialize variable
		observed_data_log_likelihood = 0
		# Pre-compute matrix level quantities
		predicted_p_mat = scipy.special.expit(np.dot(self.U, self.V) + np.dot(self.cov, self.C))
		conc_mat = np.dot(np.ones((self.N, 1)), np.asmatrix(self.conc))
		alpha_1_mat = np.dot(np.ones((self.N, 1)), np.asmatrix(self.alpha_a1))
		alpha_2_mat = np.dot(np.ones((self.N, 1)), np.asmatrix(self.alpha_a2))
		# Compute logp across all samples for each mixture component
		log_p_1 = beta_binomial_log_p(self.total_counts, self.allelic_counts, np.multiply(predicted_p_mat, conc_mat), np.multiply((1.0-predicted_p_mat), conc_mat))
		log_p_2 = beta_binomial_log_p(self.total_counts, self.allelic_counts, np.multiply((1.0-predicted_p_mat), conc_mat), np.multiply(predicted_p_mat, conc_mat))
		log_p_3 = beta_binomial_log_p(self.total_counts, self.allelic_counts, np.ones((self.N, self.T)), 7.0 + alpha_1_mat)
		log_p_4 = beta_binomial_log_p(self.total_counts, self.allelic_counts, 7.0 + alpha_2_mat, np.ones((self.N, self.T)))
		# Simple error checking to make sure no nans came about
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(log_p_2)):
			print('assumption error')
			pdb.set_trace()
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(log_p_3)):
			print('assumption error')
			pdb.set_trace()
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(log_p_4)):
			print('assumption error')
			pdb.set_trace()
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(self.total_counts)):
			print('assumption error')
			pdb.set_trace()
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(self.allelic_counts)):
			print('assumption error')
			pdb.set_trace()
		# Loop through individuals
		for individual_index in range(self.I):
			# Get indexes of samples corresponding to this individual
			sample_indices = np.where(self.z==individual_index)[0]
			num_samples = len(sample_indices)
			if num_samples == 0:
				print('assumption error: No samples for individual ' + str(individual_index))
				pdb.set_trace()
			# Temporarily replace phi with log probabilities
			indi_log_p_0 = sum_over_rows(log_p_1[sample_indices,:]) + np.log(self.w[:,0])
			indi_log_p_1 = sum_over_rows(log_p_2[sample_indices,:]) + np.log(self.w[:,1])
			indi_log_p_2 = sum_over_rows(log_p_3[sample_indices,:]) + np.log(self.w[:,2])
			indi_log_p_3 = sum_over_rows(log_p_4[sample_indices,:]) + np.log(self.w[:,3])

			for test_index in range(self.T):
				arr = [indi_log_p_0[test_index], indi_log_p_1[test_index], indi_log_p_2[test_index], indi_log_p_3[test_index]]
				indi_test_log_likelihood = scipy.special.logsumexp(arr)
				if np.isnan(indi_test_log_likelihood):
					continue
				observed_data_log_likelihood = observed_data_log_likelihood + indi_test_log_likelihood
		return observed_data_log_likelihood

	def e_step(self):
		# Pre-compute matrix level quantities
		predicted_p_mat = scipy.special.expit(np.dot(self.U, self.V) + np.dot(self.cov, self.C))
		conc_mat = np.dot(np.ones((self.N, 1)), np.asmatrix(self.conc))
		alpha_1_mat = np.dot(np.ones((self.N, 1)), np.asmatrix(self.alpha_a1))
		alpha_2_mat = np.dot(np.ones((self.N, 1)), np.asmatrix(self.alpha_a2))
		# Compute logp across all samples for each mixture component
		log_p_1 = beta_binomial_log_p(self.total_counts, self.allelic_counts, np.multiply(predicted_p_mat, conc_mat), np.multiply((1.0-predicted_p_mat), conc_mat))
		log_p_2 = beta_binomial_log_p(self.total_counts, self.allelic_counts, np.multiply((1.0-predicted_p_mat), conc_mat), np.multiply(predicted_p_mat, conc_mat))
		log_p_3 = beta_binomial_log_p(self.total_counts, self.allelic_counts, np.ones((self.N, self.T)), 7.0 + alpha_1_mat)
		log_p_4 = beta_binomial_log_p(self.total_counts, self.allelic_counts, 7.0 + alpha_2_mat, np.ones((self.N, self.T)))
		# Simple error checking to make sure no nans came about
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(log_p_2)):
			print('assumption error')
			pdb.set_trace()
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(log_p_3)):
			print('assumption error')
			pdb.set_trace()
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(log_p_4)):
			print('assumption error')
			pdb.set_trace()
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(self.total_counts)):
			print('assumption error')
			pdb.set_trace()
		if np.sum(np.isnan(log_p_1)) != np.sum(np.isnan(self.allelic_counts)):
			print('assumption error')
			pdb.set_trace()
		# Loop through individuals
		for individual_index in range(self.I):
			# Get indexes of samples corresponding to this individual
			sample_indices = np.where(self.z==individual_index)[0]
			num_samples = len(sample_indices)
			if num_samples == 0:
				print('assumption error: No samples for individual ' + str(individual_index))
				pdb.set_trace()
			# Temporarily replace phi with log probabilities
			self.phi[individual_index][:,0] = sum_over_rows(log_p_1[sample_indices,:]) + np.log(self.w[:,0])
			self.phi[individual_index][:,1] = sum_over_rows(log_p_2[sample_indices,:]) + np.log(self.w[:,1])
			self.phi[individual_index][:,2] = sum_over_rows(log_p_3[sample_indices,:]) + np.log(self.w[:,2])
			self.phi[individual_index][:,3] = sum_over_rows(log_p_4[sample_indices,:]) + np.log(self.w[:,3])
			for test_index in range(self.T):
				# Ignore nan rows
				if np.sum(np.isnan(self.phi[individual_index][test_index,:])) > 0:
					# Quick error checking
					if np.sum(np.isnan(self.phi[individual_index][test_index,:])) != 4:
						print('assumption erorr')
						pdb.set_trace()
					if sum(~np.isnan(self.allelic_counts[sample_indices,test_index])) > 0:
						print('assumption error')
						pdb.set_trace()
					if sum(~np.isnan(self.total_counts[sample_indices,test_index])) > 0:
						print('assumption error')
						pdb.set_trace()
					continue
				# Some more quick error checking 
				if sum(~np.isnan(self.allelic_counts[sample_indices,test_index])) == 0:
					print('assumption error')
					pdb.set_trace()
				if sum(~np.isnan(self.total_counts[sample_indices,test_index])) == 0:
					print('assumption error')
					pdb.set_trace()	
				if sum(~np.isnan(self.total_counts[sample_indices,test_index])) != sum(~np.isnan(self.total_counts[sample_indices,test_index])):
					print('assumption error')
					pdb.set_trace()
				# log sum exp trick to compute expectations
				self.phi[individual_index][test_index,:] = np.exp(self.phi[individual_index][test_index,:] - max(self.phi[individual_index][test_index,:]))
				self.phi[individual_index][test_index,:] = self.phi[individual_index][test_index,:]/np.sum(self.phi[individual_index][test_index,:])
	def update_U(self):
		covariate_predicted = np.dot(self.cov, self.C)
		parrallel = True
		U_update_data = []
		if parrallel == False:
			for sample_iter in range(self.N):
				U_update_data.append(outside_update_U_n(self.V, self.allelic_counts[sample_iter, :], self.total_counts[sample_iter, :], self.conc, covariate_predicted[sample_iter, :],  self.phi[self.z[sample_iter],:,:]))
		elif parrallel == True:
			U_update_data = Parallel(n_jobs=24)(delayed(outside_update_U_n)(self.V, self.allelic_counts[sample_iter, :], self.total_counts[sample_iter, :], self.conc, covariate_predicted[sample_iter, :],  self.phi[self.z[sample_iter],:,:]) for sample_iter in range(self.N))
		U_update_data = np.asarray(U_update_data)
		self.U = U_update_data

	def update_U_old(self):
		covariate_predicted = np.dot(self.cov, self.C)

		for sample_iter in range(self.N):
			individual_index = self.z[sample_iter]
			phi_full = self.phi[individual_index,:,:]
			observed_indices = np.isnan(self.allelic_counts[sample_iter,:]) == False
			# Simple error checking
			if np.sum(np.isnan(phi_full[observed_indices,:])) > 0:
				print('assumption error')
				pdb.set_trace()
			merged_intercept = covariate_predicted[sample_iter, :]
			data = dict(N=sum(observed_indices), P=self.U.shape[1], x=np.transpose(self.V[:,observed_indices]), intercept=merged_intercept[observed_indices], ys=self.allelic_counts[sample_iter, observed_indices].astype(int), ns=self.total_counts[sample_iter, observed_indices].astype(int), conc=self.conc[observed_indices], prob1=phi_full[observed_indices,0], prob2=phi_full[observed_indices,1])
			########################################
			# Run optimization
			try:
				with suppress_stdout_stderr():
					op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1)
			except RuntimeError:
				try:
					op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=2)
				except RuntimeError:
					try:
						op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=10, algorithm="BFGS")
					except RuntimeError:
						try:
							op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=11, algorithm="BFGS")
						except RuntimeError:
							try:
								op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=3)
							except RuntimeError:
								try:
									op = BB_FIXED_CONC_GLM.optimizing(data = data, verbose=False,iter=5000,seed=4)
								except RuntimeError:
									print('eror')
									pdb.set_trace()
			self.U[sample_iter, :] = op['beta']

	def update_V_and_conc_and_C(self):
		# Add column of ones (intercept to U)
		#U_new = np.hstack((X0, self.U, self.cov))
		parrallel = True
		V_update_data = []
		if parrallel == False:
			U_new = np.hstack((self.cov, self.U))
			for test_index in range(self.T):
				print(test_index)
				V_update_data.append(outside_update_V_conc_C_t(U_new, self.allelic_counts[:, test_index], self.total_counts[:, test_index], self.phi[:, test_index, :], self.z, self.concShape, self.concRate, self.cauchy_scale))
		elif parrallel == True:
			U_new = np.hstack((self.cov, self.U))
			V_update_data = Parallel(n_jobs=24)(delayed(outside_update_V_conc_C_t)(U_new, self.allelic_counts[:, test_index], self.total_counts[:, test_index], self.phi[:, test_index, :], self.z, self.concShape, self.concRate, self.cauchy_scale) for test_index in range(self.T))
		V_update_data = np.asarray(V_update_data)
		V_update_data_2 = np.transpose(V_update_data[:,3:])
		self.conc = V_update_data[:,0]
		self.alpha_a1 = V_update_data[:,1]
		self.alpha_a2 = V_update_data[:,2]
		self.C = V_update_data_2[:(self.num_cov), :]
		self.V = V_update_data_2[(self.num_cov):, :]
	def update_V(self):
		# Add column of ones (intercept to U)
		#U_new = np.hstack((X0, self.U, self.cov))
		covariate_predicted = np.dot(self.cov, self.C)
		parrallel = True
		#if self.iter > 0:
		#	parrallel = False
		V_update_data = []
		if parrallel == False:
			for test_index in range(self.T):
				print(test_index)
				V_update_data.append(outside_update_V_t(self.U, self.allelic_counts[:, test_index], self.total_counts[:, test_index], self.phi[:, test_index, :], self.z, covariate_predicted[:, test_index], self.concShape, self.concRate))
		elif parrallel == True:
			U_new = np.hstack((self.cov, self.U))
			V_update_data = Parallel(n_jobs=24)(delayed(outside_update_V_t)(self.U, self.allelic_counts[:, test_index], self.total_counts[:, test_index], self.phi[:, test_index, :], self.z, covariate_predicted[:, test_index], self.concShape, self.concRate) for test_index in range(self.T))
		V_update_data = np.asarray(V_update_data)
		self.conc = V_update_data[:,0]
		self.V = np.transpose(V_update_data[:,1:])

	def initialize_variables(self):
		self.I = len(np.unique(self.z))
		self.N = self.allelic_counts.shape[0]
		self.T = self.allelic_counts.shape[1]
		self.num_cov = self.cov.shape[1]
		# Randomly initialize 
		self.U = np.random.randn(self.N, self.K)*.01
		self.w = np.ones((self.T, 4))
		self.w[:,0] = self.w[:,0]*.49
		self.w[:,1] = self.w[:,1]*.49
		self.w[:,2] = self.w[:,2]*.01
		self.w[:,3] = self.w[:,3]*.01
		self.alpha_a1 = np.ones(self.T)
		self.alpha_a2 = np.ones(self.T)
		self.V = np.random.randn(self.K, self.T)*.01
		# Decently smart initialization of conc and C
		self.conc = np.ones(self.T)*10.0
		self.C = np.zeros((self.num_cov, self.T))
		for test_index in range(self.T):
			rats = self.allelic_counts[:,test_index]/self.total_counts[:,test_index]
			nans = np.isnan(rats)
			rats = rats[~nans]
			conc_init = np.min((1.0/np.var(rats), 1000))
			intercept_init = np.min((np.max((np.mean(rats), 1/1000.0)), (1.0 - 1.0/1000)))
			intercept_init = np.log(intercept_init/(1.0-intercept_init))
			self.C[0, test_index] = intercept_init
			self.conc[test_index] = conc_init
		# E-step stuff
		self.phi = np.ones((self.I, self.T, 4))
		start_over = False
		if start_over == True:
			root = '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/eqtl_results/ase_factorization_via_em_als_folded_beta_binomial_subsampled_high_biallelic_fraction_only_endoderm_differentiation_3_ase_factorization_no_ppca_seed_4_'
			self.U = np.loadtxt(root + '_temper_U.txt')
			self.w = np.loadtxt(root + '_temper_w.txt')
			self.alpha_a1 = np.loadtxt(root + '_temper_alpha_a1.txt')
			self.alpha_a2 = np.loadtxt(root + '_temper_alpha_a2.txt')
			self.V = np.loadtxt(root + '_temper_V.txt')
			self.conc = np.loadtxt(root + '_temper_conc.txt')
			self.C = np.loadtxt(root + '_temper_C.txt')
			self.phi = pickle.load( open( root + "temper_phi.pkl", "rb" ) )
			self.iter = np.loadtxt(root + '_iter.txt') + 0

