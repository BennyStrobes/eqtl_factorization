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

# BB_GLM = pystan.StanModel(file = "betabinomial_glm.stan")
BB_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm.pkl', 'rb'))
# BB_GLM_FIXED_CONC = pystan.StanModel(file = "betabinomial_glm_fixed_conc.stan")
BB_GLM_FIXED_CONC = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm_fixed_conc.pkl', 'rb'))

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

def regress_out_cell_line(scaled_rat, Z_mat):
	# First make 1 hot encoding of Z
	'''
	num_cell_lines = len(np.unique(Z))
	num_cells = len(Z)
	Z_mat = np.zeros((num_cells, num_cell_lines-1))
	for n in range(num_cells):
		line_index = Z[n]
		if line_index != num_cell_lines-1:
			Z_mat[n, line_index] = 1.0
	'''
	final_rat = np.copy(scaled_rat)
	N, P = scaled_rat.shape
	for column_index in range(P):
		print(column_index)
		valid_row_indices = np.isnan(scaled_rat[:,column_index]) == False
		reg = LinearRegression().fit(Z_mat[valid_row_indices,:], scaled_rat[valid_row_indices,column_index])
		predicted = reg.predict(Z_mat[valid_row_indices,:])
		final_rat[valid_row_indices,column_index] = scaled_rat[valid_row_indices,column_index] - predicted
	return final_rat

def outside_update_V_t(U_new, ys, ns, concShape, concRate, itera, conc_init, C_init, V_init):
	observed_indices = np.isnan(ys) == False
	data = dict(N=sum(observed_indices), P=U_new.shape[1], x=U_new[observed_indices,:], ys=ys[observed_indices].astype(int), ns=ns[observed_indices].astype(int), concShape=concShape, concRate=concRate)
	########################################
	# Get MOM estimates for initialization
	if itera == 0:
		rat = data['ys']/data['ns'].astype(float)
		if np.sum(np.isnan(rat)) > 0:
			pdb.set_trace()
		# moment estimator of concentration parameter
		mom_conc_init = min(1.0/np.var(rat), 1000.0)
		# thresholded moment estimator of the mean 
		m_init = min(max(np.mean(rat), 1.0/1000 ), 1.0-(1.0/1000))
		beta_init = np.zeros(data['P'])
		beta_init[0] = np.log(m_init/(1.0-m_init))
		init = dict(conc=mom_conc_init, beta=beta_init)
	else:
		init = dict(conc=conc_init, beta=np.hstack((C_init,V_init)))
	########################################
	# Run optimization
	try:
		with suppress_stdout_stderr():
			op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1, init=init)
	except RuntimeError:
		try: 
			op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1)
		except RuntimeError:
			try:
				op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1, algorithm="Newton")
			except RuntimeError:
				try:
					op = BB_GLM.optimizing(data = data, verbose=False,iter=5000,seed=2, algorithm="Newton")
				except RuntimeError:
					print('error')
					pdb.set_trace()
	return np.hstack((op['conc'],op['beta']))


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
	def __init__(self, K=5, concShape=1.0001, concRate=1.0001, max_iter=1000, output_root='temp'):
		self.max_iter = max_iter
		self.K = K
		self.concShape = concShape
		self.concRate = concRate
		self.iter = 0
		self.output_root = output_root + '_4'
	def fit(self, allelic_counts, total_counts, cov):
		""" Fit the model.
			Args:
		"""
		np.random.seed(1)
		self.allelic_counts = allelic_counts
		self.total_counts = total_counts
		self.cov = cov
		self.initialize_variables()
		for vi_iter in range(self.max_iter):
			print('ITERATION ' + str(self.iter))
			self.update_V_and_conc()
			self.update_U()
			# Save to output every five iters
			if np.mod(vi_iter, 5) == 0: 
				np.savetxt(self.output_root + '_iter.txt', np.asmatrix(vi_iter), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_U.txt', (self.U), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_V.txt', (self.V), fmt="%s", delimiter='\t')
				np.savetxt(self.output_root + '_temper_conc.txt', (self.conc), fmt="%s", delimiter='\t')
			self.iter = self.iter + 1

	def update_U(self):
		covariate_predicted = np.dot(self.cov, self.C)
		for sample_iter in range(self.N):
			observed_indices = np.isnan(self.allelic_counts[sample_iter,:]) == False
			merged_intercept = covariate_predicted[sample_iter, :]
			data = dict(N=sum(observed_indices), P=self.U.shape[1], x=np.transpose(self.V[:,observed_indices]), intercept=merged_intercept[observed_indices], ys=self.allelic_counts[sample_iter, observed_indices].astype(int), ns=self.total_counts[sample_iter, observed_indices].astype(int), conc=self.conc[observed_indices])
			########################################
			# Initialize
			if self.iter == 0:
				beta_init = np.zeros(data['P'])
			else:
				beta_init = self.U[sample_iter, :]
			init = dict(beta=beta_init)
			########################################
			# Run optimization
			try:
				with suppress_stdout_stderr():
					op = BB_GLM_FIXED_CONC.optimizing(data = data, verbose=False,iter=5000,seed=1, init=init)
			except RuntimeError:
				try:
					op = BB_GLM_FIXED_CONC.optimizing(data = data, verbose=False,iter=5000,seed=1)
				except RuntimeError:
					try:
						op = BB_GLM_FIXED_CONC.optimizing(data = data, verbose=False,iter=5000,seed=1, algorithm="Newton")
					except RuntimeError:
						try:
							op = BB_GLM_FIXED_CONC.optimizing(data = data, verbose=False,iter=5000,seed=2, algorithm="Newton")
						except RuntimeError:
							print('eror')
							pdb.set_trace()
			self.U[sample_iter, :] = op['beta']

	def update_V_and_conc(self):
		# Add column of ones (intercept to U)
		#U_new = np.hstack((X0, self.U, self.cov))
		parrallel = True
		V_update_data = []
		if parrallel == False:
			U_new = np.hstack((self.cov, self.U))
			for test_index in range(self.T):
				print(test_index)
				V_update_data.append(outside_update_V_t(U_new, self.allelic_counts[:, test_index], self.total_counts[:, test_index], self.concShape, self.concRate, self.iter, self.conc[test_index], self.C[:, test_index], self.V[:, test_index]))
		elif parrallel == True:
			U_new = np.hstack((self.cov, self.U))
			V_update_data = Parallel(n_jobs=24)(delayed(outside_update_V_t)(U_new, self.allelic_counts[:, test_index], self.total_counts[:, test_index], self.concShape, self.concRate, self.iter, self.conc[test_index], self.C[:, test_index], self.V[:, test_index]) for test_index in range(self.T))
		V_update_data = np.asarray(V_update_data)
		V_update_data_2 = np.transpose(V_update_data[:,1:])
		self.conc = V_update_data[:,0]
		self.C = V_update_data_2[:(self.num_cov), :]
		self.V = V_update_data_2[(self.num_cov):, :]

	def initialize_variables(self):
		self.N = self.allelic_counts.shape[0]
		self.T = self.allelic_counts.shape[1]
		self.num_cov = self.cov.shape[1]
		# Randomly initialize (only U matters)
		self.U = np.random.randn(self.N, self.K)
		self.V = np.random.randn(self.K, self.T) 
		self.conc = np.random.randn(self.T)
		self.C = np.random.randn(self.num_cov, self.T)
		ppca_init = True
		if ppca_init == True:
			rat = self.allelic_counts/self.total_counts
			nans = np.isnan(rat)
			scaled_rat = scale_allelic_ratios(rat)
			scaled_residual_rat = regress_out_cell_line(scaled_rat, self.cov[:,1:])
			rescaled_residual_rat = scale_allelic_ratios(scaled_residual_rat)
			ppca = PPCA()
			ppca.fit(data=np.transpose(rescaled_residual_rat), d=self.K, verbose=True, tol=1e-6)
			self.U = ppca.C/np.std(ppca.C)
