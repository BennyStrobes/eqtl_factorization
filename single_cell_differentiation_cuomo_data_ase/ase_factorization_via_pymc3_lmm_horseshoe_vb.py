import numpy as np 
import os
import sys
import pdb
import pickle
import pymc3 as pm


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
		self.output_root = output_root + '1'
	def fit(self, allelic_counts, total_counts, cov, z):
		""" Fit the model.
			Args:
		"""
		self.allelic_counts = allelic_counts
		self.total_counts = total_counts
		self.cov = cov
		self.Z = z.astype(int)
		self.initialize_variables()
		N = self.allelic_counts.shape[0]
		T = self.allelic_counts.shape[1]
		I = len(np.unique(self.Z))

		self.run_factorization(N, T, self.cov, self.Z, I, self.K, self.num_cov, self.allelic_counts, self.total_counts)
		
		#data = dict(N=N, T=T, K=self.K, num_cov=self.num_cov, cov=self.cov, ys=np.transpose(self.allelic_counts.astype(int)), ns=np.transpose(self.total_counts.astype(int)), concShape=self.concShape, concRate=self.concRate)
		#print("START")
		#aa = BB_FACTORIZATION.vb(data=data)
		#pickle.dump(aa, open(self.output_root + '_model', 'wb'))

	def run_factorization(self, N, S, X, Z, I, K, num_cov, k, n):
		# Smart initialization
		print("STARTING")
		rat = k/n
		nans = np.isnan(rat)
		conc_inits = np.zeros((1, S))
		beta_inits = np.zeros((num_cov, S))
		for index_s in range(S):
			column_rat = rat[:, index_s]
			column_nans = np.isnan(column_rat)
			valid_rat = column_rat[~column_nans]
			conc_init = min(1.0/np.var(valid_rat), 1000.0)
			m_init = min(max(np.mean(valid_rat), 1.0/1000 ), 1.0-(1.0/1000))
			conc_inits[0, index_s] = conc_init
			beta_inits[0, index_s] = np.log(m_init/(1.0-m_init))
		# Run bb-mf
		with pm.Model() as bb_glm:
			CONC = pm.HalfCauchy('CONC', beta=5, shape=(1,S), testval=conc_inits)
			BETA = pm.Normal('BETA', mu=0, tau=(1/1000000.0), shape=(S, num_cov), testval=beta_inits.T)
			#U = pm.Normal('U', mu=0, tau=(1/1.0), shape=(N, K), testval=np.random.randn(N, K))
			#U = pm.Exponential('U',lam=10.0, shape=(N, K), testval=np.abs(np.random.randn(N, K)))
			#V = pm.Normal('V', mu=0, tau=(1/100000.0), shape=(S, K), testval=np.random.randn(S, K))

			LAMBDA_U = pm.HalfCauchy('LAMBDA_U', beta=1, shape=(N, K), testval=np.ones((N,K)))
			TAU_U = pm.HalfCauchy('TAU_U', beta=1, testval=1.0)
			SIGMA_U = pm.Deterministic('SIGMA_U', TAU_U*TAU_U*LAMBDA_U*LAMBDA_U)
			U = pm.Normal('U', mu=0, sd=SIGMA_U, shape=(N, K), testval=np.random.randn(N, K))

			LAMBDA_V = pm.HalfCauchy('LAMBDA_V', beta=1, shape=(S, K), testval=np.ones((S,K)))
			TAU_V = pm.HalfCauchy('TAU_V', beta=1, testval=1.0)
			SIGMA_V = pm.Deterministic('SIGMA_V', TAU_V*TAU_V*LAMBDA_V*LAMBDA_V)
			V = pm.Normal('V', mu=0, sd=SIGMA_V, shape=(S, K), testval=np.random.randn(S, K))


			MU_A = pm.Normal("MU_A", mu=0., sd=100**2, shape=(1,S), testval=np.zeros((1,S)))
			SIGMA_A = pm.HalfCauchy("SIGMA_A", beta=5.0, shape=(1,S), testval=np.ones((1,S)))
			mu_a_mat = pm.math.dot(np.ones((I,1)), MU_A)
			sigma_a_mat = pm.math.dot(np.ones((I,1)), SIGMA_A)
			A = pm.Normal('A', mu=mu_a_mat, sigma=sigma_a_mat, shape=(I,S), testval=np.zeros((I, S)))

			p = pm.math.invlogit(pm.math.dot(X, BETA.T) + pm.math.dot(U,V.T) + A[Z,:])
			conc_mat = pm.math.dot(np.ones((N,1)), CONC)
			R = pm.BetaBinomial('like',alpha=(p*conc_mat)[~nans], beta=((1.0-p)*conc_mat)[~nans], n=n[~nans], observed=k[~nans])
			approx = pm.fit(method='advi', n=80000)
		pickle.dump(approx, open(self.output_root + '_model', 'wb'))
		approx = pickle.load( open(self.output_root + '_model', "rb" ) )
		means_dict = approx.bij.rmap(approx.params[0].eval())
		np.savetxt(self.output_root + '_temper_U.txt', (means_dict['U']), fmt="%s", delimiter='\t')
		np.savetxt(self.output_root + '_temper_V.txt', (means_dict['V'].T), fmt="%s", delimiter='\t')
		np.savetxt(self.output_root + '_temper_BETA.txt', (means_dict['BETA'].T), fmt="%s", delimiter='\t')


	def initialize_variables(self):
		self.N = self.allelic_counts.shape[0]
		self.T = self.allelic_counts.shape[1]
		self.num_cov = self.cov.shape[1]
		# Randomly initialize (only U matters)
		self.U = np.random.randn(self.N, self.K)
		self.V = np.random.randn(self.K, self.T)  # plus 1 for intercept
		self.conc = np.random.randn(self.T)
		self.C = np.random.randn(self.num_cov, self.T)

		'''
		for n in range(self.N):
			for t in range(self.T):
				if np.isnan(self.allelic_counts[n,t]) or np.isnan(self.total_counts[n,t]):
					self.allelic_counts[n,t] = 0
					self.total_counts[n,t] = 0
		'''



