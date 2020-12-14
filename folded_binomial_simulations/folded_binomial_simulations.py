import numpy as np 
import os
import sys
import pdb
import pystan
import pickle
import scipy.special


FOLDED_BINOMIAL_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/folded_binomial_glm.pkl', 'rb'))
FOLDED_BINOMIAL_INTERCEPT_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/folded_binomial_with_intercept_glm.pkl', 'rb'))
BINOMIAL_GLM = pickle.load(open('/work-zfs/abattle4/bstrober/temp/binomial_glm.pkl', 'rb'))

def simulation1():
	N = 10000  # Number of samples
	num_cov = 1  # Number of covariates
	nn = np.random.negative_binomial(1, 0.25, N) + 2  # draw total counts from negative binomial
	p = .3
	k = []
	for index in range(N):
		counts = np.random.binomial(nn[index], p)
		if np.random.random() < .5:
			counts = np.random.binomial(nn[index], 1-p)
		counts = np.min((counts, nn[index] - counts))
		k.append(counts)
	k = np.asarray(k)

	# Fit models
	data = dict(N=N, P=1, x=np.ones((N,1)), ys=k.astype(int), ns=nn.astype(int))
	op = FOLDED_BINOMIAL_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1)
	learned_p = scipy.special.expit(op['beta'])
	op2 = BINOMIAL_GLM.optimizing(data = data, verbose=False,iter=5000,seed=1)
	learned_p2 = scipy.special.expit(op2['beta'])
	print(learned_p)
	print(learned_p2)
	pdb.set_trace()









output_root = sys.argv[1]


simulation1()


