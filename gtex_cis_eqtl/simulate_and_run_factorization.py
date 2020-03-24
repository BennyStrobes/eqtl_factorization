import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_vi_spike_and_slab
import eqtl_factorization_vi_spike_and_slab_no_gaussian_loading
import pickle

def string_to_boolean(stringer):
	if stringer == "True":
		return True
	elif stringer == "False":
		return False
	else:
		print("BOOLEAN NOT RECOGNIZED")
		return

def simulate_data_for_eqtl_factorization(I, Ni, T, K, alpha_0, beta_0, a_0, b_0):
	# Simulate gammas
	gamma_U = np.random.gamma(shape=alpha_0, scale=1.0/beta_0,size=K)
	gamma_V = np.random.gamma(shape=alpha_0, scale=1.0/beta_0,size=K)
	gamma_F = np.random.gamma(shape=alpha_0,scale=1.0/beta_0)
	# Simulate thetas
	theta_U = np.random.beta(a_0, b_0, size=K)
	theta_V = np.random.beta(a_0, b_0, size=K)
	theta_F = np.random.beta(a_0, b_0)
	# Simulate other variances
	psi = np.random.gamma(shape=alpha_0, scale=1.0/beta_0, size=T)
	tau = np.random.gamma(shape=alpha_0, scale=1.0/beta_0, size=T)
	# Simulate alphas (random effects)
	alphas = np.zeros((I, T))
	for test_index in range(T):
		for individual_index in range(I):
			alphas[individual_index, test_index] = np.random.normal(loc=0, scale=np.sqrt(1.0/psi[test_index]))
	# Simulate Us
	sample_num = 0
	U = np.zeros((I*Ni, K))
	S_U = np.zeros((I*Ni,K))
	Z = []
	for individual_index in range(I):
		for ni in range(Ni):
			Z.append(str(individual_index))
			for k in range(K):
				U[sample_num, k] = np.random.normal(loc=0, scale=np.sqrt(1.0/gamma_U[k]))
				S_U[sample_num, k] = np.random.binomial(n=1,p=theta_U[k])
			sample_num = sample_num + 1
	# Simulate Vs
	V = np.zeros((K, T))
	S_V = np.zeros((K,T))
	F = np.zeros(T)
	S_F = np.zeros(T)
	for test_index in range(T):
		F[test_index] = np.random.normal(loc=0, scale=np.sqrt(1.0/gamma_F))
		S_F[test_index] = np.random.binomial(n=1,p=theta_F)
		for k in range(K):
			V[k,test_index] = np.random.normal(loc=0, scale=np.sqrt(1.0/gamma_V[k]))
			S_V[k, test_index] = np.random.binomial(n=1,p=theta_V[k])
	# Simulate Genotypes
	G = np.zeros((I*Ni, T))
	for test_number in range(T):
		sample_num = 0
		genotypes = np.random.normal(0,1,size=I)
		for individual_index in range(I):
			for ni in range(Ni):
				G[sample_num, test_number] = genotypes[individual_index]
				sample_num = sample_num + 1
	# Simulate Expression
	pdb.set_trace()
	Y = np.zeros((I*Ni, T))
	sample_num = 0
	for individual_index in range(I):
		for ni in range(Ni):
			for test_number in range(T):
				predicted_mean = 0
				for k in range(K):
					predicted_mean = predicted_mean + U[sample_num,k]*S_U[sample_num,k]*V[k, test_number]*S_V[k,test_number]
				predicted_mean = predicted_mean + F[test_number]*S_F[test_number]
				predicted_mean = predicted_mean + alphas[individual_index, test_number]
				Y[sample_num, test_number] = np.random.normal(loc=predicted_mean, scale=np.sqrt(1.0/tau[test_number]))
			sample_num = sample_num + 1
	data = {}
	data['V'] = V
	data['S_V'] = S_V
	data['U'] = U
	data['S_U'] = S_U
	data['F'] = F
	data['S_F'] = S_F
	data['psi'] = psi
	data['tau'] = tau
	data['theta_U'] = theta_U
	data['theta_V'] = theta_V
	data['theta_F'] = theta_F
	data['gamma_V'] = gamma_V
	data['gamma_U'] = gamma_U
	data['gamma_F'] = gamma_F
	data['alpha'] = alphas
	return Y, G, Z, data


def simulate_data_for_eqtl_factorization2(I, Ni, T, K):
	N = I*Ni
	U_true = np.random.randn(N, K)
	V_true = np.random.randn(K, T)
	G = np.random.randn(N, T)
	F_true = np.random.randn(1,T)


	'''
	# Simulate thetas
	theta_U = np.random.beta(a_0, b_0, size=K)


	S_U = np.zeros((N,K))
	for n in range(N):
		for k in range(K):
			S_U[n,k] = np.random.binomial(n=1,p=theta_U[k])
	'''
	S_U = np.zeros((N,K))
	for n in range(N):
		number_of_active_components = np.random.binomial(K, p=1/K,size=1)[0]
		active_components = np.random.choice(range(K), size=number_of_active_components, replace=False)
		for component_num in active_components:
			S_U[n, component_num] = 1
	theta_U = np.sum(S_U,axis=0)/(S_U.shape[0])

	Y = G*np.dot(S_U*U_true, V_true) + G*np.dot(np.ones((N,1)),F_true)
	gene_residual_sdevs = np.sqrt(np.random.exponential(size=T))*10
	for m in range(T):
		Y[:,m] = Y[:, m] + np.random.normal(0, gene_residual_sdevs[m],size=N)
	sample_num = 0
	Z = []
	for individual_index in range(I):
		for ni in range(Ni):
			Z.append(str(individual_index))

	'''
	gene_random_effects_sdevs = np.sqrt(np.random.exponential(size=T))
	alphas = np.zeros((I, T))
	for t in range(T):
		sample_num = 0
		for individual_index in range(I):
			individual_intercept = np.random.normal(0, gene_random_effects_sdevs[t])
			alphas[individual_index, t] = individual_intercept
			for ni in range(Ni):
				Y[sample_num, t] = Y[sample_num, t] + individual_intercept
				sample_num = sample_num + 1
	'''
	data = {}
	data['U'] = U_true
	data['V'] = V_true
	data['F'] = F_true
	data['resid'] = gene_residual_sdevs
	data['theta_U'] = theta_U
	data['S_U'] = S_U
	return Y, G, Z, data

#######################
# Command line args
#######################
output_root = sys.argv[1]  # Output root (where to save things)
svi_boolean = string_to_boolean(sys.argv[2]) # Whether to use SVI or not
I = int(sys.argv[3])  # Number of individuals
Ni = int(sys.argv[4])  # Number of samples per individual
T = int(sys.argv[5])  # number of tests
K = int(sys.argv[6])  # Number of latent factors
seed = int(sys.argv[7])  # random seed
parrallel_boolean = string_to_boolean(sys.argv[8])  # Whether to parrallelize loops

#########################
# Set Seed
#########################
np.random.seed(seed)

##########################
# Parameters
##########################
alpha_0 = 1e-3
beta_0 = 1e-3
a_0 = 1
b_0 = 1
max_iter=1000
gamma_v=1.0
# Only applicable to if SVI_boolean==True
sample_batch_fraction=0.3
learning_rate=0.4
forgetting_rate=0.01

##########################
# Simulate data
##########################
Y, G, z, data = simulate_data_for_eqtl_factorization2(I, Ni, T, K)

##################
# Fit eqtl factorization using home-built variational inference
##################
eqtl_vi = eqtl_factorization_vi_spike_and_slab.EQTL_FACTORIZATION_VI(K=20, alpha=alpha_0, beta=beta_0, a=a_0, b=b_0, max_iter=max_iter, gamma_v=gamma_v, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=sample_batch_fraction, learning_rate=learning_rate, forgetting_rate=forgetting_rate)
eqtl_vi.fit(G=G, Y=Y, z=z)
pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))

if eqtl_vi.SVI == True:
	np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu_full*eqtl_vi.S_U_full), fmt="%s", delimiter='\t')
else:
	np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U), fmt="%s", delimiter='\t')
np.savetxt(output_root + 'actual_U_S.txt', (data['U']*data['S_U']), fmt="%s", delimiter='\t')

