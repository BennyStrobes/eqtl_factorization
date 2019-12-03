import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_vi_shared_effect
import eqtl_factorization_vi
import eqtl_factorization_vi_ard_loadings_only
import eqtl_factorization_vi_no_ard
import eqtl_factorization_vi_theta_fixed




# Load in sample overlap data
def load_in_sample_overlap_data(sample_overlap_file):
	f = open(sample_overlap_file)
	Z = []
	Z1 = []
	temp = []
	for line in f:
		line = line.rstrip()
		Z.append([int(line)])
		Z1.append(int(line))
		temp.append(int(line))
	f.close()
	num_individuals = max(temp) + 1
	return np.asarray(Z1), int(num_individuals)


def standardize_each_column_of_matrix(G):
	num_cols = G.shape[1]
	for col_num in range(num_cols):
		G[:,col_num] = (G[:,col_num] - np.mean(G[:,col_num]))/np.std(G[:,col_num])
	return G

def subset_matrices(Y,G,Z):
	initial_N = Y.shape[0]
	initial_T = Y.shape[1]

	new_sample_indices = np.arange(0,initial_N,2)
	new_test_indices = np.arange(0, initial_T, 3)

	new_Y = Y[new_sample_indices,:][:,new_test_indices]
	new_G = G[new_sample_indices,:][:,new_test_indices]
	new_Z = Z[new_sample_indices]

	return new_Y, new_G, new_Z

def train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, bernoulli_prob):
	############################
	# Load in data
	############################
	# Load in expression data (dimension: num_samplesXnum_tests)
	Y = np.transpose(np.loadtxt(expression_training_file, delimiter='\t'))
	# Load in genotype data (dimension: num_samplesXnum_tests)
	G = np.transpose(np.loadtxt(genotype_training_file, delimiter='\t'))
	# Load in sample overlap data
	Z,  num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	#Y,G,Z = subset_matrices(Y,G,Z)

	G = standardize_each_column_of_matrix(G)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]

	# RUN MODEL
	if model_name == 'vi_shared_effect':
		eqtl_vi = eqtl_factorization_vi_shared_effect.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1, beta=1, a=1, b=1, max_iter=600, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		np.savetxt(output_root + '_U.txt', eqtl_vi.U_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', eqtl_vi.U_mu*eqtl_vi.S_U, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', eqtl_vi.F_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F_S.txt', eqtl_vi.F_mu*eqtl_vi.S_F, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', eqtl_vi.V_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V_S.txt', eqtl_vi.V_mu*eqtl_vi.S_V, fmt="%s", delimiter='\t')	
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
	elif model_name == 'vi':
		eqtl_vi = eqtl_factorization_vi.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=1, max_iter=4000, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		np.savetxt(output_root + '_U.txt', eqtl_vi.U_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', eqtl_vi.U_mu*eqtl_vi.S_U, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', eqtl_vi.V_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V_S.txt', eqtl_vi.V_mu*eqtl_vi.S_V, fmt="%s", delimiter='\t')	
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
	elif model_name == 'vi_ard_loadings_only':
		eqtl_vi = eqtl_factorization_vi_ard_loadings_only.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=500, lambda_v=500,  max_iter=12000, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		np.savetxt(output_root + '_U.txt', eqtl_vi.U_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', eqtl_vi.U_mu*eqtl_vi.S_U, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', eqtl_vi.V_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
	elif model_name == 'vi_no_ard':
		eqtl_vi = eqtl_factorization_vi_no_ard.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, p_u=bernoulli_prob, lambda_u=lasso_param, lambda_v=1, lambda_f=.0000001, max_iter=600, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		np.savetxt(output_root + '_U.txt', eqtl_vi.U_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', eqtl_vi.U_mu*eqtl_vi.S_U, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', eqtl_vi.V_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
	elif model_name == 'vi_theta_fixed':
		eqtl_vi = eqtl_factorization_vi_theta_fixed.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, p_u=bernoulli_prob, lambda_v=1, lambda_f=.0000001, max_iter=600, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		np.savetxt(output_root + '_U.txt', eqtl_vi.U_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', eqtl_vi.U_mu*eqtl_vi.S_U, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', eqtl_vi.V_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')


#######################
# Command line args
#######################
sample_overlap_file = sys.argv[1]
expression_training_file = sys.argv[2]
genotype_training_file = sys.argv[3]
expression_testing_file = sys.argv[4]
genotype_testing_file = sys.argv[5]
num_latent_factors = int(sys.argv[6])
file_stem = sys.argv[7]
eqtl_results_dir = sys.argv[8]
seed = int(sys.argv[9])
model_name = sys.argv[10]
bernoulli_prob = float(sys.argv[11])

np.random.seed(seed)
# What to save output files to
output_root = eqtl_results_dir + file_stem 

#########################
# Train model
#########################
train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, bernoulli_prob)
