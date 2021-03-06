import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_vi_spike_and_slab
import pickle




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

def string_to_boolean(stringer):
	if stringer == "True":
		return True
	elif stringer == "False":
		return False
	else:
		print("BOOLEAN NOT RECOGNIZED")
		return

def train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean):
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

	#G = standardize_each_column_of_matrix(G)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]

	if svi_boolean == True:
		max_it = 1000
	elif svi_boolean == False:
		max_it = 1000

	# RUN MODEL
	if model_name == 'eqtl_factorization_vi_spike_and_slab':
		eqtl_vi = eqtl_factorization_vi_spike_and_slab.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=1, gamma_v=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		ordered_indices = np.argsort(-eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a))
		num_indices = sum(eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a) > .05)
		ordered_filtered_indices = ordered_indices[:(num_indices)]
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full*eqtl_vi.S_U_full)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		#np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\n')


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
random_effects = string_to_boolean(sys.argv[11])
svi_boolean = string_to_boolean(sys.argv[12])
parrallel_boolean = string_to_boolean(sys.argv[13])

np.random.seed(seed)
# What to save output files to
output_root = eqtl_results_dir + file_stem 


#########################
# Train model
#########################
train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean)
