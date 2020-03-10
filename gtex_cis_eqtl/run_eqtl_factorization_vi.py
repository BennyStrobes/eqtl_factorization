import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_vi_shared_effect
import eqtl_factorization_vi_shared_effect_ard_only
import eqtl_factorization_vi_shared_effect_double_ard_only
import eqtl_factorization_vi
import eqtl_factorization_vi_shared_effect_factor_component_ard_only
import eqtl_factorization_vi_shared_effect_prior_on_loadings_only
import eqtl_factorization_vi_shared_effect_prior_on_loadings_only_special_init
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

def train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects):
	bernoulli_prob = .5
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
		#eqtl_vi = eqtl_factorization_vi_shared_effect.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=1, max_iter=600, delta_elbo_threshold=.01)
		#eqtl_vi.fit(G=G, Y=Y, z=Z)
		#pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		factor_ordering = np.where(factor_ve > .001)[0]
		U_temp = (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering]
		V_temp = (eqtl_vi.V_mu*eqtl_vi.S_V)[factor_ordering,:]

		np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering], fmt="%s", delimiter='\t')
	# RUN MODEL
	elif model_name == 'vi_shared_effect_factor_component_ard_only':
		eqtl_vi = eqtl_factorization_vi_shared_effect_factor_component_ard_only.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-6, beta=1e-6, a=1, b=10, max_iter=1000, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
	elif model_name == 'vi_shared_effect_prior_on_loadings_only':
		eqtl_vi = eqtl_factorization_vi_shared_effect_prior_on_loadings_only.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-6, beta=1e-6, a=1, b=1, gamma_v=1.0, max_iter=1000, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		#shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		#factor_ordering = np.where(factor_ve > .001)[0]
		#U_temp = (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering]
		#V_temp = (eqtl_vi.V_mu*eqtl_vi.S_V)[factor_ordering,:]
		#np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering], fmt="%s", delimiter='\t')
	elif model_name == 'vi_shared_effect_prior_on_loadings_only_special_init':
		eqtl_vi = eqtl_factorization_vi_shared_effect_prior_on_loadings_only_special_init.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=1, gamma_v=1.0, max_iter=1000, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		#shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		#factor_ordering = np.where(factor_ve > .001)[0]
		#U_temp = (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering]
		#V_temp = (eqtl_vi.V_mu*eqtl_vi.S_V)[factor_ordering,:]
		#np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering], fmt="%s", delimiter='\t')
	elif model_name == 'vi_shared_effect_ard_only':
		eqtl_vi = eqtl_factorization_vi_shared_effect_ard_only.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, max_iter=600, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
	elif model_name == 'vi_shared_effect_double_ard_only':
		#eqtl_vi = eqtl_factorization_vi_shared_effect_double_ard_only.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-10, beta=1e-10, max_iter=600, delta_elbo_threshold=.01)
		#eqtl_vi.fit(G=G, Y=Y, z=Z)
		#pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		###############
		## TEMP
		eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		factor_ordering = np.where(factor_ve > .001)[0]
		U_temp = (eqtl_vi.U_mu)[:,factor_ordering]
		V_temp = (eqtl_vi.V_mu)[factor_ordering,:]

		np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu)[:,factor_ordering], fmt="%s", delimiter='\t')

	elif model_name == 'vi':
		eqtl_vi = eqtl_factorization_vi.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=1, max_iter=1000, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		np.savetxt(output_root + '_U.txt', eqtl_vi.U_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', eqtl_vi.U_mu*eqtl_vi.S_U, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', eqtl_vi.V_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V_S.txt', eqtl_vi.V_mu*eqtl_vi.S_V, fmt="%s", delimiter='\t')	
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
	elif model_name == 'vi_ard_loadings_only':
		eqtl_vi = eqtl_factorization_vi_ard_loadings_only.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-5, beta=1e-5, a=1, b=1, lambda_v=1,  max_iter=12000, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		np.savetxt(output_root + '_U.txt', eqtl_vi.U_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', eqtl_vi.U_mu*eqtl_vi.S_U, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', eqtl_vi.V_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
	elif model_name == 'vi_no_ard':
		eqtl_vi = eqtl_factorization_vi_no_ard.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, p_u=bernoulli_prob, lambda_u=1, lambda_v=1, lambda_f=.0000001, max_iter=600, delta_elbo_threshold=.01)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		np.savetxt(output_root + '_U.txt', eqtl_vi.U_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', eqtl_vi.U_mu*eqtl_vi.S_U, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', eqtl_vi.V_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
	elif model_name == 'vi_no_ard_learn_bernoulli':
		#eqtl_vi = eqtl_factorization_vi_no_ard_learn_bernoulli.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=1, lambda_u=1, lambda_v=1, lambda_f=.0000001, max_iter=600, delta_elbo_threshold=.01)
		#eqtl_vi.fit(G=G, Y=Y, z=Z, random_effects=random_effects) 
		# Save the model
		#pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		###############
		## TEMP
		eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))

		# COMPUTE VE of each of the factors
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		# Order factors based on their theta_U
		thetas = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_a + eqtl_vi.theta_U_b)
		factor_ordering = np.argsort(-factor_ve)
		# Create list of factors that pass threshold
		num_factors_over_thresh = np.sum(factor_ve[factor_ordering] >= .001)
		filtered_factor_ordering = factor_ordering[:num_factors_over_thresh]


		# Save ELBO
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
		# Save learned model parameters without filtering factors
		np.savetxt(output_root + '_tau.txt', eqtl_vi.tau_alpha/eqtl_vi.tau_beta, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', eqtl_vi.F_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', eqtl_vi.intercept_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S_U.txt', (eqtl_vi.S_U)[:, factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu)[factor_ordering, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_theta_U.txt', (eqtl_vi.theta_U_a/(eqtl_vi.theta_U_a + eqtl_vi.theta_U_b))[factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_ve_shared_effect.txt', [shared_ve], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_ve_component_effects.txt', factor_ve[factor_ordering], fmt="%s", delimiter='\t')
		# Save learned model parameters after filtering factors
		np.savetxt(output_root + '_filtered_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,filtered_factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_filtered_S_U.txt', (eqtl_vi.S_U)[:, filtered_factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_filtered_V.txt', (eqtl_vi.V_mu)[filtered_factor_ordering, :], fmt="%s", delimiter='\t')

	elif model_name == 'vi_no_ard_learn_bernoulli_both_sides':
		#eqtl_vi = eqtl_factorization_vi_no_ard_learn_bernoulli_both_sides.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=1, lambda_u=1, lambda_v=1, lambda_f=.0000001, max_iter=600, delta_elbo_threshold=.01)
		#eqtl_vi.fit(G=G, Y=Y, z=Z)
		# Save the model
		#pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		###############
		## TEMP
		pdb.set_trace()
		eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()

		# Order factors
		factor_ordering = np.argsort(-factor_ve)
		# Create list of factors that pass threshold
		num_factors_over_thresh = np.sum(factor_ve[factor_ordering] >= .001)
		print(num_factors_over_thresh)
		filtered_factor_ordering = factor_ordering[:num_factors_over_thresh]

		# Save ELBO
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
		# Save learned model parameters without filtering factors
		np.savetxt(output_root + '_tau.txt', eqtl_vi.tau_alpha/eqtl_vi.tau_beta, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', eqtl_vi.F_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', eqtl_vi.intercept_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V_S.txt', (eqtl_vi.V_mu*eqtl_vi.S_V)[factor_ordering, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_theta_U.txt', (eqtl_vi.theta_U_a/(eqtl_vi.theta_U_a + eqtl_vi.theta_U_b))[factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_theta_V.txt', (eqtl_vi.theta_V_a/(eqtl_vi.theta_V_a + eqtl_vi.theta_V_b))[factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_ve_shared_effect.txt', [shared_ve], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_ve_component_effects.txt', factor_ve[factor_ordering], fmt="%s", delimiter='\t')
		# Save learned model parameters after filtering factors
		np.savetxt(output_root + '_filtered_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,filtered_factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_filtered_V_S.txt', (eqtl_vi.V_mu*eqtl_vi.S_V)[filtered_factor_ordering, :], fmt="%s", delimiter='\t')

	elif model_name == 'vi_no_ard_learn_bernoulli_all_sides':
		#eqtl_vi = eqtl_factorization_vi_no_ard_learn_bernoulli_all_sides.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-3, beta=1e-3, a=1, b=1, lambda_u=1, lambda_v=1, lambda_f=1, max_iter=600, delta_elbo_threshold=.01)
		#eqtl_vi.fit(G=G, Y=Y, z=Z)
		# Save the model
		#pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		###############
		## TEMP
		eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()

		# Order factors
		factor_ordering = np.argsort(-factor_ve)
		# Create list of factors that pass threshold
		num_factors_over_thresh = np.sum(factor_ve[factor_ordering] >= .001)
		filtered_factor_ordering = factor_ordering[:num_factors_over_thresh]

		# Save ELBO
		np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\t')
		# Save learned model parameters without filtering factors
		np.savetxt(output_root + '_tau.txt', eqtl_vi.tau_alpha/eqtl_vi.tau_beta, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F_S.txt', eqtl_vi.F_mu*eqtl_vi.S_F, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', eqtl_vi.intercept_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V_S.txt', (eqtl_vi.V_mu*eqtl_vi.S_V)[factor_ordering, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_theta_U.txt', (eqtl_vi.theta_U_a/(eqtl_vi.theta_U_a + eqtl_vi.theta_U_b))[factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_theta_V.txt', (eqtl_vi.theta_V_a/(eqtl_vi.theta_V_a + eqtl_vi.theta_V_b))[factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_ve_shared_effect.txt', [shared_ve], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_ve_component_effects.txt', factor_ve[factor_ordering], fmt="%s", delimiter='\t')
		# Save learned model parameters after filtering factors
		np.savetxt(output_root + '_filtered_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,filtered_factor_ordering], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_filtered_V_S.txt', (eqtl_vi.V_mu*eqtl_vi.S_V)[filtered_factor_ordering, :], fmt="%s", delimiter='\t')


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

np.random.seed(seed)
# What to save output files to
output_root = eqtl_results_dir + file_stem 

#########################
# Train model
#########################
train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects)
