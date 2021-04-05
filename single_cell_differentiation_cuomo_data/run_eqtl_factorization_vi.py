import numpy as np 
import os
import sys
import pdb
#import eqtl_factorization_vi_spike_and_slab
#import eqtl_factorization_vi_dirichlet_simplex
#import eqtl_factorization_vi_spike_and_slab_tied_residuals
#import eqtl_factorization_vi_spike_and_slab_loadings_ard_factors
#import eqtl_factorization_vi_spike_and_slab_loadings_ard_loadings
#import eqtl_factorization_vi_zero_inflated
import eqtl_factorization_vi_zero_inflated2
import eqtl_factorization_vi_with_re
import eqtl_factorization_vi_with_re_tied_variance
import eqtl_factorization_vi_loading_spike_and_slab_with_re
import eqtl_factorization_vi_factor_loading_spike_and_slab_with_re
import eqtl_factorization_vi_factor_loading_spike_and_slab_with_multiple_re
import eqtl_factorization_vi_factor_loading_spike_and_slab_with_re_and_fixed_environmental_effect
import eqtl_factorization_vi_factor_loading_spike_and_slab_with_re_and_fixed_environmental_effect_learn_cov

#import eqtl_factorization_vi_zero_inflated3
#import eqtl_factorization_als
#import eqtl_factorization_als_unconstrained
#import eqtl_factorization_als_positive_loadings
#import eqtl_factorization_als_constrained
#import eqtl_factorization_als_constrained_weighted_regularization
#import eqtl_factorization_als_gumbel_softmax_constrained
import pickle
import h5py
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt



# Load in sample overlap data
def load_in_sample_overlap_data(sample_overlap_file):
	f = open(sample_overlap_file)
	Z1 = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		#Z1.append(np.asarray(data).astype(int))
		Z1.append(int(line))
	f.close()
	return np.asarray(Z1)


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

def train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean, lasso_param, covariate_file):
	############################
	# Load in data
	############################
	print('start loading')
	# Load in expression data (dimension: num_samplesXnum_tests)
	if expression_training_file.endswith('.txt'):
		Y = np.transpose(np.loadtxt(expression_training_file, delimiter='\t'))
	elif expression_training_file.endswith('.h5'):
		Y = np.transpose(np.asarray(h5py.File(expression_training_file,'r')['data']))
	# Load in genotype data (dimension: num_samplesXnum_tests)
	if genotype_training_file.endswith('.txt'):
		G = np.transpose(np.loadtxt(genotype_training_file, delimiter='\t'))
	elif genotype_training_file.endswith('.h5'):
		G = np.transpose(np.asarray(h5py.File(genotype_training_file,'r')['data']))
	# Load in sample overlap data
	Z = load_in_sample_overlap_data(sample_overlap_file)

	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]

	max_it = 1000
	print('datA loaded')
	print(G.shape)


	############################
	# RUN MODEL
	#############################
	if model_name == 'eqtl_factorization_vi_spike_and_slab':
		eqtl_vi = eqtl_factorization_vi_spike_and_slab.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.1)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		ordered_indices = np.argsort(-eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a))
		num_indices = sum(eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a) > .025)
		ordered_filtered_indices = ordered_indices[:(num_indices)]
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full*eqtl_vi.S_U_full)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		#np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\n')
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
	if model_name == 'eqtl_factorization_vi_zero_inflated':
		eqtl_vi = eqtl_factorization_vi_zero_inflated.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=1.0, gamma_u=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_zero_inflated2':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_zero_inflated2.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=float(lasso_param), gamma_u=float(lasso_param), max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_with_re':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_with_re.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=float(lasso_param), gamma_u=float(lasso_param), max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_loading_spike_and_slab_with_re_and_fixed_environmental_effect':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_loading_spike_and_slab_with_re_and_fixed_environmental_effect.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=float(lasso_param), gamma_u=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')		
	if model_name == 'eqtl_factorization_vi_loading_spike_and_slab_with_re':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_loading_spike_and_slab_with_re.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=float(lasso_param), gamma_u=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_factor_loading_spike_and_slab_with_re':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_factor_loading_spike_and_slab_with_re.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=1.0, gamma_u=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_factor_loading_spike_and_slab_with_multiple_re':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_factor_loading_spike_and_slab_with_multiple_re.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=1.0, gamma_u=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_fixed_environmental_effect':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_factor_loading_spike_and_slab_with_re_and_fixed_environmental_effect.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=float(lasso_param), gamma_u=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_fixed_environmental_effect_learn_cov':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_factor_loading_spike_and_slab_with_re_and_fixed_environmental_effect_learn_cov.EQTL_FACTORIZATION_VI(K=num_latent_factors, J=20, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=float(lasso_param), gamma_u=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_with_re_tied_variance':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_with_re_tied_variance.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=float(lasso_param), gamma_u=float(lasso_param), max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', np.asmatrix(eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_zero_inflated3':
		cov = np.loadtxt(covariate_file, delimiter='\t')
		eqtl_vi = eqtl_factorization_vi_zero_inflated3.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=1.0, gamma_u=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2)
		eqtl_vi.fit(G=G, Y=Y, cov=cov, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_intercept.txt', (eqtl_vi.intercept_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_S.txt', (eqtl_vi.S), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_cov.txt', (eqtl_vi.cov_mu), fmt="%s", delimiter='\t')

	if model_name == 'eqtl_factorization_vi_dirichlet_simplex':
		eqtl_vi = eqtl_factorization_vi_dirichlet_simplex.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.2)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_als':
		eqtl_vi = eqtl_factorization_als.EQTL_FACTORIZATION_ALS(K=num_latent_factors, max_iter=max_it, parrallel_boolean=parrallel_boolean)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		np.savetxt(output_root + '_U.txt', (eqtl_vi.U), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V), fmt="%s", delimiter='\t')
		# shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		#if svi_boolean == False:
		#	np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		#elif svi_boolean == True:
		#	np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full), fmt="%s", delimiter='\t')
		#p.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_als_positive_loadings':
		genotype_intercept=True

		output_root = output_root + '_geno_int_' + str(genotype_intercept) + '_lasso_param_' + str(lasso_param)
		eqtl_vi = eqtl_factorization_als_positive_loadings.EQTL_FACTORIZATION_ALS_CONSTRAINED(K=num_latent_factors, genotype_intercept=genotype_intercept, lasso_param_v=lasso_param, max_iter=max_it, parrallel_boolean=parrallel_boolean)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		np.savetxt(output_root + '_U.txt', (eqtl_vi.U), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_als_unconstrained':
		genotype_intercept=True

		output_root = output_root + '_geno_int_' + str(genotype_intercept) + '_lasso_param_' + str(lasso_param)
		eqtl_vi = eqtl_factorization_als_unconstrained.EQTL_FACTORIZATION_ALS_CONSTRAINED(K=num_latent_factors, genotype_intercept=genotype_intercept, lasso_param_v=lasso_param, max_iter=max_it, parrallel_boolean=parrallel_boolean)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		np.savetxt(output_root + '_U.txt', (eqtl_vi.U), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_als_constrained':
		genotype_intercept=True

		output_root = output_root + '_geno_int_' + str(genotype_intercept) + '_lasso_param_' + str(lasso_param)
		eqtl_vi = eqtl_factorization_als_constrained.EQTL_FACTORIZATION_ALS_CONSTRAINED(K=num_latent_factors, genotype_intercept=genotype_intercept, lasso_param_v=lasso_param, max_iter=max_it, parrallel_boolean=parrallel_boolean)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		np.savetxt(output_root + '_U.txt', (eqtl_vi.U), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_als_gumbel_softmax_constrained':
		genotype_intercept=True
		temperature = 100.0
		output_root = output_root + '_geno_int_' + str(genotype_intercept) + '_lasso_param_' + str(lasso_param) + '_temp_' + str(temperature)
		eqtl_vi = eqtl_factorization_als_gumbel_softmax_constrained.EQTL_FACTORIZATION_ALS_CONSTRAINED(K=num_latent_factors, genotype_intercept=genotype_intercept, lasso_param_v=lasso_param, max_iter=max_it, parrallel_boolean=parrallel_boolean, temperature=temperature)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		np.savetxt(output_root + '_U.txt', (eqtl_vi.U), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V), fmt="%s", delimiter='\t')

	if model_name == 'eqtl_factorization_als_constrained_weighted_regularization':
		genotype_intercept=True
		output_root = output_root + '_geno_int_' + str(genotype_intercept) + '_lasso_param_' + str(lasso_param)
		eqtl_vi = eqtl_factorization_als_constrained_weighted_regularization.EQTL_FACTORIZATION_ALS_CONSTRAINED(K=num_latent_factors, genotype_intercept=genotype_intercept, lasso_param_v=lasso_param, max_iter=max_it, parrallel_boolean=parrallel_boolean)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		np.savetxt(output_root + '_U.txt', (eqtl_vi.U), fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V), fmt="%s", delimiter='\t')

	if model_name == 'eqtl_factorization_vi_spike_and_slab_tied_residuals':
		print('tied resid')
		eqtl_vi = eqtl_factorization_vi_spike_and_slab_tied_residuals.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=1, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.3)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		ordered_indices = np.argsort(-eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a))
		num_indices = sum(eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a) > .025)
		ordered_filtered_indices = ordered_indices[:(num_indices)]
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full*eqtl_vi.S_U_full)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		#np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\n')
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
	elif model_name == 'eqtl_factorization_vi_spike_and_slab_loadings_ard_factors':
		print('ard')
		eqtl_vi = eqtl_factorization_vi_spike_and_slab_loadings_ard_factors.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=1.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.1)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		ordered_indices = np.argsort(-eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a))
		num_indices = sum(eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a) > .01)
		ordered_filtered_indices = ordered_indices[:(num_indices)]
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full*eqtl_vi.S_U_full)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		#np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\n')
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
	elif model_name == 'eqtl_factorization_vi_spike_and_slab_loadings_ard_loadings':
		print('ard')
		eqtl_vi = eqtl_factorization_vi_spike_and_slab_loadings_ard_loadings.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=1e-16, beta=1e-16, a=1, b=1, gamma_v=100.0, max_iter=max_it, delta_elbo_threshold=.01, SVI=svi_boolean, parrallel_boolean=parrallel_boolean, sample_batch_fraction=.1)
		eqtl_vi.fit(G=G, Y=Y, z=Z)
		# pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
		shared_ve, factor_ve = eqtl_vi.compute_variance_explained_of_factors()
		ordered_indices = np.argsort(-eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a))
		num_indices = sum(eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a) > .01)
		ordered_filtered_indices = ordered_indices[:(num_indices)]
		if svi_boolean == False:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		elif svi_boolean == True:
			np.savetxt(output_root + '_U_S.txt', (eqtl_vi.U_mu_full*eqtl_vi.S_U_full)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + '_F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		#np.savetxt(output_root + '_elbo.txt', eqtl_vi.elbo, fmt="%s", delimiter='\n')
		pickle.dump(eqtl_vi, open(output_root + '_model', 'wb'))
		#eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))


def comparison_across_runs(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean):
	#U_full = np.loadtxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_U_S_v6.txt')
	file_seed_1 = '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_factorization_results/eqtl_factorization_single_cell_sig_tests_50_pc_data_8_factors_eqtl_factorization_als_model_False_re_False_svi_1_seed_model'
	file_seed_2 = '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_factorization_results/eqtl_factorization_single_cell_sig_tests_50_pc_data_8_factors_eqtl_factorization_als_model_False_re_False_svi_2_seed_model'
	eqtl_vi_1 = pickle.load(open(file_seed_1, 'rb'))
	eqtl_vi_2 = pickle.load(open(file_seed_2, 'rb'))
	# Are different seeds similar
	# Variance explained compared to if we fix U by cell type
	# Variance explained compared to if we fix V by known pseudobulk effect sizes
	pdb.set_trace()


def debug_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean, lasso_param):
	print('start')
	output_root = output_root + '_geno_int_' + str(True) + '_lasso_param_' + str(lasso_param)

	eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))

	pdb.set_trace()
	# Save factor PVE
	# Plot distribution of factor weights (are some stronger than others?)
	# raw_expression = np.transpose(np.asarray(h5py.File('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_input/sc_raw_expression_training_data_uncorrected_10000_bp_0.5_r_squared_pruned.h5','r')['data']))
	# raw_genotype = np.transpose(np.asarray(h5py.File('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_input/sc_genotype_training_data_uncorrected_10000_bp_0.5_r_squared_pruned.h5','r')['data']))
	#V_full = np.loadtxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_V_v6.txt')
	#U_full = np.loadtxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_U_S_v6.txt')
	#tau = np.loadtxt('/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_tau_v6.txt')
	test_names = np.loadtxt('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_input/single_cell_random_subset_sig_tests_50_pc_variant_gene_pairs_in_known_cell_types.txt',dtype=str,delimiter='\t')
	Y_full = np.transpose(np.asarray(h5py.File(expression_training_file,'r')['data']))
	G_full = np.transpose(np.asarray(h5py.File(genotype_training_file,'r')['data']))

	V_full = eqtl_vi.V[1:,:]
	U_full = eqtl_vi.U

	gene_names = test_names[1:,0]
	variant_names = test_names[1:,1]
	#ordered_indices = np.argsort(-theta)
	V = V_full
	S_U = U_full
	factor_num = 2

	factor_weights = np.abs(V[factor_num,:])
	loading_weights = (S_U)[:, factor_num]
	ordered_tests = np.argsort(-factor_weights)
	geno1 = G_full[:, ordered_tests[0]]
	expr1 = Y_full[:, ordered_tests[0]]
	counter = 1
	for test_num in ordered_tests:
		print('###########')
		print(counter)
		print(gene_names[test_num])
		print(V[factor_num,test_num])
		print(np.corrcoef(geno1,  G_full[:, test_num])[0,1])
		print(np.corrcoef(expr1,  Y_full[:, test_num])[0,1])
		pdb.set_trace()
		if counter < 5:

			fig = plt.figure()
			plt.scatter(G_full[:, test_num], loading_weights,c=Y_full[:, test_num], s=.1)
			plt.xlabel('Normalized Genotype (' + variant_names[test_num] + ')')
			plt.ylabel('Loading ' + str(factor_num))
			plt.title('Factor ' + str(factor_num) + ' / Test weight: ' + str(np.round(V[factor_num,test_num], decimals=2)))
			cbar = plt.colorbar()
			cbar.set_label(gene_names[test_num] + ' Expression')
			fig.savefig('factor_' + str(factor_num) + '_test_' + str(counter) + '_genotype_vs_loading.png')

			fig = plt.figure()
			plt.scatter(G_full[:, test_num]*loading_weights, Y_full[:, test_num], s=.1)
			plt.xlabel('Loading' + str(factor_num) + '*Normalized Genotype (' + variant_names[test_num] + ')')
			plt.ylabel(gene_names[test_num] + ' Expression')
			plt.title('Factor ' + str(factor_num) + ' / Test weight: ' + str(np.round(V[factor_num,test_num], decimals=2)))
			fig.savefig('factor_' + str(factor_num) + '_test_' + str(counter) + '_loading_times_genotype_vs_expression.png')
			
			fig = plt.figure()
	
			plt.scatter(G_full[:, test_num], Y_full[:, test_num], c=loading_weights, s=.5)
			plt.xlabel('Normalized Genotype (' + variant_names[test_num] + ')')
			plt.ylabel(gene_names[test_num] + ' Expression')
			cbar = plt.colorbar()
			cbar.set_label('Loading ' + str(factor_num))
			plt.title('Factor ' + str(factor_num) + ' / Test weight: ' + str(np.round(V[factor_num,test_num], decimals=2)))
			fig.savefig('factor_' + str(factor_num) + '_test_' + str(counter) + '_genotype_vs_expression.png')
		counter = counter + 1

	'''
	print('start3')
	nn = 100
	num_factors = eqtl_vi.U_mu_full.shape[1]
	for factor_num in range(num_factors):
		pred_0 = np.abs(np.dot(np.asmatrix((eqtl_vi.U_mu_full*eqtl_vi.S_U_full)[:,factor_num]).T, np.asmatrix(eqtl_vi.V_mu[factor_num,:])))
		flat_pred_0 = pred_0.flatten().A1
		filtered = flat_pred_0[flat_pred_0 > .5]
		sort_filt = sorted(-filtered)
		if len(sort_filt) < nn:
			print('assumption error')
			pdb.set_trace()
		thresh = -sort_filt[100]
		x_indices = np.where(pred_0 > thresh)[0]
		y_indices = np.where(pred_0 > thresh)[1]
		num_indices = len(x_indices)
		valz = []
		for index in range(num_indices):
			x_index = x_indices[index]
			y_index = y_indices[index]
			valz.append(raw_expression[x_index, y_index])
		pdb.set_trace()
		print(factor_num)
		print(valz)
	'''

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
lasso_param_v = float(sys.argv[14])
covariate_file = sys.argv[15]

np.random.seed(seed)
# What to save output files to
output_root = eqtl_results_dir + file_stem 


#########################
# Train model
#########################
train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean, lasso_param_v, covariate_file)


# debug_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean, lasso_param_v)

#comparison_across_runs(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean)
