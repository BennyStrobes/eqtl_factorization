import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_vi_spike_and_slab
import eqtl_factorization_vi_spike_and_slab_loadings_ard_factors
import pickle
import h5py
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt



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
	Z,  num_individuals = load_in_sample_overlap_data(sample_overlap_file)

	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]

	max_it = 1000
	print('datA loaded')


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
		num_indices = sum(eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a) > .05)
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


def debug_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean):
	print('start')
	eqtl_vi = pickle.load(open(output_root + '_model', 'rb'))
	# Save factor PVE
	# Plot distribution of factor weights (are some stronger than others?)
	# raw_expression = np.transpose(np.asarray(h5py.File('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_input/sc_raw_expression_training_data_uncorrected_10000_bp_0.5_r_squared_pruned.h5','r')['data']))
	# raw_genotype = np.transpose(np.asarray(h5py.File('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_input/sc_genotype_training_data_uncorrected_10000_bp_0.5_r_squared_pruned.h5','r')['data']))
	test_names = np.loadtxt('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_input/pseudobulk_sig_tests_50_pc_variant_gene_pairs_in_known_cell_types.txt',dtype=str,delimiter='\t')
	gene_names = test_names[1:,0]
	variant_names = test_names[1:,1]
	factor_num = 7
	factor_weights = np.abs(eqtl_vi.V_mu[factor_num,:])
	loading_weights = (eqtl_vi.U_mu*eqtl_vi.S_U)[:, factor_num]
	ordered_tests = np.argsort(-factor_weights)
	geno1 = eqtl_vi.Y_full[:, ordered_tests[0]]
	counter = 0
	for test_num in ordered_tests:
		print('###########')
		print(counter)
		print(gene_names[test_num])
		print(eqtl_vi.V_mu[factor_num,test_num])
		print(np.corrcoef(geno1,  eqtl_vi.Y_full[:, test_num])[0,1])
		pdb.set_trace()
		if counter < 5:

			fig = plt.figure()
			plt.scatter(eqtl_vi.G_full[:, test_num], loading_weights,c=eqtl_vi.Y_full[:, test_num], s=.1)
			plt.xlabel('Normalized Genotype (' + variant_names[test_num] + ')')
			plt.ylabel('Loading ' + str(factor_num))
			plt.title('Factor ' + str(factor_num) + ' / Test weight: ' + str(np.round(eqtl_vi.V_mu[factor_num,test_num], decimals=2)))
			cbar = plt.colorbar()
			cbar.set_label(gene_names[test_num] + ' Expression')
			fig.savefig('factor_' + str(factor_num) + '_test_' + str(counter) + '_genotype_vs_loading.png')

			fig = plt.figure()
			plt.scatter(eqtl_vi.G_full[:, test_num]*loading_weights, eqtl_vi.Y_full[:, test_num], s=.1)
			plt.xlabel('Loading' + str(factor_num) + '*Normalized Genotype (' + variant_names[test_num] + ')')
			plt.ylabel(gene_names[test_num] + ' Expression')
			plt.title('Factor ' + str(factor_num) + ' / Test weight: ' + str(np.round(eqtl_vi.V_mu[factor_num,test_num], decimals=2)))
			fig.savefig('factor_' + str(factor_num) + '_test_' + str(counter) + '_loading_times_genotype_vs_expression.png')
			
			fig = plt.figure()
			plt.scatter(eqtl_vi.G_full[:, test_num], eqtl_vi.Y_full[:, test_num], c=loading_weights, s=.1)
			plt.xlabel('Normalized Genotype (' + variant_names[test_num] + ')')
			plt.ylabel(gene_names[test_num] + ' Expression')
			cbar = plt.colorbar()
			cbar.set_label('Loading ' + str(factor_num))
			plt.title('Factor ' + str(factor_num) + ' / Test weight: ' + str(np.round(eqtl_vi.V_mu[factor_num,test_num], decimals=2)))
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

np.random.seed(seed)
# What to save output files to
output_root = eqtl_results_dir + file_stem 


#########################
# Train model
#########################
#train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean)


debug_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, output_root, model_name, random_effects, svi_boolean, parrallel_boolean)


