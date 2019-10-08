#from __future__ import absolute_import
#rom __future__ import division
#from __future__ import print_function
#from edward.models import Gamma, Poisson, Normal, PointMass, \
#    TransformedDistribution
#from edward.models import PointMass, Empirical, HalfNormal
#import edward as ed
#import tensorflow as tf
import numpy as np 
import os
import sys
import pdb
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import statsmodels.api as sm
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from joblib import Parallel, delayed
import multiprocessing


# Load in sample overlap data
def load_in_sample_overlap_data(sample_overlap_file):
	f = open(sample_overlap_file)
	Z = []
	temp = []
	for line in f:
		line = line.rstrip()
		Z.append(line)
		temp.append(int(line))
	f.close()
	num_individuals = max(temp) + 1
	return Z, int(num_individuals)

def initialize_sample_loading_matrix_with_kmeans(num_samples, num_latent_factor, Y):
	# Sample loading matrix
	U = np.zeros((num_samples, num_latent_factor))
	# Fit kmeans model
	kmeans = KMeans(n_clusters=num_latent_factor, random_state=0).fit(Y)
	# Fill in sample loading matrix with kmeans assignments
	for sample_index, kmeans_assignment in enumerate(kmeans.labels_):
		U[sample_index, kmeans_assignment] = 1.0
	return U

def update_factor_matrix_one_test(test_number, Y, G, U, Z, lasso_param):
	# Get slice of data corresponding to this test
	y_test = Y[:, test_number]
	g_test = G[:, test_number]
	# Get U scaled by genotype for this test
	U_scaled = U*g_test[:,None]

	X = np.hstack((U, U_scaled))
	# Get U for intercept terms
	'''
	# Fit linear regression model
	if lasso_param == 0:
		reg = LinearRegression().fit(U_scaled,y_test)
		params = np.hstack((reg.intercept_, reg.coef_))
	# Fit lasso model
	else:
		clf = linear_model.Lasso(alpha=lasso_param, fit_intercept=True)
		clf.fit(U_scaled, y_test)
		# Get params of fitted model
		params = np.hstack((clf.intercept_, clf.coef_))
	'''
	# Fit lasso model
	num_features = U.shape[1]
	regularization_vector = np.zeros(2*num_features)
	regularization_vector[num_features:] = regularization_vector[num_features:] + lasso_param
	model_fit = sm.OLS(y_test, X).fit_regularized(alpha = regularization_vector, L1_wt=1.0)

	#clf = linear_model.Lasso(alpha=lasso_param, fit_intercept=False)
	#clf.fit(X, y_test)
	#params = np.hstack((clf.intercept_, clf.coef_))
	#pdb.set_trace()

	return model_fit.params

# Update factor matrix (V) with linear mixed model
def update_factor_matrix(Y, G, U, Z, num_samples, num_tests, num_latent_factor, lasso_param):
	num_cores = multiprocessing.cpu_count()
	#print(num_cores)
	#V_arr = Parallel(n_jobs=num_cores)(delayed(update_factor_matrix_one_test)(test_number, Y, G, U, Z) for test_number in range(num_tests))
	genotype_intercept = 'False'
	if genotype_intercept == 'False':
		V_arr = []
		for test_number in range(num_tests):
			print(test_number)
			V_arr.append(update_factor_matrix_one_test(test_number, Y, G, U, Z, lasso_param))
		#V_arr = Parallel(n_jobs=num_cores)(delayed(update_factor_matrix_one_test)(test_number, Y, G, U, Z, penalty_vector) for test_number in range(num_tests))
		# Convert back to matrix
		V = np.transpose(np.asmatrix(V_arr))
		#intercept = np.squeeze(np.asarray(V_full[0,:]))
		#intercept_mat = (np.ones((num_samples,1))*np.asmatrix(intercept))
		#V = np.asarray(V_full[1:,:])
	return V

def update_loading_matrix_one_sample(sample_num, Y, G, V, Z, intercept, lasso_param):
	# Get slice of data corresponding to this sample
	y_test = Y[sample_num, :] - intercept
	g_test = G[sample_num, :]

	# Get V scaled by genotype for this sample
	V_scaled = np.transpose(V)*g_test[:,None]

	if lasso_param == 0:
		# Mimick being close to zero
		clf = linear_model.Lasso(alpha=0.0000000000001,positive=True, fit_intercept=False)
		clf.fit(V_scaled, y_test)
	# Fit Lasso model
	else:
		clf = linear_model.Lasso(alpha=lasso_param,positive=True, fit_intercept=False)
		#clf = linear_model.Lasso(alpha=lasso_param, fit_intercept=False)
		clf.fit(V_scaled, y_test)
	return np.asarray(clf.coef_)

# Update loading matrix (U) with l1 penalized linear model
def update_loading_matrix(Y, G, V, Z, intercept, num_samples, num_tests, num_latent_factor, lasso_param):
	# Serial version
	U = np.zeros((num_samples, num_latent_factor))
	for sample_num in range(num_samples):
		U[sample_num,:] = update_loading_matrix_one_sample(sample_num, Y, G, V, Z, intercept, lasso_param)
	return U

def initialize_sample_loading_with_residual_clustering(Y, G, Z, num_latent_factor, output_root):
	# Get dimensions of the matrix
	num_tests = Y.shape[1]
	num_samples = Y.shape[0]
	# Initialize residual matrix
	resid_matrix = np.zeros((Y.shape))

	# Loop through tests and calculate residual for each test
	for test_number in range(num_tests):
		# Extract data relevent to tests
		y_test = Y[:, test_number]
		g_test = G[:, test_number]
		# Fit linear model (eQTL) across samples
		reg = LinearRegression().fit(g_test.reshape(-1,1), y_test.reshape(-1,1))

		# Calculate residuals according to predicted model
		pred_test = reg.predict(g_test.reshape(-1,1))
		resid_test = np.squeeze(np.asarray(pred_test)) - y_test

		# Fit linear mixed model (eQTL) across samples
		#cov = np.hstack((np.ones((len(g_test), 1)), np.transpose(np.asmatrix(g_test))))
		#model = sm.MixedLM(y_test, cov, Z)
		#lmm_reg = model.fit()
		#lmm_pred_test = lmm_reg.predict(cov)

		#corry = np.corrcoef(np.squeeze(np.asarray(pred_test)), lmm_pred_test)
		#print(corry)
		resid_matrix[:, test_number] = resid_test
	# Fit kmeans model
	kmeans = KMeans(n_clusters=num_latent_factor, random_state=0).fit(np.abs(resid_matrix))
	U_init = np.zeros((num_samples, num_latent_factor))
	for sample_num, label in enumerate(kmeans.labels_):
		U_init[sample_num, label] = 1.0
	# Add small constant to each entry
	U_init = U_init + .1

	np.savetxt(output_root + '_U_init.txt', U_init, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_lm_residuals.txt', resid_matrix, fmt="%s", delimiter='\t')
	# Old code for GMM instead of KMEANS (Used kmeans because of computational efficiency)
	#GMM = GaussianMixture(n_components=num_latent_factor).fit(np.abs(resid_matrix))
	#U_init2 = GMM.predict_proba(resid_matrix) + .1
	return U_init

def standardize_each_column_of_matrix(G):
	num_cols = G.shape[1]
	for col_num in range(num_cols):
		G[:,col_num] = (G[:,col_num] - np.mean(G[:,col_num]))/np.std(G[:,col_num])
	return G

#######################
# Part 1: Train eqtl factorization model
########################
def train_eqtl_factorization_model_em_version(sample_overlap_file, expression_file, genotype_file, num_latent_factor, lasso_param_u, lasso_param_v, initialization, output_root):
	# Load in expression data (dimension: num_samplesXnum_tests)
	Y = np.transpose(np.loadtxt(expression_file, delimiter='\t'))
	# Load in genotype data (dimension: num_samplesXnum_tests)
	G = np.transpose(np.loadtxt(genotype_file, delimiter='\t'))
	# Load in sample overlap data
	Z, num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]
	# Initialize U with K-means and initialize V to zeros
	V = np.zeros((num_latent_factor*2, num_tests))
	intercept = np.zeros((num_samples,num_tests))
	if initialization == 'kmeans':
		U = initialize_sample_loading_matrix_with_kmeans(num_samples, num_latent_factor, Y)
	elif initialization == 'random':
		U = np.random.random(size=(num_samples, num_latent_factor))
	elif initialization == 'fixed':
		U = np.zeros((num_samples, num_latent_factor)) + .1
	elif initialization == 'residual_clustering':
		U = initialize_sample_loading_with_residual_clustering(Y, G, Z, num_latent_factor, output_root)
	else:
		print('ASSUMPTION ERROR: Initialization not implemented')
	print('Initialization complete.')
	#####################################
	# Start iterative optimization process
	######################################
	# Standardize
	G = standardize_each_column_of_matrix(G)
	num_iter = 1
	for itera in range(num_iter):
		print('Iteration ' + str(itera))
		# Update factor matrix (V) with linear mixed model
		old_V = V
		V = update_factor_matrix(Y, G, U, Z, num_samples, num_tests, num_latent_factor, lasso_param_v)	
		pdb.set_trace()	
		# Update loading matrix (U) with l1 penalized linear model
		old_U = U
		U = update_loading_matrix(Y, G, V, Z, intercept, num_samples, num_tests, num_latent_factor, lasso_param_u)
		frob_norm = np.linalg.norm(U - old_U)
		print('L2 norm: ' + str(frob_norm))
	# Save matrices to output
	np.savetxt(output_root + '_U.txt', U, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_V.txt', V, fmt="%s", delimiter='\t')


######################
# Command line args
######################
# File containing indexes of which rna-seq sample came from the same individuals
sample_overlap_file = sys.argv[1]
# Expression matrix used to train factor model
expression_training_file = sys.argv[2]
# Genotype matrix used to test factor model
genotype_training_file = sys.argv[3]
# Expression matrix used to test significance of factors
expression_testing_file = sys.argv[4]
# Genotype matrix used to test significance of factors
genotype_testing_file = sys.argv[5]
# Number of latent factors to model in matrix factorization
num_latent_factors = int(sys.argv[6])
# Stem to use in output files
file_stem = sys.argv[7]
# Directory to save results to
output_dir = sys.argv[8]
lasso_param_u = float(sys.argv[9])
lasso_param_v = float(sys.argv[10])
initialization = sys.argv[11]


#initialization = 'random' # kmeans, random
#genotype_intercept = 'True_penalized' # 'True', 'False', 'True_penalized'
output_root = output_dir + file_stem + '_em_model_lasso_U_' + str(lasso_param_u) + '_lasso_V_' + str(lasso_param_v) + '_initialization_' + initialization
print(output_root)
train_eqtl_factorization_model_em_version(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, lasso_param_u, lasso_param_v, initialization, output_root)

