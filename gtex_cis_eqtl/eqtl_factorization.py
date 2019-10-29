import numpy as np
import os
import sys
import pdb
import math
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import statsmodels.api as sm
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from joblib import Parallel, delayed
import multiprocessing
from scipy.optimize import minimize

import sm_source_code


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
	return Z, np.asarray(Z1), int(num_individuals)

def initialize_sample_loading_matrix_with_kmeans(num_samples, num_latent_factor, Y):
	# Sample loading matrix
	U = np.zeros((num_samples, num_latent_factor))
	# Fit kmeans model
	kmeans = KMeans(n_clusters=num_latent_factor, random_state=0).fit(Y)
	# Fill in sample loading matrix with kmeans assignments
	for sample_index, kmeans_assignment in enumerate(kmeans.labels_):
		U[sample_index, kmeans_assignment] = 1.0
	return U


def update_factor_matrix_one_test(test_number, model_name, Y, G, U, Z, lasso_param):
	# Get slice of data corresponding to this test
	y_test = Y[:, test_number]
	g_test = G[:, test_number]
	# Get U scaled by genotype for this test
	U_scaled = U*g_test[:,None]



	##################
	if model_name == 'almm':
		covariates = np.hstack((np.ones((U_scaled.shape[0],1)), U_scaled))
		alpha = np.zeros(covariates.shape[1]) + lasso_param*len(y_test)
		alpha[0] = 0.0
		model = sm_source_code.MixedLM(y_test, covariates, Z)
		if lasso_param == 0.0:
			lmm_fit =  model.fit()
		else:
			lmm_fit = model.fit_regularized(alpha = alpha)
		params = np.hstack((lmm_fit.fe_params, np.sqrt(lmm_fit.cov_re[0,0]), np.sqrt(lmm_fit.scale)))
	elif model_name == 'alm':
		clf = linear_model.ElasticNet(alpha=lasso_param, l1_ratio=.8, fit_intercept=True)
		clf.fit(U_scaled, y_test)
		params = np.hstack((clf.intercept_, clf.coef_, None, None))
	elif model_name == 'alm_genotype_intercept':
		#clf = linear_model.Lasso(alpha=lasso_param, fit_intercept=True)
		#clf.fit(U_scaled, y_test)
		#params = np.hstack((clf.intercept_, clf.coef_, None, None))
		covariates = np.hstack((np.ones((U_scaled.shape[0],1)), np.transpose(np.asmatrix(g_test)), U_scaled))
		model = sm.OLS(y_test, covariates)
		alpha_param = np.zeros(covariates.shape[1]) + lasso_param
		alpha_param[0] = 0
		alpha_param[1] = 0
		fit = model.fit_regularized(method='elastic_net', alpha=alpha_param,L1_wt=1.0)
		params = np.hstack((fit.params, None, None))
	################


	return params

# Update factor matrix (V) with linear mixed model
def update_factor_matrix(model_name, Y, G, U, Z, num_samples, num_tests, num_latent_factor, lasso_param):
	#num_cores = multiprocessing.cpu_count()
	num_cores=15
	parrallel = False

	V_arr = []
	if parrallel == True:
		V_arr = Parallel(n_jobs=num_cores)(delayed(update_factor_matrix_one_test)(test_number, model_name, Y, G, U, Z, lasso_param) for test_number in range(num_tests))
	else:
		for test_number in range(num_tests):
			V_arr.append(update_factor_matrix_one_test(test_number, model_name, Y, G, U, Z, lasso_param))
	# Convert back to matrix
	if model_name == 'alm_genotype_intercept':
		V_full = np.transpose(np.asmatrix(V_arr))
		intercept = np.squeeze(np.asarray(V_full[0,:]))
		V = np.asarray(V_full[1:(2+num_latent_factor),:])
		re_sd = np.squeeze(np.asarray(V_full[(num_latent_factor + 2),:]))
		residual_sd = np.squeeze(np.asarray(V_full[(num_latent_factor + 3),:]))
	else:
		V_full = np.transpose(np.asmatrix(V_arr))
		intercept = np.squeeze(np.asarray(V_full[0,:]))
		V = np.asarray(V_full[1:(1+num_latent_factor),:])
		re_sd = np.squeeze(np.asarray(V_full[(num_latent_factor + 1),:]))
		residual_sd = np.squeeze(np.asarray(V_full[(num_latent_factor + 2),:]))
	return V, residual_sd, re_sd, intercept



def update_loading_matrix_one_sample(sample_num, Y_scaled, G, V, lasso_param):
	# Get slice of data corresponding to this sample
	y_test = Y_scaled[sample_num, :] 
	g_test = G[sample_num, :]

	# Get V scaled by genotype for this sample
	V_scaled = np.transpose(V)*g_test[:,None]

	clf = linear_model.ElasticNet(alpha=lasso_param, l1_ratio=1.0, positive=True, fit_intercept=False)
	clf.fit(V_scaled, y_test)

	return np.asarray(clf.coef_)

# Get known covariance matrix for this individual
def extract_precision_matrix_for_individual(residual_sd, re_sd, sample_indices, num_tests):
	residual_variance = np.square(residual_sd)
	re_variance = np.square(re_sd)
	# Get number of samples this individual has
	num_samples = len(sample_indices)
	# Initialize covariance matrix
	precision = np.zeros((num_samples*num_tests, num_samples*num_tests))
	for test_number in range(num_tests):
		starting_point = test_number*num_samples
		ending_point = (test_number+1)*num_samples
		cov = np.ones((num_samples,num_samples))*re_variance[test_number]
		for pos in range(num_samples):
			cov[pos,pos] = cov[pos,pos] + residual_variance[test_number]
		precision[starting_point:ending_point,starting_point:ending_point] = np.linalg.inv(cov)
	return precision

# Get feature matrix
def get_loading_feature_matrix_for_one_individual(sample_indices, num_tests, V_t, G):
	# Num factors
	num_latent_factors = V_t.shape[1]
	# Get number of samples this individual has
	num_samples = len(sample_indices)
	# Get predictor vector (G) for this individual
	g_ind = np.transpose(G[sample_indices,:]).reshape(-1)
	# Initialize feature matrix
	feat = np.zeros((num_tests*num_samples, num_latent_factors*num_samples))
	for test_number in range(num_tests):
		for sample_num in range(num_samples):
			row_number = test_number*num_samples + sample_num
			col_start = num_latent_factors*sample_num
			col_end = (num_latent_factors)*(sample_num+1)
			feat[row_number,col_start:col_end] = V_t[test_number,:]*g_ind[row_number]
	return feat


def loading_objective_fn(X, y, beta, lasso_param, precision, nn):
	val = cp.quad_form(y - cp.matmul(X,beta), precision)
	return val/float(nn) + lasso_param*cp.pnorm(beta, p=2)**2
	#return val/nn + lasso_param*cp.norm1(beta)


def loading_matrix_one_sample(individual, Z, Y_scaled, residual_sd, re_sd, num_tests, V_t, G, lasso_param):
	# Get indixes of samples corresponding to this individual
	sample_indices = np.where(Z == individual)[0]
	# Get response vector (Y) for this individual
	# Is ordered (test1, sample1), (test1, sample2), (test1, sampleI), ..., (testT,sampleI)
	y_ind = np.transpose(Y_scaled[sample_indices,:]).reshape(-1)
	# Get known precision matrix for this individual
	precision = extract_precision_matrix_for_individual(residual_sd, re_sd, sample_indices, num_tests)
	# Get feature matrix
	X = get_loading_feature_matrix_for_one_individual(sample_indices, num_tests, V_t, G)

	num_samples_per_test = X.shape[0]

	new_U_unpenalized = np.matmul(np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(X.T, precision), X)), X.T),precision), y_ind)

	fe_params = new_U_unpenalized.copy()

	alpha = np.zeros(len(new_U_unpenalized)) + lasso_param*len(y_ind)
	ceps=1e-4
	ptol=1e-6
	maxit=200
	for itr in range(maxit):
		fe_params_s = fe_params.copy()
		for j in range(len(fe_params)):
			if abs(fe_params[j]) < ceps:
				continue
			# The residuals
			fe_params[j] = 0.
			expval = np.dot(X, fe_params)
			resid_all = y_ind - expval

			uu = np.dot(precision, X[:,j])

			aa = np.dot(uu, X[:,j])
			bb = -2 * np.dot(uu, resid_all)
			pwt1 = alpha[j]
			if bb > pwt1:
				fe_params[j] = -(bb - pwt1) / (2 * aa)
			elif bb < -pwt1:
				fe_params[j] = -(bb + pwt1) / (2 * aa)
		if np.abs(fe_params_s - fe_params).max() < ptol:
			break
	new_U = fe_params

	return new_U

# Update loading matrix (U) with l1 penalized linear model
def update_loading_matrix(model_name, Y, G, V, Z, intercept, residual_sd, re_sd, num_samples, num_tests, num_latent_factor, lasso_param):
	# Scale Y by intercept
	Y_scaled = Y - intercept
	# Loop through individuals
	individuals = np.unique(Z)
	# Get transpose of V
	V_t = np.transpose(V)
	# Parrallelization stuff
	num_cores = multiprocessing.cpu_count()
	parrallel = False
	U_indi = []

	if model_name == 'almm':
		U = np.zeros((num_samples, num_latent_factor))
		if parrallel == True:
			U_indi = Parallel(n_jobs=num_cores)(delayed(loading_matrix_one_sample)(individual, Z, Y_scaled, residual_sd, re_sd, num_tests, V_t, G, lasso_param) for individual in individuals)
		else:
			for individual in individuals:
				U_indi.append(loading_matrix_one_sample(individual, Z, Y_scaled, residual_sd, re_sd, num_tests, V_t, G, lasso_param))
		for itera, individual in enumerate(individuals):
			sample_indices = np.where(Z == individual)[0]
			U_vec = U_indi[itera]
			for i, index in enumerate(sample_indices):
				col_start = num_latent_factor*i
				col_end = num_latent_factor*(i+1)
				U[index,:] = U_vec[col_start:col_end]
	elif model_name == 'alm':
		if parrallel == True:
			U_indi = Parallel(n_jobs=num_cores)(delayed(update_loading_matrix_one_sample)(sample_num, Y_scaled, G, V, lasso_param) for sample_num in range(num_samples))
		else:
			for sample_num in range(num_samples):
				U_indi.append(update_loading_matrix_one_sample(sample_num, Y_scaled, G, V, lasso_param))
		U = np.asarray(U_indi)
	elif model_name == 'alm_genotype_intercept':
		Y_scaled = Y_scaled - (V[0,:]*G)
		if parrallel == True:
			U_indi = Parallel(n_jobs=num_cores)(delayed(update_loading_matrix_one_sample)(sample_num, Y_scaled, G, V[1:,:], lasso_param) for sample_num in range(num_samples))
		else:
			for sample_num in range(num_samples):
				U_indi.append(update_loading_matrix_one_sample(sample_num, Y_scaled, G, V[1:,:], lasso_param))
		U = np.asarray(U_indi)
	return U

def compute_rmse(resid_matrix):
	return np.sqrt(np.mean(np.square(resid_matrix)))

def compute_mse(resid_matrix):
	return np.mean(np.square(resid_matrix))

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

	# Compute RMSE on training data
	mse = compute_mse(resid_matrix)
	print('LM MSE: ' + str(mse))
	np.savetxt(output_root + '_lm_mse.txt', [mse], fmt="%s", delimiter='\t')

	return U_init

def standardize_each_column_of_matrix(G):
	num_cols = G.shape[1]
	for col_num in range(num_cols):
		G[:,col_num] = (G[:,col_num] - np.mean(G[:,col_num]))/np.std(G[:,col_num])
	return G

def make_eqtl_factorization_predictions(G, U, V, intercept):
	pred_y = np.dot(U,V)*G - intercept
	return pred_y

#######################
# Part 1: Train eqtl factorization model
########################
def train_eqtl_factorization_model(model_name, sample_overlap_file, expression_file, genotype_file, num_latent_factor, lasso_param_u, lasso_param_v, initialization, output_root):
	############################
	# Load in data
	############################
	# Load in expression data (dimension: num_samplesXnum_tests)
	Y = np.transpose(np.loadtxt(expression_file, delimiter='\t'))
	# Load in genotype data (dimension: num_samplesXnum_tests)
	G = np.transpose(np.loadtxt(genotype_file, delimiter='\t'))
	# Load in sample overlap data
	Z, Z2,  num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]

	############################
	# Initialize U matrix
	############################
	V = np.zeros((num_latent_factor, num_tests))
	if initialization == 'random':
		U = np.random.random(size=(num_samples, num_latent_factor))
	else:
		print('ASSUMPTION ERROR: Initialization not implemented')
	print('Initialization complete.')

	#####################################
	# Start iterative optimization process
	######################################
	# Standardize
	G = standardize_each_column_of_matrix(G)
	num_iter = 200
	#mse_list = []

	for itera in range(num_iter):
		print('Iteration ' + str(itera))
		
		# Update factor matrix (V)
		old_V = V
		V, residual_sd, re_sd, intercept = update_factor_matrix(model_name, Y, G, U, Z2, num_samples, num_tests, num_latent_factor, lasso_param_v)

		# Update loading matrix (U)
		old_U = U
		U = update_loading_matrix(model_name, Y, G, V, Z2, intercept, residual_sd, re_sd, num_samples, num_tests, num_latent_factor, lasso_param_u)
		frob_norm = np.linalg.norm(U - old_U)
		print('L2 norm: ' + str(frob_norm))
		print(U)
		# Make predictions
		#Y_hat = make_eqtl_factorization_predictions(G, U, V, intercept)
		#mse = compute_mse(Y_hat - Y)
		#mse_list.append(mse)
		#print('eQTL Factorization MSE: ' + str(mse))
	# Save matrices to output
	np.savetxt(output_root + '_U.txt', U, fmt="%s", delimiter='\t')
	np.savetxt(output_root + '_V.txt', V, fmt="%s", delimiter='\t')
	#np.savetxt(output_root + '_eqtl_factorization_mse.txt', np.asarray(mse_list), fmt="%s", delimiter='\n')


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
# Lasso penalty on elements of U matrix
lasso_param_u = float(sys.argv[9])
# Lasso penalty on elements of V matrix
lasso_param_v = float(sys.argv[10])
# How to initialize U matrix
initialization = sys.argv[11]
# Random seed to use for initialization
seed = int(sys.argv[12])
# Name of model to be used for optimization 
model_name = sys.argv[13]

# What to save output files to
output_root = output_dir + file_stem 

# Set random seed
np.random.seed(seed)


#########################
# Train model
#########################
train_eqtl_factorization_model(model_name, sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, lasso_param_u, lasso_param_v, initialization, output_root)

