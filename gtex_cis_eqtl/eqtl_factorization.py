#from __future__ import absolute_import
#rom __future__ import division
#from __future__ import print_function
#from edward.models import Gamma, Poisson, Normal, PointMass, \
#    TransformedDistribution
#from edward.models import PointMass, Empirical, HalfNormal
import edward as ed
import tensorflow as tf
from edward.models import PointMass, Normal
import numpy as np
import os
import sys
import pdb
import math
import edward_model
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import statsmodels.api as sm
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from joblib import Parallel, delayed
import multiprocessing
import cvxpy as cp
from scipy.optimize import minimize


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


def update_factor_matrix_one_test(test_number, Y, G, U, Z, lasso_param):
	# Get slice of data corresponding to this test
	y_test = Y[:, test_number]
	g_test = G[:, test_number]
	# Get U scaled by genotype for this test
	U_scaled = U*g_test[:,None]

	#X = np.hstack((U, U_scaled))
	# Get U for intercept terms

	# Fit linear regression model
	clf = linear_model.Lasso(alpha=lasso_param, fit_intercept=True)
	clf.fit(U_scaled, y_test)
	# Get params of fitted model
	params = np.hstack((clf.intercept_, clf.coef_))
	# Fit lasso model
	#clf = linear_model.Lasso(alpha=lasso_param, fit_intercept=False)
	#clf.fit(X, y_test)
	#params = np.hstack((clf.intercept_, clf.coef_))
	#pdb.set_trace()

	return params

# Update factor matrix (V) with linear mixed model
def update_factor_matrix(Y, G, U, Z, num_samples, num_tests, num_latent_factor, lasso_param):
	num_cores = multiprocessing.cpu_count()
	#print(num_cores)
	#V_arr = Parallel(n_jobs=num_cores)(delayed(update_factor_matrix_one_test)(test_number, Y, G, U, Z) for test_number in range(num_tests))
	genotype_intercept = 'False'
	if genotype_intercept == 'False':
		V_arr = []
		for test_number in range(num_tests):
			V_arr.append(update_factor_matrix_one_test(test_number, Y, G, U, Z, lasso_param))
		#V_arr = Parallel(n_jobs=num_cores)(delayed(update_factor_matrix_one_test)(test_number, Y, G, U, Z, penalty_vector) for test_number in range(num_tests))
		# Convert back to matrix
		V_full = np.transpose(np.asmatrix(V_arr))
		intercept = np.squeeze(np.asarray(V_full[0,:]))
		#intercept_mat = (np.ones((num_samples,1))*np.asmatrix(intercept))
		V = np.asarray(V_full[1:,:])
	return V, np.asarray(intercept)


def update_factor_matrix_edward(Y, G, U, Z, num_samples, num_tests, num_latent_factor, lasso_param_v):
	tf.reset_default_graph()
	sess = tf.InteractiveSession()
	#resid_sd = 1.0
	beta_sd = math.sqrt((1.0)/(lasso_param_v*num_samples))
	# Get number of individuals
	num_individuals = len(np.unique(Z))

	# Set up placeholders for the data inputs.
	ind_ph = tf.placeholder(tf.int32, [num_samples, 1])
	# Set up placeholders for the data inputs.
	genotype = tf.placeholder(tf.float32, [num_samples, num_tests])
	loading = tf.placeholder(tf.float32, [num_samples, num_latent_factor])

	# Set up fixed effects (intercept term)
	mu = tf.get_variable("mu", [num_tests])
	# Set up standard deviation of random effects term
	sigma_ind = tf.sqrt(tf.exp(tf.get_variable("sigma_ind", [num_tests])))
	# Set up standard deviation of residuals term
	sigma_resid = tf.sqrt(tf.exp(tf.get_variable("sigma_resid", [num_tests])))
	# Set up random effects
	eta_ind = Normal(loc=tf.zeros([num_individuals, num_tests]), scale= tf.matmul(tf.ones([num_individuals,num_tests]),tf.matrix_diag(sigma_ind)))


	V = Normal(loc=0.0, scale=beta_sd, sample_shape=[num_latent_factor, num_tests])

	yhat = (tf.multiply(genotype, tf.matmul(loading, V)) + tf.gather_nd(eta_ind, ind_ph) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))

	y = Normal(loc=yhat, scale=tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(sigma_resid)))
	#y = Normal(loc=yhat, scale=resid_sd*tf.ones([num_samples, num_tests]))
	###############################
	## Inference set up
	###############################
	q_ind_s = Normal(loc=tf.get_variable("q_ind_s/loc", [num_individuals, num_tests]),scale=tf.nn.softplus(tf.get_variable("q_ind_s/scale", [num_individuals, num_tests])))

	qV = PointMass(tf.get_variable("qV",[num_latent_factor,num_tests]))


	inference_e = ed.KLqp({eta_ind:q_ind_s}, data={y:Y, ind_ph: Z, genotype: G, loading:U, V:qV})

	inference_m = ed.MAP({V:qV},data={y:Y, ind_ph: Z, genotype: G, eta_ind:q_ind_s, loading:U})
	n_iter=70
	inference_e.initialize()
	num_m_steps_per_iter = 1
	optimizer = tf.train.AdamOptimizer(0.01, epsilon=1.0)
	inference_m.initialize(n_iter=n_iter*num_m_steps_per_iter, optimizer=optimizer)
	tf.global_variables_initializer().run()
	loss = np.empty(n_iter*num_m_steps_per_iter, dtype=np.float32)
	
	counter = 0
	for i in range(n_iter):
		info_dict_e = inference_e.update()
		for j in range(num_m_steps_per_iter):
			#print(j)
			info_dict_m = inference_m.update()
			loss[counter] = info_dict_m["loss"]
			#print(loss[counter])
			counter = counter + 1
		inference_m.print_progress(info_dict_m)
	# Get parameters to pass on
	random_effects_mean = tf.gather_nd(q_ind_s.loc.eval(), Z).eval()
	random_effects_sd = tf.gather_nd(q_ind_s.scale.eval(), Z).eval()
	test_intercept = mu.eval()
	residual_sd = sigma_resid.eval()
	re_sd = sigma_ind.eval()

	new_V = qV.eval()

	#intercept = test_intercept + random_effects_mean
	#variance = (np.square(random_effects_sd) + np.square(resid_sd))
	return new_V, residual_sd, re_sd, test_intercept

def update_loading_matrix_one_sample(sample_num, Y, G, V, Z, mu, weights, lasso_param):
	# Get slice of data corresponding to this sample
	y_test = Y[sample_num, :] - mu
	g_test = G[sample_num, :]

	# Get V scaled by genotype for this sample
	V_scaled = np.transpose(V)*g_test[:,None]

	if lasso_param == 0:
		# Mimick being close to zero
		clf = linear_model.Lasso(alpha=0.0000000000001,positive=True, fit_intercept=False)
		clf.fit(V_scaled, y_test)
	# Fit Lasso model
	else:
		model = sm.WLS(y_test, V_scaled, weights=weights)
		fit = model.fit_regularized(method='elastic_net', alpha=lasso_param, L1_wt=1.0)
		#reg = LinearRegression().fit(V_scaled, y_test, sample_weight=weights)
		#clf = linear_model.Ridge(alpha=lasso_param,positive=True, fit_intercept=False)
		#clf = linear_model.Lasso(alpha=lasso_param, fit_intercept=False)
		#clf.fit(V_scaled, y_test)
	return np.asarray(fit.params)

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
	beta = cp.Variable(X.shape[1],nonneg=True)
		#problem1 = cp.Problem(cp.Minimize(loading_objective_fn_correct(X, y_ind, beta, lasso_param, precision)))
		#problem1.solve()
	problem2 = cp.Problem(cp.Minimize(loading_objective_fn(X, y_ind, beta, lasso_param, precision, X.shape[0]/len(sample_indices))))
	problem2.solve()

	new_U = beta.value
	return new_U

# Update loading matrix (U) with l1 penalized linear model
def update_loading_matrix(Y, G, V, Z, intercept, residual_sd, re_sd, num_samples, num_tests, num_latent_factor, lasso_param):
	# Initialize output matrix
	U = np.zeros((num_samples, num_latent_factor))
	# Scale Y by intercept
	Y_scaled = Y - intercept
	# Loop through individuals
	individuals = np.unique(Z)
	# Get transpose of V
	V_t = np.transpose(V)

	U_indi = []

	for individual in individuals:
		#print(individual)
		# Get indixes of samples corresponding to this individual
		sample_indices = np.where(Z == individual)[0]
		# Get response vector (Y) for this individual
		# Is ordered (test1, sample1), (test1, sample2), (test1, sampleI), ..., (testT,sampleI)
		y_ind = np.transpose(Y_scaled[sample_indices,:]).reshape(-1)
		# Get known precision matrix for this individual
		precision = extract_precision_matrix_for_individual(residual_sd, re_sd, sample_indices, num_tests)
		# Get feature matrix
		X = get_loading_feature_matrix_for_one_individual(sample_indices, num_tests, V_t, G)

		#beta = cp.Variable(X.shape[1],nonneg=True)
		#beta = cp.Variable(X.shape[1])
		num_samples_per_test = X.shape[0]
		#problem1 = cp.Problem(cp.Minimize(loading_objective_fn_correct(X, y_ind, beta, lasso_param, precision)))
		#problem1.solve()
		#problem2 = cp.Problem(cp.Minimize(loading_objective_fn(X, y_ind, beta, lasso_param, precision, X.shape[0])))
		#problem2.solve()

		#new_U = beta.value

		new_U = np.matmul(np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(X.T, precision), X) + np.eye(X.shape[1])*lasso_param*num_samples_per_test), X.T),precision), y_ind)

		
		for i, index in enumerate(sample_indices):
			col_start = num_latent_factor*i
			col_end = num_latent_factor*(i+1)
			U[index,:] = new_U[col_start:col_end]
		#pdb.set_trace()
		# Update model
		#U[sample_num,:] = update_loading_matrix_one_sample(sample_num, Y, G, V, Z, intercept[sample_num, :], 1.0/variance[sample_num,:], lasso_param)
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
def train_eqtl_factorization_model_em_version(sample_overlap_file, expression_file, genotype_file, num_latent_factor, lasso_param_u, lasso_param_v, initialization, output_root):
	# Load in expression data (dimension: num_samplesXnum_tests)
	Y = np.transpose(np.loadtxt(expression_file, delimiter='\t'))
	# Load in genotype data (dimension: num_samplesXnum_tests)
	G = np.transpose(np.loadtxt(genotype_file, delimiter='\t'))
	# Load in sample overlap data
	Z, Z2,  num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]
	# Initialize U with K-means and initialize V to zeros
	V = np.zeros((num_latent_factor, num_tests))
	intercept = np.zeros((num_samples,num_tests))
	if initialization == 'kmeans':
		U = initialize_sample_loading_matrix_with_kmeans(num_samples, num_latent_factor, Y)
	elif initialization.startswith('random'):
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
	num_iter = 100
	mse_list = []
	for itera in range(num_iter):
		print('Iteration ' + str(itera))
		# Update factor matrix (V) with linear mixed model
		old_V = V
		#V, intercept = update_factor_matrix(Y, G, U, Z, num_samples, num_tests, num_latent_factor, lasso_param_v)
		V, residual_sd, re_sd, intercept = update_factor_matrix_edward(Y, G, U, Z, num_samples, num_tests, num_latent_factor, lasso_param_v)
		# Update loading matrix (U) with l1 penalized linear model
		old_U = U
		U = update_loading_matrix(Y, G, V, Z2, intercept, residual_sd, re_sd, num_samples, num_tests, num_latent_factor, lasso_param_u)
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
lasso_param_u = float(sys.argv[9])
lasso_param_v = float(sys.argv[10])
initialization = sys.argv[11]


#initialization = 'random' # kmeans, random
#genotype_intercept = 'True_penalized' # 'True', 'False', 'True_penalized'
output_root = output_dir + file_stem + '_em_model_lasso_U_' + str(lasso_param_u) + '_lasso_V_' + str(lasso_param_v) + '_initialization_' + initialization
print(output_root)
train_eqtl_factorization_model_em_version(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, lasso_param_u, lasso_param_v, initialization, output_root)

