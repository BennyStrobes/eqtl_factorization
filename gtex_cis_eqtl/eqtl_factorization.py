#from __future__ import absolute_import
#rom __future__ import division
#from __future__ import print_function
#from edward.models import Gamma, Poisson, Normal, PointMass, \
#    TransformedDistribution
from edward.models import PointMass, HalfNormal, Normal, Bernoulli
import edward as ed
import tensorflow as tf
import numpy as np 
import os
import sys
import pdb
from sklearn.cluster import KMeans
import statsmodels.api as sm
from sklearn import linear_model
from joblib import Parallel, delayed
import multiprocessing
#import pystan


# Load in pystan optimizizer
#UPDATE_V = pystan.StanModel(file = "update_v.stan")
#UPDATE_U = pystan.StanModel(file = "update_u.stan")




# Load in sample overlap data
def load_in_sample_overlap_data(sample_overlap_file):
	f = open(sample_overlap_file)
	Z = []
	temp = []
	for line in f:
		line = line.rstrip()
		Z.append([int(line)])
		temp.append(int(line))
	f.close()
	num_individuals = max(temp) + 1
	return Z, int(num_individuals)

# Load in sample overlap data
def load_in_sample_overlap_data_v2(sample_overlap_file):
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


# Return EDWARD lognormal variational objective
def lognormal_q(shape, name=None):
  with tf.variable_scope(name, default_name="lognormal_q"):
    min_scale = 1e-5
    loc = tf.get_variable("loc", shape)
    scale = tf.get_variable(
        "scale", shape, initializer=tf.random_normal_initializer(stddev=0.1))
    rv = TransformedDistribution(
        distribution=Normal(loc, tf.maximum(tf.nn.softplus(scale), min_scale)),
        bijector=tf.contrib.distributions.bijectors.Exp())
    return rv

def gamma_q(obj_shape, name=None):
  # Parameterize Gamma q's via shape and scale, with softplus unconstraints.
  with tf.variable_scope(name, default_name="gamma_q"):
    min_shape = 1e-3
    min_scale = 1e-5
    shape = tf.get_variable(
        "shape", obj_shape,
        initializer=tf.random_normal_initializer(mean=0.5, stddev=0.1))
    scale = tf.get_variable(
        "scale", obj_shape, initializer=tf.random_normal_initializer(stddev=0.1))
    rv = Gamma(tf.maximum(tf.nn.softplus(shape), min_shape),
               tf.maximum(1.0 / tf.nn.softplus(scale), 1.0 / min_scale))
    return rv

#######################
# Part 1: Train eqtl factorization model
########################
def train_eqtl_factorization_model(sample_overlap_file, expression_file, genotype_file, num_latent_factors, sparsity_parameter_u, sparsity_parameter_v, output_root):
	# Load in expression data
	Y = np.transpose(np.loadtxt(expression_file, delimiter='\t'))
	# Load in genotype data
	G = np.transpose(np.loadtxt(genotype_file, delimiter='\t'))
	# Load in sample overlap data
	Z, num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]

	###############################
	## MODEL
	# Y ~ 1 + (UXV)G 
	###############################
	# Set up placeholders for the data inputs.
	#ind_ph = tf.placeholder(tf.int32, [num_samples, 1])
	# Set up placeholders for the data inputs.
	genotype = tf.placeholder(tf.float32, [num_samples, num_tests])
	# Set up fixed effects (intercept term)
	mu = tf.get_variable("mu", [num_tests])
	# Set up standard deviation of random effects term
	#sigma_ind = tf.sqrt(tf.exp(tf.get_variable("sigma_ind", [num_tests])))
	# Set up standard deviation of residuals term
	sigma_resid = tf.sqrt(tf.exp(tf.get_variable("sigma_resid", [num_tests])))
	# Set up random effects
	# eta_ind = Normal(loc=tf.zeros([num_individuals, num_tests]), scale= tf.matmul(tf.ones([num_individuals,num_tests]),tf.matrix_diag(sigma_ind)))

	# higher values of sparsity parameter result in a more sparse solution
	#U = tf.nn.softplus(tf.get_variable("U", [num_samples, num_latent_factors]))
	#U = tf.get_variable("U", [num_samples, num_latent_factors])
	V = Normal(loc=0.0, scale=1.0/(2.0*sparsity_parameter_v), sample_shape=[num_latent_factors, num_tests])
	#U = HalfNormal(scale=1.0/(2.0*sparsity_parameter_u), sample_shape=[num_samples, num_latent_factors])

	U_mask = Bernoulli(logits=0.0, sample_shape=[num_samples, num_latent_factors])
	U_raw = HalfNormal(scale=1.0/(2.0*sparsity_parameter_u), sample_shape=[num_samples, num_latent_factors])
	U = tf.to_float(U_mask) * U_raw

	#U = Normal(loc=0.0, scale=1.0, sample_shape=[num_samples, num_latent_factors])
	#U = Gamma(concentration=1.0, rate=sparsity_parameter, sample_shape=[num_samples, num_latent_factors])
  	#V = Normal(loc=0.0, scale=1.0, sample_shape=[num_latent_factors, num_tests])
  	#U = ed.to_simplex(tf.get_variable("U", [num_samples, num_latent_factors-1]))

  	#U2 = tf.ones([num_samples,1])
  	#V2 = tf.get_variable("V2", [1, num_tests])
	
	#yhat = (tf.multiply(genotype, tf.matmul(U, V)) + tf.multiply(genotype, tf.matmul(U2, V2)) + tf.gather_nd(eta_ind, ind_ph) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))
	#yhat = (tf.multiply(genotype, tf.matmul(U, V)) + tf.gather_nd(eta_ind, ind_ph) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))
	yhat = (tf.multiply(genotype, tf.matmul(U, V)) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))
	y = Normal(loc=yhat, scale=tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(sigma_resid)))

	###############################
	## Inference set up
	###############################
	#q_ind_s = Normal(loc=tf.get_variable("q_ind_s/loc", [num_individuals, num_tests]),scale=tf.nn.softplus(tf.get_variable("q_ind_s/scale", [num_individuals, num_tests])))
	
	#qU = PointMass(tf.nn.softplus(tf.get_variable("qU",[num_samples,num_latent_factors])))
	qU = PointMass(tf.nn.softplus(tf.Variable(tf.random_normal([num_samples, num_latent_factors]))))
	qU_mask = PointMass(tf.Variable(tf.zeros([num_samples, num_latent_factors])))
	#qU = PointMass(tf.nn.softplus(tf.Variable(tf.zeros([num_samples,num_latent_factors]))))

  	qV = Normal(loc=tf.get_variable("qV/loc", [num_latent_factors, num_tests]),scale=tf.nn.softplus(tf.get_variable("qV/scale", [num_latent_factors, num_tests])))

	#latent_vars = {eta_ind: q_ind_s, U: qU}

	#data = {y: Y, ind_ph: Z, genotype: G}

	#inference = ed.KLqp(latent_vars, data)
	#tf.global_variables_initializer().run()
	#inference.run(n_iter=400)
	n_iter = 10
	inference_e = ed.KLqp({V:qV}, data={y:Y, genotype: G, U_raw:qU, U_mask:qU_mask})

	inference_m = ed.MAP({U_raw:qU, U_mask:qU_mask},data={y:Y, genotype: G, V: qV})

	inference_e.initialize()
	num_m_steps_per_iter = 1
	inference_m.initialize(n_iter=n_iter*num_m_steps_per_iter, optimizer="rmsprop")
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

	pdb.set_trace()
	#########################
	# Save data to output
	#########################
	output_file = output_root + '_qU_mean.txt'
	#qU_mean = np.exp(qU.distribution.loc.eval())
	qU_mean = qU.mean().eval()
	np.savetxt(output_file,qU_mean, fmt="%s",delimiter='\t')

	#output_file = output_root + '_qU_sd.txt'
	#qU_sd = tf.maximum(tf.nn.softplus(qU.distribution.scale), 1e-5).eval()
	#np.savetxt(output_file, qU_sd, fmt="%s",delimiter='\t')

	output_file = output_root + '_V.txt'
	np.savetxt(output_file, qV.mean().eval(), fmt="%s",delimiter='\t')

	output_file = output_root + '_mu.txt'
	np.savetxt(output_file, mu.eval(), fmt="%s",delimiter='\n')

	output_file = output_root + '_sigma_resid.txt'
	np.savetxt(output_file, sigma_resid.eval(), fmt="%s",delimiter='\n')


#######################
# Part 1: Train eqtl factorization model
########################
def train_eqtl_factorization_model_v2(sample_overlap_file, expression_file, genotype_file, num_latent_factors, sparsity_parameter, output_root, n_iter):
	# Load in expression data
	Y = np.transpose(np.loadtxt(expression_file, delimiter='\t'))
	# Load in genotype data
	G = np.transpose(np.loadtxt(genotype_file, delimiter='\t'))
	# Load in sample overlap data
	Z, num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]

	
	###############################
	## MODEL
	# Y ~ 1 + (UXV)G (1|individual) 
	###############################
	# Set up placeholders for the data inputs.
	ind_ph = tf.placeholder(tf.int32, [num_samples, 1])
	# Set up placeholders for the data inputs.
	genotype = tf.placeholder(tf.float32, [num_samples, num_tests])
	# Set up fixed effects (intercept term)
	mu = tf.get_variable("mu", [num_tests])
	# Set up standard deviation of random effects term
	sigma_ind = tf.sqrt(tf.exp(tf.get_variable("sigma_ind", [num_tests])))
	# Set up standard deviation of residuals term
	sigma_resid = tf.sqrt(tf.exp(tf.get_variable("sigma_resid", [num_tests])))
	# Set up random effects
	eta_ind = Normal(loc=tf.zeros([num_individuals, num_tests]), scale= tf.matmul(tf.ones([num_individuals,num_tests]),tf.matrix_diag(sigma_ind)))

	# higher values of sparsity parameter result in a more sparse solution
	#U = tf.nn.softplus(tf.get_variable("U", [num_samples, num_latent_factors]))
	#U = tf.get_variable("U", [num_samples, num_latent_factors])
	#V = tf.get_variable("V", [num_latent_factors, num_tests])
	#U = Normal(loc=0.0, scale=1.0, sample_shape=[num_samples, num_latent_factors])
	U = Gamma(concentration=1.0, rate=sparsity_parameter, sample_shape=[num_samples, num_latent_factors])
  	V = Normal(loc=0.0, scale=5.0, sample_shape=[num_latent_factors, num_tests])
  	#U = ed.to_simplex(tf.get_variable("U", [num_samples, num_latent_factors-1]))

  	#U2 = tf.ones([num_samples,1])
  	#V2 = tf.get_variable("V2", [1, num_tests])
	
	#yhat = (tf.multiply(genotype, tf.matmul(U, V)) + tf.multiply(genotype, tf.matmul(U2, V2)) + tf.gather_nd(eta_ind, ind_ph) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))
	yhat = (tf.multiply(genotype, tf.matmul(U, V)) + tf.gather_nd(eta_ind, ind_ph) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))

	y = Normal(loc=yhat, scale=tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(sigma_resid)))

	###############################
	## Inference set up
	###############################
	q_ind_s = Normal(loc=tf.get_variable("q_ind_s/loc", [num_individuals, num_tests]),scale=tf.nn.softplus(tf.get_variable("q_ind_s/scale", [num_individuals, num_tests])))
	#q_ind_s = Empirical(tf.nn.softplus(tf.Variable(tf.random_normal([num_individuals, num_tests]))))
	#qU = PointMass(tf.nn.softplus(tf.Variable(tf.random_normal([num_samples,num_latent_factors]))))
	qU = PointMass(tf.nn.softplus(tf.get_variable("qU",[num_samples,num_latent_factors])))
	#qV = Empirical(tf.nn.softplus(tf.Variable(tf.random_normal([num_latent_factors,num_tests]))))
	#qU = gamma_q(U.shape)
	#qU = Exponential(rate=tf.nn.softplus(tf.Variable(tf.random_normal([D, N]))))
	#qU = Normal(loc=tf.get_variable("qU/loc", [num_samples, num_latent_factors]),
     #         scale=tf.nn.softplus(
      #            tf.get_variable("qU/scale", [num_samples, num_latent_factors])))
  	qV = Normal(loc=tf.get_variable("qV/loc", [num_latent_factors, num_tests]),scale=tf.nn.softplus(tf.get_variable("qV/scale", [num_latent_factors, num_tests])))

	#latent_vars = {eta_ind: q_ind_s, U: qU}

	#data = {y: Y, ind_ph: Z, genotype: G}

	#inference = ed.KLqp(latent_vars, data)
	#tf.global_variables_initializer().run()
	#inference.run(n_iter=400)

	inference_e = ed.KLqp({V:qV, eta_ind:q_ind_s}, data={y:Y, ind_ph: Z, genotype: G, U:qU})

	inference_m = ed.MAP({U:qU},data={y:Y, ind_ph: Z, genotype: G, V: qV, eta_ind:q_ind_s})

	inference_e.initialize()
	num_m_steps_per_iter = 1
	inference_m.initialize(n_iter=n_iter*num_m_steps_per_iter, optimizer="rmsprop")
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

	#########################
	# Save data to output
	#########################
	output_file = output_root + '_qU_mean.txt'
	#qU_mean = np.exp(qU.distribution.loc.eval())
	qU_mean = qU.mean().eval()
	np.savetxt(output_file,qU_mean, fmt="%s",delimiter='\t')

	#output_file = output_root + '_qU_sd.txt'
	#qU_sd = tf.maximum(tf.nn.softplus(qU.distribution.scale), 1e-5).eval()
	#np.savetxt(output_file, qU_sd, fmt="%s",delimiter='\t')

	output_file = output_root + '_V.txt'
	np.savetxt(output_file, qV.mean().eval(), fmt="%s",delimiter='\t')

	output_file = output_root + '_re_mean.txt'
	np.savetxt(output_file, q_ind_s.mean().eval(), fmt="%s",delimiter='\t')

	output_file = output_root + '_mu.txt'
	np.savetxt(output_file, mu.eval(), fmt="%s",delimiter='\n')
	
	output_file = output_root + '_sigma_ind.txt'
	np.savetxt(output_file, sigma_ind.eval(), fmt="%s",delimiter='\n')

	output_file = output_root + '_sigma_resid.txt'
	np.savetxt(output_file, sigma_resid.eval(), fmt="%s",delimiter='\n')






#######################
# Part 1: Train eqtl factorization model
########################
def train_eqtl_factorization_model_v4(sample_overlap_file, expression_file, genotype_file, num_latent_factors, scale_prior, output_root, n_iter):
	# Load in expression data
	Y = np.transpose(np.loadtxt(expression_file, delimiter='\t'))
	# Load in genotype data
	G = np.transpose(np.loadtxt(genotype_file, delimiter='\t'))
	# Load in sample overlap data
	Z, num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]

	
	###############################
	## MODEL
	# Y ~ 1 + (UXV)G (1|individual) 
	###############################
	# Set up placeholders for the data inputs.
	ind_ph = tf.placeholder(tf.int32, [num_samples, 1])
	# Set up placeholders for the data inputs.
	genotype = tf.placeholder(tf.float32, [num_samples, num_tests])
	# Set up fixed effects (intercept term)
	mu = tf.get_variable("mu", [num_tests])
	# Set up standard deviation of random effects term
	sigma_ind = tf.sqrt(tf.exp(tf.get_variable("sigma_ind", [num_tests])))
	# Set up standard deviation of residuals term
	sigma_resid = tf.sqrt(tf.exp(tf.get_variable("sigma_resid", [num_tests])))
	# Set up random effects
	eta_ind = Normal(loc=tf.zeros([num_individuals, num_tests]), scale= tf.matmul(tf.ones([num_individuals,num_tests]),tf.matrix_diag(sigma_ind)))

	# higher values of sparsity parameter result in a more sparse solution
	#U = tf.nn.softplus(tf.get_variable("U", [num_samples, num_latent_factors]))
	#U = tf.get_variable("U", [num_samples, num_latent_factors])
	#V = tf.get_variable("V", [num_latent_factors, num_tests])
	#U = Normal(loc=0.0, scale=1.0, sample_shape=[num_samples, num_latent_factors])
	U = HalfNormal(scale=scale_prior, sample_shape=[num_samples, num_latent_factors])
  	V = Normal(loc=0.0, scale=5.0, sample_shape=[num_latent_factors, num_tests])
  	#U = ed.to_simplex(tf.get_variable("U", [num_samples, num_latent_factors-1]))

  	U2 = tf.ones([num_samples,1])
  	V2 = Normal(loc=0.0, scale=5.0, sample_shape=[1, num_tests])
	
	yhat = (tf.multiply(genotype, tf.matmul(U, V)) + tf.multiply(genotype, tf.matmul(U2, V2)) + tf.gather_nd(eta_ind, ind_ph) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))
	#yhat = (tf.multiply(genotype, tf.matmul(U, V)) + tf.gather_nd(eta_ind, ind_ph) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))

	y = Normal(loc=yhat, scale=tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(sigma_resid)))

	###############################
	## Inference set up
	###############################
	q_ind_s = Normal(loc=tf.get_variable("q_ind_s/loc", [num_individuals, num_tests]),scale=tf.nn.softplus(tf.get_variable("q_ind_s/scale", [num_individuals, num_tests])))
	
	#qU = PointMass(tf.nn.softplus(tf.get_variable("qU",[num_samples,num_latent_factors])))
	qU = PointMass(tf.nn.softplus(tf.Variable(tf.random_normal([num_samples,num_latent_factors]))))
	#qU = PointMass(tf.nn.softplus(tf.Variable(tf.zeros([num_samples,num_latent_factors]))))

  	qV = Normal(loc=tf.get_variable("qV/loc", [num_latent_factors, num_tests]),scale=tf.nn.softplus(tf.get_variable("qV/scale", [num_latent_factors, num_tests])))

  	qV2 = Normal(loc=tf.get_variable("qV2/loc", [1, num_tests]),scale=tf.nn.softplus(tf.get_variable("qV2/scale", [1, num_tests])))

	#latent_vars = {eta_ind: q_ind_s, U: qU}

	#data = {y: Y, ind_ph: Z, genotype: G}

	#inference = ed.KLqp(latent_vars, data)
	#tf.global_variables_initializer().run()
	#inference.run(n_iter=400)

	inference_e = ed.KLqp({V:qV, eta_ind:q_ind_s, V2:qV2}, data={y:Y, ind_ph: Z, genotype: G, U:qU})

	inference_m = ed.MAP({U:qU},data={y:Y, ind_ph: Z, genotype: G, V: qV, eta_ind:q_ind_s, V2:qV2})

	inference_e.initialize()
	num_m_steps_per_iter = 1
	inference_m.initialize(n_iter=n_iter*num_m_steps_per_iter, optimizer="rmsprop")
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

	#########################
	# Save data to output
	#########################
	output_file = output_root + '_qU_mean.txt'
	#qU_mean = np.exp(qU.distribution.loc.eval())
	qU_mean = qU.mean().eval()
	np.savetxt(output_file,qU_mean, fmt="%s",delimiter='\t')

	#output_file = output_root + '_qU_sd.txt'
	#qU_sd = tf.maximum(tf.nn.softplus(qU.distribution.scale), 1e-5).eval()
	#np.savetxt(output_file, qU_sd, fmt="%s",delimiter='\t')

	output_file = output_root + '_V.txt'
	np.savetxt(output_file, qV.mean().eval(), fmt="%s",delimiter='\t')

	output_file = output_root + '_re_mean.txt'
	np.savetxt(output_file, q_ind_s.mean().eval(), fmt="%s",delimiter='\t')

	output_file = output_root + '_mu.txt'
	np.savetxt(output_file, mu.eval(), fmt="%s",delimiter='\n')
	
	output_file = output_root + '_sigma_ind.txt'
	np.savetxt(output_file, sigma_ind.eval(), fmt="%s",delimiter='\n')

	output_file = output_root + '_sigma_resid.txt'
	np.savetxt(output_file, sigma_resid.eval(), fmt="%s",delimiter='\n')


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



	data = dict(N=len(y_test), y=y_test, X=U_scaled, K=U_scaled.shape[1], scale_beta=1.0/(2.0*lasso_param))
	op = UPDATE_V.optimizing(data=data, verbose=False, iter=5000)

	param_vec = np.hstack((op['sigma'], op['alpha'], op['beta']))

	#U_scaled_w_intercept = np.hstack((np.ones((U_scaled.shape[0],1)), U_scaled))
	# Fit model
	#model = sm.MixedLM(y_test, U_scaled_w_intercept, Z)
	#if penalty_vector[-1] != 0.0:
	#	result = model.fit_regularized(method='l1', alpha=penalty_vector)
	#	params = result.fe_params
		#clf = linear_model.Lasso(alpha=penalty_vector[-1],positive=True, fit_intercept=True)
		#clf.fit(U_scaled, y_test)
		#params = np.hstack((clf.intercept_, clf.coef_))
	#else:
	#	pdb.set_trace()
	#	model = sm.OLS(y_test, U_scaled_w_intercept)
	#	result = model.fit()
	#	params = result.params
	#return result.fe_params
	return param_vec

def update_factor_matrix_with_genotype_intercept_one_test(test_number, Y, G, U, Z, penalty_vector):
	# Get slice of data corresponding to this test
	y_test = Y[:, test_number]
	g_test = G[:, test_number]

	U2 = np.hstack((np.ones((U.shape[0],1)), U))
	# Get U scaled by genotype for this test
	U_scaled = U2*g_test[:,None]
	U_scaled_w_intercept = np.hstack((np.ones((U_scaled.shape[0],1)), U_scaled))
	# Fit model
	#model = sm.MixedLM(y_test, U_scaled_w_intercept, Z)
	model = sm.OLS(y_test, U_scaled_w_intercept)
	if penalty_vector[-1] != 0.0:
		result = model.fit_regularized(method='elastic_net', L1_wt=1.0, alpha=penalty_vector)
	else:
		result = model.fit()
	#return result.fe_params
	return result.params

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
		V_full = np.transpose(np.asmatrix(V_arr))
		sigma = np.squeeze(np.asarray(V_full[0,:]))
		intercept = np.squeeze(np.asarray(V_full[1,:]))
		#intercept_mat = (np.ones((num_samples,1))*np.asmatrix(intercept))
		V = np.asarray(V_full[2:,:])
	'''
	elif genotype_intercept == 'True':
		penalty_vector = np.ones(num_latent_factor + 2)*lasso_param
		penalty_vector[0] = 0.0
		penalty_vector[1] = 0.0
		V_arr = []
		for test_number in range(num_tests):
			V_arr.append(update_factor_matrix_with_genotype_intercept_one_test(test_number, Y, G, U, Z, penalty_vector))
		# Convert back to matrix
		V_full = np.transpose(np.asmatrix(V_arr))
		intercept = np.squeeze(np.asarray(V_full[0,:]))
		intercept_mat = (np.ones((num_samples,1))*np.asmatrix(intercept))
		V = np.asarray(V_full[2:,:])
		effect_size = np.squeeze(np.asarray(V_full[1,:]))
		intercept_mat = intercept_mat + np.dot(G, np.diag(effect_size))
	elif genotype_intercept == 'True_penalized':
		penalty_vector = np.ones(num_latent_factor + 2)*lasso_param
		penalty_vector[0] = 0.0
		V_arr = []
		for test_number in range(num_tests):
			V_arr.append(update_factor_matrix_with_genotype_intercept_one_test(test_number, Y, G, U, Z, penalty_vector))
		# Convert back to matrix
		V_full = np.transpose(np.asmatrix(V_arr))
		intercept = np.squeeze(np.asarray(V_full[0,:]))
		intercept_mat = (np.ones((num_samples,1))*np.asmatrix(intercept))
		V = np.asarray(V_full[2:,:])
		effect_size = np.squeeze(np.asarray(V_full[1,:]))
		intercept_mat = intercept_mat + np.dot(G, np.diag(effect_size))
	'''
	#else:
	#	print('Assumption error: genotype intercept of this form not implemented')
	#np.savetxt('temp_V.txt', V, fmt="%s", delimiter='\t')
	'''
	# Initialize factor matrix (V)
	#V = np.zeros((num_latent_factor, num_tests))
	#V = np.zeros(num_tests)
	# Run model for each test (independently), ie loop through tests
	for test_number in range(num_tests):
		print(test_number)
		pdb.set_trace()
		V[test_number] = update_factor_matrix_one_test(test_number, Y, G, U, Z)
	'''
	return V, np.asarray(intercept), np.asarray(sigma)

def update_loading_matrix_one_sample(sample_num, Y, G, V, Z, intercept, sigma, lasso_param):
	# Get slice of data corresponding to this sample
	y_test = Y[sample_num, :] - intercept
	g_test = G[sample_num, :]
	# Get V scaled by genotype for this sample
	V_scaled = np.transpose(V)*g_test[:,None]

	pdb.set_trace()

	data = dict(N=len(y_test), y=y_test, X=V_scaled, K=V_scaled.shape[1], scale_beta=1.0/(2.0*lasso_param), sigma=sigma)
	op = UPDATE_U.optimizing(data=data, verbose=False, iter=5000)
	
	#clf = linear_model.Lasso(alpha=lasso_param,positive=True, fit_intercept=False)
	#clf.fit(V_scaled, y_test)
	return np.asarray(op['beta'])

# Update loading matrix (U) with l1 penalized linear model
def update_loading_matrix(Y, G, V, Z, intercept, sigma, num_samples, num_tests, num_latent_factor, lasso_param):
	# Serial version
	U = np.zeros((num_samples, num_latent_factor))
	for sample_num in range(num_samples):
		U[sample_num,:] = update_loading_matrix_one_sample(sample_num, Y, G, V, Z, intercept, sigma, lasso_param)
	return U


#######################
# Part 1: Train eqtl factorization model
########################
def train_eqtl_factorization_model_em_version(sample_overlap_file, expression_file, genotype_file, num_latent_factor, lasso_param_u, lasso_param_v, initialization, output_root):
	# Load in expression data (dimension: num_samplesXnum_tests)
	Y = np.transpose(np.loadtxt(expression_file, delimiter='\t'))
	# Load in genotype data (dimension: num_samplesXnum_tests)
	G = np.transpose(np.loadtxt(genotype_file, delimiter='\t'))
	# Load in sample overlap data
	Z, num_individuals = load_in_sample_overlap_data_v2(sample_overlap_file)
	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]
	# Initialize U with K-means and initialize V to zeros
	V = np.zeros((num_latent_factor, num_tests))
	intercept = np.zeros((num_samples,num_tests))
	if initialization == 'kmeans':
		U = initialize_sample_loading_matrix_with_kmeans(num_samples, num_latent_factor, Y)
	elif initialization == 'random':
		U = np.random.random(size=(num_samples, num_latent_factor))
	elif initialization == 'fixed':
		U = np.zeros((num_samples, num_latent_factor)) + .1
	else:
		print('ASSUMPTION ERROR: Initialization not implemented')
	print('Initialization complete.')
	#####################################
	# Start iterative optimization process
	######################################
	num_iter = 70
	for itera in range(num_iter):
		print('Iteration ' + str(itera))
		# Update factor matrix (V) with linear mixed model
		old_V = V
		V, intercept, sigma = update_factor_matrix(Y, G, U, Z, num_samples, num_tests, num_latent_factor, lasso_param_v)		
		#V2 = np.loadtxt('temp_V.txt')
		# Update loading matrix (U) with l1 penalized linear model
		old_U = U
		U = update_loading_matrix(Y, G, V, Z, intercept, sigma, num_samples, num_tests, num_latent_factor, lasso_param_u)
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
output_root = output_dir + file_stem + '_edward_model_lasso_U_' + str(lasso_param_u) + '_lasso_V_' + str(lasso_param_v) + '_initialization_' + initialization
print(output_root)
train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, lasso_param_u, lasso_param_v, output_root)

