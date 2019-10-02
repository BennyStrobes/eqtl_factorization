from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from edward.models import Gamma, Poisson, Normal, PointMass, \
    TransformedDistribution
import edward as ed
import tensorflow as tf
import numpy as np 
import os
import sys
import pdb




# Load in sample overlap data
def load_in_sample_overlap_data(sample_overlap_file):
	f = open(sample_overlap_file)
	Z = []
	for line in f:
		line = line.rstrip()
		Z.append([float(line)])
	f.close()
	num_individuals = float(line) + 1
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


#######################
# Part 1: Train eqtl factorization model
########################
def train_eqtl_factorization_model(sample_overlap_file, expression_file, genotype_file, num_latent_factors, sparsity_parameter, output_root):
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
	U = Gamma(concentration=1.0, rate=sparsity_parameter, sample_shape=[num_samples, num_latent_factors])
	V = tf.get_variable("V", [num_latent_factors, num_tests])
	
	yhat = (tf.multiply(genotype, tf.matmul(U, V)) + tf.gather_nd(eta_ind, ind_ph) + tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(mu)))
	y = Normal(loc=yhat, scale=tf.matmul(tf.ones([num_samples, num_tests]), tf.matrix_diag(sigma_resid)))

	###############################
	## Inference set up
	###############################
	q_ind_s = Normal(
		loc=tf.get_variable("q_ind_s/loc", [num_individuals, num_tests]),
		scale=tf.nn.softplus(tf.get_variable("q_ind_s/scale", [num_individuals, num_tests])))
	qU = lognormal_q(U.shape)

	latent_vars = {U: qU, eta_ind: q_ind_s}

	data = {y: Y, ind_ph: Z, genotype: G}

	inference = ed.KLqp(latent_vars, data)
	tf.global_variables_initializer().run()
	inference.run(n_iter=1000)
	
	#########################
	# Save data to output
	#########################
	output_file = output_root + '_qU_mean.txt'
	qU_mean = np.exp(qU.distribution.loc.eval())
	np.savetxt(output_file, qU_mean, fmt="%s",delimiter='\t')

	output_file = output_root + '_qU_sd.txt'
	qU_sd = tf.maximum(tf.nn.softplus(qU.distribution.scale), 1e-5).eval()
	np.savetxt(output_file, qU_sd, fmt="%s",delimiter='\t')

	output_file = output_root + '_V.txt'
	np.savetxt(output_file, V.eval(), fmt="%s",delimiter='\t')

	output_file = output_root + '_re_mean.txt'
	np.savetxt(output_file, q_ind_s.mean().eval(), fmt="%s",delimiter='\t')

	output_file = output_root + '_mu.txt'
	np.savetxt(output_file, mu.eval(), fmt="%s",delimiter='\n')
	
	output_file = output_root + '_sigma_ind.txt'
	np.savetxt(output_file, sigma_ind.eval(), fmt="%s",delimiter='\n')

	output_file = output_root + '_sigma_resid.txt'
	np.savetxt(output_file, sigma_resid.eval(), fmt="%s",delimiter='\n')











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


#######################
# Part 1: Train eqtl factorization model
########################
sparsity_parameter = 1.0
output_root = output_dir + file_stem + str(sparsity_parameter) + '_sparsity_parameter'
train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, sparsity_parameter, output_root)

sparsity_parameter = 10
output_root = output_dir + file_stem + str(sparsity_parameter) + '_sparsity_parameter'
train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, sparsity_parameter, output_root)

sparsity_parameter = 100
output_root = output_dir + file_stem + str(sparsity_parameter) + '_sparsity_parameter'
train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, num_latent_factors, sparsity_parameter, output_root)