import numpy as np
import os
import sys
import pdb
import random



# Simulate a haplotype vector (across individuals) for one SNP
def simulate_haplotype(num_individuals, num_minor_alleles):
	haplotype = np.zeros(num_individuals)
	haplotype[random.sample(range(num_individuals), num_minor_alleles)] = np.ones(num_minor_alleles)
	return haplotype

# Simulate a genotype vector (across individuals) for one SNP
def simulate_genotype(maf, num_individuals):
	num_minor_alleles = int(round(maf*num_individuals))
	# Simulate each haplotype seperately
	h1 = simulate_haplotype(num_individuals, num_minor_alleles)
	h2 = simulate_haplotype(num_individuals, num_minor_alleles)
	# Genotype is then simply the sum of the two haplotypes
	genotype = h1 + h2
	return genotype

def simulate_genotype_matrix(num_individuals, num_cells_per_individual, num_tests, maf):
	# Number of samples is just the product of the number of individuals and the number of cells
	num_samples = num_individuals*num_cells_per_individual
	#Initialize genotype matrix
	genotype_matrix = np.zeros((num_samples, num_tests))
	# Loop through tests (simulate genotype for each test independently)
	for test_number in range(num_tests):
		# Simulate genotype vector for this individual
		genotype_vector = simulate_genotype(maf, num_individuals)
		# Genotype matrix is many repeats of this vector (because we have multiple samples from the same individual)
		genotype_vector_across_samples = np.repeat(genotype_vector,num_cells_per_individual)
		# Now fill in matrix for this column
		genotype_matrix[:, test_number] = genotype_vector_across_samples
	return genotype_matrix

def generate_sample_loadings(num_samples, num_latent_factors):
	# Initialize matrix 
	loadings = np.zeros((num_samples, num_latent_factors))
	# Iterate through all elements of the matrix
	for sample_num in range(num_samples):
		for latent_factor_num in range(num_latent_factors):
			# Generate 2 numbers between 0 and 1
			num1 = random.uniform(0, 1)
			num2 = random.uniform(0, 1)
			# With 50% chance make it a non-zero number between 0 and 1. Otherwise make it zero
			if num1 > .5:
				loadings[sample_num, latent_factor_num] = num2
	return loadings

def generate_factor_matrix(num_tests, num_latent_factors):
	# Initialize matrix
	factors = np.zeros((num_latent_factors, num_tests))
	# Iterate through all elements of the matrix
	for test_num in range(num_tests):
		for latent_factor_num in range(num_latent_factors):
			# Generate a random number between 0 and 1
			num1 = random.uniform(0,1)
			# 80% of time a particular (test, factor) has no effect
			if num1 < .2:
				factors[latent_factor_num, test_num] = np.random.normal()
	return factors

def simulate_expression_matrix(num_individuals, num_cells_per_individual, num_tests, num_latent_factors, genotype_matrix):
	num_samples = num_individuals*num_cells_per_individual
	# Generate loading matrix
	loading_matrix = generate_sample_loadings(num_samples, num_latent_factors)
	# Generate factor matrix
	factor_matrix = generate_factor_matrix(num_tests, num_latent_factors)

	# Initialize expression matrix
	Y = np.zeros((num_samples, num_tests))

	# Add eqtl effects to expression matrix
	R_true = np.dot(loading_matrix, factor_matrix)
	eqtl_effects = np.multiply(R_true, genotype_matrix)
	Y = Y + eqtl_effects

	# Initialize simulated matrix containing random effects for each individual
  	random_effects = np.zeros((num_individuals, num_tests))
  	Z = []
  	# Fill in matrix keeping track of which individual each sample belongs to
  	sample_num = 0
  	for indi_num in range(num_individuals):
  		for cell_num in range(num_cells_per_individual):
  			Z.append(indi_num)
  			sample_num = sample_num + 1

  	# Simulate random effect stdevs for each gene
  	gene_random_effects_sdevs =  np.sqrt(np.random.exponential(size=num_tests))
  	# Simulate resdiaul effect sdevs for each gene
  	gene_residual_sdevs = np.sqrt(np.random.exponential(size=num_tests))
  	# Simulate mean for each gene
  	gene_mean = np.random.normal(size=num_tests)
  	# Simulate gene expression data for each gene
  	for gene_num in range(num_tests):
  		# Simulate random effects for each individual
  		gene_random_effect = np.random.normal(loc=0, scale=gene_random_effects_sdevs[gene_num], size=num_individuals)
  		random_effects[:, gene_num] = gene_random_effect
  		# Loop through samples
  		for sample_num in range(num_samples):
  			# Get index of individual corresponding to this sample
  			indi_id = Z[sample_num]
  			# Get simulated predicted mean for this sample
  			predicted_mean = gene_mean[gene_num] + gene_random_effect[int(indi_id)]
  			# Randomly draw expression sample for this sample
  			Y[sample_num, gene_num] = Y[sample_num, gene_num] + np.random.normal(loc=predicted_mean, scale=gene_residual_sdevs[gene_num])
  	return Y, random_effects, gene_mean, gene_residual_sdevs, gene_random_effects_sdevs, loading_matrix, factor_matrix, Z


###################
# Command line args
####################
num_tests = int(sys.argv[1])
maf = float(sys.argv[2])
num_individuals = int(sys.argv[3])
num_cells_per_individual = int(sys.argv[4])
num_latent_factors = int(sys.argv[5])
simulated_data_dir = sys.argv[6]
file_stem = sys.argv[7]

##########################
# Simulate genotype matrix
##########################
genotype_matrix = simulate_genotype_matrix(num_individuals, num_cells_per_individual, num_tests, maf)

##########################
# Simulate gene expression
##########################
expression_matrix, random_effects, gene_mean, gene_residual_sdevs, gene_random_effects_sdevs, loading_matrix, factor_matrix, Z = simulate_expression_matrix(num_individuals, num_cells_per_individual, num_tests, num_latent_factors, genotype_matrix)

##########################
# Save matrices to output
##########################
# Genotype matrix
output_file = simulated_data_dir + file_stem + 'simulated_genotype.txt'
np.savetxt(output_file, np.transpose(genotype_matrix), fmt="%s", delimiter='\t')

# Expression matrix
output_file = simulated_data_dir + file_stem + 'simulated_expression.txt'
np.savetxt(output_file, np.transpose(expression_matrix), fmt="%s",delimiter='\t')

# Random effects
output_file = simulated_data_dir + file_stem + 'simulated_random_effects.txt'
np.savetxt(output_file, random_effects, fmt="%s",delimiter='\t')

# Gene mean
output_file = simulated_data_dir + file_stem + 'simulated_gene_mean.txt'
np.savetxt(output_file, gene_mean, fmt="%s",delimiter='\n')

# gene_residual_sdevs
output_file = simulated_data_dir + file_stem + 'simulated_gene_residual_sdevs.txt'
np.savetxt(output_file, gene_residual_sdevs, fmt="%s",delimiter='\n')

# gene random effects sdevs
output_file = simulated_data_dir + file_stem + 'simulated_gene_random_effects_sdevs.txt'
np.savetxt(output_file, gene_random_effects_sdevs, fmt="%s",delimiter='\n')

# loading matrix
output_file = simulated_data_dir + file_stem + 'simulated_loading_matrix.txt'
np.savetxt(output_file, loading_matrix, fmt="%s",delimiter='\t')

# factor matrix
output_file = simulated_data_dir + file_stem + 'simulated_factor_matrix.txt'
np.savetxt(output_file, factor_matrix, fmt="%s",delimiter='\t')

# individual id array
output_file = simulated_data_dir + file_stem + 'simulated_individual_id.txt'
np.savetxt(output_file, Z, fmt="%s", delimiter='\n')
