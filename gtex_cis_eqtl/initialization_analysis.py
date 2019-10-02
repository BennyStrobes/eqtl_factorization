import numpy as np 
import os
import sys
import pdb
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture


def initialize_sample_loading_matrix_with_kmeans(num_samples, num_latent_factor, Y):
	# Sample loading matrix
	U = np.zeros((num_samples, num_latent_factor))
	# Fit kmeans model
	kmeans = KMeans(n_clusters=num_latent_factor).fit(Y)
	# Fill in sample loading matrix with kmeans assignments
	for sample_index, kmeans_assignment in enumerate(kmeans.labels_):
		U[sample_index, kmeans_assignment] = 1.0
	return U

def initialize_sample_loading_matrix_with_gmm(num_samples, num_latent_factor, Y):
	# Sample loading matrix
	U = np.zeros((num_samples, num_latent_factor))
	# Fit kmeans model
	gmm = GaussianMixture(n_components=num_latent_factor, covariance_type='full')
	gmm_fit = gmm.fit(Y)
	return gmm_fit.predict_proba(Y)

def standardize_Y(Y):
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]
	standardize_Y = np.zeros((num_samples, num_tests))

	for test_num in range(num_tests):
		standardize_Y[:,test_num] = (Y[:, test_num] - np.mean(Y[:, test_num]))/np.std(Y[:, test_num])

	return standardize_Y

def get_tissue_names(file_name): 
	f = open(file_name)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		arr.append(line)
	return arr

def compute_kmeans_to_tissue(tissue_names, kmeans_matrix, num_latent_factors):
	dicti = {}
	num_samples = len(tissue_names)
	for sample_num in range(num_samples):
		tissue_type = tissue_names[sample_num]
		counts = kmeans_matrix[sample_num,:]
		if tissue_type not in dicti:
			dicti[tissue_type] = np.zeros(num_latent_factors)
		dicti[tissue_type] = dicti[tissue_type] + counts
	new_dicti = {}
	for tissue_type in sorted(dicti.keys()):
		new_dicti[tissue_type] = dicti[tissue_type]/np.sum(dicti[tissue_type])
		print(tissue_type)
		print(new_dicti[tissue_type])
	return new_dicti

expression_training_file = sys.argv[1]
genotype_training_file = sys.argv[2]
num_latent_factors = int(sys.argv[3])
eqtl_results_dir = sys.argv[4]
tissue_names_file = sys.argv[5]

tissue_names = get_tissue_names(tissue_names_file)

# Load in expression data (dimension: num_samplesXnum_tests)
Y = np.transpose(np.loadtxt(expression_training_file, delimiter='\t'))
# Load in genotype data (dimension: num_samplesXnum_tests)
G = np.transpose(np.loadtxt(genotype_training_file, delimiter='\t'))




num_samples = Y.shape[0]

# GMM on expression
output_file = eqtl_results_dir + 'initialization_of_gmm_on_expression.txt'
#gmm_matrix = initialize_sample_loading_matrix_with_gmm(num_samples, num_latent_factors,(Y/(G+1.0)))
#dicti = compute_kmeans_to_tissue(tissue_names, gmm_matrix, num_latent_factors)

#np.savetxt(output_file, gmm_matrix, fmt="%s",delimiter='\t')

# Kmeans on expression
output_file = eqtl_results_dir + 'initialization_of_kmeans_on_expression.txt'
pdb.set_trace()
kmeans_matrix = initialize_sample_loading_matrix_with_kmeans(num_samples, num_latent_factors, (Y/(G+2.0)))
dicti = compute_kmeans_to_tissue(tissue_names, kmeans_matrix, num_latent_factors)



np.savetxt(output_file, kmeans_matrix, fmt="%s",delimiter='\t')


