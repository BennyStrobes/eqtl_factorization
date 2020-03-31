import numpy as np 
import os
import sys
import pdb
import h5py
from sklearn.decomposition import PCA


# Print covariates to readble file
def generate_cell_covariate_file(data, cell_covariate_file):
	# Open output file handle
	t = open(cell_covariate_file, 'w')
	# print header to output file
	t.write('cell_barcode\tdisease_cov\tct_cov\tpop_cov\tind_cov\twell\tbatch_cov\tbatch\tpercent_mito\tn_counts\tSLEDAI\tBroad\tFemale\tn_genes\tPF4\tSDPR\tGNG11\tPPBP\tPC2\tPC3\tPC10\tPC12\tPC14\tPC27\tPC28\tlouvain\tleiden\tsite\tdisease_pop_site_cov\tumap_density_disease_pop_site_cov\n')
	# Iterate over cells
	for sample_num in range(len(data['obs'])):
		# Extract covariate vec for this cell
		covariate_vec = data['obs'][sample_num]
		# Error check to make sure covariate_vec has correct number of columns
		if len(covariate_vec) != 30:
			print('Error: quitting')
			sys.exit()
		# Print columns (1 by one..)
		processed_covariate_vec = []
		# Cell barcode
		processed_covariate_vec.append(covariate_vec[0])
		# disease_cov (categorical)
		processed_covariate_vec.append(data['uns']['disease_cov_categories'][covariate_vec[1]])
		# ct_cov (categorical)
		processed_covariate_vec.append(data['uns']['ct_cov_categories'][covariate_vec[2]])
		# pop_cov (categorical)
		processed_covariate_vec.append(data['uns']['pop_cov_categories'][covariate_vec[3]])
		# ind_cov (categorical)
		processed_covariate_vec.append(data['uns']['ind_cov_categories'][covariate_vec[4]])
		# well (categorical)
		processed_covariate_vec.append(data['uns']['well_categories'][covariate_vec[5]])
		# batch_cov (categorical)
		processed_covariate_vec.append(data['uns']['batch_cov_categories'][covariate_vec[6]])
		# batch (categorical)
		processed_covariate_vec.append(data['uns']['batch_cov_categories'][covariate_vec[7]])
		# percent_mito
		processed_covariate_vec.append(str(covariate_vec[8]))
		# n_counts
		processed_covariate_vec.append(str(covariate_vec[9]))
		# SLEDAI (categorical)
		processed_covariate_vec.append(data['uns']['SLEDAI_categories'][covariate_vec[10]])
		# Broad
		processed_covariate_vec.append(str(covariate_vec[11]))
		# Female
		processed_covariate_vec.append(str(covariate_vec[12]))
		# n_genes
		processed_covariate_vec.append(str(covariate_vec[13]))
		# several specific genes
		processed_covariate_vec.append(str(covariate_vec[14]))
		processed_covariate_vec.append(str(covariate_vec[15]))
		processed_covariate_vec.append(str(covariate_vec[16]))
		processed_covariate_vec.append(str(covariate_vec[17]))
		# several PCs
		processed_covariate_vec.append(str(covariate_vec[18]))
		processed_covariate_vec.append(str(covariate_vec[19]))
		processed_covariate_vec.append(str(covariate_vec[20]))
		processed_covariate_vec.append(str(covariate_vec[21]))
		processed_covariate_vec.append(str(covariate_vec[22]))
		processed_covariate_vec.append(str(covariate_vec[23]))
		processed_covariate_vec.append(str(covariate_vec[24]))
		# louvain (categorical)
		processed_covariate_vec.append(data['uns']['louvain_categories'][covariate_vec[25]])
		# leiden (categorical)
		processed_covariate_vec.append(data['uns']['leiden_categories'][covariate_vec[26]])
		# site (categorical)
		processed_covariate_vec.append(data['uns']['site_categories'][covariate_vec[27]])
		# disease_pop_site_cov
		processed_covariate_vec.append(data['uns']['disease_pop_site_cov_categories'][covariate_vec[28]])
		# umap_density_disease_pop_site_cov
		processed_covariate_vec.append(str(covariate_vec[29]))
		# Error check to make sure processed_covariate_vec has correct number of columns
		if len(processed_covariate_vec) != 30:
			print('Error: quitting')
			sys.exit()
		# Write new line to output file
		t.write('\t'.join(processed_covariate_vec) + '\n')
	# close file handle
	t.close()


# Print expression data to a readable file
def generate_expression_file(data, sc_expression_file):
	np.savetxt(sc_expression_file, data['X'][:,:], fmt="%s", delimiter='\t')

# Print gene_names to a readable file
def generate_gene_names_file(data, gene_names_file):
	# Number of genes in file
	num_genes = len(data['var'])
	# Open file handle
	t = open(gene_names_file, 'w')
	# Loop through genes
	for gene_num in range(num_genes):
		# error checking
		if len(data['var'][gene_num]) != 1:
			print('assumption error')
			pdb.set_trace()
		# Get name of gene and print it to output file
		gene_name = data['var'][gene_num][0]
		t.write(gene_name + '\n')
	t.close()

# Get cell type colors and print to file
def print_cell_type_colors_to_file(data, cell_type_colors_file):
	# Extract fields
	cell_type_names = data['uns']['ct_cov_categories'][:]
	cell_type_colors = data['uns']['ct_cov_colors'][:]
	# Open output file handle
	t = open(cell_type_colors_file, 'w')
	# print header to output file
	t.write('cell_type\thex_color\n')
	# loop through cell types
	for index, cell_type in enumerate(cell_type_names):
		t.write(cell_type + '\t' + cell_type_colors[index] + '\n')
	t.close()

# Grab UMAP scores from ye-lab data structure
def generate_umap_from_ye_lab_file(data, umap_file):
	# Open output file handle
	t = open(umap_file, 'w')
	# Get number of samples (ie cells)
	num_samples = len(data['obsm'])
	# loop through cells
	for sample_num in range(num_samples):
		# get umap scores for this sample
		umap_scores = data['obsm'][sample_num][2].astype(str)
		# write to output file
		t.write('\t'.join(umap_scores) + '\n')
	t.close()

# Grab UMAP scores from ye-lab data structure
def generate_pca_from_ye_lab_file(data, pca_file):
	# Open output file handle
	t = open(pca_file, 'w')
	# Get number of samples (ie cells)
	num_samples = len(data['obsm'])
	# loop through cells
	for sample_num in range(num_samples):
		# get umap scores for this sample
		pca_scores = data['obsm'][sample_num][0].astype(str)
		# write to output file
		t.write('\t'.join(pca_scores) + '\n')
	t.close()

def pca_understanding(X):
	num_genes = X.shape[1]
	for gene_num in range(num_genes):
	#	X[:, gene_num] = (X[:, gene_num] - np.mean(X[:, gene_num]))/np.std(X[:, gene_num])
		X[:, gene_num] = (X[:, gene_num] - np.mean(X[:, gene_num]))
	print(X.shape)
	# SVD approach
	uuu, sss, vh = np.linalg.svd(np.transpose(X), full_matrices=False)
	svd_loadings = np.transpose(vh)[:,:40]

	pdb.set_trace()

	ye = np.loadtxt('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/processed_expression/pca_scores_ye_lab.txt')


	for pc_num in range(39):
		corry = np.corrcoef(svd_loadings[:,pc_num], ye[:, pc_num])[0,1]
		print(str(pc_num) + ' : ' + str(corry))
	pdb.set_trace()

# Filter expression and covariate file to only contain SLE individuals
def filter_data_to_only_sle_individuals(cell_covariate_file, sc_expression_file, filtered_cell_covariate_file, filtered_sc_expression_file):
	# Initialize vector to keep track of which cells are from sle individuals (will be 1 if they are from, 0 if they are not)
	valid_cells = []
	# Filter covariate file first
	f = open(cell_covariate_file)
	t = open(filtered_cell_covariate_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		if data[1] != 'sle':
			valid_cells.append(0)
			continue
		valid_cells.append(1)
		t.write(line + '\n')
	f.close()
	t.close()
	# Filter expression file now
	f = open(sc_expression_file)
	t = open(filtered_sc_expression_file, 'w')
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Check if cell corresponds to sle individual
		if valid_cells[counter] == 1:
			t.write(line + '\n')
		counter = counter + 1
	f.close()
	t.close()

# Filter expression and covariate file to only contain SLE individuals
def filter_data_to_only_sle_individuals_and_randomly_subset(cell_covariate_file, sc_expression_file, filtered_cell_covariate_file, filtered_sc_expression_file):
	np.random.seed(0)
	# Initialize vector to keep track of which cells are from sle individuals (will be 1 if they are from, 0 if they are not)
	valid_cells = []
	# Filter covariate file first
	f = open(cell_covariate_file)
	t = open(filtered_cell_covariate_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		if data[1] != 'sle' or np.random.binomial(1,.1) == 0:
			valid_cells.append(0)
			continue
		valid_cells.append(1)
		t.write(line + '\n')
	f.close()
	t.close()
	# Filter expression file now
	f = open(sc_expression_file)
	t = open(filtered_sc_expression_file, 'w')
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Check if cell corresponds to sle individual
		if valid_cells[counter] == 1:
			t.write(line + '\n')
		counter = counter + 1
	f.close()
	t.close()


def standardize_gene_expression_data(input_expression_file, output_standardized_expression_file):
	X = np.loadtxt(input_expression_file)
	num_genes = X.shape[1]
	print(X.shape)
	# Loop through genes
	for gene_num in range(num_genes):
		X[:, gene_num] = (X[:, gene_num] - np.mean(X[:, gene_num]))/np.std(X[:, gene_num])
	# Save to output file
	np.savetxt(output_standardized_expression_file, X, fmt="%s", delimiter='\t')

# Generate expression PC loadings and variance explained of those expression PCs
def generate_pca_scores_and_variance_explained(filtered_standardized_sc_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file):
	# Load in data
	X = np.loadtxt(filtered_standardized_sc_expression_file)

	# Run PCA (via SVD)
	uuu, sss, vh = np.linalg.svd(np.transpose(X), full_matrices=False)
	svd_loadings = np.transpose(vh)[:,:num_pcs]

	# Save to output file
	np.savetxt(filtered_cells_pca_file, svd_loadings, fmt="%s", delimiter='\t')

	# Compute variance explained
	ve = (np.square(sss)/np.sum(np.square(sss)))[:num_pcs]
	np.savetxt(filtered_cells_pca_ve_file, ve, fmt="%s", delimiter='\n')

# Input data from Ye-Lab
input_h5py_file = sys.argv[1]
# Output directory
processed_expression_dir = sys.argv[2]
# gene annotation file
gene_annotation_file = sys.argv[3]



# Load in data
data = h5py.File(input_h5py_file, 'r')


# Print covariates to readable file
cell_covariate_file = processed_expression_dir + 'cell_covariates_ye_lab.txt'
#generate_cell_covariate_file(data, cell_covariate_file)

# Print expression data to a readable file
sc_expression_file = processed_expression_dir + 'single_cell_expression_ye_lab.txt'
#generate_expression_file(data, sc_expression_file)

# Print gene_names to a readable file
gene_names_file = processed_expression_dir + 'gene_names_ye_lab.txt'
#generate_gene_names_file(data, gene_names_file)

# Get cell type colors and print to file
cell_type_colors_file = processed_expression_dir + 'cell_type_colors_ye_lab.txt'
#print_cell_type_colors_to_file(data, cell_type_colors_file)

# Grab UMAP scores from ye-lab data structure
umap_file = processed_expression_dir + 'umap_scores_ye_lab.txt'
#generate_umap_from_ye_lab_file(data, umap_file)

# Grab UMAP scores from ye-lab data structure
# This data was genead by centering (mean 0) the expression data (sc_expression_file), but not scaling (each gene does not have SD 1)
# Also worth noting that they multiplied the loadings by the singular values
pca_file = processed_expression_dir + 'pca_scores_ye_lab.txt'
#generate_pca_from_ye_lab_file(data, pca_file)

# Filter expression and covariate file to only contain SLE individuals
filtered_cell_covariate_file = processed_expression_dir + 'cell_covariates_sle_individuals.txt'
filtered_sc_expression_file = processed_expression_dir + 'single_cell_expression_sle_individuals.txt'
#filter_data_to_only_sle_individuals(cell_covariate_file, sc_expression_file, filtered_cell_covariate_file, filtered_sc_expression_file)

# Standardize expression so each gene has mean 0 and standard deviation 1 (across individuals)
filtered_standardized_sc_expression_file = processed_expression_dir + 'single_cell_expression_sle_individuals_standardized.txt'
#standardize_gene_expression_data(filtered_sc_expression_file, filtered_standardized_sc_expression_file)

# Generate expression PC loadings and variance explained of those expression PCs
num_pcs=200
filtered_cells_pca_file = processed_expression_dir + 'pca_scores_sle_individuals.txt'
filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_explained_sle_individuals.txt'
#generate_pca_scores_and_variance_explained(filtered_standardized_sc_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file)



# Filter expression and covariate file to only contain SLE individuals
filtered_cell_covariate_file = processed_expression_dir + 'cell_covariates_sle_individuals_random_subset.txt'
filtered_sc_expression_file = processed_expression_dir + 'single_cell_expression_sle_individuals_random_subset.txt'
filter_data_to_only_sle_individuals_and_randomly_subset(cell_covariate_file, sc_expression_file, filtered_cell_covariate_file, filtered_sc_expression_file)
print('1')
# Standardize expression so each gene has mean 0 and standard deviation 1 (across individuals)
filtered_standardized_sc_expression_file = processed_expression_dir + 'single_cell_expression_sle_individuals_random_subset_standardized.txt'
standardize_gene_expression_data(filtered_sc_expression_file, filtered_standardized_sc_expression_file)
print('2')
# Generate expression PC loadings and variance explained of those expression PCs
num_pcs=200
filtered_cells_pca_file = processed_expression_dir + 'pca_scores_sle_individuals_random_subset.txt'
filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_explained_sle_individuals_random_subset.txt'
generate_pca_scores_and_variance_explained(filtered_standardized_sc_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file)
print('3')



