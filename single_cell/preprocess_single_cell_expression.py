import numpy as np 
import os
import sys
import pdb
import h5py
from sklearn.decomposition import PCA
import scanpy as sc
from anndata import AnnData




def extract_protein_coding_known_autosomal_genes(gene_struct, gene_annotation_file):
	'''
	# Simple error checking to understand gene struct
	col_names = gene_struct.columns
	ensamble_ids = gene_struct[col_names[0]]
	for col_name in col_names:
		if len(np.unique(gene_struct[col_name])) != 32738:
			print('assumption errror')
			pdb.set_trace()
		if np.array_equal(gene_struct[col_name], ensamble_ids) == False:
			print('assumption errorr')
			pdb.set_trace()
	'''
	# Make dictionary list of genes that are protein-coding, known, and autosomal
	valid_genes = {}
	f = open(gene_annotation_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Simple error check
		if len(data) != 8 and len(data) != 9:
			print('assumption errror in processing gene annotation file')
			pdb.set_trace()
		gene_id = data[5]
		# Some more error checks
		start = int(data[2])
		end = int(data[3])
		if start > end:
			print('assumption errror in processing gene annotation file')
			pdb.set_trace()
		if data[6] != 'protein_coding' or data[7] != 'KNOWN':
			continue
		if data[1] == 'chrX' or data[1] == 'chrY' or data[1] == 'chrM':
			continue
		valid_genes[data[5] + '_' + data[0].split('.')[0]] = 1
	f.close()
	# Make binary vectory corresponding to ordered list of whehter our genes are protein-coding, autosomal and known
	binary_vector = []
	ensamble_ids = gene_struct[gene_struct.columns[0]]
	gene_ids = gene_struct.index
	num_genes = len(gene_ids)
	for gene_num in range(num_genes):
		ensamble_id = ensamble_ids[gene_num]
		gene_id = gene_ids[gene_num]
		if gene_id + '_' + ensamble_id not in valid_genes:
			binary_vector.append(False)
		else:
			binary_vector.append(True)
	return np.asarray(binary_vector)




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

# Generate pseudobulk covariate file (along with ordered list of pseudobulk samples)
def generate_pseudobulk_covariate_file(adata, pseudobulk_covariate_file):
	# Loop through cell covariate file to create  list of individual-cell_type pairs (these are our samples)
	samples = {}
	num_cells = adata.obs.shape[0]
	sample_counts = {}
	for cell_num in range(num_cells):
		# Extract relevent fields
		cell_type = adata.obs.ct_cov[cell_num]
		race = adata.obs.pop_cov[cell_num]
		ind_id = adata.obs.ind_cov[cell_num]
		# Come up with sample id name
		cell_type = '_'.join(cell_type.split(' '))
		sample_name = ind_id + ':' + cell_type
		if sample_name not in samples:  # Add sample to list
			samples[sample_name] = race
		else:  # Simple error checking making sure ever individual was mapped to the same race
			if samples[sample_name] != race:
				print('assumption errro')
		if sample_name not in sample_counts:
			sample_counts[sample_name] = 0
		sample_counts[sample_name] = sample_counts[sample_name] + 1
	# Print to output file
	t = open(pseudobulk_covariate_file, 'w')
	t.write('sample_name\tct_cov\tpop_cov\tind_cov\tnum_cells\tct_cov_readable\n')
	# Loop through samples and print to output file
	for sample_name in samples.keys():
		# Extract relevent fields
		ind_cov = sample_name.split(':')[0]
		ct_cov = ' '.join(sample_name.split(':')[1].split('_'))
		pop_cov = samples[sample_name]
		# Print to output file
		t.write(sample_name + '\t' + ct_cov + '\t' + pop_cov + '\t' + ind_cov + '\t' + str(sample_counts[sample_name]) + '\t' + sample_name.split(':')[1] + '\n')
	t.close()

# Generate summed pseudobulk counts
def generate_pseudobulk_counts(adata, raw_pseudobulk_expression_file, pseudobulk_covariate_file):
	# Goal: convert adata.raw.X (cellsXgenes) to pseudobulk_matrix (pseudobulk_samplesXgenes)
	# First create mapping from pseudobulk sample name to index. Also count number of pseudobulk samples
	pseudobulk_sample_to_index = {}
	f = open(pseudobulk_covariate_file)
	head_count = 0
	num_pseudobulk_samples = 0

	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_id = data[0]
		pseudobulk_sample_to_index[sample_id] = num_pseudobulk_samples
		num_pseudobulk_samples = num_pseudobulk_samples + 1
	f.close()
	# Now create array of length number of cells where each element is mapping from cell to pseudobulk sample
	cell_to_pseudobulk_sample_index = []
	num_cells = adata.obs.shape[0]
	for cell_num in range(num_cells):
		# Extract relevent fields
		cell_type = adata.obs.ct_cov[cell_num]
		ind_id = adata.obs.ind_cov[cell_num]
		# Come up with sample id name for this cll
		cell_type = '_'.join(cell_type.split(' '))
		sample_name = ind_id + ':' + cell_type
		# Map cells sample id to an index
		pseudobulk_index = pseudobulk_sample_to_index[sample_name]
		cell_to_pseudobulk_sample_index.append(pseudobulk_index)
	cell_to_pseudobulk_sample_index = np.asarray(cell_to_pseudobulk_sample_index)
	# Now generate and fill in pseudobulk counts matrix
	num_genes = adata.raw.X.shape[1]
	pseudobulk_matrix = np.zeros((num_pseudobulk_samples, num_genes))
	# convert raw cell counts into np matrix
	raw_cell_counts = adata.raw.X.toarray()
	# Loop through cells and genes and add cell-gene count to pseudobulk count matrix
	for cell_num in range(num_cells):
		pseudobulk_sample_index = cell_to_pseudobulk_sample_index[cell_num]
		pseudobulk_matrix[pseudobulk_sample_index, :] = pseudobulk_matrix[pseudobulk_sample_index, :] + raw_cell_counts[cell_num, :]
	np.savetxt(raw_pseudobulk_expression_file, pseudobulk_matrix, fmt="%s", delimiter='\t')

# Standardize summed pseudobulk counts (samples X genes)
def standardize_pseudobulk_counts(raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file):
	# Load in raw pseodobulk expression data
	raw_pseudobulk_expression = np.loadtxt(raw_pseudobulk_expression_file)
	num_samples, num_genes = raw_pseudobulk_expression.shape

	# initialize and fill in normalized pseudobulk expression data
	normalized_pseudobulk_expression = np.zeros((num_samples, num_genes))
	for sample_num in range(num_samples):
		normalized_pseudobulk_expression[sample_num, :] = np.log(10000.0*(raw_pseudobulk_expression[sample_num,:]/np.sum(raw_pseudobulk_expression[sample_num,:])) + 1.0)
	
	# initialize and fill in standardized pseudobulk expression data
	standardized_pseudobulk_expression = np.zeros((num_samples, num_genes))
	for gene_num in range(num_genes):
		standardized_pseudobulk_expression[:, gene_num] = (normalized_pseudobulk_expression[:, gene_num] - np.mean(normalized_pseudobulk_expression[:, gene_num]))/np.std(normalized_pseudobulk_expression[:, gene_num])

	np.savetxt(standardized_pseudobulk_expression_file, standardized_pseudobulk_expression, fmt="%s", delimiter='\t')


def generate_pseudobulk_expression_data_wrapper(adata, raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file, pseudobulk_covariate_file):
	# Generate pseudobulk covariate file (along with ordered list of pseudobulk samples)
	generate_pseudobulk_covariate_file(adata, pseudobulk_covariate_file)

	# Generate summed pseudobulk counts
	generate_pseudobulk_counts(adata, raw_pseudobulk_expression_file, pseudobulk_covariate_file)

	# Standardize summed pseudobulk counts (samples X genes)
	standardize_pseudobulk_counts(raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file)

def generate_cell_type_sc_expression_data_wrapper(cell_type, adata, standardized_cell_type_expression_file, standardized_cell_type_covariate_file):
	############################
	# Filter to cells of this cell type
	###############################
	ct_data = adata[adata.obs.ct_cov == cell_type, :]

	############################
	# Save covariates to output file
	###############################
	np.savetxt(standardized_cell_type_covariate_file, ct_data.obs, fmt="%s", delimiter='\t', header='\t'.join(ct_data.obs.columns), comments='')

	######################	
	# Normalize data
	#######################
	ct_data_clean = AnnData(ct_data.raw.X.toarray())
	sc.pp.normalize_total(ct_data_clean, target_sum=1e4)
	sc.pp.log1p(ct_data_clean)
	sc.pp.scale(ct_data_clean, max_value=10)
	############################
	# Save expression to output
	###########################
	np.savetxt(standardized_cell_type_expression_file, ct_data_clean.X, fmt="%s", delimiter='\t')

#####################
# Command line args
######################
input_h5py_file = sys.argv[1]
processed_expression_dir = sys.argv[2]
gene_annotation_file = sys.argv[3]



######################
# Filtering parameters
#######################
min_genes = 400
# Min fraction of expressed cells for a gene
min_fraction_of_cells = .01
# Random subset
random_subset = False
np.random.seed(0)

'''
######################
# Load in ScanPy data
#######################
adata = sc.read_h5ad(input_h5py_file)


######################
# Filter data
#######################
# Limit to SLE individuals
adata = adata[adata.obs.disease_cov == "sle", :]

if random_subset == True:
	num_cells = adata.X.shape[0]
	adata.obs['random_subset'] = np.random.uniform(size=num_cells) < (1/8)
	adata = adata[adata.obs.random_subset == True, :]


# Standard filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Extract percent mito-counts
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

# Standard filtering
adata = adata[adata.obs.n_genes < 2500, :]
adata = adata[adata.obs.percent_mito < 0.05, :]

# Limit to protein-coding, known, autosomal genes
gene_indices = extract_protein_coding_known_autosomal_genes(adata.var, gene_annotation_file)
adata.var['protein_coding_known_autosomal'] = gene_indices
adata = adata[:, adata.var.protein_coding_known_autosomal==True]

# Filter rows and columns
print(adata.X.shape)
sc.pp.filter_cells(adata, min_genes=min_genes)
print(adata.X.shape)
sc.pp.filter_genes(adata, min_cells=(adata.X.shape[0])*min_fraction_of_cells)
print(adata.X.shape)
#######################
# Save un-normalized expression data
#######################
if random_subset == True:
	expression_output_file = processed_expression_dir + 'single_cell_raw_expression_sle_individuals_random_subset.txt'
else:
	expression_output_file = processed_expression_dir + 'single_cell_raw_expression_sle_individuals.txt'
np.savetxt(expression_output_file, adata.X.toarray(), fmt="%s", delimiter='\t')

adata.raw = adata

######################
# Normalize data
#######################
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

######################
# Run PCA on data
#######################
sc.tl.pca(adata, svd_solver='arpack')


#######################
# Save Gene IDs
#######################
if random_subset == True:
	gene_id_output_file = processed_expression_dir + 'single_cell_expression_sle_individuals_random_subset_gene_ids.txt'
else:
	gene_id_output_file = processed_expression_dir + 'single_cell_expression_sle_individuals_gene_ids.txt'
np.savetxt(gene_id_output_file, np.vstack((adata.var.index, adata.var[adata.var.columns[0]])).T, fmt="%s", delimiter='\t')

#######################
# Save expression data
#######################
if random_subset == True:
	expression_output_file = processed_expression_dir + 'single_cell_expression_sle_individuals_random_subset_standardized.txt'
else:
	expression_output_file = processed_expression_dir + 'single_cell_expression_sle_individuals_standardized.txt'
np.savetxt(expression_output_file, adata.X, fmt="%s", delimiter='\t')

#######################
# Save Covariate Info
#######################
adata.obs['cell_id'] = adata.obs.index
if random_subset == True:
	covariate_output_file = processed_expression_dir + 'cell_covariates_sle_individuals_random_subset.txt'
else:
	covariate_output_file = processed_expression_dir + 'cell_covariates_sle_individuals.txt'
np.savetxt(covariate_output_file, adata.obs, fmt="%s", delimiter='\t', header='\t'.join(adata.obs.columns), comments='')

#######################
# Save PCA Scores
#######################
num_pcs=200
if random_subset == True:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_sle_individuals_random_subset.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_explained_sle_individuals_random_subset.txt'
else:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_sle_individuals.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_explained_sle_individuals.txt'
generate_pca_scores_and_variance_explained(expression_output_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file)
#######################
# Save Output h5 file
#######################
if random_subset == True:
	h5_output_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset.h5ad'
else:
	h5_output_file = processed_expression_dir + 'scanpy_processed_single_cell_data.h5ad'
adata.write(h5_output_file)

'''


#######################
# Make cell type specific expression data
#######################
# Load in data
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset.h5ad'
else:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data.h5ad'
adata = sc.read_h5ad(processed_single_cell_h5_file)
# Get unique cell types
unique_cell_types = np.unique(adata.obs.ct_cov)
np.savetxt(processed_expression_dir + 'cell_types.txt', unique_cell_types, fmt="%s")
# Loop through cell types
for cell_type in unique_cell_types:
	print(cell_type)
	printible_cell_type = '_'.join(cell_type.split(' '))
	if random_subset == False:
		standardized_cell_type_expression_file = processed_expression_dir + printible_cell_type + '_single_cell_expression_sle_individuals_standardized.txt'
		standardized_cell_type_covariate_file = processed_expression_dir + printible_cell_type + '_cell_covariates_sle_individuals.txt'
		cell_type_pca_file = processed_expression_dir + printible_cell_type + '_pca_scores_sle_individuals.txt'
		cell_type_pca_pve_file = processed_expression_dir + printible_cell_type + '_pca_variance_explained_sle_individuals.txt'
	else:
		standardized_cell_type_expression_file = processed_expression_dir + printible_cell_type + '_single_cell_expression_sle_individuals_random_subset_standardized.txt'
		standardized_cell_type_covariate_file = processed_expression_dir + printible_cell_type + '_cell_covariates_sle_individuals_random_subset.txt'
		cell_type_pca_file = processed_expression_dir + printible_cell_type + '_pca_scores_sle_individuals_random_subset.txt'
		cell_type_pca_pve_file = processed_expression_dir + printible_cell_type + '_pca_variance_explained_sle_individuals_random_subset.txt'
	#generate_cell_type_sc_expression_data_wrapper(cell_type, adata, standardized_cell_type_expression_file, standardized_cell_type_covariate_file)
	num_pcs=200
	generate_pca_scores_and_variance_explained(standardized_cell_type_expression_file, num_pcs, cell_type_pca_file, cell_type_pca_pve_file)


'''
#######################
# Make Pseudo-bulk expression data
#######################
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset.h5ad'
	raw_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_raw_expression_sle_individuals_random_subset.txt'
	standardized_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_expression_sle_individuals_random_subset_standardized.txt'
	pseudobulk_covariate_file = processed_expression_dir + 'pseudobulk_covariates_sle_individuals_random_subset.txt'
else:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data.h5ad'
	raw_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_raw_expression_sle_individuals.txt'
	standardized_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_expression_sle_individuals_standardized.txt'
	pseudobulk_covariate_file = processed_expression_dir + 'pseudobulk_covariates_sle_individuals.txt'

adata = sc.read_h5ad(processed_single_cell_h5_file)
generate_pseudobulk_expression_data_wrapper(adata, raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file, pseudobulk_covariate_file)



#######################
# Save Pseudobulk PCA Scores
#######################
num_pcs=200
if random_subset == True:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_pseudobulk_sle_individuals_random_subset.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_pseudobulk_explained_sle_individuals_random_subset.txt'
else:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_pseudobulk_sle_individuals.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_pseudobulk_explained_sle_individuals.txt'
generate_pca_scores_and_variance_explained(standardized_pseudobulk_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file)

'''




