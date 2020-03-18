import numpy as np 
import os
import sys
import pdb
import h5py



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
# generate_cell_covariate_file(data, cell_covariate_file)

# Print expression data to a readable file
sc_expression_file = processed_expression_dir + 'single_cell_expression_ye_lab.txt'
# generate_expression_file(data, sc_expression_file)

# Print gene_names to a readable file
gene_names_file = processed_expression_dir + 'gene_names_ye_lab.txt'
#generate_gene_names_file(data, gene_names_file)

# Get cell type colors and print to file
cell_type_colors_file = processed_expression_dir + 'cell_type_colors_ye_lab.txt'
#print_cell_type_colors_to_file(data, cell_type_colors_file)

# Grab UMAP scores from ye-lab data structure
umap_file = processed_expression_dir + 'umap_scores_ye_lab.txt'
# generate_umap_from_ye_lab_file(data, umap_file)

# Grab UMAP scores from ye-lab data structure
pca_file = processed_expression_dir + 'pca_scores_ye_lab.txt'
generate_pca_from_ye_lab_file(data, pca_file)

