import numpy as np 
import os
import sys
import pdb

def standardize_expression(normalized_expression_file, standardized_gene_expression_file):
	f = open(normalized_expression_file)
	t = open(standardized_gene_expression_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(data) + '\n')
			continue
		counts = np.asarray(data[1:]).astype(float)
		standardized_counts = (counts - np.mean(counts))/np.std(counts)
		t.write(data[0] + '\t' + '\t'.join(standardized_counts.astype(str)) + '\n')
	f.close()
	t.close()

def recreate_expression(normalized_expression_file, recreated_normalized_gene_expression_file, valid_cell_indices, protein_coding_genes):
	f = open(normalized_expression_file)
	t = open(recreated_normalized_gene_expression_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split())
		if head_count == 0:
			head_count = head_count + 1
			t.write('Gene_id\t' + '\t'.join(data[valid_cell_indices]) + '\n')
			continue
		gene_id = data[0]
		ensamble_id = gene_id.split('_')[0]
		if ensamble_id not in protein_coding_genes:
			continue
		counts = data[1:]
		t.write(gene_id + '\t' + '\t'.join(counts[valid_cell_indices]) + '\n')
	t.close()
	f.close()


def fraction_of_zeros_in_each_gene(normalized_expression_file, fraction_of_zeros_in_each_gene_file):
	f = open(normalized_expression_file)
	t = open(fraction_of_zeros_in_each_gene_file, 'w')
	t.write('Gene_id\tfraction_of_zeros\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		counts = np.asarray(data[1:])
		num_samples = len(counts)
		num_zeros = len(np.where(counts == '0')[0])
		fraction_zeros = float(num_zeros)/num_samples
		t.write(gene_id + '\t' + str(fraction_zeros) + '\n')
	t.close()
	f.close()

def center_expression(normalized_expression_file, centered_normalized_gene_expression_file):
	f = open(normalized_expression_file)
	t = open(centered_normalized_gene_expression_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(data) + '\n')
			continue
		gene_id = data[0]
		counts = np.asarray(data[1:]).astype(float)
		centered_counts = counts - np.mean(counts)
		t.write(gene_id + '\t' + '\t'.join(centered_counts.astype(str)) + '\n')
	t.close()
	f.close()

# Generate expression PC loadings and variance explained of those expression PCs
def generate_pca_scores_and_variance_explained(filtered_standardized_sc_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file):
	# Load in data
	X_full = np.loadtxt(filtered_standardized_sc_expression_file, dtype=str, delimiter='\t', comments='*')
	X = np.transpose(X_full[1:,1:].astype(float))

	# Run PCA (via SVD)
	uuu, sss, vh = np.linalg.svd(np.transpose(X), full_matrices=False)
	svd_loadings = np.transpose(vh)[:,:num_pcs]

	# Save to output file
	np.savetxt(filtered_cells_pca_file, svd_loadings, fmt="%s", delimiter='\t')

	# Compute variance explained
	ve = (np.square(sss)/np.sum(np.square(sss)))[:num_pcs]
	np.savetxt(filtered_cells_pca_ve_file, ve, fmt="%s", delimiter='\n')

def recreate_cell_covariates(meta_data_file, recreated_cell_covariates_file, valid_individuals):
	f = open(meta_data_file)
	t = open(recreated_cell_covariates_file, 'w')
	valid_cells = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('cell_name_header\t' + '\t'.join(data) + '\n')
			continue
		if data[85] in valid_individuals:
			t.write('\t'.join(data) + '\n')
			valid_cells.append(True)
		else:
			valid_cells.append(False)
	f.close()
	t.close()
	return np.asarray(valid_cells)

def make_cell_to_individual_mapping(meta_data_file, cell_to_individual_mapping_file, valid_individuals):
	f = open(meta_data_file)
	t = open(cell_to_individual_mapping_file, 'w')
	head_count = 0
	t.write('cell_name\tindividual_id\n')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[85] in valid_individuals:
			valid_individuals[data[85]] = 1
			t.write(data[0] + '\t' + data[85] + '\n')
	f.close()
	t.close()


def variance_of_each_gene(recreated_normalized_gene_expression_file, variance_of_each_gene_file):
	t = open(variance_of_each_gene_file, 'w')
	t.write('Gene_id\tvariance\n')
	f = open(recreated_normalized_gene_expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		counts = np.asarray(data[1:]).astype(float)
		variance = np.var(counts)
		t.write(gene_id + '\t' + str(variance) + '\n')
	f.close()
	t.close()


def standardize_expression_for_top_nn_variable_genes(recreated_normalized_gene_expression_file, top_n_variable_genes_standardized_gene_expression_file, nn):
	variances = []
	f = open(recreated_normalized_gene_expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		counts = np.asarray(data[1:]).astype(float)
		variances.append(np.var(counts))
	f.close()
	variances = np.asarray(variances)
	variance_thresh = sorted(variances)[-nn]
	f = open(recreated_normalized_gene_expression_file)
	t = open(top_n_variable_genes_standardized_gene_expression_file, 'w')
	head_count = 0
	gene_counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		gene_var = variances[gene_counter]
		gene_counter = gene_counter + 1
		if gene_var >= variance_thresh:
			gene_id = data[0]
			counts = np.asarray(data[1:]).astype(float)
			standardized_counts = (counts - np.mean(counts))/np.std(counts)
			t.write(gene_id + '\t' + '\t'.join(standardized_counts.astype(str)) + '\n')
	f.close()
	t.close()

def compute_maf(genotype):
	af = np.sum(genotype)/(2.0*len(genotype))
	if af > .5:
		maf = 1.0 - af
	else:
		maf = af
	if maf > 0.5 or maf < 0.0:
		print('genotype assumption error')
		pdb.set_trace()
	return maf


def reformat_genotype_file(genotype_file, reformatted_genotype_file, meta_data_file):
	ordered_individuals = []
	individuals = {}
	f = open(genotype_file)
	t = open(reformatted_genotype_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(' ')
		if head_count == 0:
			head_count = head_count + 1
			t.write('variant_id')
			for indi_id in data:
				reformatted_indi_id = indi_id.split('"')[1]
				individuals[reformatted_indi_id] = 0
				ordered_individuals.append(reformatted_indi_id)
				t.write('\t' + reformatted_indi_id)
			t.write('\n')
			continue
		row_id = data[0].split('"')[1]
		chrom_num = row_id.split('--')[0]
		genotype = np.asarray(data[1:]).astype(int)
		maf = compute_maf(genotype)
		# Ignore variants with maf < .05
		if maf >= .05:
			t.write(row_id + '\t'.join(data[1:]) + '\n')
	t.close()
	return individuals

def extract_dictionary_list_of_protein_coding_ensamble_ids(gene_annotation_file):
	dicti = {}
	head_count = 0
	f = open(gene_annotation_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0].split('.')[0]
		gene_type = data[6]
		gene_known = data[7]
		if gene_type != 'protein_coding':
			continue
		if gene_known != 'KNOWN':
			continue
		dicti[ensamble_id] = 1
	f.close()
	return dicti

def generate_per_day_standardized_expression(recreated_normalized_gene_expression_file, recreated_cell_covariates_file, day, per_day_expression_file):
	# Extract binary vector corresponding to cells belonging to this day
	day_specific_cells = []
	head_count = 0
	f = open(recreated_cell_covariates_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[6] != 'day' + str(day):
			day_specific_cells.append(False)
		else:
			day_specific_cells.append(True)
	f.close()
	day_specific_cells = np.asarray(day_specific_cells)
	print(sum(day_specific_cells))
	print(len(day_specific_cells))
	# Stream gene expression file
	# Filter each row to columns that correspond with the current day
	f = open(recreated_normalized_gene_expression_file)
	t = open(per_day_expression_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split())
		if head_count == 0:
			head_count = head_count + 1
			row_name = data[0]
			cell_names = data[1:]
			t.write(row_name + '\t' + '\t'.join(cell_names[day_specific_cells]) + '\n')
			continue
		gene_id = data[0]
		counts = data[1:].astype(float)
		filtered_counts = counts[day_specific_cells]
		normalized_filtered_counts = (filtered_counts - np.mean(filtered_counts))/np.std(filtered_counts)
		t.write(gene_id + '\t' + '\t'.join(normalized_filtered_counts.astype(str)) + '\n')
	f.close()
	t.close()

def generate_per_day_cell_covariates(recreated_cell_covariates_file, day, per_day_cell_covariates_file):
	f = open(recreated_cell_covariates_file)
	t = open(per_day_cell_covariates_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		if data[6] != 'day' + str(day):
			continue
		t.write(line + '\n')
	f.close()
	t.close()

def get_dictionary_list_of_individuals_that_we_have_genotype_for(genotype_file):
	f = open(genotype_file)
	individuals = {}
	head_count = 0
	for line in f:
		if head_count == 0:
			line = line.rstrip()
			data = line.split()
			indis = data[1:]
			for indi in indis:
				individuals[indi] = 0
			head_count = head_count + 1
			continue
	f.close()
	return individuals


normalized_expression_file = sys.argv[1]
meta_data_file = sys.argv[2]
genotype_file = sys.argv[3]
gene_annotation_file = sys.argv[4]
pre_processed_data_dir = sys.argv[5]

##############################
# Extract dictionary list of protein coding ensamble ids
protein_coding_genes = extract_dictionary_list_of_protein_coding_ensamble_ids(gene_annotation_file)


##############################
# Reformat genotype data and get list of individuals that we have genotype data for
valid_individuals = get_dictionary_list_of_individuals_that_we_have_genotype_for(genotype_file)

###############################
# Create mapping from cell-id to individual id
# And filter cells to those for which we have genotype data for
cell_to_individual_mapping_file = pre_processed_data_dir + 'cell_individual_mapping.txt'
make_cell_to_individual_mapping(meta_data_file, cell_to_individual_mapping_file, valid_individuals)


###############################
# Re-create cell covariates
# And filter cells to those for which we have genotype data for
recreated_cell_covariates_file = pre_processed_data_dir + 'cell_covariates.txt'
valid_cell_indices = recreate_cell_covariates(meta_data_file, recreated_cell_covariates_file, valid_individuals)

###############################
# Re-create expression data
# Also filter to only protein coding genes
recreated_normalized_gene_expression_file = pre_processed_data_dir + 'normalized_expression_all_cells.txt'
recreate_expression(normalized_expression_file, recreated_normalized_gene_expression_file, valid_cell_indices, protein_coding_genes)


###############################
# Generate standardized gene expression data in each day independently
# And then create PCs
num_pcs=50
for day in range(4):
	print(day)
	per_day_expression_file = pre_processed_data_dir + 'standardized_normalized_per_day_' + str(day) + '_expression.txt'
	generate_per_day_standardized_expression(recreated_normalized_gene_expression_file, recreated_cell_covariates_file, day, per_day_expression_file)
	pca_loading_file = pre_processed_data_dir + 'standardized_normalized_per_day_' + str(day) + '_expression_pca_loadings.txt'
	pca_pve_file = pre_processed_data_dir + 'standardized_normalized_per_day_' + str(day) + '_expression_pca_pve.txt'
	generate_pca_scores_and_variance_explained(per_day_expression_file, num_pcs, pca_loading_file, pca_pve_file)
	per_day_cell_covariates_file =  pre_processed_data_dir + 'cell_covariates_day_' + str(day) + '.txt'
	generate_per_day_cell_covariates(recreated_cell_covariates_file, day, per_day_cell_covariates_file)


###############################
# Center expression data
centered_normalized_gene_expression_file = pre_processed_data_dir + 'centered_normalized_expression_all_cells.txt'
center_expression(recreated_normalized_gene_expression_file, centered_normalized_gene_expression_file)

###############################
# Standardize expression data
standardized_gene_expression_file = pre_processed_data_dir + 'standardized_normalized_expression_all_cells.txt'
standardize_expression(recreated_normalized_gene_expression_file, standardized_gene_expression_file)

###############################
# Generate expression data for top nn variable genes
nn = 500
top_n_variable_genes_standardized_gene_expression_file = pre_processed_data_dir + 'standardized_normalized_expression_all_cells_top_' + str(nn) + '_variable_genes.txt'
standardize_expression_for_top_nn_variable_genes(recreated_normalized_gene_expression_file, top_n_variable_genes_standardized_gene_expression_file, nn)

###############################
# Generate expression data for top nn variable genes
nn = 2000
top_n_variable_genes_standardized_gene_expression_file = pre_processed_data_dir + 'standardized_normalized_expression_all_cells_top_' + str(nn) + '_variable_genes.txt'
standardize_expression_for_top_nn_variable_genes(recreated_normalized_gene_expression_file, top_n_variable_genes_standardized_gene_expression_file, nn)

###############################
# Compute variance of each gene
variance_of_each_gene_file = pre_processed_data_dir + 'variance_of_each_gene.txt'
variance_of_each_gene(recreated_normalized_gene_expression_file, variance_of_each_gene_file)

###############################
# Compute fraction of zeros in each gene
fraction_of_zeros_in_each_gene_file = pre_processed_data_dir + 'fraction_of_zeros_in_each_gene.txt'
fraction_of_zeros_in_each_gene(recreated_normalized_gene_expression_file, fraction_of_zeros_in_each_gene_file)

###############################
# Run PCA on standardized expression data
num_pcs=50
pca_loading_file = pre_processed_data_dir + 'standardized_normalized_expression_pca_loadings.txt'
pca_pve_file = pre_processed_data_dir + 'standardized_normalized_expression_pca_pve.txt'
generate_pca_scores_and_variance_explained(standardized_gene_expression_file, num_pcs, pca_loading_file, pca_pve_file)

###############################
# Run PCA on centered expression data
num_pcs=50
pca_loading_file = pre_processed_data_dir + 'centered_normalized_expression_pca_loadings.txt'
pca_pve_file = pre_processed_data_dir + 'centered_normalized_expression_pca_pve.txt'
generate_pca_scores_and_variance_explained(centered_normalized_gene_expression_file, num_pcs, pca_loading_file, pca_pve_file)

###############################
# Run PCA on standardized expression data (top n variable genes)
num_pcs=50
nn = 500
pca_loading_file = pre_processed_data_dir + 'standardized_normalized_top_' + str(nn) + '_variable_genes_expression_pca_loadings.txt'
pca_pve_file = pre_processed_data_dir + 'standardized_normalized_top_' + str(nn) + '_variable_genes_expression_pca_pve.txt'
top_n_variable_genes_standardized_gene_expression_file = pre_processed_data_dir + 'standardized_normalized_expression_all_cells_top_' + str(nn) + '_variable_genes.txt'
generate_pca_scores_and_variance_explained(top_n_variable_genes_standardized_gene_expression_file, num_pcs, pca_loading_file, pca_pve_file)

###############################
# Run PCA on standardized expression data (top n variable genes)
num_pcs=50
nn = 2000
pca_loading_file = pre_processed_data_dir + 'standardized_normalized_top_' + str(nn) + '_variable_genes_expression_pca_loadings.txt'
pca_pve_file = pre_processed_data_dir + 'standardized_normalized_top_' + str(nn) + '_variable_genes_expression_pca_pve.txt'
top_n_variable_genes_standardized_gene_expression_file = pre_processed_data_dir + 'standardized_normalized_expression_all_cells_top_' + str(nn) + '_variable_genes.txt'
generate_pca_scores_and_variance_explained(top_n_variable_genes_standardized_gene_expression_file, num_pcs, pca_loading_file, pca_pve_file)
