from sklearn import linear_model
import numpy as np 
import os
import sys
import pdb
import gzip
import random
import pandas as pd
import rnaseqnorm


def get_tissues(file_name):
	f = open(file_name)
	arr = []
	arr2 = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		arr.append(data[0])
		arr2.append(data[1])
	f.close()
	return arr, arr2



def regress_out_covariates(expression_input_file, covariate_file, expression_output_file):
	# Load in covariates
	covariate_raw = np.transpose(np.loadtxt(covariate_file, dtype=str, delimiter='\t'))
	covariate_names = covariate_raw[1:,0]
	covariate_samples = covariate_raw[0,1:]
	covariate_mat = covariate_raw[1:,1:].astype(float)
	
	# Load in expression data
	expression_raw = np.loadtxt(expression_input_file, dtype=str, delimiter='\t')
	expression_samples = expression_raw[1:,0]
	gene_ids = expression_raw[0,1:]
	expr_mat = expression_raw[1:,1:].astype(float)
	# Simple error checking
	if np.array_equal(expression_samples, covariate_samples) == False:
		print('assumption error!')
		pdb.set_trace()

	# Initialize output matrix 
	num_samples = expr_mat.shape[0]
	num_genes = expr_mat.shape[1]

	t = open(expression_output_file, 'w')
	t.write('Gene_id\t' + '\t'.join(expression_samples) + '\n')


	model = linear_model.LinearRegression(fit_intercept=True) 
	modelfit = model.fit(np.transpose(covariate_mat),expr_mat)
	pred = modelfit.predict(np.transpose(covariate_mat))

	resid = expr_mat - pred
	for gene_number in range(num_genes):
		# print(np.std(resid[:,gene_number]))
		#residual_expression = regress_out_covariates_for_one_gene(expr_mat[:,gene_number], covariate_mat)
		gene_id = gene_ids[gene_number]
		t.write(gene_id + '\t' + '\t'.join(resid[:,gene_number].astype(str)) + '\n')
	t.close()
	'''
	# Open expression input and output files
	f = gzip.open(expression_input_file)
	t = open(expression_output_file, 'w')
	# Stream input file line by line
	head_count = 0  # to identify header line
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			# check to make sure expression samples are in the same order as covariate sampels
			if np.array_equal(data[4:],covariate_samples) == False:
				print('Assumption error: samples not alligned')
			# Print header to output file
			t.write('\t'.join(data[3:]) + '\n')
			continue
		gene_id = data[3]
		expression = np.asarray(data[4:]).astype(float)
		residual_expression = regress_out_covariates_for_one_gene(expression, covariate_mat)
		t.write(gene_id + '\t' + '\t'.join(residual_expression.astype(str)) + '\n')
	f.close()
	t.close()
	'''

def get_sample_names(tissues, gtex_expression_dir, sample_name_file):
	# Initialize arr
	samples = []
	# get samples in tisssue
	for tissue in tissues:
		expression_file = gtex_expression_dir + tissue + '.v8.normalized_expression.bed.gz'
		f = gzip.open(expression_file)
		# sample names are only in the header
		head_count = 0
		for line in f:
			if head_count == 0:
				line = line.rstrip()
				data = line.split()
				head_count = head_count +1
				for indi_id in data[4:]:
					samples.append(indi_id + ':' + tissue)
				continue
			break
		f.close()
	# Get mapping from sample_name to index
	sample_to_index = {}
	# Print to output file
	t = open(sample_name_file,'w')
	for i,sample in enumerate(samples):
		t.write(sample + '\n')
		sample_to_index[sample] = i
	t.close()
	return samples, sample_to_index

# Extract array of tests to use
# Extract a dictionary for variants that maps from each variant in tests to an initialized array of length== length(sample_names)
# Extract a dictionary for genes that maps from each genes in tests to an initialized array of length== length(sample_names)
def extract_tests(tissues, gtex_egene_dir, num_samples, genes_tested_in_all_tissues, valid_variants):
	#num_tests_per_tissue = int(round(25000.0/len(tissues)))
	num_tests_per_tissue = 1000
	# Initailize output objects
	tests = []
	variants = {}
	genes = {}
	# Loop through tissues
	for tissue in tissues:
		possible_tests = []
		# Get egene file for this tissue
		egene_file = gtex_egene_dir + tissue + '.v8.egenes.txt'
		# Stream egene file for this tissue
		f = open(egene_file)
		head_count = 0
		data = 'hi'
		for line in f:
			line = line.rstrip()
			prev_data = data
			data = line.split('\t')
			if len(data) != 33:
				continue
			# skip header
			if head_count == 0:
				head_count = head_count + 1
				header = data
				continue
			ensamble_id = data[0]
			# Skip genes not tested in all tissues
			if ensamble_id not in genes_tested_in_all_tissues:
				continue
			snp_id = data[11]
			# Skip variant if not in our list
			if snp_id not in valid_variants:
				continue
			# limit to autosomal chromosomes
			if snp_id.split('_')[0] == 'chrX' or snp_id.split('_')[0] == 'chrY':
				continue
			rs_id = data[18]
			qval = float(data[28])
			if qval < .05:
				possible_tests.append(ensamble_id + ':' + snp_id)
		f.close()
		if len(possible_tests) < num_tests_per_tissue:
			print('assumption error')
			pdb.set_trace()
		tests_in_this_tissue = random.sample(possible_tests, num_tests_per_tissue)
		for test in tests_in_this_tissue:
			ensamble_id = test.split(':')[0]
			snp_id = test.split(':')[1]
			tests.append(ensamble_id + ':' + snp_id)
			variants[snp_id] = np.zeros(num_samples)
			genes[ensamble_id] = np.zeros(num_samples)
	return np.unique(tests), variants, genes

# Print test names to output file
def print_test_names(tests, test_name_file):
	f = open(test_name_file, 'w')
	for test in tests:
		f.write(test + '\n')
	f.close()

# We are going to limit analysis to genes tested in all tissues
def get_genes_tested_in_all_tissues(tissues, gtex_expression_dir):
	# Initialize dictionaries to keep track of which genes were in all tissues
	gene_counts = {}
	valid_genes = {}
	# For each tissue, keep track fo which genes were used
	for tissue in tissues:
		gene_file_name = gtex_expression_dir + tissue + '.v8.normalized_expression.bed.gz'
		head_count = 0
		f = gzip.open(gene_file_name)
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			ensamble_id = data[3]
			# Add gene to dictionary if not observed
			if ensamble_id not in gene_counts:
				gene_counts[ensamble_id] = 0
			# Keep track of how many tissues this gene was observed in 
			gene_counts[ensamble_id] = gene_counts[ensamble_id] + 1
		f.close()
	# Loop through all observed ensamble ids
	# Only take genes that are expressed in all tissues
	for ensamble_id in gene_counts.keys():
		if gene_counts[ensamble_id] == len(tissues):
			valid_genes[ensamble_id] = 1
	return valid_genes

# Fill in 'genes' dictionary with expression values from each tissue
def add_expression_values_to_data_structure(expr_file, sample_to_index, genes):
	f = open(expr_file)
	# to identify header
	head_count = 0
	# Stream file
	for line in f:
		line = line.rstrip()
		data = line.split()
		# header
		if head_count == 0:
			head_count = head_count + 1
			# Creating mapping from index to sample name
			mapping = {}
			for i, indi_id in enumerate(data[1:]):
				mapping[i] = indi_id
			continue
		ensamble_id = data[0]
		# Only consider lines where ensamble_id is in genes dictionary
		if ensamble_id not in genes:
			continue
		expr_vec = np.asarray(data[1:]).astype(float)
		for index, expr in enumerate(expr_vec):
			sample_name = mapping[index] 
			genes[ensamble_id][sample_to_index[sample_name]] = expr
	f.close()
	return genes

# Fill in 'genes' dictionary with expression values from each tissue
def add_expression_values_to_data_structure_t(expr_file, sample_to_index, genes):
	aa = np.loadtxt(expr_file,dtype=str,delimiter='\t')
	np.savetxt('temp.txt', np.transpose(aa), fmt="%s", delimiter='\t')
	f = open('temp.txt')
	# to identify header
	head_count = 0
	# Stream file
	for line in f:
		line = line.rstrip()
		data = line.split()
		# header
		if head_count == 0:
			head_count = head_count + 1
			# Creating mapping from index to sample name
			mapping = {}
			for i, indi_id in enumerate(data[1:]):
				mapping[i] = indi_id
			continue
		ensamble_id = data[0]
		# Only consider lines where ensamble_id is in genes dictionary
		if ensamble_id not in genes:
			continue
		expr_vec = np.asarray(data[1:]).astype(float)
		for index, expr in enumerate(expr_vec):
			sample_name = mapping[index] 
			genes[ensamble_id][sample_to_index[sample_name]] = expr
	f.close()
	return genes


# Fill in 'variants' dictionary with genotype values
def add_genotype_values_to_data_structure(variants, sample_names, gtex_genotype_dir):
	used_variants = {}
	# loop through chromosomes
	for chrom_num in range(1,23):
		genotype_file = gtex_genotype_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + str(chrom_num) + '_dosage_MAF_05.txt'
		head_count = 0
		# Stream genotype file for this chromosome
		f = open(genotype_file)
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				# Header contains indi_ids
				indi_ids = data
				# Create mapping from column index to sample inidices
				mapping = {}
				for index, indi_id in enumerate(indi_ids):
					valid = False
					arr = []
					for index2, sample_name in enumerate(sample_names):
						if sample_name.split(':')[0] == indi_id:
							valid = True
							arr.append(index2)
					if valid == True:
						mapping[index] = np.asarray(arr)
				continue
			snp_id = data[0]
			# skip lines where variant not in 'variants' dictionary
			if snp_id not in variants:
				continue
			try:
				genotype_vec = np.asarray(data[1:]).astype(float)
			except:
				for index, genotype in enumerate(data[1:]):
					if index not in mapping or genotype == '-':
						continue
					variants[snp_id][mapping[index]] = float(genotype)
				continue
			used_variants[snp_id] = 1
			for index, genotype in enumerate(genotype_vec):
				if index not in mapping:
					continue
				variants[snp_id][mapping[index]] = genotype
		f.close()
	return variants


# Fill in 'variants' dictionary with genotype values
def get_variants_we_have_genotype_data_for( gtex_genotype_dir):
	used_variants = {}
	# loop through chromosomes
	for chrom_num in range(1,23):
		genotype_file = gtex_genotype_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + str(chrom_num) + '_dosage_MAF_05.txt'
		head_count = 0
		# Stream genotype file for this chromosome
		f = open(genotype_file)
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			snp_id = data[0]
			valid = True
			for ele in data[1:]:
				if ele == '-':
					valid=False
			if valid == True:
				used_variants[snp_id] = data[0]
		f.close()
	return used_variants

def print_big_matrix(output_file, tests, data_struct, index):
	t = open(output_file,'w')
	for test in tests:
		name = test.split(':')[index]
		t.write('\t'.join(data_struct[name].astype(str)) + '\n')
	t.close()

def print_individual_id_file_and_tissue_file(sample_names, output_file, output_tissue_file):
	# Create mapping for indi_id to index
	counter = 0
	mapping = {}
	for sample_name in sample_names:
		indi_id = sample_name.split(':')[0]
		if indi_id not in mapping:
			mapping[indi_id] = counter
			counter = counter  + 1
	t = open(output_file, 'w')
	for sample_name in sample_names:
		indi_id = sample_name.split(':')[0]
		index = mapping[indi_id]
		t.write(str(index) + '\n')
	t.close()
	t = open(output_tissue_file, 'w')
	for sample_name in sample_names:
		tissue_name = sample_name.split(':')[1]
		t.write(tissue_name + '\n')
	t.close()


def remove_genes_with_zero_expression(tpm_matrix, ordered_genes):
	valid_indices = []
	for gene_number in range(len(ordered_genes)):
		if np.sum(tpm_matrix[:,gene_number]) != 0:
			valid_indices.append(gene_number)
	valid_indices = np.asarray(valid_indices)
	return tpm_matrix[:,valid_indices], np.asarray(ordered_genes)[valid_indices]

# Generate TPM expression matrix
def generate_tpm_expression_matrix(tissues, tissues_alt, ordered_sample_names, genes_tested_in_all_tissues, gtex_tpm_dir, output_file):
	# Generate mapping from sample_name to row number
	sample_name_to_row_number = {}
	for index, sample_name in enumerate(ordered_sample_names):
		sample_name_to_row_number[sample_name] = index
	# Generate mapping from gene to column number
	ordered_genes = sorted(genes_tested_in_all_tissues.keys())
	gene_to_column_number = {}
	for index, gene in enumerate(ordered_genes):
		gene_to_column_number[gene] = index
	# Num samples and number of genes
	num_samples = len(sample_name_to_row_number)
	num_genes = len(ordered_genes)
	# Initialize tpm matrix
	tpm_matrix = np.zeros((num_samples, num_genes))
	counter = 0
	used = {}
	# Loop through each tissue and fill in the tpm matrix
	for tissue_index, tissue in enumerate(tissues):
		tissue_alt = tissues_alt[tissue_index]
		# Stream tpm file in this tissue
		tpm_file = gtex_tpm_dir + tissue_alt + '.txt'
		head_count = 0
		f = open(tpm_file)
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Header
			if head_count == 0:
				head_count = head_count + 1
				# Get ordered list of sample names
				sample_names = []
				for sample_name_temp in data[2:]:
					sample_name = sample_name_temp.split('-')[0] + '-' + sample_name_temp.split('-')[1] + ':' + tissue
					sample_names.append(sample_name)
				continue
			# Get relevent fields from line
			ensamble_id = data[0]
			# Skip genes we are not interested in
			if ensamble_id not in gene_to_column_number:
				continue
			tpm_counts = np.asarray(data[2:]).astype(float)
			# Loop through samples and add sample/tpm count to tpm_matrix
			for index, sample_name in enumerate(sample_names):
				# Ignore samples that aren't in our list
				if sample_name not in sample_name_to_row_number:
					continue
				# Get row corresponding to this smample
				row_index = sample_name_to_row_number[sample_name]
				# Get column corresponding to this gene
				column_index = gene_to_column_number[ensamble_id]

				ele_name = str(row_index) + '_' + str(column_index)

				# Add tpm count to matrix
				tpm_matrix[row_index, column_index] = tpm_counts[index]
				counter = counter + 1
		f.close()
	# Remove genes with zero expression across all samples
	filtered_tpm_matrix, filtered_orderd_genes = remove_genes_with_zero_expression(tpm_matrix, ordered_genes)
	log_filtered_tpm_matrix = np.log2(filtered_tpm_matrix + 1.0)
	# Print to output file
	t = open(output_file, 'w')
	# print header
	t.write('GeneId\t' + '\t'.join(filtered_orderd_genes) + '\n')
	for sample_num, sample_name in enumerate(ordered_sample_names):
		tpm_expr = log_filtered_tpm_matrix[sample_num,:].astype(str)
		t.write(sample_name + '\t' + '\t'.join(tpm_expr) + '\n')
	t.close()


def standardize_expression(tpm_expression_matrix_file, standardized_tpm_expression_matrix_file):
	tpm_full = np.loadtxt(tpm_expression_matrix_file, dtype=str,delimiter='\t')
	tpm = tpm_full[1:,1:].astype(float)
	samples = tpm_full[1:,0]
	genes = tpm_full[0,1:]
	# Quantile normalize the samples
	df = pd.DataFrame(np.transpose(tpm))
	#rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
	#temp_out = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
	#tpm_quantile_normalized = np.transpose(np.asarray(temp_out))
	temp_out = rnaseqnorm.normalize_quantiles(df)
	norm_df = rnaseqnorm.inverse_normal_transform(temp_out)
	standardized_tpm = np.transpose(np.asarray(norm_df))

	###
	#tpm_quantile_normalized = np.transpose(np.asarray(temp_out))
	###

	# Standardize the genes
	#num_genes = tpm_quantile_normalized.shape[1]
	#num_samples = tpm_quantile_normalized.shape[0]

	####
	#standardized_tpm = np.zeros((num_samples, num_genes))
	#for gene_num in range(num_genes):
	#	standardized_tpm[:,gene_num] = (tpm_quantile_normalized[:, gene_num] - np.mean(tpm_quantile_normalized[:, gene_num]))/np.std(tpm_quantile_normalized[:, gene_num])
	####
	# Print to output file
	t = open(standardized_tpm_expression_matrix_file, 'w')
	# print header
	t.write('GeneId\t' + '\t'.join(genes) + '\n')
	for sample_num, sample_name in enumerate(samples):
		#expr = tpm_quantile_normalized[sample_num, :].astype(str)
		###
		expr = standardized_tpm[sample_num, :].astype(str)
		###
		t.write(sample_name + '\t' + '\t'.join(expr) + '\n')
	t.close()

def extract_covariates(expr_file, output_file, num_expression_pcs):
	# Load in expression data
	expr_full = np.loadtxt(expr_file, dtype=str,delimiter='\t')
	expr = expr_full[1:,1:].astype(float)
	samples = expr_full[1:,0]
	genes = expr_full[0,1:]

	# Compute gene expression pcs on expression data
	uuu, sss, vh = np.linalg.svd(np.transpose(expr))
	expr_pc_loadings = np.transpose(vh)[:,:num_expression_pcs]
	# print to output file 
	pc_names = []
	for pc_num in range(num_expression_pcs):
		pc_names.append('PC' + str(pc_num))
	t = open(output_file, 'w')
	t.write('Sample_name\t' + '\t'.join(pc_names) + '\n')
	for sample_num, sample_name in enumerate(samples):
		t.write(sample_name + '\t' + '\t'.join(expr_pc_loadings[sample_num,:].astype(str)) + '\n')
	t.close()


def get_genes_we_have_expression_data_for(file_name):
	genes = {}
	head_count = 0
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes[data[0]] = 1
	return genes



tissues_file = sys.argv[1]
gtex_expression_dir = sys.argv[2]
gtex_tpm_dir = sys.argv[3]
gtex_covariate_dir = sys.argv[4]
gtex_genotype_dir = sys.argv[5]
gtex_egene_dir = sys.argv[6]
processed_data_dir = sys.argv[7]

tissues, tissues_alt = get_tissues(tissues_file)


# Extract file of sample names
sample_name_file = processed_data_dir + 'sample_names.txt'
sample_names, sample_to_index = get_sample_names(tissues, gtex_expression_dir, sample_name_file)

print_individual_id_file_and_tissue_file(sample_names, processed_data_dir + 'individual_id.txt', processed_data_dir + 'sample_tissue_names.txt')

# We are going to limit analysis to genes tested in all tissues
genes_tested_in_all_tissues = get_genes_tested_in_all_tissues(tissues, gtex_expression_dir)

# Generate TPM expression matrix
tpm_expression_matrix_file = processed_data_dir + 'cross_tissue_tpm.txt'
generate_tpm_expression_matrix(tissues, tissues_alt, sample_names, genes_tested_in_all_tissues, gtex_tpm_dir, tpm_expression_matrix_file)

# Quantile normalize and standardize TPM expression matrix
standardized_tpm_expression_matrix_file = processed_data_dir + 'cross_tissue_tpm_standardized.txt'
standardize_expression(tpm_expression_matrix_file, standardized_tpm_expression_matrix_file)





# Extract covariates (expression pcs)
num_expression_pcs = 50
covariate_file = processed_data_dir + 'covariates.txt'
extract_covariates(standardized_tpm_expression_matrix_file, covariate_file, num_expression_pcs)

# Regress out covariates
residual_standardized_tpm_expression_matrix_file = processed_data_dir + 'residual_cross_tissue_tpm_standardized.txt'
regress_out_covariates(standardized_tpm_expression_matrix_file, covariate_file, residual_standardized_tpm_expression_matrix_file)

# Limit to genes in our analysis
valid_genes = get_genes_we_have_expression_data_for(residual_standardized_tpm_expression_matrix_file)

# Limit to variants we have genoytpe data for
valid_variants = get_variants_we_have_genotype_data_for(gtex_genotype_dir)



# Extract array of tests to use
# Extract a dictionary for variants that maps from each variant in tests to an initialized array of length== length(sample_names)
# Extract a dictionary for genes that maps from each genes in tests to an initialized array of length== length(sample_names)
tests, variants, genes = extract_tests(tissues, gtex_egene_dir, len(sample_names), valid_genes, valid_variants)
genes_uncorrected = genes.copy()


# Print test names to output file
test_name_file = processed_data_dir + 'test_names.txt'
print_test_names(tests, test_name_file)


# Fill in 'genes' dictionary with expression values from each tissue
genes = add_expression_values_to_data_structure(residual_standardized_tpm_expression_matrix_file, sample_to_index, genes)


# Fill in 'variants' dictionary with genotype values
variants = add_genotype_values_to_data_structure(variants, sample_names, gtex_genotype_dir)


print_big_matrix(processed_data_dir + 'expr.txt', tests, genes, 0)

print_big_matrix(processed_data_dir + 'genotype.txt', tests, variants, 1)




# Do the same for un-corrected gene expression
genes_uncorrected = add_expression_values_to_data_structure_t(standardized_tpm_expression_matrix_file, sample_to_index, genes_uncorrected)
print_big_matrix(processed_data_dir + 'expr_uncorrected.txt', tests, genes_uncorrected, 0)

















################################
# OLD (NOT USED CURRENTLY)
################################

'''
# For each tissue, regress out covariates, and then standardized expression data
for tissue in tissues:
	print(tissue)
	expression_input_file = gtex_expression_dir + tissue + '.v8.normalized_expression.bed.gz'
	expression_output_file = processed_data_dir + tissue + '.v8.normalized_expression_covariate_regressed.bed'
	covariate_file= gtex_covariate_dir + tissue + '.v8.covariates.txt'
	regress_out_covariates(expression_input_file, covariate_file, expression_output_file)
'''
