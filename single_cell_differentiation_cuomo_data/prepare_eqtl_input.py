import numpy as np 
import os
import sys
import pdb
from sklearn import linear_model
import h5py
import statsmodels.api as sm




# Get mapping from genes to (chrom_num, position)
def get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes):
	# Convert gene array to dictionary
	gene_mapping = {}
	for gene in genes:
		gene_mapping[gene] = [0,0]
	# Fill in gene positions
	f = open(gene_annotation_file)
	used_genes = {}
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
		if gene_id not in gene_mapping:
			continue
		used_genes[gene_id] = 1
		# Some more error checks
		start = int(data[2])
		end = int(data[3])
		if start > end:
			print('assumption errror in processing gene annotation file')
			pdb.set_trace()
		if data[6] != 'protein_coding' or data[7] != 'KNOWN':
			continue
		if data[1] == 'chrX' or data[1] == 'chrY':
			continue
		# Extract relevent info on gene: Chrom num and TSS
		chrom_num = int(data[1].split('hr')[1])
		strand = data[4]
		gene_id = data[5]
		if strand == '+':
			tss = start
		elif strand == '-':
			tss = end
		else:
			print('assumption error while processing gene annotation file')
		# Add info to gene dictionary
		if gene_mapping[gene_id][0] != 0:
			if gene_mapping[gene_id][0] != chrom_num and gene_mapping[gene_id][1] != tss:
				gene_mapping.pop(gene_id)
		else:
			gene_mapping[gene_id] = (chrom_num, tss)
	f.close()
	print(len(used_genes))
	return gene_mapping

# Create array where each element is a BP in this chromosome
# 'Null' if no genes in distance BP of gene
# Otherwise is a list of gene names
def create_gene_chromsome(chrom_num, gene_mapping, distance):
	# Initialize chromosome
	chromosome = ['Null']*250000000
	# Loop through genes
	for gene_id in gene_mapping.keys():
		# Extract dictionary on position of gene
		gene_tuple = gene_mapping[gene_id]
		gene_chrom_num = gene_tuple[0]
		gene_tss = gene_tuple[1]
		# Make sure gene is on correct chromosome
		if gene_chrom_num != chrom_num:
			continue
		# Iterate across positions in a distance KB window around gene and add gene to chromosome
		for pos in range((gene_tss-distance),(gene_tss+distance+1)):
			if chromosome[pos] == 'Null':
				chromosome[pos] = gene_id
			else:
				chromosome[pos] = chromosome[pos] + ':' + gene_id
	return chromosome

def get_maf(genotype):
	af = np.sum(genotype)/(2.0*len(genotype))
	if af > .5:
		maf = 1.0 - af
	else:
		maf = af
	if maf > .5 or maf < 0.0:
		print('assumption error')
		pdb.set_trace()
	return maf

########################
# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
########################
def extract_variant_gene_pairs_for_eqtl_testing(gene_file, gene_annotation_file, distance, genotype_data_dir, variant_gene_pair_file):
	# Extract gene list
	genes = np.loadtxt(gene_file, delimiter='\t',dtype=str)[:,0]
	print(len(genes))
	# Get mapping from genes to (chrom_num, position)
	gene_mapping = get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes)
	# Open file handle to output file containing variant-gene pair tests and print header
	t = open(variant_gene_pair_file, 'w')
	t.write('Gene_id\tvariant_id\tchrom_num\tgene_tss\tvariant_position\n')

	# Fill in file containing lists of variant gene pairs for each chromosome iteratively
	used_genes = {}
	for chrom_num in range(1,23):
		print(chrom_num)
		# Create array where each element is a BP in this chromosome
		# 'Null' if no genes in distance BP of gene
		# Otherwise is a list of gene names
		chromosome = create_gene_chromsome(chrom_num, gene_mapping, distance)
		# Now loop through variants on this chromosome
		genotype_file = genotype_data_dir + 'chr' + str(chrom_num) + '.genotypes.matrix.eqtl.txt'
		f = open(genotype_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			if len(data) != 120:
				print('assumption error!')
			variant_id = data[0]
			variant_chrom = int(variant_id.split(':')[0])
			variant_pos = int(variant_id.split(':')[1])
			# Simple error check
			if variant_chrom != chrom_num:
				print('assumption error')
				pdb.set_trace()
			# No genes within 10KB of variant
			if chromosome[variant_pos] == 'Null':
				continue
			genotype = np.asarray(data[1:]).astype(float)
			maf = get_maf(genotype)
			if maf < .05:
				print('skipped variant')
				continue
			# List of genes that variant maps to
			mapping_genes = chromosome[variant_pos].split(':')
			for gene_id in mapping_genes:
				# THIS IS A VARIANT-GENE PAIR WE WILL TEST
				# PRINT TO OUTPUT
				used_genes[gene_id] = 1
				t.write(gene_id + '\t' + variant_id + '\t' + str(chrom_num) + '\t' + str(gene_mapping[gene_id][1]) + '\t' + str(variant_pos) + '\n')
		f.close()
	t.close()

# Helper function to load expression file
def load_expression_file(expression_file):
	arr = []
	f = open(expression_file)
	counter=0
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split()).astype(float)
		arr.append(data)
		counter = counter + 1
	expr = np.asmatrix(arr)
	return expr

# First extract list of variants
def extract_variants_from_variant_gene_pair_file(variant_gene_pair_file):
	f = open(variant_gene_pair_file)
	variants = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		variants[data[1]] = 1
	f.close()
	return variants

# Create mapping from variants in variat_list to genotype vectors
def create_mapping_from_variants_to_genotype(variant_list, genotype_file):
	variants = {}
	f = open(genotype_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[0]
		# Limit to variants in variant_list
		if variant_id not in variant_list:
			continue
		variants[variant_id] = np.asarray(data[1:]).astype(float)
	f.close()
	if len(variants) != len(variant_list):
		print('assumption error')
		pdb.set_trace()
	return variants

# Create list of gene_variant pairs after pruning
def create_mapping_from_gene_to_variants_in_variant_gene_pair_file(variant_gene_pair_file):
	f = open(variant_gene_pair_file)
	genes = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		variant_id = data[1]
		if gene_id not in genes:
			genes[gene_id] = []
		genes[gene_id].append(variant_id)
	f.close()
	return genes

def get_list_of_eligable_variants(gene_mapped_variants, gene_mapped_ld_pruned_variants, variants, r_squared_threshold):
	# Initialize list of eligable variants
	eligable_variants = []
	for variant in gene_mapped_variants:
		variant_passes = True
		for variant_in_set in gene_mapped_ld_pruned_variants:
			r_squared = np.square(np.corrcoef(variants[variant], variants[variant_in_set])[0,1])
			if r_squared >= r_squared_threshold:
				variant_passes = False
		if variant_passes == True:
			eligable_variants.append(variant)
	return eligable_variants

def get_list_of_pruned_gene_variant_pairs(gene_to_variant_mapping_unpruned, variants, r_squared_threshold):
	# Initialize list
	pruned_gene_variant_pairs = {}
	# Loop through genes
	for i, gene_id in enumerate(gene_to_variant_mapping_unpruned.keys()):
		# Get list of all variants mapped to this gene
		gene_mapped_variants = gene_to_variant_mapping_unpruned[gene_id]
		# Initialize list of ld-pruned variants mapped to this gene
		gene_mapped_ld_pruned_variants = []
		# randomly select first variant to include
		gene_mapped_ld_pruned_variants.append(np.random.choice(gene_mapped_variants))
		# boolean variable to keep track of if there are more_variants variants to add
		more_variants = True 

		# While loop where each iteration adds a variant until it can't
		while more_variants == True:
			# Get list of variants not in high ld with any variant in gene_mapped_ld_pruned_variants
			eligable_variants = get_list_of_eligable_variants(gene_mapped_variants, gene_mapped_ld_pruned_variants, variants, r_squared_threshold)
			# End the loop
			if len(eligable_variants) == 0:
				more_variants = False
			else:
				gene_mapped_ld_pruned_variants.append(np.random.choice(eligable_variants))
		# Add variant gene pairs to a larger list
		for variant in gene_mapped_ld_pruned_variants:
			test_id = gene_id + '_' + variant 
			pruned_gene_variant_pairs[test_id] = 1
	return pruned_gene_variant_pairs



########################
# Step 2: LD prune the above file in each gene (ie limit to only independent snps per gene)
########################
def ld_prune_variant_gene_pair_file(variant_gene_pair_file, ld_pruned_variant_gene_pair_file, r_squared_threshold, genotype_file, random_seed):
	np.random.seed(random_seed)
	# First extract list of variants
	variant_list = extract_variants_from_variant_gene_pair_file(variant_gene_pair_file)
	# Second extract mapping from genes to array of variants
	gene_to_variant_mapping_unpruned = create_mapping_from_gene_to_variants_in_variant_gene_pair_file(variant_gene_pair_file)
	# Create mapping from variants in variat_list to genotype vectors
	variants = create_mapping_from_variants_to_genotype(variant_list, genotype_file)
	# Create list of gene_variant pairs after pruning
	gene_variant_pairs_pruned = get_list_of_pruned_gene_variant_pairs(gene_to_variant_mapping_unpruned, variants, r_squared_threshold)
	
	# Print to output file
	f = open(variant_gene_pair_file)
	t = open(ld_pruned_variant_gene_pair_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		test_id = data[0] + '_' + data[1]
		if test_id not in gene_variant_pairs_pruned:
			continue
		t.write(line + '\n')
	f.close()
	t.close()

def get_ordered_list_of_gene_names(ld_pruned_variant_gene_pair_file):
	gene_names = []
	f = open(ld_pruned_variant_gene_pair_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_names.append(data[0])
	return np.asarray(gene_names)

def get_ordered_list_of_variant_names(ld_pruned_variant_gene_pair_file):
	variant_names = []
	dicti = {}
	f = open(ld_pruned_variant_gene_pair_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_names.append(data[1])
		dicti[data[1]] = 1
	return np.asarray(variant_names), dicti

def create_mapping_from_gene_names_to_expression_vectors(sc_expression_file):
	f = open(sc_expression_file)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0].split('_')[0]
		gene_expression = np.asarray(data[1:])
		mapping[gene_id] = gene_expression
	f.close()
	'''
	# Get gene names
	gene_names = np.loadtxt(gene_names_file,dtype=str, delimiter='\t')[:,0]
	print(gene_names)
	# Load in expression matrix (every column corresponds to a gene)
	expression_matrix = np.loadtxt(sc_expression_file, dtype=str, delimiter='\t')
	# Create mapping
	mapping = {}
	for index, gene_name in enumerate(gene_names):
		mapping[gene_name] = expression_matrix[:, index]
	'''
	return mapping

def generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, sc_expression_file, single_cell_expression_eqtl_traing_data_file):
	# Get ordered list of gene names (this will be the order that the output file will be saved in)
	ordered_gene_names = get_ordered_list_of_gene_names(ld_pruned_variant_gene_pair_file)
	# Create mapping from gene names to expression vectors
	gene_name_to_expression_vector = create_mapping_from_gene_names_to_expression_vectors(sc_expression_file)
	# print to output file
	t = open(single_cell_expression_eqtl_traing_data_file, 'w')
	for gene_name in ordered_gene_names:
		# Use map to get expression vector corresponding to this gene
		expression_vector = gene_name_to_expression_vector[gene_name]
		# Print to output file
		t.write('\t'.join(expression_vector) + '\n')
	t.close()


def extract_covariates_from_cell_level_info_file(cell_level_info_file):
	f = open(cell_level_info_file)
	batch_dicti = {}
	well_dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		batch_dicti[data[5]] = 1
	f.close()
	num_batches = len(batch_dicti)
	batch_mapping = {}
	counter = 0
	for batch_name in batch_dicti.keys():
		vec = np.zeros(num_batches - 1)
		if counter < num_batches -1:
			vec[counter] = 1.0
		batch_mapping[batch_name] = vec
		counter = counter + 1

	f = open(cell_level_info_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			print(data)
			continue
		line_arr = []
		if data[2] == 'WHITE':
			line_arr.append(0.0)
		else:
			line_arr.append(1.0)
		line_arr.append(float(data[8]))
		line_arr.append(float(data[9]))
		line_arr.append(float(data[10]))
		line_arr.append(batch_mapping[data[5]])
		line_arr = np.hstack(line_arr)
		arr.append(line_arr)
	f.close()
	return np.asarray(arr)

def regress_out_covariates(expression_input_file, covariate_file, expression_output_file, num_pcs):
	# Load in covariates
	covariate_raw = np.loadtxt(covariate_file, dtype=str, delimiter='\t')

	pc_covariate_mat = covariate_raw[:,:num_pcs].astype(float)
	
	#cell_level_covariate_mat = extract_covariates_from_cell_level_info_file(cell_level_info_file)
	covariate_mat = pc_covariate_mat
	#covariate_mat = np.hstack((pc_covariate_mat, cell_level_covariate_mat))
	f = open(expression_input_file)
	t = open(expression_output_file, 'w')

	counter = 0
	for line in f:
		if np.mod(counter, 50) == 0:
			print(counter)
		counter = counter + 1
		line = line.rstrip()
		data = line.split()
		# Extract expression (each line corresponds to a gene)
		expr = np.asarray(data).astype(float)
		# Fit linear model of covariates onto expression
		model = linear_model.LinearRegression(fit_intercept=True) 
		modelfit = model.fit(covariate_mat, expr)
		# Get predicted expression according to the LM
		pred = modelfit.predict(covariate_mat)
		# Get residual expression
		resid = expr - pred
		# Print residual expression to output file
		t.write('\t'.join(resid.astype(str)) + '\n')
	f.close()
	t.close()

	# Load in expression data
	#expression_raw = np.loadtxt(expression_input_file, dtype=str, delimiter='\t')
	#expression_samples = expression_raw[1:,0]
	#gene_ids = expression_raw[0,1:]
	# This is num_samplesXnum_genes
	#expr_mat = np.transpose(expression_raw.astype(float))


	# Initialize output matrix 
	#num_samples = expr_mat.shape[0]
	#num_genes = expr_mat.shape[1]

	#t = open(expression_output_file, 'w')
	#t.write('Gene_id\t' + '\t'.join(expression_samples) + '\n')


	#model = linear_model.LinearRegression(fit_intercept=True) 
	#modelfit = model.fit(np.transpose(covariate_mat),expr_mat)
	#pred = modelfit.predict(np.transpose(covariate_mat))

	#resid = expr_mat - pred
	#for gene_number in range(num_genes):
		# print(np.std(resid[:,gene_number]))
		#residual_expression = regress_out_covariates_for_one_gene(expr_mat[:,gene_number], covariate_mat)
		#gene_id = gene_ids[gene_number]
		#t.write('\t'.join(resid[:,gene_number].astype(str)) + '\n')
	#t.close()

def get_cell_level_ordered_individaul_array(cell_level_info_file):
	array = []
	head_count = 0
	f = open(cell_level_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		array.append(data[85])
	f.close()
	return np.asarray(array)

def get_genotype_level_ordered_individual_array(genotype_file):
	f = open(genotype_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			line = line.rstrip()
			data = line.split()
			indi = data[1:]
			head_count = head_count + 1
			continue
	f.close()
	return np.asarray(indi)

def construct_standardized_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file):
	variant_names, variant_list = get_ordered_list_of_variant_names(ld_pruned_variant_gene_pair_file)

	ordered_individuals_cell_level = get_cell_level_ordered_individaul_array(cell_level_info_file)
	
	ordered_individuals_genotype_level = get_genotype_level_ordered_individual_array(genotype_file)

	# Create mapping from variant_id to genotype vectors
	variant_to_genotype = create_mapping_from_variants_to_genotype(variant_list, genotype_file)
	
	# Create array mapping from genotype-level array to single cell-level array
	mapping_array = []
	converter = {}
	for i, indi in enumerate(ordered_individuals_genotype_level):
		converter[indi] = i 
	for indi in ordered_individuals_cell_level:
		mapping_array.append(converter[indi])
	mapping_array = np.asarray(mapping_array)
	# Simple error checking
	if np.array_equal(ordered_individuals_cell_level, ordered_individuals_genotype_level[mapping_array]) == False:
		print('assumption error')
		pdb.set_trace()

	# print new genotype file
	t = open(single_cell_genotype_eqtl_training_data_file, 'w')
	for variant_id in variant_names:
		genotype_vector = (variant_to_genotype[variant_id]).astype(str)
		cell_level_genotype_vector = genotype_vector[mapping_array]
		cell_level_genotype_vector = cell_level_genotype_vector.astype(float)
		standardized_cell_level_genotype_vector = (cell_level_genotype_vector - np.mean(cell_level_genotype_vector))/np.std(cell_level_genotype_vector)
		t.write('\t'.join(standardized_cell_level_genotype_vector.astype(str)) + '\n')
	t.close()


def construct_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file):
	variant_names, variant_list = get_ordered_list_of_variant_names(ld_pruned_variant_gene_pair_file)

	ordered_individuals_cell_level = get_cell_level_ordered_individaul_array(cell_level_info_file)
	
	ordered_individuals_genotype_level = get_genotype_level_ordered_individual_array(genotype_file)

	# Create mapping from variant_id to genotype vectors
	variant_to_genotype = create_mapping_from_variants_to_genotype(variant_list, genotype_file)
	
	# Create array mapping from genotype-level array to single cell-level array
	mapping_array = []
	converter = {}
	for i, indi in enumerate(ordered_individuals_genotype_level):
		converter[indi] = i 
	for indi in ordered_individuals_cell_level:
		mapping_array.append(converter[indi])
	mapping_array = np.asarray(mapping_array)
	# Simple error checking
	if np.array_equal(ordered_individuals_cell_level, ordered_individuals_genotype_level[mapping_array]) == False:
		print('assumption error')
		pdb.set_trace()

	# print new genotype file
	t = open(single_cell_genotype_eqtl_training_data_file, 'w')
	for variant_id in variant_names:
		genotype_vector = (variant_to_genotype[variant_id]).astype(str)
		cell_level_genotype_vector = genotype_vector[mapping_array]
		t.write('\t'.join(cell_level_genotype_vector) + '\n')
	t.close()

# Step 6: Generate individual id file (z matrix in eqtl factorization)
def generate_individual_id_file(cell_level_info_file, single_cell_individual_id_file):
	# Get array of length number of cells and each element corresponds to the individual id corresponding to that cell
	indi_ids = get_cell_level_ordered_individaul_array(cell_level_info_file)
	# Creating mapping from individual id to a number
	unique_indi_ids = np.unique(indi_ids)
	indi_to_number = {}
	for i, indi_id in enumerate(unique_indi_ids):
		indi_to_number[indi_id] = i
	# print to output file
	t = open(single_cell_individual_id_file, 'w')
	for indi_id in indi_ids:
		number = indi_to_number[indi_id]
		t.write(str(number) + '\n')
	t.close()

def save_as_h5_file(input_text_file):
	f = open(input_text_file)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		arr.append(data)
	f.close()
	mat = np.asarray(arr).astype(float)
	# save as h5 file
	h5_file_name = input_text_file.split('.tx')[0] + '.h5'
	h5f = h5py.File(h5_file_name, 'w')
	h5f.create_dataset('data', data=mat)
	h5f.close()

def extract_nominal_sig_variant_gene_pairs_from_known_cell_types(single_cell_eqtl_dir, nominal_p, variant_gene_pair_file, expression_file):
	# Extract dictionary list of gene names
	genes = {}
	f = open(expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes[data[0].split('_')[0]] = 1
	f.close()
	# Loop through each time step
	t = open(variant_gene_pair_file, 'w')
	t.write('Gene_id\tvariant_id\tchrom_num\tgene_tss\tvariant_position\n')
	used_variant_gene_pairs = {}
	#for time_step in range(3,-1,-1):
	for time_step in range(4):
		counter = 0
		cell_type_sig_file = single_cell_eqtl_dir + 'sc_per_time_step_eqtl_analysis_' + str(time_step) + '_day_10_pcs_all_variant_gene_pairs_merged.txt' 
		f = open(cell_type_sig_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			test_name = data[0] + '_' + data[1]
			if data[0] in genes and float(data[7]) < nominal_p and test_name not in used_variant_gene_pairs:
				t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[3] + '\t' + data[4] + '\n')
				used_variant_gene_pairs[test_name] = 1
				counter = counter + 1
		print(str(time_step) + '\t' + str(counter))
		f.close()
	t.close()

def extract_sig_variant_gene_pairs_from_known_cell_types(single_cell_eqtl_dir, ld_pruned_variant_gene_pair_file, gene_file, random_seed):
	np.random.seed(random_seed)
	# First extract genes we have expression for
	genes = {}
	f = open(gene_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0].split('_')[0]
		genes[gene_id] = 1
	f.close()
	# Re-order cell types
	tests = {}
	used_genes = {}
	used_snps = {}
	t = open(ld_pruned_variant_gene_pair_file, 'w')
	t.write('Gene_id\tvariant_id\tchrom_num\tgene_tss\tvariant_position\n')
	# Add test keys to tests one cell type at a time
	for day in range(3,-1,-1):
		counter = 0
		cell_type_sig_file = single_cell_eqtl_dir + 'sc_per_time_step_eqtl_analysis_' + str(day) + '_day_10_pcs_all_variant_gene_pairs_multiple_testing_bf_bh_0.1_fdr_.txt'
		f = open(cell_type_sig_file)
		arr = []
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			arr.append(data)
		f.close()
		new_mat = np.asarray(arr)
		nrows = new_mat.shape[0]
		for row_num in range(nrows):
			gene_id = new_mat[row_num, 0]
			variant_id = new_mat[row_num, 1]
			test_name = gene_id + '_' + variant_id
			if test_name in tests:
				continue
			if gene_id in used_genes:
				continue
			if variant_id in used_snps:
				continue
			if gene_id not in genes:
				continue
			tests[test_name] = 1
			used_genes[gene_id] = 1
			used_snps[variant_id] = 1
			counter = counter + 1
			t.write('\t'.join(new_mat[row_num, :5]) + '\n')
	t.close()

def center_data(input_file, output_file):
	f = open(input_file)
	t = open(output_file, 'w')
	for line in f:
		line = line.rstrip()
		data = line.split()
		expr = np.asarray(data).astype(float)
		centered_expr = expr - np.mean(expr)
		t.write('\t'.join(centered_expr.astype(str)) + '\n')
	f.close()
	t.close()

def zero_center_expression_data(corrected_expression_file, raw_expression_file, output_file):
	t = open(output_file, 'w')
	f = open(corrected_expression_file)
	g = open(raw_expression_file)
	for line in f:
		line = line.rstrip()
		corrected_expr = np.asarray(line.split('\t')).astype(float)
		raw_expr = np.asarray(g.next().rstrip().split('\t'))
		zero_indices = (raw_expr.astype(float)) == 0.0
		if sum(zero_indices) == 0:
			expr_zeros_at_zero = corrected_expr - 0
		else:
			zero_mean = np.mean(corrected_expr[zero_indices])
			expr_zeros_at_zero = corrected_expr - zero_mean
		t.write('\t'.join(expr_zeros_at_zero.astype(str)) + '\n')
	t.close()
	f.close()
	f.close()

def subset_covariates(covariate_file, subset_covariate_file, num_cov):
	f = open(covariate_file)
	t = open(subset_covariate_file, 'w')
	for line in f:
		line = line.rstrip()
		data = line.split()
		t.write('\t'.join(data[:num_cov]) + '\n')
	f.close()
	t.close()



def prepare_eqtl_factorization_files_wrapper(output_root, gene_annotation_file, distance, genotype_data_dir, gene_file,raw_expression_file, expression_file, r_squared_threshold, max_variants_per_gene, num_pcs, covariate_file, cell_level_info_file, known_cell_type_file, pseudobulk_eqtl_dir, random_seed):
	'''
	########################
	# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
	########################
	# Output file containing list of variant gene pairs
	variant_gene_pair_file = output_root + '_variant_gene_pairs.txt'
	extract_variant_gene_pairs_for_eqtl_testing(gene_file, gene_annotation_file, distance, genotype_data_dir, variant_gene_pair_file)


	########################
	# Step 2: LD prune the above file in each gene (ie limit to only independent snps per gene)
	########################
	# Output file containing list of (pruned) variant gene pairs
	ld_pruned_variant_gene_pair_file = output_root + '_variant_gene_pairs_r_squared_pruned.txt'
	ld_prune_variant_gene_pair_file(variant_gene_pair_file, ld_pruned_variant_gene_pair_file, r_squared_threshold, genotype_data_dir, max_variants_per_gene, random_seed)
	'''
	ld_pruned_variant_gene_pair_file = output_root + '_variant_gene_pairs_in_known_cell_types.txt'
	num_tests_per_cell_type = 800
	extract_sig_variant_gene_pairs_from_known_cell_types(known_cell_type_file, single_cell_eqtl_dir, num_tests_per_cell_type, ld_pruned_variant_gene_pair_file, gene_file, random_seed)


	########################
	# Step 3: Generate raw expression matrix
	########################
	# Output file
	single_raw_cell_expression_eqtl_traing_data_file = output_root + '_raw_expression_training_data_uncorrected_r_squared_pruned.txt'
	generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, raw_expression_file, gene_file, single_raw_cell_expression_eqtl_traing_data_file)
	save_as_h5_file(single_raw_cell_expression_eqtl_traing_data_file)

	########################
	# Step 3: Generate expression matrix
	########################
	# Output file
	single_cell_expression_eqtl_traing_data_file = output_root + '_expression_training_data_uncorrected_r_squared_pruned.txt'
	generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, expression_file, gene_file, single_cell_expression_eqtl_traing_data_file)

	########################
	# Step 3.5: Generate centered expression matrix
	########################
	single_cell_centered_expression_eqtl_traing_data_file = output_root + '_expression_training_data_centered_uncorrected_r_squared_pruned.txt'
	center_data(single_cell_expression_eqtl_traing_data_file, single_cell_centered_expression_eqtl_traing_data_file)
	save_as_h5_file(single_cell_centered_expression_eqtl_traing_data_file)

	########################
	# Step 3b: Get subset of covariates
	########################
	num_cov = 10
	subset_covariate_file = output_root + '_covariate_subset_' + str(num_cov) + '.txt'
	subset_covariates(covariate_file, subset_covariate_file, num_cov)


	########################
	# Step 3b: center zeros at zero
	########################
	# Residual expression file
	single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file = output_root + '_expression_training_data_uncorrected_zero_centered_r_squared_pruned.txt'
	zero_center_expression_data(single_cell_expression_eqtl_traing_data_file, single_raw_cell_expression_eqtl_traing_data_file, single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file)
	save_as_h5_file(single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file)

	########################
	# Step 4: Generate residual expression
	########################
	# Residual expression file
	single_cell_corrected_expression_eqtl_traing_data_file = output_root + '_expression_training_data_corrected_r_squared_pruned.txt'
	regress_out_covariates(single_cell_expression_eqtl_traing_data_file, covariate_file, single_cell_corrected_expression_eqtl_traing_data_file, num_pcs)
	save_as_h5_file(single_cell_corrected_expression_eqtl_traing_data_file)
	

	########################
	# Step 4b: center zeros at zero
	########################
	# Residual expression file
	single_cell_corrected_zero_centered_expression_eqtl_traing_data_file = output_root + '_expression_training_data_corrected_zero_centered_r_squared_pruned.txt'
	zero_center_expression_data(single_cell_corrected_expression_eqtl_traing_data_file, single_raw_cell_expression_eqtl_traing_data_file, single_cell_corrected_zero_centered_expression_eqtl_traing_data_file)
	save_as_h5_file(single_cell_corrected_zero_centered_expression_eqtl_traing_data_file)

	########################
	# Step 5: Generate Genotype matrix
	########################
	# Output file
	single_cell_genotype_eqtl_training_data_file = output_root + '_genotype_training_data_uncorrected_r_squared_pruned.txt'
	construct_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_data_dir, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)
	save_as_h5_file(single_cell_genotype_eqtl_training_data_file)

	########################
	# Step 6: Generate residual genotype
	########################
	# corrected expression file (output file)
	single_cell_corrected_genotype_eqtl_traing_data_file = output_root + '_genotype_training_data_corrected_r_squared_pruned.txt'
	#regress_out_covariates(single_cell_genotype_eqtl_training_data_file, covariate_file, single_cell_corrected_genotype_eqtl_traing_data_file, num_pcs)
	#save_as_h5_file(single_cell_corrected_genotype_eqtl_traing_data_file)

	########################
	# Step 5: Generate Genotype matrix
	########################
	# Output file
	single_cell_genotype_eqtl_training_data_file = output_root + '_standardized_genotype_training_data_uncorrected_r_squared_pruned.txt'
	construct_standardized_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_data_dir, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)
	save_as_h5_file(single_cell_genotype_eqtl_training_data_file)

	########################
	# Step 6: Generate residual genotype
	########################
	# corrected expression file (output file)
	single_cell_corrected_genotype_eqtl_traing_data_file = output_root + '_standardized_genotype_training_data_corrected_r_squared_pruned.txt'
	regress_out_covariates(single_cell_genotype_eqtl_training_data_file, covariate_file, single_cell_corrected_genotype_eqtl_traing_data_file, num_pcs)
	save_as_h5_file(single_cell_corrected_genotype_eqtl_traing_data_file)

	########################
	# Step 7: Generate individual id file (z matrix in eqtl factorization)
	########################
	# Output file
	single_cell_individual_id_file = output_root + '_individual_id.txt'
	generate_individual_id_file(cell_level_info_file, single_cell_individual_id_file)




def prepare_eqtl_factorization_files_wrapper_v2(output_root, gene_annotation_file, distance, genotype_file, raw_expression_file, expression_file, nominal_p, r_squared_threshold, num_pcs, covariate_file, cell_level_info_file, single_cell_eqtl_dir, random_seed):
	########################
	# Step 1: Extract variant gene pairs that reached nominal significance in cell type specific eqtl analyis (ignore hla genes)
	########################
	variant_gene_pair_file = output_root + '_nominal_sig_variant_gene_pairs_in_known_time_steps.txt'
	extract_nominal_sig_variant_gene_pairs_from_known_cell_types(single_cell_eqtl_dir, nominal_p, variant_gene_pair_file, expression_file)

	########################
	# Step 2: LD prune the above file in each gene (ie limit to only independent snps per gene)
	########################
	# Output file containing list of (pruned) variant gene pairs
	ld_pruned_variant_gene_pair_file = output_root + '_nominal_sig_variant_gene_pairs_in_known_cell_types_r_squared_pruned.txt'
	ld_prune_variant_gene_pair_file(variant_gene_pair_file, ld_pruned_variant_gene_pair_file, r_squared_threshold, genotype_file, random_seed)

	########################
	# Step 3: Generate raw expression matrix
	########################
	# Output file
	single_raw_cell_expression_eqtl_traing_data_file = output_root + '_raw_expression_training_data_uncorrected_r_squared_pruned.txt'
	generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, raw_expression_file, single_raw_cell_expression_eqtl_traing_data_file)
	save_as_h5_file(single_raw_cell_expression_eqtl_traing_data_file)


	########################
	# Step 4: Generate normalized expression matrix
	########################
	# Output file
	single_cell_expression_eqtl_traing_data_file = output_root + '_expression_training_data_uncorrected_r_squared_pruned.txt'
	generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, expression_file, single_cell_expression_eqtl_traing_data_file)

	
	########################
	# Step 5: center zeros at zero
	########################
	# Residual expression file
	single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file = output_root + '_expression_training_data_uncorrected_zero_centered_r_squared_pruned.txt'
	zero_center_expression_data(single_cell_expression_eqtl_traing_data_file, single_raw_cell_expression_eqtl_traing_data_file, single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file)
	save_as_h5_file(single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file)

	
	########################
	# Step 6: Get subset of covariates
	########################
	subset_covariate_file = output_root + '_covariate_subset_' + str(num_pcs) + '.txt'
	subset_covariates(covariate_file, subset_covariate_file, num_pcs)	

	########################
	# Step 7: Generate Genotype matrix
	########################
	# Output file
	single_cell_genotype_eqtl_training_data_file = output_root + '_genotype_training_data_uncorrected_r_squared_pruned.txt'
	construct_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)
	save_as_h5_file(single_cell_genotype_eqtl_training_data_file)

	########################
	# Step 8: Generate standardized Genotype matrix
	########################
	# Output file
	single_cell_genotype_eqtl_training_data_file = output_root + '_standardized_genotype_training_data_uncorrected_r_squared_pruned.txt'
	construct_standardized_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)
	save_as_h5_file(single_cell_genotype_eqtl_training_data_file)

	########################
	# Step 9: Generate individual id file (z matrix in eqtl factorization)
	########################
	# Output file
	single_cell_individual_id_file = output_root + '_individual_id.txt'
	generate_individual_id_file(cell_level_info_file, single_cell_individual_id_file)

def prepare_eqtl_factorization_files_wrapper_v3(output_root, gene_annotation_file, distance, genotype_file, raw_expression_file, expression_file, num_pcs, covariate_file, cell_level_info_file, single_cell_eqtl_dir, random_seed):
	########################
	# Step 1 and step 2: Extract variant gene pairs that reached nominal significance in cell type specific eqtl analyis (ignore hla genes)
	########################
	ld_pruned_variant_gene_pair_file = output_root + '_variant_gene_pairs_in_known_days.txt'
	extract_sig_variant_gene_pairs_from_known_cell_types(single_cell_eqtl_dir, ld_pruned_variant_gene_pair_file, expression_file, random_seed)

	########################
	# Step 3: Generate raw expression matrix
	########################
	# Output file
	single_raw_cell_expression_eqtl_traing_data_file = output_root + '_raw_expression_training_data_uncorrected_r_squared_pruned.txt'
	generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, raw_expression_file, single_raw_cell_expression_eqtl_traing_data_file)
	save_as_h5_file(single_raw_cell_expression_eqtl_traing_data_file)


	########################
	# Step 4: Generate normalized expression matrix
	########################
	# Output file
	single_cell_expression_eqtl_traing_data_file = output_root + '_expression_training_data_uncorrected_r_squared_pruned.txt'
	generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, expression_file, single_cell_expression_eqtl_traing_data_file)

	
	########################
	# Step 5: center zeros at zero
	########################
	# Residual expression file
	single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file = output_root + '_expression_training_data_uncorrected_zero_centered_r_squared_pruned.txt'
	zero_center_expression_data(single_cell_expression_eqtl_traing_data_file, single_raw_cell_expression_eqtl_traing_data_file, single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file)
	save_as_h5_file(single_cell_uncorrected_zero_centered_expression_eqtl_traing_data_file)

	
	########################
	# Step 6: Get subset of covariates
	########################
	subset_covariate_file = output_root + '_covariate_subset_' + str(num_pcs) + '.txt'
	subset_covariates(covariate_file, subset_covariate_file, num_pcs)	

	########################
	# Step 7: Generate Genotype matrix
	########################
	# Output file
	single_cell_genotype_eqtl_training_data_file = output_root + '_genotype_training_data_uncorrected_r_squared_pruned.txt'
	construct_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)
	save_as_h5_file(single_cell_genotype_eqtl_training_data_file)

	########################
	# Step 8: Generate standardized Genotype matrix
	########################
	# Output file
	single_cell_genotype_eqtl_training_data_file = output_root + '_standardized_genotype_training_data_uncorrected_r_squared_pruned.txt'
	construct_standardized_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)
	save_as_h5_file(single_cell_genotype_eqtl_training_data_file)

	########################
	# Step 9: Generate individual id file (z matrix in eqtl factorization)
	########################
	# Output file
	single_cell_individual_id_file = output_root + '_individual_id.txt'
	generate_individual_id_file(cell_level_info_file, single_cell_individual_id_file)
	
######################
# Command line args
######################
gene_annotation_file = sys.argv[1]
pre_processed_data_dir = sys.argv[2]
eqtl_factorization_input_dir = sys.argv[3]
single_cell_eqtl_dir = sys.argv[4]


################
# single cell randm subset
#############
# Variant must be within $distance BP from TSS of gene
distance = 25000

# Input file containing expression
raw_expression_file = pre_processed_data_dir + 'normalized_expression_all_cells.txt'
# Input file containing expression
expression_file = pre_processed_data_dir + 'standardized_normalized_expression_all_cells.txt'
# Covariate file
covariate_file = pre_processed_data_dir + 'standardized_normalized_expression_pca_loadings.txt'
# File containing mapping from cell index to individual id
cell_level_info_file = pre_processed_data_dir + 'cell_covariates.txt'
# Genotype file
genotype_file = pre_processed_data_dir + 'genotype_mean_inputed.txt'
# Nominal p-value threshold
nominal_p=.001
# Only allow snps with r_squared threshold less than this
r_squared_threshold=0.2
# Number of PCs to use
num_pcs = 10
# random seed used for ld pruning
random_seed=1
# Output root
output_root = eqtl_factorization_input_dir + 'single_cell_sig_tests_nominal_p_' + str(nominal_p) + '_r_squared_thresh_' + str(r_squared_threshold)


#prepare_eqtl_factorization_files_wrapper_v2(output_root, gene_annotation_file, distance, genotype_file, raw_expression_file, expression_file, nominal_p, r_squared_threshold, num_pcs, covariate_file, cell_level_info_file, single_cell_eqtl_dir, random_seed)



################
# single cell randm subset
#############
# Variant must be within $distance BP from TSS of gene

# Input file containing expression
raw_expression_file = pre_processed_data_dir + 'normalized_expression_all_cells.txt'
# Input file containing expression
expression_file = pre_processed_data_dir + 'standardized_normalized_expression_all_cells.txt'
# Covariate file
covariate_file = pre_processed_data_dir + 'standardized_normalized_expression_pca_loadings.txt'
# File containing mapping from cell index to individual id
cell_level_info_file = pre_processed_data_dir + 'cell_covariates.txt'
# Genotype file
genotype_file = pre_processed_data_dir + 'genotype_mean_inputed.txt'
# Number of PCs to use
num_pcs = 10
# random seed used for ld pruning
random_seed=1
# Output root
output_root = eqtl_factorization_input_dir + 'single_cell_sig_tests_in_each_day'


prepare_eqtl_factorization_files_wrapper_v3(output_root, gene_annotation_file, distance, genotype_file, raw_expression_file, expression_file, num_pcs, covariate_file, cell_level_info_file, single_cell_eqtl_dir, random_seed)








