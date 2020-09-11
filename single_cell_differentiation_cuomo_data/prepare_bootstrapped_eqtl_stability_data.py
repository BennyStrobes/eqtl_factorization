import numpy as np 
import os
import sys
import pdb


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
		gene_id = data[0].split('.')[0]
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
		if data[1] == 'chrX' or data[1] == 'chrY' or data[1] == 'chrM':
			continue
		# Extract relevent info on gene: Chrom num and TSS
		chrom_num = int(data[1].split('hr')[1])
		strand = data[4]
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



########################
# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
########################
def extract_variant_gene_pairs_for_eqtl_testing(gene_annotation_file, expression_file, distance, genotype_file, variant_gene_pair_file):
	# Extract gene list
	genes = []
	f = open(expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0].split('_')[0]
		genes.append(gene_name)
	f.close()
	genes = np.asarray(genes)
	# Get mapping from genes to (chrom_num, position)
	gene_mapping = get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes)
	# Open file handle to output file containing variant-gene pair tests and print header
	t = open(variant_gene_pair_file, 'w')
	t.write('Gene_id\tvariant_id\tchrom_num\tgene_tss\tvariant_position\tmaf\n')

	# Fill in file containing lists of variant gene pairs for each chromosome iteratively
	used_genes = {}
	for chrom_num in range(1,23):
		print(chrom_num)
		# Create array where each element is a BP in this chromosome
		# 'Null' if no genes in distance BP of gene
		# Otherwise is a list of gene names
		chromosome = create_gene_chromsome(chrom_num, gene_mapping, distance)
		# Now loop through variants on this chromosome
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
			variant_chrom = int(variant_id.split(':')[0])
			variant_pos = int(variant_id.split(':')[1])
			# Simple error check
			if variant_chrom != chrom_num:
				continue
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
				t.write(gene_id + '\t' + variant_id + '\t' + str(chrom_num) + '\t' + str(gene_mapping[gene_id][1]) + '\t' + str(variant_pos) + '\t' + str(maf) + '\n')
		f.close()
	t.close()

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
	return mapping

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

def subset_covariates(covariate_file, subset_covariate_file, num_cov):
	f = open(covariate_file)
	t = open(subset_covariate_file, 'w')
	for line in f:
		line = line.rstrip()
		data = line.split()
		t.write('\t'.join(data[:num_cov]) + '\n')
	f.close()
	t.close()

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

def get_cell_level_ordered_experimental_array(cell_level_info_file):
	array = []
	head_count = 0
	f = open(cell_level_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		array.append(data[9])
	f.close()
	return np.asarray(array)

def get_cell_level_ordered_plate_array(cell_level_info_file):
	array = []
	head_count = 0
	f = open(cell_level_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		array.append(data[59])
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
		cell_level_genotype_vector = cell_level_genotype_vector.astype(float)
		#standardized_cell_level_genotype_vector = (cell_level_genotype_vector - np.mean(cell_level_genotype_vector))/np.std(cell_level_genotype_vector)
		t.write('\t'.join(cell_level_genotype_vector.astype(str)) + '\n')
	t.close()

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

# Step 6: Generate individual id file (z matrix in eqtl factorization)
def generate_experiment_id_file(cell_level_info_file, single_cell_individual_id_file):
	# Get array of length number of cells and each element corresponds to the individual id corresponding to that cell
	indi_ids = get_cell_level_ordered_experimental_array(cell_level_info_file)
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


# Step 6: Generate individual id file (z matrix in eqtl factorization)
def generate_plate_id_file(cell_level_info_file, single_cell_individual_id_file):
	# Get array of length number of cells and each element corresponds to the individual id corresponding to that cell
	indi_ids = get_cell_level_ordered_plate_array(cell_level_info_file)
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

# Step 6: Generate individual id file (z matrix in eqtl factorization)
def generate_environmental_variable_file(cell_level_info_file, single_cell_environment_file):
	# Get array of length number of cells and each element corresponds to the individual id corresponding to that cell
	t = open(single_cell_environment_file, 'w')
	f = open(cell_level_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		t.write(data[93] + '\n')
	f.close()
	t.close()

def prepare_bootstrapped_eqtl_stability_input_data(gene_annotation_file, expression_file, covariate_file, cell_level_info_file, genotype_file, num_pcs, distance, random_seed, output_root):
	np.random.seed(random_seed)
	########################
	# Step 1: Generate input variant gene pairs
	########################
	variant_gene_pair_file = output_root + 'variant_gene_pairs.txt'
	#extract_variant_gene_pairs_for_eqtl_testing(gene_annotation_file, expression_file, distance, genotype_file, variant_gene_pair_file)

	########################
	# Step 2: Generate normalized expression matrix
	########################
	# Output file
	single_cell_expression_eqtl_traing_data_file = output_root + 'expression.txt'
	#generate_single_cell_expression_eqtl_training_data(variant_gene_pair_file, expression_file, single_cell_expression_eqtl_traing_data_file)

	########################
	# Step 3: Get subset of covariates
	########################
	subset_covariate_file = output_root + 'covariate_subset_' + str(num_pcs) + '.txt'
	#subset_covariates(covariate_file, subset_covariate_file, num_pcs)	

	########################
	# Step 4: Generate standardized Genotype matrix
	########################
	# Output file
	single_cell_genotype_eqtl_training_data_file = output_root + 'standardized_genotype.txt'
	#construct_standardized_genotype_matrix(variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)

	########################
	# Step 4: Generate Genotype matrix
	########################
	# Output file
	single_cell_genotype_eqtl_training_data_file = output_root + 'genotype.txt'
	#construct_genotype_matrix(variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)

	########################
	# Step 5: Generate individual id file (z matrix in eqtl factorization)
	########################
	# Output file
	single_cell_individual_id_file = output_root + 'individual_id.txt'
	#generate_individual_id_file(cell_level_info_file, single_cell_individual_id_file)

	########################
	# Step 5: Generate experiment id (z matrix in eqtl factorization)
	########################
	# Output file
	single_cell_experiment_id_file = output_root + 'experiment_id.txt'
	#generate_experiment_id_file(cell_level_info_file, single_cell_experiment_id_file)

	########################
	# Step 5: Generate plate id (z matrix in eqtl factorization)
	########################
	# Output file
	single_cell_plate_id_file = output_root + 'plate_id.txt'
	#generate_plate_id_file(cell_level_info_file, single_cell_plate_id_file)

	########################
	# Step 5: Generate environmental variable file 
	########################
	# Output file
	single_cell_environment_file = output_root + 'pseudotime.txt'
	generate_environmental_variable_file(cell_level_info_file, single_cell_environment_file)


gene_annotation_file = sys.argv[1]
pre_processed_data_dir = sys.argv[2]
bootstrapped_eqtl_stability_dir = sys.argv[3]




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
# Distance
distance=250000
# Output root
output_root = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_input_data_cis_window_' + str(distance) + '_'

prepare_bootstrapped_eqtl_stability_input_data(gene_annotation_file, expression_file, covariate_file, cell_level_info_file, genotype_file, num_pcs, distance, random_seed, output_root)
