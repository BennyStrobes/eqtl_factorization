import numpy as np 
import os
import sys
import pdb
from sklearn import linear_model





# Get mapping from genes to (chrom_num, position)
def get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes):
	# Convert gene array to dictionary
	gene_mapping = {}
	for gene in genes:
		gene_mapping[gene] = [0,0]
	# Fill in gene positions
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
		if gene_id not in gene_mapping:
			continue
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
	genes = np.loadtxt(gene_file, delimiter='\n',dtype=str)
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
	print(len(used_genes))

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
def create_mapping_from_variants_to_genotype(variant_list, genotype_data_dir):
	variants = {}
	for chrom_num in range(1,23):
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
			# Simple error checking
			if len(data) != 120:
				print('assumption error!')
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
def ld_prune_variant_gene_pair_file(variant_gene_pair_file, ld_pruned_variant_gene_pair_file, r_squared_threshold, genotype_data_dir, random_seed):
	np.random.seed(random_seed)
	# First extract list of variants
	variant_list = extract_variants_from_variant_gene_pair_file(variant_gene_pair_file)
	# Second extract mapping from genes to array of variants
	gene_to_variant_mapping_unpruned = create_mapping_from_gene_to_variants_in_variant_gene_pair_file(variant_gene_pair_file)
	# Create mapping from variants in variat_list to genotype vectors
	variants = create_mapping_from_variants_to_genotype(variant_list, genotype_data_dir)
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

def create_mapping_from_gene_names_to_expression_vectors(sc_expression_file, gene_names_file):
	# Get gene names
	gene_names = np.loadtxt(gene_names_file,dtype=str)
	# Load in expression matrix (every column corresponds to a gene)
	expression_matrix = np.loadtxt(sc_expression_file, dtype=str, delimiter='\t')
	# Create mapping
	mapping = {}
	for index, gene_name in enumerate(gene_names):
		mapping[gene_name] = expression_matrix[:, index]
	return mapping

def generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, sc_expression_file, gene_names_file, single_cell_expression_eqtl_traing_data_file):
	# Get ordered list of gene names (this will be the order that the output file will be saved in)
	ordered_gene_names = get_ordered_list_of_gene_names(ld_pruned_variant_gene_pair_file)
	# Create mapping from gene names to expression vectors
	gene_name_to_expression_vector = create_mapping_from_gene_names_to_expression_vectors(sc_expression_file, gene_names_file)
	# print to output file
	t = open(single_cell_expression_eqtl_traing_data_file, 'w')
	for gene_name in ordered_gene_names:
		# Use map to get expression vector corresponding to this gene
		expression_vector = gene_name_to_expression_vector[gene_name]
		# Print to output file
		t.write('\t'.join(expression_vector) + '\n')
	t.close()


def regress_out_covariates(expression_input_file, covariate_file, expression_output_file, num_pcs):
	# Load in covariates
	covariate_raw = np.loadtxt(covariate_file, dtype=str, delimiter='\t')

	covariate_mat = covariate_raw[:,:num_pcs].astype(float)
	

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
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		array.append(data[4])
	f.close()
	return np.asarray(array)

def get_genotype_level_ordered_individual_array(genotype_data_dir):
	chrom_num = 1
	genotype_file = genotype_data_dir + 'chr' + str(chrom_num) + '.genotypes.matrix.eqtl.txt'
	f = open(genotype_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			line = line.rstrip()
			data = line.split()
			indi = data
			head_count = head_count + 1
			continue
	f.close()
	'''
	# Errror checking making sure other chromosmes have same order
	for chrom_num in range(1,23):
		genotype_file = genotype_data_dir + 'chr' + str(chrom_num) + '.genotypes.matrix.eqtl.txt'
		f = open(genotype_file)
		head_count = 0
		for line in f:
			if head_count == 0:
				head_count = head_count + 1
				line = line.rstrip()
				data = line.split()
				if np.array_equal(np.asarray(data), np.asarray(indi)) == False:
					print('assumption erro')
					pdb.set_trace()
				continue
		f.close()
	'''
	return np.asarray(indi)


def construct_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_data_dir, cell_level_info_file, single_cell_genotype_eqtl_training_data_file):
	variant_names, variant_list = get_ordered_list_of_variant_names(ld_pruned_variant_gene_pair_file)

	ordered_individuals_cell_level = get_cell_level_ordered_individaul_array(cell_level_info_file)
	
	ordered_individuals_genotype_level = get_genotype_level_ordered_individual_array(genotype_data_dir)

	# Create mapping from variant_id to genotype vectors
	variant_to_genotype = create_mapping_from_variants_to_genotype(variant_list, genotype_data_dir)
	
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


######################
# Command line args
######################
gene_annotation_file = sys.argv[1]
processed_expression_dir = sys.argv[2]
eqtl_input_dir = sys.argv[3]
genotype_data_dir = sys.argv[4]



########################
# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
########################
# Variant must be within $distance BP from TSS of gene
distance = 10000
# File containing availible genes
gene_file = processed_expression_dir + 'gene_names_ye_lab.txt'
# Output file containing list of variant gene pairs
variant_gene_pair_file = eqtl_input_dir + 'variant_gene_pairs_' + str(distance) + '_bp.txt'

#extract_variant_gene_pairs_for_eqtl_testing(gene_file, gene_annotation_file, distance, genotype_data_dir, variant_gene_pair_file)


########################
# Step 2: LD prune the above file in each gene (ie limit to only independent snps per gene)
########################
# random seed used for ld pruning
random_seed=1
# Only allow snps with r_squared threshold less than this
r_squared_threshold=0.7
# Input file containing list of (un-pruned) variant gene pairs
variant_gene_pair_file = eqtl_input_dir + 'variant_gene_pairs_' + str(distance) + '_bp.txt'
# Output file containing list of (pruned) variant gene pairs
ld_pruned_variant_gene_pair_file = eqtl_input_dir + 'variant_gene_pairs_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'

#ld_prune_variant_gene_pair_file(variant_gene_pair_file, ld_pruned_variant_gene_pair_file, r_squared_threshold, genotype_data_dir, random_seed)

# Only allow snps with r_squared threshold less than this
r_squared_threshold=0.5
# Input file containing list of (un-pruned) variant gene pairs
variant_gene_pair_file = eqtl_input_dir + 'variant_gene_pairs_' + str(distance) + '_bp.txt'
# Output file containing list of (pruned) variant gene pairs
ld_pruned_variant_gene_pair_file = eqtl_input_dir + 'variant_gene_pairs_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'

#ld_prune_variant_gene_pair_file(variant_gene_pair_file, ld_pruned_variant_gene_pair_file, r_squared_threshold, genotype_data_dir, random_seed)




########################
# Step 3: Generate expression matrix
########################
distance = 10000
r_squared_threshold=0.5
# Output file
single_cell_expression_eqtl_traing_data_file = eqtl_input_dir + 'sc_expression_training_data_uncorrected_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'
# Input file 
sc_expression_file = processed_expression_dir + 'single_cell_expression_sle_individuals_standardized.txt'
# Input gene list
gene_name_file = processed_expression_dir + 'gene_names_ye_lab.txt'
# Input test names
ld_pruned_variant_gene_pair_file = eqtl_input_dir + 'variant_gene_pairs_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'

#generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, sc_expression_file, gene_name_file, single_cell_expression_eqtl_traing_data_file)


########################
# Step 3: Generate residual expression
########################
# num_pcs
num_pcs = 50
# Input file
single_cell_expression_eqtl_traing_data_file = eqtl_input_dir + 'sc_expression_training_data_uncorrected_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'
# Output file
single_cell_corrected_expression_eqtl_traing_data_file = eqtl_input_dir + 'sc_corrected_expression_training_data_uncorrected_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'
# Covariate file
covariate_file = processed_expression_dir + 'pca_scores_sle_individuals.txt'

regress_out_covariates(single_cell_expression_eqtl_traing_data_file, covariate_file, single_cell_corrected_expression_eqtl_traing_data_file, num_pcs)


########################
# Step 4: Generate Genotype matrix
########################
# File containing mapping from cell index to individual id
cell_level_info_file = processed_expression_dir + 'cell_covariates_sle_individuals.txt'
# Output file
single_cell_genotype_eqtl_training_data_file = eqtl_input_dir + 'sc_genotype_training_data_uncorrected_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'

# construct_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_data_dir, cell_level_info_file, single_cell_genotype_eqtl_training_data_file)


########################
# Step 5: Generate residual genotype
########################
# num_pcs
num_pcs = 50
# Input file
single_cell_genotype_eqtl_training_data_file = eqtl_input_dir + 'sc_genotype_training_data_uncorrected_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'
# Output file
single_cell_corrected_genotype_eqtl_traing_data_file = eqtl_input_dir + 'sc_corrected_genotype_training_data_uncorrected_' + str(distance) + '_bp_' + str(r_squared_threshold) + '_r_squared_pruned.txt'
# Covariate file
covariate_file = processed_expression_dir + 'pca_scores_sle_individuals.txt'

regress_out_covariates(single_cell_genotype_eqtl_training_data_file, covariate_file, single_cell_corrected_genotype_eqtl_traing_data_file, num_pcs)











# STEP 6 Make z matrix
