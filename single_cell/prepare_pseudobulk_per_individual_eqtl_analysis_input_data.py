import numpy as np 
import os
import sys
import pdb

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


########################
# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
########################
def extract_variant_gene_pairs_for_eqtl_testing(gene_file, gene_annotation_file, distance, genotype_data_dir, variant_gene_pair_file):
	# Extract gene list
	genes = np.loadtxt(gene_file, delimiter='\t',dtype=str)[:,0]
	# Get mapping from genes to (chrom_num, position)
	gene_mapping = get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes)
	# Open file handle to output file containing variant-gene pair tests and print header
	t = open(variant_gene_pair_file, 'w')
	t.write('Gene_id\tvariant_id\tchrom_num\tgene_tss\tvariant_position\n')

	# Fill in file containing lists of variant gene pairs for each chromosome iteratively
	used_genes = {}
	for chrom_num in range(1,23):
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
	gene_names = np.loadtxt(gene_names_file,dtype=str, delimiter='\t')[:,0]
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
		array.append(data[2])
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
	return np.asarray(indi)

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


def construct_covariate_matrix(pseudobulk_covariate_file, eqtl_covariate_file, num_pcs):
	pcs = np.loadtxt(pseudobulk_covariate_file)
	np.savetxt(eqtl_covariate_file, pcs[:, :num_pcs], fmt="%s", delimiter='\t')

def generate_eqtl_input_files(genotype_data_dir, pseudobulk_expression_file, pseudobulk_covariate_file, sample_covariate_file, eqtl_variant_gene_pairs_file, eqtl_expression_file, eqtl_genotype_file0, eqtl_covariate_file, distance, gene_annotation_file, gene_id_file, num_pcs):
	########################
	# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
	########################
	extract_variant_gene_pairs_for_eqtl_testing(gene_id_file, gene_annotation_file, distance, genotype_data_dir, eqtl_variant_gene_pairs_file)

	########################
	# Step 2: Generate expression matrix
	########################
	generate_single_cell_expression_eqtl_training_data(eqtl_variant_gene_pairs_file, pseudobulk_expression_file, gene_id_file, eqtl_expression_file)

	########################
	# Step 3: Generate Genotype matrix
	#########################
	construct_genotype_matrix(eqtl_variant_gene_pairs_file, genotype_data_dir, sample_covariate_file, eqtl_genotype_file)

	########################
	# Step 4: Generate Covariate matrix
	#########################
	construct_covariate_matrix(pseudobulk_covariate_file, eqtl_covariate_file, num_pcs)


###################
# Command Line args
###################
processed_expression_dir = sys.argv[1]  # Input dir
gene_annotation_file = sys.argv[2]
genotype_data_dir = sys.argv[3]
pseudobulk_eqtl_dir = sys.argv[4]  # Output dir



###################
# Input files
###################
# Raw pseudobulk expression file
pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_per_individual_expression_sle_individuals_standardized.txt'
# Pseudobulk sample covariate file (in same order as raw_pseudobulk_expression_file)
pseudobulk_covariate_file = processed_expression_dir + 'pca_scores_pseudobulk_per_individual_sle_individuals.txt'
# Pseudobulk gene names (in same order as raw_pseudobulk_expression_file)
gene_id_file = processed_expression_dir + 'single_cell_expression_sle_individuals_min_expressed_cells_0.05_log_transform_transform_regress_out_batch_True_gene_ids.txt'
# sample covariate file
sample_covariate_file = processed_expression_dir + 'pseudobulk_per_individual_covariates_sle_individuals.txt'




###################
# Generate eqtl input files
###################
distance=10000
num_pcs=15
# Output files
eqtl_variant_gene_pairs_file = pseudobulk_eqtl_dir  + 'pseudobulk_per_individual_eqtl_input_variant_gene_pairs.txt'
eqtl_expression_file = pseudobulk_eqtl_dir + 'pseudobulk_per_individual_eqtl_input_expression.txt'
eqtl_raw_expression_file = pseudobulk_eqtl_dir + 'pseudobulk_per_individual_eqtl_input_raw_expression.txt'
eqtl_genotype_file = pseudobulk_eqtl_dir + 'pseudobulk_per_individual_eqtl_input_genotype.txt'
eqtl_covariate_file = pseudobulk_eqtl_dir + 'pseudobulk_per_individual_pca_scores.txt'
generate_eqtl_input_files(genotype_data_dir, pseudobulk_expression_file, pseudobulk_covariate_file, sample_covariate_file, eqtl_variant_gene_pairs_file, eqtl_expression_file, eqtl_genotype_file, eqtl_covariate_file, distance, gene_annotation_file, gene_id_file, num_pcs)
