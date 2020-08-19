import numpy as np 
import os
import sys
import pdb

def create_gene_names_file(expression_file, gene_names_file):
	f = open(expression_file)
	t = open(gene_names_file, 'w')
	t.write('ensamble_id\tgene_id\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_full = data[0]
		ensamble_id = gene_full.split('_')[0]
		gene_id = gene_full.split('_')[1]
		t.write(ensamble_id + '\t' + gene_id + '\n')
	t.close()
	f.close()

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
		ensamble_id = data[0].split('.')[0]
		if ensamble_id not in gene_mapping:
			continue
		chrom_num = data[1]
		if chrom_num == 'chrX' or chrom_num == 'chrY' or chrom_num == 'chrM':
			continue
		# Simple error checking
		start = int(data[2])
		end = int(data[3])
		if start > end:
			print('assumption eroror')
			pdb.set_trace()
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
		if gene_mapping[ensamble_id][0] != 0:
			print('miss')
			pdb.set_trace()
			if gene_mapping[ensamble_id][0] != chrom_num and gene_mapping[ensamble_id][1] != tss:
				gene_mapping.pop(ensamble_id)
		else:
			gene_mapping[ensamble_id] = (chrom_num, tss)
	f.close()
	return gene_mapping

# Create array where each element is a BP in this chromosome
# 'Null' if no genes in distance BP of gene
# Otherwise is a list of gene names
def create_gene_chromsome(chrom_num, gene_mapping, distance):
	# Initialize chromosome
	chromosome = ['Null']*400000000
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
def extract_variant_gene_pairs_for_eqtl_testing(gene_file, gene_annotation_file, distance, genotype_file, variant_gene_pair_file):
	# Extract gene list
	genes = np.loadtxt(gene_file, delimiter='\t',dtype=str)[:,0]
	# Get mapping from genes to (chrom_num, position)
	gene_mapping = get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes)
	# Open file handle to output file containing variant-gene pair tests and print header
	t = open(variant_gene_pair_file, 'w')
	t.write('Gene_id\tvariant_id\tchrom_num\tgene_tss\tvariant_position\n')
	# Fill in file containing lists of variant gene pairs for each chromosome iteratively
	for chrom_num in range(1,23):
		print(chrom_num)
		# Create array where each element is a BP in this chromosome
		# 'Null' if no genes in distance BP of gene
		# Otherwise is a list of gene names
		chromosome = create_gene_chromsome(chrom_num, gene_mapping, distance)
		# Now loop through variants on this chromosome
		f = open(genotype_file)
		skipped = 0
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			if len(data) != 105:
				print('assumption error!')
				pdb.set_trace()
			variant_id = data[0]
			string_pos = variant_id.split('--')[1]
			if len(string_pos.split('-')) > 1:
				continue
			variant_chrom = int(variant_id.split('--')[0])
			variant_pos = int(variant_id.split('--')[1])
			# Ignore variants from other chromosomes
			if variant_chrom != chrom_num:
				continue
			# No genes within 10KB of variant
			if variant_pos > 400000000:
				skipped = skipped + 1
				pdb.set_trace()
				continue
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
				t.write(gene_id + '\t' + variant_id + '\t' + str(chrom_num) + '\t' + str(gene_mapping[gene_id][1]) + '\t' + str(variant_pos) + '\n')
		f.close()
		print(skipped)
	t.close()
	print(skipped)


def generate_cell_type_eqtl_input_files(day, genotype_file, gene_names_file, cell_type_sc_expression_file, cell_type_sc_sample_covariate_file, cell_type_eqtl_variant_gene_pairs_file, cell_type_eqtl_expression_file, cell_type_eqtl_genotype_file, distance, gene_annotation_file, cell_type_sample_overlap_file):
	########################
	# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
	########################
	extract_variant_gene_pairs_for_eqtl_testing(gene_names_file, gene_annotation_file, distance, genotype_file, cell_type_eqtl_variant_gene_pairs_file)



#####################
# Command line args
######################
pre_processed_data_dir = sys.argv[1]
gene_annotation_file = sys.argv[2]
per_time_step_eqtl_input_data_dir = sys.argv[3]


################
gene_names_file = per_time_step_eqtl_input_data_dir + 'gene_names.txt'
expression_file = pre_processed_data_dir + 'standardized_normalized_expression_all_cells.txt'
#create_gene_names_file(expression_file, gene_names_file)




###################
# For each day generate eqtl input files
###################
genotype_file = pre_processed_data_dir + 'genotype.txt'
distance=10000
#for day in range(4):
for day in range(1):
	print(day)
	# Input files
	cell_type_sc_expression_file = pre_processed_data_dir + 'standardized_normalized_per_day_' + str(day) + '_expression.txt'
	cell_type_sc_sample_covariate_file = pre_processed_data_dir + 'cell_covariates_day_' + str(day) + '.txt'
	# Output files
	cell_type_eqtl_variant_gene_pairs_file = per_time_step_eqtl_input_data_dir + 'day_' + str(day) + '_eqtl_input_variant_gene_pairs.txt'
	cell_type_eqtl_expression_file = per_time_step_eqtl_input_data_dir + 'day_' + str(day) + '_eqtl_input_expression.txt'
	cell_type_eqtl_genotype_file = per_time_step_eqtl_input_data_dir + 'day_' + str(day) + '_eqtl_input_genotype.txt'
	cell_type_sample_overlap_file = per_time_step_eqtl_input_data_dir +  'day_' + str(day) + '_eqtl_input_sample_overlap.txt'
	generate_cell_type_eqtl_input_files(day, genotype_file, gene_names_file, cell_type_sc_expression_file, cell_type_sc_sample_covariate_file, cell_type_eqtl_variant_gene_pairs_file, cell_type_eqtl_expression_file, cell_type_eqtl_genotype_file, distance, gene_annotation_file, cell_type_sample_overlap_file)