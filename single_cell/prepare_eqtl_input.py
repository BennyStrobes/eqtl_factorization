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
	genes = np.loadtxt(gene_file, delimiter='\n',dtype=str)
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
			# List of genes that variant maps to
			mapping_genes = chromosome[variant_pos].split(':')
			for gene_id in mapping_genes:
				# THIS IS A VARIANT-GENE PAIR WE WILL TEST
				# PRINT TO OUTPUT
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

########################
# Step 2: Run eQTL analysis across all cells
########################
def run_eqtl_analysis_across_all_cells(expression_file, cell_file, gene_file, variant_gene_pair_file, eqtl_output_file):
	# Extract gene list
	genes = np.loadtxt(gene_file, delimiter='\n',dtype=str)
	# Load in residual expression
	#expr = load_expression_file(expression_file)
	pdb.set_trace()


######################
# Command line args
######################
gene_annotation_file = sys.argv[1]
processed_expression_dir = sys.argv[2]
eqtl_input_dir = sys.argv[3]
genotype_data_dir = sys.argv[4]


# Normalization method used in preprocesssing
normalization_method = 'sctransform_with_covariates_and_individual'

########################
# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
########################
# Variant must be within $distance BP from TSS of gene
distance = 10000
# File containing availible genes
gene_file = processed_expression_dir + normalization_method + '_gene_names.txt'
# Output file containing list of variant gene pairs
variant_gene_pair_file = eqtl_input_dir + 'variant_gene_pairs_' + str(distance) + '_bp.txt'

#extract_variant_gene_pairs_for_eqtl_testing(gene_file, gene_annotation_file, distance, genotype_data_dir, variant_gene_pair_file)


########################
# Step 2: Run eQTL analysis across all cells
########################
# Expression file
expression_file = processed_expression_dir + normalization_method + '_residual_expression.txt'
# Cell file: File containing column names of expression file
cell_file = processed_expression_dir + normalization_method + '_cell_info.txt'
# Output file containing list of variant gene pairs and their corresponding p-values
eqtl_output_file = eqtl_input_dir + 'variant_gene_pairs_' + str(distance) + '_bp_eqtl_results_across_all_cells.txt'

run_eqtl_analysis_across_all_cells(expression_file, cell_file, gene_file, variant_gene_pair_file, eqtl_output_file)
