import numpy as np 
import os
import sys
import pdb

def get_pseudobulk_sample_names(sample_names_file):
	f = open(sample_names_file)
	head_count = 0
	sample_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_names.append(data[3])
	f.close()
	return np.asarray(sample_names)

def get_sc_sample_names(sample_names_file):
	f = open(sample_names_file)
	head_count = 0
	sample_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_names.append(data[3])
	f.close()
	return np.asarray(sample_names)

def extract_variant_name_corresponding_to_this_test(pseudobulk_variant_gene_pair_file, test_num):
	f = open(pseudobulk_variant_gene_pair_file)
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if counter == test_num:
			variant_id = data[1]
		counter = counter + 1
	f.close()
	return variant_id

def get_processed_genotype(pseudobulk_genotype_file, test_num):
	f = open(pseudobulk_genotype_file)
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if counter == test_num:
			genotype = np.asarray(data).astype(float)
		counter = counter + 1
	f.close()
	return genotype

def get_genotype_and_ordered_individuals(genotype_file, variant_id):
	f = open(genotype_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			if len(data) != 119:
				print('assumption error!')
				pdb.set_trace()
			ordered_individuals = np.asarray(data)
			continue
		if len(data) != 120:
			print('assumption error!')
			pdb.set_trace()
		line_variant_id = data[0]
		if line_variant_id != variant_id:
			continue
		genotype = np.asarray(data[1:]).astype(float)
	f.close()
	return genotype, ordered_individuals


def check_sc_genotype(sc_variant_gene_pair_file, sc_genotype_file, sc_sample_names_file, genotype_data_dir, test_num):
	indi_names = get_sc_sample_names(sc_sample_names_file)

	variant_id = extract_variant_name_corresponding_to_this_test(sc_variant_gene_pair_file, test_num)

	chrom_num = variant_id.split(':')[0]

	processed_genotype = get_processed_genotype(sc_genotype_file, test_num)

	genotype_file = genotype_data_dir + 'chr' + str(chrom_num) + '.genotypes.matrix.eqtl.txt'
	raw_genotype, raw_ordered_individuals = get_genotype_and_ordered_individuals(genotype_file, variant_id)
	
	reordering = []
	for indi_name in indi_names:
		for i, raw_indi in enumerate(raw_ordered_individuals):
			if raw_indi == indi_name:
				reordering.append(i)
	if np.array_equal(np.asarray(raw_ordered_individuals[reordering]), indi_names) == False:
		print('assumption error')
		pdb.set_trace()
	if np.array_equal(np.asarray(raw_genotype[reordering]), processed_genotype) == False:
		print('assumption error')
		pdb.set_trace()

def check_pseudobulk_genotype(pseudobulk_variant_gene_pair_file, pseudobulk_genotype_file, sample_names_file, genotype_data_dir, test_num):
	indi_names = get_pseudobulk_sample_names(sample_names_file)

	variant_id = extract_variant_name_corresponding_to_this_test(pseudobulk_variant_gene_pair_file, test_num)
	chrom_num = variant_id.split(':')[0]

	processed_genotype = get_processed_genotype(pseudobulk_genotype_file, test_num)

	genotype_file = genotype_data_dir + 'chr' + str(chrom_num) + '.genotypes.matrix.eqtl.txt'
	raw_genotype, raw_ordered_individuals = get_genotype_and_ordered_individuals(genotype_file, variant_id)
	

	reordering = []
	for indi_name in indi_names:
		for i, raw_indi in enumerate(raw_ordered_individuals):
			if raw_indi == indi_name:
				reordering.append(i)

	if np.array_equal(np.asarray(raw_ordered_individuals[reordering]), indi_names) == False:
		print('assumption error')
		pdb.set_trace()
	if np.array_equal(np.asarray(raw_genotype[reordering]), processed_genotype) == False:
		print('assumption error')
		pdb.set_trace()





genotype_data_dir = sys.argv[1]
data_dir = sys.argv[2]
#chr5.genotypes.matrix.eqtl.txt


pseudobulk_variant_gene_pair_file = data_dir + 'B_cells_pseudobulk_eqtl_variant_gene_pairs.txt'
pseudobulk_genotype_file = data_dir + 'B_cells_pseudobulk_genotype.txt'
sample_names_file = data_dir + 'B_cells_pseudobulk_sample_covariates.txt'

sc_variant_gene_pair_file = data_dir + 'B_cells_sc_eqtl_variant_gene_pairs.txt'
sc_genotype_file = data_dir + 'B_cells_sc_genotype.txt'
sc_sample_names_file = data_dir + 'B_cells_sc_cell_covariates.txt'

test_num = 0
check_sc_genotype(sc_variant_gene_pair_file, sc_genotype_file, sc_sample_names_file, genotype_data_dir, test_num)
check_pseudobulk_genotype(pseudobulk_variant_gene_pair_file, pseudobulk_genotype_file, sample_names_file, genotype_data_dir, test_num)





test_num = 1
check_sc_genotype(sc_variant_gene_pair_file, sc_genotype_file, sc_sample_names_file, genotype_data_dir, test_num)
check_pseudobulk_genotype(pseudobulk_variant_gene_pair_file, pseudobulk_genotype_file, sample_names_file, genotype_data_dir, test_num)


test_num = 2
check_sc_genotype(sc_variant_gene_pair_file, sc_genotype_file, sc_sample_names_file, genotype_data_dir, test_num)
check_pseudobulk_genotype(pseudobulk_variant_gene_pair_file, pseudobulk_genotype_file, sample_names_file, genotype_data_dir, test_num)





test_num = 45
check_sc_genotype(sc_variant_gene_pair_file, sc_genotype_file, sc_sample_names_file, genotype_data_dir, test_num)
check_pseudobulk_genotype(pseudobulk_variant_gene_pair_file, pseudobulk_genotype_file, sample_names_file, genotype_data_dir, test_num)
