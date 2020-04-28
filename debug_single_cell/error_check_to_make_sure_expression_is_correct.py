import numpy as np 
import os
import sys
import pdb
import scanpy as sc


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

def get_sc_cell_names(sc_sample_names_file):
	f = open(sc_sample_names_file)
	head_count = 0
	sample_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[10]
		sample_names.append(cell_id)
	f.close()
	return np.asarray(sample_names)

def extract_gene_name_corresponding_to_this_test(pseudobulk_variant_gene_pair_file, test_num):
	f = open(pseudobulk_variant_gene_pair_file)
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if counter == test_num:
			gene_id = data[0]
		counter = counter + 1
	f.close()
	return gene_id

def get_processed_expression(pseudobulk_genotype_file, test_num):
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

def regenerate_pseudobulk(adata, cell_type, gene_id, indi_names):
	readable_ct = ' '.join(cell_type.split('_'))
	counts = []
	gene_arr = np.where(adata.var.index == gene_id)[0]
	if len(gene_arr) != 1:
		print('assumption error')
		pdb.set_trace()
	gene_index = np.where(adata.var.index == gene_id)[0][0]
	for indi in indi_names:
		indices = (adata.obs.ind_cov == indi) & (adata.obs.ct_cov == readable_ct)
		positions = np.where(indices == True)[0]
		counts.append(sum(adata.raw.X[positions, gene_index].toarray())[0])
	return np.asarray(counts)

def regenerate_sc(adata, gene_id, cell_names):
	gene_arr = np.where(adata.var.index == gene_id)[0]
	if len(gene_arr) != 1:
		print('assumption error!')
		pdb.set_trace()
	gene_index = np.where(adata.var.index == gene_id)[0][0]
	mapping = {}
	for i, cell_id in enumerate(adata.obs.cell_id):
		mapping[cell_id] = i
	counts = []
	counter = 0
	for cell_name in cell_names:
		#cell_pos_arr = np.where(adata.obs.cell_id == cell_name)[0]
		#if len(cell_pos_arr) != 1:
		#	print('assumption error')
		#	pdb.set_trace()
		#cell_pos = cell_pos_arr[0]
		cell_pos = mapping[cell_name]
		counts.append(adata.raw.X[cell_pos, gene_index])
	return np.asarray(counts)


def check_pseudobulk_expression(pseudobulk_variant_gene_pair_file, pseudobulk_expression_file, sample_names_file, test_num, adata, cell_type):
	indi_names = get_pseudobulk_sample_names(sample_names_file)

	gene_id = extract_gene_name_corresponding_to_this_test(pseudobulk_variant_gene_pair_file, test_num)

	processed_expression = get_processed_expression(pseudobulk_expression_file, test_num)


	regenerated_expression = regenerate_pseudobulk(adata, cell_type, gene_id, indi_names)

	if np.array_equal(regenerated_expression, processed_expression) == False:
		print('assumption error!')
		pdb.set_trace()

def check_sc_expression(sc_variant_gene_pair_file, sc_expression_file, sc_sample_names_file, test_num, adata, cell_type):
	cell_names = get_sc_cell_names(sc_sample_names_file)

	gene_id = extract_gene_name_corresponding_to_this_test(sc_variant_gene_pair_file, test_num)

	processed_expression = get_processed_expression(sc_expression_file, test_num)

	regenerated_expression = regenerate_sc(adata, gene_id, cell_names)


	if np.array_equal(regenerated_expression, processed_expression) == False:
		print('assumption error!')
		pdb.set_trace()




data_dir = sys.argv[1]
expression_file = sys.argv[2]

cell_type = 'B_cells'
pseudobulk_variant_gene_pair_file = data_dir + 'B_cells_pseudobulk_eqtl_variant_gene_pairs.txt'
pseudobulk_expression_file = data_dir + 'B_cells_pseudobulk_raw_expression.txt'
sample_names_file = data_dir + 'B_cells_pseudobulk_sample_covariates.txt'

sc_variant_gene_pair_file = data_dir + 'B_cells_sc_eqtl_variant_gene_pairs.txt'
sc_expression_file = data_dir + 'B_cells_sc_raw_expression.txt'
sc_sample_names_file = data_dir + 'B_cells_sc_cell_covariates.txt'

adata = sc.read_h5ad(expression_file)
#adata = []
print('loaded')

test_num = 0
check_sc_expression(sc_variant_gene_pair_file, sc_expression_file, sc_sample_names_file, test_num, adata, cell_type)
check_pseudobulk_expression(pseudobulk_variant_gene_pair_file, pseudobulk_expression_file, sample_names_file, test_num, adata, cell_type)


test_num = 1
check_sc_expression(sc_variant_gene_pair_file, sc_expression_file, sc_sample_names_file, test_num, adata, cell_type)
check_pseudobulk_expression(pseudobulk_variant_gene_pair_file, pseudobulk_expression_file, sample_names_file, test_num, adata, cell_type)


test_num = 18
check_sc_expression(sc_variant_gene_pair_file, sc_expression_file, sc_sample_names_file, test_num, adata, cell_type)
check_pseudobulk_expression(pseudobulk_variant_gene_pair_file, pseudobulk_expression_file, sample_names_file, test_num, adata, cell_type)


test_num = 201
check_sc_expression(sc_variant_gene_pair_file, sc_expression_file, sc_sample_names_file, test_num, adata, cell_type)
check_pseudobulk_expression(pseudobulk_variant_gene_pair_file, pseudobulk_expression_file, sample_names_file, test_num, adata, cell_type)

test_num = 107
check_sc_expression(sc_variant_gene_pair_file, sc_expression_file, sc_sample_names_file, test_num, adata, cell_type)
check_pseudobulk_expression(pseudobulk_variant_gene_pair_file, pseudobulk_expression_file, sample_names_file, test_num, adata, cell_type)