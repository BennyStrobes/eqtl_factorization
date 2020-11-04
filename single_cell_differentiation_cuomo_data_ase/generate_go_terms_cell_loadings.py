import numpy as np 
import os
import sys
import pdb
import gzip
import h5py



def load_in_expression_and_save_as_h5_file(standardized_total_expression_file, h5_expression_file):
	f = open(standardized_total_expression_file)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		arr.append(data)
	f.close()
	mat = np.asarray(arr)
	h5f = h5py.File(h5_expression_file, 'w')
	h5f.create_dataset('data', data=mat)
	h5f.close()

def get_gene_symbols(gene_ids_full):
	gene_symbols = []
	for ele in gene_ids_full:
		info = ele.split('_')
		if len(info) != 2:
			print('erororo')
			pdb.set_trace()
		gene_symbols.append(info[1])
	return np.asarray(gene_symbols)



######################
# Command Line Args
######################
standardized_total_expression_file = sys.argv[1]
go_terms_file = sys.argv[2]
gencode_gene_annotation_file = sys.argv[3]
processed_data_dir = sys.argv[4]

h5_expression_file = processed_data_dir + 'standardized_cell_expression.h5'
#load_in_expression_and_save_as_h5_file(standardized_total_expression_file, h5_expression_file)
exp_full = np.asarray(h5py.File(h5_expression_file,'r')['data'])
cells = exp_full[0,1:]
gene_ids_full = exp_full[1:,0]
exp = exp_full[1:,1:].astype(float)
gene_symbols = get_gene_symbols(gene_ids_full)

# Open output file handle
go_term_cell_loading_output_file = processed_data_dir + 'go_terms_cell_loadings.txt'
t = open(go_term_cell_loading_output_file,'w')
# print header
t.write('go_term\t' + '\t'.join(cells) + '\n')

f = open(go_terms_file)
count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	go_term = data[0]
	gene_set = data[2:]
	indices = []
	for gene in gene_set:
		positions = np.where(gene_symbols==gene)[0]
		if len(positions) == 0:
			continue
		elif len(positions) == 1:
			indices.append(positions[0])
		else:
			print('assumptoiner oeroor')
			pdb.set_trace()
	# get indices of genes
	indices = np.asarray(indices)
	if len(indices) < 2:
		continue
	# run pca on gene subset
	uuu, sss, vh = np.linalg.svd(exp[indices,:], full_matrices=False)
	pc1 = np.transpose(vh)[:,0]
	t.write(go_term + '\t' + '\t'.join(pc1.astype(str)) + '\n')
	count = count + 1
	print(count)
f.close()
t.close()



#exp = np.loadtxt(standardized_total_expression_file, dtype=str, delimiter='\t', comments='*')
#pdb.set_trace()