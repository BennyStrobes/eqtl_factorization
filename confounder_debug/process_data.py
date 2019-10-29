import numpy as np 
import os
import sys
import pdb
from sklearn import linear_model


def regress_out_covs(raw_gtex_data, pc_file, residual_expression_file):
	expr_full = np.loadtxt(raw_gtex_data, dtype=str)
	cov_full = np.loadtxt(pc_file, dtype=str)
	expr = expr_full[1:,4:].astype(float)
	cov = cov_full[1:,1:].astype(float)
	num_genes = expr.shape[0]

	t = open(residual_expression_file,'w')
	t.write('\t'.join(expr_full[0,:]) + '\n')

	for gene_num in range(num_genes):
		expr_vec = expr[gene_num,:]
		model = linear_model.LinearRegression(fit_intercept=True) 
		# fit model
		modelfit = model.fit(np.transpose(cov),np.transpose(np.asmatrix(expr_vec)))
		beta = modelfit.coef_
		residual_expression = np.zeros(len(expr_vec))

		for i,ele in enumerate(expr_vec):
			residual_expression[i] = ele
		for i,val in enumerate(beta[0]):
			residual_expression = residual_expression - val*cov[i,:]
		# residual_expression = (residual_expression - np.mean(residual_expression))/np.std(residual_expression)

		t.write('\t'.join(expr_full[gene_num+1,0:4]) + '\t' + '\t'.join(residual_expression.astype(str)) + '\n')
	t.close()



def filter_genes(input_file, gene_subset, output_file):
	genes = {}
	count = 0
	f = open(gene_subset)
	for line in f:
		line = line.rstrip()
		data = line.split()
		genes[line.split(':')[0]] = 1
		count = count +1
	f.close()
	f = open(input_file)
	t = open(output_file, 'w')
	head_count = 0

	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(data[3:]) + '\n')
			continue
		ensamble = data[3]
		if ensamble not in genes:
			continue
		t.write('\t'.join(data[3:]) + '\n')
	t.close()
	f.close()

def filter_genes2(input_file, gene_subset, output_file):
	genes = {}
	count = 0
	f = open(gene_subset)
	for line in f:
		line = line.rstrip()
		data = line.split()
		genes[line.split(':')[0]] = 1
		count = count +1
	f.close()
	f = open(input_file)
	t = open(output_file, 'w')
	head_count = 0

	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(data) + '\n')
			continue
		ensamble = data[0]
		if ensamble not in genes:
			continue
		t.write('\t'.join(data) + '\n')
	t.close()
	f.close()

def reorder_samples_according_to_extraction_type(input_file, cov_file, output_file):
	covs = np.loadtxt(cov_file, dtype=str, delimiter='\t')
	cohort = covs[1:,0]
	sample_names = np.loadtxt('/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl_ancestry/processed_data/Adipose_Subcutaneous_sample_names.txt', dtype=str)

	new_sample_order = []
	for i, val in enumerate(sample_names):
		indi_cohort = cohort[i]
		if indi_cohort == 'Postmortem':
			new_sample_order.append(val)
	for i, val in enumerate(sample_names):
		indi_cohort = cohort[i]
		if indi_cohort != 'Postmortem':
			new_sample_order.append(val)

	f = open(input_file)
	t = open(output_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			indi_ids = data[1:]
			indices = []
			for sampler in new_sample_order:
				for i, indi in enumerate(indi_ids):
					if indi == sampler:
						indices.append(i)
			#t.write(data[0] + '\t' + '\t'.join(np.asarray(indi_ids)[indices]) + '\n')
			continue		
		t.write( '\t'.join(np.asarray(data[1:])[indices]) + '\n')
	f.close()
	t.close()


raw_gtex_data = sys.argv[1]
cov_file = sys.argv[2]
pc_file = sys.argv[3]
gene_subset = sys.argv[4]
output_dir = sys.argv[5]



residual_expression_file = output_dir + 'residual_expression.txt'
#regress_out_covs(raw_gtex_data, pc_file, residual_expression_file)


raw_gene_filtered = output_dir + 'raw_gene_filtered.txt'
filter_genes(raw_gtex_data, gene_subset, raw_gene_filtered)

resid_gene_filtered = output_dir + 'resid_gene_filtered.txt'
filter_genes(residual_expression_file, gene_subset, resid_gene_filtered)

my_data = '/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl_ancestry/processed_data/' + 'Adipose_Subcutaneous_residual_cross_tissue_tpm_standardized.txt'
my_data_filtered = output_dir + 'personal_resid_gene_filtered.txt'
filter_genes2(my_data, gene_subset, my_data_filtered)


reorder_samples_according_to_extraction_type(my_data_filtered, cov_file, output_dir + 'personal_resid_gene_filtered_reordered.txt')

reorder_samples_according_to_extraction_type(raw_gene_filtered, cov_file, output_dir + 'raw_gene_filtered_reordered.txt')

reorder_samples_according_to_extraction_type(resid_gene_filtered, cov_file, output_dir + 'resid_gene_filtered_reordered.txt')