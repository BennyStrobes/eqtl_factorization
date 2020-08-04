import numpy as np 
import os
import sys
import pdb
import pandas as pd
from sklearn import linear_model
from pymer4.models import Lmer
import statsmodels.api as sm
import time
import scipy.stats



# Generate eqtl effect size (beta) and p-value for one test
def run_eqtl_one_test_lm(expression, genotype, covariates):
	# Covariate matrix
	X = np.vstack((genotype,covariates.T)).T
	# Add intercept
	X2 = sm.add_constant(X)

	est = sm.OLS(expression, X2)
	est2 = est.fit()

	beta = est2.params[1]
	standard_error = est2.bse[1]
	pvalue = est2.pvalues[1]
	residual_scale = est2.scale
	ll_full = est2.llf


	# FIT NULL MODEL
	X3 = sm.add_constant(covariates)
	model = sm.OLS(expression, X3)
	fit = model.fit()
	ll_null = fit.llf


	pseudo_r_squared = 1.0 - (ll_full/ll_null)

	return beta, standard_error, pvalue, residual_scale, pseudo_r_squared

# Generate eqtl effect size (beta) and p-value for one test
def run_eqtl_one_test_lmm(expression, genotype, covariates, groups):
	num_cov = covariates.shape[1]
	# Covariate matrix
	X = np.vstack((expression, groups, genotype,covariates.T)).T
	# Create column names
	cov_names = ['cov'+str(i) for i in range(num_cov)]
	col_names = ['y', 'group', 'g'] + cov_names

	# Make df
	df = pd.DataFrame(X, columns=col_names)
	# Make formula for LMM
	if num_cov > 0:
		formula = 'y ~ g + ' + ' + '.join(cov_names) + ' + (1 | group)'
	else:
		formula = 'y ~ g + ' + '(1 | group)'

	model = Lmer(formula, data=df)
	model.fit()
	
	beta = model.coefs['Estimate'][1]
	standard_error = model.coefs['SE'][1]
	pvalue = model.coefs['P-val'][1]
	#t_value = fit['T-stat'][1]
	#normal_approx_p = 2.0*(1.0 - scipy.stats.norm.cdf(abs(t_value)))
	residual_scale = model.ranef_var.Std[1]
	return beta, standard_error, pvalue, residual_scale



# Generate eqtl effect size (beta) and p-value for one test
def run_eqtl_one_test2(expression, genotype, covariates):
	# Fit linear model of covariates onto expression
	model = linear_model.LinearRegression(fit_intercept=True) 
	modelfit = model.fit(covariates, expression)
	# Get predicted expression according to the LM
	pred = modelfit.predict(covariates)
	# Get residual expression
	resid_expression = expression - pred

	# Fit linear model of covariates onto expression
	model = linear_model.LinearRegression(fit_intercept=True) 
	modelfit = model.fit(covariates, genotype)
	# Get predicted expression according to the LM
	pred = modelfit.predict(covariates)
	# Get residual expression
	resid_genotype = genotype - pred
	#X2 = sm.add_constant(resid_genotype)
	est = sm.OLS(resid_expression, resid_genotype)
	est2 = est.fit()
	beta = est2.params[0]
	standard_error = est2.bse[0]
	return beta, standard_error

# Run eQTL analysis
def eqtl_analysis(covariate_file, test_names_file, expression_file, genotype_file, sample_overlap_file, num_pcs, variant_gene_pairs_eqtl_results_file, start_number, end_number):
	# Load in covariates (fixed across all tests)
	covariates = np.loadtxt(covariate_file)
	covariates = covariates[:,:num_pcs]
	# Load in individual names (for re term)
	individuals = np.loadtxt(sample_overlap_file)
	# Open up file handles
	test_name_handle = open(test_names_file)
	expression_handle = open(expression_file)
	genotype_handle = open(genotype_file)
	# Output file handle
	t = open(variant_gene_pairs_eqtl_results_file, 'w')

	# Loop through tests
	head_count = 0  # Used to skip header
	counter = 0
	for line in test_name_handle:
		test_name_line  = line.rstrip()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			t.write(test_name_line + '\tbeta\tstandard_error\tpvalue\tresidual_scale\tbeta_lm\tstandard_error_lm\tpvalue_lm\tresidual_scale_lm\tpseudo_r_squared_lm\n')
			continue
		if counter < start_number or counter > end_number:
			counter = counter + 1
			temp_e = expression_handle.readline()
			temp_g = genotype_handle.readline()
			continue
		print(counter)
		counter = counter + 1
		expression = np.asarray(expression_handle.readline().rstrip().split('\t')).astype(float)
		genotype = np.asarray(genotype_handle.readline().rstrip().split('\t')).astype(float)
		beta, std_err, pvalue, residual_scale = run_eqtl_one_test_lmm(expression, genotype, covariates, individuals)
		beta_lm, std_err_lm, pvalue_lm, residual_scale_lm, pseudo_r_squared_lm = run_eqtl_one_test_lm(expression, genotype, covariates)
		t.write(test_name_line + '\t' + str(beta) + '\t' + str(std_err) + '\t' + str(pvalue) + '\t' + str(residual_scale) + '\t' + str(beta_lm) + '\t' + str(std_err_lm) + '\t' + str(pvalue_lm) + '\t' + str(residual_scale_lm) + '\t' + str(pseudo_r_squared_lm) + '\n')
		if np.mod(counter, 20) == 0:
			t.flush()
	t.close()
	test_name_handle.close()
	expression_handle.close()
	genotype_handle.close()

def bf_fdr_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, fdr_thresh):
	f = open(variant_gene_pairs_eqtl_results_file)
	t = open(multple_testing_correction_results_file, 'w')
	head_count = 0
	genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\tnum_snps_in_gene\tfdr\n')
			continue
		gene_id = data[0]
		variant_id = data[1]
		pvalue = float(data[7])
		if gene_id not in genes:
			genes[gene_id] = (variant_id, pvalue, 1, line)
		else:
			old_pvalue = genes[gene_id][1]
			old_count = genes[gene_id][2]
			if pvalue <= old_pvalue:
				genes[gene_id] = (variant_id, pvalue, old_count+1, line)
			else:
				genes[gene_id] = (genes[gene_id][0], genes[gene_id][1], old_count+1, genes[gene_id][3])
	f.close()
	# Loop through genes and do BF correction
	bf_gene_array = []
	for gene in genes.keys():
		lead_variant = genes[gene][0]
		lead_nominal_pvalue = genes[gene][1]
		num_variants_at_gene = genes[gene][2]
		test_line = genes[gene][3]
		bf_corrected_pvalue = lead_nominal_pvalue*num_variants_at_gene
		if bf_corrected_pvalue > 1.0:
			bf_corrected_pvalue = 1.0
		bf_gene_array.append((bf_corrected_pvalue, lead_variant, gene, num_variants_at_gene, test_line))
	sorted_bf_gene_array = sorted(bf_gene_array, key=lambda tup: tup[0])
	# BH correction
	kk = 1
	num_genes = len(sorted_bf_gene_array)
	sig = True
	for gene_tuple in sorted_bf_gene_array:
		bf_pvalue = gene_tuple[0]
		fdr = num_genes*bf_pvalue/kk 
		kk = kk + 1
		if fdr > fdr_thresh:
			sig = False
		if sig == True:
			t.write(gene_tuple[4] + '\t' + str(gene_tuple[3]) + '\t' + str(fdr) + '\n')
	t.close()

def bf_top_nn_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, nn_thresh):
	f = open(variant_gene_pairs_eqtl_results_file)
	t = open(multple_testing_correction_results_file, 'w')
	head_count = 0
	genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\tnum_snps_in_gene\tfdr\n')
			continue
		gene_id = data[0]
		variant_id = data[1]
		pvalue = float(data[7])
		if gene_id not in genes:
			genes[gene_id] = (variant_id, pvalue, 1, line)
		else:
			old_pvalue = genes[gene_id][1]
			old_count = genes[gene_id][2]
			if pvalue <= old_pvalue:
				genes[gene_id] = (variant_id, pvalue, old_count+1, line)
			else:
				genes[gene_id] = (genes[gene_id][0], genes[gene_id][1], old_count+1, genes[gene_id][3])
	f.close()
	# Loop through genes and do BF correction
	bf_gene_array = []
	for gene in genes.keys():
		lead_variant = genes[gene][0]
		lead_nominal_pvalue = genes[gene][1]
		num_variants_at_gene = genes[gene][2]
		test_line = genes[gene][3]
		bf_corrected_pvalue = lead_nominal_pvalue*num_variants_at_gene
		if bf_corrected_pvalue > 1.0:
			bf_corrected_pvalue = 1.0
		bf_gene_array.append((bf_corrected_pvalue, lead_variant, gene, num_variants_at_gene, test_line))
	sorted_bf_gene_array = sorted(bf_gene_array, key=lambda tup: tup[0])
	# BH correction
	kk = 1
	num_genes = len(sorted_bf_gene_array)
	sig = True
	for gene_tuple in sorted_bf_gene_array:
		bf_pvalue = gene_tuple[0]
		fdr = num_genes*bf_pvalue/kk 
		kk = kk + 1
		if kk > (nn_thresh+1):
			sig = False
		if sig == True:
			t.write(gene_tuple[4] + '\t' + str(gene_tuple[3]) + '\t' + str(fdr) + '\n')
	t.close()


def get_number_of_tests(file_name):
	f = open(file_name)
	count = -1
	for line in f:
		line = line.rstrip()
		count = count + 1
	f.close()
	return count

# For parallelization purposes
def parallelization_start_and_end(num_tasks, job_number, total_jobs):
	tasks_per_job = int(num_tasks/total_jobs) + 1
	start_task = job_number*tasks_per_job
	end_task = (job_number + 1)*tasks_per_job -1 
	return start_task, end_task

#####################
# Command line args
#####################
expression_file = sys.argv[1]  # Input file of dimension num_testsXnum_samples
genotype_file = sys.argv[2]  # Input file of dimension num_testsXnum_samples
test_names_file = sys.argv[3]
covariate_file = sys.argv[4]  # Input file of dimension num_samplesXnum_covariates
sample_overlap_file = sys.argv[5]  # File containing info on which samples came from which individuals (for random effects)
num_pcs = int(sys.argv[6])  # Integer containing number of pcs (covariates to include)
output_root = sys.argv[7]  # output root
job_number = int(sys.argv[8])
total_jobs = int(sys.argv[9])

#For parallelization purposes
number_of_tests = get_number_of_tests(test_names_file)
start_number, end_number = parallelization_start_and_end(number_of_tests, job_number, total_jobs)
print(end_number-start_number)

####################
# Run eQTL analysis
####################
# Output file
variant_gene_pairs_eqtl_results_file = output_root + 'all_variant_gene_pairs_' + str(job_number) + '_' + str(total_jobs) + '.txt'
eqtl_analysis(covariate_file, test_names_file, expression_file, genotype_file, sample_overlap_file, num_pcs, variant_gene_pairs_eqtl_results_file, start_number, end_number)

####################
# Multiple-testing correction
####################
# Output file
fdr_thresh=.05
multple_testing_correction_results_file = output_root + 'multiple_testing_bf_bh_' + str(fdr_thresh) + '_fdr_' + '.txt'
# Perform bonferonni correction at variant level (per gene) and then BH at gene level
# bf_fdr_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, fdr_thresh)


####################
# Multiple-testing correction
####################
# Output file
fdr_thresh=.1
multple_testing_correction_results_file = output_root + 'multiple_testing_bf_bh_' + str(fdr_thresh) + '_fdr_' + '.txt'
# Perform bonferonni correction at variant level (per gene) and then BH at gene level
# bf_fdr_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, fdr_thresh)


####################
# Top nn genes
####################
# Output file
nn_thresh= 800
multple_testing_correction_results_file = output_root + 'multiple_testing_bf_top_' + str(nn_thresh) + '_genes.txt'
# Perform bonferonni correction at variant level (per gene) and then BH at gene level
# bf_top_nn_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, nn_thresh)