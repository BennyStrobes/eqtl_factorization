import numpy as np 
import os
import sys
import pdb
import pandas as pd
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from pymer4.models import Lmer
# import statsmodels.api as sm
import time
import scipy.stats
#from statsmodels.stats.diagnostic import het_white
from statsmodels.stats.diagnostic import het_breuschpagan

def run_eqtl_one_test_lm(expression, genotype, covariates):
	X = np.vstack((genotype, covariates.T)).T
	reg = LinearRegression().fit(X, expression)
	beta = reg.coef_[0]
	return beta

def run_eqtl_on_residual_expression_one_test_lm(expression, genotype):
	X = np.vstack((genotype))
	reg = LinearRegression().fit(X, expression)
	beta = reg.coef_[0]
	return beta

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
	#residual_scale = model.ranef_var.Std[1]
	return beta, standard_error, pvalue

# Generate eqtl effect size (beta) and p-value for one test
def run_dynamic_eqtl_one_test_lmm(expression, genotype, covariates, groups, environmental_variable):
	num_cov = covariates.shape[1]
	# Covariate matrix
	X = np.vstack((expression, groups, genotype, environmental_variable, environmental_variable*genotype, covariates.T)).T
	# Create column names
	cov_names = ['cov'+str(i) for i in range(num_cov)]
	col_names = ['y', 'group', 'g', 'e', 'gXe'] + cov_names

	# Make df
	df = pd.DataFrame(X, columns=col_names)
	# Make formula for LMM
	if num_cov > 0:
		formula = 'y ~ g + e + gXe + ' + ' + '.join(cov_names) + ' + (1 | group)'
	else:
		formula = 'y ~ g + e + gXe + ' + '(1 | group)'

	model = Lmer(formula, data=df)
	model.fit()
	
	beta = model.coefs['Estimate'][3]
	standard_error = model.coefs['SE'][3]
	pvalue = model.coefs['P-val'][3]
	#t_value = fit['T-stat'][1]
	#normal_approx_p = 2.0*(1.0 - scipy.stats.norm.cdf(abs(t_value)))
	#residual_scale = model.ranef_var.Std[1]
	return pvalue

# Generate eqtl effect size (beta) and p-value for one test
def run_dynamic_eqtl_one_test_lm(expression, genotype, covariates, groups, environmental_variable):
	num_cov = covariates.shape[1]
	# Covariate matrix
	X = np.vstack((genotype, environmental_variable, environmental_variable*genotype, covariates.T)).T

	#reg = LinearRegression().fit(X, expression)
	#beta = reg.coef_[2]
	X2 = sm.add_constant(X)
	est = sm.MixedLM(endog=expression, exog=X2, groups=groups).fit()

	#t_value = fit['T-stat'][1]
	#normal_approx_p = 2.0*(1.0 - scipy.stats.norm.cdf(abs(t_value)))
	#residual_scale = model.ranef_var.Std[1]
	return est.pvalues[3]

def get_bootstrapped_indices(individuals, individual_to_cells, sampling_fraction):
#	num_sampled_cells = np.ceil(len(individuals)*sampling_fraction)
#	new_indices = np.random.choice(range(len(individuals)),size=int(num_sampled_cells), replace=False)
#	return new_indices
	new_indices = []
	for cell_arr in individual_to_cells:
		num_sampled_cells = np.ceil(len(cell_arr)*sampling_fraction)
		temp_indices = np.random.choice(cell_arr,size=int(num_sampled_cells), replace=False)
		new_indices.append(temp_indices)
	return np.hstack(new_indices)


def run_bootstrapped_eqtl_lmm_stability_one_test(expression, genotype, covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction):
	num_cov = covariates.shape[1]

	# Covariate matrix
	X = np.vstack((expression, individuals, genotype,covariates.T)).T
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

	bootstrapped_betas = []
	for bootstrap_num in range(num_bootstraps):
		print(bootstrap_num)
		indices = get_bootstrapped_indices(individuals, individual_to_cells, sampling_fraction)
		model = Lmer(formula, data=df.iloc[indices,:])
		model.fit()
		bootstrapped_beta = model.coefs['Estimate'][1]
		#bootstrapped_beta, bootstrapped_std_err, bootstrapped_pvalue = run_eqtl_one_test_lmm(expression[indices], genotype[indices], covariates[indices,:], individuals[indices])
		bootstrapped_betas.append(bootstrapped_beta)
	return np.asarray(bootstrapped_betas)

def run_bootstrapped_eqtl_stability_one_test(expression, genotype, covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction):
	
	bootstrapped_betas = []
	for bootstrap_num in range(num_bootstraps):
		indices = get_bootstrapped_indices(individuals, individual_to_cells, sampling_fraction)
		bootstrapped_beta = run_eqtl_one_test_lm(expression[indices], genotype[indices], covariates[indices])
		bootstrapped_betas.append(bootstrapped_beta)
	return np.asarray(bootstrapped_betas)


def regress_out_covariates(expression, covariates):
	X = covariates
	reg = LinearRegression().fit(X, expression)
	residual_expression = np.copy(expression)
	num_cov = X.shape[1]
	# compute residuals
	for cov_num in range(num_cov):
		residual_expression = residual_expression - reg.coef_[cov_num]*X[:,cov_num]
	return residual_expression


def run_bootstrapped_eqtl_stability_with_residuals_one_test(expression, genotype, covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction, seed):
	np.random.seed(seed)
	residual_expression = regress_out_covariates(expression, covariates)
	residual_genotype = regress_out_covariates(genotype, covariates)
	bootstrapped_betas = []
	for bootstrap_num in range(num_bootstraps):
		indices = get_bootstrapped_indices(individuals, individual_to_cells, sampling_fraction)
		bootstrapped_beta = run_eqtl_on_residual_expression_one_test_lm(residual_expression[indices], residual_genotype[indices])
		bootstrapped_betas.append(bootstrapped_beta)
	return np.asarray(bootstrapped_betas)


def run_bootstrapped_eqtl_stability_with_residuals_one_test_v2(expression, genotype, covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction, seed):
	np.random.seed(seed)
	#residual_expression = regress_out_covariates(expression, covariates)
	#residual_genotype = regress_out_covariates(genotype, covariates)

	# Covariate matrix
	num_cov = covariates.shape[1]
	X = np.vstack((expression, individuals.astype(str), genotype,covariates.T)).T
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
	eqtl_pvalue = model.coefs['P-val'][1]
	bp_test = het_breuschpagan(model.residuals, np.vstack(genotype))
	pdb.set_trace()
	#X2 = sm.add_constant(X)
	#reg = LinearRegression().fit(X, expression)
	#est = sm.MixedLM(endog=expression, exog=X2, groups=individuals).fit()
	#est = sm.OLS(expression,X2).fit()
	#eqtl_pvalue = est.pvalues[1]
	#bp_test = het_breuschpagan(est.resid, np.vstack(genotype))
	#bp_test = het_breuschpagan(est.resid, X)
	#white_test = het_white(est.resid,X)
	#print(white_test)
	#model = ols(expression, X)
	#for bootstrap_num in range(num_bootstraps):
	#	indices = get_bootstrapped_indices(individuals, individual_to_cells, sampling_fraction)
	#	bootstrapped_beta = run_eqtl_on_residual_expression_one_test_lm(residual_expression[indices], genotype[indices])
		#bootstrapped_betas.append(bootstrapped_beta)
		#bootstrapped_perm_beta = run_eqtl_on_residual_expression_one_test_lm(residual_expression[indices], np.random.permutation(genotype[indices]))
		#bootstrapped_perm_betas.append(bootstrapped_perm_beta)
	#print(np.max(bootstrapped_betas) - np.min(bootstrapped_betas))
	#print(np.max(bootstrapped_perm_betas) - np.min(bootstrapped_perm_betas))
	#print(np.var(bootstrapped_betas))
	#print(np.var(bootstrapped_perm_betas))
	#print(np.mean(bootstrapped_betas))
	#print(np.mean(bootstrapped_perm_betas))
	return eqtl_pvalue, bp_test[3]



def get_individual_to_cells(individuals):
	arr = []
	unique_indi = np.unique(individuals)
	for indi in unique_indi:
		individual_indices = np.where(individuals == indi)[0]
		arr.append(individual_indices)
	return arr

def compute_maf(genotype):
	af = np.sum(genotype)/(2.0*len(genotype))
	if af > .5:
		maf = 1.0 - af
	else:
		maf = af
	pdb.set_trace()
	return maf

def permute_genotype_by_individual(genotype, individual_to_cells):
	individual_genotype_arr = []
	for individual_cell_arr in individual_to_cells:
		individual_genotype = genotype[individual_cell_arr[0]]
		individual_genotype_arr.append(individual_genotype)
	permuted_individual_genotype_arr = np.random.permutation(individual_genotype_arr)
	permuted_cell_genotype = np.copy(genotype)
	for i, individual_cell_arr in enumerate(individual_to_cells):
		permuted_cell_genotype[individual_cell_arr] = permuted_individual_genotype_arr[i]
	return permuted_cell_genotype

def run_bootstrapped_eqtl_stability(covariate_file, test_names_file, expression_file, genotype_file, sample_overlap_file, environmental_variable_file, num_bootstraps, sampling_fraction, output_file, method):
	# Load in covariates (fixed across all tests)
	covariates = np.loadtxt(covariate_file)
	# Load in individual names (for re term)
	individuals = np.loadtxt(sample_overlap_file)
	# Create list where each element is a vector corresponding to all cells from that individual
	individual_to_cells = get_individual_to_cells(individuals)
	# Load in environmental variable (for dynamic eqtl calling)
	environmental_variable = np.loadtxt(environmental_variable_file)
	# Open up file handles
	test_name_handle = open(test_names_file)
	expression_handle = open(expression_file)
	genotype_handle = open(genotype_file)
	# Output file handle
	t = open(output_file, 'w')

	# Loop through tests
	head_count = 0  # Used to skip header
	counter = 0
	for line in test_name_handle:
		test_name_line  = line.rstrip()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			t.write(test_name_line + '\tdynamic_eqtl_pvalue\teqtl_pvalue\thet_eqtl_pvalue\n')
			continue
		print(counter)
		counter = counter + 1
		expression = np.asarray(expression_handle.readline().rstrip().split('\t')).astype(float)
		genotype = np.asarray(genotype_handle.readline().rstrip().split('\t')).astype(float)
		if method == 'eqtl_beta_variance_relative_to_permutation':
			dynamic_eqtl_pvalue = run_dynamic_eqtl_one_test_lmm(expression, genotype, covariates, individuals, environmental_variable)
			eqtl_pvalue, het_test_pvalue = run_bootstrapped_eqtl_stability_with_residuals_one_test_v2(expression, genotype, covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction, counter)
			t.write(test_name_line + '\t' + str(dynamic_eqtl_pvalue) + '\t' + str(eqtl_pvalue) + '\t' + str(het_test_pvalue) + '\n')
			print(str(dynamic_eqtl_pvalue) + '\t' + str(eqtl_pvalue) + '\t' + str(het_test_pvalue) + '\n')
				#perm_bootstrapped_betas = run_bootstrapped_eqtl_stability_with_residuals_one_test(expression, np.random.permutation(genotype), covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction, counter)
		# Dynamic eQTL analysis (assuming known pseudotime)
		#dynamic_eqtl_beta = run_dynamic_eqtl_one_test_lm(expression, genotype, covariates, individuals, environmental_variable)

		# Bootstrapped eqtl stability analysis
		#bootstrapped_betas = run_bootstrapped_eqtl_lmm_stability_one_test(expression, genotype, covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction)
		#bootstrapped_betas = run_bootstrapped_eqtl_stability_one_test(expression, genotype, covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction)
		#bootstrapped_betas = run_bootstrapped_eqtl_stability_with_residuals_one_test(expression, genotype, covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction, counter)
		#bootstrapped_var = np.var(bootstrapped_betas)

		#eqtl_beta, eqtl_standard_error, eqtl_pvalue = run_eqtl_one_test_lmm(expression, genotype, covariates, individuals)
		#print('#########')
		#print(str(eqtl_beta) + '\t' + str(eqtl_standard_error) + '\t' + str(eqtl_pvalue))

		#emperical_counter = 0
		#num_perms = 10
		#for perm_num in range(num_perms):
		#	np.random.seed()
		#	permuted_genotype = permute_genotype_by_individual(genotype, individual_to_cells)
		#	perm_eqtl_beta, perm_eqtl_standard_error, perm_eqtl_pvalue = run_eqtl_one_test_lmm(expression, permuted_genotype, covariates, individuals)
		#	print(str(perm_eqtl_beta) + '\t' + str(perm_eqtl_standard_error) + '\t' + str(perm_eqtl_pvalue))
			#permuted_genotype = permute_genotype_by_individual(genotype, individual_to_cells)
			#perm_bootstrapped_betas = run_bootstrapped_eqtl_stability_with_residuals_one_test(expression, np.random.permutation(genotype), covariates, individuals, individual_to_cells, num_bootstraps, sampling_fraction, counter)
			#perm_bootstrapped_var = np.var(perm_bootstrapped_betas)
		#	if perm_eqtl_standard_error >= eqtl_standard_error:
		#		emperical_counter = emperical_counter + 1
		#variance_pval = float(emperical_counter)/float(num_perms)

		#print(test_name_line + '\t' + str(dynamic_eqtl_beta) + '\t' + str(variance_pval) + '\n')
		# Print to output file
		#t.write(test_name_line + '\t' + str(dynamic_eqtl_beta) + '\t' + str(variance_pval) + '\n')
		#t2.write(test_name_line + '\t' + ','.join(bootstrapped_betas.astype(str)) + '\n')
		if np.mod(counter, 20) == 0:
			t.flush()
		#	t2.flush()	
	t.close()




bootstrapped_eqtl_stability_dir = sys.argv[1]
method = sys.argv[2]


################
# Input data files
covariate_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_input_data_cis_window_250000_covariate_subset_10.txt'
test_names_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_input_data_cis_window_250000_variant_gene_pairs.txt'
expression_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_input_data_cis_window_250000_expression.txt'
genotype_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_input_data_cis_window_250000_standardized_genotype.txt'
sample_overlap_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_input_data_cis_window_250000_individual_id.txt'
environmental_variable_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_input_data_cis_window_250000_environmental_variable.txt'


#################
# Model parameters
num_bootstraps = 500
sampling_fraction = .05

################
# Output data file
if method == 'residual_variance_relative_to_permutation':
	output_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_' + method + '_results.txt'
elif method == 'eqtl_beta_variance_relative_to_permutation':
	output_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_' + method + '_' + str(num_bootstraps) + '_' + str(sampling_fraction) + '_results.txt'
#output_file = bootstrapped_eqtl_stability_dir + 'bootstrapped_eqtl_stability_' + str(num_bootstraps) + '_samples_' + str(sampling_fraction) + '_sampling_fraction_results.txt'

################
# Run analysis
run_bootstrapped_eqtl_stability(covariate_file, test_names_file, expression_file, genotype_file, sample_overlap_file, environmental_variable_file, num_bootstraps, sampling_fraction, output_file, method)
