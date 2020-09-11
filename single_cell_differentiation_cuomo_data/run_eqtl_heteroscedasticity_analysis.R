args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)

get_dynamic_eqtl_lmm_pvalue <- function(expr, geno, covariates, environmental_variable, groups) {
	fit_full <- lmer(expr ~ geno + environmental_variable + covariates + environmental_variable:geno + (1|groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + environmental_variable + covariates + (1|groups), REML=FALSE)
	lrt <- anova(fit_null,fit_full)
	pvalue <- lrt[[8]][2]
	return(pvalue)
}

get_non_linear_dynamic_eqtl_lmm_pvalue <- function(expr, geno, covariates, environmental_variable, groups) {
	squared_environmental_variable = environmental_variable*environmental_variable
	fit_full <- lmer(expr ~ geno + environmental_variable + squared_environmental_variable + covariates + environmental_variable:geno + squared_environmental_variable:geno + (1|groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + environmental_variable + squared_environmental_variable + covariates + (1|groups), REML=FALSE)
	lrt <- anova(fit_null,fit_full)
	obj <- lrtest(fit_null, fit_full)
	pvalue <- obj[[5]][2]
	return(pvalue)
}

get_het_eqtl_lmm_pvalue <- function(expr, geno, covariates, groups) {
	fit_full <- lmer(expr ~ geno + covariates + (1|groups), REML=FALSE)
	fit_null <- lmer(expr ~ covariates + (1|groups), REML=FALSE)
	lrt <- anova(fit_null,fit_full)
	pvalue <- lrt[[8]][2]

	mod <- lmer(expr ~ geno + covariates + (1|groups))

	# print(summary(mod))

	vc <- VarCorr(mod)
	vc_df = as.data.frame(vc)
	random_effect_variance = vc_df$vcov[1]
	residual_variance = vc_df$vcov[2]
	random_effect_variance_fraction = random_effect_variance/(random_effect_variance + residual_variance)


	#p0 <- predict(mod)
	#p1 <- predict(mod, re.form=NA)

	#print(head(residuals(mod)))
	#print(head(expr-p0))
	#print(head(expr-p1))

	bp = bptest(residuals(mod) ~ geno)
	#bp = bptest((expr-p1) ~ geno)
	het_pvalue = bp$p.value

	return(list(het_pvalue=het_pvalue, random_effect_variance_fraction=random_effect_variance_fraction, eqtl_pvalue=pvalue))
}

get_maf <- function(geno) {
	af = sum(geno)/(2.0*length(geno))
	if (af > .5) {
		maf = 1.0-af
	} else {
		maf = af
	}
	return(maf)
}


#####################
# Command line args
#####################
bootstrapped_eqtl_stability_dir <- args[1]

################
# Input data files
covariate_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_covariate_subset_10.txt")
test_names_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_variant_gene_pairs.txt")
expression_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_expression.txt")
genotype_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_genotype.txt")
sample_overlap_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_individual_id.txt")
environmental_variable_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_environmental_variable.txt")


covariates <- as.matrix(read.table(covariate_file, header=FALSE))
environmental_variable = read.table(environmental_variable_file, header=FALSE)$V1
groups <- read.table(sample_overlap_file, header=FALSE)$V1


output_file <- paste0(bootstrapped_eqtl_stability_dir, "output2.txt")
sink(output_file)

# Stream files
stop = FALSE
count = 0

f_test = file(test_names_file, "r")
f_expr = file(expression_file, "r")
f_geno = file(genotype_file, "r")


line_test = readLines(f_test, n=1)
line_test = readLines(f_test, n=1)
line_expr = readLines(f_expr, n=1)
line_geno = readLines(f_geno, n=1)

while(!stop) {
	# Unpack data
	expr = as.numeric(strsplit(line_expr,'\t')[[1]])
	geno = as.numeric(strsplit(line_geno,'\t')[[1]])
	# Run eqtl analysis
	dynamic_eqtl_pvalue = get_dynamic_eqtl_lmm_pvalue(expr, geno, covariates, environmental_variable, groups)
	non_linear_dynamic_eqtl_pvalue = get_non_linear_dynamic_eqtl_lmm_pvalue(expr, geno, covariates, environmental_variable, groups)


	eqtl_info = get_het_eqtl_lmm_pvalue(expr, geno, covariates, groups)
	het_pvalue = eqtl_info$het_pvalue
	eqtl_pvalue = eqtl_info$eqtl_pvalue
	maf = get_maf(geno)
	random_effect_variance_fraction = eqtl_info$random_effect_variance_fraction
	#print(paste0(maf, "     ", dynamic_eqtl_pvalue, "     ", eqtl_pvalue, "     ", random_effect_variance_fraction, "     ", het_pvalue))

    new_line <- paste0(line_test, "\t", maf, "\t", dynamic_eqtl_pvalue, "\t", non_linear_dynamic_eqtl_pvalue, "\t", eqtl_pvalue, "\t", random_effect_variance_fraction, "\t", het_pvalue, "\n")
    cat(new_line)

	# Get info for next line
	count = count + 1
	line_test = readLines(f_test, n=1)
	line_expr = readLines(f_expr, n=1)
	line_geno = readLines(f_geno, n=1)
	if (length(line_test) == 0) {
		stop=TRUE
		close(f_test)
		close(f_expr)
		close(f_geno)
	}
}
sink()
