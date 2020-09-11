args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)
library(ggplot2)
library(reshape2)
library(cowplot)
options(bitmapType='cairo')

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}

get_dynamic_eqtl_lmm_pvalue <- function(expr, geno, covariates, groups, experiment_groups, plate_groups, pseudotime) {
	fit_full <- lmer(expr ~ geno + covariates + pseudotime + pseudotime:geno + (1|groups) + (1|plate_groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + pseudotime + (1|groups) + (1|plate_groups), REML=FALSE)
	lrt <- anova(fit_null,fit_full)
	pvalue <- lrt[[8]][2]

	return(pvalue)
}




get_eqtl_lmm_pvalue <- function(expr, geno, covariates, groups, experiment_groups) {
	fit_full <- lmer(expr ~ geno + covariates + (1|groups) + (1|experiment_groups), REML=FALSE)
	fit_null <- lmer(expr ~ covariates + (1|groups)+ (1|experiment_groups), REML=FALSE)
	lrt <- anova(fit_null,fit_full)
	pvalue <- lrt[[8]][2]
	return(pvalue)
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
total_lines <- as.numeric(args[2])
job_number = as.numeric(args[3])
num_jobs = as.numeric(args[4])


# Determine number of lines each parrallelized job will complete
lines_per_job = ceiling(total_lines/num_jobs)
start_num = job_number*lines_per_job
end_num = (job_number + 1)*lines_per_job


################
# Input data files
covariate_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_covariate_subset_10.txt")
test_names_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_variant_gene_pairs.txt")
expression_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_expression.txt")
genotype_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_genotype.txt")
sample_overlap_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_individual_id.txt")
experiment_overlap_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_experiment_id.txt")
plate_overlap_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_plate_id.txt")
pseudotime_file = paste0(bootstrapped_eqtl_stability_dir, "bootstrapped_eqtl_stability_input_data_cis_window_250000_pseudotime.txt")

covariates <- as.matrix(read.table(covariate_file, header=FALSE))
groups <- read.table(sample_overlap_file, header=FALSE)$V1
experiment_groups <- read.table(experiment_overlap_file, header=FALSE)$V1
plate_groups <- read.table(plate_overlap_file, header=FALSE)$V1
pseudotime <- read.table(pseudotime_file, header=FALSE)$V1


output_file <- paste0(bootstrapped_eqtl_stability_dir, "dynamic_eqtl_results_", job_number, "_", num_jobs, ".txt")
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
	if (count >= start_num & count < end_num) {
	# Unpack data
		expr = as.numeric(strsplit(line_expr,'\t')[[1]])
		geno = as.numeric(strsplit(line_geno,'\t')[[1]])
		# Run eqtl analysis
		line_info <- strsplit(line_test,'\t')[[1]]
		# output_fig <- paste0(bootstrapped_eqtl_stability_dir, "covariate_modulated_viz_", line_info[1], "_", line_info[2], ".pdf")
		dynamic_eqtl_pvalue = get_dynamic_eqtl_lmm_pvalue(expr, geno, covariates, groups, experiment_groups, plate_groups, pseudotime)
		eqtl_pvalue = get_eqtl_lmm_pvalue(expr, geno, covariates, groups, experiment_groups)
		maf = get_maf(geno)
		new_line <- paste0(line_test, "\t", maf, "\t", dynamic_eqtl_pvalue, "\t", eqtl_pvalue, "\n")
    	cat(new_line)
	}
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
