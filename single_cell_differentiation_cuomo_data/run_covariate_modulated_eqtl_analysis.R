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

make_covariate_modulated_eqtl_plot <- function(expr, geno, covariates, groups, experiment_groups, plate_groups, covariate_num) {
	#covariate_num <- 1
	num_cells = length(covariates[,covariate_num])
	temp <- data.frame(cov=covariates[,covariate_num],geno=geno, expr=expr, quartile=rep(NA, num_cells))
	temp$quartile <- with(temp, cut(cov, 
                                breaks=quantile(cov, probs=seq(0,1, by=0.1), na.rm=TRUE), 
                                include.lowest=TRUE))

	p <- ggplot(temp, aes(x=quartile, y=expr, fill=factor(geno))) +
  		geom_boxplot() +
  		figure_theme() + 
  		labs(x=paste0("Covariate ", covariate_num),y="Expression", color="Genotype") 
  	return(p)
}

get_covariate_modulated_eqtl_lmm_pvalue <- function(expr, geno, covariates, groups, experiment_groups, plate_groups, output_fig) {
	#num_cov = dim(covariates)[2]
	#pvalue_arr <- c()
	#squared_covariates = covariates*covariates
	#fit_null <- lmer(expr ~ geno + covariates + squared_covariates + (1|groups) + (1|plate_groups), REML=FALSE)
	#for (cov_num in 1:num_cov) {
	#	fit_full <- lmer(expr ~ geno + covariates + squared_covariates + (covariates[,cov_num]):geno + (1|groups) + (1|plate_groups), REML=FALSE)
	#	lrt <- anova(fit_null,fit_full)
	#	pvalue <- lrt[[8]][2]
	#	pvalue_arr <- c(pvalue_arr, pvalue)
	#}
	#print(length(unique(paste0(groups, "_", experiment_groups))))
	plate_indi_groups = factor(paste0(groups, "_", plate_groups))
	experiment_indi_groups = factor(paste0(groups, "_", experiment_groups))

	#print(experiment_indi_groups)
	if (FALSE) {
	fit_full <- lmer(expr ~ geno + covariates + covariates:geno + (1|groups) + (1|plate_groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + (1|groups) + (1|plate_groups), REML=FALSE)
	t_stats = as.numeric(coef(summary(fit_full))[,"t value"])
	t_stats = t_stats[13:length(t_stats)]
	lrt <- anova(fit_null,fit_full)
	#print(fit_full)
	#print(fit_null)
	pvalue_2 <- lrt[[8]][2]

	fit_full <- lmer(expr ~ geno + covariates + covariates:geno + (1|groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + (1|groups), REML=FALSE)
	t_stats = as.numeric(coef(summary(fit_full))[,"t value"])
	t_stats = t_stats[13:length(t_stats)]
	lrt <- anova(fit_null,fit_full)
	#print(fit_full)
	#print(fit_null)
	pvalue_1 <- lrt[[8]][2]

	fit_full <- lmer(expr ~ geno + covariates + covariates:geno + (1|groups) + (1|plate_groups) + (1|experiment_groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + (1|groups) + (1|plate_groups) + (1|experiment_groups), REML=FALSE)
	t_stats = as.numeric(coef(summary(fit_full))[,"t value"])
	t_stats = t_stats[13:length(t_stats)]
	lrt <- anova(fit_null,fit_full)
	#print(fit_full)
	#print(fit_null)
	pvalue_3<- lrt[[8]][2]

	squared_covariates <- covariates^2

	fit_full <- lmer(expr ~ geno + covariates + squared_covariates + covariates:geno + (1|groups) + (1|plate_groups) + (1|experiment_groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + squared_covariates + (1|groups) + (1|plate_groups) + (1|experiment_groups), REML=FALSE)
	t_stats = as.numeric(coef(summary(fit_full))[,"t value"])
	t_stats = t_stats[13:length(t_stats)]
	lrt <- anova(fit_null,fit_full)
	#print(fit_full)
	#print(fit_null)
	pvalue_4<- lrt[[8]][2]

	fit_full <- lmer(expr ~ geno + covariates + squared_covariates + covariates:geno + (1|groups) + (1|plate_groups) + (1|experiment_groups) + (1|plate_indi_groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + squared_covariates + (1|groups) + (1|plate_groups) + (1|experiment_groups) + (1|plate_indi_groups), REML=FALSE)
	t_stats = as.numeric(coef(summary(fit_full))[,"t value"])
	t_stats = t_stats[13:length(t_stats)]
	lrt <- anova(fit_null,fit_full)
	#print(fit_full)
	#print(fit_null)
	pvalue_5<- lrt[[8]][2]
	}

	fit_full <- lmer(expr ~ geno + covariates + covariates:geno + (1|groups) + (1 | plate_groups) + (1 | plate_indi_groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + (1|groups) + (1 | plate_groups) + (1 | plate_indi_groups), REML=FALSE)
	t_stats = as.numeric(coef(summary(fit_full))[,"t value"])
	t_stats = t_stats[13:length(t_stats)]
	lrt <- anova(fit_null,fit_full)
	#print(fit_full)
	#print(fit_null)
	pvalue_6<- lrt[[8]][2]

	if (FALSE) {
	fit_full <- lmer(expr ~ geno + covariates + covariates:geno + (1|groups) + (1 + geno | plate_groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + (1|groups) + (1 + geno| plate_groups), REML=FALSE)
	t_stats = as.numeric(coef(summary(fit_full))[,"t value"])
	t_stats = t_stats[13:length(t_stats)]
	lrt <- anova(fit_null,fit_full)
	#print(fit_full)
	#print(fit_null)
	pvalue_7<- lrt[[8]][2]
	}

	#cat(paste0(pvalue_1, "\t", pvalue_2, "\t", pvalue_3, "\t", pvalue_4, "\t", pvalue_5, "\t", pvalue_6, "\t", pvalue_7,"\n"))
	#print(paste0(which.min(pvalue_arr), "\t", min(pvalue_arr), "\t", pvalue))

	pvalue <- pvalue_6
	print(pvalue)
	if (pvalue < 1e-5) {
		covariate_num <- which.max(abs(t_stats))
		plotter <- make_covariate_modulated_eqtl_plot(expr, geno, covariates, groups, experiment_groups, plate_groups, covariate_num)
		ggsave(plotter, file=output_fig, width=7.2, height=5, units="in")
	}


	#t_stat_string = paste(as.character(t_stats), collapse=",")
	#res <- list(pvalue=pvalue, t_stat_string=t_stat_string)
	#print(paste(as.character(pvalue_arr), collapse=","))
	#print((pvalue))
	#print(t_stats[13:length(t_stats)])
	#num_cov = dim(covariates)[2]
	#coefs <- c()
	#for (cov_num in 1:num_cov) {
#		breaks = cut(covariates[, cov_num], breaks=10)
	#	names = levels(breaks)
	#	summer = 0
	#	devs = c()
	#	for (name_iter in 1:length(names)) {
	#		name = names[name_iter]
	#		indices = breaks==name
	#		sdev = sd(geno[indices])
	#		devs <- c(devs, sdev)
	#	}
	#	deciles = 1:length(names)
	#	lin_mod <- lm(devs ~ deciles)
		#coefs <- c(coefs, )
	#	coefs <- c(coefs, as.numeric(coef(summary(lin_mod))[,"t value"])[2])
	#}
	#print(pvalue)
	#print(t_stats)
	#print(coefs)
	#print(cor(abs(t_stats), abs(coefs)))


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

get_eqtl_lmm_pvalue <- function(expr, geno, covariates, groups) {
	fit_full <- lmer(expr ~ geno + covariates + (1|groups), REML=FALSE)
	fit_null <- lmer(expr ~ covariates + (1|groups), REML=FALSE)
	lrt <- anova(fit_null,fit_full)
	pvalue <- lrt[[8]][2]
	return(list(eqtl_pvalue=pvalue))
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



covariates <- as.matrix(read.table(covariate_file, header=FALSE))
groups <- read.table(sample_overlap_file, header=FALSE)$V1
experiment_groups <- read.table(experiment_overlap_file, header=FALSE)$V1
plate_groups <- read.table(plate_overlap_file, header=FALSE)$V1


output_file <- paste0(bootstrapped_eqtl_stability_dir, "covariate_modulated_eqtl_results3_", job_number, "_", num_jobs, ".txt")
#sink(output_file)

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
		output_fig <- paste0(bootstrapped_eqtl_stability_dir, "covariate_modulated_viz_", line_info[1], "_", line_info[2], ".pdf")
		covariate_modulated_eqtl_info = get_covariate_modulated_eqtl_lmm_pvalue(expr, geno, covariates, groups, experiment_groups, plate_groups, output_fig)
		eqtl_pvalue = get_eqtl_lmm_pvalue(expr, geno, covariates, groups)
		maf = get_maf(geno)
		#new_line <- paste0(line_test, "\t", maf, "\t", covariate_modulated_eqtl_info$pvalue, "\t", eqtl_pvalue, "\t", covariate_modulated_eqtl_info$t_stat_string, "\n")
    	#cat(new_line)
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
#sink()
