args = commandArgs(trailingOnly=TRUE)
library(flexmix)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
#library(rhdf5)

options(bitmapType = 'cairo', device = 'pdf')

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


generate_eqtl_mixture_input <- function(G, Y) {
	num_tests = dim(G)[1]
	df = data.frame(G1=as.numeric(G[1,]), Y1=as.numeric(Y[1,]))
	for (test_num in 2:num_tests) {
		df[paste0("G", test_num)] = as.numeric(G[test_num,])
		df[paste0("Y", test_num)] = as.numeric(Y[test_num,])
	}
	return(df)
}

run_flexmix <- function(eqtl_mixture_input_df, num_tests, num_latent_factors) {
	model_list = list()
	print(num_latent_factors)
	for (test_num in 1:num_tests) {
		#model_list[[test_num]] = FLXMRglm(as.formula(paste0("Y", test_num, " ~ G", test_num)))
		model_list[[test_num]] = FLXMRglm(as.formula(paste0("Y", test_num, " ~ -1 + G", test_num)))
	}

	# Run model

	res <- flexmix(formula = ~ -1, data=eqtl_mixture_input_df, k=num_latent_factors, model=model_list)
	return(res)
}

make_mixture_model_loading_proportion_plot <- function(loadings, covariates) {
	K <- dim(loadings)[2]
	cell_types <- c()
	factor_num <- c()
	counts <- c()

	ct_vec <- as.character(covariates$day)
	unique_cell_types <- as.character(unique(ct_vec))
	num_cell_types <- length(unique_cell_types)


	for (k in 1:K) {
		#factor_indices = loadings[, k] > .5
		#print(length)

		for (cell_type_num in 1:num_cell_types) {
			cell_type <- as.character(unique_cell_types[cell_type_num])

			tot <- sum(ct_vec[loadings[,k] > .5] == cell_type)
			print(tot)

			cell_types <- c(cell_types, cell_type)
			factor_num <- c(factor_num, k)
			counts <- c(counts, tot)
		}
	}

	df <- data.frame(cell_type=factor(cell_types), factor_num=factor(factor_num), counts=counts)
	print(df)
	plotter <- ggplot(df,aes(x=factor_num, y=counts, fill=cell_type)) + 
    	geom_bar(stat="identity", position="fill") +
    	labs(x="Factor num",y="Num samples",fill="") + 
    	gtex_v8_figure_theme() +
    	theme(legend.position="bottom")
    return(plotter)

}

make_mixture_model_real_valued_covariate_stratefied_by_mixture <- function(loadings, covariate, covariate_name) {
	mixture_assignments <- c()
	num_samples <- dim(loadings)[1]
	K <- dim(loadings)[2]

	for (sample_num in 1:num_samples) {
		factor_int <- which.max(loadings[sample_num, ])
		mixture_assignments <- c(mixture_assignments, paste0("mixture_", factor_int))
	}
	df <- data.frame(mixture_assignments=factor(mixture_assignments), covariate=covariate)
	boxplot <- ggplot(df, aes(x=mixture_assignments, y=covariate)) + geom_boxplot() +
				gtex_v8_figure_theme() + 
	        	labs(x="", y = covariate_name, fill=covariate_name) +
	        	theme(legend.position="bottom") +
	        	guides(colour = guide_legend(override.aes = list(size=2))) +
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2)))
	return(boxplot)
}

genotype_training_file = args[1]
expression_training_file = args[2]
sample_overlap_file = args[3]
num_latent_factors = as.numeric(args[4])
cell_covariates_file = args[5]
output_root = args[6]
if (FALSE) {
print(num_latent_factors)
print(genotype_training_file)
print(expression_training_file)

G = read.table(genotype_training_file, header=FALSE)
print(dim(G))
Y = read.table(expression_training_file, header=FALSE)
num_tests <- dim(G)[1]
num_samples <- dim(G)[2]

#subsampled_indices = sample(1:num_samples, size=500, replace=FALSE)
#G <- G[, subsampled_indices]
#Y <- Y[, subsampled_indices]
#num_tests <- dim(G)[1]
#num_samples <- dim(G)[2]

#print(dim(Y))

G <- as.matrix(G)
Y <- as.matrix(Y)

eqtl_mixture_input_df <- generate_eqtl_mixture_input(G, Y)
print("START RUN")
flex_result <- run_flexmix(eqtl_mixture_input_df, num_tests, num_latent_factors)

output_file <- paste0(output_root, "flexmix_model.RDS")
saveRDS(flex_result, file=output_file)


output_file <- paste0(output_root, "flexmix_model.RDS")
flex_result <- readRDS(output_file)

posterior_probs <- posterior(flex_result)
}

output_file <- paste0(output_root, "flexmix_model_posterior_prob.txt")

posterior_probs <- read.table(output_file)
print(summary(posterior_probs))

if (FALSE) {
write.table(posterior_probs, file=output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

print(summary(posterior(flex_result)))
}


print(output_file)
print(cell_covariates_file)
covariates <- read.table(cell_covariates_file, comment.char="*", header=TRUE, sep="\t")


output_file <- paste0(output_root, "eqtl_mixture_loading_day_proportions.pdf")
plotter <- make_mixture_model_loading_proportion_plot(posterior_probs, covariates)
ggsave(plotter, file=output_file, width=7.2, height=6, units="in")

output_file <- paste0(output_root, "eqtl_mixture_pseudotime_boxplot.pdf")
plotter <- make_mixture_model_real_valued_covariate_stratefied_by_mixture(posterior_probs, covariates$pseudotime, "Pseudotime")
ggsave(plotter, file=output_file, width=7.2, height=6, units="in")
