args = commandArgs(trailingOnly=TRUE)
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




make_umap_loading_scatter_plot_colored_by_categorical_variable <- function(covariates, umap_loadings, covariate_name) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=factor(covariates))
	#valid_indices <- df$loading_1 < 10 & df$loading_2 < 10
	#df <- df[valid_indices, ]
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.001) +
	           gtex_v8_figure_theme() + 
	           guides(colour = guide_legend(override.aes = list(size=2))) +
	           labs(x="UMAP 1", y = "UMAP 2", color=covariate_name) + 
	           guides(colour=guide_legend(nrow=3,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.position="bottom") + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}


make_umap_loading_scatter_plot_colored_by_real_valued_variable <- function(covariates, umap_loadings, covariate_name) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=covariates)
	#valid_indices <- df$loading_1 < 10 & df$loading_2 < 10
	#df <- df[valid_indices, ]


	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.01) +
	           gtex_v8_figure_theme() + 
	           labs(x="UMAP 1", y = "UMAP 2", color=covariate_name) + 
	           scale_color_gradient(low="pink",high="blue") +
	           theme(legend.position="bottom") + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}

loading_scatter_plot_colored_by_categorical_covariate <- function(covariates, loadings, dim1, dim2, covariate_name) {
	df <- data.frame(loading_1=loadings[,dim1], loading_2=loadings[,dim2], covariate=factor(covariates))
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.001) +
	           gtex_v8_figure_theme() + 
	           labs(x="Loading 1", y = "Loading 2", color=covariate_name) + 
	           theme(legend.position="bottom") + 
	           guides(colour=guide_legend(nrow=3,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}





make_loading_boxplot_plot_by_categorical_covariate <- function(covariates, loadings, covariate_name) {
	loading_vec <- c()
	covariate_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		covariate_vec <- c(covariate_vec, as.character(covariates))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates)))
	}


	df <- data.frame(loading=loading_vec, covariate=factor(covariate_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))



	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=covariate)) + geom_boxplot(outlier.size = .00001) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill=covariate_name) +
	        	theme(legend.position="bottom") +
	        	guides(colour = guide_legend(override.aes = list(size=2))) +
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2)))
	   

	return(boxplot)
}

make_loading_boxplot_for_one_factor_by_categorical_covariate <- function(covariate, loadings, factor_number, covariate_name) {
	df <- data.frame(loading=loadings, covariate=factor(covariate))

	boxplot <- ggplot(df, aes(x=covariate, y=loading, fill=covariate)) + geom_boxplot(outlier.size = .00001) +
				gtex_v8_figure_theme() + 
	        	labs(x="", y = paste0("Sample loading (", factor_number,")"), fill="") +
	        	theme(legend.position="none") +
	        	guides(colour = guide_legend(override.aes = list(size=2))) +
	        	theme(axis.text.x=element_blank()) + 
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) +
	           	geom_hline(yintercept=0)


}


make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate <- function(covariate, loadings, covariate_name) {
	loading_vec <- c()
	covariate_vec <- c()
	num_factors <- dim(loadings)[2]
	print(num_factors)

	factor_number <- 1
	factor_1_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 2
	factor_2_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)


	factor_number <- 3
	factor_3_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 4
	factor_4_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 5
	factor_5_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)


	combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, factor_4_boxplot, factor_5_boxplot, ncol=1)
	#combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, factor_4_boxplot, ncol=1)

	return(combined)
}

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
	# Covariates columns to consider
	#print(summary(covariates[, valid_covariates]))
	# Remove unimportant columns
	covariates$experiment_day = factor(paste0(covariates$experiment, "_", covariates$day))

    loadings <- as.matrix(loadings)
    #valid_covariates <- 2:85
    #covs <- covariates[,valid_covariates]
    valid_covariates <- c(6, 7, 8, 10, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 36, 37, 38, 40, 41, 43, 44, 45,46, 47, 48, 49, 51, 52, 53, 54, 59, 60, 69, 71, 87, 95, 96, 97, 98)
 	covariate_type <- c("num", "cat", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "cat", "num", "num", "num", "num", "num", "num", "num")

 	cov_names <- colnames(covariates)[valid_covariates]
 	num <- length(cov_names)

 	print(cov_names)

  	covs <- covariates[,valid_covariates]
	#covs <- covariates


    # Initialize PVE heatmap
    factor_colnames <- paste0("Factor", 1:(dim(loadings)[2]))
    factor_rownames <- colnames(covs)
    pve_map <- matrix(0, dim(covs)[2] + 2, dim(loadings)[2])
    colnames(pve_map) <- factor_colnames
    rownames(pve_map) <- c(colnames(covs), "experiment+day", "experiment+pseudotime")


    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(loadings)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
    	pc_vec <- loadings[,num_pc]
        for (num_cov in 1:num_covs) {
            cov_vec <- covs[,num_cov]
            if (covariate_type[num_cov] == "cat") {
            #print(cov_vec[1:10])
            	lin_model <- lm(pc_vec ~ factor(cov_vec))
        	} else {
        		lin_model <- lm(pc_vec ~ cov_vec)
        	}
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
        lin_model <- lm(pc_vec ~ factor(covs[,4]) + covs[,2])
        pve_map[(num_covs+1), num_pc] <- summary(lin_model)$adj.r.squared

       	lin_model <- lm(pc_vec ~ factor(covs[,4]) + covs[,38])
        pve_map[(num_covs+2), num_pc] <- summary(lin_model)$adj.r.squared
    }
    
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "Loading","PVE")

    melted_mat$Covariate = factor(melted_mat$Covariate, levels=rownames(pve_map)[ord])
    melted_mat$Loading = factor(melted_mat$Loading, levels=factor_colnames)
	 #  Use factors to represent covariate and pc name
   	# melted_mat$Covariate 
    # melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    #melted_mat$PC <- substr(as.character(melted_mat$PC),3,5)
    #melted_mat$PC <- factor(melted_mat$PC, levels=paste0("", 1:(length(unique(melted_mat$PC)))))

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=Loading)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + labs(y="",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)

}


make_mixture_model_loading_proportion_plot <- function(loadings, covariates, title) {
	K <- dim(loadings)[2]
	cell_types <- c()
	factor_num <- c()
	counts <- c()

	ct_vec <- as.character(covariates$ct_cov)
	unique_cell_types <- as.character(unique(ct_vec))
	num_cell_types <- length(unique_cell_types)

	print(head(ct_vec[loadings[,1] > .5]))
	print(head(ct_vec[loadings$V1 > .5]))

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
    	labs(x="Factor num",y="Num samples",fill="", title=title) + 
    	gtex_v8_figure_theme() +
    	theme(legend.position="bottom")
    return(plotter)

}


make_factor_distribution_histogram_for_one_factor <- function(eqtl_factor, title, max_value) {
	df = data.frame(eqtl_factor=abs(eqtl_factor))
	histo <- ggplot(df, aes(x=eqtl_factor)) + geom_histogram(breaks=seq(0,max_value,.05)) +
			gtex_v8_figure_theme() +
			labs(x="Absolute factor value", title=title) 
	return(histo)

}

make_factor_distribution_histograms <- function(factors) {
	maxy = max(abs(factors))
	maxy = 3
	f_1_histo <- make_factor_distribution_histogram_for_one_factor(as.numeric(factors[1,]), paste0("Factor 1"), maxy)
	f_2_histo <- make_factor_distribution_histogram_for_one_factor(as.numeric(factors[2,]), paste0("Factor 2"), maxy)
	f_3_histo <- make_factor_distribution_histogram_for_one_factor(as.numeric(factors[3,]), paste0("Factor 3"), maxy)
	f_4_histo <- make_factor_distribution_histogram_for_one_factor(as.numeric(factors[4,]), paste0("Factor 4"), maxy)
	f_5_histo <- make_factor_distribution_histogram_for_one_factor(as.numeric(factors[5,]), paste0("Factor 5"), maxy)

	combined <- plot_grid(f_1_histo, f_2_histo, f_3_histo, f_4_histo, f_5_histo, ncol=1)

	return(combined)
}

make_factor_distribution_violins <- function(factors) {
	num_factors <- dim(factors)[1]
	
	factor_names <- c()
	factor_values <- c()
	for (factor_num in 1:num_factors) {
		factor_values <- c(factor_values, abs(as.numeric(factors[factor_num,])))
		factor_names <- c(factor_names, rep(paste0("factor",factor_num), length(as.numeric(factors[factor_num,]))))
	}
	df <- data.frame(factor_value=factor_values, factor_name=factor(factor_names))

	p <- ggplot(df, aes(factor_name, factor_value)) + geom_violin() +
		gtex_v8_figure_theme() +
		labs(y="Absolute value of Factor")
	return(p)

} 

make_factor_residual_variance_scatterplot <- function(factors, taus) {
	num_tests <- dim(factors)[2]
	max_factors <- c()
	for (test_num in 1:num_tests) {
		max_factors <- c(max_factors, max(abs(factors[,test_num])))
	}
	df <- data.frame(factor_value=max_factors, residual_precision=taus)

	p <- ggplot(df, aes(x=factor_value, y=residual_precision)) + geom_point() +
		gtex_v8_figure_theme() +
		labs(y="Residual Precision of test", x="Max factor value of test")
	return(p)

}

compute_pve_of_eqtl_factors <- function(loadings, factors, taus) {
	pve_vec <- c()
	s_k_vec <- c()

	K <- dim(loadings)[2]
	N <- dim(loadings)[1]
	T <- dim(factors)[2]
	print(T)

	for (k in 1:K) {
		aa = as.matrix(loadings[, k]) %*% as.matrix(factors[k,])
		s_k = sum(aa*aa)
		s_k_vec <- c(s_k_vec, s_k)
	}
	if (length(taus) ==1) {
		tau_term <- sum(1/taus)*N*T
	} else {
		tau_term <- sum(1/taus)*N
	}

	for (k in 1:K) {
		pve <- s_k_vec[k]/(sum(s_k_vec) + tau_term)
		pve_vec <- c(pve_vec, pve)
	}


	return(pve_vec)

}

make_pseudotime_factor_scatter <- function(pseudotime, loading, day) {
	df <- data.frame(loading=loading, pseudotime=pseudotime, day=factor(day))
	print(cor(pseudotime, loading))
	p <- ggplot(df, aes(x=loading, y=pseudotime, colour=day)) + geom_point(size=.0001) +
		gtex_v8_figure_theme() +
		labs(y="Pseudotime", x="Loading")
	return(p)
}

make_loading_covariate_scatter_for_each_factor <- function(expression_pc, loading, factor_number, x_axis_label) {
	df <- data.frame(expression_pc=expression_pc, loading=loading)
	print(factor_number)
	print(cor(expression_pc, loading)*cor(expression_pc, loading))

	plotter <- ggplot(df, aes(x=expression_pc, y=loading)) + 
	           geom_point(size=.001, alpha=.35) +
	           gtex_v8_figure_theme() + 
	           labs(x=x_axis_label, y = paste0("Loading ", factor_number)) + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) +
	           geom_smooth()
	return(plotter)

}


loading_covariate_scatter_with_row_for_every_factor <- function(expression_pc, loadings, x_axis_label) {
	factor_number <- 1
	factor_1_scatterplot <- make_loading_covariate_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)


	factor_number <- 2
	factor_2_scatterplot <- make_loading_covariate_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)


	factor_number <- 3
	factor_3_scatterplot <- make_loading_covariate_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)

	factor_number <- 4
	factor_4_scatterplot <- make_loading_covariate_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)

	factor_number <- 5
	ffactor_5_scatterplot <- make_loading_covariate_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)



	combined <- plot_grid(factor_1_scatterplot, factor_2_scatterplot, factor_3_scatterplot, factor_4_scatterplot, ffactor_5_scatterplot, ncol=1)


	return(combined)
}

make_loading_boxplot_for_one_factor_by_2_categorical_covariates <- function(covariate1, covariate2, loadings, factor_number, covariate_name) {
	df <- data.frame(loading=loadings, covariate1=factor(covariate1), covariate2=factor(covariate2))

	boxplot <- ggplot(df, aes(x=covariate1, y=loading, fill=covariate2)) + geom_boxplot(outlier.size = .00001) +
				gtex_v8_figure_theme() + 
	        	labs(x="", y = paste0("Sample loading (", factor_number,")"), fill="") +
	        	theme(legend.position="none") +
	        	guides(colour = guide_legend(override.aes = list(size=2))) +
	        	theme(axis.text.x=element_blank()) + 
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) +
	           	geom_hline(yintercept=0)

}



make_loading_boxplot_plot_with_row_for_every_factor_by_2_categorical_covariates <- function(covariate1, covariate2, loadings, covariate_name) {
	loading_vec <- c()
	covariate_vec <- c()
	num_factors <- dim(loadings)[2]
	print(num_factors)

	factor_number <- 1
	factor_1_boxplot <- make_loading_boxplot_for_one_factor_by_2_categorical_covariates(covariate1, covariate2, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 2
	factor_2_boxplot <- make_loading_boxplot_for_one_factor_by_2_categorical_covariates(covariate1, covariate2, loadings[,factor_number], factor_number, covariate_name)


	factor_number <- 3
	factor_3_boxplot <- make_loading_boxplot_for_one_factor_by_2_categorical_covariates(covariate1, covariate2, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 4
	factor_4_boxplot <- make_loading_boxplot_for_one_factor_by_2_categorical_covariates(covariate1, covariate2, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 5
	factor_5_boxplot <- make_loading_boxplot_for_one_factor_by_2_categorical_covariates(covariate1, covariate2, loadings[,factor_number], factor_number, covariate_name)

	#factor_number <- 4
	#factor_4_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	#factor_number <- 5
	#factor_5_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)


	#combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, factor_4_boxplot, factor_5_boxplot, ncol=1)
	combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, factor_4_boxplot, factor_5_boxplot, ncol=1)

	return(combined)
}


#############################
# Command line args
#############################
pre_processed_data_dir <- args[1]
eqtl_factorization_input_dir <- args[2]
eqtl_factorization_results_dir <- args[3]
eqtl_visualization_dir <- args[4]


eqtl_factorization_stem <- "/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/eqtl_factorization_results/eqtl_factorization_single_cell_sig_dynamic_eqtls_min_expressed_cells_5_factors_eqtl_factorization_vi_fixed_environmental_effect_learn_cov_model_False_re_False_svi_0_seed_1_lasso_param_temper_"
print(eqtl_factorization_stem)
#eqtl_factorization_stem <- "/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/eqtl_factorization_results/eqtl_factorization_single_cell_sig_tests_50_pc_min_expressed_cells_5_factors_eqtl_factorization_vi_with_re_model_False_re_False_svi_0_seed_100_lasso_param_temper_"
#eqtl_factorization_stem <- "/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/eqtl_factorization_results/eqtl_factorization_single_cell_sig_tests_50_pc_min_expressed_cells_5_factors_eqtl_factorization_vi_zero_inflated2_model_False_re_False_svi_0_seed_100_lasso_param_temper_"
output_stem = "sig_dynamic_eqtls_top_2000_re_model_spike_and_slab_and_environment_fixed_effect_learn_cov_20_"
#output_stem = "re_model_"
#output_stem = "no_re_model_"

eqtl_visualization_dir <- paste0(eqtl_visualization_dir, output_stem)



#genotype_file <- "/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data/eqtl_factorization_input/single_cell_sig_tests_in_each_day_standardized_genotype_training_data_uncorrected_r_squared_pruned.txt"
#G = read.table(genotype_file, header=FALSE)


# Input files
covariate_file <- paste0(pre_processed_data_dir, "cell_covariates_appeneded.txt")
# Load in data
covariates <- read.table(covariate_file, comment.char="*", header=TRUE, sep="\t")

model_covariates_file <- paste0(eqtl_factorization_input_dir, "single_cell_sig_covariate_modulated_eqtls_top_2000_covariate_subset_10.txt")
expression_pcs = read.table(model_covariates_file, header=FALSE, sep="\t")
print(summary(expression_pcs))

eqtl_factorization_loading_file <- paste0(eqtl_factorization_stem, "U_S.txt")
loadings <- read.table(eqtl_factorization_loading_file, header=FALSE)

eqtl_factorization_factor_file <- paste0(eqtl_factorization_stem, "V.txt")
factors <- read.table(eqtl_factorization_factor_file, header=FALSE)

eqtl_factorization_tau_file <- paste0(eqtl_factorization_stem, "tau.txt")
taus <- read.table(eqtl_factorization_tau_file, header=FALSE)$V1


######################################
# Scatterplot of pseudotime with 1st factor
#######################################
output_file <- paste0(eqtl_visualization_dir, "pseudotime_factor1_scatterplot.pdf")
#pseudotime_factor1_scatter <- make_pseudotime_factor_scatter(covariates$princ_curve_scaled01, loadings[,1], covariates$day)
#ggsave(pseudotime_factor1_scatter, file=output_file, width=7.2, height=6.0, units="in")


# Compute PVE of eqtl factors
pve_of_eqtl_factors <- compute_pve_of_eqtl_factors(loadings, factors, taus)
print(pve_of_eqtl_factors)

print("UMAP START")
umap_loadings = umap(loadings)$layout
saveRDS( umap_loadings, paste0(eqtl_visualization_dir, "umap_loadings.rds"))
#umap_loadings <- readRDS(paste0(eqtl_visualization_dir, "umap_loadings.rds"))
print("UMAP DONE")



######################################
# Make Scatter plot of max factor value for test and residual variance
#######################################
output_file <- paste0(eqtl_visualization_dir, "factor_value_residual_variance_scatterplot.pdf")
scatterplot <- make_factor_residual_variance_scatterplot(factors, taus) 
ggsave(scatterplot, file=output_file, width=7.2, height=5, units="in")

######################################
# Make Vionlin plot of factor values for each factor
#######################################
output_file <- paste0(eqtl_visualization_dir, "factor_distribution_violin.pdf")
violin_histogram <- make_factor_distribution_violins(factors) 
ggsave(violin_histogram, file=output_file, width=7.2, height=7, units="in")

######################################
# Make Histogram of factor values for each factor
#######################################
output_file <- paste0(eqtl_visualization_dir, "factor_distribution_histograms.pdf")
factor_histogram <- make_factor_distribution_histograms(factors) 
ggsave(factor_histogram, file=output_file, width=7.2, height=7, units="in")



######################################
# Make Heatmap correlating known covariates with Expression PCs
#######################################
output_file <- paste0(eqtl_visualization_dir, "expression_pc_covariate_heatmap.pdf")
#heatmap <- make_covariate_loading_correlation_heatmap(covariates, expression_pcs) 
#ggsave(heatmap, file=output_file, width=10.2, height=5.5, units="in")


######################################
# Make Heatmap correlating known covariates with eqtl factorization loadings
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_covariate_heatmap.pdf")
heatmap <- make_covariate_loading_correlation_heatmap(covariates, loadings) 
ggsave(heatmap, file=output_file, width=10.2, height=5.5, units="in")




######################################
# Make loading boxplot colored by cell type
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_boxplot_colored_by_day.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$day, loadings, "Differentiation day")
ggsave(boxplot, file=output_file, width=10.2, height=5.5, units="in")
print("done")


output_file <- paste0(eqtl_visualization_dir, "boxplot_colored_by_indi_id.pdf")
indi_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$donor_long_id, loadings, "Individual")
ggsave(indi_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")


output_file <- paste0(eqtl_visualization_dir, "boxplot_colored_by_experiment_day.pdf")
exp_day_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(factor(paste0(covariates$experiment,"_",covariates$day)), loadings, "Experiment_day")
ggsave(exp_day_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")


output_file <- paste0(eqtl_visualization_dir, "boxplot_colored_by_experiment.pdf")
experiment_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$experiment, loadings, "Experiment")
ggsave(experiment_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")

output_file <- paste0(eqtl_visualization_dir, "boxplot_colored_by_plate_id.pdf")
plate_id_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(factor(covariates$plate_id), loadings, "plate_id")
ggsave(plate_id_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")


output_file <- paste0(eqtl_visualization_dir, "boxplot_colored_by_experiment_day2.pdf")
exp_day_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_2_categorical_covariates(factor(paste0(covariates$experiment,"_",covariates$day)), factor(covariates$day), loadings, "Experiment_day")
ggsave(exp_day_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known percent mito
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_expr_PC1.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$princ_curve_scaled01, umap_loadings, "Expr PC1")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_differentiation_day.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$day, umap_loadings, "Differentiation day")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_indi_id.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$donor_long_id, umap_loadings, "Individual")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_donor.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$donor, umap_loadings, "Donor")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_log_total_counts.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$log10_total_counts, umap_loadings, "log_total_counts")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_pseudotime.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$pseudotime, umap_loadings, "Pseudotime")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_G2_M_transition.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$G2_M_transition, umap_loadings, "G2_M_transition")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_sterol_biosynthesis.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$sterol_biosynthesis, umap_loadings, "sterol_biosynthesis")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_G1_S_transition.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$G1_S_transition, umap_loadings, "G1_S_transition")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_differentiation.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$differentiation, umap_loadings, "differentiation")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_respiration.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$respiration, umap_loadings, "respiration")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_loading_1.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(loadings[,1], umap_loadings, "loading_1")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_loading_2.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(loadings[,2], umap_loadings, "loading_2")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_loading_3.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(loadings[,3], umap_loadings, "loading_3")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_loading_4.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(loadings[,4], umap_loadings, "loading_4")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_loading_5.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(loadings[,5], umap_loadings, "loading_5")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_well_id.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$well_id, umap_loadings, "Well ID")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_plate_id.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(factor(covariates$plate_id), umap_loadings, "Plate ID")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_experiment.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(factor(covariates$experiment), umap_loadings, "Experiment")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Make loading scatterplot with row for every factor colored by a real valued covariate
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_pseudotime_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_covariate_scatter_with_row_for_every_factor(covariates$pseudotime , loadings, "Pseudotime")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by a real valued covariate
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_G2_M_transition_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_covariate_scatter_with_row_for_every_factor(covariates$G2_M_transition , loadings, "G2_M_transition")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by a real valued covariate
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_sterol_biosynthesis_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_covariate_scatter_with_row_for_every_factor(covariates$sterol_biosynthesis , loadings, "sterol_biosynthesis")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")


######################################
# Make loading scatterplot with row for every factor colored by a real valued covariate
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_G1_S_transition_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_covariate_scatter_with_row_for_every_factor(covariates$G1_S_transition , loadings, "G1_S_transition")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by a real valued covariate
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_differentiation_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_covariate_scatter_with_row_for_every_factor(covariates$differentiation , loadings, "differentiation")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")


######################################
# Make loading scatterplot with row for every factor colored by a real valued covariate
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_respiration_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_covariate_scatter_with_row_for_every_factor(covariates$respiration , loadings, "respiration")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

