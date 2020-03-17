args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')


make_loading_scatter_plot_for_fixed_dimensions <- function(loadings, tissues, factor_1, factor_2, tissue_colors) {
	unique_tissues = unique(tissues)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}


	df <- data.frame(loading_1=loadings[,factor_1], loading_2=loadings[,factor_2], tissue=factor(tissues))
	plotter <- ggplot(df) + 
	           geom_point( aes(x=loading_1, y=loading_2, color=tissue),size=.01) +
	           gtex_v8_figure_theme() + 
	           scale_color_manual(values=colors) +
	           labs(x=paste0("Loading ", factor_1), y = paste0("Loading ", factor_2), color="") +
	           theme(legend.position="none")
	return(plotter)
}

make_loading_scatter_plot_for_fixed_dimensions_legend <- function(loadings, tissues, factor_1, factor_2, tissue_colors) {
	unique_tissues = unique(tissues)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}

	df <- data.frame(loading_1=loadings[,factor_1], loading_2=loadings[,factor_2], tissue=factor(tissues))
	plotter <- ggplot(df) + 
	           geom_point( aes(x=loading_1, y=loading_2, color=tissue),size=.01) +
	           gtex_v8_figure_theme() + 
	           scale_color_manual(values=colors) +
	           labs(x=paste0("Loading ", factor_1), y = paste0("Loading ", factor_2), color="") +
	           guides(colour=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.position="bottom")

	return(get_legend(plotter))
}

make_loading_scatter_plot <- function(tissues, sample_covariate_file, tissue_colors, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")
	race_factor = factor(as.character(as.numeric(covariates$race==3)))
	num_samples <- dim(loadings)[1]
	num_factors <- dim(loadings)[2]
	if (num_factors == 2) {
	unique_tissues = unique(tissues)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}
	df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues))
	plotter <- ggplot(df) + 
	           geom_point( aes(x=loading_1, y=loading_2, color=tissue),size=.005) +
	           gtex_v8_figure_theme() + 
	           scale_color_manual(values=colors) + 
	           labs(x="Loading 1", y = "Loading 2", color="")
	} else if (num_factors == 3) {
		plot_arr <- c()
		list_counter <- 1
		for (factor_1 in 1:num_factors) {
			for (factor_2 in 1:num_factors) {
				temp_plot <- make_loading_scatter_plot_for_fixed_dimensions(loadings, tissues, factor_1, factor_2, tissue_colors)
				if (factor_1 == factor_2) {
					plot_arr[[list_counter]] = NULL
				} else {
					plot_arr[[list_counter]] <- temp_plot
				}
				
				list_counter <- list_counter + 1

			}
		}

		legend <- make_loading_scatter_plot_for_fixed_dimensions_legend(loadings, tissues, factor_1, factor_2, tissue_colors)
		
		plotter <- plot_grid(plot_grid(plot_arr[[1]], plot_arr[[2]], plot_arr[[3]], plot_arr[[4]], plot_arr[[5]], plot_arr[[6]], plot_arr[[7]], plot_arr[[8]], ncol=num_factors), legend, ncol=1,rel_heights=c(1,.15))
	} else if (num_factors == 4) {
		plot_arr <- c()
		list_counter <- 1
		for (factor_1 in 1:num_factors) {
			for (factor_2 in 1:num_factors) {
				temp_plot <- make_loading_scatter_plot_for_fixed_dimensions(loadings, tissues, factor_1, factor_2, tissue_colors)
				if (factor_1 == factor_2) {
					plot_arr[[list_counter]] = NULL
				} else {
					plot_arr[[list_counter]] <- temp_plot
				}
				
				list_counter <- list_counter + 1

			}
		}

		legend <- make_loading_scatter_plot_for_fixed_dimensions_legend(loadings, tissues, factor_1, factor_2, tissue_colors)
		
		plotter <- plot_grid(plot_grid(plot_arr[[2]], plot_arr[[3]], plot_arr[[4]], plot_arr[[6]], plot_arr[[7]], plot_arr[[8]], NULL, plot_arr[[11]], plot_arr[[12]], ncol=num_factors-1), legend, ncol=1,rel_heights=c(1,.15))

	} else if (num_factors == 5) {
		plot_arr <- c()
		list_counter <- 1
		for (factor_1 in 1:num_factors) {
			for (factor_2 in 1:num_factors) {
				temp_plot <- make_loading_scatter_plot_for_fixed_dimensions(loadings, tissues, factor_1, factor_2, tissue_colors)
				if (factor_1 == factor_2) {
					plot_arr[[list_counter]] = NULL
				} else {
					plot_arr[[list_counter]] <- temp_plot
				}
				
				list_counter <- list_counter + 1

			}
		}

		legend <- make_loading_scatter_plot_for_fixed_dimensions_legend(loadings, tissues, factor_1, factor_2, tissue_colors)
		
		plotter <- plot_grid(plot_grid(plot_arr[[1]], plot_arr[[2]], plot_arr[[3]], plot_arr[[4]], plot_arr[[5]], plot_arr[[6]], plot_arr[[7]], plot_arr[[8]], plot_arr[[9]], plot_arr[[10]], plot_arr[[11]], plot_arr[[12]], plot_arr[[13]], plot_arr[[14]], plot_arr[[15]], plot_arr[[16]],plot_arr[[17]], plot_arr[[18]], plot_arr[[19]], plot_arr[[20]], plot_arr[[21]], plot_arr[[22]], plot_arr[[23]], plot_arr[[24]], ncol=num_factors), legend, ncol=1,rel_heights=c(1,.15))

	}

	return(plotter)
}

make_umap_loading_scatter_plot <- function(tissues, tissue_colors, sample_covariate_file, loading_file) {
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")
	loadings <- read.table(loading_file, header=FALSE)

	unique_tissues = unique(tissues)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}
	# 0.034993038124563364, -0.02809246824633377, 0.06527282210541974
	umap_loadings = umap(loadings)$layout

	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], tissue=factor(tissues), race=factor(as.character(as.numeric(covariates$race==3))))
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=tissue, shape=race)) +
	           scale_color_manual(values=colors) + 
	           scale_shape_manual(values = c(0,4)) +
	           gtex_v8_figure_theme() + 
	           guides(colour = guide_legend(override.aes = list(size=2))) +
	           labs(x="UMAP 1", y = "UMAP 2", color="", shape="Race") + 
	           guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.position="bottom")
	return(plotter)
}


make_loading_boxplot_plot_by_race <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	race_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		race_vec <- c(race_vec, as.character(as.numeric(covariates$race==3)))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$race)))
	}


	df <- data.frame(loading=loading_vec, race=factor(race_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=race)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known race") +
	        	theme(legend.position="bottom")

	return(boxplot)
}

make_loading_boxplot_plot_by_cm_proportion <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	race_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		race_vec <- c(race_vec, as.character(as.numeric(covariates$cm_cell_type_composition < .5)))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$cm_cell_type_composition)))
	}


	df <- data.frame(loading=loading_vec, race=factor(race_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=race)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known Cell type composition") +
	        	theme(legend.position="bottom")

	return(boxplot)
}

make_loading_boxplot_plot_by_sex <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	sex_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		sex_vec <- c(sex_vec, as.character(as.numeric(covariates$sex==2)))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$sex)))
	}


	df <- data.frame(loading=loading_vec, sex=factor(sex_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=sex)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known sex") +
	        	theme(legend.position="bottom")

	return(boxplot)
}

make_loading_boxplot_plot_by_cohort <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	cohort_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		cohort_vec <- c(cohort_vec, as.character(covariates$cohort=="Postmortem"))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$cohort)))
	}


	df <- data.frame(loading=loading_vec, cohort=factor(cohort_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=cohort)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known cohort") +
	        	theme(legend.position="bottom")

	return(boxplot)
}

make_loading_boxplot_plot_by_age <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	race_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		age_vec <- c(race_vec, covariates$age)
		print(factor_number)
		print(cor(loadings[,factor_number],age_vec))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$age)))
	}

	if (FALSE) {
	df <- data.frame(loading=loading_vec, sex=factor(sex_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=sex)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known sex") +
	        	theme(legend.position="bottom")

	return(boxplot)
	}
}

make_loading_boxplot_plot_by_tissue <- function(tissues,tissue_colors, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))

	loading_vec <- c()
	tissue_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		tissue_vec <- c(tissue_vec, as.character(tissues))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(tissues)))
	}


	df <- data.frame(loading=loading_vec, tissue=factor(tissue_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	unique_tissues = unique(df$tissue)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=tissue)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
				scale_fill_manual(values=colors) + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known tissue") +
	        	theme(legend.position="bottom") +
	        	guides(colour = guide_legend(override.aes = list(size=2))) +
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2)))
	   

	return(boxplot)
}


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

get_tissue_names <- function(sample_file_name) {
	aa <- read.table(sample_file_name)
	vecy = as.character(aa$V1)
	tissues <- c()
	for (iter in 1:length(vecy)) {
		tissue <- strsplit(vecy[iter],":")[[1]][2]
		tissues <- c(tissues, tissue)
	}
	return(tissues)
}

 make_residual_clustering_heatmap <- function(lm_residual_file, tissue_names) {
 	# Get matrix of dimension num_genesXnum_indi
 	residual_mat <- abs(t(read.table(lm_residual_file)))
 	#residual_mat <- residual_mat[1:100, 1:200]
 	residual_mat[residual_mat > 3] = 3.0
 	num_genes <- dim(residual_mat)[1]
 	num_indi <- dim(residual_mat)[2]
	colnames(residual_mat) = paste0("indi_", 1:num_indi)
	rownames(residual_mat) = paste0("gene_", 1:num_genes)

	ord <- hclust( dist(scale(residual_mat), method = "euclidean"), method = "ward.D" )$order

	melted_mat <- melt(residual_mat)
    colnames(melted_mat) <- c("gene", "sample", "error")

    #melted_mat$gene <- factor(paste0("gene_", 1:num_genes), levels = paste0("gene_", 1:num_genes)[ord])
    melted_mat$gene <- factor(melted_mat$gene, levels = paste0("gene_", 1:num_genes)[ord])

    melted_mat$sample <- factor(melted_mat$sample, levels = paste0("indi_", 1:num_indi))


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=gene, y=sample)) + geom_tile(aes(fill=error)) + 
		gtex_v8_figure_theme() +
   		labs(y="RNA-seq Sample", x="Variant-gene", fill="Absolute lm residual") +
   		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="pink",high="blue") +
   		scale_x_discrete(expand=c(0,0)) +
 		scale_y_discrete(expand=c(0,0))
	return(heatmap)
 }

  make_abs_expr_clustering_heatmap <- function(expr_file, tissue_names) {
  	print(expr_file)
 	# Get matrix of dimension num_genesXnum_indi
 	residual_mat <- abs(t(t(read.table(expr_file))))
 	residual_mat[residual_mat > 3] = 3.0
 	#residual_mat[residual_mat < -3] = -3.0
 	num_genes <- dim(residual_mat)[1]
 	num_indi <- dim(residual_mat)[2]


	colnames(residual_mat) = paste0("indi_", 1:num_indi)
	rownames(residual_mat) = paste0("gene_", 1:num_genes)


	ord <- hclust( dist(scale(residual_mat), method = "euclidean"), method = "ward.D" )$order
	melted_mat <- melt(residual_mat)
    colnames(melted_mat) <- c("gene", "sample", "error")

    #melted_mat$gene <- factor(paste0("gene_", 1:num_genes), levels = paste0("gene_", 1:num_genes)[ord])
    melted_mat$gene <- factor(melted_mat$gene, levels = paste0("gene_", 1:num_genes)[ord])

    melted_mat$sample <- factor(melted_mat$sample, levels = paste0("indi_", 1:num_indi))


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=gene, y=sample)) + geom_tile(aes(fill=error)) + 
		gtex_v8_figure_theme() +
   		labs(y="RNA-seq Sample", x="Variant-gene", fill="Abs Expr") +
   		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="pink",high="blue") +
   		scale_x_discrete(expand=c(0,0)) +
 		scale_y_discrete(expand=c(0,0))
	return(heatmap)
 }

 make_expr_clustering_heatmap <- function(expr_file, tissue_names) {
  	print(expr_file)
 	# Get matrix of dimension num_genesXnum_indi
 	residual_mat <- t(t(read.table(expr_file)))
 	residual_mat[residual_mat > 3] = 3.0
 	residual_mat[residual_mat < -3] = -3.0
 	num_genes <- dim(residual_mat)[1]
 	num_indi <- dim(residual_mat)[2]


	colnames(residual_mat) = paste0("indi_", 1:num_indi)
	rownames(residual_mat) = paste0("gene_", 1:num_genes)


	ord <- hclust( dist(scale(residual_mat), method = "euclidean"), method = "ward.D" )$order
	melted_mat <- melt(residual_mat)
    colnames(melted_mat) <- c("gene", "sample", "error")

    #melted_mat$gene <- factor(paste0("gene_", 1:num_genes), levels = paste0("gene_", 1:num_genes)[ord])
    melted_mat$gene <- factor(melted_mat$gene, levels = paste0("gene_", 1:num_genes)[ord])

    melted_mat$sample <- factor(melted_mat$sample, levels = paste0("indi_", 1:num_indi))


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=gene, y=sample)) + geom_tile(aes(fill=error)) + 
		gtex_v8_figure_theme() +
   		labs(y="RNA-seq Sample", x="Variant-gene", fill="Expr") +
   		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="pink",high="blue") +
   		scale_x_discrete(expand=c(0,0)) +
 		scale_y_discrete(expand=c(0,0))
	return(heatmap)
 }

elbo_line_plot_for_various_random_initializations <- function(eqtl_results_dir, model_stem, num_seeds, hot_start) {
	elbo_arr <- c()
	seed_id_arr <- c()
	iter_arr <- c()
	for (seed_number in 0:(num_seeds-1)) {
		seed_elbo_file <- paste0(eqtl_results_dir, model_stem , seed_number, "_seed_elbo.txt")
		temp_data <- read.table(seed_elbo_file, header=FALSE)
		elbo = temp_data[,1]
		if (hot_start == "True") {
			elbo = elbo[5:length(elbo)]
		}
		elbo_arr <- c(elbo_arr, elbo)
		seed_id_arr <- c(seed_id_arr, rep(paste0("seed ", seed_number), length(elbo)))
		iter_arr <- c(iter_arr, 1:length(elbo))
	}
	df <- data.frame(elbo=elbo_arr, iterations=iter_arr, seed=factor(seed_id_arr))


	p<-ggplot(df, aes(x=iterations, y=elbo, group=seed)) +
  		geom_line(aes(color=seed))+
		gtex_v8_figure_theme() +
		labs(x="Variational inference iteration", y="ELBO", color="")
	return(p)
}

make_absolute_effect_size_boxplot <- function(effect_size_file, tissue_colors, factor_file, num_tests) {
	effect_sizes <- read.table(effect_size_file, header=TRUE)
	all_factors <- read.table(factor_file, header=FALSE)
	num_tissues = dim(effect_sizes)[2]
	num_factors <- dim(all_factors)[1]

	effect_size_arr <- c()
	factor_num_arr <- c()
	tissue_name_arr <- c()

	unique_tissues = unique(colnames(effect_sizes))
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}

	for (factor_num in 1:num_factors) {
		# Find indices of top n factors
		sorted_factor <- sort(abs(all_factors[factor_num,]), decreasing=TRUE)
		min_val = as.numeric(sorted_factor[num_tests])
		indices = abs(all_factors[factor_num,]) >= min_val
	
		# Put info into data frame
		for (tissue_num in 1:num_tissues) {
			tissue_name <- colnames(effect_sizes)[tissue_num]
			tissue_effect_sizes <- effect_sizes[indices, tissue_num]
			effect_size_arr <- c(effect_size_arr, abs(tissue_effect_sizes))
			tissue_name_arr <- c(tissue_name_arr, rep(tissue_name, length(tissue_effect_sizes)))
			factor_num_arr <- c(factor_num_arr, rep(factor_num, length(tissue_effect_sizes)))
		}

	}
	df <- data.frame(effect_size=effect_size_arr, tissue=factor(tissue_name_arr), latent_factor=factor(factor_num_arr))
	boxplot <- ggplot(df, aes(x=latent_factor, y=effect_size, fill=tissue)) + geom_violin() +
				gtex_v8_figure_theme() + 
				scale_fill_manual(values=colors) + 
	        	labs(x="Latent factor", y = "Absolute effect size", fill="Known tissue") +
	        	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) +
	        	theme(legend.position="bottom")

	return(boxplot)

}


processed_data_dir <- args[1]
eqtl_results_dir <- args[2]
visualization_dir <- args[3]
tissue_colors_file <- args[4]


# Read in tissue colors and names
tissue_colors = read.table(tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")
# slight mislabeling
for (tiss_num in 1:length(tissue_colors$tissue_id)) {
	if (tissue_colors$tissue_id[tiss_num] == "Brain_Spinal_cord_cervical_c1") {
		tissue_colors$tissue_id[tiss_num] = "Brain_Spinal_cord_cervical_c.1"
	}
	if (tissue_colors$tissue_id[tiss_num] == "Cells_EBVtransformed_lymphocytes") {
		tissue_colors$tissue_id[tiss_num] = "Cells_EBV.transformed_lymphocytes"
	}
}

############################
# Model Specification
############################
model_name <- "eqtl_factorization_vi_spike_and_slab"
num_factors <- "25"
num_tissues <- "4"
random_effects_bool <- "False"
model_stem <- paste0("eqtl_factorization_tissues_subset_", num_tissues, "_gtex_data_", num_factors, "_factors_", model_name, "_model_", random_effects_bool, "_re_")
seed_number=5
seed_model_stem <- paste0(model_stem, seed_number, "_seed_")
loading_file <- paste0(eqtl_results_dir, seed_model_stem, "U_S.txt")

tissue_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_sample_names.txt")
sample_covariate_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_sample_covariates.txt")
effect_size_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_test_effect_sizes.txt")
tissue_names <- get_tissue_names(tissue_file)

######################
# Make box plot for each Race, showing loading distributions
output_file <- paste0(visualization_dir, seed_model_stem, "race_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_race(sample_covariate_file, loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each tissue, showing loading distributions
output_file <- paste0(visualization_dir, seed_model_stem, "tissue_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_tissue(tissue_names, tissue_colors, loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")
#####################
# Run Umap on loadings. Plot Umap loadings in scatter plot color by observed tissue type
output_file <- paste0(visualization_dir, seed_model_stem, "umap_loading_scatter.pdf")
umap_scatter <- make_umap_loading_scatter_plot(tissue_names, tissue_colors, sample_covariate_file, loading_file)
ggsave(umap_scatter, file=output_file, width=7.2*1.5, height=5.5*1.5, units="in")

######################
# Make scatter plot where each sample is a point, x and y axis are factor loadings, and points are colored by their tissue type
output_file <- paste0(visualization_dir, seed_model_stem, "loading_scatter.pdf")
scatter <- make_loading_scatter_plot(tissue_names,sample_covariate_file, tissue_colors, loading_file)
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")


if (FALSE) {
############################
# Model Specification
############################
model_name <- "vi_prior_on_loadings_only_special_init"
num_factors <- "15"
num_tissues <- "4"
random_effects_bool <- "False"
model_stem <- paste0("eqtl_factorization_tissues_subset_", num_tissues, "_gtex_data_", num_factors, "_factors_", model_name, "_model_", random_effects_bool, "_re_")
seed_number=1
seed_model_stem <- paste0(model_stem, seed_number, "_seed_")
loading_file <- paste0(eqtl_results_dir, seed_model_stem, "U_S.txt")

tissue_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_sample_names.txt")
sample_covariate_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_sample_covariates.txt")
effect_size_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_test_effect_sizes.txt")
tissue_names <- get_tissue_names(tissue_file)

######################
# Make box plot for each Race, showing loading distributions
output_file <- paste0(visualization_dir, seed_model_stem, "race_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_race(sample_covariate_file, loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each sex, showing loading distributions
output_file <- paste0(visualization_dir, seed_model_stem, "sex_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_sex(sample_covariate_file, loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each sex, showing loading distributions
output_file <- paste0(visualization_dir, seed_model_stem, "cohort_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_cohort(sample_covariate_file, loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each tissue, showing loading distributions
output_file <- paste0(visualization_dir, seed_model_stem, "tissue_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_tissue(tissue_names, tissue_colors, loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")
#####################
# Run Umap on loadings. Plot Umap loadings in scatter plot color by observed tissue type
output_file <- paste0(visualization_dir, seed_model_stem, "umap_loading_scatter.pdf")
umap_scatter <- make_umap_loading_scatter_plot(tissue_names, tissue_colors, sample_covariate_file, loading_file)
ggsave(umap_scatter, file=output_file, width=7.2*1.5, height=5.5*1.5, units="in")

}




if (FALSE) {
############################
# Model Specification
############################
model_name <- "vi_no_ard_learn_bernoulli"
num_factors <- "15"
num_tissues <- "4"
random_effects_bool <- "False"
model_stem <- paste0("eqtl_factorization_tissues_subset_", num_tissues, "_gtex_data_", num_factors, "_factors_", model_name, "_model_", random_effects_bool, "_re_")

############################
#Parameters
############################
num_seeds <- 6





#####################
# Make line plot showing ELBO over iterations for various random initializations
######################
hot_start="False"
output_file <- paste0(visualization_dir, model_stem, "_elbo_line_plot_hot_start_", hot_start, ".pdf")
elbo_line_plot <- elbo_line_plot_for_various_random_initializations(eqtl_results_dir, model_stem, num_seeds, hot_start)
ggsave(elbo_line_plot, file=output_file, width=12.2, height=5.5, units="in")

hot_start="True"
output_file <- paste0(visualization_dir, model_stem, "_elbo_line_plot_hot_start_", hot_start, ".pdf")
elbo_line_plot <- elbo_line_plot_for_various_random_initializations(eqtl_results_dir, model_stem, num_seeds, hot_start)
ggsave(elbo_line_plot, file=output_file, width=12.2, height=5.5, units="in")



######################
# Run remaining analysis with specified seed
######################
seed_number=5
seed_model_stem <- paste0(model_stem, seed_number, "_seed_")
loading_file <- paste0(eqtl_results_dir, seed_model_stem, "filtered_U_S.txt")
factor_file <- paste0(eqtl_results_dir, seed_model_stem, "filtered_V.txt")

tissue_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_sample_names.txt")
sample_covariate_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_sample_covariates.txt")
effect_size_file <- paste0(processed_data_dir, "tissues_subset_", num_tissues, "_test_effect_sizes.txt")
tissue_names <- get_tissue_names(tissue_file)



######################
# Make box plot showing absolute effect sizes in each tissue of genes loaded by each factor
num_tests <- 20
output_file <- paste0(visualization_dir, seed_model_stem, "absolute_effect_size_of_", num_tests, "_highest_loaded_factors_boxplot.pdf")
boxplot <- make_absolute_effect_size_boxplot(effect_size_file, tissue_colors, factor_file, num_tests)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each Race, showing loading distributions
output_file <- paste0(visualization_dir, seed_model_stem, "race_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_race(sample_covariate_file, loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each tissue, showing loading distributions
output_file <- paste0(visualization_dir, seed_model_stem, "tissue_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_tissue(tissue_names, tissue_colors, loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")
#####################
# Run Umap on loadings. Plot Umap loadings in scatter plot color by observed tissue type
output_file <- paste0(visualization_dir, seed_model_stem, "umap_loading_scatter.pdf")
umap_scatter <- make_umap_loading_scatter_plot(tissue_names, tissue_colors, sample_covariate_file, loading_file)
ggsave(umap_scatter, file=output_file, width=7.2*1.5, height=5.5*1.5, units="in")

######################
# Make scatter plot where each sample is a point, x and y axis are factor loadings, and points are colored by their tissue type
output_file <- paste0(visualization_dir, seed_model_stem, "loading_scatter.pdf")
scatter <- make_loading_scatter_plot(tissue_names,tissue_colors, loading_file)
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")

}

