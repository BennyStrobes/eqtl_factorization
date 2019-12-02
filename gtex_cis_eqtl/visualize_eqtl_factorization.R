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

make_loading_scatter_plot <- function(tissues, tissue_colors, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
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
		
		plotter <- plot_grid(plot_grid(plot_arr[[1]], plot_arr[[2]], plot_arr[[3]], plot_arr[[4]], plot_arr[[5]], plot_arr[[6]], plot_arr[[7]], plot_arr[[8]], plot_arr[[9]], plot_arr[[10]], plot_arr[[11]], plot_arr[[12]], plot_arr[[13]], plot_arr[[14]], plot_arr[[15]], ncol=num_factors), legend, ncol=1,rel_heights=c(1,.15))

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

make_umap_loading_scatter_plot <- function(tissues, tissue_colors, loading_file) {
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

	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], tissue=factor(tissues))
	plotter <- ggplot(df) + 
	           geom_point( aes(x=loading_1, y=loading_2, color=tissue),size=.01) +
	           scale_color_manual(values=colors) + 
	           gtex_v8_figure_theme() + 
	           guides(colour = guide_legend(override.aes = list(size=2))) +
	           labs(x="UMAP 1", y = "UMAP 2", color="") 
	return(plotter)
}


make_loading_boxplot_plot <- function(tissues,tissue_colors, loading_file) {
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
				gtex_v8_figure_theme() + ylim(0,5) +
				scale_fill_manual(values=colors) + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known tissue") +
	        	theme(legend.position="bottom")

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



num_tissue=4
num_factors=4
seed=0
loading_file <- paste0(eqtl_results_dir, "eqtl_factorization_tissues_subset_", num_tissue, "_gtex_data_4_factors_vi_model_",seed,"_seed_U_S.txt")
tissue_file <- paste0(processed_data_dir, "tissues_subset_", num_tissue, "_sample_names.txt")
tissue_names <- get_tissue_names(tissue_file)

######################
# Make box plot for each tissue, showing loading distributions
output_file <- paste0(visualization_dir,"eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_vi_model_",seed,"_seed_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot(tissue_names, tissue_colors, loading_file)
ggsave(boxplot, file=output_file, width=12.2, height=5.5, units="in")
#####################
# Run Umap on loadings. Plot Umap loadings in scatter plot color by observed tissue type
output_file <- paste0(visualization_dir,"eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_vi_model_", seed, "_seed_umap_loading_scatter.pdf")
umap_scatter <- make_umap_loading_scatter_plot(tissue_names, tissue_colors, loading_file)
ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make scatter plot where each sample is a point, x and y axis are factor loadings, and points are colored by their tissue type
output_file <- paste0(visualization_dir, "eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_vi_model_", seed, "_seed_loading_scatter.pdf")
scatter <- make_loading_scatter_plot(tissue_names,tissue_colors, loading_file)
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")


if (FALSE) {
lasso_param_us = c("0.0001")
initializations = c("random_0")
num_factor_arr = c(4)
num_tissues = c(4)
for (lasso_param_u_iter in 1:length(lasso_param_us)) {
	for (initialization_iter in 1:length(initializations)) {
		for (num_factor_iter in 1:length(num_factor_arr)) {
			for (num_tissue_iter in 1:length(num_tissues)) {

				lasso_param_u <- lasso_param_us[lasso_param_u_iter]
				lasso_param_v <-  lasso_param_us[lasso_param_u_iter]
				initialization <- initializations[initialization_iter]
				num_factors <- num_factor_arr[num_factor_iter]
				num_tissue <- num_tissues[num_tissue_iter]
				loading_file <- paste0(eqtl_results_dir, "eqtl_factorization_tissues_subset_", num_tissue, "_gtex_data_", num_factors, "_factors_alm_lasso_U_", lasso_param_u, "_lasso_V_",lasso_param_v, "_initialization_", initialization, "_U.txt")

				tissue_file <- paste0(processed_data_dir, "tissues_subset_", num_tissue, "_sample_names.txt")
				tissue_names <- get_tissue_names(tissue_file)


				######################
				# Make heatmap showing residual clustering (heatmap of dimension num_genesXnum_samples)
				if (initialization == "residual_clustering") {
					# residual clustering data
					lm_residual_file <- paste0(eqtl_results_dir, "eqtl_factorization_tissues_subset_", num_tissue, "_gtex_data_", num_factors, "_factors_em_model_lasso_U_", lasso_param_u, "_lasso_V_", lasso_param_v, "_initialization_", initialization, "_lm_residuals.txt")
					output_file <- paste0(visualization_dir,"eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_em_model_lasso_U_", lasso_param_u, "_lasso_V_",lasso_param_v, "_initialization_", initialization, "_residual_clustering_heatmap.pdf")
					#heatmap <- make_residual_clustering_heatmap(lm_residual_file, tissue_names)
					#ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")

					expr_file <- paste0(processed_data_dir, "tissues_subset_", num_tissue, "_expr2.txt")
					output_file <- paste0(visualization_dir,"eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_em_model_lasso_U_", lasso_param_u, "_lasso_V_",lasso_param_v, "_initialization_", initialization, "_abs_expr_clustering_heatmap.pdf")
					#heatmap <- make_abs_expr_clustering_heatmap(expr_file, tissue_names)
					#ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")

					expr_file <- paste0(processed_data_dir, "tissues_subset_", num_tissue, "_expr2.txt")
					output_file <- paste0(visualization_dir,"eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_em_model_lasso_U_", lasso_param_u, "_lasso_V_",lasso_param_v, "_initialization_", initialization, "_expr_clustering_heatmap.pdf")
					#heatmap <- make_expr_clustering_heatmap(expr_file, tissue_names)
					#ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")


				}

				######################
				# Make box plot for each tissue, showing loading distributions
				output_file <- paste0(visualization_dir,"eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_alm_lasso_U_", lasso_param_u, "_lasso_V_",lasso_param_v, "_initialization_", initialization, "_loading_boxplot.pdf")
				boxplot <- make_loading_boxplot_plot(tissue_names, tissue_colors, loading_file)
				ggsave(boxplot, file=output_file, width=12.2, height=5.5, units="in")
				#####################
				# Run Umap on loadings. Plot Umap loadings in scatter plot color by observed tissue type
				output_file <- paste0(visualization_dir,"eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_alm_lasso_U_", lasso_param_u, "_lasso_V_",lasso_param_v, "_initialization_", initialization, "_umap_loading_scatter.pdf")
				umap_scatter <- make_umap_loading_scatter_plot(tissue_names, tissue_colors, loading_file)
				ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")

				######################
				# Make scatter plot where each sample is a point, x and y axis are factor loadings, and points are colored by their tissue type
				output_file <- paste0(visualization_dir, "eqtl_factorization_of_", num_tissue, "_tissues_with_", num_factors, "_factors_alm_lasso_U_", lasso_param_u, "_lasso_V_",lasso_param_v, "_initialization_", initialization, "_loading_scatter.pdf")
				scatter <- make_loading_scatter_plot(tissue_names,tissue_colors, loading_file)
				ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")


			}
		}
	}

}




##### KMEANS
output_file <- paste0(visualization_dir, "kmeans_expression_init_boxplot.pdf")
#boxplot <- make_loading_boxplot_plot(tissue_file, paste0(eqtl_results_dir, "initialization_of_kmeans_on_expression.txt"))
#ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")
#loading_file <- paste0(eqtl_results_dir, "eqtl_factorization_3_factors_em_model_lasso_U_0.01_lasso_V_0.0_initialization_random_genotype_intercept_True_U.txt")
#print(loading_file)


######################
# Make box plot for each tissue, showing loading distributions
output_file <- paste0(visualization_dir, "loading_boxplot for each tissue.pdf")
#boxplot <- make_loading_boxplot_plot(tissue_file, loading_file)
#ggsave(boxplot, file=output_file, width=12.2, height=5.5, units="in")

######################
# Make scatter plot where each sample is a point, x and y axis are factor loadings, and points are colored by their tissue type
#output_file <- paste0(visualization_dir, "loading_scatter_colored_by_tissue_type.pdf")
#scatter <- make_loading_scatter_plot(tissue_file, loading_file)
#ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")
}