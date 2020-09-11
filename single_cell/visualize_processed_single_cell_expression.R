library(ggplot2)
library(reshape2)
library(Matrix)
library(cowplot)
options(bitmapType='cairo')
options(warn=1)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}




make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors <- function(categorical_variable, dim_reduction_1, dim_reduction_2, categorical_label, x_axis_label, y_axis_label) {
	# Put everthing in compact data frame
	unique_categories = as.character(unique(categorical_variable))

	df <- data.frame(categorical_variable=factor(categorical_variable, levels=unique_categories), dim_reduction_1=dim_reduction_1, dim_reduction_2=dim_reduction_2)
	

	scatter <- ggplot(df, aes(x=dim_reduction_1, y=dim_reduction_2, color=categorical_variable)) +
  				geom_point(size=.01) +
  				figure_theme() + 
  				labs(x=x_axis_label,y=y_axis_label, color=categorical_label) +
  				guides(colour = guide_legend(override.aes = list(size=3)))
  	
  	
  	return(scatter)
}




make_pc_variance_explained_line_plot <- function(variance_explained, num_pcs) {
	variance_explained <- variance_explained[1:num_pcs]
	df <- data.frame(variance_explained = variance_explained, pc_num = 1:num_pcs)

	# PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=variance_explained)) +
                geom_line() +
                geom_point() +
                ylim(0,max(variance_explained) + .002) + 
                scale_x_continuous(breaks=seq(0,(num_pcs-1),5)) +
                labs(x = "PC number", y = "Variance Explained") + 
                figure_theme() 

    return(line_plot)
}




#########################
# Command line args
##########################
args <- commandArgs(TRUE)
processed_expression_dir <- args[1]  # Input dir
visualize_processed_expression_dir <- args[2]  # Output Dir



##########################
# Load in data
##########################
# Load in Covariates
filtered_covariate_file <- paste0(processed_expression_dir, "cell_covariates_sle_individuals_min_expressed_cells_0.05_log_transform_transform.txt")
filtered_covariate_data <- read.table(filtered_covariate_file, header=TRUE, sep="\t")

# Load in PCA Scores
filtered_pca_file <- paste0(processed_expression_dir, "pca_scores_sle_individuals_min_expressed_cells_0.05_log_transform_transform.txt")
filtered_pca_data <- read.table(filtered_pca_file, header=FALSE, sep="\t")

# Load in PCA PVE
pca_pve_file <- paste0(processed_expression_dir, "pca_variance_explained_sle_individuals_min_expressed_cells_0.05_log_transform_transform.txt")
pca_pve <- read.table(pca_pve_file, header=FALSE, sep="\t")




##########################
# Make UMAP Plot colored by cell type
##########################
#umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(covariate_data$ct_cov, umap_data[,1], umap_data[,2], "Cell Type", "UMAP1", "UMAP2", cell_type_colors_data)
#output_file <- paste0(visualize_processed_expression_dir, "umap_scatter_colored_by_cell_type.pdf")
#ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

##########################
# Make PCA Plot colored by cell type
##########################
pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$ct_cov, filtered_pca_data[,1], filtered_pca_data[,2], "Cell Type", "PC1", "PC2")
output_file <- paste0(visualize_processed_expression_dir, "sle_individuals_filtered_pca_1_2_scatter_colored_by_cell_type.pdf")
ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

#pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(covariate_data$ct_cov, pca_data[,1], pca_data[,2], "Cell Type", "PC1", "PC2", cell_type_colors_data)
#output_file <- paste0(visualize_processed_expression_dir, "pca_1_2_scatter_colored_by_cell_type.pdf")
#ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

#pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$ct_cov, filtered_pca_data[,2], filtered_pca_data[,3], "Cell Type", "PC2", "PC3", cell_type_colors_data)
#output_file <- paste0(visualize_processed_expression_dir, "sle_individuals_filtered_pca_2_3_scatter_colored_by_cell_type.pdf")
#ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

#pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(covariate_data$ct_cov, pca_data[,2], pca_data[,3], "Cell Type", "PC2", "PC3", cell_type_colors_data)
#output_file <- paste0(visualize_processed_expression_dir, "pca_2_3_scatter_colored_by_cell_type.pdf")
#ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

#pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$ct_cov, filtered_pca_data[,3], filtered_pca_data[,4], "Cell Type", "PC3", "PC4", cell_type_colors_data)
#output_file <- paste0(visualize_processed_expression_dir, "sle_individuals_filtered_pca_3_4_scatter_colored_by_cell_type.pdf")
#ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

#pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(covariate_data$ct_cov, pca_data[,3], pca_data[,4], "Cell Type", "PC3", "PC4", cell_type_colors_data)
#output_file <- paste0(visualize_processed_expression_dir, "pca_3_4_scatter_colored_by_cell_type.pdf")
#ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

##########################
# Make PCA PVE line plot
##########################
num_pcs <- 50
output_file <- paste0(visualize_processed_expression_dir, "sle_individuals_pca_variance_explained_", num_pcs, "_pcs_line_plot.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pca_pve[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")

num_pcs <- 100
#output_file <- paste0(visualize_processed_expression_dir, "sle_individuals_pca_variance_explained_", num_pcs, "_pcs_line_plot.pdf")
#ve_line_plot <- make_pc_variance_explained_line_plot(pca_pve[,1], num_pcs)
#ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")
