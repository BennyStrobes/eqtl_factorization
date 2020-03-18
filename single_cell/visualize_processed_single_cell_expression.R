library(ggplot2)
library(reshape2)
library(Matrix)
library(cowplot)
options(bitmapType='cairo')
options(warn=1)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}




make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors <- function(categorical_variable, dim_reduction_1, dim_reduction_2, categorical_label, x_axis_label, y_axis_label, color_data) {
	# Put everthing in compact data frame
	unique_categories = as.character(unique(categorical_variable))
	colors <- c()
	for (category_iter in 1:length(unique_categories)) {
		category <- unique_categories[category_iter]
		hex = color_data$hex_color[color_data$cell_type == category]
		colors <- c(colors, paste0(hex))
	}
	df <- data.frame(categorical_variable=factor(categorical_variable, levels=unique_categories), dim_reduction_1=dim_reduction_1, dim_reduction_2=dim_reduction_2)
	

	scatter <- ggplot(df, aes(x=dim_reduction_1, y=dim_reduction_2, color=categorical_variable)) +
  				geom_point(size=.01) +
  				figure_theme() + 
  				labs(x=x_axis_label,y=y_axis_label, color=categorical_label) +
  				scale_color_manual(values=colors) +
  				guides(colour = guide_legend(override.aes = list(size=3)))
  	
  	
  	return(scatter)
}









#########################
# Command line args
##########################
args <- commandArgs(TRUE)
processed_expression_dir <- args[1]  # Input dir
visualize_processed_expression_dir <- args[2]  # Output Dir


print(processed_expression_dir)
print(visualize_processed_expression_dir)

##########################
# Load in data
##########################
# Load in Covariates
covariate_file <- paste0(processed_expression_dir, "cell_covariates_ye_lab.txt")
covariate_data <- read.table(covariate_file, header=TRUE, sep="\t")
# Load in UMAP Scores
umap_file <- paste0(processed_expression_dir, "umap_scores_ye_lab.txt")
umap_data <- read.table(umap_file, header=FALSE, sep="\t")
# Load in PCA Scores
pca_file <- paste0(processed_expression_dir, "pca_scores_ye_lab.txt")
pca_data <- read.table(pca_file, header=FALSE, sep="\t")
# Load in Cell type colors
cell_type_colors_file <- paste0(processed_expression_dir, "cell_type_colors_ye_lab.txt")
cell_type_colors_data <- read.table(cell_type_colors_file, header=TRUE, sep="\t", comment.char="")



##########################
# Make UMAP Plot colored by cell type
##########################
umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(covariate_data$ct_cov, umap_data[,1], umap_data[,2], "Cell Type", "UMAP1", "UMAP2", cell_type_colors_data)
output_file <- paste0(visualize_processed_expression_dir, "umap_scatter_colored_by_cell_type.pdf")
ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

##########################
# Make PCA Plot colored by cell type
##########################
pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(covariate_data$ct_cov, pca_data[,1], pca_data[,2], "Cell Type", "PC1", "PC2", cell_type_colors_data)
output_file <- paste0(visualize_processed_expression_dir, "pca_1_2_scatter_colored_by_cell_type.pdf")
ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")


pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(covariate_data$ct_cov, pca_data[,2], pca_data[,3], "Cell Type", "PC2", "PC3", cell_type_colors_data)
output_file <- paste0(visualize_processed_expression_dir, "pca_2_3_scatter_colored_by_cell_type.pdf")
ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(covariate_data$ct_cov, pca_data[,3], pca_data[,4], "Cell Type", "PC3", "PC4", cell_type_colors_data)
output_file <- paste0(visualize_processed_expression_dir, "pca_3_4_scatter_colored_by_cell_type.pdf")
ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")


