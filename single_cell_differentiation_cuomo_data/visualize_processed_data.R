library(ggplot2)
library(reshape2)
library(Matrix)
library(cowplot)
options(bitmapType='cairo')
options(warn=1)

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}


make_dimensionality_reduction_scatter_colored_by_categorical_variable <- function(categorical_variable, dim_reduction_1, dim_reduction_2, categorical_label, x_axis_label, y_axis_label) {

	df <- data.frame(categorical_variable=factor(categorical_variable), dim_reduction_1=dim_reduction_1, dim_reduction_2=dim_reduction_2)
	

	scatter <- ggplot(df, aes(x=dim_reduction_1, y=dim_reduction_2, color=categorical_variable)) +
  				geom_point(size=.01) +
  				figure_theme() + 
  				labs(x=x_axis_label,y=y_axis_label, color=categorical_label) +
  				guides(colour = guide_legend(override.aes = list(size=3))) +
          theme(legend.position="none")
  	
  	
  	return(scatter)
}

make_gene_feature_scatterplot <- function(variance, fraction_of_zeros) {
	df <- data.frame(variance=variance, fraction_of_zeros=fraction_of_zeros)
	scatter <- ggplot(df, aes(x=variance, y=fraction_of_zeros)) +
  				geom_point(size=.01) +
  				figure_theme() + 
  				labs(x="Gene Variance",y="Fraction of zeros in each gene") +
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


make_fraction_of_zeros_histogram <- function(fraction_of_zeros) {
	df <- data.frame(fraction_of_zeros)
	histo <- ggplot(df, aes(x=fraction_of_zeros))+
  			geom_histogram(color="darkblue", fill="lightblue") + 
  			labs(x="Fraction of zeros / Gene") + 
  			figure_theme()
  	return(histo)
}

make_num_cells_per_plate_histogram <- function(covariate_data) {
  unique_plates = unique(covariate_data$plate_id)
  num_cells = c()
  for (plate_iter in 1:length(unique_plates)) {
    plate_name = unique_plates[plate_iter]
    cell_indices = covariate_data$plate_id == plate_name
    total = sum(cell_indices)
    num_cells <- c(num_cells, total)
  }

  df <- data.frame(num_cells)
  histo <- ggplot(df, aes(x=num_cells))+
        geom_histogram(color="darkblue", fill="lightblue") + 
        labs(x="Num cells / plate") + 
        figure_theme()
  return(histo)
}


make_num_experiments_per_plate_histogram <- function(covariate_data) {
  unique_plates = unique(covariate_data$plate_id)
  num_experiments = c()
  for (plate_iter in 1:length(unique_plates)) {
    plate_name = unique_plates[plate_iter]
    cell_indices = covariate_data$plate_id == plate_name
    total <- length(unique(covariate_data$experiment[cell_indices]))
    num_experiments <- c(num_experiments, total)
  }
  df <- data.frame(num_experiments)
  histo <- ggplot(df, aes(x=num_experiments))+
        geom_histogram(color="darkblue", fill="lightblue") + 
        labs(x="Num experiments / plate") + 
        figure_theme()
  return(histo)
}


make_num_days_per_plate_histogram <- function(covariate_data) {
  unique_plates = unique(covariate_data$plate_id)
  num_days = c()
  for (plate_iter in 1:length(unique_plates)) {
    plate_name = unique_plates[plate_iter]
    cell_indices = covariate_data$plate_id == plate_name
    total <- length(unique(covariate_data$day[cell_indices]))
    num_days <- c(num_days, total)
  }
  df <- data.frame(num_days)
  histo <- ggplot(df, aes(x=num_days))+
        geom_histogram(color="darkblue", fill="lightblue") + 
        labs(x="Num days / plate") + 
        figure_theme()
  return(histo)
}

make_num_plates_per_experiment_histogram <- function(covariate_data) {
  unique_experiments = unique(covariate_data$experiment)
  num_plates = c()
  for (experiment_iter in 1:length(unique_experiments)) {
    experiment_name = unique_experiments[experiment_iter]
    cell_indices = covariate_data$experiment == experiment_name
    total <- length(unique(covariate_data$plate_id[cell_indices]))
    num_plates <- c(num_plates, total)
  }
  df <- data.frame(num_plates)
  print(num_plates)
  histo <- ggplot(df, aes(x=num_plates))+
        geom_histogram(color="darkblue", fill="lightblue") + 
        labs(x="Num plates / experiment") + 
        figure_theme()
  return(histo)
}


debug <- function(covariate_data) {
  unique_plates = unique(covariate_data$plate_id)
  new_names <- paste0(covariate_data$experiment, " ", covariate_data$day)
  print(length(unique(new_names)))
  unique_experiments = unique(covariate_data$experiment)
  total = c()
  for (experiment_iter in 1:length(unique_experiments)) {
    experiment_name = unique_experiments[experiment_iter]
    for (day in 0:3) {
      cell_indices = (covariate_data$day == paste0("day",day)) & (covariate_data$experiment == experiment_name)
      total <- c(total, sum(cell_indices))

      if (length(unique(covariate_data$plate_id[cell_indices])) > 1) {
        print(summary(covariate_data[cell_indices,]))
      }
    }
  }

}


#########################
# Command line args
##########################
args <- commandArgs(TRUE)
processed_data_dir <- args[1]  # Input dir
visualize_processed_data_dir <- args[2]  # Output Dir




##########################
# Load in data
##########################
# Load in Covariates
cell_covariate_file <- paste0(processed_data_dir, "cell_covariates.txt")
covariate_data <- read.table(cell_covariate_file, sep="\t", header=TRUE, comment.char = "*")

# Load in PCS
pca_loadings_top_500_vg_file <- paste0(processed_data_dir, "standardized_normalized_top_500_variable_genes_expression_pca_loadings.txt")
pca_loadings_top_500_vg <- read.table(pca_loadings_top_500_vg_file, header=FALSE, sep="\t")

pca_loadings_top_2000_vg_file <- paste0(processed_data_dir, "standardized_normalized_top_2000_variable_genes_expression_pca_loadings.txt")
pca_loadings_top_2000_vg <- read.table(pca_loadings_top_2000_vg_file, header=FALSE, sep="\t")

pca_loadings_file <- paste0(processed_data_dir, "standardized_normalized_expression_pca_loadings.txt")
pca_loadings <- read.table(pca_loadings_file, header=FALSE, sep="\t")

# Load in PCA PVE
pca_pve_file <- paste0(processed_data_dir, "standardized_normalized_expression_pca_pve.txt")
pca_pve <- read.table(pca_pve_file, header=FALSE, sep="\t")

pca_pve_top_500_vg_file <- paste0(processed_data_dir, "standardized_normalized_top_500_variable_genes_expression_pca_pve.txt")
pca_pve_top_500_vg <- read.table(pca_pve_top_500_vg_file, header=FALSE, sep="\t")

pca_pve_top_2000_vg_file <- paste0(processed_data_dir, "standardized_normalized_top_2000_variable_genes_expression_pca_pve.txt")
pca_pve_top_2000_vg <- read.table(pca_pve_top_2000_vg_file, header=FALSE, sep="\t")


# Load in file containing variance of each gene
gene_variance_file <- paste0(processed_data_dir, "variance_of_each_gene.txt")
gene_variance <- read.table(gene_variance_file, header=TRUE, sep="\t")

# Load in file fraction of cells expressed in each gene
fraction_of_cells_expressed_file <- paste0(processed_data_dir, "fraction_of_zeros_in_each_gene.txt")
fraction_of_cells_expressed <- read.table(fraction_of_cells_expressed_file, header=TRUE, sep="\t")


debug(covariate_data)

##########################
# Make histogram showing number of plates per experiment
##########################
num_plates_per_experiment_histogram <- make_num_plates_per_experiment_histogram(covariate_data)
output_file <- paste0(visualize_processed_data_dir, "num_plates_per_experiment_histogram.pdf")
ggsave(num_plates_per_experiment_histogram, file=output_file, width=7.2, height=5, units="in")

##########################
# Make histogram showing number of days per plate
##########################
num_days_per_plate_histogram <- make_num_days_per_plate_histogram(covariate_data)
output_file <- paste0(visualize_processed_data_dir, "num_days_per_plate_histogram.pdf")
ggsave(num_days_per_plate_histogram, file=output_file, width=7.2, height=5, units="in")

##########################
# Make histogram showing number of experiments per plate
##########################
num_experiments_per_plate_histogram <- make_num_experiments_per_plate_histogram(covariate_data)
output_file <- paste0(visualize_processed_data_dir, "num_experiments_per_plate_histogram.pdf")
ggsave(num_experiments_per_plate_histogram, file=output_file, width=7.2, height=5, units="in")

##########################
# Make histogram showing number of cells per plate
##########################
num_cells_per_plate_histogram <- make_num_cells_per_plate_histogram(covariate_data)
output_file <- paste0(visualize_processed_data_dir, "num_cells_per_plate_histogram.pdf")
ggsave(num_cells_per_plate_histogram, file=output_file, width=7.2, height=5, units="in")


if (FALSE) {
##########################
# Make Scatter plot correlating gene variance with fraction of cells expressed by gene
##########################
gene_feature_scatter_plot <- make_gene_feature_scatterplot(gene_variance$variance, fraction_of_cells_expressed$fraction_of_zeros)
output_file <- paste0(visualize_processed_data_dir, "scatter_correlating_gene_variance_with_fraction_of_cells_expressed.pdf")
ggsave(gene_feature_scatter_plot, file=output_file, width=7.2, height=5, units="in")


##########################
# Make Histogram showing fraction of zeros in each gene
##########################
fraction_of_zeros_histogram <- make_fraction_of_zeros_histogram(fraction_of_cells_expressed$fraction_of_zeros)
output_file <- paste0(visualize_processed_data_dir, "fraction_of_zeros_histogram.pdf")
ggsave(fraction_of_zeros_histogram, file=output_file, width=7.2, height=5, units="in")


##########################
# Make PCA Plot colored by differentiation
##########################
pc_scatter_colored_by_day <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(covariate_data$day, pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_differentiation_day.pdf")
ggsave(pc_scatter_colored_by_day, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_plate_id <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$plate_id), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_plate_id.pdf")
ggsave(pc_scatter_colored_by_plate_id, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_experiment <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$experiment), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_experiment.pdf")
ggsave(pc_scatter_colored_by_experiment, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_donor <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$donor_long_id), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_donor.pdf")
ggsave(pc_scatter_colored_by_donor, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_day <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(covariate_data$day, pca_loadings_top_500_vg[,1], pca_loadings_top_500_vg[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_top_500_vg_colored_by_differentiation_day.pdf")
ggsave(pc_scatter_colored_by_day, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_day <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(covariate_data$day, pca_loadings_top_2000_vg[,1], pca_loadings_top_2000_vg[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_top_2000_vg_colored_by_differentiation_day.pdf")
ggsave(pc_scatter_colored_by_day, file=output_file, width=7.2, height=5, units="in")



##########################
# Make PCA PVE line plot
##########################
num_pcs <- 50
output_file <- paste0(visualize_processed_data_dir, "pca_pve_line_plot.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pca_pve[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")

num_pcs <- 50
output_file <- paste0(visualize_processed_data_dir, "pca_pve_top_500_vg_line_plot.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pca_pve_top_500_vg[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")

num_pcs <- 50
output_file <- paste0(visualize_processed_data_dir, "pca_pve_top_2000_vg_line_plot.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pca_pve_top_2000_vg[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")
}
