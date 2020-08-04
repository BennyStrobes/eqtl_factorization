args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_scatter_plot_with_color <- function(x_data, y_data, color_data, x_axis_label, y_axis_label, color_axis_label) {
	df <- data.frame(x_data=x_data, y_data=y_data, color_data=color_data)

	plotter <- ggplot(df) + 
	           geom_point(aes(x=x_data, y=y_data,color=color_data)) +
	           gtex_v8_figure_theme() + 
	           labs(x=x_axis_label, y = y_axis_label, color=color_axis_label) +
	           geom_abline(size=.1, color="grey")
	return(plotter)
}


make_scatter_plot <- function(x_data, y_data, x_axis_label, y_axis_label) {
	df <- data.frame(x_data=x_data, y_data=y_data)

	plotter <- ggplot(df) + 
	           geom_point(aes(x=x_data, y=y_data)) +
	           gtex_v8_figure_theme() + 
	           labs(x=x_axis_label, y = y_axis_label) +
	           geom_abline(size=.1, color="grey")
	return(plotter)
}





dir_name <- args[1]

r_squared_rand_1_file <- paste0(dir_name, "r_squared_run_1.txt")
r_squared_rand_2_file <- paste0(dir_name, "r_squared_run_2.txt")
r_squared_cell_type_loading_file <- paste0(dir_name, "r_squared_run_cell_type_loading.txt")
percent_expressed_cells_file <- paste0(dir_name, "percent_cells_expresssed_per_test.txt")

r_squared_rand_1 <- read.table(r_squared_rand_1_file, header=FALSE)
r_squared_rand_2 <- read.table(r_squared_rand_2_file, header=FALSE)
r_squared_cell_type_loading <- read.table(r_squared_cell_type_loading_file, header=FALSE)
percent_expressed_cells <- read.table(percent_expressed_cells_file, header=FALSE)

#####################################
# Correlation scatter plot of R^2 of each gene for two independent runs
#####################################
scatter <- make_scatter_plot_with_color(r_squared_rand_1$V1, r_squared_rand_2$V1, percent_expressed_cells$V1, "R^2 of test in eqtl factorization (rand seed 1)", "R^2 of test in eqtl factorization (rand seed 2)", "fraction cells expressed")
output_file <- paste0(dir_name, "eqtl_factorization_rand_1_vs_eqtl_factorization_rand_2_r_squared_scatter.pdf")
ggsave(scatter, file=output_file, width=7.2, height=6, units="in")


#####################################
# Correlation scatter plot of R^2 of each gene for two independent runs
#####################################
scatter <- make_scatter_plot_with_color(r_squared_rand_1$V1, r_squared_cell_type_loading$V1, percent_expressed_cells$V1, "R^2 of test in eqtl factorization (rand seed 1)", "R^2 of test in eqtl factorization (loading fixed to cell type)", "fraction cells expressed")
output_file <- paste0(dir_name, "eqtl_factorization_rand_1_vs_eqtl_factorization_cell_type_loading_r_squared_scatter.pdf")
ggsave(scatter, file=output_file, width=7.2, height=6, units="in")


#####################################
# Correlation scatter plot of R^2 with fraction of expressed cells
#####################################
scatter <- make_scatter_plot(r_squared_cell_type_loading$V1, percent_expressed_cells$V1, "R^2 of test in eqtl factorization (loading fixed to cell type)", "fraction cells expressed")
output_file <- paste0(dir_name, "eqtl_factorization_cell_type_loading_r_squared_vs_fraction_expressed_cells_scatter.pdf")
ggsave(scatter, file=output_file, width=7.2, height=6, units="in")
