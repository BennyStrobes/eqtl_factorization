args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(PRROC)
library(cowplot)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')




make_loading_scatter_plot <- function(tissue_file, loading_file) {
	tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	plotter <- ggplot(df) + 
	           geom_point( aes(x=loading_1, y=loading_2, color=tissue),size=.01) +
	           gtex_v8_figure_theme() + 
	           labs(x="Loading 1", y = "Loading 2", color="")
	return(plotter)
}

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=14), text = element_text(size=14),axis.text=element_text(size=13), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=14), legend.title = element_text(size=13)))
}

compare_loading_scatter_plot_to_true <- function(loading_true, loading) {
	num_samp <- dim(loading)[1]
	true <- c(loading_true$V1, loading_true$V2, loading_true$V3)
	learned <- c(loading$V3, loading$V1, loading$V2)
	names <- c(rep("loading1", num_samp), rep("loading2", num_samp), rep("loading3", num_samp))
	df <- data.frame(true=true, learned=learned, name=names)
	plotter <- ggplot(df) + 
	           geom_point( aes(x=true, y=learned, color=name),size=.01) +
	           gtex_v8_figure_theme() + 
	           labs(x="True loading", y = "Learned loading", color="")
	return(plotter)
}

compare_random_effects_scatter_plot_to_true <- function(loading_true, loading) {
	df <- data.frame(true=loading_true, learned=loading)
	plotter <- ggplot(df) + 
	           geom_point( aes(x=true, y=learned),size=.01) +
	           gtex_v8_figure_theme() + 
	           labs(x="True random effect", y = "Learned random effect", color="")
	return(plotter)
}

simulated_data_dir <- args[1]
eqtl_results_dir <- args[2]
visualization_dir <- args[3]
file_stem <- args[4]

# Compare random effects
true_random_effect_file <- paste0(simulated_data_dir, file_stem, "simulated_random_effects.txt")
learned_random_effect_file <- paste0(eqtl_results_dir, file_stem, "1.0_sparsity_parameter_re_mean.txt")
true_random_effects <- unlist(read.table(true_random_effect_file, header=FALSE), use.names=FALSE)
learned_random_effects <- unlist(read.table(learned_random_effect_file, header=FALSE), use.names=FALSE)
output_file <- paste0(visualization_dir, file_stem, "random_effect_comparison_scatter.pdf")
mysample <- sample(length(true_random_effects), 10000, replace=FALSE)
random_effects_scatter <- compare_random_effects_scatter_plot_to_true(true_random_effects[mysample], learned_random_effects[mysample])
ggsave(random_effects_scatter, file=output_file, width=7.2, height=5.5, units="in")

# Compare loadings matrices
loading_true_file <- paste0(simulated_data_dir, file_stem, "simulated_loading_matrix.txt")
loading_learned_file <- paste0(eqtl_results_dir, file_stem, "1.0_sparsity_parameter_qU_mean.txt")
loading_true <- read.table(loading_true_file, header=FALSE)
loading <- read.table(loading_learned_file, header=FALSE)
output_file <- paste0(visualization_dir, file_stem, "loading_comparison_scatter.pdf")
loading_scatter <- compare_loading_scatter_plot_to_true(loading_true, loading)
ggsave(loading_scatter, file=output_file, width=7.2, height=5.5, units="in")