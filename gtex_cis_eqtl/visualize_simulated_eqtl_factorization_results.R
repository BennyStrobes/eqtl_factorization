args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')



figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



#Make heatmap showing correlation of loadings between models
make_loading_correlation_heatmap <- function(model1_loadings, model2_loadings, x_axis_label, y_axis_label) {
	num_dim_x = dim(model1_loadings)[2]
	num_dim_y = dim(model2_loadings)[2]
	
	corr_matrix = matrix(0, num_dim_x, num_dim_y)

	for (x_index in 1:num_dim_x) {
		for (y_index in 1:num_dim_y) {
			corr_matrix[x_index, y_index] <- abs(cor(model1_loadings[,x_index], model2_loadings[, y_index]))
		}
	}
	
	melted_mat <- melt(corr_matrix)
    colnames(melted_mat) <- c("model1", "model2", "correlation")


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=model1, y=model2)) + geom_tile(aes(fill=correlation)) + 
		figure_theme() +
   		labs(y=y_axis_label, x=x_axis_label, fill="Absolute\nPearson correlation") +
   		#theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		#theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="white",high="blue") +
   		scale_x_continuous(breaks=1:num_dim_x, labels=1:num_dim_x) +
   		scale_y_continuous(breaks=1:num_dim_y, labels=1:num_dim_y)
	return(heatmap)

}








######################
# Command line args
######################
eqtl_results_dir <- args[1]
visualization_dir <- args[2]





#########################
# Load in specific model
#########################
num_individuals="1000"
num_samples_per_individual="1"
num_tests="1000"
num_latent_factors="7"
seed="0"

# File stems
svi="False"
model_stem_no_svi=paste0(eqtl_results_dir, num_individuals, "_individuals_", num_samples_per_individual, "_samples_per_individual_", num_tests, "_tests_", num_latent_factors, "_latent_factors_", svi, "_svi_boolean_", seed, "_seed_")
svi="True"
model_stem_svi=paste0(eqtl_results_dir, num_individuals, "_individuals_", num_samples_per_individual, "_samples_per_individual_", num_tests, "_tests_", num_latent_factors, "_latent_factors_", svi, "_svi_boolean_", seed, "_seed_")

# Load in loading matrices
no_svi_loadings <- read.table(paste0(model_stem_no_svi, "U_S.txt"), header=FALSE)
svi_loadings <- read.table(paste0(model_stem_svi, "U_S.txt"), header=FALSE)
simulated_loadings <- read.table(paste0(model_stem_no_svi, "actual_U_S.txt"), header=FALSE)


######################################
# Make heatmap showing correlation of loadings between models
######################################
# SVI vs no SVI
output_file <- paste0(visualization_dir, "svi_vs_no_svi_loading_correlation_heatmap.pdf")
svi_vs_no_svi_heatmap <- make_loading_correlation_heatmap(no_svi_loadings, svi_loadings, "Loadings (standard VI)", "Loadings (Stochastic VI)")
ggsave(svi_vs_no_svi_heatmap, file=output_file, width=7.2, height=5.5, units="in")


# Simulated (actual) vs SVI
output_file <- paste0(visualization_dir, "svi_vs_actual_loading_correlation_heatmap.pdf")
svi_vs_simulated_heatmap <- make_loading_correlation_heatmap(simulated_loadings, svi_loadings, "Loadings (simulated)", "Loadings (Stochastic VI)")
ggsave(svi_vs_simulated_heatmap, file=output_file, width=7.2, height=5.5, units="in")



