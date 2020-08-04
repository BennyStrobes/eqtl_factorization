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




make_umap_loading_scatter_plot_colored_by_categorical_variable <- function(covariates, umap_loadings, covariate_name) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=factor(covariates))
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
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) 

}


make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate <- function(covariate, loadings, covariate_name) {
	loading_vec <- c()
	covariate_vec <- c()
	num_factors <- dim(loadings)[2]

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

	factor_number <- 6
	factor_6_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 7
	factor_7_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 8
	factor_8_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 9
	factor_9_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 10
	factor_10_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 11
	factor_11_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 12
	factor_12_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 13
	factor_13_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)


	combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, factor_4_boxplot, factor_5_boxplot, factor_6_boxplot, factor_7_boxplot, factor_8_boxplot, factor_9_boxplot, factor_10_boxplot, factor_11_boxplot, factor_12_boxplot, factor_13_boxplot, ncol=1)

	return(combined)
}

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
	# Covariates columns to consider
	valid_covariates <- c(2, 4, 5, 6, 9, 10, 11)
	#print(summary(covariates[, valid_covariates]))
	# Remove unimportant columns
    loadings <- as.matrix(loadings)
    covs <- covariates[,valid_covariates]



    # Initialize PVE heatmap
    factor_colnames <- paste0("Factor", 1:(dim(loadings)[2]))
    factor_rownames <- colnames(covs)
    pve_map <- matrix(0, dim(covs)[2], dim(loadings)[2])
    colnames(pve_map) <- factor_colnames
    rownames(pve_map) <- colnames(covs)


    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(loadings)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- loadings[,num_pc]
            cov_vec <- covs[,num_cov]
            lin_model <- lm(pc_vec ~ cov_vec)
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    print(pve_map)
    
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

#############################
# Command line args
#############################
processed_expression_dir <- args[1]
eqtl_input_dir <- args[2]
eqtl_results_dir <- args[3]
eqtl_mixture_results_dir <- args[4]
eqtl_visualization_dir <- args[5]


if (FALSE) {
# Input files
covariate_file <- paste0(processed_expression_dir, "pseudobulk_covariates_sle_individuals.txt")
covariates <- read.table(covariate_file, header=TRUE, sep="\t")

print(covariate_file)
type = "pseudobulk"
num_components <- "3"
loading_file <- paste0(eqtl_mixture_results_dir, type, "_mixture_", num_components, "_components_flexmix_model_posterior_prob.txt")
print(loading_file)
loadings <- read.table(loading_file, header=FALSE)

output_file <- paste0(eqtl_visualization_dir, "eqtl_mixture_", type, "_", num_components, "_loading_cell_type_proportions.pdf")
plotter <- make_mixture_model_loading_proportion_plot(loadings, covariates, paste0("Pseudo-bulk /", num_components, " components"))
ggsave(plotter, file=output_file, width=7.2, height=6, units="in")


type = "pseudobulk"
num_components <- "5"
loading_file <- paste0(eqtl_mixture_results_dir, type, "_mixture_", num_components, "_components_flexmix_model_posterior_prob.txt")
print(loading_file)
loadings <- read.table(loading_file, header=FALSE)

output_file <- paste0(eqtl_visualization_dir, "eqtl_mixture_", type, "_", num_components, "_loading_cell_type_proportions.pdf")
plotter <- make_mixture_model_loading_proportion_plot(loadings, covariates, paste0("Pseudo-bulk /", num_components, " components"))
ggsave(plotter, file=output_file, width=7.2, height=6, units="in")



type = "pseudobulk"
num_components <- "7"
loading_file <- paste0(eqtl_mixture_results_dir, type, "_mixture_", num_components, "_components_flexmix_model_posterior_prob.txt")
print(loading_file)
loadings <- read.table(loading_file, header=FALSE)

output_file <- paste0(eqtl_visualization_dir, "eqtl_mixture_", type, "_", num_components, "_loading_cell_type_proportions.pdf")
plotter <- make_mixture_model_loading_proportion_plot(loadings, covariates, paste0("Pseudo-bulk /", num_components, " components"))
ggsave(plotter, file=output_file, width=7.2, height=6, units="in")

# Input files
covariate_file <- paste0(processed_expression_dir, "cell_covariates_sle_individuals_random_subset.txt")
covariates <- read.table(covariate_file, header=TRUE, sep="\t")


type = "sc_random_subset"
num_components <- "3"
loading_file <- paste0(eqtl_mixture_results_dir, type, "_mixture_", num_components, "_components_flexmix_model_posterior_prob.txt")
print(loading_file)
loadings <- read.table(loading_file, header=FALSE)

output_file <- paste0(eqtl_visualization_dir, "eqtl_mixture_", type, "_", num_components, "_loading_cell_type_proportions.pdf")
plotter <- make_mixture_model_loading_proportion_plot(loadings, covariates, paste0("single-cell /", num_components, " components"))
ggsave(plotter, file=output_file, width=7.2, height=6, units="in")



type = "sc_random_subset"
num_components <- "5"
loading_file <- paste0(eqtl_mixture_results_dir, type, "_mixture_", num_components, "_components_flexmix_model_posterior_prob.txt")
#print(loading_file)
loadings <- read.table(loading_file, header=FALSE)

output_file <- paste0(eqtl_visualization_dir, "eqtl_mixture_", type, "_", num_components, "_loading_cell_type_proportions.pdf")
plotter <- make_mixture_model_loading_proportion_plot(loadings, covariates, paste0("single-cell /", num_components, " components"))
ggsave(plotter, file=output_file, width=7.2, height=6, units="in")



}


print("Hello")

frac_expressed = "0.1"
# Input files
covariate_file <- paste0(processed_expression_dir, "cell_covariates_sle_individuals_random_subset_min_expressed_cells_", frac_expressed, ".txt")
#covariate_file <- paste0(processed_expression_dir, "cell_covariates_sle_individuals.txt")
#covariate_file <- paste0(processed_expression_dir, "pseudobulk_covariates_sle_individuals.txt")

#eqtl_factorization_loading_file <- paste0(eqtl_results_dir, "eqtl_factorization_pseudobulk_data_20_factors_eqtl_factorization_vi_spike_and_slab_model_True_re_False_svi_0_seed_U_S.txt")
#eqtl_factorization_loading_file <- "/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_U_als.txt"
#eqtl_factorization_loading_file <- paste0(eqtl_results_dir, "eqtl_factorization_single_cell_sig_tests_50_pc_min_expressed_cells_", frac_expressed, "_data_uncorrected_genotype_5_factors_eqtl_factorization_als_model_False_re_False_svi_0_seed_U.txt")
eqtl_factorization_loading_file <- "/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/temp_output/temper_U_als.txt"

# Load in data
covariates <- read.table(covariate_file, header=TRUE, sep="\t")
loadings <- read.table(eqtl_factorization_loading_file, header=FALSE)
# Filter loadings file

#loadings <- loadings[,good_loadings]

# Create UMAP factors
umap_loadings = umap(loadings)$layout
saveRDS( umap_loadings, "umap_loadings.rds")
print("UMAP DONE")
#umap_loadings <- readRDS("umap_loadings.rds")


eqtl_visualization_dir <- paste0(eqtl_visualization_dir, frac_expressed, "_")
######################################
# Make scatter-plot of loading1 vs loading2 colored by known cell type
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading1_vs_loading2_colored_by_known_cell_type.pdf")
#scatter <- loading_scatter_plot_colored_by_categorical_covariate(covariates$ct_cov, loadings, 1, 2, "Known cell type")
#ggsave(scatter, file=output_file, width=7.2, height=6, units="in")

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
output_file <- paste0(eqtl_visualization_dir, "covariate_loading_correlation_heatmap.pdf")
#heatmap <- make_covariate_loading_correlation_heatmap(covariates, loadings)
#ggsave(heatmap, file=output_file, width=7.2, height=6, units="in")

######################################
# Make loading boxplot with row for every factor colored by cell type
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_cell_type.pdf")
#boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$ct_cov, loadings, "Known cell type")
#ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")

######################################
# Make loading boxplot with row for every factor colored by individual
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_individual.pdf")
#boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$ind_cov, loadings, "Individual")
#ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")

######################################
# Make loading boxplot with row for every factor colored by batch
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_batch.pdf")
#boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$batch_cov, loadings, "Batch")
#ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")

######################################
# Make loading boxplot colored by cell type
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_boxplot_colored_by_cell_type.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$ct_cov, loadings, "Known cell type")
ggsave(boxplot, file=output_file, width=10.2, height=5.5, units="in")
print("done")
######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(eqtl_visualization_dir, "loading_boxplot_colored_by_ancestry.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$pop_cov, loadings, "Known ancestry")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell type
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_cell_type.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$ct_cov, umap_loadings, "Known cell type")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known ancestry
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_ancestry.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$pop_cov, umap_loadings, "Known ancestry")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by male/femaile
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_gender.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$Female, umap_loadings, "Known female/male")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by batch
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_batch.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$batch_cov, umap_loadings, "Known batch")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by individual
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_individual.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$ind_cov, umap_loadings, "Known individual")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by known percent mito
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_percent_mito.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$percent_mito, umap_loadings, "Percent mito")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_cell_counts.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$n_counts, umap_loadings, "Cell counts")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by known number of genes
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_number_of_genes.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$n_genes, umap_loadings, "Number of genes")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

