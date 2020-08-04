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


make_number_of_egenes_per_cell_type_bar_plot <- function(cell_types, pseudobulk_eqtl_dir, processed_expression_dir) {
	cell_type_arr <- c()
	egene_arr <- c()
	num_cell_arr <- c()
	covariates <- read.table(paste0(processed_expression_dir, "pseudobulk_covariates_sle_individuals.txt"), header=TRUE,sep="\t")

	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		egene_file <- paste0(pseudobulk_eqtl_dir, cell_type, "_pseudobulk_eqtl_analysis_multiple_testing_bf_bh_0.05_fdr_.txt")
		egene_data <- read.table(egene_file, header=TRUE)
		num_egenes <- dim(egene_data)[1]

		num_cells = sum(covariates$num_cells[covariates$ct_cov_readable==cell_type])

		cell_type_arr <- c(cell_type_arr, cell_type)
		egene_arr <- c(egene_arr, num_egenes)
		num_cell_arr <- c(num_cell_arr, num_cells)
	}


	df <- data.frame(cell_type=factor(cell_type_arr), num_egene=egene_arr, num_cells=num_cell_arr)
	print(df[order(df$num_cells),])
	df2 <- df[order(df$num_cells),]
	print(df2)
	#ordered_cell_types = as.character(df[order(df$num_cells),]$cell_type)
	#rint(ordered_cell_types)

	#df2 = data.frame(cell_type=factor(cell_type_arr, levels=ordered_cell_types), num_egene=egene_arr, num_cells=num_cell_arr)
	#rint(df2)

	p<-ggplot(data=df2, aes(x=reorder(cell_type,num_cells), y=num_egene)) +
  		geom_bar(stat="identity") +
  		coord_flip() + 
  		gtex_v8_figure_theme() + 
  		labs(y="Number of pseudo-bulk eGenes (FDR < .05)", x="")
  	return(p)
}

summary_stat_scatter_plot <- function(title, pseudobulk_stat, sc_stat, cov, x_axis_label, y_axis_label, color_label) {
	df <- data.frame(pseudo=pseudobulk_stat, sc=sc_stat, cov=cov)

	missing <- !is.na(df$sc) 

	corry = cor(df$pseudo[missing], df$sc[missing])
	# Basic scatter plot
	scatter <- ggplot(df, aes(x=pseudo, y=sc, color=cov)) + geom_point(size=.1) +
		labs(x=x_axis_label, y=y_axis_label, color=color_label) + 
		labs(title=paste0(title, " / Pearson rho: ", corry)) +
		figure_theme() + 
		geom_hline(yintercept=0, color = "grey", size=.1) + 
		geom_vline(xintercept=0, color = "grey", size=.1) + 
		geom_abline(color="grey", size=.1)

	return(scatter)
}

pvalue_scatter_plot <- function(title, pseudobulk_stat, sc_stat, cov, x_axis_label, y_axis_label, color_label) {
	pseudobulk_stat[pseudobulk_stat < 0.0] = 0.0
	sc_stat[sc_stat < 0.0] = 0.0
	df <- data.frame(pseudo=pseudobulk_stat, sc=sc_stat, cov=cov)
	missing <- !is.na(df$sc) 
	corry = cor(df$pseudo[missing], df$sc[missing])
	# Basic scatter plot
	scatter <- ggplot(df, aes(x=pseudo, y=sc, color=cov)) + geom_point(size=.1) +
		labs(x=x_axis_label, y=y_axis_label, color=color_label) + 
		labs(title=paste0(title, " / Pearson rho: ", corry)) +
		figure_theme() + 
		geom_abline(color = "grey", size=.1)

	return(scatter)
}




processed_expression_dir <- args[1]
pseudobulk_eqtl_dir <- args[2]
single_cell_eqtl_dir <- args[3]
visualize_pseudobulk_eqtl_dir <- args[4]


cell_type_file <- paste0(pseudobulk_eqtl_dir, "cell_types.txt")
cell_types <- as.character(read.table(cell_type_file)$V1)


######################################################
# Make Bar plot showing number of eGenes per cell type
######################################################
if (FALSE) {
output_file <- paste0(visualize_pseudobulk_eqtl_dir, "number_egenes_per_cell_type_bar_plot.pdf")
bar_plot <- make_number_of_egenes_per_cell_type_bar_plot(cell_types, pseudobulk_eqtl_dir, processed_expression_dir)
ggsave(bar_plot, file=output_file, width=7.2, height=6.0, units="in")
}




cell_type <- "B_cells"
pseudobulk_eqtl_file <- paste0(pseudobulk_eqtl_dir, cell_type, "_pseudobulk_eqtl_analysis_all_variant_gene_pairs.txt")
sc_eqtl_0_pc_file <- paste0(single_cell_eqtl_dir, cell_type, "_sc_eqtl_analysis_0_pcs_all_variant_gene_pairs_merged.txt")
sc_eqtl_25_pc_file <- paste0(single_cell_eqtl_dir, cell_type, "_sc_eqtl_analysis_25_pcs_all_variant_gene_pairs_merged.txt")
sc_nb_eqtl_25_pc_file <- paste0(single_cell_eqtl_dir, cell_type, "_sc_nb_eqtl_analysis_25_pcs_all_variant_gene_pairs_merged.txt")

test_info_file <- paste0(single_cell_eqtl_dir, cell_type, "_eqtl_input_test_info.txt")

ct_pseudobulk_eqtls <- read.table(pseudobulk_eqtl_file, header=TRUE)
ct_sc_0_pc_eqtls <- read.table(sc_eqtl_0_pc_file, header=TRUE)
ct_sc_25_pc_eqtls <- read.table(sc_eqtl_25_pc_file, header=TRUE)
ct_nb_sc_25_pc_eqtls <- read.table(sc_nb_eqtl_25_pc_file, header=TRUE)

ct_nb_sc_25_pc_eqtls$beta[ct_nb_sc_25_pc_eqtls$beta > 1.5] = 1.5
ct_nb_sc_25_pc_eqtls$beta[ct_nb_sc_25_pc_eqtls$beta < -1.5] = -1.5
ct_nb_sc_25_pc_eqtls$pseudo_r_squared[ct_nb_sc_25_pc_eqtls$pseudo_r_squared < 0] = 0
ct_nb_sc_25_pc_eqtls$pseudo_r_squared[ct_nb_sc_25_pc_eqtls$pseudo_r_squared > 1] = 1

ct_test_info <- read.table(test_info_file, header=TRUE)


#########################################
# Pseudo-bulk LMM beta comparison scatter
########################################
output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_beta_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 0 PCs")
scatter_plot <- summary_stat_scatter_plot(title, ct_pseudobulk_eqtls$beta, ct_sc_0_pc_eqtls$beta, ct_test_info$percent_expressed_cells, "pseudobulk beta", "sc beta", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_beta_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 0 PCs")
scatter_plot <- summary_stat_scatter_plot(title, ct_pseudobulk_eqtls$beta, ct_sc_0_pc_eqtls$beta, log(ct_test_info$total_counts), "pseudobulk beta", "sc beta", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_25_pc_beta_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- summary_stat_scatter_plot(title, ct_pseudobulk_eqtls$beta, ct_sc_25_pc_eqtls$beta, ct_test_info$percent_expressed_cells, "pseudobulk beta", "sc beta", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_25_pc_beta_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- summary_stat_scatter_plot(title, ct_pseudobulk_eqtls$beta, ct_sc_25_pc_eqtls$beta, log(ct_test_info$total_counts), "pseudobulk beta", "sc beta", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

#########################################
# Pseudo-bulk LMM pvalue comparison scatter
########################################
output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_pvalue_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 0 PCs")
scatter_plot <- pvalue_scatter_plot(title, -log10(ct_pseudobulk_eqtls$pvalue + 1e-15), -log10(ct_sc_0_pc_eqtls$pvalue + 1e-15), ct_test_info$percent_expressed_cells, "pseudobulk -log10(pvalue)", "sc -log10(pvalue)", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_pvalue_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 0 PCs")
scatter_plot <- pvalue_scatter_plot(title, -log10(ct_pseudobulk_eqtls$pvalue + 1e-15), -log10(ct_sc_0_pc_eqtls$pvalue + 1e-15), log(ct_test_info$total_counts), "pseudobulk -log10(pvalue)", "sc -log10(pvalue)", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_25_pc_pvalue_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- pvalue_scatter_plot(title, -log10(ct_pseudobulk_eqtls$pvalue + 1e-15), -log10(ct_sc_25_pc_eqtls$pvalue + 1e-15), ct_test_info$percent_expressed_cells, "pseudobulk -log10(pvalue)", "sc -log10(pvalue)", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_25_pc_pvalue_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- pvalue_scatter_plot(title, -log10(ct_pseudobulk_eqtls$pvalue + 1e-15), -log10(ct_sc_25_pc_eqtls$pvalue + 1e-15), log(ct_test_info$total_counts), "pseudobulk -log10(pvalue)", "sc -log10(pvalue)", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")


#########################################
# Pseudo-bulk-NB beta comparison scatter
########################################
output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_nb_25_pc_beta_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- summary_stat_scatter_plot(title, ct_pseudobulk_eqtls$beta, ct_nb_sc_25_pc_eqtls$beta, ct_test_info$percent_expressed_cells, "pseudobulk beta", "sc nb beta", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_nb_25_pc_beta_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- summary_stat_scatter_plot(title, ct_pseudobulk_eqtls$beta, ct_nb_sc_25_pc_eqtls$beta, log(ct_test_info$total_counts), "pseudobulk beta", "sc nb beta", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")


#########################################
# sc_gaussian-sc_NB beta comparison scatter
########################################
output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_sc_gaussian_sc_nb_25_pc_beta_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- summary_stat_scatter_plot(title, ct_sc_25_pc_eqtls$beta, ct_nb_sc_25_pc_eqtls$beta, ct_test_info$percent_expressed_cells, "sc gaussian beta", "sc nb beta", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_sc_gaussian_sc_nb_25_pc_beta_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- summary_stat_scatter_plot(title, ct_sc_25_pc_eqtls$beta, ct_nb_sc_25_pc_eqtls$beta, log(ct_test_info$total_counts), "sc gaussian beta", "sc nb beta", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")



if (FALSE) {
#########################################
# pseudobulk-sc_gaussian pseudo-r^2 comparison scatter
########################################
output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_gaussian_25_pc_pseudo_rsquared_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- pvalue_scatter_plot(title, ct_pseudobulk_eqtls$pseudo_r_squared, ct_sc_25_pc_eqtls$pseudo_r_squared_lm, ct_test_info$percent_expressed_cells, "pseudobulk pseudo R^2", "sc gaussian pseudo R^2", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_gaussian_25_pc_pseudo_rsquared_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- pvalue_scatter_plot(title, ct_pseudobulk_eqtls$pseudo_r_squared, ct_sc_25_pc_eqtls$pseudo_r_squared_lm, log(ct_test_info$total_counts), "pseudobulk pseudo R^2", "sc gaussian pseudo R^2", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

#########################################
# pseudobulk-sc_nb pseudo-r^2 comparison scatter
########################################
output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_nb_25_pc_pseudo_rsquared_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- pvalue_scatter_plot(title, ct_pseudobulk_eqtls$pseudo_r_squared, ct_nb_sc_25_pc_eqtls$pseudo_r_squared, ct_test_info$percent_expressed_cells, "pseudobulk pseudo R^2", "sc nb pseudo R^2", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_pseudobulk_sc_nb_25_pc_pseudo_rsquared_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- pvalue_scatter_plot(title, ct_pseudobulk_eqtls$pseudo_r_squared, ct_nb_sc_25_pc_eqtls$pseudo_r_squared, log(ct_test_info$total_counts), "pseudobulk pseudo R^2", "sc nb pseudo R^2", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")


#########################################
# sc_gaussian-sc_NB pseudo-r^2 comparison scatter
########################################
output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_sc_gaussian_sc_nb_25_pc_pseudo_rsquared_comparison_scatter_colored_by_percent_expressed.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- pvalue_scatter_plot(title, ct_sc_25_pc_eqtls$pseudo_r_squared_lm, ct_nb_sc_25_pc_eqtls$pseudo_r_squared, ct_test_info$percent_expressed_cells, "sc gaussian pseudo R^2", "sc nb pseudo R^2", "Fraction\nexpressed cells")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualize_pseudobulk_eqtl_dir, cell_type, "_sc_gaussian_sc_nb_25_pc_pseudo_rsquared_comparison_scatter_colored_by_total_counts.pdf")
title <- paste0(cell_type, " 25 PCs")
scatter_plot <- pvalue_scatter_plot(title, ct_sc_25_pc_eqtls$pseudo_r_squared_lm, ct_nb_sc_25_pc_eqtls$pseudo_r_squared, log(ct_test_info$total_counts), "sc gaussian pseudo R^2", "sc nb pseudo R^2", "log(test depth)")
ggsave(scatter_plot, file=output_file, width=7.2, height=6.0, units="in")
}
