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



processed_expression_dir <- args[1]
pseudobulk_eqtl_dir <- args[2]
visualize_pseudobulk_eqtl_dir <- args[3]


cell_type_file <- paste0(pseudobulk_eqtl_dir, "cell_types.txt")
cell_types <- as.character(read.table(cell_type_file)$V1)




######################################################
# Make Bar plot showing number of eGenes per cell type
######################################################
output_file <- paste0(visualize_pseudobulk_eqtl_dir, "number_egenes_per_cell_type_bar_plot.pdf")
bar_plot <- make_number_of_egenes_per_cell_type_bar_plot(cell_types, pseudobulk_eqtl_dir, processed_expression_dir)
ggsave(bar_plot, file=output_file, width=7.2, height=6.0, units="in")


