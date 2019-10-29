args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(PRROC)
library(cowplot)
library(umap)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

 make_expr_clustering_heatmap <- function(expr_file) {
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
   		labs(y="RNA-seq Sample", x="Variant-gene", fill="expr") +
   		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="pink",high="blue") +
   		scale_x_discrete(expand=c(0,0)) +
 		scale_y_discrete(expand=c(0,0))
	return(heatmap)
 }


 make_abs_expr_clustering_heatmap <- function(expr_file) {
 	# Get matrix of dimension num_genesXnum_indi
 	residual_mat <- abs(t(t(read.table(expr_file))))
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
   		labs(y="RNA-seq Sample", x="Variant-gene", fill="Abs expr") +
   		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="pink",high="blue") +
   		scale_x_discrete(expand=c(0,0)) +
 		scale_y_discrete(expand=c(0,0))
	return(heatmap)
 }



output_dir <- args[1]



raw_expression_file <- paste0(output_dir, "raw_gene_filtered_reordered.txt")


residual_expression_file <- paste0(output_dir, "resid_gene_filtered_reordered.txt")


personal_expression_file <- paste0(output_dir, "personal_resid_gene_filtered_reordered.txt")

output_file <- paste0(output_dir, "raw_expr_real_valued_heatmap.pdf")
heatmap <- make_expr_clustering_heatmap(raw_expression_file)
ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")


output_file <- paste0(output_dir, "resid_expr_real_valued_heatmap.pdf")
heatmap <- make_expr_clustering_heatmap(residual_expression_file)
ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")

output_file <- paste0(output_dir, "personal_resid_expr_real_valued_heatmap.pdf")
heatmap <- make_expr_clustering_heatmap(personal_expression_file)
ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")

print("DONE")


output_file <- paste0(output_dir, "raw_expr_abs_valued_heatmap.pdf")
heatmap <- make_abs_expr_clustering_heatmap(raw_expression_file)
ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")


output_file <- paste0(output_dir, "resid_expr_abs_valued_heatmap.pdf")
heatmap <- make_abs_expr_clustering_heatmap(residual_expression_file)
ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")


output_file <- paste0(output_dir, "personal_resid_expr_abs_valued_heatmap.pdf")
heatmap <- make_abs_expr_clustering_heatmap(personal_expression_file)
ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")