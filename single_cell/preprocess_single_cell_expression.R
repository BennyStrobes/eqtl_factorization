library(dplyr)
library(Seurat)
library(ggplot2)
library(reshape2)
library(feather)
library(sctransform)
library(Matrix)
library(cowplot)
library(future)
options(bitmapType='cairo')
options(warn=1)


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}


concatenate_umi_counts_across_all_cells <- function(raw_umi_count_dir, umi_count_file, cell_classes) {
	listy = c()
	for (class_iter in 1:length(cell_classes)) {
		class_name <- cell_classes[class_iter]
		input_umi_file <- paste0(raw_umi_count_dir, "umi_counts_", class_name, ".rds")
		temp <- readRDS(input_umi_file)
		listy <- c(listy, temp)
	}
	concatenated_umi_counts <- do.call(cbind, listy)
	saveRDS(concatenated_umi_counts, umi_count_file)
}

concatenate_meta_data_across_all_cells <- function(meta_data_dir, meta_data_file, cell_classes, umi_counts) {
	for (class_iter in 1:length(cell_classes)) {
		class_name <- cell_classes[class_iter]
		input_meta_data_file <- paste0(meta_data_dir, "metadata_", class_name, ".txt_per_cell.txt")
		meta_data <- read.delim(input_meta_data_file, header = T, row.names = 1)
		if (class_iter==1) {
			concatenated_mat = meta_data
		} else {
			concatenated_mat = rbind(concatenated_mat, meta_data)
		}
	}
	concatenated_meta_data <- concatenated_mat[match(colnames(umi_counts), rownames(concatenated_mat)), ]
	saveRDS(concatenated_meta_data, meta_data_file)
}


make_boxplot_of_cells_per_individual <- function(meta_data) {
	indi_arr <- c()
	cell_type_arr <- c()
	num_cell_arr <- c()

	cell_types <- as.character(unique(meta_data$ct_cov))
	individuals <- as.character(unique(meta_data$ind_cov))

	for (cell_type_iter in 1:length(cell_types)) {
		for (indi_iter in 1:length(individuals)) {
			cell_type <- cell_types[cell_type_iter]
			individual <- individuals[indi_iter]
			num_cells <- sum(as.character(meta_data$ct_cov) == cell_type & as.character(meta_data$ind_cov) == individual)

			indi_arr <- c(indi_arr, individual)
			cell_type_arr <- c(cell_type_arr, cell_type)
			num_cell_arr <- c(num_cell_arr, num_cells)
		}
	}
	df <- data.frame(individual=factor(indi_arr), cell_type=factor(cell_type_arr), number_of_cells=num_cell_arr)

	plotter <- ggplot(df,aes(x=cell_type, y=number_of_cells, fill=cell_type)) + 
    	geom_boxplot() + 
    	labs(x="Cell Type",y="Number of cells / individual") + 
    	gtex_v8_figure_theme() +
    	theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
    	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
    return(plotter)

}

make_histogram_of_log_mean_expression_per_gene <- function(umi_counts) {
	average_counts_per_gene <- rowMeans(umi_counts)
	df <- data.frame(average_counts_per_gene=log(average_counts_per_gene))
	p <- ggplot(df, aes(x=average_counts_per_gene)) + geom_histogram(color="darkblue", fill="lightblue") + gtex_v8_figure_theme() +
	labs(x="log(Average counts) / gene")
	return(p)
}

make_histogram_of_number_of_expressed_cells_per_gene <- function(umi_counts) {
	average_counts_per_gene <- rowSums(umi_counts>0.0)
	df <- data.frame(average_counts_per_gene=(average_counts_per_gene))
	p <- ggplot(df, aes(x=average_counts_per_gene)) + geom_histogram(color="darkgreen", fill="lightgreen") + gtex_v8_figure_theme() +
	labs(x="Number of expressed (count > 1) cells / gene)")
	return(p)
}

make_histogram_of_log_number_of_expressed_cells_per_gene <- function(umi_counts) {
	average_counts_per_gene <- rowSums(umi_counts>0.0)
	df <- data.frame(average_counts_per_gene=log(average_counts_per_gene))
	p <- ggplot(df, aes(x=average_counts_per_gene)) + geom_histogram(color="darkgreen", fill="lightgreen") + gtex_v8_figure_theme() +
	labs(x="log(Number of expressed (count > 1) cells) / gene)")
	return(p)
}

make_histogram_of_genes_expression <- function(expr_vec) {
	df <- data.frame(expr=expr_vec)
	p <- ggplot(df, aes(x=expr)) + geom_histogram(color="darkgreen", fill="lightgreen") + gtex_v8_figure_theme() +
		labs(x="Normalized expression")
	return(p)

}

make_histogram_of_genes_expression_across_genes <- function(expr) {
	set.seed(5)
	p1 <- make_histogram_of_genes_expression(expr[sample(1:10000, 1),])
	p2 <- make_histogram_of_genes_expression(expr[sample(1:10000, 1),])
	p3 <- make_histogram_of_genes_expression(expr[sample(1:10000, 1),])
	p4 <- make_histogram_of_genes_expression(expr[sample(1:10000, 1),])
	p5 <- make_histogram_of_genes_expression(expr[sample(1:10000, 1),])
	p6 <- make_histogram_of_genes_expression(expr[sample(1:10000, 1),])
	p7 <- make_histogram_of_genes_expression(expr[sample(1:10000, 1),])
	p8 <- make_histogram_of_genes_expression(expr[sample(1:10000, 1),])


	combined <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol=1)
	return(combined)
}

#########################
# Filter genes to:
# 1. protein_coding
# 2. gene_status=='KNOWN'
# 3. Autosomal
##########################
filter_genes <- function(counts, gene_annotation_file) {
	# Load in gene annotation file
	gene_anno <- read.table(gene_annotation_file,header=TRUE,sep='\t')
	# Get indices of gene annotation file that pass our filters
	filter =as.character(gene_anno$chr) != 'chrX' & as.character(gene_anno$chr) != 'chrY' & as.character(gene_anno$chr) != 'chrM' & as.character(gene_anno$gene_status) == "KNOWN" & as.character(gene_anno$gene_type) == "protein_coding"
	# Get gene ids that pass our filter
	genes_passed_filter <- gene_anno[filter,]$gene_name
	# Get row indices of counts data that are genes that pass our filters
	indices = rownames(counts) %in% as.character(genes_passed_filter)
	return(counts[indices,])
}

#########################
# Filter Individuals to those with disease_cov=="sle"
##########################
filter_to_individuals_with_sle <- function(umi_counts, meta_data) {
	# Get indices of cells
	indices <- as.character(meta_data$disease_cov) == "sle"
}

#########################
# Make heatmap showing correlation of loadings withknown covariates
##########################
make_pca_covariate_correlation_heatmap <- function(pcs, covariates) {
	################
	# Fair amount of manual processing
	################
	# Filter to covariates that could possible have impact on expression
	covs <- covariates[,c(2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18)]

    # Initialize PVE heatmap
    pve_map <- matrix(0, dim(covs)[2], dim(pcs)[2])
    colnames(pve_map) <- colnames(pcs)
    rownames(pve_map) <- colnames(covs)
    
    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(pcs)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- pcs[,num_pc]
            cov_vec <- covs[,num_cov]
            lin_model <- lm(pc_vec ~ cov_vec)
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "PC","PVE")

	melted_mat$Covariate <- factor(chartr("_", " ",melted_mat$Covariate), levels = chartr("_", " ",rownames(pve_map)[ord]))
	melted_mat$PC <- factor(melted_mat$PC)
	print(head(melted_mat))

    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=PC)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + labs(y="PC number",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 

    return(heatmap)
}

# Make plot showing variance explained of first n pcs
plot_pca_variance_explained <- function(variance_explained, N) {

    # Merge global vectors into a data frame
    df <- data.frame(variance_explained = variance_explained[1:N], pc_num = 1:N)

    # PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=variance_explained)) +
                geom_line() +
                geom_point() +
                ylim(0,max(variance_explained) + .01) + 
                scale_x_continuous(breaks=1:N) +
                labs(x = "PC Number", y = "Variance Explained") + 
                gtex_v8_figure_theme()
                 #theme(plot.title=element_text(size=8,face="plain"),text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8), axis.text.x = element_text(vjust=.5)) 

    # SAVE PLOT
    return(line_plot)
}

# Helper method to save DGE data structure to tab-deliminated text file
save_python_style_matrix <- function(counts, output_file, row_label_name) {
    #  Convert DGE data structure to matrix
    temp_mat <- as.matrix(counts)

    #  Edit colnames to include a header over the row-labels.
    revised_column_names <- colnames(temp_mat)
    revised_column_names[1] <- paste0(row_label_name,"\t",revised_column_names[1])

    write.table(temp_mat, output_file, quote=FALSE,col.names = revised_column_names, sep="\t")

}


#########################
# Load in command line args
##########################
args <- commandArgs(TRUE)
raw_umi_count_dir <- args[1]
meta_data_dir <- args[2]
processed_expression_dir <- args[3]
visualize_processed_expression_dir <- args[4]
normalization_method <- args[5]
gene_annotation_file <- args[6]


#########################
# Concatenate data across cell types
##########################
# Ordered list of cell classes
cell_classes <- c("B_cells", "CD14+_Monocytes", "CD4_T_cells", "CD8_T_cells", "Dendritic_cells", "FCGR3A+_Monocytes", "Megakaryocytes", "NK_cells")
# Generate sparse file containing umi counts from each cell types
umi_count_file <- paste0(processed_expression_dir, "raw_umi_all_cells.rds")
#concatenate_umi_counts_across_all_cells(raw_umi_count_dir, umi_count_file, cell_classes)
umi_counts <- readRDS(umi_count_file)
# Generate metadata file across all cells
meta_data_file <- paste0(processed_expression_dir, "meta_data_all_cells.rds")
#concatenate_meta_data_across_all_cells(meta_data_dir, meta_data_file, cell_classes, umi_counts)
meta_data <- readRDS(meta_data_file)

#########################
# Parameters
##########################
# Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
min.cells = 300
# Include cells where at least this many features are detected.
min.genes = 400
# Filter out cells with this amount of mitochondrial levels
mito_thresh = 5

#########################
# Filter Individuals to those with disease_cov=="sle"
##########################
indices <- as.character(meta_data$disease_cov) == "sle"
umi_counts <- umi_counts[,indices]
meta_data <- meta_data[indices,]


#########################
# Make boxplot of number of cells per individual, stratified by cell type
##########################
output_file <- paste0(visualize_processed_expression_dir, "boxplot_of_number_of_cells_per_individual_stratified_by_cell_type.pdf")
boxplot <- make_boxplot_of_cells_per_individual(meta_data)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")



#########################
# Generate Seurat Object
##########################
expr_seurat <- CreateSeuratObject(counts = umi_counts, meta.data = meta_data, min.cells = min.cells, min.features = min.genes)
print(dim(expr_seurat))
expr_seurat[["percent.mt"]] <- PercentageFeatureSet(expr_seurat, pattern = "^MT-")


#########################
# Visualize QC metrics as a violin plot
##########################
p <- VlnPlot(expr_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0.0)
ggsave(p, file=paste0(visualize_processed_expression_dir, "nfeature_violin_plot.pdf"), width=7.2, height=5.5, units="in")

#########################
# Further filter cells
##########################
expr_seurat <- subset(expr_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < mito_thresh)

#########################
# Make histogram showing number of expressed cells per gene
##########################
output_file <- paste0(visualize_processed_expression_dir, "histogram_showing_number_of_expressed_cells_per_gene.pdf")
filtered_umi_counts <- expr_seurat[["RNA"]]@data
hist <- make_histogram_of_number_of_expressed_cells_per_gene(filtered_umi_counts)
ggsave(hist, file=output_file, width=7.2, height=2.5, units="in")

output_file <- paste0(visualize_processed_expression_dir, "histogram_showing_log_number_of_expressed_cells_per_gene.pdf")
filtered_umi_counts <- expr_seurat[["RNA"]]@data
hist <- make_histogram_of_log_number_of_expressed_cells_per_gene(filtered_umi_counts)
ggsave(hist, file=output_file, width=7.2, height=2.5, units="in")

#########################
# Make histogram showing log mean expression across cells for a given
##########################
output_file <- paste0(visualize_processed_expression_dir, "histogram_showing_mean_expression_per_gene.pdf")
filtered_umi_counts <- expr_seurat[["RNA"]]@data
hist <- make_histogram_of_log_mean_expression_per_gene(filtered_umi_counts)
ggsave(hist, file=output_file, width=7.2, height=2.5, units="in")



#########################
# log normalize the data
##########################
if (normalization_method == "log") {
	expr_seurat <- NormalizeData(expr_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
	expr_seurat <- FindVariableFeatures(expr_seurat, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(expr_seurat)
	expr_seurat <- ScaleData(expr_seurat, features = all.genes)
	scaled_data <- expr_seurat[["RNA"]]@scale.data
}
if (normalization_method == "sctransform") {
	options(future.globals.maxSize = 1000 * 1024^2)
	expr_seurat <- SCTransform(expr_seurat, return.only.var.genes = FALSE, verbose = FALSE)
	scaled_data <- expr_seurat[["SCT"]]@scale.data
}
if (normalization_method == "log_with_covariates") {
	#options(future.globals.maxSize = 100000 * 1024^2)
	expr_seurat <- NormalizeData(expr_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
	expr_seurat <- FindVariableFeatures(expr_seurat, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(expr_seurat)
	#plan(strategy = "multicore", workers = 6)
	expr_seurat <- ScaleData(expr_seurat, features = all.genes, vars.to.regress = c("batch", "percent.mt", "disease_cov", "pop_cov"))
	scaled_data <- expr_seurat[["RNA"]]@scale.data
}
if (normalization_method == "sctransform_with_covariates") {
	options(future.globals.maxSize = 1000 * 1024^2)
	expr_seurat <- SCTransform(expr_seurat, return.only.var.genes = FALSE, verbose = FALSE, vars.to.regress = c("batch", "batch_cov", "pop_cov", "well", "percent.mt"))
	scaled_data <- expr_seurat[["SCT"]]@scale.data
}
if (normalization_method == "sctransform_with_covariates_and_individual") {
	options(future.globals.maxSize = 1000 * 1024^2)
	expr_seurat <- SCTransform(expr_seurat, return.only.var.genes = FALSE, verbose = FALSE, vars.to.regress = c("batch", "batch_cov", "pop_cov", "well", "ind_cov", "percent.mt"))
	scaled_data <- expr_seurat[["SCT"]]@scale.data
}
#############################
# Save scaled data to output
#############################
saveRDS(expr_seurat, paste0(processed_expression_dir, "seurat_", normalization_method, "_normalized_object.rds"))


#########################
# Run pca
##########################
expr_seurat <- RunPCA(expr_seurat, features = VariableFeatures(object = expr_seurat))

pca_cell_type <- DimPlot(expr_seurat, group.by="ct_cov", reduction = "pca", pt.size = .01)
pca_batch <- DimPlot(expr_seurat, group.by="batch", reduction = "pca", pt.size = .01)
pca_population <- DimPlot(expr_seurat, group.by="pop_cov", reduction = "pca", pt.size = .01)
pca_batch_cov <- DimPlot(expr_seurat, group.by="batch_cov", reduction = "pca", pt.size = .01)
pca_batch_well <- DimPlot(expr_seurat, group.by="well", reduction = "pca", pt.size = .01)
#pca_lib_size <- FeaturePlot(expr_seurat_no_regression, features="lib_size", reduction.use="pca")
#pca_sparsity <- FeaturePlot(expr_seurat_no_regression, features="sparsity", reduction.use="pca")

p <- CombinePlots(plots = list(pca_cell_type, pca_batch, pca_population, pca_batch_cov, pca_batch_well), ncol=2)
ggsave(p, file=paste0(visualize_processed_expression_dir, normalization_method, "_pca_plot.pdf"), width=13.2, height=7.2, units="in")

#########################
# Make heatmap showing correlation of loadings withknown covariates
##########################
pca_embeddings = Embeddings(expr_seurat, reduction = "pca")[, 1:10]
covariates = expr_seurat@meta.data
pca_covariate_correlation_heatmap <- make_pca_covariate_correlation_heatmap(pca_embeddings, covariates)
ggsave(pca_covariate_correlation_heatmap, file=paste0(visualize_processed_expression_dir, normalization_method, "_pca_covariate_correlation.pdf"), width=7.2, height=7.2, units="in")

#########################
# Run UMAP
##########################
expr_seurat <- RunUMAP(expr_seurat, dims = 1:10)

umap_cell_type <- DimPlot(expr_seurat, group.by="ct_cov", reduction = "umap", pt.size = .01)
umap_batch <- DimPlot(expr_seurat, group.by="batch", reduction = "umap", pt.size = .01)
umap_population <- DimPlot(expr_seurat, group.by="pop_cov", reduction = "umap", pt.size = .01)
umap_batch_cov <- DimPlot(expr_seurat, group.by="batch_cov", reduction = "umap", pt.size = .01)
umap_well <- DimPlot(expr_seurat, group.by="well", reduction = "umap", pt.size = .01)


umap_lib_size <- FeaturePlot(expr_seurat, features="lib_size")


p <- CombinePlots(plots = list(umap_cell_type, umap_batch, umap_population, umap_batch_cov, umap_well, umap_lib_size), ncol=2)
ggsave(p, file=paste0(visualize_processed_expression_dir, normalization_method, "_umap_plot.pdf"), width=13.2, height=7.2, units="in")


#############
# TEMPER
#############
#print("STARTING")
#expr_seurat <- readRDS(paste0(processed_expression_dir, "seurat_", normalization_method, "_normalized_object.rds"))
#scaled_data <- expr_seurat[["SCT"]]@scale.data
#print(dim(scaled_data))
#############
# END TEMPER
#############



################
# Remove non-protein coding genes
#################
scaled_data_gene_filtered = filter_genes(scaled_data, gene_annotation_file)


################
# Get expression PCs
#################
svd1 <- svd(scaled_data_gene_filtered)
variance_explained <- (svd1$d^2/sum(svd1$d^2))
loadings <- svd1$v[,1:10]

################
# Make VE plot
#################
ve_plot <- plot_pca_variance_explained(variance_explained,20)
ggsave(ve_plot,file=paste0(visualize_processed_expression_dir, normalization_method, "pc_pve_protein_coding_genes.pdf"),width=7.2, height=5, units="in")



################
# Regress Out PCs to get residual expression
#################
lm.fit <- lm(t(scaled_data_gene_filtered) ~ loadings)
residual_expression <- t(resid(lm.fit))


#################
# Save to output files
##################
# Save gene names
gene_names = rownames(residual_expression)
write.table(gene_names, paste0(processed_expression_dir, normalization_method, "_gene_names.txt"), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\n")
# Save expr
write.table(residual_expression, paste0(processed_expression_dir, normalization_method, "_residual_expression.txt"), row.names=FALSE,col.names=FALSE,sep="\t")
save_python_style_matrix(expr_seurat@meta.data, paste0(processed_expression_dir, normalization_method, "_cell_info.txt"), "cell_bar_code")

#####################################################
# OLD
#####################################################



#########################
# Filter genes to:
# 1. protein_coding
# 2. gene_status=='KNOWN'
# 3. Autosomal
##########################
#umi_counts = filter_genes(raw_umi_counts, gene_annotation_file)





#############################
# Save scaled data to output
#############################
#expr_seurat <- readRDS(paste0(processed_expression_dir, "seurat_", normalization_method, "_normalized_object.rds"))
#output_file <-paste0(visualize_processed_expression_dir, normalization_method, "_expression_across_cells.pdf")
#hist_plot <- make_histogram_of_genes_expression_across_genes(scaled_data)
#ggsave(hist_plot, file=output_file, width=7.2, height=12, units="in")

