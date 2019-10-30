library(dplyr)
library(Seurat)
library(ggplot2)
library(feather)
library(sctransform)
library(Matrix)
options(bitmapType='cairo')



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


#########################
# Load in command line args
##########################
args <- commandArgs(TRUE)
raw_umi_count_dir <- args[1]
meta_data_dir <- args[2]
processed_expression_dir <- args[3]
visualize_processed_expression_dir <- args[4]



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
min.cells = 100
# Include cells where at least this many features are detected.
min.genes = 600


#########################
# Generate Seurat Object
##########################
expr_seurat <- CreateSeuratObject(counts = umi_counts, meta.data = meta_data, min.cells = min.cells, min.features = min.genes)

expr_seurat[["percent.mt"]] <- PercentageFeatureSet(expr_seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
p <- VlnPlot(expr_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(p, file=paste0(visualize_processed_expression_dir, "nfeature_violin_plot.png"), width=7.2, height=5.5, units="in")


