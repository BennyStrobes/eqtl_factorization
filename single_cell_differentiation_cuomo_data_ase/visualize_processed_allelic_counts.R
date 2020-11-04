args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')




generate_pseudotime_af_heatmap <- function(allelic_fraction_pseutome_binned_file) {
	af <- read.table(allelic_fraction_pseutome_binned_file, header=TRUE, sep="\t")
	rowz = af[,1]
	af <- af[,2:(dim(af)[2])]
	af_mat = as.matrix(af)

	rownames(af_mat) =rowz
	colz <- colnames(af_mat)

	scaled_af_mat <- t(scale(t(af_mat)))


	ord <- hclust( dist(scaled_af_mat, method = "euclidean"), method = "ward.D" )$order

	melted_mat <- melt(scaled_af_mat)
	colnames(melted_mat) <- c("Exonic_site", "pseudotime","allelic_fraction")

	melted_mat$Exonic_site = factor(melted_mat$Exonic_site, levels=rownames(scaled_af_mat)[ord])
	melted_mat$pseudotime = factor(melted_mat$pseudotime, levels=colz)

	heatmap <- ggplot(data=melted_mat, aes(x=pseudotime, y=Exonic_site)) + geom_tile(aes(fill=allelic_fraction)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
	heatmap <- heatmap + labs(y="Exonic site", x="Pseudotime",fill="Standardized\nAllelic Fraction")
	heatmap <- heatmap + theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
	heatmap <- heatmap + theme(axis.text.x=element_blank())
	heatmap <- heatmap + theme(axis.text.y=element_blank())
	return(heatmap)
}

generate_U_af_heatmap <- function(allelic_fraction_pseutome_binned_file, component_num) {
	af <- read.table(allelic_fraction_pseutome_binned_file, header=TRUE, sep="\t")
	rowz = af[,1]
	af <- af[,2:(dim(af)[2])]
	af_mat = as.matrix(af)

	rownames(af_mat) =rowz
	colz <- colnames(af_mat)

	scaled_af_mat <- t(scale(t(af_mat)))


	ord <- hclust( dist(scaled_af_mat, method = "euclidean"), method = "ward.D" )$order

	melted_mat <- melt(scaled_af_mat)
	colnames(melted_mat) <- c("Exonic_site", "pseudotime","allelic_fraction")

	melted_mat$Exonic_site = factor(melted_mat$Exonic_site, levels=rownames(scaled_af_mat)[ord])
	melted_mat$pseudotime = factor(melted_mat$pseudotime, levels=colz)

	heatmap <- ggplot(data=melted_mat, aes(x=pseudotime, y=Exonic_site)) + geom_tile(aes(fill=allelic_fraction)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
	heatmap <- heatmap + labs(y="Exonic site", x=paste0("loading ", component_num),fill="Standardized\nAllelic Fraction")
	heatmap <- heatmap + theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
	heatmap <- heatmap + theme(axis.text.x=element_blank())
	heatmap <- heatmap + theme(axis.text.y=element_blank())
	return(heatmap)
}

generate_experiment_af_heatmap <- function(allelic_fraction_pseutome_binned_file) {
	af <- read.table(allelic_fraction_pseutome_binned_file, header=TRUE, sep="\t")
	rowz = af[,1]
	af <- af[,2:(dim(af)[2])]
	af_mat = as.matrix(af)

	rownames(af_mat) =rowz
	colz <- colnames(af_mat)

	scaled_af_mat <- t(scale(t(af_mat)))


	ord <- hclust( dist(scaled_af_mat, method = "euclidean"), method = "ward.D" )$order

	melted_mat <- melt(scaled_af_mat)
	colnames(melted_mat) <- c("Exonic_site", "pseudotime","allelic_fraction")

	melted_mat$Exonic_site = factor(melted_mat$Exonic_site, levels=rownames(scaled_af_mat)[ord])
	melted_mat$pseudotime = factor(melted_mat$pseudotime, levels=colz)

	heatmap <- ggplot(data=melted_mat, aes(x=pseudotime, y=Exonic_site)) + geom_tile(aes(fill=allelic_fraction)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
	heatmap <- heatmap + labs(y="Exonic site", x="Pseudotime",fill="Standardized\nAllelic Fraction")
	heatmap <- heatmap + theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
	heatmap <- heatmap + theme(axis.text.y=element_blank())
	return(heatmap)
}

generate_plate_af_heatmap <- function(allelic_fraction_pseutome_binned_file) {
	af <- read.table(allelic_fraction_pseutome_binned_file, header=TRUE, sep="\t")
	rowz = af[,1]
	af <- af[,2:(dim(af)[2])]
	af_mat = as.matrix(af)

	rownames(af_mat) =rowz
	colz <- colnames(af_mat)

	scaled_af_mat <- t(scale(t(af_mat)))


	ord <- hclust( dist(scaled_af_mat, method = "euclidean"), method = "ward.D" )$order

	melted_mat <- melt(scaled_af_mat)
	colnames(melted_mat) <- c("Exonic_site", "pseudotime","allelic_fraction")

	melted_mat$Exonic_site = factor(melted_mat$Exonic_site, levels=rownames(scaled_af_mat)[ord])
	melted_mat$pseudotime = factor(melted_mat$pseudotime, levels=colz)

	heatmap <- ggplot(data=melted_mat, aes(x=pseudotime, y=Exonic_site)) + geom_tile(aes(fill=allelic_fraction)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
	heatmap <- heatmap + labs(y="Exonic site", x="Plate id",fill="Standardized\nAllelic Fraction")
	heatmap <- heatmap + theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
	heatmap <- heatmap + theme(axis.text.y=element_blank()) + theme(axis.text.x=element_blank())
	return(heatmap)
}


generate_experiment_af_heatmap_cluster_experiments <- function(allelic_fraction_pseutome_binned_file) {
	af <- read.table(allelic_fraction_pseutome_binned_file, header=TRUE, sep="\t")
	rowz = af[,1]
	af <- af[,2:(dim(af)[2])]
	af_mat = as.matrix(af)

	rownames(af_mat) =rowz
	colz <- colnames(af_mat)

	scaled_af_mat <- t(scale(t(af_mat)))


	ord <- hclust( dist(scaled_af_mat, method = "euclidean"), method = "ward.D" )$order
	ord2 <- hclust( dist(t(scaled_af_mat), method = "euclidean"), method = "ward.D" )$order


	melted_mat <- melt(scaled_af_mat)
	colnames(melted_mat) <- c("Exonic_site", "pseudotime","allelic_fraction")

	melted_mat$Exonic_site = factor(melted_mat$Exonic_site, levels=rownames(scaled_af_mat)[ord])
	melted_mat$pseudotime = factor(melted_mat$pseudotime, levels=colz[ord2])

	heatmap <- ggplot(data=melted_mat, aes(x=pseudotime, y=Exonic_site)) + geom_tile(aes(fill=allelic_fraction)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
	heatmap <- heatmap + labs(y="Exonic site", x="Pseudotime",fill="Standardized\nAllelic Fraction")
	heatmap <- heatmap + theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
	heatmap <- heatmap + theme(axis.text.y=element_blank())
	return(heatmap)
}

make_fraction_biallelic_boxplot_by_experiment <- function(cell_info) {
	p<-ggplot(cell_info, aes(x=experiment, y=fraction_biallelic, color=experiment)) +
  		geom_boxplot() +
  		labs(x="Experiment", y="Fraction of biallelic sites",color="") + 
  		theme(legend.position = "none") +
  		theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
  	return(p)
}

make_fraction_biallelic_boxplot_by_plate <- function(cell_info) {
	p<-ggplot(cell_info, aes(x=factor(plate_id), y=fraction_biallelic, color=factor(plate_id))) +
  		geom_boxplot() +
  		labs(x="Plate ID", y="Fraction of biallelic sites",color="") + 
  		theme(legend.position = "none") +
  		theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) +
  		theme(axis.text.x=element_blank())
  	return(p)
}

make_fraction_biallelic_histogram <- function(cell_info) {
	p<-ggplot(cell_info, aes(x=fraction_biallelic)) +
  		geom_histogram() +
  		labs(x="Biallelic Fraction") + 
  		theme(legend.position = "none") +
  		theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
  		return(p)
}

make_fraction_biallelic_scatterplot <- function(cell_info) {
	p<-ggplot(cell_info, aes(x=number_biallelic_sites, y=num_expressed_sites)) +
  		geom_point(size=.001) +
  		labs(x="Number of biallelic sites", y="Number of sites") + 
  		theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
  		return(p)
}

make_fraction_biallelic_scatterplot_colored_by_experiment <- function(cell_info) {
	p<-ggplot(cell_info, aes(x=number_biallelic_sites, y=num_expressed_sites, color=experiment)) +
  		geom_point(size=.001) +
  		labs(x="Number of biallelic sites", y="Number of sites") + 
  		theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
  		return(p)
}

make_fraction_biallelic_scatterplot_colored_by_plate <- function(cell_info) {
	p<-ggplot(cell_info, aes(x=number_biallelic_sites, y=num_expressed_sites, color=factor(plate_id))) +
  		geom_point(size=.001) +
  		labs(x="Number of biallelic sites", y="Number of sites") + 
  		theme(legend.position="none") + 
  		theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
  		return(p)
}

make_fraction_biallelic_expression_pc1_scatterplot <- function(cell_info) {
	cell_info2 <- cell_info[cell_info$fraction_biallelic > .85, ]
	p<-ggplot(cell_info2, aes(x=fraction_biallelic, y=PC1_top500hvgs)) +
  		geom_point(size=.001) +
  		geom_smooth() +
  		labs(x="Biallelic fraction", y="Expression PC1") + 
  		theme(text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
  		return(p)
}

processed_data_dir = args[1]



cell_info_file <- paste0(processed_data_dir, "cell_info_after_filtering_0.3_0.5.txt")
allelic_fraction_pseutome_binned_file <- paste0(processed_data_dir, "ase_40_binned_by_pseudotime_allelic_fraction.txt")
allelic_fraction_experiment_binned_file <- paste0(processed_data_dir, "ase_binned_by_experiment_allelic_fraction.txt")
allelic_fraction_plate_id_binned_file <- paste0(processed_data_dir, "ase_binned_by_plate_id_allelic_fraction.txt")
allelic_fraction_U0_binned_file <- paste0(processed_data_dir, "ase_30_binned_by_U_0_allelic_fraction.txt")
allelic_fraction_U1_binned_file <- paste0(processed_data_dir, "ase_30_binned_by_U_1_allelic_fraction.txt")
allelic_fraction_U2_binned_file <- paste0(processed_data_dir, "ase_30_binned_by_U_2_allelic_fraction.txt")


cell_info <- read.table(cell_info_file, header=TRUE, sep="\t", comment.char="",quote="")

if (FALSE) {
binomial_p_mle = sum(cell_info$number_biallelic_sites)/sum(cell_info$num_expressed_sites)

biallelic_sites = cell_info$number_biallelic_sites
expressed_sites = cell_info$num_expressed_sites
num_sites = length(biallelic_sites)
for (site_num in 1:num_sites) {
	aa = binom.test(biallelic_sites[site_num], expressed_sites[site_num], p = binomial_p_mle, alternative = "two.sided",conf.level = 0.95)
	print(paste0(cell_info$fraction_biallelic[site_num], "\t", aa$p.value))
}
}

biallelic_scatter <- make_fraction_biallelic_expression_pc1_scatterplot(cell_info)
ggsave(biallelic_scatter, file=paste0(processed_data_dir, "fraction_biallelic_expression_pc1_scatter.pdf"), width=7.2, height=5.0, units="in")


biallelic_scatter <- make_fraction_biallelic_scatterplot(cell_info)
ggsave(biallelic_scatter, file=paste0(processed_data_dir, "fraction_biallelic_scatter.pdf"), width=7.2, height=5.0, units="in")


biallelic_scatter <- make_fraction_biallelic_scatterplot_colored_by_experiment(cell_info)
ggsave(biallelic_scatter, file=paste0(processed_data_dir, "fraction_biallelic_scatter_colored_by_experiment.pdf"), width=7.2, height=5.0, units="in")

biallelic_scatter <- make_fraction_biallelic_scatterplot_colored_by_plate(cell_info)
ggsave(biallelic_scatter, file=paste0(processed_data_dir, "fraction_biallelic_scatter_colored_by_plate_id.pdf"), width=7.2, height=5.0, units="in")



biallelic_histogram <- make_fraction_biallelic_histogram(cell_info)
ggsave(biallelic_histogram, file=paste0(processed_data_dir, "fraction_biallelic_histogram.pdf"), width=7.2, height=5.0, units="in")




experiment_fraction_biallelic_boxplot <- make_fraction_biallelic_boxplot_by_experiment(cell_info)
ggsave(experiment_fraction_biallelic_boxplot, file=paste0(processed_data_dir, "fraction_biallelic_boxplot_by_experiment.pdf"), width=7.2, height=5.0, units="in")

plate_fraction_biallelic_boxplot <- make_fraction_biallelic_boxplot_by_plate(cell_info)
ggsave(plate_fraction_biallelic_boxplot, file=paste0(processed_data_dir, "fraction_biallelic_boxplot_by_plate.pdf"), width=7.2, height=5.0, units="in")


U_0_heatmap <- generate_U_af_heatmap(allelic_fraction_U0_binned_file, 1)
ggsave(U_0_heatmap, file=paste0(processed_data_dir, "af_U1_binned_heatmap.pdf"), width=7.2, height=5.0, units="in")

U_1_heatmap <- generate_U_af_heatmap(allelic_fraction_U1_binned_file, 2)
ggsave(U_1_heatmap, file=paste0(processed_data_dir, "af_U2_binned_heatmap.pdf"), width=7.2, height=5.0, units="in")

U_2_heatmap <- generate_U_af_heatmap(allelic_fraction_U2_binned_file, 3)
ggsave(U_2_heatmap, file=paste0(processed_data_dir, "af_U3_binned_heatmap.pdf"), width=7.2, height=5.0, units="in")

combined_U_heatmap <- plot_grid(U_0_heatmap, U_1_heatmap, U_2_heatmap, ncol=1)
ggsave(combined_U_heatmap, file=paste0(processed_data_dir, "af_U_binned_heatmap_merged.pdf"), width=7.2, height=8.0, units="in")


plate_heatmap <- generate_plate_af_heatmap(allelic_fraction_plate_id_binned_file)
ggsave(plate_heatmap, file=paste0(processed_data_dir, "af_plate_binned_heatmap.pdf"), width=7.2, height=5.0, units="in")

if (FALSE) {
pseudotime_heatmap <- generate_pseudotime_af_heatmap(allelic_fraction_pseutome_binned_file)
ggsave(pseudotime_heatmap, file=paste0(processed_data_dir, "af_pseudotime_binned_heatmap.pdf"), width=7.2, height=5.0, units="in")

experiment_heatmap <- generate_experiment_af_heatmap(allelic_fraction_experiment_binned_file)
ggsave(experiment_heatmap, file=paste0(processed_data_dir, "af_experiment_binned_heatmap.pdf"), width=7.2, height=5.0, units="in")

experiment_heatmap <- generate_experiment_af_heatmap_cluster_experiments(allelic_fraction_experiment_binned_file)
ggsave(experiment_heatmap, file=paste0(processed_data_dir, "af_experiment_binned_clustered_heatmap.pdf"), width=7.2, height=5.0, units="in")
}
