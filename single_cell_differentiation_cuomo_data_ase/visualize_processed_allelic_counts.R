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


pseudotime_binned_file <- paste0(processed_data_dir,"ase_30_binned_by_pseudotime_allelic_fraction.txt")

pseudotime_heatmap <- generate_pseudotime_af_heatmap(pseudotime_binned_file)
ggsave(pseudotime_heatmap, file=paste0(processed_data_dir, "pseudotime_heatmap.pdf"), width=7.2, height=5.0, units="in")


binomial_allelic_fraction_U0_binned_file <- paste0(processed_data_dir, "ase_binomial_mean_30_binned_by_mixture_U_0_allelic_fraction.txt")
binomial_allelic_fraction_U1_binned_file <- paste0(processed_data_dir, "ase_binomial_mean_30_binned_by_mixture_U_1_allelic_fraction.txt")
binomial_allelic_fraction_U2_binned_file <- paste0(processed_data_dir, "ase_binomial_mean_30_binned_by_mixture_U_2_allelic_fraction.txt")


if (FALSE) {
U_0_heatmap <- generate_U_af_heatmap(binomial_allelic_fraction_U0_binned_file, 1)
U_1_heatmap <- generate_U_af_heatmap(binomial_allelic_fraction_U1_binned_file, 2)
U_2_heatmap <- generate_U_af_heatmap(binomial_allelic_fraction_U2_binned_file, 3)
combined_U_heatmap <- plot_grid(U_0_heatmap, U_1_heatmap, U_2_heatmap, ncol=1)
ggsave(combined_U_heatmap, file=paste0(processed_data_dir, "af_binomial_fraction_U_binned_heatmap_merged.pdf"), width=7.2, height=8.0, units="in")
print('done')
folded_binomial_allelic_fraction_U0_binned_file <- paste0(processed_data_dir, "ase_folded_binomial_mean_30_binned_by_mixture_U_0_allelic_fraction.txt")
folded_binomial_allelic_fraction_U1_binned_file <- paste0(processed_data_dir, "ase_folded_binomial_mean_30_binned_by_mixture_U_1_allelic_fraction.txt")
folded_binomial_allelic_fraction_U2_binned_file <- paste0(processed_data_dir, "ase_folded_binomial_mean_30_binned_by_mixture_U_2_allelic_fraction.txt")

U_0_heatmap <- generate_U_af_heatmap(folded_binomial_allelic_fraction_U0_binned_file, 1)
U_1_heatmap <- generate_U_af_heatmap(folded_binomial_allelic_fraction_U1_binned_file, 2)
U_2_heatmap <- generate_U_af_heatmap(folded_binomial_allelic_fraction_U2_binned_file, 3)
combined_U_heatmap <- plot_grid(U_0_heatmap, U_1_heatmap, U_2_heatmap, ncol=1)
ggsave(combined_U_heatmap, file=paste0(processed_data_dir, "af_folded_binomial_fraction_U_binned_heatmap_merged.pdf"), width=7.2, height=8.0, units="in")
print('done')

}
beta_binomial_allelic_fraction_U0_binned_file <- paste0(processed_data_dir, "ase_beta_binomial_mean_30_binned_by_mixture_U_0_allelic_fraction.txt")
beta_binomial_allelic_fraction_U1_binned_file <- paste0(processed_data_dir, "ase_beta_binomial_mean_30_binned_by_mixture_U_1_allelic_fraction.txt")
beta_binomial_allelic_fraction_U2_binned_file <- paste0(processed_data_dir, "ase_beta_binomial_mean_30_binned_by_mixture_U_2_allelic_fraction.txt")


U_0_heatmap <- generate_U_af_heatmap(beta_binomial_allelic_fraction_U0_binned_file, 1)
U_1_heatmap <- generate_U_af_heatmap(beta_binomial_allelic_fraction_U1_binned_file, 2)
U_2_heatmap <- generate_U_af_heatmap(beta_binomial_allelic_fraction_U2_binned_file, 3)
combined_U_heatmap <- plot_grid(U_0_heatmap, U_1_heatmap, U_2_heatmap, ncol=1)
ggsave(combined_U_heatmap, file=paste0(processed_data_dir, "af_beta_binomial_fraction_U_binned_heatmap_merged.pdf"), width=7.2, height=8.0, units="in")
print('done')

folded_beta_binomial_allelic_fraction_U0_binned_file <- paste0(processed_data_dir, "ase_folded_beta_binomial_mean_30_binned_by_mixture_U_0_allelic_fraction.txt")
folded_beta_binomial_allelic_fraction_U1_binned_file <- paste0(processed_data_dir, "ase_folded_beta_binomial_mean_30_binned_by_mixture_U_1_allelic_fraction.txt")
folded_beta_binomial_allelic_fraction_U2_binned_file <- paste0(processed_data_dir, "ase_folded_beta_binomial_mean_30_binned_by_mixture_U_2_allelic_fraction.txt")

U_0_heatmap <- generate_U_af_heatmap(folded_beta_binomial_allelic_fraction_U0_binned_file, 1)
U_1_heatmap <- generate_U_af_heatmap(folded_beta_binomial_allelic_fraction_U1_binned_file, 2)
U_2_heatmap <- generate_U_af_heatmap(folded_beta_binomial_allelic_fraction_U2_binned_file, 3)
combined_U_heatmap <- plot_grid(U_0_heatmap, U_1_heatmap, U_2_heatmap, ncol=1)
ggsave(combined_U_heatmap, file=paste0(processed_data_dir, "af_folded_beta_binomial_fraction_U_binned_heatmap_merged.pdf"), width=7.2, height=8.0, units="in")

print('done')


