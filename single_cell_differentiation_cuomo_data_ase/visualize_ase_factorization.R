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


make_loading_expression_pc_scatter_for_each_factor <- function(expression_pc, loading, factor_number, x_axis_label) {
	df <- data.frame(expression_pc=expression_pc, loading=loading)
	print(factor_number)
	print(cor(expression_pc, loading)*cor(expression_pc, loading))

	plotter <- ggplot(df, aes(x=expression_pc, y=loading)) + 
	           geom_point(size=.001, alpha=.35) +
	           gtex_v8_figure_theme() + 
	           labs(x=x_axis_label, y = paste0("Loading ", factor_number)) + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) +
	           geom_smooth()
	return(plotter)

}

make_loading_cov_cov_scatter_colored_by_loadings_for_each_factor <- function(cov1, cov2, loading, factor_number, x_axis_label, y_axis_label) {
	df <- data.frame(cov1=cov1, cov2=cov2, loading=loading)
	plotter <- ggplot(df, aes(x=cov1, y=cov2, color=loading)) + 
	           geom_point(size=.001) +
	           gtex_v8_figure_theme() + 
	           labs(x=x_axis_label, y = y_axis_label) + 
	           scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) 
	return(plotter)

}

loading_expression_pc1_scatter_with_row_for_every_factor <- function(expression_pc, loadings, x_axis_label) {
	factor_number <- 1
	factor_1_scatterplot <- make_loading_expression_pc_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)


	factor_number <- 2
	factor_2_scatterplot <- make_loading_expression_pc_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)


	factor_number <- 3
	factor_3_scatterplot <- make_loading_expression_pc_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)

	factor_number <- 4
	#factor_4_scatterplot <- make_loading_expression_pc_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)

	factor_number <- 5
	#factor_5_scatterplot <- make_loading_expression_pc_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)



	combined <- plot_grid(factor_1_scatterplot, factor_2_scatterplot, factor_3_scatterplot, ncol=1)


	return(combined)
}

loading_cov_cov_scatter_colored_by_loadings_with_row_for_every_factor <- function(cov1, cov2, loadings, x_axis_label, y_axis_label) {
	factor_number <- 1
	factor_1_scatterplot <- make_loading_cov_cov_scatter_colored_by_loadings_for_each_factor(cov1, cov2, loadings[, factor_number], factor_number, x_axis_label, y_axis_label)


	factor_number <- 2
	factor_2_scatterplot <- make_loading_cov_cov_scatter_colored_by_loadings_for_each_factor(cov1, cov2, loadings[, factor_number], factor_number, x_axis_label, y_axis_label)


	factor_number <- 3
	factor_3_scatterplot <- make_loading_cov_cov_scatter_colored_by_loadings_for_each_factor(cov1, cov2, loadings[, factor_number], factor_number, x_axis_label, y_axis_label)

	factor_number <- 4
	#factor_4_scatterplot <- make_loading_expression_pc_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)

	factor_number <- 5
	#factor_5_scatterplot <- make_loading_expression_pc_scatter_for_each_factor(expression_pc, loadings[, factor_number], factor_number, x_axis_label)



	combined <- plot_grid(factor_1_scatterplot, factor_2_scatterplot, factor_3_scatterplot, ncol=1)


	return(combined)

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

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_go_term_loading_loading_correlation_heatmap <- function(covariates, loadings) {
	# Covariates columns to consider
	#print(summary(covariates[, valid_covariates]))
	# Remove unimportant columns

    loadings <- as.matrix(loadings)
    print(summary(covariates))
    #valid_covariates <- 2:85
    #covs <- covariates[,valid_covariates]
    #valid_covariates <- c(6, 7, 8, 10, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 36, 37, 38, 40, 41, 43, 44, 45,46, 47, 48, 49, 51, 52, 53, 54, 59, 60, 69, 71, 87, 95, 96, 97,100)
 	#covariate_type <- c("num", "cat", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "cat", "num", "num", "num", "num", "num", "num", "num")

 	cov_names <- colnames(covariates)
 	num <- length(cov_names)



  	covs <- as.matrix(covariates)
	#covs <- covariates


    # Initialize PVE heatmap
    factor_colnames <- paste0("Factor", 1:(dim(loadings)[2]))
    factor_rownames <- colnames(covs)
    pve_map <- matrix(0, dim(covs)[2], dim(loadings)[2])
    print(dim(pve_map))
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
    print(which(colnames(covs) == "AEROBIC_RESPIRATION"))
    index = which(colnames(covs) == "AEROBIC_RESPIRATION")
    print(pve_map[index,])
    print(which.max(pve_map[,1]))
    print(pve_map[which.max(pve_map[,1]),])
    print(which.max(pve_map[,2]))
    print(pve_map[which.max(pve_map[,2]),])
    print(which.max(pve_map[,3]))
    print(pve_map[which.max(pve_map[,3]),])


    
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
    heatmap <- heatmap + labs(y="", x="GO Term Loading",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    heatmap <- heatmap + theme(axis.text.x=element_blank())
    # Save File
    return(heatmap)

}


######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
	# Covariates columns to consider
	#print(summary(covariates[, valid_covariates]))
	# Remove unimportant columns
	covariates$experiment_day = factor(paste0(covariates$experiment, "_", covariates$day))
	covariates$frac_alt_reads = covariates$n_alt_reads/covariates$n_total_reads
    loadings <- as.matrix(loadings)
    print(summary(covariates))
    #valid_covariates <- 2:85
    #covs <- covariates[,valid_covariates]
    valid_covariates <- c(6, 7, 8, 10, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 36, 37, 38, 40, 41, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,55, 56, 57, 58, 59, 60, 69, 71, 77, 78, 79, 80, 81,83, 87, 95, 96, 97,98,99,100, 101,102,103, 104, 105, 106, 107, 108, 109, 110, 111, 112)
 	covariate_type <- c("num", "cat", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "cat", "num", "num", "num", "num", "num", "num", "num", "cat", "num", "num", "num", "num", "num", "num", "cat","num","num","num", "num", "num", "num", "num", "num", "num", "num", "cat", "num")

 	num_cov = length(valid_covariates)
 	for (iter in 1:num_cov) {
 		print(paste0(colnames(covariates)[valid_covariates[iter]], " ", covariate_type[iter]))
 	}
 	print(length(valid_covariates))
 	print(length(covariate_type))
 	print(colnames(covariates))
 	cov_names <- colnames(covariates)[valid_covariates]
 	num <- length(cov_names)
 	print(cov_names)



  	covs <- covariates[,valid_covariates]
	#covs <- covariates


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
            if (covariate_type[num_cov] == "cat") {
            #print(cov_vec[1:10])
            	lin_model <- lm(pc_vec ~ factor(cov_vec))
        	} else {
        		lin_model <- lm(pc_vec ~ cov_vec)
        	}
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    
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

make_loading_boxplot_for_one_factor_by_categorical_covariate <- function(covariate, loadings, factor_number, covariate_name) {
	df <- data.frame(loading=loadings, covariate=factor(covariate))

	boxplot <- ggplot(df, aes(x=covariate, y=loading, fill=covariate)) + geom_boxplot(outlier.size = .00001) +
				gtex_v8_figure_theme() + 
	        	labs(x="", y = paste0("Sample loading (", factor_number,")"), fill="") +
	        	theme(legend.position="none") +
	        	guides(colour = guide_legend(override.aes = list(size=2))) +
	        	theme(axis.text.x=element_blank()) + 
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) +
	           	geom_hline(yintercept=0)


}
make_loading_boxplot_for_one_factor_by_2_categorical_covariates <- function(covariate1, covariate2, loadings, factor_number, covariate_name) {
	df <- data.frame(loading=loadings, covariate1=factor(covariate1), covariate2=factor(covariate2))

	boxplot <- ggplot(df, aes(x=covariate1, y=loading, fill=covariate2)) + geom_boxplot(outlier.size = .00001) +
				gtex_v8_figure_theme() + 
	        	labs(x="", y = paste0("Sample loading (", factor_number,")"), fill="") +
	        	theme(legend.position="none") +
	        	guides(colour = guide_legend(override.aes = list(size=2))) +
	        	theme(axis.text.x=element_blank()) + 
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) +
	           	geom_hline(yintercept=0)


}


make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate <- function(covariate, loadings, covariate_name) {
	loading_vec <- c()
	covariate_vec <- c()
	num_factors <- dim(loadings)[2]
	print(num_factors)

	factor_number <- 1
	factor_1_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 2
	factor_2_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)


	factor_number <- 3
	factor_3_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 4
	#factor_4_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 5
	#factor_5_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	#factor_number <- 4
	#factor_4_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	#factor_number <- 5
	#factor_5_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)


	#combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, factor_4_boxplot, factor_5_boxplot, ncol=1)
	combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, ncol=1)

	return(combined)
}

make_loading_boxplot_plot_with_row_for_every_factor_by_2_categorical_covariates <- function(covariate1, covariate2, loadings, covariate_name) {
	loading_vec <- c()
	covariate_vec <- c()
	num_factors <- dim(loadings)[2]
	print(num_factors)

	factor_number <- 1
	factor_1_boxplot <- make_loading_boxplot_for_one_factor_by_2_categorical_covariates(covariate1, covariate2, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 2
	factor_2_boxplot <- make_loading_boxplot_for_one_factor_by_2_categorical_covariates(covariate1, covariate2, loadings[,factor_number], factor_number, covariate_name)


	factor_number <- 3
	factor_3_boxplot <- make_loading_boxplot_for_one_factor_by_2_categorical_covariates(covariate1, covariate2, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 4
	#factor_4_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	factor_number <- 5
	#factor_5_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	#factor_number <- 4
	#factor_4_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)

	#factor_number <- 5
	#factor_5_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_number], factor_number, covariate_name)


	#combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, factor_4_boxplot, factor_5_boxplot, ncol=1)
	combined <- plot_grid(factor_1_boxplot, factor_2_boxplot, factor_3_boxplot, ncol=1)

	return(combined)
}


#############################
# Command line args
#############################
eqtl_results_dir <- args[1]
processed_data_dir <- args[2]
covariate_file <- args[3]
go_loadings_file <- args[4]
visualization_dir <- args[5]



# Input files
#covariate_file <- paste0(processed_data_dir, "annotated_sample_names.txt")
covariates <- read.table(covariate_file, header=TRUE, sep="\t", comment.char="",quote="")


#go_loadings <- read.table(go_loadings_file, header=TRUE, sep="\t", comment.char="", quote="")

loading_file <- paste0(eqtl_results_dir, "ase_factorization_via_pca_non_min_counts_regress_out_cell_line_subsampled_high_biallelic_fraction_only_endoderm_differentiation_3_ase_factorization3_temper_U.txt")
loadings <- read.table(loading_file, header=FALSE)


output_stem <- "endoderm_differentiation_pca_non_min_k_3_"

# Create UMAP factors
umap_loadings = umap(loadings)$layout
saveRDS( umap_loadings, "umap_loadings.rds")
#print("UMAP DONE")
#umap_loadings <- readRDS("umap_loadings.rds")


######################################
# Visualize UMAP scatter plot colored by tissue type
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_tissue_type.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$tissue_id, umap_loadings, "Known tissue type")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")



######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
output_file <- paste0(visualization_dir, output_stem, "covariate_loading_correlation_heatmap.pdf")
#heatmap <- make_covariate_loading_correlation_heatmap(covariates, loadings)
#ggsave(heatmap, file=output_file, width=18.2, height=10, units="in")


######################################
# Visualize UMAP scatter plot colored by day
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_day.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$day, umap_loadings, "Day")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "boxplot_colored_by_plate_id.pdf")
plate_id_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(factor(covariates$plate_id), loadings, "plate_id")
ggsave(plate_id_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "boxplot_colored_by_plate_id.pdf")
plate_id_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(factor(covariates$plate_id), loadings, "plate_id")
ggsave(plate_id_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "boxplot_colored_by_indi_id.pdf")
indi_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$donor_long_id, loadings, "Individual")
ggsave(indi_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "boxplot_colored_by_experiment_day.pdf")
exp_day_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_2_categorical_covariates(factor(paste0(covariates$experiment,"_",covariates$day)), factor(covariates$day), loadings, "Experiment_day")
ggsave(exp_day_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "boxplot_colored_by_experiment_day2.pdf")
exp_day_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_2_categorical_covariates(factor(paste0(covariates$experiment,"_",covariates$day)), factor(covariates$experiment), loadings, "Experiment_day")
ggsave(exp_day_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "boxplot_colored_by_experiment.pdf")
experiment_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$experiment, loadings, "Experiment")
ggsave(experiment_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "boxplot_colored_by_day.pdf")
experiment_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$day, loadings, "Day")
ggsave(experiment_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")


output_file <- paste0(visualization_dir, output_stem, "boxplot_colored_by_cell_cycle_state.pdf")
cc_categorical_boxplot = make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$cc_phase, loadings, "Cell cycle state")
ggsave(cc_categorical_boxplot, file=output_file, width=7.2, height=7.0, units="in")


######################################
# Make Heatmap correlating known covariates with eqtl factorization loadings
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_covariate_heatmap.pdf")
heatmap <- make_covariate_loading_correlation_heatmap(covariates, loadings) 
ggsave(heatmap, file=output_file, width=10.2, height=5.5, units="in")


######################################
# Make Heatmap correlating known covariates with eqtl factorization loadings
#######################################
output_file <- paste0(visualization_dir, "loading_go_term_loading_heatmap.pdf")
#heatmap <- make_go_term_loading_loading_correlation_heatmap(go_loadings[,2:(dim(go_loadings)[2])], loadings) 
#ggsave(heatmap, file=output_file, width=25.2, height=5.5, units="in")
print("DONE")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_expression_S_score.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$S_score, umap_loadings, "S-score")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by G2_M_transition
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_G2_M_transition.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$G2_M_transition, umap_loadings, "G2_M_transition")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_expression_G2M_score.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$G2M_score, umap_loadings, "G2M-score")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_expression_pc1.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$PC1_top500hvgs, umap_loadings, "expression pc1")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_fraction_biallelic.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$fraction_biallelic, umap_loadings, "fraction biallelic")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_log10_total_features.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$log10_total_features, umap_loadings, "log10(total features)")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_log10_total_counts.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$log10_total_counts, umap_loadings, "log10_total_counts")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_pct_counts_top_200_features_endogenous.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$pct_counts_top_200_features_endogenous, umap_loadings, "pct_counts_top_200_features_endogenous")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_pct_counts_top_500_features_endogenous.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$pct_counts_top_500_features_endogenous, umap_loadings, "pct_counts_top_500_features_endogenous")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_differentiation.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$differentiation, umap_loadings, "differentiation")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
aa = covariates$pct_counts_top_500_features_endogenous
print(summary(aa))
aa[aa < 52] = 52
aa[aa > 65] = 65
print(summary(aa))
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_pct_counts_top_500_features_endogenous_thresh.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(aa, umap_loadings, "pct_counts_top_500_features_endogenous")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_pct_counts_top_200_features.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$pct_counts_top_200_features, umap_loadings, "pct_counts_top_200_features")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_pct_counts_top_500_features.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$pct_counts_top_500_features, umap_loadings, "pct_counts_top_500_features")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell line
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_cell_cycle_phase.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$cc_phase, umap_loadings, "Cell cycle")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_experiment.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$experiment, umap_loadings, "Experiment")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell line
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_cell_line.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$donor_long_id, umap_loadings, "cell_line")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell line
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_plate.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$plate_id, umap_loadings, "cell_line")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")




######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "pseudotime_G2_M_transition_scatter_with_row_for_every_factor_colored_by_loadings.pdf")
scatterplot <- loading_cov_cov_scatter_colored_by_loadings_with_row_for_every_factor(covariates$pseudo, covariates$G2_M_transition, loadings, "Pseudotime", "G2_M_transition")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "pseudotime_G1_S_transition_scatter_with_row_for_every_factor_colored_by_loadings.pdf")
scatterplot <- loading_cov_cov_scatter_colored_by_loadings_with_row_for_every_factor(covariates$pseudo, covariates$G1_S_transition, loadings, "Pseudotime", "G1_S_transition")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "pseudotime_differentiation_scatter_with_row_for_every_factor_colored_by_loadings.pdf")
scatterplot <- loading_cov_cov_scatter_colored_by_loadings_with_row_for_every_factor(covariates$pseudo, covariates$differentiation, loadings, "Pseudotime", "Differentiation")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "pseudotime_respiration_scatter_with_row_for_every_factor_colored_by_loadings.pdf")
scatterplot <- loading_cov_cov_scatter_colored_by_loadings_with_row_for_every_factor(covariates$pseudo, covariates$respiration, loadings, "Pseudotime", "respiration")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "pseudotime_sterol_biosynthesis_scatter_with_row_for_every_factor_colored_by_loadings.pdf")
scatterplot <- loading_cov_cov_scatter_colored_by_loadings_with_row_for_every_factor(covariates$pseudo, covariates$sterol_biosynthesis, loadings, "Pseudotime", "sterol_biosynthesis")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_expression_pc1_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$PC1_top500hvgs , loadings, "Expression PC1")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_pseudotime2_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$pseudotime2 , loadings, "Pseudotime")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_pseudoXG2_M_transition_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$pseudo*covariates$G2_M_transition, loadings, "pseudotimeXG2_M_transition")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "loading_pseudoXsterol_biosynthesis_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$pseudo*covariates$sterol_biosynthesis, loadings, "pseudotimeXsterol_biosynthesis")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

output_file <- paste0(visualization_dir, output_stem, "loading_G2_M_transitionXsterol_biosynthesis_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$G2_M_transition*covariates$sterol_biosynthesis, loadings, "G2_M_transitionXsterol_biosynthesis")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_pseudo_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$pseudo , loadings, "Pseudotime")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_fraction_biallelic_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$fraction_biallelic , loadings, "Fraction biallelic")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")


######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_log10_total_feature_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$log10_total_features , loadings, "log10(total features)")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_G2_M_transition_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$G2_M_transition, loadings, "G2_M_transition")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_differentiation_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$differentiation, loadings, "differentiation")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")


######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_fraction_biallelic_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$fraction_biallelic, loadings, "fraction_biallelic")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_number_biallelic_sites_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$number_biallelic_sites, loadings, "number_biallelic_sites")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_global_allelic_fraction_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$global_allelic_fraction, loadings, "global_allelic_fraction")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_G1_S_transition_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$G1_S_transition, loadings, "G1_S_transition")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_respiration_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$respiration, loadings, "respiration")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_sterol_biosynthesis_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$sterol_biosynthesis, loadings, "sterol_biosynthesis")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_S_score_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$S_score , loadings, "S_score")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_G2M_score_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$G2M_score , loadings, "G2M_score")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_pct_counts_top_200_features_endogenous_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$pct_counts_top_200_features_endogenous, loadings, "pct_counts_top_200_features_endogenous")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_pct_counts_top_500_features_endogenous_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$pct_counts_top_500_features_endogenous, loadings, "pct_counts_top_500_features_endogenous")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_pct_counts_top_200_features_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$pct_counts_top_200_features, loadings, "pct_counts_top_200_features")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make loading scatterplot with row for every factor colored by expression PC1
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_pct_counts_top_500_features_scatter_with_row_for_every_factor.pdf")
scatterplot <- loading_expression_pc1_scatter_with_row_for_every_factor(covariates$pct_counts_top_500_features, loadings, "pct_counts_top_500_features_endogenous")
ggsave(scatterplot, file=output_file, width=7.2, height=9.0, units="in")


































######################################
# Make loading boxplot with row for every factor colored by day
#######################################
output_file <- paste0(visualization_dir, output_stem, "loading_boxplot_with_row_for_every_factor_colored_by_day.pdf")
#boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(factor(covariates$day), covariates$day, loadings, "Day")
#ggsave(boxplot, file=output_file, width=7.2, height=6.5, units="in")

######################################
# Visualize UMAP scatter plot colored by COHORT
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_COHORT.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$COHORT, umap_loadings, "COHORT")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by DTHLUCOD
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_DTHLUCOD.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$DTHLUCOD, umap_loadings, "Cause of death")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by MHBLDDNDR
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_MHBLDDNDR.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$MHBLDDNDR, umap_loadings, "Blood donor")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")



######################################
# Visualize UMAP scatter plot colored by DTHPLCE
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_DTHPLCE.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$DTHPLCE, umap_loadings, "Death Place")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")



######################################
# Visualize UMAP scatter plot colored by TRISCHD
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_TRISCHD.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$TRISCHD, umap_loadings, "TRISCHD")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by Race
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_race.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$RACE, umap_loadings, "Race")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by Height
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_height.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$HGHT, umap_loadings, "Height")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by Age
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_Age.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$AGE, umap_loadings, "AGE")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by BMI
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_BMI.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$BMI, umap_loadings, "BMI")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by fibroblase counts
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_endothelial_cell_percentage.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$endothelial.cell, umap_loadings, "Endothelial cell %")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by fibroblase counts
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_fibroblast_cell_percentage.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$fibroblast, umap_loadings, "Fibroblast cell %")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by endocardial counts
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_endocardial_cell_percentage.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$endocardial.cell, umap_loadings, "Endocardial cell %")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cardiomyocyte counts
#######################################
output_file <- paste0(visualization_dir, output_stem, "umap_loading_scatter_colored_by_cardiac_muscle_cell_percentage.pdf")
#umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$cardiac.muscle.cell, umap_loadings, "Cardiac muscle cell %")
#ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


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





print("Hello")

frac_expressed = "0.1"
# Input files
covariate_file <- paste0(processed_expression_dir, "cell_covariates_sle_individuals_random_subset_min_expressed_cells_0.05_log_transform_transform.txt")
#covariate_file <- paste0(processed_expression_dir, "cell_covariates_sle_individuals.txt")
#covariate_file <- paste0(processed_expression_dir, "pseudobulk_covariates_sle_individuals.txt")

#eqtl_factorization_loading_file <- paste0(eqtl_results_dir, "eqtl_factorization_pseudobulk_data_20_factors_eqtl_factorization_vi_spike_and_slab_model_True_re_False_svi_0_seed_U_S.txt")
#eqtl_factorization_loading_file <- "/home-1/bstrobe1@jhu.edu/work/ben/temp/temper_U_als.txt"
#eqtl_factorization_loading_file <- paste0(eqtl_results_dir, "eqtl_factorization_single_cell_sig_tests_50_pc_min_expressed_cells_", frac_expressed, "_data_uncorrected_genotype_5_factors_eqtl_factorization_als_model_False_re_False_svi_0_seed_U.txt")
eqtl_factorization_loading_file <- "/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/eqtl_factorization_results/temper2_U_S.txt"

# Load in data
covariates <- read.table(covariate_file, header=TRUE, sep="\t")
loadings <- read.table(eqtl_factorization_loading_file, header=FALSE)
# Filter loadings file

#loadings <- loadings[,good_loadings]

# Create UMAP factors
#umap_loadings = umap(loadings)$layout
#saveRDS( umap_loadings, "umap_loadings.rds")
print("UMAP DONE")
umap_loadings <- readRDS("umap_loadings.rds")


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
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$batch_cov, umap_loadings, "Known batch")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

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
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$percent_mito, umap_loadings, "Percent mito")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by known cell counts
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_cell_counts.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$n_counts, umap_loadings, "Cell counts")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by known number of genes
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_number_of_genes.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$n_genes, umap_loadings, "Number of genes")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
}

