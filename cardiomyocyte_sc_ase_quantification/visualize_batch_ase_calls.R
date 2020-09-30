args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(cowplot)
library(RColorBrewer)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



make_number_of_sites_histo <- function(num_sites, x_axis_label) {
	df <- data.frame(num_sites=num_sites)

	histo <- ggplot(df, aes(x=num_sites))+
 		geom_histogram(color="darkblue", fill="lightblue") +
 		figure_theme() + 
 		labs(y="Number of cells", x=x_axis_label)

 	return(histo)
}

make_lib_size_num_sites_scatter <- function(lib_size, num_sites, x_label, y_label) {
	df <- data.frame(num_sites=num_sites, library_size=lib_size)
	scatter <- ggplot(df, aes(x=library_size, y=num_sites)) + geom_point() +
		figure_theme() + 
		labs(x=x_label, y=y_label)
	return(scatter)
}




ase_read_count_dir <- args[1]
batch_name <- args[2]


ase_read_count_info_file <- paste0(ase_read_count_dir, batch_name, "ase_counts_summary.txt")
ase_read_count_info <- read.table(ase_read_count_info_file, header=TRUE)

# Lib size-num biallelic sites scatter
output_file <- paste0(ase_read_count_dir, batch_name, "_lib_size_num_biallelic_sites_scatter.pdf")
plotter <- make_lib_size_num_sites_scatter(ase_read_count_info$library_size, ase_read_count_info$num_biallelic_sites, "Library size", "Number of Biallelic sites")
ggsave(plotter, file=output_file, width=7.2, height=6.0, units="in")



# Make num-sites histogram
output_file <- paste0(ase_read_count_dir, batch_name, "_num_sites_with_ge_3_reads.pdf")
plotter <- make_number_of_sites_histo(ase_read_count_info$num_sites_ge_3, "Number of sites with ge 3 reads")
ggsave(plotter, file=output_file, width=7.2, height=6.0, units="in")

# Make num-sites histogram
output_file <- paste0(ase_read_count_dir, batch_name, "_num_sites_with_ge_4_reads.pdf")
plotter <- make_number_of_sites_histo(ase_read_count_info$num_sites_ge_4, "Number of sites with ge 4 reads")
ggsave(plotter, file=output_file, width=7.2, height=6.0, units="in")

# Make num-sites histogram
output_file <- paste0(ase_read_count_dir, batch_name, "_num_sites_with_ge_5_reads.pdf")
plotter <- make_number_of_sites_histo(ase_read_count_info$num_sites_ge_5, "Number of sites with ge 5 reads")
ggsave(plotter, file=output_file, width=7.2, height=6.0, units="in")

# Make num-sites histogram
output_file <- paste0(ase_read_count_dir, batch_name, "_num_biallelic_sites.pdf")
plotter <- make_number_of_sites_histo(ase_read_count_info$num_biallelic_sites, "Number of bi-allelic sites")
ggsave(plotter, file=output_file, width=7.2, height=6.0, units="in")
