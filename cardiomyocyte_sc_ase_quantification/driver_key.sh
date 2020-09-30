#!/bin/bash -l

#SBATCH
#SBATCH --time=14:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=10GB
#SBATCH --nodes=1



#########################
# Output directories
#########################
# Root directory for simulation analysis
root_dir="/project2/gilad/bstrober/ase/cardiomyocyte_sc_ase_quantification/"

# Directory containing genotype info
input_data_dir=$root_dir"input_data/"

# Directory containing genotype info
genotype_dir=$root_dir"genotype/"

# Directory containing ASE read counts
ase_read_count_dir=$root_dir"ase_read_counts/"


#########################
# Input Data
#########################
# File containing list of input directories (one for each batch, lane)
bam_input_directories_file=$input_data_dir"bam_input_directories.txt"

# VCF File
vcf_file="/project2/gilad/bstrober/ase/YRI_input_data/YRI_genotype_hg38_liftover.vcf.gz"

# Master barcode file
master_barcode_file="/project2/gilad/reem/singlecellCM/sampleandcol_metadata_allcells_xfinalprb_postfilter.txt"

# Genome build
genome_build_file="/project2/gilad/bstrober/ase/YRI_input_data/hg38.fa"

# GATK jar file
gatk_jar="/home/bstrober/ipsc_differentiation/preprocess/v5/scripts/GenomeAnalysisTK.jar"

# PICARD Jar file
picard_jar="/home/bstrober/ipsc_differentiation/preprocess/v5/scripts/picard.jar"



bam_input_directory="/project2/gilad/reem/singlecellCM/round1/xfinalmerged/CD2/CD2col3/output/"
batch_name="E1CD2col3"
sh align_one_batch.sh $batch_name $bam_input_directory $vcf_file $master_barcode_file $genome_build_file $gatk_jar $picard_jar $ase_read_count_dir



