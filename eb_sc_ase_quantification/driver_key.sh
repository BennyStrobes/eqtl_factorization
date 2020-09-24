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
root_dir="/project2/gilad/bstrober/ase/eb_sc_ase_quantification/"

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


bam_input_directory="/project2/gilad/kenneth/HumanChimpEBsOrtho/OrthoCellranger/Batch1/Batch1_Lane1_human/"
batch_name="batch1_lane1"
sh align_one_batch.sh $batch_name $bam_input_directory $ase_read_count_dir