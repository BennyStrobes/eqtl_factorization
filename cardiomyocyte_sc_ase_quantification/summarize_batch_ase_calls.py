import numpy as np 
import os
import sys
import pdb



def get_number_of_sites(full_file_name, min_number_reads):
	f = open(full_file_name)
	num_sites = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		total_reads = int(data[7])
		if total_reads >= min_number_reads:
			num_sites = num_sites + 1
	f.close()
	return num_sites


def get_number_of_biallelic_sites(full_file_name):
	f = open(full_file_name)
	num_sites = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		a1_reads = int(data[5])
		a2_reads = int(data[6])
		if a1_reads > 0 and a2_reads > 0:
			num_sites = num_sites + 1
	f.close()
	return num_sites

def get_cell_to_lib_size_mapping(input_file):
	f = open(input_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[3]
		lib_size = data[2]
		dicti[cell_id] = lib_size
	f.close()
	return dicti


ase_read_count_dir = sys.argv[1]
batch_name = sys.argv[2]


output_file = ase_read_count_dir + batch_name + 'ase_counts_summary.txt'

t = open(output_file, 'w')
# print header
t.write('cell_id\tnum_sites_ge_1\tnum_sites_ge_2\tnum_sites_ge_3\tnum_sites_ge_4\tnum_sites_ge_5\tnum_biallelic_sites\tlibrary_size\n')


cell_to_lib_size_mapping = get_cell_to_lib_size_mapping('/project2/gilad/bstrober/ase/cardiomyocyte_sc_ase_quantification/input_data/josh_libsizes.tsv')

# Loop through ase count for each cell
for file_name in os.listdir(ase_read_count_dir):
	if file_name.endswith('_read_counter_table.txt') == False:
		continue
	full_file_name = ase_read_count_dir + file_name
	# Get cell id
	cell_id = file_name.split('_')[0] + '_' + file_name.split('_')[1]
	# Extract info from data
	num_sites_1 = get_number_of_sites(full_file_name, 1)
	num_sites_2 = get_number_of_sites(full_file_name, 2)
	num_sites_3 = get_number_of_sites(full_file_name, 3)
	num_sites_4 = get_number_of_sites(full_file_name, 4)
	num_sites_5 = get_number_of_sites(full_file_name, 5)
	num_biallelic_sites = get_number_of_biallelic_sites(full_file_name)
	if cell_id not in cell_to_lib_size_mapping:
		lib_size = 'NA'
	else:
		lib_size = cell_to_lib_size_mapping[cell_id]
	# Print to output
	t.write(cell_id + '\t' + str(num_sites_1) + '\t' + str(num_sites_2) + '\t' + str(num_sites_3) + '\t' + str(num_sites_4) + '\t' + str(num_sites_5) + '\t' + str(num_biallelic_sites) + '\t' + str(lib_size) + '\n')
t.close()