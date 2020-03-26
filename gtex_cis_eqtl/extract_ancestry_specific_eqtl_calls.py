import numpy as np 
import os
import sys
import pdb
import pyarrow
import pandas as pd



def get_tissue_names(tissues_file):
	f = open(tissues_file)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		arr.append(data[0])
	f.close()
	return np.asarray(arr)


def get_test_names(test_names_file):
	f = open(test_names_file)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split(':')
		arr.append(data)
	f.close()
	return arr

def fill_in_effect_sizes(effect_sizes, test_names, input_data_root):
	for chrom_num in range(1,23):
		print(chrom_num)
		chrom_string = 'chr' + str(chrom_num) + '_'
		input_file = input_data_root + 'chr' + str(chrom_num) + '.parquet'
		df = pd.read_parquet(input_file, engine='pyarrow')
		for index, test_arr in enumerate(test_names):
			ensamble_id = test_arr[0]
			variant_id = test_arr[1]
			if variant_id.startswith(chrom_string) == False:
				continue
			try:
				slope = df[(df['phenotype_id'] == ensamble_id) & (df['variant_id'] == variant_id)]['slope']
				effect_sizes[index] = str(np.asarray(slope)[0])
			except:
				print("exception occured")
	return effect_sizes

def extract_ancestry_specific_eqtl_effect_sizes_for_a_specific_tissue(tissue_name, output_file, test_names, european_ancestry_gtex_eqtl_dir, african_ancestry_gtex_eqtl_dir):
	aa_effect_sizes = ['NaN']*len(test_names)
	ea_effect_sizes = ['NaN']*len(test_names)

	ea_effect_sizes = fill_in_effect_sizes(ea_effect_sizes, test_names, european_ancestry_gtex_eqtl_dir + tissue_name + '.v8.EUR.allpairs.')
	aa_effect_sizes = fill_in_effect_sizes(aa_effect_sizes, test_names, african_ancestry_gtex_eqtl_dir + tissue_name + '.v8.AFR.eqtl_allpairs.')
	

	t = open(output_file, 'w')
	t.write('ensamble_id\tvariant_id\tea_effect_size\taa_effect_size\n')
	for index, test_arr in enumerate(test_names):
		ensamble_id = test_arr[0]
		variant_id = test_arr[1]
		t.write(ensamble_id + '\t' + variant_id + '\t' + ea_effect_sizes[index] + '\t' + aa_effect_sizes[index] + '\n')
	t.close()


test_names_file = sys.argv[1]
tissues_file = sys.argv[2]
european_ancestry_gtex_eqtl_dir = sys.argv[3]
african_ancestry_gtex_eqtl_dir = sys.argv[4]
output_root = sys.argv[5]

tissue_names = get_tissue_names(tissues_file)

test_names = get_test_names(test_names_file)


for tissue_name in tissue_names:
	print(tissue_name)
	output_file = output_root + tissue_name + '_effect_sizes.txt'
	extract_ancestry_specific_eqtl_effect_sizes_for_a_specific_tissue(tissue_name, output_file, test_names, european_ancestry_gtex_eqtl_dir, african_ancestry_gtex_eqtl_dir)