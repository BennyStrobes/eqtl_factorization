import numpy as np 
import os
import sys
import pdb


def get_ase_sites_dictionary(ase_counts_file):
	f = open(ase_counts_file)
	ase_sites = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		full_id = data[0]
		short_id = full_id.split('_')[0] + ':' + (full_id.split('_')[1])
		ase_sites[short_id] = 0.0
	f.close()
	return ase_sites

def get_ase_site_to_genotype_vec_mapping(ase_sites, genotype_file):
	f = open(genotype_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 106:
			print('assumptoine eroro')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			ordered_individuals = np.asarray(data[1:])
			continue
		site_id = data[0]
		if site_id not in ase_sites:
			continue
		genotype_vec = np.asarray(data[1:]).astype(float)
		dicti[site_id] = genotype_vec
	f.close()
	return dicti, ordered_individuals


ase_counts_file = sys.argv[1]
cell_info_file = sys.argv[2]
genotype_file = sys.argv[3]
processed_data_dir = sys.argv[4]

thresh = 10

# First get names of sites with ase data for
ase_sites = get_ase_sites_dictionary(ase_counts_file)
# Then create mapping from site (ie variant id) to genotype vector across individuals
site_to_genotype_vec, ordered_individuals = get_ase_site_to_genotype_vec_mapping(ase_sites, genotype_file)
# Then for each cell, compute fraction of expressed sites (>=10 reads) where each individual gets it "Right"

head_count = 0
counter = 0
f = open(ase_counts_file)
for line in f:
	print(counter)
	counter = counter + 1
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		cell_ids = np.asarray(data[1:])
		numerz = np.zeros((len(ordered_individuals), len(cell_ids)))
		denerz = np.zeros((len(ordered_individuals), len(cell_ids)))
		continue
	full_id = data[0]
	short_id = full_id.split('_')[0] + ':' + (full_id.split('_')[1])
	if short_id not in site_to_genotype_vec:
		continue
	allelic_counts = np.asarray(data[1:])
	genotype_vec = site_to_genotype_vec[short_id]

	for cell_index, allelic_count in enumerate(allelic_counts):
		if allelic_count == 'NA':
			continue
		numer = int(allelic_count.split('/')[0])
		dener = int(allelic_count.split('/')[1])
		if dener < thresh:
			continue
		for indi_index, individual in enumerate(ordered_individuals):
			indi_genotype = np.round(genotype_vec[indi_index])
			denerz[indi_index, cell_index] = denerz[indi_index, cell_index] + 1
			# If monoallelic
			if numer == dener or numer == 0:
				# If homozygous at snp
				if indi_genotype == 0.0 or indi_genotype == 2.0:
					numerz[indi_index, cell_index] = numerz[indi_index, cell_index] + 1
			else:
				if indi_genotype == 1.0:
					numerz[indi_index, cell_index] = numerz[indi_index, cell_index] + 1
f.close()

numerz2 = np.transpose(numerz)
denerz2 = np.transpose(denerz)
t = open(processed_data_dir + 'ase_contamination_check_counts.txt','w')
t.write('cell_id\t' + '\t'.join(ordered_individuals) + '\n')
for cell_index, cell_id in enumerate(cell_ids):
	t.write(cell_id)
	for individual_index in range(len(ordered_individuals)):
		t.write('\t' + str(numerz2[cell_index, individual_index]) + '/' + str(denerz2[cell_index, individual_index]))
	t.write('\n')
t.close()

frac = numerz2/denerz2.astype(float)
pdb.set_trace()



