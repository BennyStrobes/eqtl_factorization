import numpy as np 
import os
import sys
import pdb


def get_ordered_cell_names_from_ase_file(ase_file):
	f = open(ase_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			cells = np.asarray(data[1:])
			continue
		break
	return cells

def get_ordered_pc1(cell_covariates_file, ase_cells):
	f = open(cell_covariates_file)
	head_count = 0
	counter = 0
	pcs = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if ase_cells[counter] != data[0]:
			print('assumption error')
			pdb.set_trace()
		counter = counter + 1
		pcs.append(data[89])
	return np.asarray(pcs).astype(float)

def equal_bin(N, m):
	sep = (N.size/float(m))*np.arange(1,m+1)
	idx = sep.searchsorted(np.arange(N.size))
	return idx[N.argsort().argsort()]

def split_ase_file_into_bins(ase_file, num_bins, cell_bins, allelic_counts_output_file):
	t = open(allelic_counts_output_file, 'w')
	f = open(ase_file)
	header = []
	for bin_num in range(num_bins):
		header.append('bin_' + str(bin_num))
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(data[0] + '\t' + '\t'.join(header) + '\n')
			continue
		exonic_site = data[0]
		counts = np.asarray(data[1:])
		allelic_counts = []
		for bin_num in range(num_bins):
			numerz = 0
			denerz = 0
			bin_indices = np.where(cell_bins==bin_num)[0]
			for index in bin_indices:
				if counts[index] != 'NA':
					numer = int(counts[index].split('/')[0])
					dener = int(counts[index].split('/')[1])
					corrected_numer = np.min((numer, dener-numer))
					numerz = numerz + corrected_numer 
					denerz = denerz + dener
			if denerz == 0:
				allelic_counts.append('NA')
				print('miss')
			else:
				allelic_counts.append(str(numerz) + '/' + str(denerz))
		t.write(exonic_site + '\t' + '\t'.join(allelic_counts) + '\n')
	f.close()
	t.close()

def compute_allelic_fractions(allelic_counts_output_file, allelic_fraction_output_file):
	f = open(allelic_counts_output_file)
	t = open(allelic_fraction_output_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		allelic_fractions = []
		exonic_site = data[0]
		allelic_counts = data[1:]
		for allelic_count in allelic_counts:
			if allelic_count == 'NA':
				allelic_fractions.append('NA')
			else:
				numer = float(allelic_count.split('/')[0])
				dener = float(allelic_count.split('/')[1])
				allelic_fractions.append(numer/dener)
		t.write(exonic_site + '\t' + '\t'.join(np.asarray(allelic_fractions).astype(str)) + '\n')
	f.close()
	t.close()

ase_file = sys.argv[1]
cell_covariates_file = sys.argv[2]
num_bins = int(sys.argv[3])
loading_file = sys.argv[4]
output_root = sys.argv[5]
distribution = sys.argv[6]


U = np.loadtxt(loading_file)
K = U.shape[1]
ase_cells = get_ordered_cell_names_from_ase_file(ase_file)

#pc1 = get_ordered_pc1(cell_covariates_file, ase_cells)
for k in range(K):
	U_k = U[:,k]

	cell_bins = equal_bin(U_k, num_bins)

	allelic_counts_output_file = output_root + str(k) + '_' + 'allelic_counts.txt'
	split_ase_file_into_bins(ase_file, num_bins, cell_bins, allelic_counts_output_file)

	allelic_fraction_output_file = output_root + str(k) + '_' + 'allelic_fraction.txt'
	compute_allelic_fractions(allelic_counts_output_file, allelic_fraction_output_file)

