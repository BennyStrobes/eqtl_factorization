import numpy as np 
import os
import sys
import pdb
import pystan
import pickle
import scipy.special


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
			numerz = []
			denerz = []
			bin_indices = np.where(cell_bins==bin_num)[0]
			for index in bin_indices:
				if counts[index] != 'NA':
					numer = int(counts[index].split('/')[0])
					dener = int(counts[index].split('/')[1])
					if dener == 1:
						continue
					corrected_numer = np.min((numer, dener-numer))
					numerz.append(corrected_numer)
					denerz.append(dener)
					#numerz = numerz + corrected_numer 
					#denerz = denerz + dener
			if len(denerz) == 0:
				allelic_counts.append('NA')
				print('miss')
			else:
				ratios = []
				for i, numer in enumerate(numerz):
					dener = denerz[i]
					if dener == 0:
						pdb.set_trace()
					ratios.append(str(numer) + '/' + str(dener))
				ratios = np.asarray(ratios)
				allelic_counts.append(','.join(ratios))
		t.write(exonic_site + '\t' + '\t'.join(allelic_counts) + '\n')
	f.close()
	t.close()

def compute_allelic_mean(allelic_counts_output_file, allelic_fraction_output_file):
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
				ratios = allelic_count.split(',')
				numerz = []
				denerz = []
				for ratio in ratios:
					numer = ratio.split('/')[0]
					dener = ratio.split('/')[1]
					numerz.append(int(numer))
					denerz.append(int(dener))
				x = np.ones((len(numerz), 1))
				data = dict(N=len(numerz), P=x.shape[1], x=x, ys=np.asarray(numerz), ns=np.asarray(denerz),concShape=1.0001, concRate=1e-4)
				op = MODEL.optimizing(data = data, verbose=False,iter=5000,seed=1)
				mean = scipy.special.expit(op['beta'])
				#numer = float(allelic_count.split('/')[0])
				#dener = float(allelic_count.split('/')[1])
				allelic_fractions.append(mean)
		t.write(exonic_site + '\t' + '\t'.join(np.asarray(allelic_fractions).astype(str)) + '\n')
	f.close()
	t.close()

ase_file = sys.argv[1]
cell_covariates_file = sys.argv[2]
num_bins = int(sys.argv[3])
loading_file = sys.argv[4]
output_root = sys.argv[5]
distribution = sys.argv[6]

if distribution == 'binomial':
	MODEL = pickle.load(open('/work-zfs/abattle4/bstrober/temp/binomial_glm.pkl', 'rb'))
	MODEL_ALT = pickle.load(open('/work-zfs/abattle4/bstrober/temp/folded_binomial_glm.pkl', 'rb'))
elif distribution == 'folded_binomial':
	MODEL = pickle.load(open('/work-zfs/abattle4/bstrober/temp/folded_binomial_glm.pkl', 'rb'))
elif distribution == 'folded_beta_binomial':
	MODEL = pickle.load(open('/work-zfs/abattle4/bstrober/temp/folded_beta_binomial_glm.pkl', 'rb'))
elif distribution == 'beta_binomial':
	MODEL = pickle.load(open('/work-zfs/abattle4/bstrober/temp/betabinomial_glm.pkl', 'rb'))

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
	compute_allelic_mean(allelic_counts_output_file, allelic_fraction_output_file)

