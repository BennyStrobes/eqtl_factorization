import numpy as np 
import os
import sys
import pdb


def summarize_ase_counts_file(counts_file, summary_file):
	f = open(counts_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			num_cells = len(data[1:])
			cell_names = np.asarray(data[1:])
			ge_3 = np.zeros(num_cells)
			ge_3_biallelic = np.zeros(num_cells)
			ge_5 = np.zeros(num_cells)
			ge_10 = np.zeros(num_cells)
			ge_15 = np.zeros(num_cells)
			continue
		exonic_site = data[0]
		counts = data[1:]
		counter = counter + 1
		for i, count in enumerate(counts):
			if count == 'NA':
				continue
			alt_count = int(count.split('/')[0])
			total_count = int(count.split('/')[1])
			if total_count >= 3:
				ge_3[i] = ge_3[i] + 1
				if alt_count != total_count and (total_count - alt_count) != total_count:
					ge_3_biallelic[i] = ge_3_biallelic[i] + 1
			if total_count >= 5:
				ge_5[i] = ge_5[i] + 1
			if total_count >= 10:
				ge_10[i] = ge_10[i] + 1
			if total_count >= 15:
				ge_15[i] = ge_15[i] + 1
	f.close()
	t = open(summary_file, 'w')
	# header
	t.write('cell_id\tge_3\tge_3_biallelic\tge_5\tge_10\tge_15\n')
	for i in range(num_cells):
		t.write(cell_names[i] + '\t')
		t.write(str(ge_3[i]) + '\t')
		t.write(str(ge_3_biallelic[i]) + '\t')
		t.write(str(ge_5[i]) + '\t')
		t.write(str(ge_10[i]) + '\t')
		t.write(str(ge_15[i]) + '\n')
	t.close()

processed_data_dir = sys.argv[1]



raw_counts_file = processed_data_dir + 'raw_ase_counts.txt'
raw_counts_summary_file = processed_data_dir + 'raw_ase_counts_summary.txt'
summarize_ase_counts_file(raw_counts_file, raw_counts_summary_file)