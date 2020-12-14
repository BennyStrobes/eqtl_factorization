import numpy as np 
import os
import sys
import pdb
import gzip

def extract_ordered_cells(annotated_samples_file):
	f = open(annotated_samples_file)
	head_count = 0
	ordered_cells = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[0]
		ordered_cells.append(cell_id)
	f.close()
	ordered_cells = np.asarray(ordered_cells)
	cell_mapping = {}
	for index, cell_id in enumerate(ordered_cells):
		cell_mapping[cell_id] = index
	return ordered_cells, cell_mapping

def fill_in_dictionary_with_ase_counts(ase_dicti, alt_count_file, total_count_file, cell_mapping, num_cells):
	f = open(alt_count_file)
	g = open(total_count_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data_alt = line.split('\t')
		data_total = g.readline().rstrip().split('\t')
		if head_count == 0:
			head_count = head_count + 1
			cell_names = np.asarray(data_alt[1:])
			if np.array_equal(cell_names, np.asarray(data_total[1:])) == False:
				print('assumption error')
				pdb.set_trace()
			'''
			indices = []
			valid = []
			for cell_name in cell_names:
				if cell_name in cell_mapping:
					indices.append(cell_mapping[cell_name])
					valid.append(True)
				else:
					valid.append(False)
			indices = np.asarray(indices)
			valid = np.asarray(valid)
			'''
			continue
		if data_alt[0] != data_alt[0]:
			print('assumptino eroror')
			pdb.set_trace()
		if len(data_alt) != len(data_total):
			print('assumpotinoi eroror')
			pdb.set_trace()
		site_id = data_alt[0]
		if site_id not in ase_dicti:
			ase_vec = ['NA']*num_cells
		else:
			ase_vec = ase_dicti[site_id]
		alt_counts = data_alt[1:]
		total_counts = data_total[1:]

		for i, alt_count in enumerate(alt_counts):
			total_count = total_counts[i]
			if total_count == '':
				continue
			cell_name = cell_names[i]
			if cell_name not in cell_mapping:
				continue
			master_index = cell_mapping[cell_name]
			ase_vec[master_index] = alt_count + '/' + total_count
		ase_dicti[site_id] = ase_vec
	return ase_dicti


def generate_raw_ase_counts_file(ordered_cells, cell_mapping, ase_input_data_dir, raw_ase_counts_file):
	# First initialize dictionary to keep track of data
	# Key is exonic site. Value is vector of length num cells orderedd by ordered cells
	ase_dicti = {}
	# Loop through ase files (one per individual)
	counter = 0
	for file_name in os.listdir(ase_input_data_dir):
		if file_name.endswith('.altcount.tsv') == False:
			continue
		counter = counter + 1
		print(file_name + '\t' + str(counter))
		# Each individual has an 'alt count file' and a 'total count file'
		info = file_name.split('.')
		alt_count_file = ase_input_data_dir + file_name
		total_count_file = ase_input_data_dir + info[0] + '.' + info[1] + '.' + info[2] + '.totalcount.tsv'
		ase_dicti = fill_in_dictionary_with_ase_counts(ase_dicti, alt_count_file, total_count_file, cell_mapping, len(ordered_cells))
	# Print to output file
	t = open(raw_ase_counts_file,'w')
	# header 
	t.write('exonic_site\t' + '\t'.join(ordered_cells) + '\n')
	for het_site_id in ase_dicti.keys():
		allelic_counts = ase_dicti[het_site_id]
		t.write(het_site_id + '\t' + '\t'.join(allelic_counts) + '\n')
	t.close()

def compute_fraction_monoallelic(het_counts):
	total_het = len(het_counts)
	total_monoallelic = 0
	for het_count in het_counts:
		ref_count = int(het_count.split('/')[0])
		total_count = int(het_count.split('/')[1])
		if total_count < 3 or min(ref_count/float(total_count), (total_count-ref_count)/float(total_count)) < .01:
			total_monoallelic = total_monoallelic + 1
	return float(total_monoallelic)/total_het

def filter_ase_matrix(unfiltered_ase_file, filtered_ase_file, fraction_samples, fraction_monoallelic_thresh):
	f = open(unfiltered_ase_file)
	t = open(filtered_ase_file, 'w')
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		counter = counter + 1
		if np.mod(counter,10000) == 0:
			print(counter)
		site_id = data[0]

		counts = np.asarray(data[1:])
		num_samples = len(counts)
		het_indis = np.where(counts != 'NA')[0]

		if len(het_indis)/float(num_samples) < fraction_samples:
			continue
		fraction_monoallelic = compute_fraction_monoallelic(counts[het_indis])
		if fraction_monoallelic > fraction_monoallelic_thresh:
			continue
		t.write(line + '\n')
	f.close()
	t.close()

def get_gene_name(stringer):
    gene_name = 'nuller'
    fields = stringer.split(';')
    for field in fields:
        field_info = field.split(' ')
        if field_info[0] == 'gene_id':
            gene_name = field_info[1].split('"')[1]
    if gene_name == 'nuller':
        print('get_gene_name assumption erroro')
    if gene_name.startswith('ENSG') == False:
        print('get_gene_name assumption erroro')
    return gene_name

def get_gene_type(stringer):
    gene_name = 'nuller'
    fields = stringer.split('; ')
    for field in fields:
        field_info = field.split(' ')
        if field_info[0] == 'gene_type':
            gene_name = field_info[1].split('"')[1]
    if gene_name == 'nuller':
        print('get_gene_type assumption erroro')
    return gene_name

def get_gene_status(stringer):
    gene_name = 'nuller'
    fields = stringer.split('; ')
    for field in fields:
        field_info = field.split(' ')
        if field_info[0] == 'gene_status':
            gene_name = field_info[1].split('"')[1]
    if gene_name == 'nuller':
        print('get_gene_type assumption erroro')
    return gene_name


#Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
def add_gene_to_chromosome_object(chromosome, start, end, gene_name):
    for pos in range(start, end+1):
        if chromosome[pos] == 'NULL':  # No genes already mapped to this position
            chromosome[pos] = gene_name
        else:  # at least one gene already mapped to this position
            chromosome[pos] = chromosome[pos] + ',' + gene_name
    return chromosome

def make_chromosome_with_exon_positions(gencode_gene_annotation_file, chrom_num):
	f = gzip.open(gencode_gene_annotation_file)
	chromosome = ['NULL']*259250621
	types = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip headers
		if line.startswith('##'):
			continue
		if len(data) != 9:
			print('assumption erorro!')
			pdb.set_trace()
		# Extract relevent fields
		line_chrom_num = data[0]
		gene_part = data[2]
		gene_name = get_gene_name(data[8])
		gene_type = get_gene_type(data[8])
		gene_status = get_gene_status(data[8])
		start = int(data[3])  # 5' gene start site (this is the min of all UTR,exons,etc)
		end = int(data[4])  # 3' gene end site (this is the max of all UTR,exons,etc)
		if line_chrom_num != 'chr' + chrom_num:  # limit to chromosome of interest
			continue
		if gene_part != 'exon':  # ignore other parts of the gene as 'gene' encomposes everything
			continue
		if gene_type != 'protein_coding':
			continue
		if gene_status != 'KNOWN':
			continue
		chromosome = add_gene_to_chromosome_object(chromosome, start, end, gene_name)
	f.close()
	return chromosome

def map_to_genes(ase_site_file, gene_mapped_ase_site_file, gencode_gene_annotation_file):
	t = open(gene_mapped_ase_site_file, 'w')
	for chrom_num in range(1, 23):
		print(chrom_num)
		chromosome = make_chromosome_with_exon_positions(gencode_gene_annotation_file, str(chrom_num))
		f = open(ase_site_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				if chrom_num == 1:
					t.write(line + '\n')
				continue
			site_id = data[0]
			site_info = site_id.split('_')
			site_chrom = site_info[0]
			if site_chrom != str(chrom_num):
				continue
			site_pos = int(site_info[1])
			associated_genes = chromosome[site_pos]
			if associated_genes == 'NULL':
				t.write(site_id + '_' + 'NULL' + '\t' + '\t'.join(data[1:]) + '\n')
			else:
				unique_genes = np.unique(associated_genes.split(','))
				if len(unique_genes) == 1:
					t.write(site_id + '_' + unique_genes[0] + '\t' + '\t'.join(data[1:]) + '\n')
				else:
					t.write(site_id + '_' + 'NULL' + '\t' + '\t'.join(data[1:]) + '\n')
		f.close()
	t.close()
def compute_depth(counts):
	depth = 0
	for count in counts:
		depth = depth + int(count.split('/')[1])
	return depth

def filter_ase_matrix_to_one_site_per_gene(filtered_ase_file, filtered_one_site_per_gene_ase_file):
	# Create mapping from gene id to best site and total number of reads at that site
	f = open(filtered_ase_file)
	head_count = 0
	gene_to_site = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		full_site_id = data[0]
		site_id = full_site_id.split('_')[0] + '_' + full_site_id.split('_')[1] + '_' + full_site_id.split('_')[2] + '_' + full_site_id.split('_')[3]
		gene_id = full_site_id.split('_')[4]
		counts = np.asarray(data[1:])
		het_indis = np.where(counts != 'NA')[0]
		depth = compute_depth(counts[het_indis])
		if gene_id not in gene_to_site:
			gene_to_site[gene_id] = (site_id, depth)
		else:
			old_tuple = gene_to_site[gene_id]
			old_site_id = old_tuple[0]
			old_depth = old_tuple[1]
			if depth > old_depth:
				gene_to_site[gene_id] = (site_id, depth)
	f.close()
	# Create dictionary list of site_genes
	site_genes = {}
	for gene_id in gene_to_site.keys():
		site_id = gene_to_site[gene_id][0]
		site_genes[site_id + '_' + gene_id] = 1
		head_count = 0
	f = open(filtered_ase_file)
	t = open(filtered_one_site_per_gene_ase_file, 'w')
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		if data[0] not in site_genes:
			continue
		t.write(line + '\n')
	f.close()
	t.close()

def remove_zero_cells(filtered_ase_gene_mapped_one_per_gene_file, final_ase_file):
	head_count = 0
	f = open(filtered_ase_gene_mapped_one_per_gene_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			cell_ids = data[1:]
			totals = np.zeros(len(cell_ids))
			continue
		counts = np.asarray(data[1:])
		for i, count in enumerate(counts):
			if count != 'NA':
				num_reads = float(count.split('/')[1])
				if num_reads >= 5:
					totals[i] = totals[i] + 1
	f.close()
	f = open(filtered_ase_gene_mapped_one_per_gene_file)
	t = open(final_ase_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		t.write(data[0])
		for i, ele in enumerate(data[1:]):
			if totals[i] > 0.5:
				t.write('\t' + ele)
		t.write('\n')
	f.close()
	t.close()

def filter_cell_info_file(annotated_samples_file, filtered_cell_info_file, final_ase_file):
	f = open(final_ase_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			num_cells = len(data[1:])
			num_expressed = np.zeros(num_cells)
			num_biallelic = np.zeros(num_cells)
			numerator = np.zeros(num_cells)
			denomenator = np.zeros(num_cells)
			continue
		for i, stringer in enumerate(data[1:]):
			if stringer == 'NA':
				continue
			numer = int(stringer.split('/')[0])
			dener = int(stringer.split('/')[1])
			numer = np.min((numer, dener-numer))
			numerator[i] = numerator[i] + numer
			denomenator[i] = denomenator[i] + dener
			if dener < 10:
				continue
			num_expressed[i] = num_expressed[i] + 1
			true_numer = np.min((numer, dener-numer))
			if true_numer > 0:
				num_biallelic[i] = num_biallelic[i] + 1
	f.close()
	frac = num_biallelic.astype(float)/num_expressed
	frac2 = numerator.astype(float)/denomenator
	valid_cells = {}
	f = open(final_ase_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			if len(data[1:]) != len(frac):
				print('assumption error')
			for i, cell_name in enumerate(data[1:]):
				valid_cells[cell_name] = (num_biallelic[i], num_expressed[i], frac[i], frac2[i])
			continue
		break
	f.close()

	f = open(annotated_samples_file)
	t = open(filtered_cell_info_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'number_biallelic_sites\tnum_expressed_sites' '\t' + 'fraction_biallelic' + '\t' + 'global_allelic_fraction' + '\n')
			continue
		cell_name = data[0]
		if cell_name in valid_cells:
			t.write(line + '\t' + str(valid_cells[cell_name][0]) + '\t' + str(valid_cells[cell_name][1]) + '\t' + str(valid_cells[cell_name][2])+ '\t' + str(valid_cells[cell_name][3]) + '\n')
	f.close()
	t.close()

def generate_sample_overlap_file(filtered_cell_info_file, cell_overlap_file, category_name):
	f = open(filtered_cell_info_file)
	head_count = 0
	donors = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			category_index = np.where(np.asarray(data) == category_name)[0][0]
			continue
		donor_id = data[category_index]
		donors.append(donor_id)
	f.close()
	unique_donors = np.unique(donors)
	donor_to_index = {}
	for i, donor in enumerate(unique_donors):
		donor_to_index[donor] = i
	f = open(filtered_cell_info_file)
	t = open(cell_overlap_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		donor_id = data[category_index]
		index = donor_to_index[donor_id]
		t.write(str(index) + '\n')
	t.close()
	f.close()

def subsample_cells(final_ase_file, subsampled_ase_file, subsampling_fraction):
	f = open(final_ase_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			num_cells = len(data[1:])
			continue
		break
	f.close()
	num_subsampled_cells = np.floor(num_cells*subsampling_fraction).astype(int)
	indices = np.asarray(sorted(np.random.choice(range(num_cells), size=num_subsampled_cells, replace=False)))
	f = open(final_ase_file)
	t = open(subsampled_ase_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		t.write(data[0] + '\t' + '\t'.join(np.asarray(data[1:])[indices]) + '\n')
	f.close()
	t.close()

def remove_low_biallelic_fraction_cells(filtered_cell_info_file, final_ase_file, ase_high_biallelic_fraction_file):
	valid_cells = {}
	all_cells = {}
	f = open(filtered_cell_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[0]
		num_biallelic_sites = float(data[97])
		fraction_biallelic_sites = float(data[99])
		if num_biallelic_sites >= 50.0 and fraction_biallelic_sites > .85:
			valid_cells[cell_id] = 0
		all_cells[cell_id] = 0
	f.close()
	f = open(final_ase_file)
	t = open(ase_high_biallelic_fraction_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = data[0]
			cell_ids = np.asarray(data[1:])
			valid_indices = []
			for i, cell_id in enumerate(cell_ids):
				if cell_id in valid_cells:
					valid_indices.append(i)
			valid_indices = np.asarray(valid_indices)
			t.write(header + '\t' + '\t'.join(cell_ids[valid_indices]) + '\n')
			for cell_id in cell_ids[valid_indices]:
				if cell_id not in valid_cells:
					print('assumption error')
					pdb.set_trace()
			continue
		header = data[0]
		counts = np.asarray(data[1:])
		t.write(header + '\t' + '\t'.join(counts[valid_indices]) + '\n')
	f.close()
	t.close()

def filter_go_term_loading_file(go_term_loading_file, filtered_go_term_loading_file, final_ase_file):
	valid_cells = {}
	f = open(final_ase_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			ordered_cells = np.asarray(data[1:])
			for i, cell_name in enumerate(data[1:]):
				valid_cells[cell_name] = 0
			continue
		break
	f.close()
	go_data = np.transpose(np.loadtxt(go_term_loading_file, dtype=str,delimiter='\t', comments='*'))
	go_cells = go_data[1:,0]
	go_terms = go_data[0,1:]
	go_loadings = go_data[1:,1:]
	t = open(filtered_go_term_loading_file, 'w')
	t.write('cell_id\t' + '\t'.join(go_terms) + '\n')
	for i, cell_name in enumerate(go_cells):
		if cell_name in valid_cells:
			loadings = go_loadings[i,:]
			t.write(cell_name + '\t' + '\t'.join(loadings) + '\n')
	t.close()

def add_cell_cycle_info_to_annotated_samples_file(annotated_samples_file, cell_cycle_file, cell_state_file, cc_annotated_samples_file):
	f = open(cell_cycle_file)
	head_count = 0
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[3]
		mapping[cell_id] = np.asarray(data[:3])
	f.close()
	f = open(cell_state_file)
	mapping2 = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[0]
		mapping2[cell_id] = np.asarray(data[1:])
	f.close()
	f = open(annotated_samples_file)
	t = open(cc_annotated_samples_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'S_score\tG2M_score\tcc_phase\trespiration\tdifferentiation\tG1_S_transition\tsterol_biosynthesis\tG2_M_transition\tpseudotime2\n')
			continue
		cell_id = data[0]
		info = mapping[cell_id]
		info2 = mapping2[cell_id]
		t.write(line + '\t' + '\t'.join(info) + '\t' + '\t'.join(info2) + '\n')
	f.close()
	t.close()


# Directory containing input allelic counts for each individual
ase_input_data_dir = sys.argv[1]
# File containing list of cells (and annotations describing those cells) used by cuomo et al. analysis
annotated_samples_file = sys.argv[2]
# Output directory
processed_data_dir = sys.argv[3]
# Genecode gene annotation file
gencode_gene_annotation_file = sys.argv[4]
# cell cycle file
cell_cycle_file = sys.argv[5]
# Cell State file
cell_state_file = sys.argv[6]


go_term_loading_file = processed_data_dir + 'go_terms_cell_loadings.txt'

# First create ordered list of cells (ordered_cells) and mapping from cell to cell_index (cell_mapping)
ordered_cells, cell_mapping = extract_ordered_cells(annotated_samples_file)
num_cells = len(ordered_cells)

# update annotated samples file with cell cycle info
cc_annotated_samples_file = processed_data_dir + 'annotated_samples_with_cell_cycle.txt'
add_cell_cycle_info_to_annotated_samples_file(annotated_samples_file, cell_cycle_file, cell_state_file, cc_annotated_samples_file)


# Rows are exonic sites, columns are cells (according to ordered_cells)
# 'NA' if cell is not heterozygous at site or cell doesn't have expression at site
# otherwise ref_count/total_count
raw_ase_counts_file = processed_data_dir + 'raw_ase_counts.txt'
#generate_raw_ase_counts_file(ordered_cells, cell_mapping, ase_input_data_dir, raw_ase_counts_file)


# Require site to have at least min_samples het individuals
fraction_samples = .3
fraction_monoallelic = .5

filtered_ase_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '.txt'
#filter_ase_matrix(raw_ase_counts_file, filtered_ase_file, fraction_samples, fraction_monoallelic)

# map to gene
filtered_ase_gene_mapped_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_gene_mapped.txt'
#map_to_genes(filtered_ase_file, filtered_ase_gene_mapped_file, gencode_gene_annotation_file)

filtered_ase_gene_mapped_one_per_gene_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_gene_mapped_one_per_gene.txt'
#filter_ase_matrix_to_one_site_per_gene(filtered_ase_gene_mapped_file, filtered_ase_gene_mapped_one_per_gene_file)


final_ase_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_final.txt'
#remove_zero_cells(filtered_ase_gene_mapped_one_per_gene_file, final_ase_file)


filtered_cell_info_file = processed_data_dir + 'cell_info_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '.txt'
filter_cell_info_file(cc_annotated_samples_file, filtered_cell_info_file, final_ase_file)

filtered_go_term_loading_file = processed_data_dir + 'go_terms_cell_loadings_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '.txt'
#filter_go_term_loading_file(go_term_loading_file, filtered_go_term_loading_file, final_ase_file)

cell_overlap_file = processed_data_dir + 'cell_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '.txt'
#generate_sample_overlap_file(filtered_cell_info_file, cell_overlap_file, 'donor_long_id')

plate_overlap_file = processed_data_dir + 'cell_plate_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '.txt'
#generate_sample_overlap_file(filtered_cell_info_file, plate_overlap_file, 'plate_id')

experiment_overlap_file = processed_data_dir + 'cell_experiment_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '.txt'
#generate_sample_overlap_file(filtered_cell_info_file, experiment_overlap_file, 'experiment')


############################################
# Subsample
############################################
final_ase_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_final.txt'
subsampled_ase_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_final_subsampled.txt'
subsampling_fraction = .1
#subsample_cells(final_ase_file, subsampled_ase_file, subsampling_fraction)

filtered_cell_info_subsampled_file = processed_data_dir + 'cell_info_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_subsampled.txt'
filter_cell_info_file(cc_annotated_samples_file, filtered_cell_info_subsampled_file, subsampled_ase_file)

cell_overlap_subsampled_file = processed_data_dir + 'cell_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_subsampled.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, cell_overlap_subsampled_file, 'donor_long_id')

filtered_go_term_loading_file = processed_data_dir + 'go_terms_cell_loadings_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_subsampled.txt'
#filter_go_term_loading_file(go_term_loading_file, filtered_go_term_loading_file, subsampled_ase_file)


plate_overlap_subsampled_file = processed_data_dir + 'cell_plate_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_subsampled.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, plate_overlap_subsampled_file, 'plate_id')

experiment_overlap_subsampled_file = processed_data_dir + 'cell_experiment_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_subsampled.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, experiment_overlap_subsampled_file, 'experiment')

############################################
# Remove cells with low biallelic fraction
############################################
ase_high_biallelic_fraction_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_final_high_biallelic_fraction_only.txt'
#remove_low_biallelic_fraction_cells(filtered_cell_info_file, final_ase_file, ase_high_biallelic_fraction_file)

filtered_cell_info_high_biallelic_fraction_file = processed_data_dir + 'cell_info_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only.txt'
filter_cell_info_file(cc_annotated_samples_file, filtered_cell_info_high_biallelic_fraction_file, ase_high_biallelic_fraction_file)

filtered_go_term_loading_file = processed_data_dir + 'go_terms_cell_loadings_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only.txt'
#filter_go_term_loading_file(go_term_loading_file, filtered_go_term_loading_file, ase_high_biallelic_fraction_file)

cell_overlap_high_biallelic_fraction_file = processed_data_dir + 'cell_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only.txt'
#generate_sample_overlap_file(filtered_cell_info_high_biallelic_fraction_file, cell_overlap_high_biallelic_fraction_file, 'donor_long_id')

plate_overlap_high_biallelic_fraction_file = processed_data_dir + 'cell_plate_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only.txt'
#generate_sample_overlap_file(filtered_cell_info_high_biallelic_fraction_file, plate_overlap_high_biallelic_fraction_file, 'plate_id')

experiment_overlap_high_biallelic_fraction_file = processed_data_dir + 'cell_experiment_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only.txt'
#generate_sample_overlap_file(filtered_cell_info_high_biallelic_fraction_file, experiment_overlap_high_biallelic_fraction_file, 'experiment')



############################################
# Subsample high BF cells
############################################
subsampled_ase_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_final_high_biallelic_fraction_only_subsampled.txt'
subsampling_fraction = .15
#subsample_cells(ase_high_biallelic_fraction_file, subsampled_ase_file, subsampling_fraction)

filtered_cell_info_subsampled_file = processed_data_dir + 'cell_info_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled.txt'
filter_cell_info_file(cc_annotated_samples_file, filtered_cell_info_subsampled_file, subsampled_ase_file)

filtered_go_term_loading_subsampled_file = processed_data_dir + 'go_terms_cell_loadings_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled.txt'
#filter_go_term_loading_file(go_term_loading_file, filtered_go_term_loading_subsampled_file, subsampled_ase_file)

cell_overlap_subsampled_file = processed_data_dir + 'cell_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, cell_overlap_subsampled_file, 'donor_long_id')

plate_overlap_subsampled_file = processed_data_dir + 'cell_plate_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, plate_overlap_subsampled_file, 'plate_id')

experiment_overlap_subsampled_file = processed_data_dir + 'cell_experiment_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, experiment_overlap_subsampled_file, 'experiment')





############################################
# Subsample high BF cells
############################################
subsampled_ase_file = processed_data_dir + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_final_high_biallelic_fraction_only_subsampled_small.txt'
subsampling_fraction = .05
#subsample_cells(ase_high_biallelic_fraction_file, subsampled_ase_file, subsampling_fraction)

filtered_cell_info_subsampled_file = processed_data_dir + 'cell_info_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled_small.txt'
filter_cell_info_file(cc_annotated_samples_file, filtered_cell_info_subsampled_file, subsampled_ase_file)

filtered_go_term_loading_subsampled_file = processed_data_dir + 'go_terms_cell_loadings_after_filtering_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled_small.txt'
#filter_go_term_loading_file(go_term_loading_file, filtered_go_term_loading_subsampled_file, subsampled_ase_file)

cell_overlap_subsampled_file = processed_data_dir + 'cell_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled_small.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, cell_overlap_subsampled_file, 'donor_long_id')

plate_overlap_subsampled_file = processed_data_dir + 'cell_plate_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled_small.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, plate_overlap_subsampled_file, 'plate_id')

experiment_overlap_subsampled_file = processed_data_dir + 'cell_experiment_overlap_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_high_biallelic_fraction_only_subsampled_small.txt'
#generate_sample_overlap_file(filtered_cell_info_subsampled_file, experiment_overlap_subsampled_file, 'experiment')
