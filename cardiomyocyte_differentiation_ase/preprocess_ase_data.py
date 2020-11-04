import numpy as np 
import os
import sys
import pdb


def extract_tissue_names(tissues_file):
	 dicti = {}
	 dicti2 = {}
	 f = open(tissues_file)
	 head_count = 0
	 for line in f:
	 	line = line.rstrip()
	 	data = line.split('\t')
	 	if head_count == 0:
	 		head_count = head_count + 1
	 		continue
	 	dicti[data[2]] = data[0]
	 return dicti

def generate_raw_ase_sample_names(ase_input_data_dir, tissue_names, unfiltered_samples_file):
	# Create dictionary list of sample names
	sample_names = {}
	# loop through individuals (each individual has its own ase file)
	count = 0
	for individual_ase_file in os.listdir(ase_input_data_dir):
		# Error checking
		if individual_ase_file.endswith('.tsv') == False:
			continue
		print(count)
		count = count + 1
		# Open file
		f = open(ase_input_data_dir + individual_ase_file)
		head_count = 0
		site_filtering = {}
		for line in f:
			line = line.rstrip()
			data = line.split()
			# error checking
			if len(data) != 24:
				print('assumption error')
				pdb.set_trace()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			site_id = data[2]
			sample_id = data[5]
			tissue_id = data[7]
			# Skip lines corresponding to samples not in relevent tissues
			if tissue_id not in tissue_names:
				continue
			gene_id = data[20]
			if gene_id == 'NA':
				continue
			low_mapability = data[21]
			mapping_bias_sim = data[22]
			genotype_warning = data[23]
			if low_mapability == '1' or mapping_bias_sim == '1' or genotype_warning == '1':
				continue
			sample_names[sample_id] = tissue_id
		f.close()
	t = open(unfiltered_samples_file, 'w')
	for sample_name in sample_names.keys():
		t.write(sample_name + '\t' + sample_names[sample_name] + '\n')
	t.close()

def extract_ordered_samples(unfiltered_samples_file):
	f = open(unfiltered_samples_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[1])
	return arr

def create_het_site_obj(het_prob_file):
	f = open(het_prob_file)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if line.startswith('##'):
			continue
		if head_count == 0:
			head_count = head_count + 1
			indis = np.asarray(data[9:])
			continue
		rs_id = data[2]
		if rs_id.startswith('rs') == False:
			continue
		het_probs = np.asarray(data[9:]).astype(float)
		if rs_id in dicti:
			print('assumptio eoror')
			pdb.set_trace()
		dicti[rs_id] = {}
		for i, het_prob in enumerate(het_probs):
			indi = indis[i]
			if het_prob > .9:
				dicti[rs_id][indi] = 1
	f.close()
	return dicti

def generate_raw_ase_matrix(ase_input_data_dir, het_prob_file, unfiltered_samples_file, unfiltered_ase_file):
	# Extract ordered array of rna-seq samples
	samples_arr = extract_ordered_samples(unfiltered_samples_file)
	num_samples = len(samples_arr)
	# Create mapping from sample_id to sample_position
	sample_id_to_position = {}
	for i, sample_id in enumerate(samples_arr):
		sample_id_to_position[sample_id] = i
	# Create data structure containing info on whether a site, sample is heterozygous
	het_site_obj = create_het_site_obj(het_prob_file)

	# Initialize mapping from site id to a vector of length num samples where vector contains ase information
	sites = {}
	site_to_gene = {}
	# loop through individuals (each individual has its own ase file)
	count = 0
	for sample_id in samples_arr:
		individual_ase_file = sample_id + '_merged_wasp_corrected.txt'
		indi_id = sample_id.split('_')[0]
		count = count + 1
		# Open file
		f = open(ase_input_data_dir + individual_ase_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			# error checking
			if len(data) != 13:
				print('assumption error')
				pdb.set_trace()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			site_id = data[2]
			if site_id == '.':
				continue
			# Skip lines corresponding to samples not in relevent tissues
			ref_count = data[5]
			total_count = data[7]
			if site_id not in sites:
				sites[site_id] = ['NA']*num_samples
			sample_position = sample_id_to_position[sample_id]
			het_site_dicti = het_site_obj[site_id]
			if int(total_count) == 0:
				continue
			if indi_id in het_site_dicti and len(het_site_dicti) > 2:
				sites[site_id][sample_position] = ref_count + '/' + total_count
		f.close()
	t = open(unfiltered_ase_file,'w')
	t.write('site_id\t' + '\t'.join(np.asarray(samples_arr)) + '\n')
	for site_id in sorted(sites.keys()):
		t.write(site_id + '\t' + '\t'.join(sites[site_id]) + '\n')
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
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
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
		site_id = full_site_id.split(':')[0]
		gene_id = full_site_id.split(':')[1]
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
		site_genes[site_id + ':' + gene_id] = 1
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


def make_min_allele_ase_file(filtered_one_site_per_gene_ase_file, filtered_one_site_per_gene_min_allele_ase_file):
	f = open(filtered_one_site_per_gene_ase_file)
	t = open(filtered_one_site_per_gene_min_allele_ase_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		arr = []
		for ele in data[1:]:
			if ele == 'NA':
				arr.append(ele)
			else:
				ref_count = int(ele.split('/')[0])
				tot_count = int(ele.split('/')[1])
				min_count = min(ref_count, tot_count-ref_count)
				new_string = str(min_count) + '/' + str(tot_count)
				arr.append(new_string)
		t.write(data[0] + '\t' + '\t'.join(arr) + '\n')
	f.close()
	t.close()

def annotate_sample_ids(unfiltered_samples_file, annotated_samples_file, gtex_individual_information_file, cell_type_decomposition_hlv_file, tissue_names):
	f = open(gtex_individual_information_file)
	subject_id_to_covariates = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			subject_header = np.asarray(data[1:])
			continue
		subject_id = data[0]
		if len(data[1:]) != 174:
			anno_arr = data[1:]
			anno_arr.append('NaN')
			anno_arr.append('NaN')
		else:
			anno_arr = data[1:]
		subject_id_to_covariates[subject_id] = np.asarray(anno_arr)
	f.close()
	f = open(cell_type_decomposition_hlv_file)
	subject_id_to_hlv_cell_type_decomp = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		if head_count == 0:
			head_count = head_count + 1
			hlv_cell_type_decomp_header = np.asarray(data[1:5])
			continue
		subject_id = data[0]
		subject_id_to_hlv_cell_type_decomp[subject_id] = np.asarray(data[1:5])
	f.close()
	f = open(unfiltered_samples_file)
	t = open(annotated_samples_file, 'w')
	t.write('sample_id\ttissue_id\t' + '\t'.join(subject_header) + '\t' + '\t'.join(hlv_cell_type_decomp_header) + '\n')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		subject_id = data[0].split('-')[0] + '-' + data[0].split('-')[1]
		t.write(data[0] + '\t' + tissue_names[data[1]] + '\t')
		t.write('\t'.join(subject_id_to_covariates[subject_id]) + '\t')
		if subject_id in subject_id_to_hlv_cell_type_decomp:
			t.write('\t'.join(subject_id_to_hlv_cell_type_decomp[subject_id]) + '\n')
		else:
			t.write('NaN\tNaN\tNaN\tNaN\n')
	f.close()
	t.close()
	f = open(annotated_samples_file)
	print(annotated_samples_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 180:
			print(data)
			print('assumption errror')
	f.close()

def make_covariate_file_based_on_sample_repeat_structure(annotated_samples_file, covariate_file):
	f = open(annotated_samples_file)
	indis = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi = data[1].split('_')[0]
		indis[indi] = 1
	f.close()
	indi_arr = []
	for indi in indis.keys():
		indi_arr.append(indi)
	num_indis = len(indi_arr)
	used_indi_arr = indi_arr[:(num_indis-1)]
	num_used_indis = len(used_indi_arr)
	indi_to_pos = {}
	for i, indi in enumerate(used_indi_arr):
		indi_to_pos[indi] = i
	f = open(annotated_samples_file)
	t = open(covariate_file, 'w')
	head_count = 0 
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		covs = np.zeros(num_used_indis)
		indi = data[1].split('_')[0]
		if indi in indi_to_pos:
			pos = indi_to_pos[indi]
			covs[pos] = 1
		t.write('\t'.join(covs.astype(str)) + '\n')
	t.close()
	f.close()

def make_sample_overlap_file(annotated_samples_file, covariate_file):
	f = open(annotated_samples_file)
	indis = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		indi = data[1].split('_')[0]
		indis[indi] = 1
	f.close()
	indi_arr = []
	for indi in indis.keys():
		indi_arr.append(indi)
	num_indis = len(indi_arr)
	used_indi_arr = np.copy(indi_arr)
	num_used_indis = len(used_indi_arr)
	indi_to_pos = {}
	for i, indi in enumerate(used_indi_arr):
		indi_to_pos[indi] = i
	f = open(annotated_samples_file)
	t = open(covariate_file, 'w')
	head_count = 0 
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		covs = np.zeros(num_used_indis)
		indi = data[1].split('_')[0]
		pos = indi_to_pos[indi]
		t.write(str(pos) + '\n')
	t.close()
	f.close()

def add_pcs_to_annotated_samples_file(annotated_samples_file, pca_file, updated_annotated_samples_file):
	f = open(pca_file)
	sample_names = []
	pcs = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_names.append(data[0])
		pcs.append(data[1:3])
	f.close()
	f = open(annotated_samples_file)
	t = open(updated_annotated_samples_file, 'w')
	head_count = 0
	count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'expression_pc1' + '\t' + 'expression_pc2' + '\n')
			continue
		if data[1] != sample_names[count]:
			print('assumption error')
			pdb.set_trace()
		t.write(line + '\t' + '\t'.join(pcs[count]) + '\n')
		count = count + 1
	f.close()
	t.close()


ase_input_data_dir = sys.argv[1]
het_prob_file = sys.argv[2]
annotated_samples_file = sys.argv[3]
pca_file = sys.argv[4]
output_root = sys.argv[5]

# Add PCs to annotated samples file
updated_annotated_samples_file = output_root + 'annotated_samples.txt'
#add_pcs_to_annotated_samples_file(annotated_samples_file, pca_file, updated_annotated_samples_file)


# Generate unfiltered ase matrix 
# columns are rna-seq samples
# rows are het sites
unfiltered_ase_file = output_root + 'raw_ase_counts.txt'
#generate_raw_ase_matrix(ase_input_data_dir, het_prob_file, annotated_samples_file, unfiltered_ase_file)


# Generate Filtered ase matrix 
# columns are rna-seq samples
# rows are het sites

# Require site to have at least min_samples het individuals
fraction_samples = .35
fraction_monoallelic = .5

filtered_ase_file = output_root + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '.txt'
#filter_ase_matrix(unfiltered_ase_file, filtered_ase_file, fraction_samples, fraction_monoallelic)



filtered_min_allele_ase_file = output_root + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_min_allele.txt'
#make_min_allele_ase_file(filtered_ase_file, filtered_min_allele_ase_file)



covariate_file = output_root + 'covariates.txt'
# make_covariate_file_based_on_sample_repeat_structure(annotated_samples_file, covariate_file)




sample_overlap_file = output_root + 'sample_overlap.txt'
make_sample_overlap_file(annotated_samples_file, sample_overlap_file)








#filtered_one_site_per_gene_ase_file = output_root + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_one_site_per_gene.txt'
#filter_ase_matrix_to_one_site_per_gene(filtered_ase_file, filtered_one_site_per_gene_ase_file)

