import numpy as np 
import os
import sys
import pdb
import gzip


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
	for line in f:
		line = line.rstrip()
		data = line.split()
		arr.append(data[0])
	return arr


def compute_fraction_monoallelic(het_counts):
	total_het = len(het_counts)
	total_monoallelic = 0
	for het_count in het_counts:
		ref_count = int(het_count.split('/')[0])
		total_count = int(het_count.split('/')[1])
		if total_count < 3 or min(ref_count/float(total_count), (total_count-ref_count)/float(total_count)) < .01:
			total_monoallelic = total_monoallelic + 1
	return float(total_monoallelic)/total_het

def filter_ase_matrix(unfiltered_ase_file, filtered_ase_file, fraction_samples, fraction_monoallelic):
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
		gene_id = site_id.split(':')[1]
		counts = np.asarray(data[1:])
		num_samples = len(counts)
		het_indis = np.where(counts != 'NA')[0]
		if gene_id == 'NULL':
			continue
		if len(het_indis)/float(num_samples) < fraction_samples:
			continue
		fraction_monoallelic = compute_fraction_monoallelic(counts[het_indis])
		if fraction_monoallelic > .5:
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

#Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
def add_gene_to_chromosome_object(chromosome, start, end, gene_name):
    for pos in range(start, end+1):
        if chromosome[pos] == 'NULL':  # No genes already mapped to this position
            chromosome[pos] = gene_name
        else:  # at least one gene already mapped to this position
            chromosome[pos] = chromosome[pos] + ',' + gene_name
    return chromosome

def map_het_sites_to_genes(ase_site_file, gene_mapped_ase_site_file, gencode_gene_annotation_file):
	t = open(gene_mapped_ase_site_file, 'w')
	t.write('chr\tasePos\tGene_id\n')
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
				continue
			site_chrom = data[0]
			if site_chrom != str(chrom_num):
				continue
			site_pos = int(data[1])
			associated_genes = chromosome[site_pos]
			if associated_genes == 'NULL':
				t.write(data[0] + '\t' + data[1] + '\t' 'NULL\n')
			else:
				unique_genes = np.unique(associated_genes.split(','))
				if len(unique_genes) == 1:
					t.write(data[0] + '\t' + data[1] + '\t' + unique_genes[0] + '\n')
				else:
					t.write(data[0] + '\t' + data[1] + '\t' 'NULL\n')
		f.close()
	t.close()

def generate_raw_ase_matrix(ase_input_data_dir, gene_mapped_ase_site_file, unfiltered_ase_file):
	# First create mapping from site id to gene
	site_to_gene_mapping = {}
	f = open(gene_mapped_ase_site_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		site_id = 'chr' + data[0] + '_' + data[1]
		gene_id = data[2]
		if site_id in site_to_gene_mapping:
			print('assumption error')
			pdb.set_trace()
		site_to_gene_mapping[site_id] = gene_id
	f.close()
	# Next get ordered sites
	site_file = ase_input_data_dir + 'ase_locus_annot.txt'
	ordered_sites = []
	f = open(site_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		site_id = 'chr' + data[0] + '_' + data[1]
		gene_id = site_to_gene_mapping[site_id]
		full_site_id = site_id + ':' + gene_id
		ordered_sites.append(full_site_id)
	f.close()
	# Load in input ase data
	ref_file = ase_input_data_dir + 'ref.txt'
	alt_file = ase_input_data_dir + 'alt.txt'
	het_file = ase_input_data_dir + 'het.txt'

	ref_full = np.loadtxt(ref_file, dtype=str)
	alt_full = np.loadtxt(alt_file, dtype=str)
	het_full = np.loadtxt(het_file, dtype=str)
	ref_counts = np.transpose(ref_full[:, 1:].astype(int))
	alt_counts = np.transpose(alt_full[:, 1:].astype(int))
	het_counts = np.transpose(het_full[:, 1:].astype(int))
	sample_names = ref_full[:, 0]
	sample_names2 = alt_full[:, 0]
	if np.array_equal(sample_names, sample_names2) == False:
		print('assumptino error')
		pdb.set_trace()
	num_sites = ref_counts.shape[0]
	num_samples = ref_counts.shape[1]
	# Print to output file
	t = open(unfiltered_ase_file, 'w')
	# print header
	t.write('site_id\t' + '\t'.join(sample_names) + '\n')
	for site_num in range(num_sites):
		site_id = ordered_sites[site_num]
		site_arr = []
		for sample_num in range(num_samples):
			ref_instance = ref_counts[site_num, sample_num]
			alt_instance = alt_counts[site_num, sample_num]
			het_instance = het_counts[site_num, sample_num]
			total_instance = ref_instance + alt_instance
			if total_instance == 0 or het_instance == 0:
				site_arr.append('NA')
			else:
				stringer = str(ref_instance) + '/' + str(total_instance)
				site_arr.append(stringer)
		site_arr = np.asarray(site_arr)
		t.write(site_id + '\t' + '\t'.join(site_arr) + '\n')
	t.close()


def generate_annotated_sample_names(unfiltered_ase_file, dgn_technical_covariates, dgn_biological_covariates, annotated_sample_names_file):
	ase_sample_names = []
	f = open(unfiltered_ase_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			ase_sample_names = np.asarray(data[1:])
			continue
	f.close()
	f = open(dgn_technical_covariates)
	sample_to_technical_cov = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 36:
			print('assumptionerroro')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			technical_covariate_header = np.asarray(data[1:])
			continue
		sample_id = data[0]
		tech_cov = np.asarray(data[1:])
		sample_to_technical_cov[sample_id] = tech_cov
	f.close()
	f = open(dgn_biological_covariates)
	sample_to_biological_cov = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 40:
			print('assumtpoineroror')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			biological_covariate_header = np.asarray(data[1:])
			continue
		sample_id = data[0]
		biological_cov = np.asarray(data[1:])
		sample_to_biological_cov[sample_id] = biological_cov
	f.close()
	t = open(annotated_sample_names_file, 'w')
	# Print header
	t.write('sample_id\t' + '\t'.join(technical_covariate_header) + '\t' + '\t'.join(biological_covariate_header) + '\n')
	for sample_id in ase_sample_names:
		t.write(sample_id + '\t' + '\t'.join(sample_to_technical_cov[sample_id]) + '\t' + '\t'.join(sample_to_biological_cov[sample_id]) + '\n')
	t.close()



ase_input_data_dir = sys.argv[1]
dgn_technical_covariates = sys.argv[2]
dgn_biological_covariates = sys.argv[3]
gencode_gene_annotation_file = sys.argv[4]
output_root = sys.argv[5]

# Map ASE Heterozygous sites to genes
ase_site_file = ase_input_data_dir + 'ase_locus_annot.txt'
gene_mapped_ase_site_file = output_root + 'ase_locus_annot_gene_mapped.txt'
# map_het_sites_to_genes(ase_site_file, gene_mapped_ase_site_file, gencode_gene_annotation_file)



# Generate unfiltered ase matrix 
# columns are rna-seq samples
# rows are het sites
unfiltered_ase_file = output_root + 'raw_ase_counts.txt'
#generate_raw_ase_matrix(ase_input_data_dir, gene_mapped_ase_site_file, unfiltered_ase_file)


# Generate annotated sample names
annotated_sample_names_file = output_root + 'annotated_sample_names.txt'
generate_annotated_sample_names(unfiltered_ase_file, dgn_technical_covariates, dgn_biological_covariates, annotated_sample_names_file)


# Generate Filtered ase matrix 
# columns are rna-seq samples
# rows are het sites

# Require site to have at least min_samples het individuals
fraction_samples = .35
fraction_monoallelic = .5

filtered_ase_file = output_root + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '.txt'
#filter_ase_matrix(unfiltered_ase_file, filtered_ase_file, fraction_samples, fraction_monoallelic)



filtered_one_site_per_gene_ase_file = output_root + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_one_site_per_gene.txt'
#filter_ase_matrix_to_one_site_per_gene(filtered_ase_file, filtered_one_site_per_gene_ase_file)


filtered_one_site_per_gene_min_allele_ase_file = output_root + 'filtered_ase_counts_' + str(fraction_samples) + '_' + str(fraction_monoallelic) + '_one_site_per_gene_min_allele.txt'
#make_min_allele_ase_file(filtered_one_site_per_gene_ase_file, filtered_one_site_per_gene_min_allele_ase_file)

