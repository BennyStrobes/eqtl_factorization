import numpy as np 
import os
import sys
import pdb


def extract_individuals_from_file(input_individual_file):
	f = open(input_individual_file)
	indis = []
	for line in f:
		line = line.rstrip()
		indis.append(line)
	return np.asarray(indis)

def extract_snp_names_from_file(input_snp_name_file, chrom_string):
	f = open(input_snp_name_file)
	snp_names = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if data[0] != chrom_string:
			print('assumption error!')
			pdb.set_trace()
		pos = int(data[1])
		if pos < 0:
			print('assumption error')
			pdb.set_trace()
		if pos > 250000000:
			print('assumption error')
			pdb.set_trace()
		snp_name = data[0] + ':' + data[1]
		snp_names.append(snp_name)
	return np.asarray(snp_names)

def generate_chromosome_specific_genotype(chrom_num, input_genotype_file, input_individual_file, input_snp_name_file, output_genotype_file):
	# Extract input data
	indis = extract_individuals_from_file(input_individual_file)
	snp_names = extract_snp_names_from_file(input_snp_name_file, str(chrom_num))
	genotype = np.transpose(np.loadtxt(input_genotype_file))
	genotype = genotype[1:,:]
	# Open output file 
	t = open(output_genotype_file, 'w')
	# Print header of output file
	t.write('snp_id\t' + '\t'.join(indis) + '\n')
	# Print each snp to outputfile
	num_snps = len(snp_names)
	for snp_num in range(num_snps):
		snp_name = snp_names[snp_num]
		genotype_vec = genotype[snp_num, :]
		t.write(snp_name + '\t' + '\t'.join(genotype_vec.astype(str)) + '\n')
	t.close()


#####################
# Command Line Args
#####################
genotype_dir = sys.argv[1]
pre_processed_data_dir = sys.argv[2]


for chrom_num in range(1,23):
	print(chrom_num)
	input_root = genotype_dir + 'chr' + str(chrom_num) + '/'
	input_genotype_file = input_root + 'combined.chr' + str(chrom_num) + '.common.012'
	input_individual_file = input_root + 'combined.chr' + str(chrom_num) + '.common.012.indv'
	input_snp_name_file = input_root + 'combined.chr' + str(chrom_num) + '.common.012.pos'
	output_genotype_file = pre_processed_data_dir + 'genotype_chr_' + str(chrom_num) + '.txt'
	generate_chromosome_specific_genotype(chrom_num, input_genotype_file, input_individual_file, input_snp_name_file, output_genotype_file)


