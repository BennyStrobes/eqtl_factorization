import numpy as np 
import os
import sys
import pdb


def make_sure_files_exist(output_root, total_jobs, suffix):
	booly = True
	for job_number in range(total_jobs):
		file_name = output_root + str(job_number) + '_' + str(total_jobs) + suffix
		if os.path.isfile(file_name) == False:
			print(file_name)
			booly = False
	return booly

def merge_parallelized_results(output_root, suffix, total_jobs):
	to_run = make_sure_files_exist(output_root, total_jobs, suffix)
	if to_run == False:
		print('Missing required input files. Please re-evaluate :)')
		return
	# Open output (merged result) file handle
	t = open(output_root + 'merged' + suffix, 'w')
	# Loop through parrallelized jobs
	for job_number in range(total_jobs):
		file_name = output_root + str(job_number) + '_' + str(total_jobs) + suffix
		# Open file for one job
		f = open(file_name)
		# To identify header
		head_count = 0
		# Stream file from one job
		for line in f:
			line = line.rstrip()
			# HEADER
			if head_count == 0:
				head_count = head_count + 1
				# Print header if this the first job
				if job_number == 0:
					t.write(line + '\n')
				continue
			# Standard line
			t.write(line + '\n')
		f.close()
		# Delete file from single job
		os.system ('rm ' + file_name)
	t.close()



output_root = sys.argv[1]
total_jobs = int(sys.argv[2])



merge_parallelized_results(output_root, ".txt", total_jobs)
