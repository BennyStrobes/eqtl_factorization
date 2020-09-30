import numpy as np 
import os
import sys
import pdb




batch_name = sys.argv[1]
master_barcode_file = sys.argv[2]
batch_barcode_file = sys.argv[3]



batch_barcodes = {}

f = open(master_barcode_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if len(data) != 7:
		print('assumption eroror')
		pdb.set_trace()
	if head_count == 0:
		head_count = head_count + 1
		continue
	barcode_info = data[0].split('_')
	if barcode_info[1] == batch_name:
		batch_barcodes[barcode_info[2]] = 1
f.close()
t = open(batch_barcode_file,'w')
for barcode in batch_barcodes.keys():
	t.write(barcode + '\n')
t.close()
