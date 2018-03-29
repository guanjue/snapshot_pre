import os
import numpy as np
from subprocess import call
################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################
### ideas_matrix2ideas_bed
def index_label2meansig(bed_signal_matrix_file, bed_sig_col, significant_index_label, significant_index_label_col, outputname):
	### read od_index_label matrix & significant label
	bed_signal_matrix = read2d_array(bed_signal_matrix_file, str)
	significant_index_label = read2d_array(significant_index_label, str)
	### signif index_label to dict
	significant_index_label_dict = {}
	for records in significant_index_label:
		label_tmp = records[significant_index_label_col-1]
		significant_index_label_dict[label_tmp] = 0
	### get insignif label
	insignif_label = ''
	print(bed_signal_matrix[0][-1])
	eg = bed_signal_matrix[0][-1].split('_')
	for i in range(0,len(eg)-1):
		insignif_label = insignif_label + 'X_'
	insignif_label = insignif_label + 'X'
	significant_index_label_dict[insignif_label] = 0
	### od to signif
	bed_signal_matrix_dict = {}
	for records in bed_signal_matrix:
		index_set_label = records[-1]
		#print(index_set_label)
		if index_set_label in significant_index_label_dict:
			if index_set_label in bed_signal_matrix_dict:
				bed_signal_matrix_dict[index_set_label].append( records[(bed_sig_col-1):len(records)-1] )
			else:
				bed_signal_matrix_dict[index_set_label] = [ records[(bed_sig_col-1):len(records)-1] ]
		else:
			index_set_label = insignif_label
			if index_set_label in bed_signal_matrix_dict:
				bed_signal_matrix_dict[index_set_label].append( records[(bed_sig_col-1):len(records)-1] )
			else:
				bed_signal_matrix_dict[index_set_label] = [ records[(bed_sig_col-1):len(records)-1] ]			
	### write output 
	r1=open(outputname,'w')
	for records in bed_signal_matrix_dict:
		sig_matrix = np.array(bed_signal_matrix_dict[records], dtype=float)
		#print(sig_matrix.shape)
		sig_matrix_mean = np.mean(sig_matrix, axis=0)
		print(sig_matrix_mean)
		r1.write(str(records)+'\t')
		for i in range(0,len(sig_matrix_mean)-1):
			r1.write(str(sig_matrix_mean[i])+'\t')
		r1.write(str(sig_matrix_mean[len(sig_matrix_mean)-1])+'\n')
	r1.close()

	call('sort -k1,1 ' + outputname + ' > ' + outputname+'.meansig.indexsort.txt', shell=True)
	call('rm ' + outputname, shell=True)
############################################################################
#time python index_label2meansig.py -i B_SPL.signal_indexlabel_matrix.signif.txt -a 4 -s all5cell.binary_matrix.png.binary.counts.thresh.bed -c 2 -o B_SPL.meansig.signif.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:a:s:c:o:")
	except getopt.GetoptError:
		print 'time python index_label2meansig.py -i bed_signal_matrix_file -a bed_sig_col -s significant_index_label -c significant_index_label_col -o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python index_label2meansig.py -i od_index_label_table -a bed_sig_col -s significant_index_label -c significant_index_label_col -o outputname'
			sys.exit()
		elif opt=="-i":
			bed_signal_matrix_file=str(arg.strip())
		elif opt=="-a":
			bed_sig_col=int(arg.strip())
		elif opt=="-s":
			significant_index_label=str(arg.strip())		
		elif opt=="-c":
			significant_index_label_col=int(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())

	index_label2meansig(bed_signal_matrix_file, bed_sig_col, significant_index_label, significant_index_label_col, outputname)

if __name__=="__main__":
	main(sys.argv[1:])