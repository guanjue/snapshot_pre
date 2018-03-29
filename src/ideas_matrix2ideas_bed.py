import os
import numpy as np
################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	### get colnames
	colnames = data.readline()
	colnames = [x.strip() for x in colnames.split(' ')]
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(' ')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return {'data0': data0, 'colnames': colnames}

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
def ideas_matrix2ideas_bed(input_matrix, bed_start_col, bed_end_col, sig_start_col, sig_end_col, outputname):
	### read input matrix and colnames
	input_matrix = read2d_array(input_matrix, str)
	### get signal matrix and colnames
	input_matrix_bedsig = input_matrix['data0']
	input_matrix_colnames = input_matrix['colnames'][(sig_start_col-1):(sig_end_col)]
	### bed matrix & signal matrix
	input_matrix_bed = input_matrix_bedsig[:, (bed_start_col-1):(bed_end_col)]
	input_matrix_sig = input_matrix_bedsig[:, (sig_start_col-1):(sig_end_col)]
	### 
	for i in range(0,len(input_matrix_colnames)):
		tmp_sig = input_matrix_sig[:,i]
		tmp_sig = tmp_sig.reshape(tmp_sig.shape[0],1)
		ct_matrix_bed_sig = np.concatenate((input_matrix_bed, tmp_sig), axis = 1)
		celltypename = input_matrix_colnames[i]
		write2d_array(ct_matrix_bed_sig, outputname+'.'+celltypename+'.bed')

############################################################################
#time python ideas_matrix2ideas_bed.py -i tmp.txt -a 2 -b 4 -c 5 -d 24 -o tmp

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:a:b:c:d:o:")
	except getopt.GetoptError:
		print 'time python ideas_matrix2ideas_bed.py -i input_matrix -a bed_start_col -b bed_end_col -c sig_start_col -d sig_end_col -o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python ideas_matrix2ideas_bed.py -i input_matrix -a bed_start_col -b bed_end_col -c sig_start_col -d sig_end_col -o outputname'
			sys.exit()
		elif opt=="-i":
			input_matrix=str(arg.strip())
		elif opt=="-a":
			bed_start_col=int(arg.strip())
		elif opt=="-b":
			bed_end_col=int(arg.strip())		
		elif opt=="-c":
			sig_start_col=int(arg.strip())
		elif opt=="-d":
			sig_end_col=int(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())

	ideas_matrix2ideas_bed(input_matrix, bed_start_col, bed_end_col, sig_start_col, sig_end_col, outputname)

if __name__=="__main__":
	main(sys.argv[1:])