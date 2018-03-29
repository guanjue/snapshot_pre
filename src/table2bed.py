import os
import numpy as np

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

def table2bed(input_table_file, bedcol, binary_col_st, binary_col_ed, sep_sign):
	### read input
	data_all = read2d_array(input_table_file, str)
	header_all = data_all[0,:]
	bed = data_all[1:,0:bedcol]
	binary = data_all[1:,binary_col_st-1:binary_col_ed]
	header_binary = header_all[binary_col_st-1:binary_col_ed]

	ct_dict = {}
	ct_array = []
	for i in range(0,header_binary.shape[0]):
		ct = header_binary[i].split(sep_sign)[0]
		print(ct)
		if not (ct in ct_dict):
			ct_dict[ct] = [ binary[:,i] ]
			ct_array.append(ct)
		else:
			ct_dict[ct].append(binary[:,i])

	ct_binary_dict = {}
	for ct in ct_array:
		binary_matrix = np.array(ct_dict[ct], dtype='float')
		binary_matrix_rowsum = np.sum(binary_matrix, axis=0)
		### get pk in the ct
		binary_matrix_used = binary_matrix_rowsum != 0
		### merge bed and binary label
		print(bed.shape)
		print(np.reshape(binary_matrix_rowsum, (binary_matrix_rowsum.shape[0],1)).shape)
		bed_tmp_all = np.concatenate((bed, np.reshape(binary_matrix_rowsum, (binary_matrix_rowsum.shape[0],1))), axis=1)
		### only keep pk
		bed_tmp_ct = bed_tmp_all[binary_matrix_used,:]
		### add to dict
		ct_binary_dict[ct] = bed_tmp_ct

	### write output
	for ct in ct_array:
		bed = ct_binary_dict[ct]
		write2d_array(bed,ct+'.pk.bed')


############################################################################
#time python table2bed.py -i tmp.txt -b 4 -s 36 -s 66 -p _

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:b:s:e:p:")
	except getopt.GetoptError:
		print 'time python table2bed.py -i input_table_file -b 4 -s binary_start_col -e binary_end_col -p sep_sign'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python table2bed.py -i input_table_file -b 4 -s binary_start_col -e binary_end_col -p sep_sign'
			sys.exit()
		elif opt=="-i":
			input_table_file=str(arg.strip())
		elif opt=="-b":
			bedcol=int(arg.strip())
		elif opt=="-s":
			binary_col_st=int(arg.strip())
		elif opt=="-e":
			binary_col_ed=int(arg.strip())		
		elif opt=="-p":
			sep_sign=str(arg.strip())	

	table2bed(input_table_file, bedcol, binary_col_st, binary_col_ed, sep_sign)

if __name__=="__main__":
	main(sys.argv[1:])

