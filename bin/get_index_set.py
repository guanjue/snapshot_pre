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

################################################################################################
### write index sets 
def get_index_set(inputmatrix, given_column_order, given_signal_level, sorted_index_output, sorted_index_set_output):
	### read given cell type order
	column_oder = read2d_array(given_column_order,str)
	order = np.array(np.transpose(column_oder[:,0]), dtype=int)
	header = np.array([column_oder[:,1]])
	print('cell type order:')
	print(header)

	### read given sorting signal level
	signal_level_range = read2d_array(given_signal_level, str)
	signal_level_range = np.array(signal_level_range[1:,0:3], dtype=float)
	print('signal level range:')
	print(signal_level_range)

	### read primary matrix
	prime_mark0 = read2d_array(inputmatrix,str)
	prime_mark_id = np.transpose([prime_mark0[1:,0]])

	### sort cell types by given cell type order
	prime_mark_info_sorted = np.array(prime_mark0[1:,order], dtype=float)

	### change signal to user defined signal levels
	for levels in signal_level_range:
		prime_mark_info_sorted[((prime_mark_info_sorted>=levels[1]) * (prime_mark_info_sorted<levels[2]))] = levels[0]

	### sort rows based on index pattern
	for i in range(0,prime_mark_info_sorted.shape[1])[::-1]:
		### get rows order based on each column at a time
		row_order = prime_mark_info_sorted[:,i].argsort(kind='mergesort')
		### sort pattern matrix
		prime_mark_info_sorted = prime_mark_info_sorted[row_order,:]
		### sort DNA region id
		prime_mark_id = prime_mark_id[row_order,:]

	### add header
	prime_mark_id_sorted = np.concatenate((np.array([['name']]), prime_mark_id),axis = 0)
	prime_mark_info_sorted_info = np.concatenate((header, prime_mark_info_sorted), axis = 0)

	### add row names
	prime_mark_info_sorted_matrix = np.concatenate((prime_mark_id_sorted, prime_mark_info_sorted_info), axis = 1)

	### write the all sorted DNA region binary pattern (index)
	write2d_array(prime_mark_info_sorted_matrix, sorted_index_output+'.all.txt')

	### keep the peaks that are called in both replicates in at least one cell type
	### split header from the rest of the matrix (so the matrix can be converted to float)
	prime_mark_info_sorted_matrix_reliable = [prime_mark_info_sorted_matrix[0]]
	for i in range(1, prime_mark_info_sorted_matrix.shape[0]):
		### get pattern & 2float
		pattern = np.array(prime_mark_info_sorted_matrix[i,1:], dtype=float)
		### check if called in both replicates in at least one cell type
		info = np.sum(pattern)
		if info > 0:
			prime_mark_info_sorted_matrix_reliable.append(prime_mark_info_sorted_matrix[i,:])
	write2d_array(prime_mark_info_sorted_matrix_reliable, sorted_index_output)

	### get index sets
	index_set_dict = {}
	index_set = []
	for index in prime_mark_info_sorted:
		index_merge = ''
		for i in index:
			index_merge = index_merge+'_'+str(i)
		### get the number of DNA regions in index set
		if not (index_merge in index_set_dict):
			index_set_dict[index_merge] = 1
			index_set.append(index_merge)
		else:
			index_set_dict[index_merge] = index_set_dict[index_merge]+1

	### write index set
	def write_index_set(output_name, index_set_dict, index_set_array, header):
		result = open(output_name,'w')
		### write header
		result.write('name'+'\t')
		for i in range(0,header.shape[1]-1):
			result.write(header[0,i]+'\t')
		result.write(str(header[0,header.shape[1]-1])+'\n')

		### write index set info
		for pattern in index_set_array:
			pattern_num = index_set_dict[pattern]
			binary_label = pattern.split('_')
			if np.sum(np.array(binary_label[1:], dtype=float)) > 0: ### reomve the empty index set (all 0)
				result.write(pattern[1:]+'\t')
				### label times the number of DNA regions
				for i in range(1, len(binary_label)-1):
					result.write(str(float(binary_label[i]) * pattern_num)+'\t') 
				result.write(str(float(binary_label[len(binary_label)-1]) * pattern_num)+'\n')
		result.close()

	write_index_set(sorted_index_set_output, index_set_dict, index_set, header)

############################################################################
#time python get_index_set.py -i celltype.binary_pattern.txt -r celltype.order.txt -l signal_level_range.txt -f celltype.binary_pattern.sorted.txt -s celltype.index_set.sorted.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:r:l:f:s:")
	except getopt.GetoptError:
		print 'time python get_index_set.py -i inputmatrix -r celltype_order -l signal_level_range -f index_output -s index_set_output'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python get_index_set.py -i inputmatrix -r celltype_order -l signal_level_range -f index_output -s index_set_output'
			sys.exit()
		elif opt=="-i":
			inputmatrix=str(arg.strip())
		elif opt=="-r":
			given_column_order=str(arg.strip())
		elif opt=="-l":
			given_signal_level=str(arg.strip())		
		elif opt=="-f":
			sorted_index_output=str(arg.strip())
		elif opt=="-s":
			sorted_index_set_output=str(arg.strip())

	get_index_set(inputmatrix, given_column_order, given_signal_level, sorted_index_output, sorted_index_set_output)

if __name__=="__main__":
	main(sys.argv[1:])
