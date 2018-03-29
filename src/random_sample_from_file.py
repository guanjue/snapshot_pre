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
### plot violin
def random_sample_from_file(filename, random_seed, sample_num, outputname):
	### read input file
	data = read2d_array(filename, str)

	### random sampling
	np.random.seed(random_seed)
	idx = np.random.randint(low = 0, high = data.shape[0], size=sample_num)
	data_sample = data[idx,:]

	### write output file
	write2d_array(data_sample, outputname)

############################################################################
#time python random_sample_from_file.py -i 200_noblack.11_22_2017.bed -o 200_noblack.sample.bed -s 2017 -n 200000

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:s:n:")
	except getopt.GetoptError:
		print 'time python index_label2meansig.py -i input_file_list -o outputname -s random_seed -n random_sampling_number'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python index_label2meansig.py -i input_file_list -o outputname -s random_seed -n random_sampling_number'
			sys.exit()
		elif opt=="-i":
			filename=str(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())
		elif opt=="-s":
			random_seed=int(arg.strip())
		elif opt=="-n":
			sample_num=int(arg.strip())
	### run random sampling
	random_sample_from_file(filename, random_seed, sample_num, outputname)

if __name__=="__main__":
	main(sys.argv[1:])

