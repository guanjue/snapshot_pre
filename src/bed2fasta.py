import os
import numpy as np

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used, split_mark):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(split_mark)]
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


def bed2fasta(input_bed_file, sequence_col, outputname):
	data = read2d_array(input_bed_file, str, '\t')
	seq = data[:,sequence_col-1]
	bed = data[:,0:3]
	fasta_array = []
	for b,s in zip(bed, seq):
		name = '>'+b[0]+':'+b[1]+'-'+b[2]
		fasta_seq = s.replace(',','')
		fasta_array.append([name])
		fasta_array.append([fasta_seq])

	write2d_array(fasta_array, outputname)

############################################################################
#time python index_label2meansig.py -i B_SPL.signal_indexlabel_matrix.signif.txt -a 4 -s all5cell.binary_matrix.png.binary.counts.thresh.bed -c 2 -o B_SPL.meansig.signif.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:c:o:")
	except getopt.GetoptError:
		print 'time python bed2fasta.py -i input_bed_file -c sequence_col -o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python bed2fasta.py -i input_bed_file -c sequence_col -o outputname'
			sys.exit()
		elif opt=="-i":
			input_bed_file=str(arg.strip())
		elif opt=="-c":
			sequence_col=int(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())

	bed2fasta(input_bed_file, sequence_col, outputname)

if __name__=="__main__":
	main(sys.argv[1:])