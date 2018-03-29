import os
import numpy as np
from subprocess import call
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
### plot violin
def plot_violin(input_file_list, outputname, log2, small_num, lowerlim, upperlim):
	input_file_list = open(input_file_list, 'r')
	signal_track_list = []
	filename_list = []
	for records in input_file_list:
		filename = records.split()[0]
		filename_list.append(records.split()[1])
		### read file
		signal_track = read2d_array(filename, float)
		### select data by lower & upper limit
		signal_track = signal_track[signal_track[:,0]>=lowerlim,:]
		signal_track = signal_track[signal_track[:,0]<=upperlim,:]
		### log2 transform
		if log2=='T':
			signal_track = np.log2(signal_track + small_num)
		print(signal_track[:,0])
		signal_track_list.append(signal_track[:,0])

	### plot violin plot
	pos = range(1,(len(filename_list))*5,5)
	print(len(pos))
	print(len(signal_track_list))
	print('plot violinplot of index:' + outputname)
	plt.figure(figsize=(10,10), dpi=200)
	plt.violinplot(signal_track_list, pos, points=50, widths=4, showmeans=True, showextrema=False, showmedians=False)
	plt.xticks(pos, filename_list, rotation='vertical')
	plt.savefig(outputname + '.violin.pdf')

############################################################################
#time python plot_violin.py -i atac_list.txt -o atac_list -l T -s 2 -l 4 -u 100

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:l:s:d:u:")
	except getopt.GetoptError:
		print 'time python index_label2meansig.py -i input_file_list -o outputname -l log2 -s small_num -d lowerlim -u upperlim'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python index_label2meansig.py -i input_file_list -o outputname -l log2 -s small_num -d lowerlim -u upperlim'
			sys.exit()
		elif opt=="-i":
			input_file_list=str(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())
		elif opt=="-l":
			log2=str(arg.strip())
		elif opt=="-s":
			small_num=float(arg.strip())		
		elif opt=="-d":
			lowerlim=float(arg.strip())		
		elif opt=="-u":
			upperlim=float(arg.strip())		


	plot_violin(input_file_list, outputname, log2, small_num, lowerlim, upperlim)

if __name__=="__main__":
	main(sys.argv[1:])

