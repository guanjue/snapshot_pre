import os
import numpy as np
from subprocess import call
from scipy import stats
from scipy.stats import t, nbinom, poisson
from collections import Counter
from sklearn.cluster import AgglomerativeClustering, KMeans
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
### get convert bedtools window output to matrix of pk and intersect ideas label info
def ideas_label_info(input_bedtools_window, id_col, lb_col, pk_col, ideas_col):
	data_info0=open(input_bedtools_window, 'r')
	### read DNA region orders
	data_info=[]
	for records in data_info0:
		tmp=[x.strip() for x in records.split('\t')]
		### get intersect region; midpoint dist; TF peak length
		if ((int(tmp[ideas_col-1]) - int(tmp[pk_col-1]))>=0) and ((int(tmp[ideas_col]) - int(tmp[pk_col]))<=0) :
			### IDEAS Bin >= pk region
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[ideas_col])-int(tmp[ideas_col-1]), (float(tmp[ideas_col])+float(tmp[ideas_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[ideas_col])-int(tmp[ideas_col-1]) ]
		elif ((int(tmp[ideas_col-1]) - int(tmp[pk_col-1]))<0) and ((int(tmp[ideas_col]) - int(tmp[pk_col]))>0) :
			### IDEAS Bin < pk region
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[pk_col])-int(tmp[pk_col-1]), (float(tmp[ideas_col])+float(tmp[ideas_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[ideas_col])-int(tmp[ideas_col-1]) ]
		elif ((int(tmp[ideas_col-1]) - int(tmp[pk_col-1]))<0) and ((int(tmp[ideas_col]) - int(tmp[pk_col]))<=0) :
			### IDEAS Bin upstream < pk region upstream & IDEAS Bin downstream <= pk region downstream
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[ideas_col])-int(tmp[pk_col-1]), (float(tmp[ideas_col])+float(tmp[ideas_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[ideas_col])-int(tmp[ideas_col-1]) ]
		elif ((int(tmp[ideas_col-1]) - int(tmp[pk_col-1]))>=0) and ((int(tmp[ideas_col]) - int(tmp[pk_col]))>0) :
			### IDEAS Bin upstream >= pk region upstream & IDEAS Bin downstream > pk region downstream
			tmp_vec = [ tmp[id_col-1], tmp[lb_col-1], int(tmp[pk_col])-int(tmp[ideas_col-1]), (float(tmp[ideas_col])+float(tmp[ideas_col-1])-float(tmp[pk_col])-float(tmp[pk_col-1]))/2, int(tmp[ideas_col])-int(tmp[ideas_col-1]) ]
		data_info.append(tmp_vec)
	data_info0.close()
	###### return output dict
	return (data_info)

################################################################################################
### get peak's ideas labels
def get_cRE_ideas_state(data_info_matrix, id_col, lb_col, cover_col, middist_col, ideaslen, bed_od_file, bed_od_idcol, outputname):
	### read DNA region IDEAS state info matrix
	pk_id_list = []
	data_ideas1={}
	data_ideas1_maxcover={} ### coverage size
	data_ideas1_middist={} ### midpoint dist
	data_ideas1_statelen={} ### ideas state len
	### initialize problem counter
	k=0
	for info in data_info_matrix:
		pk_id = info[id_col-1]
		### creat pk_id_list for keeping the id order in output
		pk_id_list.append(pk_id)
		if not (pk_id in data_ideas1):
			data_ideas1[pk_id] = info[lb_col-1]
			data_ideas1_maxcover[pk_id] = info[cover_col-1]
			data_ideas1_middist[pk_id] = info[middist_col-1]
			data_ideas1_statelen[pk_id] = info[ideaslen-1]
		elif info[cover_col-1] > data_ideas1_maxcover[pk_id]:
			### if multiple cover; select the highest covering state
			data_ideas1[pk_id] = info[lb_col-1]
			data_ideas1_maxcover[pk_id] = info[cover_col-1]
			data_ideas1_middist[pk_id] = info[middist_col-1]
			data_ideas1_statelen[pk_id] = info[ideaslen-1]
		elif info[cover_col-1] == data_ideas1_maxcover[pk_id]: ### if 2 states cover the same region with same length
			if info[middist_col-1] < data_ideas1_middist[pk_id]: 
				### if cover the same; check mid point distance
				data_ideas1[pk_id] = info[lb_col-1]
				data_ideas1_maxcover[pk_id] = info[cover_col-1]
				data_ideas1_middist[pk_id] = info[middist_col-1]
				data_ideas1_statelen[pk_id] = info[ideaslen-1]
			elif info[middist_col-1] == data_ideas1_middist[pk_id]: ### if 2 states cover the same region with same length; with same midpoint dist
				if info[ideaslen-1] < data_ideas1_statelen[pk_id]:
					### if cover same & mid point distance same; check state len 
					data_ideas1[pk_id] = info[lb_col-1]
					data_ideas1_maxcover[pk_id] = info[cover_col-1]
					data_ideas1_middist[pk_id] = info[middist_col-1]
					data_ideas1_statelen[pk_id] = info[ideaslen-1]
				else: ### if 2 states cover the same region with same length; with same midpoint dist; with same state length ...give attention!
					k=k+1
					print('problem!')
					print(k)
	### read original bed file to get the pk id list
	bed_od_file=open(bed_od_file,'r')
	bed_od_id_list = []
	for records in bed_od_file:
		bed_od_id_list.append(records.split()[bed_od_idcol-1])
	bed_od_file.close()
	### write ideas label output
	result=open(outputname,'w')
	for pkid in bed_od_id_list:
		if pkid in data_ideas1:
			tmp=data_ideas1[pkid]
			result.write(tmp+'\n')
		else:
			tmp=records
			result.write('NA'+'\n')
	result.close()

################################################################################################
### get index/signal matrix
def get_mark_matrix(peak_bed, peak_bed_colnum, mark_list, output_file, signal_col, method, sort):
	### sort input bed files
	sort_bed_file = peak_bed + '.sort.bed'
	call('gsort -k1,1 -V -s -k2,2n ' + peak_bed + ' > ' + sort_bed_file, shell=True)
	call('cp ' + sort_bed_file + ' ' + output_file, shell=True)
	### generate index mark matrix
	mark_list_vec = open(mark_list, 'r')
	celltype_list = []
	for mark_bed in mark_list_vec:	
		tmp = [x.strip() for x in mark_bed.split('\t')]
		### read bianry label file list
		mark_bed_file = tmp[0]
		print(mark_bed_file)
		### add cell type name to cell type list
		celltype_list.append(tmp[1])
		### sort bianry label bed files
		if sort == 'T':
			call('gsort -k1,1 -V -s -k2,2n ' + mark_bed_file + ' > ' + mark_bed_file+'.sort.bed', shell=True)
		else:
			call('cp ' + mark_bed_file + ' ' + mark_bed_file+'.sort.bed', shell=True)
		### use bedtools to generate the index/signal matrix
		if method == 'intersect':
			### used bedtools intersect to get the binary label of each peak
			call('bedtools intersect -c -a ' + sort_bed_file + ' -b ' + mark_bed_file+'.sort.bed' + ' > ' + mark_bed_file+'.tmp01.txt', shell=True)
		elif method == 'map':
			### used bedtools map to get the average signal of each peak
			call('bedtools map -c ' + str(signal_col) + ' -null 0 -F 0.5 -o mean -a ' + sort_bed_file + ' -b ' + mark_bed_file+'.sort.bed' + ' > ' + mark_bed_file+'.tmp01.txt', shell=True)
		elif method == 'window':
			### used bedtools map to get the average signal of each peak
			call('bedtools window -a ' + sort_bed_file + ' -b ' + mark_bed_file+'.sort.bed' + ' -w 0 > ' + mark_bed_file+'.tmp01.txt', shell=True)
			### convert bedtools window output to matrix of pk and intersect ideas label info (intersect region; midpoint dist; TF peak length)
			data_info_matrix = ideas_label_info(mark_bed_file+'.tmp01.txt', 4, 8, 2, 6)
			### get peak's ideas labels based on intersect region; midpoint dist; TF peak length
			get_cRE_ideas_state(data_info_matrix, 1, 2, 3, 4, 5, sort_bed_file, 4, mark_bed_file+'.tmp01.txt')
		### cut the map number column
		call('cut -f'+ str(peak_bed_colnum+1) +" -d$'\t' " + mark_bed_file+'.tmp01.txt' + ' > ' + mark_bed_file+'.tmp02.txt', shell=True)
		### cbind to matrix
		call('paste ' + output_file + ' ' + mark_bed_file+'.tmp02.txt' + ' > ' + output_file+'.tmp.txt' + ' && mv ' + output_file+'.tmp.txt ' + output_file, shell=True)
		### remove tmp files
		call('rm ' + mark_bed_file+'.tmp01.txt' + ' ' + mark_bed_file+'.tmp02.txt' + ' ' + mark_bed_file+'.sort.bed', shell=True)
	mark_list_vec.close()

################################################################################################
### convert index matrix to index counts dict
def index_matrix2index_count_dict(index_matrix, index_matrix_start_col):
	### binary matrix to index X_X_X
	index_vector = []
	for vec in index_matrix:
		index_tmp = ''
		for i in range(index_matrix_start_col-1, len(vec)-1):
			index_tmp = index_tmp + vec[i] + '_'
		index_tmp = index_tmp + vec[len(vec)-1]
		index_vector.append( index_tmp )
	### index_vector 2 np array
	index_vector = np.array(index_vector)
	### index peak counts (dict)
	index_uniq_count_dict = Counter(index_vector)
	### index peak counts dict 2 index peak counts np array
	index_vector_count_vec = []
	for index in index_uniq_count_dict:
		index_vector_count_vec.append(index_uniq_count_dict[ index ])
	index_vector_count_vec = np.array(index_vector_count_vec)
	###### return output dict
	return {'index_vector': index_vector.reshape(index_vector.shape[0], 1), 'index_uniq_count_dict': index_uniq_count_dict, 'index_vector_count_vec': index_vector_count_vec}

################################################################################################
### get index set counts threshold
def index_count_thresh(count_vec, thesh):
	mean = np.mean(count_vec)
	std = np.std(count_vec)
	n = len(count_vec)
	t_stat = t.ppf(thesh, n)
	index_count_thresh = mean + t_stat * std / np.sqrt(n)
	###### return output dict
	return(index_count_thresh)

################################################################################################
### use while loop to select the threshold of index set counts
def select_index_set_counts_thresh(index_matrix, index_matrix_start_col, siglevel_counts):
	index_count_dict = index_matrix2index_count_dict(index_matrix, index_matrix_start_col)
	#print(index_count_dict)
	index_vector = index_count_dict['index_vector']
	index_uniq_count_dict = index_count_dict['index_uniq_count_dict']
	index_vector_count_vec = index_count_dict['index_vector_count_vec']
	print('calculating index peak counts thresh...')
	### initalize thresh 
	index_vector_count_vec = np.log2(index_vector_count_vec)
	index_vector_count_vec_select = index_vector_count_vec

	index_count_thresh_1 = max(index_vector_count_vec)
	index_count_thresh_2 = 0
	i = 1
	#while index_count_thresh_1 != index_count_thresh_2:
	for i in range(0,1):
		print('select index peak counts thresh - round: ' + str(i))
		#print('index_count_thresh_1: ' + str(index_count_thresh_1))
		#print('index_count_thresh_2: ' + str(index_count_thresh_2))
		i = i+1
		### use  significant level to find insignificant index counts 
		index_count_thresh_1 = index_count_thresh(index_vector_count_vec_select, siglevel_counts)
		index_vector_count_vec_select = index_vector_count_vec_select[index_vector_count_vec_select<index_count_thresh_1]
		index_count_thresh_2 = 2**index_count_thresh(index_vector_count_vec_select, siglevel_counts)
		print('index_count_thresh_1: ' + str(index_count_thresh_1))
		print('index_count_thresh_2: ' + str(index_count_thresh_2))
	index_count_thresh_2 = 200
	print('select index peak counts thresh: ' + str(index_count_thresh_2))
	### extract insignificant index names
	insig_index = []
	for index in index_uniq_count_dict:
		if index_uniq_count_dict[ index ] < index_count_thresh_2:
			insig_index.append( index )
	###### return output dict
	return { 'index_vector': index_vector, 'insig_index': insig_index, 'index_count_thresh': index_count_thresh_2, 'index_vector_count_vec': index_vector_count_vec }

################################################################################################
### calculating multiple variable norm density score
def mvn_density_score(signal_matrix_od, signal_matrix_start_col, log_signal, small_value, qda_round, index_vector, insig_index, scale):
	print('calculating multiple variable norm density score...')
	### initalize signal matrix for each index (dict)
	index_signal_matrix_dict = {}
	### initalize uniq index vector
	uniq_index = []
	### extract bed files
	signal_matrix_bed = signal_matrix_od[:,range(0, signal_matrix_start_col-1)]
	### extract signal matrix
	signal_matrix = signal_matrix_od[:,range(signal_matrix_start_col-1,signal_matrix_od.shape[1])]
	### convert string matrix to float matrix
	signal_matrix = signal_matrix.astype(float)
	### log transform the data
	if log_signal == 'T':
		signal_matrix = np.log2(signal_matrix+small_value)
	### scale
	if scale == 'T':
		signal_matrix_mean = np.mean(signal_matrix, axis=0)
		signal_matrix_std = np.std(signal_matrix, axis=0)
		signal_matrix = (signal_matrix -  signal_matrix_mean) / signal_matrix_std
	print('check scale...')
	print(np.mean(signal_matrix, axis=0))
	print(np.std(signal_matrix, axis=0))
	###############
	### QDA start
	for l in range(0, qda_round):
		print('Round: '+ str(l))
		### initialize the index name vector for all peaks
		if l==0:
			### in the 1st round, we use the orginal index to label the peaks
			index_vector_loop = index_vector
			index_name_vec_index_set = index_vector.reshape(index_vector.shape[0],1)
		else:
			### after 1st round, we use previous round generated index to label the peaks
			index_vector_loop = index_name_vec
		### convert index to index and insignificant index set
		print('filter insig_index...')
		if l == qda_round-1:
			index_vector_filter = []
		for i in range(0, len(index_vector_loop)):
			index = index_vector_loop[i, 0]
			signal_vector = signal_matrix[i,:]
			### if not in the insig_index vector, we use the index as one index set
			if not (index in insig_index):
				if index in index_signal_matrix_dict:
					index_signal_matrix_dict[ index ].append( signal_vector )
				else:
					uniq_index.append( index )
					index_signal_matrix_dict[ index ] = [ signal_vector ]
			### if in the insig_index vector, we use the insig_index ('X_X_X...') as one index set
			else:
				### replace index by X_...
				index_list = index.split('_')
				index = ''
				for i in range(0, len(index_list)-1):
					index = index + 'X_'
				index = index + 'X'
				### add to 
				if index in index_signal_matrix_dict:
					index_signal_matrix_dict[ index ].append( signal_vector )
				else:
					### keep order of index set
					uniq_index.append( index )
					index_signal_matrix_dict[ index ] = [ signal_vector ]
			### in last round, we append to index_vector_filter for original index vector relabeling
			if l == qda_round-1:
				index_vector_filter.append(index)
		print('calculating mean vector and cov matrix...')
		### initialize cov dict & mean dict
		index_signal_cov_dict = {}
		index_signal_mean_dict = {}
		index_set_signal_matrix_dict = {}
		### initialize mean matrix & counts matrix for index set output
		index_set_signal_mean_matrix = []
		index_set_peak_counts_matrix = []
		### loop uniq index
		for index in uniq_index:
			### extract index signal matrix
			one_index_matrix = np.array(index_signal_matrix_dict[ index ], dtype = float)
			#print(str(index)+': '+str(one_index_matrix.shape[0]))
			### calculating index matrix covariance matrix & mean vector
			one_index_matrix_cov = np.cov(one_index_matrix, rowvar = False)
			one_index_matrix_mean = np.mean(one_index_matrix, axis = 0)
			### append to covariance matrix dict & mean vector dict & index set matrix
			index_signal_cov_dict[ index ] = one_index_matrix_cov
			index_signal_mean_dict[ index ] = one_index_matrix_mean
			index_set_signal_matrix_dict[ index ] = one_index_matrix
			index_set_signal_mean_matrix.append(one_index_matrix_mean)
			index_set_peak_counts_matrix.append(one_index_matrix.shape[0])

		print('calculating Quadratic Scores...')
		if l==0:
			index_p_vec_index_set = np.empty((0, 1), float)
		score_i_exp_matrix = np.empty((signal_matrix.shape[0], 0), float)
		for index in uniq_index:
			### extract index signal matrix of one index set
			cov_i = index_signal_cov_dict[ index ]
			cov_i_inverse = np.linalg.inv(cov_i)
			mean_i = index_signal_mean_dict[ index ]
			signal_matrix_i = index_set_signal_matrix_dict[ index ]
			if l == 0:
				### calculate log scale score (for Just index set itself)
				d_index_set = np.sum(- 0.5 * np.log( abs(cov_i) ))
				score_i_index_set = d_index_set - 0.5 * np.sum( np.dot((signal_matrix_i-mean_i), cov_i_inverse) * (signal_matrix_i-mean_i), axis = 1 )
				score_i_exp_index_set = np.exp(score_i_index_set).reshape((score_i_index_set.shape[0],1))
				index_p_vec_index_set = np.concatenate((index_p_vec_index_set, score_i_exp_index_set), axis=0)
				#index_p_vec_index_set = np.array(index_p_vec_index_set).reshape(len(index_p_vec_index_set),1)
				### add od index set index name * nrow
				index_name_vec_i = np.array([[ index ]]*signal_matrix_i.shape[0])
				
			### calculate log scale score (for entire signal_matrix)
			d = np.sum(- 0.5 * np.log( abs(cov_i) ))
			score_i = d - 0.5 * np.sum( np.dot((signal_matrix-mean_i), cov_i_inverse) * (signal_matrix-mean_i), axis = 1 )
			### convert to exp scale
			score_i_exp = np.exp(score_i).reshape((score_i.shape[0],1))
			### cbind exp_scale score to score_i_exp_matrix
			score_i_exp_matrix = np.concatenate((score_i_exp_matrix, score_i_exp), axis=1)
		print('calculating Quadratic Scores...DONE')
		print('calculating max p and index...') 
		index_name_vec = []
		index_p_vec = []
		for p_vec in score_i_exp_matrix:
			index_i = uniq_index[np.argmax(p_vec)]
			index_name_vec.append(index_i)
			p_max = np.max(p_vec) / np.sum(p_vec)
			index_p_vec.append(p_max)
		print('calculating max p and index...DONE') 
		index_p_vec = np.array(index_p_vec).reshape(len(index_p_vec),1)
		index_name_vec = np.array(index_name_vec).reshape(len(index_name_vec),1)	
	### return all objects
	return { 'signal_matrix_bed': signal_matrix_bed, 'index_name_vec': index_name_vec, 'index_p_vec': index_p_vec, 'index_name_vec_index_set': index_name_vec_index_set, 'index_p_vec_index_set': index_p_vec_index_set, 'signal_matrix': signal_matrix, 'uniq_index': uniq_index, 'index_set_peak_counts_matrix': index_set_peak_counts_matrix, 'index_set_signal_mean_matrix': index_set_signal_mean_matrix, 'index_vector_filter': index_vector_filter }

################################################################################################
### vector most frequent element
def frequent(element_vector):
	most_freq_element = Counter(element_vector).most_common(1)[0][0]
	###### return output dict
	return most_freq_element

################################################################################################
### vector most enriched element
def ideas_enrich(indexset_vector, totalvector):
	### get ( index set & all ) ideas state counts
	indexset_counts_dict = Counter(indexset_vector)
	indexset_counts = []
	all_counts_dict = Counter(totalvector)
	all_counts = []
	ideas_state_list = []
	for records in indexset_counts_dict:
		indexset_counts.append(indexset_counts_dict[records])
		all_counts.append(all_counts_dict[records])
		### save the ideas state label order
		ideas_state_list.append(records)
	indexset_counts = np.array(indexset_counts, dtype = float)
	all_counts = np.array(all_counts, dtype = float)
	ideas_state_list = np.array(ideas_state_list)
	### calculate enrichment 
	smallnum = 100
	enrichment = (indexset_counts + smallnum) / (all_counts + smallnum)
	### get the most enriched ideas state
	most_enriched_element = ideas_state_list[np.argmax(enrichment)]
	return most_enriched_element

################################################################################################
### Shannon Entropy
def shannon_entropy(element_vector):
	### get element counts
	counts_dict = Counter(element_vector)
	counts = []
	for records in counts_dict:
		counts.append(counts_dict[records])
	counts = np.array(counts)
	### get entropy
	total = np.sum(counts)
	probability = counts / float(total)
	probability = probability[probability!=0]
	shannon_entropy = -np.sum(probability * np.log2(probability))
	###### return output dict
	return shannon_entropy

################################################################################################
### column based calculation
def matrix_col_cal(matrix, function, para=None):
	### get total column number
	colnum = matrix.shape[1]
	### for loop columns
	col_score_list = []
	for i in range(0, colnum):
		### extract column vector
		col_vec = matrix[:,i]
		### calculation
		if para is None:
			col_vec_score = function(col_vec)
		elif para.shape[1]==colnum:
			col_vec_score = function(col_vec, para[:,i])
		else:
			col_vec_score = function(col_vec, para)
		col_score_list.append(col_vec_score)
	col_score_list = np.array(col_score_list)
	###### return output dict
	return col_score_list

################################################################################################
### index set sth some score matrix
def index_set_score(index_name_vec, index_p_vec, sth_matrix_file, sth_start_col, uniq_index, method, smallnum, preorder, output_filename):
	### read sth matrix file
	sth_matrix = read2d_array(sth_matrix_file, 'str')
	### rm unused cols
	sth_matrix = sth_matrix[:, (sth_start_col-1):]
	### merge index id, index set p-value, sth matrix
	#print(index_name_vec.shape)
	#print(index_p_vec.shape)
	#print(sth_matrix.shape)
	sth_matrix_indexed = np.concatenate((index_name_vec, index_p_vec, sth_matrix), axis = 1)
	### get index set enriched score
	sth_enriched_dict = {}
	for records in sth_matrix_indexed:
		if not (records[0] in sth_enriched_dict):
			sth_enriched_dict[records[0]] = [ records[2:] ]
		else:
			sth_enriched_dict[records[0]].append( records[2:] )
	### get the index set score of sth
	index_set_sth_matrix = []
	index_set_list = []
	for index_set in uniq_index:
		### original index include all possible index (index number threshold will filter some index)
		if index_set in sth_enriched_dict:
			### keep the index set index order
			index_set_list.append(index_set)
			### column sum/mean
			if method == 'sum':
				matrix = np.array(sth_enriched_dict[ index_set ], dtype = float)
				matrix_col = np.sum(matrix, axis=0)
			elif method == 'mean':
				matrix = np.array(sth_enriched_dict[ index_set ], dtype = float)
				matrix_col = np.mean(matrix, axis=0)
			elif method == 'mostfreq':
				matrix = np.array(sth_enriched_dict[ index_set ], dtype = str)
				matrix_col = matrix_col_cal(matrix, frequent)
			elif method == 'sh':
				matrix = np.array(sth_enriched_dict[ index_set ], dtype = str)
				matrix_col = matrix_col_cal(matrix, shannon_entropy)
			elif method == 'mostenrich':
				matrix = np.array(sth_enriched_dict[ index_set ], dtype = str)
				matrix_col = matrix_col_cal(matrix, ideas_enrich, sth_matrix)
			elif method == 'enrich':
				matrix = np.array(sth_enriched_dict[ index_set ], dtype = float)
				### get expected counts by random chance
				sth_matrix = np.array(sth_matrix, dtype = float)
				col_total = np.sum(sth_matrix, axis = 0)
				pk_total = sth_matrix.shape[0]
				pk_index = matrix.shape[0]
				random_exp = col_total * pk_index / pk_total
				### get enrichment score
				matrix_col = (np.sum(matrix, axis=0) + smallnum) / (random_exp+ smallnum)
			index_set_sth_matrix.append( matrix_col )
	### convert to np array for concatenate (cbind)
	index_set_sth_matrix = np.array(index_set_sth_matrix)
	### cbind index_set index id & signal matrix
	index_set_list = np.array(index_set_list)
	index_set_list = index_set_list.reshape(index_set_list.shape[0], 1)
	index_set_sth_matrix = np.concatenate((index_set_list, index_set_sth_matrix), axis = 1)
	###### write output
	if preorder =='F':
		### indexed sth matrix
		write2d_array(sth_matrix_indexed, output_filename+'.indexed.txt')
		call('sort -k1,1 -k2,2rn ' + output_filename+'.indexed.txt' + ' > ' + output_filename+'.indexed.sort.txt', shell=True)
		call('rm ' + output_filename+'.indexed.txt', shell=True)
		### index set of sth 
		write2d_array(index_set_sth_matrix, output_filename+'.index_set.txt')
		call('sort -k1,1 ' + output_filename+'.index_set.txt' + ' > ' + output_filename+'.index_set.sort.txt', shell=True)
		call('rm ' + output_filename+'.index_set.txt', shell=True)

	elif preorder =='T':
		### use negative binomial to get new index
		m0 = sth_matrix_indexed[:,1:]
		m0 = np.array(m0, dtype=float)
		for i in range(0, 10):
			if i==0:
				m1 = m0
			else:
				m1 = m1[m1<=thresh]
			sigmean = np.mean(m1)
			sigvar = np.var(m1)
			sigprob = sigmean / sigvar
			### set threshold for sigprob
			if sigprob <= 0.1:
				### to avoid var too high relative to mean
				sigprob = 0.1
				sigsize = sigmean * sigprob / (1-sigprob)
				### get nb value
				thresh1 = nbinom.ppf(0.95, sigsize, sigprob)
			elif sigprob >= 1:
				### to avoid mean too high relative to var
				print('Should not use negative binomial!!!')
				if i ==0:
					thresh = poisson.ppf(0.95, sigmean)
				break

			if i >=1:
				if thresh1 == thresh:
					print('converged!!!')
					break
				else:
					thresh = thresh1
			else:
				thresh = thresh1
			print('round('+str(i)+'): '+str(thresh))

		print('for index set matrix...')
		### add new index to index_set matrix
		new_index_set_vector = []
		for info in index_set_sth_matrix:
			info_sig = info[1:]
			new_index = ''
			for i in info_sig:
				if float(i) > thresh:
					new_index = new_index+'1'
				else:
					new_index = new_index+'0'
			new_index_set_vector.append(new_index)
		new_index_set_vector = np.array(new_index_set_vector)
		print('sort by ni col & remove first ni col...')
		print(new_index_set_vector[np.argsort(new_index_set_vector)])
		new_index_set_order = np.argsort(new_index_set_vector)
		index_set_sth_matrix_ni_ordered = index_set_sth_matrix[new_index_set_order,:]
		### save order vector
		index_set_ni_sort = index_set_sth_matrix_ni_ordered[:,0].reshape(index_set_sth_matrix_ni_ordered.shape[0],1)
		### write output
		write2d_array(index_set_sth_matrix_ni_ordered, output_filename+'.index_set.sort.txt')
		write2d_array(index_set_ni_sort, output_filename+'.index_set_ni_sorted.txt')

		print('for index matrix...')
		### add new index to index matrix
		new_index_vector = []
		for info in sth_matrix_indexed:
			info_sig = info[1:]
			new_index = ''
			for i in info_sig:
				if float(i) > thresh:
					new_index = new_index+'1'
				else:
					new_index = new_index+'0'
			new_index_vector.append(new_index)
		new_index_vector = np.array(new_index_vector)
		print('sort by ni col & remove first ni col...')
		new_index_order = np.argsort(new_index_vector)
		sth_matrix_indexed_ni_ordered = sth_matrix_indexed[new_index_order,:]
		### get index set counts
		print('get index set counts...')
		unique, counts = np.unique(sth_matrix_indexed_ni_ordered[:,0], return_counts=True)
		unique_counts =  np.asarray((unique, counts)).T
		write2d_array(unique_counts, output_filename+'.index_set.count.txt')
		print(unique_counts)
		### save order vector
		index_ni_sort = sth_matrix_indexed_ni_ordered[:,0].reshape(sth_matrix_indexed_ni_ordered.shape[0],1)
		### write output
		write2d_array(index_ni_sort, output_filename+'.index.sort.txt')

	elif preorder =='kmeans':
		### heirarchical clustering
		n_cluster = index_set_sth_matrix.shape[0]
		#model = AgglomerativeClustering(linkage='average',n_clusters=n_cluster, compute_full_tree=True)
		model = KMeans(n_clusters=5, random_state=0)
		model.fit(index_set_sth_matrix[:,1:])
		hclust_label = model.labels_
		### hclust label sort
		index_set_ni_dict = {}
		for index, ni in zip(index_set_sth_matrix, hclust_label):
			index_set_ni_dict[index[0]] = ni
		### reorder index matrix and index set matrix
		print('reorder...')
		print('for index set matrix...')
		sth_matrix_index_set_ni_order = []
		for index_set_sth_matrix_indexed in index_set_sth_matrix:
			### read index
			index_set_tmp = index_set_sth_matrix_indexed[0]
			### index to hcluster labels
			sth_matrix_index_set_ni_order.append(index_set_ni_dict[index_set_tmp])
		### list to np.array
		sth_matrix_index_set_ni_order = np.array(sth_matrix_index_set_ni_order)
		print('sort by ni col & remove first ni col...')
		index_set_ni_order = np.argsort(sth_matrix_index_set_ni_order)
		index_set_sth_matrix_ni_ordered = index_set_sth_matrix[index_set_ni_order,:]
		### save order vector
		index_set_ni_sort = index_set_sth_matrix_ni_ordered[:,0].reshape(index_set_sth_matrix_ni_ordered.shape[0],1)
		### write output
		write2d_array(index_set_sth_matrix_ni_ordered, output_filename+'.index_set.sort.txt')
		write2d_array(index_set_ni_sort, output_filename+'.index_set_ni_sorted.txt')
		#call('sort -k1,1 ' + output_filename+'.index_set.txt ' + ' > ' + output_filename+'.index_set.sort.txt', shell=True)
		#call('rm ' + output_filename+'.index_set.txt', shell=True)

		print('for index matrix...')
		sth_matrix_index_ni_order = []
		for index in sth_matrix_indexed:
			### read index
			index_tmp = index[0]
			### index to hcluster labels
			sth_matrix_index_ni_order.append(index_set_ni_dict[index_tmp])
		### list to np.array
		sth_matrix_index_ni_order = np.array(sth_matrix_index_ni_order)
		print(sth_matrix_index_ni_order.shape)
		print('sort by ni col & remove first ni col...')
		index_ni_order = np.argsort(sth_matrix_index_ni_order)
		index_sth_matrix_ni_ordered = sth_matrix_indexed[index_ni_order,:]
		### write output
		write2d_array(index_sth_matrix_ni_ordered, output_filename+'.indexed.sort.txt')
		#call('sort -k1,1 -k3,3rn ' + output_filename+'.indexed.txt ' + ' > ' + output_filename+'.indexed.sort.txt', shell=True)
		#call('rm ' + output_filename+'.indexed.txt', shell=True)

	else:
		### read index order
		index_set_ni_ordered = read2d_array(preorder, str)
		### hclust order to dict
		index_set_ni_dict = {}
		for i in range(0, index_set_ni_ordered.shape[0]):
			index_set_tmp = index_set_ni_ordered[i,0]
			index_set_ni_dict[index_set_tmp] = i

		### reorder index matrix and index set matrix
		print('reorder...')
		print('for index set matrix...')
		sth_matrix_index_set_ni_order = []
		for index_set_sth_matrix_indexed in index_set_sth_matrix:
			### read index
			index_set_tmp = index_set_sth_matrix_indexed[0]
			### index to hcluster labels
			sth_matrix_index_set_ni_order.append(index_set_ni_dict[index_set_tmp])
		### list to np.array
		sth_matrix_index_set_ni_order = np.array(sth_matrix_index_set_ni_order)
		print('sort by ni col & remove first ni col...')
		index_set_ni_order = np.argsort(sth_matrix_index_set_ni_order)
		index_set_sth_matrix_ni_ordered = index_set_sth_matrix[index_set_ni_order,:]
		### save order vector
		index_set_ni_sort = index_set_sth_matrix_ni_ordered[:,0].reshape(index_set_sth_matrix_ni_ordered.shape[0],1)
		### write output
		write2d_array(index_set_sth_matrix_ni_ordered, output_filename+'.index_set.sort.txt')

		print('for index matrix...')
		sth_matrix_index_ni_order = []
		for index in sth_matrix_indexed:
			### read index
			index_tmp = index[0]
			### index to hcluster labels
			sth_matrix_index_ni_order.append(index_set_ni_dict[index_tmp])
		### list to np.array
		sth_matrix_index_ni_order = np.array(sth_matrix_index_ni_order)
		print('sort by ni col & remove first ni col...')
		index_ni_order = np.argsort(sth_matrix_index_ni_order)
		index_sth_matrix_ni_ordered = sth_matrix_indexed[index_ni_order,:]
		### write output
		write2d_array(index_sth_matrix_ni_ordered, output_filename+'.indexed.sort.txt')		


################################################################################################

################################################################################################
### get Binary index matrix
peak_bed = 'homerTable3.peaks.filtered.interval.bed'
peak_bed_colnum = 4
mark_list_index = 'peak_list.txt'
output_file_index = 'homerTable3.peaks.filtered.interval.bed.index.matrix.txt'
signal_col = 'N/A'
method = 'intersect'
sort_sigbed = 'T'
print('get binary matrix...')
#get_mark_matrix(peak_bed, peak_bed_colnum, mark_list_index, output_file_index, signal_col, method, sort_sigbed)

### get signal matrix
peak_bed = 'homerTable3.peaks.filtered.interval.bed'
peak_bed_colnum = 4
mark_list_signal = 'signal_list.txt'
output_file_signal = 'homerTable3.peaks.filtered.interval.bed.signal.matrix.txt'
signal_col = 5
method = 'map'
sort_sigbed = 'F'
print('get signal matrix...')
#get_mark_matrix(peak_bed, peak_bed_colnum, mark_list_signal, output_file_signal, signal_col, method, sort_sigbed)

### get ideas label matrix
peak_bed = 'homerTable3.peaks.filtered.interval.bed'
peak_bed_colnum = 0
mark_list_ideas = 'ideas_list.txt'
output_file_ideas = 'homerTable3.peaks.filtered.interval.bed.ideas.matrix.txt'
signal_col = 'N/A'
method = 'window'
sort_sigbed = 'F'
print('get ideas matrix...')
#get_mark_matrix(peak_bed, peak_bed_colnum, mark_list_ideas, output_file_ideas, signal_col, method, sort_sigbed)

### get TF ChIP-seq matrix
peak_bed = 'homerTable3.peaks.filtered.interval.bed'
peak_bed_colnum = 4
mark_list_chip = 'chip_list.txt'
output_file_chip = 'homerTable3.peaks.filtered.interval.bed.chip.matrix.txt'
signal_col = 'N/A'
method = 'intersect'
sort_sigbed = 'T'
print('get chip matrix...')
#get_mark_matrix(peak_bed, peak_bed_colnum, mark_list_chip, output_file_chip, signal_col, method, sort_sigbed)


### Multi-variable norm p-value (QDA)
index_matrix_start_col = 5
signal_matrix_start_col = 5
siglevel_counts = 0.95
small_value = 1
log_signal = 'T'
qda_round = 50
bins_folder = '/Volumes/MAC_Data/data/labs/zhang_lab/01projects/CD_viewer/bin/'
index_matrix = read2d_array('homerTable3.peaks.filtered.interval.bed.index.matrix.txt', 'str')
signal_matrix_od = read2d_array('homerTable3.peaks.filtered.interval.bed.signal.matrix.txt', 'str')
scale = 'F'

### use while loop to select the threshold of index set counts
index_matrix = read2d_array('homerTable3.peaks.filtered.interval.bed.index.matrix.txt', 'str')
insig_index_dict = select_index_set_counts_thresh(index_matrix, index_matrix_start_col, siglevel_counts)
index_vector = insig_index_dict['index_vector']
insig_index = insig_index_dict['insig_index']
index_count_thresh_2 = insig_index_dict['index_count_thresh']
index_vector_count_vec = insig_index_dict['index_vector_count_vec']

### calculating multiple variable norm density score
mvn_density_score_dict = mvn_density_score(signal_matrix_od, signal_matrix_start_col, log_signal, small_value, qda_round, index_vector, insig_index, scale)
signal_matrix_bed = mvn_density_score_dict['signal_matrix_bed']
index_name_vec = mvn_density_score_dict['index_name_vec']
index_p_vec = mvn_density_score_dict['index_p_vec']
index_name_vec_index_set = mvn_density_score_dict['index_name_vec_index_set']
index_p_vec_index_set = mvn_density_score_dict['index_p_vec_index_set']
signal_matrix = mvn_density_score_dict['signal_matrix']
uniq_index = mvn_density_score_dict['uniq_index']
index_set_peak_counts_matrix = mvn_density_score_dict['index_set_peak_counts_matrix']
index_set_signal_mean_matrix = mvn_density_score_dict['index_set_signal_mean_matrix']
index_vector_filter = mvn_density_score_dict['index_vector_filter']

print('check!!!check!!!check!!!check!!!check!!!check!!!')
print(signal_matrix_bed[0:10,:])
print('check!!!check!!!check!!!check!!!check!!!check!!!')
print(index_name_vec_index_set[0:10,:])
print(signal_matrix[0:10,:])
print('check!!!check!!!check!!!check!!!check!!!check!!!')
print(index_p_vec_index_set[0:10,:])
print('check!!!check!!!check!!!check!!!check!!!check!!!')

#######################

print('write signal mean matrix...')
output_file_signal_index_set = output_file_signal+'.index_set.txt'
#index_set_score(index_name_vec_index_set, index_p_vec_index_set, output_file_signal, 5, uniq_index, 'mean', 0,  'F', output_file_signal_index_set)

print('write binary sum matrix...')
output_file_index_index_set = output_file_index+'.index_set.txt'
#index_set_score(index_name_vec_index_set, index_p_vec_index_set, output_file_index, 5, uniq_index, 'sum', 0, 'F', output_file_index_index_set)

print('write pval mean matrix...')
output_file_pval = 'homerTable3.peaks.filtered.interval.bed.pval.matrix.txt'+'.index_set.txt'
p_matrix_index_set = np.concatenate((index_p_vec_index_set, index_p_vec_index_set), axis = 1)
write2d_array( p_matrix_index_set, output_file_pval)
output_file_pval_index_set = output_file_pval+'.index_set.txt'
#index_set_score(index_name_vec_index_set, index_p_vec_index_set, output_file_pval, 1, uniq_index, 'mean', 0, 'F', output_file_pval_index_set)

print('write mvn index matrix...')
output_file_index_mvn = 'homerTable3.peaks.filtered.interval.bed.mvn_index.matrix.txt'+'.index_set.txt'
index_mvn = []
for records in index_name_vec_index_set:
	tmp0 = records[0].split('_')
	### replace X by 1
	tmp = []
	for index in tmp0:
		if index == 'X':
			tmp.append('1')
		else:
			tmp.append(index)
	index_mvn.append(tmp)
index_mvn = np.array(index_mvn)
write2d_array(index_mvn, output_file_index_mvn)

print('write mvn binary sum matrix...')
output_file_index_mvn_index_set = output_file_index_mvn+'.index_set.txt'
#index_set_score(index_name_vec_index_set, index_p_vec_index_set, output_file_index_mvn, 1, uniq_index, 'sum', 0, 'F', output_file_index_mvn_index_set)

print('write chip enrichment matrix...')
smallnum = 50
output_file_chip_index_set = output_file_chip+'.index_set.txt'
#index_set_score(index_name_vec_index_set, index_p_vec_index_set, output_file_chip, 5, uniq_index, 'enrich', smallnum, 'F', output_file_chip_index_set)

print('write ideas most frequent label matrix...')
output_file_ideas_index_set_freq = output_file_ideas+'.freq'+'.index_set.txt'
#index_set_score(index_name_vec_index_set, index_p_vec_index_set, output_file_ideas, 5, uniq_index, 'mostfreq', 0, 'F', output_file_ideas_index_set_freq)
print('write ideas most enriched label matrix...')
output_file_ideas_index_set_enrich = output_file_ideas+'.enrich'+'.index_set.txt'
#index_set_score(index_name_vec_index_set, index_p_vec_index_set, output_file_ideas, 5, uniq_index, 'mostenrich', 0, 'F', output_file_ideas_index_set_enrich)
print('write ideas Shanno Entropy matrix...')
output_file_ideas_index_set_sh = output_file_ideas+'.sh'+'.index_set.txt'
#index_set_score(index_name_vec_index_set, index_p_vec_index_set, output_file_ideas, 5, uniq_index, 'sh', 0, 'F', output_file_ideas_index_set_sh)



#######################

mark_list_ideas = 'ideas_list.txt'
ideas_range_color_file = 'ideas_range_color.txt'
output_file_ideas_index_set = 'homerTable3.peaks.filtered.interval.bed.ideas.matrix.txt'
ideas_index_matrix_start_col = 3
ideas_index_set_matrix_start_col = 2
ideas_sh_index_set_matrix_start_col = 2

ideas_log2_transform = 'F'
ideas_log2_transform_add_smallnum = 0
ideas_index_boarder_color = 'NA'
ideas_index_set_boarder_color = 'gray'

print('use rect to plot ideas heatmap...')
#call('time Rscript ' + bins_folder + 'plot_rect.R ' + output_file_ideas_index_set_enrich+'.indexed.sort.txt' + ' ' + output_file_ideas_index_set_enrich+'.indexed.sort.txt' + '.png ' + mark_list_ideas + ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_matrix_start_col) + ' ' + ideas_index_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)
#call('time Rscript ' + bins_folder + 'plot_rect.R ' + output_file_ideas_index_set_enrich+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_enrich+'.index_set.sort.txt' + '.png ' + mark_list_ideas+ ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + ideas_index_set_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)
#call('time Rscript ' + bins_folder + 'plot_rect.R ' + output_file_ideas_index_set_freq+'.indexed.sort.txt' + ' ' + output_file_ideas_index_set_freq+'.indexed.sort.txt' + '.png ' + mark_list_ideas + ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_matrix_start_col) + ' ' + ideas_index_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)
#call('time Rscript ' + bins_folder + 'plot_rect.R ' + output_file_ideas_index_set_freq+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_freq+'.index_set.sort.txt' + '.png ' + mark_list_ideas+ ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + ideas_index_set_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)

print('use rect to plot ideas heatmap...with Shannon Entropy')
#call('time Rscript ' + bins_folder + 'plot_rect_filter.R ' + output_file_ideas_index_set_enrich+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_sh+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_enrich+'.sh.index_set.sort.txt' + '.png ' + mark_list_ideas+ ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + str(ideas_sh_index_set_matrix_start_col) + ' ' + ideas_index_set_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)
#call('time Rscript ' + bins_folder + 'plot_rect_filter.R ' + output_file_ideas_index_set_freq+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_sh+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_freq+'.sh.index_set.sort.txt' + '.png ' + mark_list_ideas+ ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + str(ideas_sh_index_set_matrix_start_col) + ' ' + ideas_index_set_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)

signal_high_color = 'red'
signal_low_color = 'white'
signal_log2_transform = 'F'
signal_log2_transform_add_smallnum = 0.001
signal_index_matrix_start_col = 3
signal_index_set_matrix_start_col = 2

print('use pheatmap to plot signal index & index set heatmap...')
#call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_signal_index_set+'.index_set.sort.txt' + ' ' + output_file_signal_index_set+'.index_set.sort.txt' + '.png ' + mark_list_signal + ' ' + str(signal_index_set_matrix_start_col) + ' ' + signal_high_color + ' ' + signal_low_color + ' ' + signal_log2_transform + ' ' + str(signal_log2_transform_add_smallnum), shell=True)
#call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_signal_index_set+'.indexed.sort.txt' + ' ' + output_file_signal_index_set+'.indexed.sort.txt' + '.png ' + mark_list_signal + ' ' + str(signal_index_matrix_start_col) + ' ' + signal_high_color + ' ' + signal_low_color + ' ' + signal_log2_transform + ' ' + str(signal_log2_transform_add_smallnum), shell=True)


chip_high_color = 'blue'
chip_low_color = 'white'
chip_log2_transform = 'F'
chip_log2_transform_add_smallnum = 0.001
chip_index_matrix_start_col = 3
chip_index_set_matrix_start_col = 2

print('use pheatmap to plot signal index & index set heatmap...')
#call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_chip_index_set+'.index_set.sort.txt' + ' ' + output_file_chip_index_set+'.index_set.sort.txt' + '.png ' + mark_list_chip + ' ' + str(chip_index_set_matrix_start_col) + ' ' + chip_high_color + ' ' + chip_low_color + ' ' + chip_log2_transform + ' ' + str(chip_log2_transform_add_smallnum), shell=True)
#call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_chip_index_set+'.indexed.sort.txt' + ' ' + output_file_chip_index_set+'.indexed.sort.txt' + '.png ' + mark_list_chip + ' ' + str(chip_index_matrix_start_col) + ' ' + chip_high_color + ' ' + chip_low_color + ' ' + chip_log2_transform + ' ' + str(chip_log2_transform_add_smallnum), shell=True)


index_high_color = 'black'
index_low_color = 'white'
index_log2_transform = 'T'
index_log2_transform_add_smallnum = 0.001
index_index_matrix_start_col = 3
index_index_set_matrix_start_col = 2

print('use pheatmap to plot signal index & index set heatmap...')
#call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_index_mvn_index_set+'.index_set.sort.txt' + ' ' + output_file_index_mvn_index_set+'.index_set.sort.txt' + '.png ' + mark_list_index + ' ' + str(index_index_set_matrix_start_col) + ' ' + index_high_color + ' ' + index_low_color + ' ' + index_log2_transform + ' ' + str(index_log2_transform_add_smallnum), shell=True)
#call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_index_mvn_index_set+'.indexed.sort.txt' + ' ' + output_file_index_mvn_index_set+'.indexed.sort.txt' + '.png ' + mark_list_index + ' ' + str(index_index_matrix_start_col) + ' ' + index_high_color + ' ' + index_low_color + ' ' + index_log2_transform + ' ' + str(index_log2_transform_add_smallnum), shell=True)


pval_high_color = 'orange'
pval_low_color = 'white'
pval_log2_transform = 'T'
pval_log2_transform_add_smallnum = 0.001
pval_index_matrix_start_col = 3
pval_index_set_matrix_start_col = 2

print('use pheatmap to plot p-value index & index set heatmap...')
#call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_pval_index_set+'.index_set.sort.txt' + ' ' + output_file_pval_index_set+'.index_set.sort.txt' + '.png ' + mark_list_index + ' ' + str(pval_index_set_matrix_start_col) + ' ' + pval_high_color + ' ' + pval_low_color + ' ' + pval_log2_transform + ' ' + str(pval_log2_transform_add_smallnum), shell=True)
#call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_pval_index_set+'.indexed.sort.txt' + ' ' + output_file_pval_index_set+'.indexed.sort.txt' + '.png ' + mark_list_index + ' ' + str(pval_index_matrix_start_col) + ' ' + pval_high_color + ' ' + pval_low_color + ' ' + pval_log2_transform + ' ' + str(pval_log2_transform_add_smallnum), shell=True)




#######################
#######################
#######################
#######################
#######################
#######################

print('write signal mean matrix...')
output_file_signal_index_set = output_file_signal
index_set_score(index_name_vec, index_p_vec, output_file_signal, 5, uniq_index, 'mean', 0, 'T', output_file_signal_index_set)

print('write binary sum matrix...')
output_file_index_index_set = output_file_index
index_set_score(index_name_vec, index_p_vec, output_file_index, 5, uniq_index, 'sum', 0, output_file_signal_index_set+'.index_set_ni_sorted.txt', output_file_index_index_set)

print('write pval mean matrix...')
output_file_pval = 'homerTable3.peaks.filtered.interval.bed.pval.matrix.txt'
p_matrix = np.concatenate((index_p_vec, index_p_vec), axis = 1)
write2d_array( p_matrix, output_file_pval)
output_file_pval_index_set = output_file_pval
index_set_score(index_name_vec, index_p_vec, output_file_pval, 1, uniq_index, 'mean', 0, output_file_signal_index_set+'.index_set_ni_sorted.txt', output_file_pval_index_set)

print('write chip enrichment matrix...')
smallnum = 50
output_file_chip_index_set = output_file_chip
index_set_score(index_name_vec, index_p_vec, output_file_chip, 5, uniq_index, 'enrich', smallnum, output_file_signal_index_set+'.index_set_ni_sorted.txt', output_file_chip_index_set)

print('write ideas most frequent label matrix...')
output_file_ideas_index_set_freq = output_file_ideas+'.freq'
index_set_score(index_name_vec, index_p_vec, output_file_ideas, 5, uniq_index, 'mostfreq', 0, output_file_signal_index_set+'.index_set_ni_sorted.txt', output_file_ideas_index_set_freq)
print('write ideas most enriched label matrix...')
output_file_ideas_index_set_enrich = output_file_ideas+'.enrich'
index_set_score(index_name_vec, index_p_vec, output_file_ideas, 5, uniq_index, 'mostenrich', 0, output_file_signal_index_set+'.index_set_ni_sorted.txt', output_file_ideas_index_set_enrich)
print('write ideas Shanno Entropy matrix...')
output_file_ideas_index_set_sh = output_file_ideas+'.sh'
index_set_score(index_name_vec, index_p_vec, output_file_ideas, 5, uniq_index, 'sh', 0, output_file_signal_index_set+'.index_set_ni_sorted.txt', output_file_ideas_index_set_sh)

print('write mvn index matrix...')
output_file_index_mvn = 'homerTable3.peaks.filtered.interval.bed.mvn_index.matrix.txt'
index_mvn = []
for records in index_name_vec:
	tmp0 = records[0].split('_')
	### replace X by 1
	tmp = []
	for index in tmp0:
		if index == 'X':
			tmp.append('1')
		else:
			tmp.append(index)
	index_mvn.append(tmp)
index_mvn = np.array(index_mvn)
write2d_array(index_mvn, output_file_index_mvn)

print('write mvn binary sum matrix...')
output_file_index_mvn_index_set = output_file_index_mvn
index_set_score(index_name_vec, index_p_vec, output_file_index_mvn, 1, uniq_index, 'sum', 0, output_file_signal_index_set+'.index_set_ni_sorted.txt', output_file_index_mvn_index_set)



#######################



signal_high_color = 'red'
signal_low_color = 'white'
signal_log2_transform = 'F'
signal_log2_transform_add_smallnum = 0.001
signal_index_matrix_start_col = 3
signal_index_set_matrix_start_col = 2

print('use pheatmap to plot signal index & index set heatmap...')
call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_signal_index_set+'.index_set.sort.txt' + ' ' + output_file_signal_index_set+'.index_set.sort.txt' + '.png ' + mark_list_signal + ' ' + str(signal_index_set_matrix_start_col) + ' ' + signal_high_color + ' ' + signal_low_color + ' ' + signal_log2_transform + ' ' + str(signal_log2_transform_add_smallnum), shell=True)
call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_signal_index_set+'.indexed.sort.txt' + ' ' + output_file_signal_index_set+'.indexed.sort.txt' + '.png ' + mark_list_signal + ' ' + str(signal_index_matrix_start_col) + ' ' + signal_high_color + ' ' + signal_low_color + ' ' + signal_log2_transform + ' ' + str(signal_log2_transform_add_smallnum), shell=True)


index_high_color = 'black'
index_low_color = 'white'
index_log2_transform = 'T'
index_log2_transform_add_smallnum = 0.001
index_index_matrix_start_col = 3
index_index_set_matrix_start_col = 2

print('use pheatmap to plot signal index & index set heatmap...')
call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_index_mvn_index_set+'.index_set.sort.txt' + ' ' + output_file_index_mvn_index_set+'.index_set.sort.txt' + '.png ' + mark_list_index + ' ' + str(index_index_set_matrix_start_col) + ' ' + index_high_color + ' ' + index_low_color + ' ' + index_log2_transform + ' ' + str(index_log2_transform_add_smallnum), shell=True)
call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_index_mvn_index_set+'.indexed.sort.txt' + ' ' + output_file_index_mvn_index_set+'.indexed.sort.txt' + '.png ' + mark_list_index + ' ' + str(index_index_matrix_start_col) + ' ' + index_high_color + ' ' + index_low_color + ' ' + index_log2_transform + ' ' + str(index_log2_transform_add_smallnum), shell=True)


mark_list_ideas = 'ideas_list.txt'
ideas_range_color_file = 'ideas_range_color.txt'
output_file_ideas_index_set = 'homerTable3.peaks.filtered.interval.bed.ideas.matrix.txt'
ideas_index_matrix_start_col = 3
ideas_index_set_matrix_start_col = 2
ideas_sh_index_set_matrix_start_col = 2

ideas_log2_transform = 'F'
ideas_log2_transform_add_smallnum = 0
ideas_index_boarder_color = 'NA'
ideas_index_set_boarder_color = 'gray'

print('use rect to plot ideas heatmap...')
call('time Rscript ' + bins_folder + 'plot_rect.R ' + output_file_ideas_index_set_enrich+'.indexed.sort.txt' + ' ' + output_file_ideas_index_set_enrich+'.indexed.sort.txt' + '.png ' + mark_list_ideas + ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_matrix_start_col) + ' ' + ideas_index_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)
call('time Rscript ' + bins_folder + 'plot_rect.R ' + output_file_ideas_index_set_enrich+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_enrich+'.index_set.sort.txt' + '.png ' + mark_list_ideas+ ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + ideas_index_set_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)
call('time Rscript ' + bins_folder + 'plot_rect.R ' + output_file_ideas_index_set_freq+'.indexed.sort.txt' + ' ' + output_file_ideas_index_set_freq+'.indexed.sort.txt' + '.png ' + mark_list_ideas + ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_matrix_start_col) + ' ' + ideas_index_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)
call('time Rscript ' + bins_folder + 'plot_rect.R ' + output_file_ideas_index_set_freq+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_freq+'.index_set.sort.txt' + '.png ' + mark_list_ideas+ ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + ideas_index_set_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)

print('use rect to plot ideas heatmap...with Shannon Entropy')
call('time Rscript ' + bins_folder + 'plot_rect_filter.R ' + output_file_ideas_index_set_enrich+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_sh+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_enrich+'.sh.index_set.sort.txt' + '.png ' + mark_list_ideas+ ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + str(ideas_sh_index_set_matrix_start_col) + ' ' + ideas_index_set_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)
call('time Rscript ' + bins_folder + 'plot_rect_filter.R ' + output_file_ideas_index_set_freq+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_sh+'.index_set.sort.txt' + ' ' + output_file_ideas_index_set_freq+'.sh.index_set.sort.txt' + '.png ' + mark_list_ideas+ ' ' + str(ideas_range_color_file) + ' ' + str(ideas_index_set_matrix_start_col) + ' ' + str(ideas_sh_index_set_matrix_start_col) + ' ' + ideas_index_set_boarder_color + ' ' + ideas_log2_transform + ' ' + str(ideas_log2_transform_add_smallnum), shell=True)



chip_high_color = 'blue'
chip_low_color = 'white'
chip_log2_transform = 'F'
chip_log2_transform_add_smallnum = 0.001
chip_index_matrix_start_col = 3
chip_index_set_matrix_start_col = 2

print('use pheatmap to plot signal index & index set heatmap...')
call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_chip_index_set+'.index_set.sort.txt' + ' ' + output_file_chip_index_set+'.index_set.sort.txt' + '.png ' + mark_list_chip + ' ' + str(chip_index_set_matrix_start_col) + ' ' + chip_high_color + ' ' + chip_low_color + ' ' + chip_log2_transform + ' ' + str(chip_log2_transform_add_smallnum), shell=True)
call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_chip_index_set+'.indexed.sort.txt' + ' ' + output_file_chip_index_set+'.indexed.sort.txt' + '.png ' + mark_list_chip + ' ' + str(chip_index_matrix_start_col) + ' ' + chip_high_color + ' ' + chip_low_color + ' ' + chip_log2_transform + ' ' + str(chip_log2_transform_add_smallnum), shell=True)



pval_high_color = 'orange'
pval_low_color = 'white'
pval_log2_transform = 'T'
pval_log2_transform_add_smallnum = 0.001
pval_index_matrix_start_col = 3
pval_index_set_matrix_start_col = 2

print('use pheatmap to plot p-value index & index set heatmap...')
call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_pval_index_set+'.index_set.sort.txt' + ' ' + output_file_pval_index_set+'.index_set.sort.txt' + '.png ' + mark_list_index + ' ' + str(pval_index_set_matrix_start_col) + ' ' + pval_high_color + ' ' + pval_low_color + ' ' + pval_log2_transform + ' ' + str(pval_log2_transform_add_smallnum), shell=True)
call('time Rscript ' + bins_folder + 'plot_pheatmap.R ' + output_file_pval_index_set+'.indexed.sort.txt' + ' ' + output_file_pval_index_set+'.indexed.sort.txt' + '.png ' + mark_list_index + ' ' + str(pval_index_matrix_start_col) + ' ' + pval_high_color + ' ' + pval_low_color + ' ' + pval_log2_transform + ' ' + str(pval_log2_transform_add_smallnum), shell=True)



