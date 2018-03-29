library(networkD3, lib.loc='/storage/home/gzx103/work/r_package/')
#library(data.tree)
library(igraph, lib.loc='/storage/home/gzx103/work/r_package/')

####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
cd_tree = args[2]
signal_input_list = args[3]
signal_matrix_start_col = args[4]
high_color = args[5]
low_color = args[6]
log2 = args[7]
smallnum = as.numeric(args[8])

#################################################### 
############ read input files
####################################################
### read signal matrix file
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE))
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
### get index_set name
index_set_name = signal_matrix_od[,1]
### convert to numeric matrix
class(signal_matrix) = 'numeric'

###### read colnames file
colname_file = read.table(signal_input_list, header=F)
colname = colname_file[,2]

### read cell development tree file
tree = read.table(cd_tree, header = F, sep=',')
tree.df = as.data.frame(tree)
colnames(tree.df) = c('Node.1', 'Node.2')

source_node = unique(tree.df[,1])
target_node = unique(tree.df[,2])
all_source_node = setdiff(source_node, target_node)
all_target_node = setdiff(target_node, source_node)

### get color list
color_number = 300
color_range = max(signal_matrix) - min(signal_matrix)

my_colorbar=colorRampPalette(c(low_color, high_color))(n = color_number+1)

### plot trees
for (i in seq(1,dim(signal_matrix)[1])){
	### convert mean signal vector to color vector
	values = signal_matrix[i,]
	### get color list
	values_id = as.integer( (values - min(signal_matrix))/color_range * color_number )
	### value to color (10X)
	value_col = my_colorbar[values_id]

	### get tree
	tree.igraph = graph.data.frame(tree.df, directed=TRUE)
	tree_names = V(tree.igraph)$name
	### sort colnames by tree nodes id
	match_id = match(tree_names, colname)
	V(tree.igraph)$color = value_col[match_id]
	V(tree.igraph)$size = 25

	### correlation score
	path_sig_cor_vector = c()
	path_lw_cor = 0
	dff_sum_all = 0
	j=1
	for (s in all_source_node){
		for (t in all_target_node){
			path = unlist(all_simple_paths(tree.igraph, s, t))
			print(path)
			path_names = names(path)
			print(path_names)
			print(colname)
			### sort colnames by tree nodes id
			match_id = match(path_names, colname)
			print(match_id)
			path_signal = signal_matrix[i,][match_id]
			path_signal = path_signal[!is.na(path_signal)]
			print(path_signal)
			if (length(path_signal)>2){
				dff_sum_all = dff_sum_all + diff(path_signal)
				path_sig_cor = cor(path_signal, seq(1,length(path_signal)), method = 'pearson')
				path_sig_cor_vector[j] = length(path_signal) * path_sig_cor
				path_lw_cor = path_lw_cor + length(path_signal) * abs(path_sig_cor)
				print(path_sig_cor)
				j = j+1		
			}
		}
	}
	print(path_sig_cor_vector)
	max_cor = max(abs(path_sig_cor_vector))
	print('max path correlation')
	print(max_cor)

	print('length weighted path correlation')
	print(path_lw_cor)

	print('diff sum')
	print(dff_sum_all)

	### plot tree
	png(paste(toString(as.integer(dff_sum_all*1000)), '.', toString(i-1), '.', signal_input_list, index_set_name[i], '.tree.png', sep = ''), width = 1200, height = 1200)
	plot(tree.igraph, layout = layout_as_tree(tree.igraph, root=c(1)))
	dev.off()
}

