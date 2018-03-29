library(pheatmap)
### get parameters
args = commandArgs(trailingOnly=TRUE)
index_set_inputfile = args[1]
index_inputfile = args[2]

filter_index_set_output = args[3]

index_set_all_heatmap = args[4]
index_set_thresh_heatmap = args[5]
index_all_heatmap = args[6]

threshold = as.numeric(args[7])
color = args[8]

### set heatmap colors
my_colorbar=colorRampPalette(c('white',color))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

### read index set matrix
read_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

###########
data_index_set_od = read_matrix(index_set_inputfile)
### log2 scale
data_index_set = log2(data_index_set_od+1)
### plot index set heatmap
pheatmap(data_index_set, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_all_heatmap)
### plot index set heatmap (filter by threshold)
data_index_set_filtered = data_index_set[apply(data_index_set, 1, max) > log2(threshold), ]
pheatmap(data_index_set_filtered, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_thresh_heatmap)

###### write filtered index_set matrix
data_index_set_filtered_od = data_index_set_od[apply(data_index_set, 1, max) > log2(threshold), ]
write.table(data_index_set_filtered_od, file=filter_index_set_output, quote=FALSE, sep='\t', row.names = TRUE, col.names = NA)

###########
### read index matrix
data_index = read_matrix(index_inputfile)
### plot index heatmap
pheatmap(data_index, color=my_colorbar, cluster_cols = FALSE, cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_all_heatmap)
###########

# Rscript plot_index_set_module.R celltype.index_set.sorted.txt celltype.binary_pattern.sorted.txt black 200 index_set_all.pdf index_set_thresh.pdf index.png