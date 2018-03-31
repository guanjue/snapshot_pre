####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
output_filename = args[2]
signal_input_list = args[3]
signal_matrix_start_col = args[4]
high_color = args[5]
low_color = args[6]
log2 = args[7]
smallnum = as.numeric(args[8])

####################################################
### use pheatmap to plot heatmaps
color_heatmap = function(color_matrix, high_color, low_color, format, outputname){
	library(pheatmap)
	format(outputname, width = 1000, height = 1000) ### output name
	par(mar=c(5,0.5,0.5,0.5)) ### set heatmap margins
	### plot pheatmap
	my_colorbar=colorRampPalette(c(low_color, high_color))(n = 128)
	col_breaks = c(seq(0, 2000,length=33))
	pheatmap(color_matrix, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE)
	dev.off()
}

#################################################### 
############ read input files
####################################################
### read signal matrix file
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE))
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
### convert to numeric matrix
class(signal_matrix) = 'numeric'

### log2 transform
if (log2=='T'){
	signal_matrix = log2(signal_matrix+smallnum)
}

###### read colnames file
colname_file = read.table(signal_input_list, header=F)
### add colnames
if (dim(signal_matrix)[2]==length(colname_file[,2])){
	colnames(signal_matrix) = colname_file[,2]
}

format = png
color_heatmap(signal_matrix, high_color, low_color, format, output_filename)


