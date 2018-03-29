####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
counts_table = args[1]
output_filename = args[2]
high_color = args[3]
low_color = args[4]
log2 = args[5]
smallnum = as.numeric(args[6])

####################################################
### use pheatmap to plot heatmaps
color_heatmap = function(color_matrix, high_color, low_color, format, outputname){
	library(pheatmap)
	format(outputname, width = 1000, height = 1000) ### output name
	par() ### set heatmap margins
	### plot pheatmap
	my_colorbar=colorRampPalette(c(low_color, high_color))(n = 128)
	col_breaks = c(seq(0, 2000,length=33))
	pheatmap(color_matrix, color=my_colorbar, cellnote=color_matrix, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE)
	dev.off()
}

#################################################### 
############ read input files
####################################################
### read signal matrix file
counts_table_od = as.matrix(read.table(counts_table, header=FALSE))
### extract signal matrix without info
counts_table = cbind(counts_table_od[,2], counts_table_od[,2])
### convert to numeric matrix
class(counts_table) = 'numeric'

### log2 transform
if (log2=='T'){
	counts_table = log2(counts_table+smallnum)
}


format = png
color_heatmap(counts_table, high_color, low_color, format, output_filename)


