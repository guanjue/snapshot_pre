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

####################################################
color_heatmap = function(color_matrix, outputname, format, border_color){
	format(outputname) ### output name
	par(mar=c(5,0.5,0.5,0.5)) ### set heatmap margins
	colbin_len = 10 ### column bin size
	rowbin_len = 10 ### row bin size
	### plot areas
	plot(c(0, dim(color_matrix)[2]*colbin_len), c(0, dim(color_matrix)[1]*rowbin_len), xaxt = "n", yaxt = "n", xaxs="i", yaxs="i", type = "n", xlab = "", ylab = "",main = "")
	### add color matrix colname as heatmap colname
	axis(1, c(1 : dim(color_matrix)[2])*colbin_len-0.5*colbin_len, colnames(color_matrix), las = 2, col.axis = "black", tick=FALSE)
	### use for loop to add rectangle with different color
	for (coln in c(1 : dim(color_matrix)[2])){ ### loop columns
		for (rown in c(1 : dim(color_matrix)[1])){ ### loop rows
			### add rectangle
			rect( (coln-1)*colbin_len, (rown-1)*rowbin_len, coln*colbin_len, rown*rowbin_len, col = color_matrix[rown,coln], border=border_color, lwd = 0 )
		}
	}
	dev.off()
}
####################################################

### set heatmap colors
color_info = read.table(color, header=TRUE)
### extract rgb colors
color_vector = t(apply(color_info, 1, FUN=function(x) as.numeric(unlist(strsplit(x[4], ','))) ))

### read index set matrix
read_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=TRUE))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

###########
data_index_set_od = read_matrix(index_set_inputfile)
### log2 scale
data_index_set = log2(data_index_set_od+1)

max(data_index_set)
### plot index set heatmap
#color_heatmap(data_index_set_colormatrix, index_set_all_heatmap)
#pheatmap(data_index_set, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_all_heatmap)
### plot index set heatmap (filter by threshold)
data_index_set_filtered = data_index_set[apply(data_index_set, 1, max) > log2(threshold), ]


#pheatmap(data_index_set_filtered, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE, filename = index_set_thresh_heatmap)

data_index_set_filtered_col_matrix = data_index_set_filtered
#print(head(data_index_set_filtered))
print(color_vector[1,])
data_index_set_filtered_col_matrix[data_index_set_filtered==0] = rgb(color_vector[1,1],color_vector[1,2],color_vector[1,3])
data_index_set_filtered_col_matrix[data_index_set_filtered!=0] = rgb(color_vector[2,1],color_vector[2,2],color_vector[2,3])

color_heatmap(data_index_set_filtered_col_matrix, paste(index_set_thresh_heatmap,'.pdf',sep=''), pdf, 'gray')

###### write filtered index_set matrix
data_index_set_filtered_od = data_index_set_od[apply(data_index_set, 1, max) > log2(threshold), ]
write.table(data_index_set_filtered_od, file=filter_index_set_output, quote=FALSE, sep='\t', row.names = TRUE, col.names = NA)

###########
### read index matrix
data_index = read_matrix(index_inputfile)

data_index_colormatrix = data_index
for (i in c(1:dim(color_info)[1])){
	data_index_colormatrix[data_index == color_info[i,1]] = rgb(color_vector[i,1],color_vector[i,2],color_vector[i,3])
}

### plot index heatmap
color_heatmap(data_index_colormatrix, index_all_heatmap, png, 'gray')
###########





# Rscript plot_index_set_module.R celltype.index_set.sorted.txt celltype.binary_pattern.sorted.txt black 200 index_set_all.pdf index_set_thresh.pdf index.png


