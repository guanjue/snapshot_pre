args = commandArgs(trailingOnly=TRUE)
color_table = args[1]
outputname = args[2]


color_heatmap = function(color_matrix, outputname, format, border_color){
	format(outputname, width = 1000, height = 1000) ### output name
	par(mar=c(5,0.5,0.5,0.5)) ### set heatmap margins
	colbin_len = 10 ### column bin size
	rowbin_len = 10 ### row bin size
	### row reverse
	color_matrix = color_matrix[nrow(color_matrix):1,]
	### plot areas
	plot(c(0, dim(color_matrix)[2]*colbin_len), c(0, dim(color_matrix)[1]*rowbin_len), xaxt = "n", yaxt = "n", xaxs="i", yaxs="i", type = "n", xlab = "", ylab = "",main = "")
	### add color matrix colname as heatmap colname
	axis(1, c(1 : dim(color_matrix)[2])*colbin_len-0.5*colbin_len, colnames(color_matrix), las = 2, col.axis = "black", tick=FALSE)
	### use for loop to add rectangle with different color
	for (coln in c(1 : dim(color_matrix)[2])){ ### loop columns
		for (rown in c(1 : dim(color_matrix)[1])){ ### loop rows
			### add rectangle
			rect( (coln-1)*colbin_len, (rown-1)*rowbin_len, coln*colbin_len, rown*rowbin_len, col = color_matrix[rown, coln], border=border_color, lwd = 0 )
		}
	}
	dev.off()
}


#################################################### 
############ read input files
#################################################### 
###### read signal matrix 
heatmap_save_type = png
heatmap_boarder_col = 'black'
### read signal matrix file
signal_matrix_od = as.matrix(read.table(color_table, header=FALSE))
### extract signal matrix without info
signal_matrix = signal_matrix_od[,2]
print(signal_matrix)
signal_matrix_color = NULL
for (rgb_num in signal_matrix){
	color_vector_tmp = unlist(strsplit(rgb_num, ','))
	class(color_vector_tmp) = 'numeric'
	color_vector_tmp = color_vector_tmp / 255
	color_tmp = rgb(color_vector_tmp[1], color_vector_tmp[2], color_vector_tmp[3])
	signal_matrix_color = rbind(signal_matrix_color,c(color_tmp, color_tmp))
}
rownames(signal_matrix_color) = signal_matrix_od[,1]
####################################################
###### plot heatmap
color_heatmap(signal_matrix_color, outputname, heatmap_save_type, heatmap_boarder_col)

# Rscript /Volumes/MAC_Data/data/labs/zhang_lab/01projects/CD_viewer/bin/plot_color_table.R cell_merged_state17.signal.matrix.txt.sort.rgb.bed.color.table.txt cell_merged_state17.color.png

