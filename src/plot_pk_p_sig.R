library(pheatmap)
### get parameters
args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]

p_col = as.numeric(args[2])
p_color = args[3]
sig_col_start = as.numeric(args[4])
sig_color = args[5]

signal_list = args[6]
label_list_col = as.numeric(args[7])

outputname = args[8]

####################################################
data_index_set_od = as.matrix(read.table(inputfile, header=FALSE))
data_index_sig_matrix = data_index_set_od[,sig_col_start:dim(data_index_set_od)[2]]
data_index_p_matrix = cbind(data_index_set_od[,p_col], data_index_set_od[,p_col])

### add colnames
signal_list_matrix = as.matrix(read.table(signal_list, header=FALSE))[,2]
colnames(data_index_sig_matrix) = signal_list_matrix
### for p-value matrix use the longest string as colnames
ct_len = lapply(signal_list_matrix, function(x) length(unlist(strsplit(x, split=''))) )
colnames(data_index_p_matrix) = c(signal_list_matrix[which.max(ct_len)], signal_list_matrix[which.max(ct_len)])

### convert to numeric matrix
class(data_index_sig_matrix) = "numeric" 
class(data_index_p_matrix) = "numeric" 

### plot pheatmap
print('plot p-value heatmap...')
my_colorbar=colorRampPalette(c('white', p_color))(n = 128)
col_breaks = c(seq(0, 2000,length=33))
png(paste(outputname,'.p.png', sep=''), width = 1000, height = 1000)
pheatmap(data_index_p_matrix, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE)
dev.off()

print('plot signal matrix heatmap...')
my_colorbar=colorRampPalette(c('white', sig_color))(n = 128)
col_breaks = c(seq(0, 2000,length=33))
png(paste(outputname,'.sig.png', sep=''), width = 1000, height = 1000)
pheatmap(data_index_sig_matrix, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE)
dev.off()