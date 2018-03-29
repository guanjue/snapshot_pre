####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
binary_matrix_file = args[1]
output_file = args[2]
uselog2 = args[3]
smallnum = as.numeric(args[4])
limline = as.numeric(args[5])

### read file
binary_matrix = as.matrix(read.table(binary_matrix_file, header=FALSE))
index_set_table = table(binary_matrix)
index_set_num = index_set_table

index_set_table_thresh = index_set_table[index_set_table>=limline]

if (uselog2=='T'){
	index_set_num = log2(index_set_num+smallnum)
}

png(output_file)
plot(density(index_set_num))
if (uselog2=='T'){
	abline(v=log2(limline), col='red')
} else {
	abline(v=(limline), col='red')
}
dev.off()



write.table(index_set_table_thresh, paste(output_file, '.binary.counts.thresh.bed', sep=''), sep='\t', quote = FALSE, row.names = TRUE, col.names = FALSE)
