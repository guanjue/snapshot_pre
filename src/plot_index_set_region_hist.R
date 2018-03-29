### get parameters
args = commandArgs(trailingOnly=TRUE)
index_set_inputfile = args[1]
index_set_outfile = args[2]

### read index set matrix
read_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	rownames(data_index_set) = data_index_set[,1]
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(data_index_set)
}

###########
data_index_set = read_matrix(index_set_inputfile)

DNA_region_num = apply(data_index_set, 1, max) 

hist_table=table(DNA_region_num)

print(min(hist_table))
print(mean(hist_table))
print(median(hist_table))
print(max(hist_table))

h=hist(DNA_region_num,breaks=1000,plot=F)
#print((h$counts))
h$counts=log10(h$counts+1)
#print((h$counts))
pdf(index_set_outfile)
plot(h)
dev.off()

# Rscript plot_index_set_region_hist.R celltype.index_set.sorted.txt

