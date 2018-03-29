### get parameters
args = commandArgs(trailingOnly=TRUE)
index_set_counts_inputfile = args[1]
index_set_counts_hist_outfile = args[2]

###########
data_index_set_counts = read.table(index_set_counts_inputfile)
data_index_set_counts_vec = data_index_set_counts[1:dim(data_index_set_counts)[1]-1,]
print(summary(data_index_set_counts))

h=hist(data_index_set_counts_vec,breaks=1000,plot=F)
#print((h$counts))
h$counts=log10(h$counts+1)
#print((h$counts))
png(index_set_counts_hist_outfile)
plot(h)
abline(v=data_index_set_counts[dim(data_index_set_counts)[1],1], col = 'blue', lty=2)
dev.off()

