### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix = args[1]
count_threshold = as.numeric(args[2])


index = read.table(index_matrix, header = F)

index = apply(index, 1, function(x) paste(x[5:dim(index)[2]], collapse='_'))
#print(head(index))
index_count = table(index)
#print(head(index_count))

png('density_index_count.png', width = 8, height = 8, units = 'in', res = 300)
plot(density(log2(index_count)), main='density plot of index counts (log2 scale)')
abline(v=log2(count_threshold), col='red', lwd=1.5, lty=2)
dev.off()

