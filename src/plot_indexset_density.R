data = read.table('atac_20cell.bed.index.matrix.txt', header = F)

data_label = data[,5:dim(data)[2]]

data_index = apply(data_label, 1, function(x) paste(x, collapse='_'))

data_index_table = table(data_index)

smallnum = 0

limit = 200

uselog2 = 'T'

if (uselog2=='T'){
	sig = log2(data_index_table+smallnum)
}

png('index_set_densit.png')
plot(density(sig))
abline(v=log2(limit), col='red')
dev.off()



