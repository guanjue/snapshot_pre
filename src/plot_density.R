####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
bed_file = args[1]
output_file = args[2]
uselog2 = args[3]
smallnum = as.numeric(args[4])

### read file
bed = as.matrix(read.table(bed_file, header=FALSE))
sig = as.numeric(bed[,4])

if (uselog2=='T'){
	sig = log2(sig+smallnum)
}

png(output_file)
plot(density(sig))
abline(v=log2(3), col='red')
dev.off()