
library(dplyr)
library(stringr)

#read-in gtf file
gtf = read.table("/Users/ryan.gorzek/Downloads/Monodelphis_domestica.ASM229v1.112.filtered.gtf", skip = 5, 
                 sep = "\t", header = F)
gtf$gene = apply(gtf, 1, function(x) strsplit(strsplit(x[9], ";")[[1]][1], " ")[[1]][2]) #make column for gene name
# gtf.lnc = gtf[str_detect(gtf$V9, "lncRNA") & gtf$V3 == "gene",]
# gtf.lnc$V3 <- "CDS"
# gtf = rbind(gtf, gtf.lnc)

gtf.gene = gtf[gtf$V3 == "gene",] #subset entries called 'gene'
gtf.gene = gtf.gene[,c("V1", "V3", "gene", "V4", "V5", "V7")]
gtf.gene = gtf.gene[order(gtf.gene$V1, gtf.gene$V4),] #order by chromosome and start

gtf.bed = gtf.gene[,c("V1", "V4", "V5", "gene")]
gtf.bed$score = ""
gtf.bed$strand = gtf.gene$V7
write.table(gtf.bed, file = "/Users/ryan.gorzek/Downloads/Monodelphis_domestica.ASM229v1.112.filtered.bed", quote = F, sep = "\t", row.names = F, col.names = F)
