
library(dplyr)
library(stringr)

#read-in gtf file
gtf = read.table("/Users/ryan.gorzek/Downloads/Monodelphis_domestica.ASM229v1.112.filtered.gtf", skip = 5, 
                 sep = "\t", header = F)
gtf$gene = apply(gtf, 1, function(x) strsplit(strsplit(x[9], ";")[[1]][1], " ")[[1]][2]) #make column for gene name
gtf.lnc = gtf[str_detect(gtf$V9, "lncRNA") & gtf$V3 == "gene",]
gtf.lnc$V3 <- "CDS"
gtf = rbind(gtf, gtf.lnc)

gtf.gene = gtf[gtf$V3 == "gene",] #subset entries called 'gene'
gtf.gene = gtf.gene[,c("V1", "V3", "gene", "V4", "V5", "V7")]
gtf.gene = gtf.gene[order(gtf.gene$V1, gtf.gene$V4),] #order by chromosome and start

#find genes overlapping with upstream gene
x = c()
for (i in c(2:nrow(gtf.gene)))
{
  print(i)
  x = c(x, ifelse(gtf.gene[i,]$V4 > gtf.gene[i-1,]$V4 && 
                  gtf.gene[i,]$V4 < gtf.gene[i-1,]$V5 && 
                  gtf.gene[i,]$V7 == gtf.gene[i-1,]$V7, 
                  "T", "F"))
}
x = c("F", x)
gtf.gene$overlapLeft = x
rm(x)

#find genes overlapping with downstream gene
x = c()
for (i in c(1:(nrow(gtf.gene)-1)))
{
  print(i)
  x = c(x, ifelse(gtf.gene[i,]$V5 > gtf.gene[i+1,]$V4 && 
                  gtf.gene[i,]$V5 < gtf.gene[i+1,]$V5 && 
                  gtf.gene[i,]$V7 == gtf.gene[i-1,]$V7, 
                  "T", "F"))
}
x = c(x, "F")
gtf.gene$overlapRight = x
rm(x)

gtf.CDS = gtf[gtf$V3 == "CDS",] #subset entries called 'CDS'
#modify CDS file to make one entry per gene with smallest starts & largest stops
gtf.CDS = gtf.CDS %>% group_by(V3, gene) %>% summarize(V4 = min(V4), V5 = max(V5)) 

merge = merge(gtf.gene, gtf.CDS, by.x = "gene", by.y = "gene") #merge gene and CDS dataframes
#reduce to 'normal chromosomes'
#rest = gtf[(gtf$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "X", "MT")) == FALSE,]
#merge = merge[merge$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "X", "MT"),]

#nrow(merge[merge$V4.y - merge$V4.x == 0,]) #number of genes with CDS start = gene start
#nrow(merge[merge$V5.y - merge$V5.x == 0,]) #number of genes with CDS stop = gene stop

gtf.bed = gtf.gene[gtf.gene$gene %in% merge$gene, c("V1", "V4", "V5", "gene")]
gtf.bed$score = ""
gtf.bed$strand = "+" # gtf.gene$V7
write.table(gtf.bed, file = "/Users/ryan.gorzek/Downloads/opossum_gene.bed", quote = F, sep = "\t", row.names = F, col.names = F)

# #read in left and right distance to closest gene output from bedtools
temp.lr = read.table("/Users/ryan.gorzek/Downloads/dists.bed", sep = "\t", header = F)
# temp = read.table("/users/saumyajain/Downloads/bed_dist_left.txt", sep = "\t", header = F)
temp = temp.lr[temp.lr$V13 < 0,]
temp = temp[, c("V4", "V13")]
temp = temp %>% group_by(V4) %>% summarise(V13 = min(V13))
colnames(temp) = c("gene", "dist.Left")
merge = merge(merge, temp, by.x = "gene", by.y = "gene")
merge$dist.Left = merge$dist.Left*-1
# temp = read.table("/users/saumyajain/Downloads/bed_dist_right.txt", sep = "\t", header = F)
temp = temp.lr[temp.lr$V13 > 0,]
temp = temp[, c("V4", "V13")]
temp = temp %>% group_by(V4) %>% summarise(V13 = min(V13))
colnames(temp) = c("gene", "dist.Right")
merge = merge(merge, temp, by.x = "gene", by.y = "gene")
rm(temp)
merge$dist.Left = ifelse(merge$overlapLeft == "Y", 0, merge$dist.Left)
merge$dist.Right = ifelse(merge$overlapRight == "Y", 0, merge$dist.Right)

gtf$modify.Left = ifelse(gtf$gene %in% merge[abs(merge$V4.x - merge$V4.y) < 10000,]$gene, "Y", "N")
gtf$modify.Right = ifelse(gtf$gene %in% merge[abs(merge$V5.x - merge$V5.y) < 10000,]$gene, "Y", "N")

gtf.modified = gtf[gtf$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "X", "MT"),] %>% 
               group_by(gene) %>% 
               mutate(minV4 = min(V4), maxV5 = max(V5))
temp = merge[, c("gene", "dist.Left", "dist.Right")]

gtf.modified = merge(gtf.modified, temp, by.x = "gene", by.y = "gene", all.x = T)
gtf.modified$dist.Left = as.numeric(gtf.modified$dist.Left)
gtf.modified$dist.Right = as.numeric(gtf.modified$dist.Right)

gtf.modified$dist.Left = apply(gtf.modified, 1, function(x) min(5000, abs(as.numeric(x[15]))))
gtf.modified$dist.Right = apply(gtf.modified, 1, function(x) min(5000, abs(as.numeric(x[16]))))
gtf.modified[is.na(gtf.modified)] = 0
gtf.modified$modify.Left = ifelse((gtf.modified$V7 == "-" & gtf.modified$modify.Left == "Y"), "Y", "N")
gtf.modified$modify.Right = ifelse((gtf.modified$V7 == "+" & gtf.modified$modify.Right == "Y"), "Y", "N")

x.utr = filter(gtf.modified, modify.Left == "Y", V3 == "transcript", V4 == minV4)
x.utr = rbind(x.utr, filter(gtf.modified, modify.Right == "Y", V3 == "transcript", V5 == maxV5))
x.utr$V3 = "three_prime_utr"
x.utr$V4 = ifelse(x.utr$V7 == "+", x.utr$V5+1, x.utr$V4)
x.utr$V5 = ifelse(x.utr$V7 == "-", x.utr$V4-1, x.utr$V5)
x.utr$V5 = floor(ifelse(x.utr$V7 == "+", x.utr$V5+x.utr$dist.Right, x.utr$V5))
x.utr$V4 = floor(ifelse(x.utr$V7 == "-", x.utr$V4 - x.utr$dist.Left, x.utr$V4))

x.exon <- x.utr
x.exon$V3 <- "exon"

gtf.modified$modify.Left = ifelse(gtf.modified$V3 %in% c("gene", "transcript") & 
                                    gtf.modified$modify.Left == "Y" & 
                                    gtf.modified$V4 == gtf.modified$minV4, 
                                  "Y", "N")
gtf.modified$modify.Right = ifelse(gtf.modified$V3 %in% c("gene", "transcript") & 
                                    gtf.modified$modify.Right == "Y" & 
                                    gtf.modified$V5 == gtf.modified$maxV5, 
                                  "Y", "N")

gtf.modified$V4 = floor(ifelse(gtf.modified$modify.Left == "Y", 
                               gtf.modified$V4 - gtf.modified$dist.Left, 
                               gtf.modified$V4))

gtf.modified$V5 = floor(ifelse(gtf.modified$modify.Right == "Y", 
                               gtf.modified$V5 + gtf.modified$dist.Right, 
                               gtf.modified$V5))

n.unique = length(unique(gtf.modified[gtf.modified$modify.Left == "Y" | gtf.modified$modify.Right == "Y",]$gene))
print(paste0("Modifying ", n.unique, " genes..."))

gtf.modified = rbind(gtf.modified, x.utr, x.exon)
gtf.modified = gtf.modified[order(gtf.modified$V1, gtf.modified$V4),]
gtf.final = gtf.modified[, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")]
#gtf.final.final = rbind(gtf.final, rest[, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")])

write.table(gtf.final, file = "/Users/ryan.gorzek/Downloads/Monodelphis_domestica.ASM229v1.112.filtered.lt10kb.ex5kb.gtf", 
            quote = F, row.names = F, col.names = F, sep = "\t")
