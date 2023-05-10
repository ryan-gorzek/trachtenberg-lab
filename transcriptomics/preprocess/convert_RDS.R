
library(Seurat)
library(Matrix)

rds = readRDS("C:/Users/TLab/Downloads/Mouse_M1_10xV3_Matrix.RDS")
writeMM(rds, file="C:/Users/TLab/Downloads/Mouse_M1_10xV3_Matrix.csv")
write.table(names[1], file='C:/Users/TLab/Downloads/Mouse_M1_10xV3_Matrix_genes.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(names[2], file='C:/Users/TLab/Downloads/Mouse_M1_10xV3_Matrix_barcodes.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
