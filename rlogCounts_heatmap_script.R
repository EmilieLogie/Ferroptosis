
###### Make Heatmap of RSL3 genes #######


library("dplyr")
library("pheatmap")
library("RColorBrewer")

# read & convert data

DEG_RF <- read.csv("DESeq2_results_group-RSL3vsWA.txt"  , header = T, row.names = 1)
Genes <-row.names(DEG_RF)
DEG_RF_n <- cbind(Genes, DEG_RF)
rownames(DEG_RF_n) <- c()

significant_genes <-subset(DEG_RF_n, padj<.01 & abs(log2FoldChange)>1.5)


counts <-read.csv("rlogCountRSL3vsWA.txt", header=T, row.names=1, sep=",")
Genes <-row.names(counts)
counts_n <- cbind(Genes, counts)
rownames(counts_n) <- c()


total <- merge(significant_genes, counts_n, by="Genes")
#total_n <- total[,c(1,8:37)]

matrix_total <- as.matrix(total[,-1])
rownames(matrix_total) <- total[,1]
str(matrix_total)

matrix_total_order <- matrix_total[order(matrix_total[,6]),]
matrix_final <- matrix_total_order[,c(7:36)]
matrix_final_n <- matrix_final[c(1:120),c(1,5,7, 11,13,16,19,20,23,27,28,30)]



# heatmap
pheatmap(matrix_final_n, color=c(brewer.pal(7, "PuOr")), fontsize_row = 5, fontsize_col = 8)
 

