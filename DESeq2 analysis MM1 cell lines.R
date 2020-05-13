###### RNASeq analysis ferroptosis #####


##### Install packages #####

#DESEQ2 package


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library ("DESeq2")


# Other useful packages

library("reshape2")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("pheatmap")
library("genefilter")


##### Read in the data #####


#Read in data

count_data <- read.table(file="all_readcounts.tsv", header=T, sep="\t", row.names = 1)
count_data <- count_data[-c(1:4),]
count_data <- count_data[,-(31:32)] 


col_data <- read.table(file="MetaData_RNAseq_ferroptosis.txt", header=T, sep="\t")

col_data$cell_line <- as.factor(col_data$cell_line)
col_data$treatment <- as.factor(col_data$treatment)
col_data$replicate <- as.factor(col_data$replicate)


# Inspect the data

head(col_data)
str(col_data)
head(count_data)
str(count_data)



##### DESeq2 analysis #####

# read in the data as DESeq data file

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~treatment + cell_line)

dds


# explore this data with clustering or pca

rld = rlog(dds)
distsRL <- dist(t(assay(rld))) 

topVarianceGenes <- head(order(-rowVars(assay(rld))), 500)
matrix <- assay(rld)[topVarianceGenes,]
matrix <- matrix - rowMeans(matrix)
df <- as.data.frame(colData(rld)[,c("cell_line")])
rownames(df) <-colnames(matrix)


pheatmap(matrix, annotation_col=df)



mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(treatment)) 
colnames(mat) = NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(mat,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=colors)


plotPCA(rld, intgroup=c("cell_line"))
plotPCA(rld, intgroup=c("treatment"))
plotPCA(rld, intgroup=c("replicate"))
plotPCA(rld, intgroup =c("treatment", "cell_line"))


plotCounts(dds, "CDKN2A", intgroup = "treatment", normalized = TRUE, transform = FALSE, main = "Counts CDKN2A gene", xlab = "treatment", returnData = FALSE, replaced = FALSE, col = dds$cell_line)


#DEG analysis


dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds, test="Wald")

counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(counts), file = "/Users/emili/OneDrive/Documenten/Emilie/RNAseq ferroptosis ANALYSIS/DESeq2_normalized_readcountsMM1.txt")
write.csv(assay(rld), file="rlogCountMM1.txt", sep="/t")

res <- results(dds,name=, contrast = c("cell_line", "MM1R", "MM1S"))
res

res <- res[order(res$padj),]
mcols(res, use.names=TRUE)
write.csv(as.data.frame(res), file = "/Users/emili/OneDrive/Documenten/Emilie/RNAseq ferroptosis ANALYSIS/DESeq2_results_MM1RvsMM1S_cell_line2.txt")



#MA plot

plotMA(res, ylim=c(-7,7))


# Volcano


data_wa <- read.table(file="DESeq2_results_MM1RvsMM1S_cell_line2.txt", header = T, sep = ",", dec=".")

colnames(data_wa) <-c("Genes", "baseMean", "log2FoldChange", "lfcSE","stat", "pvalue", "padj")
rownames(data_wa) <- data_wa$Genes

par(mfrow=c(1,2))

plot(data_wa$log2FoldChange, -log10(data_wa$padj),
     pch=16,
     xlab="logFC",
     ylab="-log10(FDR)",
     main = "MM1R vs MM1S",
     xlim = c(-10,10),
     col=ifelse(data_wa$padj < 0.05 & abs(data_wa$log2FoldChange) > 1, "deeppink", "black"))

legend(-10, 200, legend = c("logFC >1 & padj <0,05", "not significant"), col = c("deeppink", "black"), pch=16)



# Volcano plot 1.b
DEgenes_DESeq <- data_wa[which(abs(data_wa$log2FoldChange) > 4 & data_wa$padj < 0.00001),]
NODEgenes <- data_wa
plot(DEgenes_DESeq$log2FoldChange, -log10(DEgenes_DESeq$padj),
     pch=16,
     xlab="logFC",
     ylab="-log10(FDR)",
     main = "MM1R vs MM1S",
     xlim = c(-10,10),
     ylim = c(0,240),
     col= "deeppink")

text(DEgenes_DESeq$log2FoldChange, -log10(DEgenes_DESeq$padj), labels= ifelse(DEgenes_DESeq$padj < 0.05 & abs(DEgenes_DESeq$log2FoldChange) > 1,rownames(DEgenes_DESeq), "NA"), cex=0.7, pos=3)
points(data_wa$log2FoldChange, -log10(data_wa$padj))

legend(-10, 220, legend = c("logFC >4 & padj <0,00001", "not significant"), col = c("deeppink", "black"), pch=c(16,1))




# Heatmap 2

alpha <- 0.001
epsilon <- 1
gene.kept <- rownames(res)[res$padj< 0.01 & !is.na(res$padj) & abs(res$log2FoldChange)>1]

gene.kept <- gene.kept[1:250]  #TOP250 DEGs

count.table.kept <- log2(counts + epsilon)[gene.kept,]
count.table.kept.UT <- count.table.kept[,c(2,3,4,6,8,9)]
colnames(count.table.kept.UT) <- c("MM1S", "MM1S", "MM1R", "MM1R", "MM1S", "MM1R")
dim(count.table.kept)
colfunc <- colorRampPalette(c("yellow", "blue"))
colCols <- ifelse(colnames(count.table.kept.UT)>"MM1S", "purple",
                           ifelse(colnames(count.table.kept.UT)>"MM1R","gray85", "gray44"))
heatmap.2(as.matrix(count.table.kept.UT), col=colfunc,
          scale="row", trace = "none", density="none", cexCol = 0.7, ColSideColors = colCols, labCol = T, labRow = "")




hclust = function(x) hclust(x, method = "average"),
distfun= function(x) as.dist((1-cor(t(x)))/2),

labRow = "",
cexCol=0.7,
na.rm=T,
ColSideColors = colCols, labCol = F, Rowv=T)
