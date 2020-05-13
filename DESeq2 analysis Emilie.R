###### RNASeq analysis ferroptosis #####

#DESEQ2 package


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library ("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("S4Vectors")
library("S4Vectors")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DelayedArray")
library("DelayedArray")

# Other useful packages

library("reshape2")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("pheatmap")
library("genefilter")
library("dplyr")
library("ggrepel")

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

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~cell_line + treatment)

dds


# explore this data with clustering or pca

rld <- rlog(dds)

distsRL <- dist(t(assay(rld))) 

topVarianceGenes <- head(order(-rowVars(assay(rld))), 275)
matrix <- assay(rld)[topVarianceGenes,]
matrix <- matrix - rowMeans(matrix)
df <- as.data.frame(colData(rld)[,c("treatment", "cell_line")])
rownames(df) <-colnames(matrix)

pheatmap(matrix, annotation_col=df)



mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(treatment)) 
colnames(mat) = NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
pheatmap(mat,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=colors)


plotPCA(rld, intgroup=c("cell_line"))
plotPCA(rld, intgroup=c("treatment"))
plotPCA(rld, intgroup=c("replicate"))
plotPCA(rld, intgroup =c("treatment", "cell_line"))


plot <-plotCounts(dds, "NFKBIA", intgroup = "treatment", normalized = TRUE, transform = FALSE, xlab = "treatment", returnData = TRUE, replaced = FALSE, col = dds$cell_line)
ggplot(plot, aes(x=treatment, y=count, color = dds$cell_line)) + geom_point(position=position_jitter(w = 0.1,h = 0)) + theme_bw() + ggtitle("Flt3") + theme(plot.title = element_text(hjust = 0.5))

#DEG analysis


dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds, test="Wald")

counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(counts), file = "/Users/ELogie/Desktop/Emilie/Experimenten/Emilie/RNA seq/Ferroptose/RNAseq ferroptosis ANALYSIS/RESULTS_WA/DESeq2_normalized_readcountsWA.txt")
write.csv(assay(rld), file="rlogCountRSL3vsWA.txt", sep="/t")

res <- results(dds,name=, contrast = c("treatment", "RSL3", "untreated"))
res

res <- res[order(res$padj),]
mcols(res, use.names=TRUE)
write.csv(as.data.frame(res), file = "/Users/emili/OneDrive/Documenten/Emilie/RNAseq ferroptosis ANALYSIS/DESeq2_results_group-RSL3vsWA.txt")


results = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res))
head(results)

DEgenes_DESeq <- results[which(abs(results$log2FoldChange) > 1.0 & results$padj < 0.05),]
write.csv(as.data.frame(DEgenes_DESeq), file = "/Users/emili/OneDrive/Documenten/Emilie/RNAseq ferroptosis ANALYSIS/DESeq2_results_group-RSL3vsWA_filtered.txt")




##### Data Visualization ####


# MA plot

plotMA(res, ylim=c(-7,7))





# Volcano plot 1 # bijzetten script Jöran


data_wa <- read.table(file="DESeq2_results_group-RSL3vsWA.txt", header = T, sep = ",", dec=".")


colnames(data_wa) <-c("Genes", "baseMean", "log2FoldChange", "lfcSE","stat", "pvalue", "padj")

rownames(data_wa) <- data_wa$Genes

par(mfrow=c(1,2))

plot(data_wa$log2FoldChange, -log10(data_wa$padj),
     pch=16,
     xlab="logFC",
     ylab="-log10(FDR)",
     main = "RSL3 vs WA",
     xlim = c(-6,15),
     ylim = c(0,50),
     col=ifelse(data_wa$padj < 0.05 & abs(data_wa$log2FoldChange) > 1, "darkgoldenrod", "black"))

legend(6, 45, legend = c("logFC >1 & padj <0,05", "not significant"), col = c("darkgoldenrod", "black"), pch=16)

with(subset(data_wa, ))

text(data_wa$log2FoldChange, -log10(data_wa$padj), labels= ifelse(data_wa$padj < 0.05 & abs(data_wa$log2FoldChange) > 1,rownames(data_wa), "NA"), cex=0.7, pos=3)


text(data_wa$log2FoldChange, -log10(data_wa$padj), labels= data_wa$Genes, cex=0.7, pos=3)




# Volcano plot 1.b
DEgenes_DESeq <- data_wa[which(abs(data_wa$log2FoldChange) > 1.5 & data_wa$padj < 0.0001),]
NODEgenes <- data_wa
plot(DEgenes_DESeq$log2FoldChange, -log10(DEgenes_DESeq$padj),
     pch=16,
     xlab="logFC",
     ylab="-log10(FDR)",
     main = "RSL3 vs WA",
     xlim = c(-6,15),
     ylim = c(0,50),
     col= "darkgoldenrod")

text(DEgenes_DESeq$log2FoldChange, -log10(DEgenes_DESeq$padj), labels= ifelse(DEgenes_DESeq$padj < 0.05 & abs(DEgenes_DESeq$log2FoldChange) > 1,rownames(DEgenes_DESeq), "NA"), cex=0.7, pos=3)
text(DEgenes_DESeq$log2FoldChange, -log10(DEgenes_DESeq$padj), labels= DEgenes_DESeq$Genes, cex=0.7, pos=3)
points(data_wa$log2FoldChange, -log10(data_wa$padj))

legend(6, 45, legend = c("logFC >1.5 & padj <0,0001", "not significant"), col = c("darkgoldenrod", "black"), pch=c(16,1))



# Volcano plot 2

p = ggplot2::ggplot(results, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("blue", "black")) +
  ggplot2::ggtitle("Volcano Plot of DESeq2 analysis")

p + ggrepel::geom_text_repel(data=results[1:10, ], ggplot2::aes(label=rownames(results[1:10, ])))

p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))




# Volcano plot 3

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=c("MT2A", "MT1G", "MT1X", "PDGFB", "BCL2L11", "MT1F", "SIK1B", "NR4A3"), cex=.8))



# explore this data with clustering or pca

rld = rlog(dds, blind=F)
distsRL <- dist(t(assay(rld))) 

topVarianceGenes <- head(order(-rowVars(assay(rld))), 100)
matrix <- assay(rld)[topVarianceGenes,]
matrix <- matrix - rowMeans(matrix)
df <- as.data.frame(colData(rld)[,c("replicate","cell_line")])
rownames(df) <-colnames(matrix)

matrix_new <- matrix[,c(1:9, 13, 23, 30)]
matrix_RSL3 <- matrix_new[,c(1,5,7,10,11,12)]

pheatmap(matrix_RSL3, annotation_col=df, color=brewer.pal(7, "YlGnBu"), scale = "none", cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean")

#Remark: heatmap above does not yet filter cell line specific effects: better to make heatmap from output DeSeq analysis

sampleDistMatrix <- as.matrix(distsRL)
rownames(sampleDistMatrix) <- paste(rld$treatment,sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=distsRL,
         clustering_distance_cols=distsRL,
         col=colors)




