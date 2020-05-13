
###### Venn Diagram ######

# Load library
library(VennDiagram)



# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)


# Read data
data_RSL3 <- read.csv("DESeq2_results_group-RSL3vsuntreated_filtered_genes.txt", header = F)
list_RSL3 <- as.vector(t(data_RSL3))
data_WA <- read.csv("DESeq2_results_group-WAvsuntreated_filtered_Genes.txt", header=F)
list_WA <- as.vector(t(data_WA))


# Chart

myCol= c("#ffbc42","#218380")
venn.plot <-venn.diagram(list(list_RSL3, list_WA), filename='Overlap_RSL3_WA.png', category.names=c("RSL3", "WA"), output=T, 
                         lwd=2, lty='blank',fill=myCol, cex=0.6, fontfamily="sans", imagetype = "png", height=700, width=700, resolution =300, compression = "lzw",
                         cat.cex=0.6, cat.fontface="bold", cat.default.pos="outer", cat.fontfamily="sans")

Common_genes<- intersect(list_RSL3, list_WA)
Common_genes <- as.data.frame(Common_genes, row.names=NULL)
write.csv(Common_genes, file="Overlap_WA_RSL3.txt",quote=F, row.names = F)
