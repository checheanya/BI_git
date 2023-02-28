# Installing the packages -----------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "DESeq2"))
BiocManager::install("apeglm")
BiocManager::install('EnhancedVolcano')

install.packages('gplots')

library(gplots)
library(apeglm)
library(DESeq2)
library(edgeR)
library(scales)
library(EnhancedVolcano)
library(dplyr)
library("ggplot2") 
library("ggrepel")


# Analysing data --------------------------------------------------------------

counts <- read.table("simple_counts.txt", header = TRUE, row.names = 1, check.names = FALSE)
groups <- read.table("ann.txt", header = TRUE)
expr_group <- "fermentation_30_min"
ctrl_group <- "control_0_min"


groups$Group <- relevel(factor(groups$Group), ref = ctrl_group)
rownames(groups) <- groups$Sample
groups$Sample <- NULL

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = groups,
                              design = ~ Group)

dds <- DESeq(dds)

res <- lfcShrink(dds, coef = paste("Group_", expr_group, "_vs_", ctrl_group, sep=""), type = "apeglm")
res <- res[order(res$padj),]

# saving the table with data
write.table(res, file = "results.tsv", sep = "\t", quote = FALSE)


# Plotting --------------------------------------------------------------------

# OPTION 1 +++++++++++++++++++++++

#Adjusted P values
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~p~adj)))

with(subset(res, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#with(subset(res, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(res), res$padj<0.001 & abs(res$log2FoldChange)>5), cex=0.8, pos=3))

# Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(res$pvalue[res$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

# OPTION 2  +++++++++++++++++++++++

df <- as.data.frame(res)
names <- rownames(res)
names_norm <- c()
for(name in names){
  namen <- substr(name, 6, nchar(name))
  names_norm <- append(names_norm, namen)
}

df$genes <- names_norm
# Will have different colors depending on significance
mutated_df <- mutate(df, significance = ifelse(df$padj < 0.01, "p_adj < 0.01", "Not significant"))
# convert the rownames to a column
input <- cbind(gene=rownames(mutated_df), mutated_df) 

#volcanoplot with log2Foldchange versus pvalue
volc <- ggplot(input, aes(log2FoldChange, -log10(pvalue)))+ 
  geom_point(aes(col=significance)) + #add points colored by significance
  scale_color_manual(values=c("black", "red")) + 
  ggtitle("Volcano plot for the fermentation data") #e.g. 'Volcanoplot DESeq2'

volc + geom_text_repel(data = head(input, 20), aes(label=gene)) #adding text for the top 20 genes
#ggsave("Volcanoplot.jpeg", device="jpeg") #In case you want to easily save to disk
volc


# OPTION 3  +++++++++++++++++++++++

EnhancedVolcano(res,
                lab = names_norm,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Difference before and after 30mins of fermentation',
                subtitle = '',
                pCutoff = 10e-32,
                FCcutoff = 5,
                pointSize = 2.0,
                labSize = 4.0)










# Scripts from teachers -------------------------------------------------------

# Set up the conditions based on the experimental setup.
cond_1 = rep("control", 2)
cond_2 = rep("30mins", 2)

# Read the data from the standard input.
countData = read.table("simple_counts.txt", header=TRUE, sep="\t", row.names=1 )

# Build the dataframe from the conditions
samples = names(countData)
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

#Set the reference to be compared
dds$condition = relevel(dds$condition,"control")

# Run deseq
dds = DESeq(dds)

# Format the results.
res = results(dds)

# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj, -log2FoldChange)), ]

# Turn it into a dataframe to have proper column names.
sorted.df = data.frame("id"=rownames(sorted),sorted)

# Write the table out.
write.table(sorted.df, file="result2.txt", sep="\t", col.names=NA, quote=FALSE)

# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)

# Save the normalize data matrix.
write.table(dt, file="norm-matrix-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)



# Script from teachers for plotting -------------------------------------------

# Read normalized counts
data = read.table("norm-matrix-deseq2.txt", header=T, sep="\t", as.is=TRUE)

gene = data[,1]
vals = as.matrix(data[,2:ncol(data)])

# Adds a little noise to each element
# To avoid the clusteing function failing on zero
# variance datalines.
vals = jitter(vals, factor = 1, amount=0.00001)


# Calculate zscore
score = NULL
for (i in 1:nrow(vals)) {
  row=vals[i,]
  zscore=(row-mean(row))/sd(row)
  score =rbind(score,zscore)
}

row.names(score) = gene
zscore=score

# Generate heatmap
mat = as.matrix(zscore)

# Opent the drawing device.
pdf('output.pdf')

colors = colorRampPalette(c("green","black","red"),space="rgb")(256)
heatmap.2(mat,col=colors,density.info="none",trace="none", margins=c(14,14),lhei=c(1,5))

invisible(dev.off())

s