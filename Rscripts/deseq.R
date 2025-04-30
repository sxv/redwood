#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

# DESEQ2
library("DESeq2")
read_counts <- read.csv("/Users/x2a/helix/udh/udh_dcis_nl_n28.reads.csv", sep=",", row.names=1)
design <- read.csv("/Users/x2a/helix/udh/udh_dcis_nl_n28.design.csv", header=TRUE)

design$patient <- as.factor(design$patient)

design$dx <- as.factor(design$dx)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~ patient + dx)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, contrast=c("dx", "UDH", "DCIS"), alpha=.05)
res <- res[which(res$padj < 0.05),]
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/udh/UDHvDCIS.csv")


design$udh <- as.factor(design$udh)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~ patient + udh)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, contrast=c("udh", "dx1", "dx0"), alpha=.05)
res <- res[which(res$padj < 0.05),]
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/udh/UDHvALL.csv")

design$dcis <- as.factor(design$dcis)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~patient + dcis)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, contrast=c("dcis", "dx1", "dx0"), alpha=.05)
res <- res[which(res$padj < 0.05),]
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/udh/DCISvALL.csv")

design$nl <- as.factor(design$nl)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~ patient + nl)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, contrast=c("nl", "dx1", "dx0"), alpha=.05)
res <- res[which(res$padj < 0.05),]
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/udh/NLvALL.csv")

design$udhvnl <- as.factor(design$udhvnl)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~ udhvnl)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, contrast=c("udhvnl", "dx1", "dx0"), alpha=.05)
res <- res[which(res$padj < 0.05),]
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/udh/UDHvNL.csv")

design$udhvnl <- as.factor(design$udhvnl)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~ patient + udhvnl)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, contrast=c("udhvnl", "dx1", "dx0"), alpha=.05)
res <- res[which(res$padj < 0.05),]
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/udh/UDHvNLwithP.csv")
