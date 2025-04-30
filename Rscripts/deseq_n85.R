# DESEQ2
library("DESeq2")
read_counts <- read.csv("/Users/x2a/helix/bph2023/bph_str_n85_reads.csv", sep=",", row.names=1)
design <- read.csv("/Users/x2a/helix/bph2023/bph_str_n85_design.csv", header=TRUE)


design$patient <- as.factor(design$patient)
design$nod1rank2cluster <- as.factor(design$nod1rank2cluster)

design$nod1rank2cluster <- as.factor(design$nod1rank2cluster)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~ patient + nod1rank2cluster)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, contrast=c("nod1rank2cluster", "1", "2"), alpha=.05)
res <- res[which(res$padj < 0.05),]
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/bph2023/nod1/deseq_str_n85_nod1rank2cluster1v2.csv")

