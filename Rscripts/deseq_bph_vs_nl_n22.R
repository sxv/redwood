# Load DESeq2
library("DESeq2")

# Load data
read_counts <- read.csv("/Users/x2a/helix/bph2025/2025-03-05_BPH_Str.n22.reads.csv", sep=",", check.names=FALSE)
rownames(read_counts) <- make.unique(read_counts[, 1], sep="_")  # Make first column unique
read_counts <- read_counts[, -1]  # Remove the first column since it's now row names


design <- read.csv("/Users/x2a/helix/bph2025/2025-03-05_BPH_Str.n22.design.csv", header=TRUE)

design$patient <- factor(design$patient)  # Ensure it's a factor
design$dx <- factor(design$dx)  # Ensure it's a factor

# Create DESeqDataSet (simplified design)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~ patient + dx)


# Estimate size factors and run DESeq
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)

# Extract results for the BPH vs NL comparison
res <- results(dds, contrast=c("dx", "BPH", "NL"), alpha=0.05)

# Filter significant results (adjusted p-value < 0.05)
res <- res[which(res$padj < 0.05),]

# Export results
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/bph2025/deseq_bph_vs_nl_patient_blocked_n22.csv")

