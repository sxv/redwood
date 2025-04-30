# Load DESeq2
library("DESeq2")

# Load data
read_counts <- read.csv("/Users/x2a/helix/ozge/dod_n30_reads.csv", sep=",", row.names=1)
design <- read.csv("/Users/x2a/helix/ozge/dod_n30_design.csv", header=TRUE)

design$group <- factor(design$group)  # Ensure it's a factor
# Ensure `group` has the reference level set to "Control"
design$group <- relevel(design$group, ref = "Control")

# Create DESeqDataSet (simplified design)
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = design, design = ~ group)

# Estimate size factors and run DESeq
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)

# Extract results for the Autism vs Control comparison
res <- results(dds, contrast=c("group", "Autism", "Control"), alpha=0.05)

# Filter significant results (adjusted p-value < 0.05)
res <- res[which(res$padj < 0.05),]

# Export results
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/ozge/deseq_group_autism_vs_control.csv")














# Load libraries
library("DESeq2")
library("dplyr")
library(tibble)

# Load data
read_counts <- read.csv("/Users/x2a/helix/ozge/dod_n30_reads.csv", sep=",", row.names=1)
design <- read.csv("/Users/x2a/helix/ozge/dod_n30_design.csv", header=TRUE)

# Ensure the columns of read_counts correspond to rows in the design file
if (!all(colnames(read_counts) == design$sample)) {
  stop("Column names of read_counts do not match design$sample.")
}

# Add `sub1` column to `read_counts` (repeat `sub1` for each sample column)
read_counts_long <- as.data.frame(t(read_counts))
read_counts_long$sub1 <- design$sub1

# Average counts by `sub1`
averaged_counts <- read_counts_long %>%
  group_by(sub1) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames(var = "sub1")

# Transpose back to the original format (genes as rows, samples as columns)
averaged_counts <- as.data.frame(t(averaged_counts))

# Simplify design matrix to one row per patient
design <- design %>%
  group_by(sub1) %>%
  summarise(group = first(group))

# Convert `group` to factor
design$group <- factor(design$group)
design$group <- relevel(design$group, ref = "Control")

# round them to integer
averaged_counts <- round(as.matrix(averaged_counts))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = as.matrix(averaged_counts), colData = design, design = ~ group)

# Estimate size factors and run DESeq
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=TRUE)

# Extract results for the Autism vs Control comparison
res <- results(dds, contrast=c("group", "Autism", "Control"), alpha=0.05)

# Filter significant results (adjusted p-value < 0.05)
res <- res[which(res$padj < 0.05),]

# Export results
write.csv(as.data.frame(res[order(res$padj),]), file="/Users/x2a/helix/ozge/deseq_group_autism_vs_control_averaged.csv")
