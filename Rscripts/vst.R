library(DESeq2)

# LOAD raw reads and design file
edata <- read.delim("/Users/x2a/helix/comet/n1457/comet_reads_n1457_v33.csv", sep=",", row.names=1, check.names=FALSE)
design.file <- read.csv('/Users/x2a/helix/comet/n1457/comet_design_n1457.csv', sep=",", header=TRUE)

dds <- DESeqDataSetFromMatrix(countData = edata, colData = design.file, design = ~1 )
dds <- DESeq(dds, parallel=TRUE)

# alternative to vst
#normalized_counts <- counts(dds, normalized=TRUE)

# unneceessary step as this is performed by vst()
#dds <- estimateSizeFactors(dds)

vsd <- varianceStabilizingTransformation(dds)
write.csv(round(assay(vsd), 2), file="/Users/x2a/helix/comet/n1457/comet_reads_n1457_v33.VST.csv")


# Filter the design file for samples with dx_new == "DCIS"
dds_filtered <- dds[, design.file$dx_new == "DCIS"]

# Ensure the colData is updated after subsetting
dds_filtered <- dds_filtered[rowSums(counts(dds_filtered)) > 0, ] # Remove rows with all zeros
dds_filtered <- DESeq(dds_filtered, parallel=TRUE)

# Apply the variance stabilizing transformation to the filtered dataset
vsd_filtered <- varianceStabilizingTransformation(dds_filtered)

# Write the transformed counts to a CSV file
write.csv(round(assay(vsd_filtered), 2), file="/Users/x2a/helix/comet/n1457/comet/comet_reads_n1457_v33_.VST.DCIS_n765.csv")
