library(Seurat)
library(sctransform)
library(ggplot2)

#BiocManager::install("glmGamPoi")

expmatrix <- read.csv("/Users/x2a/helix/udh/udh_dcis_nl_n28.reads.csv", sep=",", row.names=1, check.names=FALSE)
expmatrix <- read.csv("/Users/x2a/helix/tbcrc/tbcrc_rna3_n216.reads.csv", sep=",", row.names=1, check.names=FALSE)

expm <- CreateSeuratObject(counts = expmatrix)

# Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
expm <- SCTransform(object = expm, verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
expm <- RunPCA(object = expm,  npcs = 20, verbose = FALSE)
expm <- RunUMAP(object = expm, dims = 1:20, n.neighbors=15, verbose = FALSE)
expm <- FindNeighbors(object = expm, dims = 1:20, verbose = FALSE)
expm <- FindClusters(object = expm, verbose = FALSE)

# write out the cluster data to csv
#write.csv(as.data.frame(expm@meta.data), file="/Users/xv/helix/bph2023/bph_only_n160_clusters.csv")

# import metadata for use in plotting
metadata <- read.csv("/Users/x2a/helix/udh/udh_dcis_nl_n28.design.csv", header=TRUE)
metadata <- read.csv("/Users/x2a/helix/tbcrc/tbcrc_rna3_n216.design.csv", header=TRUE)

expm@meta.data$tissue <- metadata$t
expm@meta.data$dx <- metadata$dx
expm@meta.data$rna3 <- metadata$rna3
expm@meta.data$rna3_her2 <- metadata$rna3_her2

# plot umap, color coded by cluster, subtype etc
#png("/Users/xv/helix/bph2023/str_n67_12v345POS_umap_nod.png", 900, 600)

#DimPlot(object = expm, , pt.size=1.2)
#DimPlot(object = expm, , group.by="batch", pt.size=1.2)
#DimPlot(object = expm, , group.by="tissue", pt.size=1.2)
DimPlot(object = expm, , group.by="dx", pt.size=1.2)
DimPlot(object = expm, , group.by="tissue", pt.size=1.2)
DimPlot(object = expm, , group.by="rna3_her2", pt.size=1.2)

#DimPlot(object = expm, , group.by="umap", pt.size=1.2)

# with labels
#plot <- DimPlot(object = expm, , group.by="batch", pt.size = 5)
#LabelPoints(plot = plot, NULL, repel = TRUE, max.overlaps=500)
#graphics.off()

#FeaturePlot(expm, features = c("ESR1", "ERBB2", "MYC", "TP53", "PIK3CA"), cols=c("red", "blue") )
#FeaturePlot(expm, features = c("PIK3CA"))

# find and export all markers
#expm.markers <- FindAllMarkers(expm, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(as.data.frame(expm.markers), file="/Users/xv/helix/en_dcis_pairs/umap_markers_n42.csv")
