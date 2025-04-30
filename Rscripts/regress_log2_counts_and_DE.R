#load data
COMET_design = read.csv('~/Google Drive/My Drive/DCIS_progression_classifier/old_pipeline_comet_design_n1393.csv')
COMET_data = read.csv('~/Google Drive/My Drive/DCIS_progression_classifier/old_pipeline_comet_reads_n1393_v33.csv')

rownames(COMET_data) = COMET_data$gene_name
COMET_data$gene_name = NULL
all(COMET_design$sample_id == colnames(COMET_data))

sum(rownames(COMET_data) %in% gene_ids)

#create Seurat object
library(Seurat)
COMET = CreateSeuratObject(counts = Matrix::Matrix(as.matrix(COMET_data), sparse = TRUE), meta.data = COMET_design)

table(COMET$source)

#remove controls
COMET_1 = subset(COMET, subset = source == 'Control', invert = T) 

table(COMET_1$source)

#normalize and pereprocess
COMET_1$log2counts = log2(COMET_1$nCount_RNA + 1)
COMET_1 <- NormalizeData(COMET_1)
COMET_1 <- FindVariableFeatures(COMET_1)
COMET_1 <- ScaleData(COMET_1, vars.to.regress = "log2counts")
COMET_1 <- RunPCA(COMET_1, features = VariableFeatures(object = COMET_1))
COMET_1 <- FindNeighbors(COMET_1, dims = 1:30, reduction = "pca")
COMET_1 <- FindClusters(COMET_1, resolution = 2, verbose = FALSE)
COMET_1 <- RunUMAP(COMET_1, dims = 1:30, reduction = "pca")

#plot UMAP
DimPlot(COMET_1, reduction = "umap", label = TRUE)
ggsave('~/Google Drive/My Drive/DCIS_progression_classifier/plots/by_seurat_clusters.png')

DimPlot(COMET_1, reduction = "umap", label = TRUE, group.by = 'type')
ggsave('~/Google Drive/My Drive/DCIS_progression_classifier/plots/by_type.png')

FeaturePlot(COMET_1, 'log2counts', cols = topo.colors(10))
FeaturePlot(COMET_1, 'log2counts')
ggsave('~/Google Drive/My Drive/DCIS_progression_classifier/plots/by_log2_counts.png')


#differential expression
markers <- FindAllMarkers(COMET_1, only.pos = TRUE, 
                          min.pct = 0.25, logfc.threshold = 0.25,
                          test.use = "DESeq2")

#select 15 top DE genes
library(magrittr)
library(dplyr)
top15 = markers %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC)

DoHeatmap(COMET_1, features = top15$gene) + NoLegend()

tail(top15)
#save data
write.csv(markers, '~/Google Drive/My Drive/DCIS_progression_classifier/merkers_DESeq2_all_samples.csv', row.names = F)
saveRDS(COMET_1, '~/Google Drive/My Drive/DCIS_progression_classifier/COMET_object_no_controls_library_size_regressed.RDS')


levels(COMET_1) = rev(c("NL",
                       "BBD",
                       "ADH",
                       "DCIS",
                       "ALH",
                       "LCIS"))

#plots
DotPlot(COMET_1, features = c('KIT','CD163','CSF1R','FOLR2','LYVE1','MARCO','F13A1', 'SELENOP','CFD',
                              'JCHAIN', 'IGHA1','IGLC2','IGKC','IGHM','IGLC3','IGHA2',
                                'SPI1','FCGR3A', 'FCGBP','CX3CR1','C3', 
                                'ERBB2','CD68','CSF1','CXCL9','CXCL11','CCR7','SPP1', 'CD19','MMP9','IL4I1',
                                'NLRP3', 'IL1B','CCL19','CXCL13')) + 
  scale_colour_gradient2(
    low = "blue",
    mid = "white",
    high = "red") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("") 

Idents(COMET_1) = COMET_1$type
Idents(COMET_1) = COMET_1$seurat_clusters

#subset DCIS
DCIS = subset(COMET_1, subset = type == 'DCIS')

#subcluster DCIS
DCIS <- NormalizeData(DCIS)
DCIS <- FindVariableFeatures(DCIS)
DCIS <- ScaleData(DCIS, vars.to.regress = "log2counts")
DCIS <- RunPCA(DCIS, features = VariableFeatures(object = DCIS))
DCIS <- FindNeighbors(DCIS, dims = 1:30, reduction = "pca")
DCIS <- FindClusters(DCIS, resolution = 2, verbose = FALSE)
DCIS <- RunUMAP(DCIS, dims = 1:30, reduction = "pca")

DimPlot(DCIS, reduction = "umap", label = TRUE)
DimPlot(DCIS, reduction = "umap", label = TRUE, group.by = 'type')
FeaturePlot(DCIS, 'log2counts', cols = topo.colors(10))

markers <- FindAllMarkers(DCIS, only.pos = TRUE, 
                          min.pct = 0.25, logfc.threshold = 0.25,
                          test.use = "MAST") #MAST, DESeq2

top20 = markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

DotPlot(DCIS, features = c('KIT','CD163','CSF1R','FOLR2','LYVE1','MARCO','F13A1', 'SELENOP','CFD',
                              'JCHAIN', 'IGHA1','IGLC2','IGKC','IGHM','IGLC3','IGHA2',
                              'SPI1','FCGR3A', 'FCGBP','CX3CR1','C3', 
                              'ERBB2','CD68','CSF1','CCR2','CXCL9','CXCL11','CCR7','SPP1', 'CD19','MMP9','IL4I1',
                              'NLRP3', 'IL1B','CCL19','CXCL13')) + 
  scale_colour_gradient2(
    low = "blue",
    mid = "white",
    high = "red") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("") 

FeaturePlot(DCIS, c('ERBB2','CXCL11','CXCL9','IL4I1','CCR2','CCR7','CCL19','SPP1'))
FeaturePlot(DCIS, c('ERBB2','CCR2','CXCL9'))

###Classify
library(DESeq2)

apply_variance_stabilizing_transformation <-
  function(rna_matrix) {
    # convert rna matrix to DESeq data set
    dds <- DESeqDataSetFromMatrix(
      countData = rna_matrix,
      data.frame(cond=rep(1,ncol(rna_matrix))),
      design=~1
    )
    # estimate library size
    dds <- estimateSizeFactors(dds)
    # extracts counts normalized by size factors
    #rna_matrix_norm <- counts(dds,normalized=T)
    # apply variance stablizing transformation to counts normalized by size factors
    vsd <-varianceStabilizingTransformation(dds)
    # extract matrix of transformed values
    rna_matrix_ep <- assay(vsd)
    return(rna_matrix_ep)
  }

genes <- read.csv('~/Google Drive/My Drive/DCIS_progression_classifier/TableS3.csv',
                  as.is = TRUE,
                  row.names = 1,
                  header = TRUE
)

gene_ids <- rownames(genes)

COMET_data_1 = COMET_data[,COMET_design$source != 'Control']

example_raw_data = COMET_data_1

rrna_matrix_1 <-
  apply_variance_stabilizing_transformation(example_raw_data)

rrna_matrix_2 <- t(scale(t(rrna_matrix_1)))
rrna_matrix_2[1:10,1:10]



sum(rownames(rrna_matrix_2) %in% gene_ids) #812
gene_ids[!gene_ids %in% rownames(rrna_matrix_2)]
rahbt_data <- as.data.frame(t(rrna_matrix_2[rownames(rrna_matrix_2) %in% gene_ids,]))
# predict COMET recurrence probabilities
#load model
load("~/Library/CloudStorage/GoogleDrive-mmatusia@stanford.edu/My Drive/DCIS_progression_classifier/rf_fit.RData")

###just predict by Magda
predictions_prob <- predict(rf_fit, rahbt_data, type ='prob')
rahbt_data[1:10,1:10]
COMET_predictions = data.frame(sample = rownames(rahbt_data), predictions_prob)
COMET_predictions[1:10,1:3]
colnames(COMET_predictions)[2:3] = c('no_recurrence','recurrence')
write.csv(COMET_predictions, "~/Google Drive/My Drive/DCIS_progression_classifier/COMET_predictions_MM_based_on_812gene_classifier.csv", row.names = F)

all(COMET_1$sample_id == COMET_predictions$sample)
COMET_1$classfier_Nrec = COMET_predictions$no_recurrence
COMET_1$classfier_rec = COMET_predictions$recurrence
###

#######DCIS only prediction
COMET_DCIS = COMET_data[,COMET_design$type == 'DCIS']

example_raw_data = COMET_DCIS

rrna_matrix_1 <-
  apply_variance_stabilizing_transformation(example_raw_data)

rrna_matrix_2 <- t(scale(t(rrna_matrix_1)))
rrna_matrix_2[1:10,1:10]

sum(rownames(rrna_matrix_2) %in% gene_ids) #812
rahbt_data <- as.data.frame(t(rrna_matrix_2[rownames(rrna_matrix_2) %in% gene_ids,]))
# predict COMET recurrence probabilities
#load model
load("~/Library/CloudStorage/GoogleDrive-mmatusia@stanford.edu/My Drive/DCIS_progression_classifier/rf_fit.RData")

###just predict by Magda
predictions_prob <- predict(rf_fit, rahbt_data, type ='prob')
rahbt_data[1:10,1:10]
COMET_predictions = data.frame(sample = rownames(rahbt_data), predictions_prob)
COMET_predictions[1:10,1:3]
colnames(COMET_predictions)[2:3] = c('no_recurrence','recurrence')

all(COMET_1$sample_id[COMET_1$type == 'DCIS'] == COMET_predictions$sample)
COMET_1$DCIS_only_classfier_Nrec = NA
COMET_1$DCIS_only_classfier_Nrec[COMET_1$type == 'DCIS'] = COMET_predictions$no_recurrence

FeaturePlot(COMET_1, features = 'DCIS_only_classfier_Nrec') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave('~/Google Drive/My Drive/DCIS_progression_classifier/plots/Classifier_on_DCIS_only_Nonrecurrence_prob.png')

COMET_1$DCIS_only_classfier_rec = NA
COMET_1$DCIS_only_classfier_rec[COMET_1$type == 'DCIS'] = COMET_predictions$recurrence

FeaturePlot(COMET_1, features = 'DCIS_only_classfier_rec') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave('~/Google Drive/My Drive/DCIS_progression_classifier/plots/Classifier_on_DCIS_only_recurrence_prob.png')

#######

library(RColorBrewer)

FeaturePlot(COMET_1, features = 'classfier_rec') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave('~/Google Drive/My Drive/DCIS_progression_classifier/plots/Classifier_on_allSamples_recurrence_prob.png')

FeaturePlot(COMET_1, features = 'classfier_Nrec') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave('~/Google Drive/My Drive/DCIS_progression_classifier/plots/Classifier_on_allSamples_NONrecurrence_prob.png')

DimPlot(COMET_1, group.by = 'type')
DimPlot(COMET_1, label = 10)
FeaturePlot(COMET_1, features = c('log2counts'))
FeaturePlot(COMET_1, features = c('LYVE1','FCGBP','SPP1','CXCL9'))
FeaturePlot(COMET_1, features = c('CCR2'))

FeaturePlot(COMET_1, features = c('ERBB2','CXCL9','CXCL11','CCR7'))

FeaturePlot(COMET_1, features = 'classfier_Nrec') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#boxplot prob progression by seurat cluster, by type
annot = data.frame(names = rownames(rrna_matrix_1)) #gene names

library(AnnotationDbi)
library(org.Hs.eg.db)

annot$EntrezGene.ID = mapIds(org.Hs.eg.db,
                             keys= as.character(annot$names),
                             column=c("ENTREZID"),
                             keytype = "SYMBOL",
                             multiVals = "first")

data = t(rrna_matrix_1)
dim(data)

PAM50Preds<-molecular.subtyping(sbt.model = "AIMS",data=data,
                                annot=annot,do.mapping=F)

all(rownames(PAM50Preds$subtype) == COMET_1$sample_id)
COMET_1$AIMS_prediction = PAM50Preds$subtype

DimPlot(COMET_1, group.by = 'AIMS_prediction')
ggsave('~/Google Drive/My Drive/DCIS_progression_classifier/plots/AIMS_molecular_subtype_classifier.png')

ggplot(COMET_1@meta.data, 
       aes(seurat_clusters, classfier_rec, color = seurat_clusters)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.4) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 15)) 


COMET_1$type = factor(COMET_1$type, levels = c('NL','BBD','ADH','DCIS','ALH','LCIS'))
ggplot(COMET_1@meta.data, 
       aes(type, classfier_rec, color = type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.4) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 15)) 


COMET_1$case = sapply(COMET_1$sample_id, function(x) strsplit(x,"_")[[1]][1])
COMET_1$case[COMET_1$classfier_rec > 0.5]
