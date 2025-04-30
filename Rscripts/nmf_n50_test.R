### RNA NMF #######################################################################################
# Run NMF to cluster DCIS RNA samples

### PREAMBLE ######################################################################################
# load libraries
#library(argparse)
#library(viridis)
#library(DESeq2)
#BiocManager::install("fgsea")
#library(fgsea)
#library(tidyr)
library(NMF)
#library(rms)

rna_matrix_3<-read.csv("./bph_epi_nod_n50_reads.VST.csv", sep=",", row.names=1)

edata <- t(scale(t(rna_matrix_3),center = FALSE))

var.data <- apply(edata, 1, var)
vfilter <- 0
vdata <- edata[which(var.data > vfilter ),]
#dim(vdata)

setwd("./")

# run NMF from 2 to 10 factors, 30 times
# takes time and memory
res <- nmf(vdata, 2:10, nrun=30, method="brunet", .opt="vp10")


save(res, file = './bph_epi_nod_n50_reads.VST.RData')

