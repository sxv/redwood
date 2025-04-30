### RNA NMF #######################################################################################
# Run NMF to cluster DCIS RNA samples

### PREAMBLE ######################################################################################
# load libraries
library(NMF)


rna_matrix_3<-read.csv("/media/data/sxv/nmf/bcrf_rahbt_dcis_n142.reads.VST.csv", sep=",", row.names=1)

edata <- t(scale(t(rna_matrix_3),center = FALSE))

var.data <- apply(edata, 1, var)
vfilter <- 0
vdata <- edata[which(var.data > vfilter ),]
#dim(vdata)

setwd("/media/data/sxv/nmf/")

# run NMF from 2 to 10 factors, 30 times
# takes time and memory
res <- nmf(vdata, 2:10, nrun=30,.opt="vp20")


save(res, file = './bcrf_rahbt_dcis_n142.RData')

