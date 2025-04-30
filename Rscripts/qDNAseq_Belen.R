#!/usr/bin/env Rscript

# Quantitative DNA sequencing for chromosomal aberrations. The genome is divided into non-overlapping fixed-sized bins, 
# number of sequence reads in each counted, adjusted with a simultaneous two-dimensional loess correction for sequence mappability 
# and GC content, and filtered to remove spurious regions in the genome.
# Downstream steps of segmentation and calling are also implemented via packages DNAcopy and CGHcall, respectively.

## Author: Aziz Khan <azizk@stanford.edu>
## Version 1.0

##### Run this script from the terminal

### bam directory lab@path-liza2:/media/stroma-common/Belen/LVI_hg38_trimmed
### qdnaseq results: lab@path-liza2:/media/drive3/belen/qdnaseq.lvi_trimmed_2.15.21


### for b in *bam; do Rscript <script R> -f ${b} -n ${b}; done

start.time <- Sys.time()


library("optparse")
library(tidyverse)
library(Biobase) 
library("QDNAseq")
library(CGHbase)
future::plan("multicore")
library(data.table)
library(ACE)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Input BAM file or path", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./", 
              help="output path name [default= %default]", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="sample1", 
              help="sample bam file name/id [default= %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="hg38", 
              help="Genome (hg19, hg38) [default= %default]", metavar="character"),
  make_option(c("-m", "--minmapq"), type="integer", default=37, 
              help="Minimum quality score [default= %default]", metavar="integer"),## where is this coming from
  make_option(c("-b", "--bin"), type="integer", default=50,
              help="bin size (kb). Use 5, 10, 15, 30, 50, 100, 200, 500, 1000 [default= %default]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt$out = "/media/drive3/belen/qdnaseq.lvi_trimmed_2.15.21"
opt$name <- gsub("\\.sort.marked_duplicates.bam$", "", opt$name)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input bam file).n", call.=FALSE)
}

sink(paste0(opt$out,opt$name,"_",opt$bin,"_log.txt"))

if (opt$genome == 'hg38'){
  ## https://github.com/asntech/QDNAseq.hg38
  library("QDNAseq.hg38")
  ##
  ##Install the QDNAseq.hg38 package using remotes
  #remotes::install_github("asntech/QDNAseq.hg38@main")
  ##or devtools
  # devtools::install_github("asntech/QDNAseq.hg38@main")
  
  bins <- getBinAnnotations(binSize=opt$bin, genome="hg38")
}else{
  # library("QDNAseq.hg19")
  # bins <- getBinAnnotations(binSize=opt$bin)
}

outputPath = paste0(opt$out,opt$name,"_",opt$bin,"_")

probability_matrix <- function(cgh) {
  probs = matrix(0, dim(cgh)[1], dim(cgh)[2])
  colnames(probs) = colnames(cgh)
  rownames(probs) = rownames(cgh)
  
  for (i in 1:dim(cgh)[1]){
    for (j in 1:dim(cgh)[2]){
      if (calls(cgh)[i,j] == 2){
        probs[i,j] = probamp(cgh)[i,j]
      } else if (calls(cgh)[i,j] == 1){
        probs[i,j] = probgain(cgh)[i,j]
      } else if (calls(cgh)[i,j] == -1){
        probs[i,j] = probloss(cgh)[i,j]
      } else if (calls(cgh)[i,j] == -2){
        probs[i,j] = probdloss(cgh)[i,j]
      }
    }
  }
  
  probs
}

getSegTable <- function(x)
{
  dat<-x
  sn<-assayDataElement(dat,"segmented")
  fd <- fData(dat)
  fd$use -> use
  fdfiltfull<-fd[use,]
  sn<-sn[use,]
  segTable<-c()
  for(c in unique(fdfiltfull$chromosome))
  {
    snfilt<-sn[fdfiltfull$chromosome==c]
    fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
    sn.rle<-rle(snfilt)
    starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
    ends <- cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from <- fdfilt$start[starts[s]]
      to <- fdfilt$end[ends[s]]
      segValue <- sn.rle$value[s]
      c(fdfilt$chromosome[starts[s]], from, to, segValue)
    }) -> segtmp
    segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
    segTable<-rbind(segTable,segTableRaw)
  }
  colnames(segTable) <- c("chrom", "start", "end", "segVal")
  return(segTable)
}

if (file_test("-f", opt$file)){
  
  bams <- c(opt$file)
}else {
  bams <- paste0(opt$file, list.files(opt$file, pattern=".bam", recursive=F))
  
}

if (file_test("-f", paste0(opt$out,opt$name,"_",opt$bin,"_copyNumbersSegmented.rds"))){
  copyNumbersSegmented <- readRDS(paste0(opt$out,opt$name,"_",opt$bin,"_copyNumbersSegmented.rds"))
} else {
  
  if (file_test("-f", paste0(opt$out,opt$name,"_",opt$bin,"_readCountsFilteredCorrected.Rdata"))){
    copyNumbers <- readRDS(paste0(opt$out,opt$name,"_",opt$bin,"_readCountsFilteredCorrected.Rdata"))
  } else {
    
    readCounts <- binReadCounts(bins,  bamfiles=bams, minMapq=opt$minmapq)
    saveRDS(readCounts, file=paste0(opt$out,opt$name,"_",opt$bin,"_readCounts.Rdata"))
    
    pdf(paste0(opt$out,opt$name,"_",opt$bin,"_rawprofile.pdf"), width=20, height=8, useDingbats=F)
    ## Plot the readcounts with filtered reads highlighted.
    plot(readCounts, logTransform=TRUE)
    highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
    dev.off()
    
    ##Apply QDNAseq filters.
    readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE) ##default function excludes sex chrs
    # fix for issue with copy number segmentation arising from zero count bins
    # plot median read counts as a function of GC content and mappability as an isobar plot
    pdf(paste(outputPath,"isobar.pdf",sep=''))
    isobarPlot(readCountsFiltered)
    dev.off()
    
    ##Calculate CG correction.
    readCountsFiltered <- estimateCorrection(readCountsFiltered)
    pdf(paste(outputPath,"noise.pdf",sep=''))
    noisePlot(readCountsFiltered)
    # noise plot showing relationship between the observed standard deviation in the data
    # and its read depth
    dev.off()
    
    ##Apply GC correction.
    copyNumbers <- correctBins(readCountsFiltered)
    
    saveRDS(copyNumbers, file=paste0(opt$out,opt$name,"_",opt$bin,"_readCountsFilteredCorrected.Rdata"))
  }
  
  
  ##Normalise and smooth outliers.
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  
  pdf(paste0(opt$out,opt$name,"_",opt$bin,"_copyNumbersSmooth.pdf"), width=20, height=8, useDingbats=F)
  ##Plot the smoothed copy-number.
  plot(copyNumbersSmooth)
  dev.off()
  
  exportBins(copyNumbersSmooth, file=paste0(opt$out,opt$name,"_",opt$bin,"_CNA_smooth.txt"))
  exportBins(copyNumbersSmooth, file=paste0(opt$out,opt$name,"_",opt$bin,"_CNA_smooth.igv"), format="igv")
  exportBins(copyNumbersSmooth, file=paste0(opt$out,opt$name,"_",opt$bin,"_CNA_smooth.bed"), format="bed")
  
  # output raw and fitted read counts
  features <- fData(copyNumbersSmooth) %>%
    as.data.frame %>%
    rownames_to_column(var = "location") %>%
    transmute(location, chrom = chromosome, start = as.integer(start), end = as.integer(end))
  
  rawReadCounts <- assayData(copyNumbersSmooth)$counts %>%
    as.data.frame %>%
    rownames_to_column(var = "location")
  features %>%
    dplyr::left_join(rawReadCounts, by = "location") %>%
    write_tsv(paste(outputPath,"rawReadCounts.txt",sep=''))
  
  fittedReadCounts <- assayData(copyNumbersSmooth)$fit %>%
    as.data.frame %>%
    rownames_to_column(var = "location") %>%
    mutate_if(is.numeric, list(~round(., digits = 3)))
  features %>%
    dplyr::left_join(fittedReadCounts, by = "location") %>%
    write_tsv(paste(outputPath,"fittedReadCounts.txt",sep=''))
  
  ## Segment the copy-number profile. 
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
  #copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  
  pdf(paste0(opt$out,opt$name,"_",opt$bin,"_CN_segmented.pdf"), width=20, height=8, useDingbats=F)
  ##Plot the segmented profile.
  plot(copyNumbersSegmented)
  
  pData(copyNumbersSegmented) %>%
    rownames_to_column(var = "sample") %>%
    dplyr::select(id = name, everything()) %>%
    write_tsv(paste0(opt$out,opt$name,"_",opt$bin,"_SegmentedSummary.txt"))
  
  saveRDS(copyNumbersSegmented, file=paste0(opt$out,opt$name,"_",opt$bin,"_copyNumbersSegmented.rds"))
}

## Run ACE
if (!dir.exists(paste0(opt$out,opt$name,"_",opt$bin,"_copyNumbersSegmented/"))){
  ploidies = c(2,3,4)
  runACE(inputdir = opt$out , outputdir= opt$out , filetype = 'rds', genome = opt$genome,
         binsizes=50, ploidies = ploidies, imagetype = 'pdf', method = 'RMSE', penalty = 0, 
         cap = 12, bottom = 0, trncname = ".recal", printsummaries = TRUE, 
         autopick = TRUE)
}

# runACE(inputdir = './', outputdir='./', filetype = 'rds', genome = opt$genome,
#   binsizes=50, ploidies = ploidies, imagetype = 'pdf', method = 'RMSE', penalty = 0, 
#   cap = 12, bottom = 0, trncname = ".recal", printsummaries = TRUE, 
#   autopick = TRUE)
# }

##QDNAseq use CGHcall for CN calls
#CGHcall(inputSegmented, prior = "auto", nclass = 5, organism = "human", cellularity=1, robustsig="yes", nsegfit=3000, maxnumseg=100, minlsforfit=0.5, build="GRCh37",ncpus=1)

##get cellularity from ACE
fitpicker_2N <- read.delim(paste0(opt$out,opt$name,"_",opt$bin,"_copyNumbersSegmented/2N/fitpicker_2N.tsv"))
fitpicker_3N <- read.delim(paste0(opt$out,opt$name,"_",opt$bin,"_copyNumbersSegmented/3N/fitpicker_3N.tsv"))
fitpicker_4N <- read.delim(paste0(opt$out,opt$name,"_",opt$bin,"_copyNumbersSegmented/4N/fitpicker_4N.tsv"))

##pick max cellularity from 2n,3n,4n
cellularity <- max(fitpicker_2N$likely_fit,fitpicker_3N$likely_fit,fitpicker_4N$likely_fit)

## Call aberrations from segmented copy number data
# opt$out = "/Volumes/GoogleDrive/My Drive/Postdoc/BCRF/QDNAseq_results_noChrX_10.13.21/"

if (opt$genome == 'hg38'){
  copyNumbersCalled <- callBins(copyNumbersSegmented, cellularity=cellularity, build="GRCh38")
}else{
  # copyNumbersCalled <- callBins(copyNumbersSegmented, cellularity=cellularity, build="GRCh37")
}

exportBins(copyNumbersCalled, file=paste0(opt$out,opt$name,"_",opt$bin,"_CNA_called.vcf"), format="vcf")
exportBins(copyNumbersCalled, file=paste0(opt$out,opt$name,"_",opt$bin,"_CNA_called.seg"), format="seg")

## Plot final profile.
plot(copyNumbersCalled)
dev.off()

pdf(paste0(opt$out,opt$name,"_",opt$bin,"_frequency_plot.pdf"), width=20, height=8, useDingbats=F)
frequencyPlot(copyNumbersCalled)
dev.off()

cgh <- makeCgh(copyNumbersCalled)
saveRDS(copyNumbersSegmented, file=paste0(opt$out,opt$name,"_",opt$bin,"_cgh_object.Rdata"))
write.table(probability_matrix(cgh), paste0(opt$out,opt$name,"_",opt$bin,"_probability_matrix.txt"), sep='\t', quote=FALSE, col.names=NA)
write.table(copynumber(cgh), paste0(opt$out,opt$name,"_",opt$bin,"_copynumber.txt"), sep='\t', quote=FALSE, col.names=NA)
write.table(calls(cgh), paste0(opt$out,opt$name,"_",opt$bin,"_calls.txt"), sep='\t', quote=FALSE, col.names=NA)
write.table(segmented(cgh), paste0(opt$out,opt$name,"_",opt$bin,"_segmented.txt"), sep='\t', quote=FALSE, col.names=NA)

## Output a table of all segments with a probability of loss is greater than 0.99
# filteredCN <- copyNumbersCalled[fData(copyNumbersCalled)$use,]
# regions_of_loss <- filteredCN[assayDataElement(filteredCN,"probloss")>0.99,]
# relative_loss <- getSegTable(regions_of_loss)
# write.table(relative_loss, file=paste0(opt$out,opt$name,"_",opt$bin,"_segments_relative_loss.txt"), quote = FALSE, row.names = FALSE)

#system(paste0("mv ",tail(stringr::str_split(opt$name, "/")[[1]],1),"  ", opt$out))
#system(paste0("mv ",opt$name,"_",opt$bin,"_CNA_called.seg ", opt$out))


# runACE('./', outputdir='./', filetype = 'rds', genome = 'hg38',
#   binsizes=50, ploidies = 3, imagetype = 'pdf', method = 'RMSE', penalty = 0.5, 
#   cap = 12, bottom = 0, trncname = FALSE, printsummaries = TRUE, 
#   autopick = FALSE)

sink()

######Segmentation using copynumber R package #####
#library(copynumber)



end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

rm(list = ls())

# filesstrings::move_files(files = list.files("/Volumes/GoogleDrive/My Drive/Postdoc/BCRF/QDNAseq_results_noChrX_10.13.21"
#                                             ,pattern = ".rds"), 
#                          destinations = "/Volumes/GoogleDrive/My Drive/Postdoc/BCRF/QDNAseq_results_noChrX_10.13.21/done/")
# 
# list.files(opt$out)

