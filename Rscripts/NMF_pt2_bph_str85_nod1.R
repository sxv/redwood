### RNA NMF #######################################################################################
# Run NMF to cluster DCIS RNA samples

### PREAMBLE ######################################################################################
# load libraries
library(argparse)
library(viridis)
library(DESeq2)
#BiocManager::install("fgsea")
library(fgsea)
library(tidyr)
library(NMF)
library(rms)

load("/Users/x2a/helix/comet/NMF_comet_dcis_adh.RData")

design = read.csv('/Users/x2a/helix/comet/n1457/comet_design_n1457.csv', row.names=1)
design$group <- factor(design$dx_new)
head(design)

setwd("~/helix/comet/nmf_2025_02/")

if(requireNamespace("Biobase", quietly=TRUE)){
  # only compute the scores
  s <- featureScore()
  summary(s)
  # compute the scores and characterize each metagene
  s <- extractFeatures(res_dcis)
  str(s)
}

if(requireNamespace("Biobase", quietly=TRUE)){
  # only compute the scores
  s <- featureScore(res_dcis)
  summary(s)
  # compute the scores and characterize each metagene
  s <- extractFeatures(res_dcis)
  str(s)
}


# Function to extract and save genes per cluster
extract_genes_per_cluster <- function(nmf_fit, output_file) {
  # Extract the basis matrix
  basis_matrix <- basis(nmf_fit)
  
  # Find the cluster with the maximum weight for each gene
  gene_clusters <- apply(basis_matrix, 1, which.max)
  
  # Create a data frame with gene names and their corresponding clusters
  gene_cluster_df <- data.frame(
    Gene = rownames(basis_matrix),
    Cluster = gene_clusters
  )
  
  # Save the data frame to a CSV file
  write.csv(gene_cluster_df, file = output_file, row.names = FALSE, quote = FALSE)
  
  # Return the data frame for inspection
  return(gene_cluster_df)
}

# Example for rank 8
gene_cluster_mapping <- extract_genes_per_cluster(res_dcis$fit$`8`, "NMF_genes_per_cluster_rank8.csv")
print(head(gene_cluster_mapping))

############# JAS code below: #########################

# here we evealate the clusters quality to pick the best number of clusters
jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_rank_survey.jpg",width = 1724, height=1052)
plot(res_dcis) 
dev.off()

# same here, we evaluate how robust are the clusters
jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_2-10_hires.jpg",width = 6800, height=4202)
consensusmap(res_dcis)
dev.off()


# Loop through k2 to k10
for (k in 3:3) {
  # Extract sample names from the consensus matrix for the current k
  nmf_samples <- colnames(res_dcis$consensus[[as.character(k)]])
  
  # Filter the design file to keep only the matching samples
  design_filtered <- design[rownames(design) %in% nmf_samples, , drop=FALSE]
  
  # Reorder the design file rows to match the NMF consensus matrix order
  design_filtered <- design_filtered[nmf_samples, , drop=FALSE]
  design_filtered$group <- factor(design_filtered$group, levels = c("ADH", "DCIS"))
  
  # Ensure that design_filtered$group is used as annotations for columns
  annCol <- data.frame(Group = design_filtered$group)
  
  # Map the group to colors (blue for ADH and red for DCIS)
  annColors <- list(Group = c("ADH" = "blue", "DCIS" = "red"))
  
  # Now plot the consensus map with the annotations
  jpeg(file=paste0("/Users/x2a/helix/comet/nmf_2025_01/NMF_consensus_k", k, "_annotated.jpg"), width = 3000, height = 2000)
  consensusmap(res_dcis$consensus[[as.character(k)]], 
               annCol = annCol, 
               annColors = annColors, 
               labCol = 'sample', 
               main = paste('Cluster Stability for k =', k),
               sub = 'Consensus matrix with group annotations')
  dev.off()
}

for (k in 2:9) {
  # Load the required CSV with sample-cluster data for the chosen rank (e.g., k = 3)
  rank_sample_data <- read.csv(paste0("/Users/x2a/helix/comet/nmf_2025_01/NMF_sample_clusters/NMF_samples_per_cluster_rank", k, ".csv"))
  
  # Prepare column annotations based on the cluster labels from the data
  annCol <- data.frame(Cluster = factor(rank_sample_data$Cluster))
  rownames(annCol) <- rank_sample_data$Sample
  
  # Optional: Set colors for the clusters
  annColors <- list(Cluster = rainbow(length(unique(rank_sample_data$Cluster))))
  
  # Get the consensus matrix for the given rank (e.g., k = 3)
  consensus_matrix <- res_dcis$consensus[[as.character(k)]]  # Adjust for other k values as needed
  
  # Now plot the consensus map with the sample cluster annotations
  jpeg(file=paste0("/Users/x2a/helix/comet/nmf_2025_01/nmf_clustered_k", k, ".jpg"), width = 3000, height = 2000)
  consensusmap(consensus_matrix, 
               annCol = annCol, 
               annColors = annColors, 
               labCol = 'sample', 
               main = paste('Cluster Stability for k =', k),
               sub = 'Consensus matrix with sample cluster annotations')
  dev.off()
}




jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_k2_hires.jpg",width = 3000, height=2000)
consensusmap(res_dcis$consensus$`2`)
dev.off()

jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_k3_hires.jpg",width = 3000, height=2000)
consensusmap(res_dcis$consensus$`3`)
dev.off()

jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_k4_hires.jpg",width = 3000, height=2000)
consensusmap(res_dcis$consensus$`4`)
dev.off()
jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_k5_hires.jpg",width = 3000, height=2000)
consensusmap(res_dcis$consensus$`5`)
dev.off()
jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_k6_hires.jpg",width = 3000, height=2000)
consensusmap(res_dcis$consensus$`6`)
dev.off()
jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_k7_hires.jpg",width = 3000, height=2000)
consensusmap(res_dcis$consensus$`7`)
dev.off()
jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_k8_hires.jpg",width = 3000, height=2000)
consensusmap(res_dcis$consensus$`8`)
dev.off()
jpeg(file="/Users/x2a/helix/bph2023/nod1/NMF_bph_str85_nod1_100k_consensus_k9_hires.jpg",width = 3000, height=2000)
consensusmap(res_dcis$consensus$`9`)
dev.off()

for (k in 2:9) {
  write.csv(res_dcis$consensus$`8`, file=paste0("NMF_dcis_adh_n863_k",k,".csv"), quote=FALSE)
}

# here we get the cluster labels
k_2 = predict(res_dcis$fit$`2`,"consensus")
k_3 = predict(res_dcis$fit$`3`,"consensus")
k_4 = predict(res_dcis$fit$`4`,"consensus")
k_5 = predict(res_dcis$fit$`5`,"consensus")
k_6 = predict(res_dcis$fit$`6`,"consensus")
k_7 = predict(res_dcis$fit$`7`,"consensus")
k_8 = predict(res_dcis$fit$`8`,"consensus")
k_9 = predict(res_dcis$fit$`9`,"consensus")

#plot the silhuette plot (avoid solutions with negative values)
jpeg(file="NMF_bph_str85_nod1_100k_silhouette_hires.jpg",width = 10000, height=7500)
par(mfrow=c(2,4))
plot(silhouette(res_dcis$fit$`2`))
plot(silhouette(res_dcis$fit$`3`))
plot(silhouette(res_dcis$fit$`4`))
plot(silhouette(res_dcis$fit$`5`))
plot(silhouette(res_dcis$fit$`6`))
plot(silhouette(res_dcis$fit$`7`))
plot(silhouette(res_dcis$fit$`8`))
plot(silhouette(res_dcis$fit$`9`))
dev.off()



# Function to extract genes, clusters, and weights for a given rank
extract_genes_all_ranks <- function(nmf_result, ranks, output_dir) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each rank and extract information
  all_ranks_data <- data.frame()
  for (rank in ranks) {
    cat("Processing rank:", rank, "\n")
    
    # Extract the basis matrix for the current rank
    nmf_fit <- nmf_result$fit[[as.character(rank)]]
    basis_matrix <- basis(nmf_fit)
    
    # Find the cluster with the maximum weight for each gene
    max_weights <- apply(basis_matrix, 1, max)
    gene_clusters <- apply(basis_matrix, 1, which.max)
    
    # Create a data frame for the current rank
    rank_data <- data.frame(
      Gene = rownames(basis_matrix),
      Cluster = gene_clusters,
      Max_Weight = max_weights,
      Rank = rank
    )
    
    # Save the data for the current rank to a CSV file
    output_file <- file.path(output_dir, paste0("NMF_genes_weights_rank", rank, ".csv"))
    write.csv(rank_data, file = output_file, row.names = FALSE, quote = FALSE)
    
    # Combine with all ranks data
    all_ranks_data <- rbind(all_ranks_data, rank_data)
  }
  
  # Save combined data for all ranks
  combined_file <- file.path(output_dir, "NMF_genes_weights_all_ranks.csv")
  write.csv(all_ranks_data, file = combined_file, row.names = FALSE, quote = FALSE)
  
  return(all_ranks_data)
}

# Specify ranks and output directory
ranks <- 2:9
output_dir <- "./NMF_gene_weights_TEST"

# Extract and save data for all ranks
all_ranks_data <- extract_genes_all_ranks(res_dcis, ranks, output_dir)
print(head(all_ranks_data))





# Function to extract samples and clusters for all ranks
extract_samples_all_ranks <- function(nmf_result, ranks, output_dir) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each rank and extract sample-cluster associations
  all_ranks_sample_data <- data.frame()
  for (rank in ranks) {
    cat("Processing rank:", rank, "\n")
    
    # Extract the coefficient matrix for the current rank
    nmf_fit <- nmf_result$fit[[as.character(rank)]]
    coef_matrix <- coef(nmf_fit)
    
    # Find the cluster with the maximum coefficient for each sample
    sample_clusters <- apply(coef_matrix, 2, which.max)
    
    # Create a data frame for the current rank
    rank_sample_data <- data.frame(
      Sample = colnames(coef_matrix),
      Cluster = sample_clusters,
      Max_Coefficient = apply(coef_matrix, 2, max),
      Rank = rank
    )
    
    # Save the data for the current rank to a CSV file
    output_file <- file.path(output_dir, paste0("NMF_samples_per_cluster_rank", rank, ".csv"))
    write.csv(rank_sample_data, file = output_file, row.names = FALSE, quote = FALSE)
    
    # Combine with all ranks data
    all_ranks_sample_data <- rbind(all_ranks_sample_data, rank_sample_data)
  }
  
  # Save combined data for all ranks
  combined_file <- file.path(output_dir, "NMF_samples_per_cluster_all_ranks.csv")
  write.csv(all_ranks_sample_data, file = combined_file, row.names = FALSE, quote = FALSE)
  
  return(all_ranks_sample_data)
}

# Specify ranks and output directory
ranks <- 2:9
output_dir <- "./NMF_sample_clusters"

# Extract and save data for all ranks
all_ranks_sample_data <- extract_samples_all_ranks(res_dcis, ranks, output_dir)
print(head(all_ranks_sample_data))
