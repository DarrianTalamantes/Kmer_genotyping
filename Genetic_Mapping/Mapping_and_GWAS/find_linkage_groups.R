#!/usr/bin/env Rscript

# Author: Darrian Talamantes
# Email: darrianrtalamantes6@gmail.com

# Purpose: This script will grab a hapmat file and create linkage groups from the 
# kmers that it provides

# Load libraries
library(reshape2)
library(tidyverse)
library(vegan)
library(ggtree)

############## Script Begin ####################
# Loading in data and setting variables
# This code reads a hapmap file, gets rid of metadata, transposes it
cross <- read.table(
  "~/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Genotype_Files/301x313_hapmap.txt", 
  quote="", sep ="\t", row.names = 1, fill = TRUE, comment.char ="", header = TRUE )
cross <- subset(cross, select = -c(1:10))
cross.2 <- data.frame(t(cross))


# Trying to do clustering with my actual data
cross_matrix <- as.matrix(cross)

cross_matrix_subset <- as.matrix(cross[1:1000,])

cross_matrix_t <- as.matrix(cross.2)

# Function makes the A's and C's into numerical data
letter2Number <- function(matrix){
  matrix2 <- matrix
  for (row in 1:nrow(matrix)){
    for (col in 1:ncol(matrix)){
      if (matrix[row,col] == "A"){
        matrix2[row,col] = 1
      }
      else{
        matrix2[row,col] = 0
      }
    }
  } 
  return(matrix2)
}

# cross_matrix_num <- letter2Number(cross_matrix)
cross_matrix_t_num <- letter2Number(cross_matrix_t)
cross_matrix_subset_num <- letter2Number(cross_matrix_subset)

# Making distance matrix
cross_distances <- dist(cross_matrix_t_num)
kmer_distances <- dist(cross_matrix_subset)
kmer_distances_subset <- dist(cross_matrix_subset_num)

kmer_clusters_subset <- hclust(kmer_distances_subset)
kmer_clusters <- hclust(kmer_distances)
cross_clsuters <- hclust(cross_distances)

kmer_clusters_subset$dist.method
# tree graphs
ggtree(cross_clsuters) +
  geom_tiplab()
ggtree(kmer_clusters_subset) 


# Examples
# NMDS example
set.seed(2)
community_matrix=matrix(
  sample(1:100,300,replace=T),nrow=10,
  dimnames=list(paste("community",1:10,sep=""),paste("sp",1:30,sep="")))

example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                     k=2) # The number of reduced dimensions
stressplot(example_NMDS)

# Hclust example
head(iris)
distances <- dist(iris[, 3:4])
clusters <- hclust(distances)
plot(clusters)

#ggtree examples
tree <- rtree(50)
ggtree(tree)
ggtree(tree, layout="circular")
ggtree(tree, branch.length='none')

# example vegclust

vegclust




