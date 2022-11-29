#!/usr/bin/env Rscript

# Author: Darrian Talamantes
# Email: darrianrtalamantes6@gmail.com

# Purpose: This script will grab a hapmat file and create linkage groups from the 
# kmers that it provides

# Load libraries
library(reshape2)
library(tidyverse)
library(vegan)


############## Script Begin ####################
# Loading in data and setting variables
# This code reads a hapmap file, gets rid of metadata, transposes it
cross <- read.table(
  "/home/drt83172/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Genotype_Files/301x313_hapmap.txt", 
  quote="", sep ="\t", row.names = 1, fill = TRUE, comment.char ="", header = TRUE )
cross <- subset(cross, select = -c(1:10))
cross.2 <- data.frame(t(cross))


# Trying to do clustering with my actual data
cross_matrix <- as.matrix(cross.2)
for (row in 1:nrow(cross_matrix)){
  for (col in 1:ncol(cross_matrix)){
    if (cross_matrix[row,col] == "A"){
      cross_matrix[row,col] = 1
    }
    else{
      cross_matrix[row,col] = 0
    }
  }
}


cross_distances <- dist(cross_matrix)
head(cross_distances)
cross_clsuters <- hclust(cross_distances)
plot(cross_clsuters, xaxt='n', ann=FALSE)

ggtree(cross_clsuters)


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

