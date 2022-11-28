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
# This code loads a hapmat file but loses the column names. Turns rs# into rs and uses it as row names
cross <- read.table(
  "/home/drt83172/Documents/Tall_fescue/Kmer_Genotyping/Kmer_genotyping/Genetic_Mapping/Data/Genotype_Files/312x314_hapmap.txt", 
  quote="", sep ="\t", row.names = 1, fill = TRUE, comment.char ="", header = TRUE )






