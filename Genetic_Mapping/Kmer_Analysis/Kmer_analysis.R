#!/usr/bin/env Rscript
######################### Checking for args ###################################
args = commandArgs(trailingOnly=TRUE)
if (length(args)<=2) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)}
######################### Loading Libraries ####################################
#library(ggplot2)
#library(dplyr)
library(reshape2)
library(tidyverse)
library(vegan)
######################### Importing files into big data table###################
parent_files <- list.files(include.dirs =TRUE, args[1])
# parent_files <- list.files(include.dirs =TRUE, "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/Parents")
number_of_parents <- length(parent_files)
# N=8 #activate for testing
for (N in 1:number_of_parents){
  parent = paste(args[1], parent_files[N], sep = "/")
  # parent = paste("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/Parents", parent_files[N], sep = "/")
  name = parent_files[N]
  if (N == 1 ){
    print(parent)
    all_parents <- read.table(parent, sep = "\t" )
    all_parents <- all_parents %>% rename(KMER = V1)
    all_parents <- all_parents %>% rename(!!name := V2)
    all_parents <- all_parents %>% remove_rownames %>% column_to_rownames(var="KMER")
    all_parents <- as.data.frame(t(all_parents))
    }
  else{
    print(parent)
    temp_file <- read.table(parent, sep = "\t" )
    temp_file <- temp_file %>% rename(KMER = V1)
    temp_file <- temp_file %>% rename(!!name := V2)
    temp_file <- temp_file %>% remove_rownames %>% column_to_rownames(var="KMER")
    temp_file <- as.data.frame(t(temp_file))
    all_parents <- rbind(all_parents,temp_file)
  }
}
# write.table(all_parents, file = "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/R_parents.txt")
write.table(all_parents, file = args[3])

progeny_files <- list.files(include.dirs =TRUE, args[2])
# progeny_files <- list.files(include.dirs =TRUE, "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/Progeny")
number_of_progeny <- length(progeny_files)
# N=8 #activate for testing
for (N in 1:number_of_progeny){
  progeny = paste(args[2], progeny_files[N], sep = "/")
  # progeny = paste("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/Progeny", progeny_files[N], sep = "/")
  name = progeny_files[N]
  if (N == 1 ){
    print(progeny)
    all_progeny <- read.table(progeny, sep = "\t" )
    all_progeny <- all_progeny %>% rename(KMER = V1)
    all_progeny <- all_progeny %>% rename(!!name := V2)
    all_progeny <- all_progeny %>% remove_rownames %>% column_to_rownames(var="KMER")
    all_progeny <- as.data.frame(t(all_progeny))
  }
  else{
    print(progeny)
    temp_file <- read.table(progeny, sep = "\t" )
    temp_file <- temp_file %>% rename(KMER = V1)
    temp_file <- temp_file %>% rename(!!name := V2)
    temp_file <- temp_file %>% remove_rownames %>% column_to_rownames(var="KMER")
    temp_file <- as.data.frame(t(temp_file))
    all_progeny <- rbind(all_progeny,temp_file)
  }
}

# write.table(all_progeny, file = "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/R_progeny.txt")
write.table(all_progeny, file = args[4])


