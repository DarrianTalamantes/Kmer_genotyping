#!/usr/bin/env Rscript
library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)
library(car)

# Argument inputs
args = commandArgs(trailingOnly=TRUE)


# input variables will need to be a number for the seed and all input files. 
Scores <- read.table ("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Score_table.csv", sep = ",", header = TRUE, row.names = 1)
half_key_progeny <- read.table ("/home/drt83172/Documents/Tall_fescue/half_key_parents.txt", sep = "\t", header = TRUE)
half_key_parents <- read.table ("/home/drt83172/Documents/Tall_fescue/half_key_progeny.txt", sep = "\t", header = TRUE)
file_to_feild_key <- read.table("/home/drt83172/Documents/Tall_fescue/progeny_key.csv", sep = ",", header = TRUE)

# # Organizing data 
key <- cbind(half_key_progeny,half_key_parents)
colnames(file_to_feild_key)[1] <- "Progeny_feild"
file_to_feild_key <- tibble::column_to_rownames(file_to_feild_key, "FileName")

# # Make column of mean values for each progeny across all parents. 
# # Make standard deviation values for each progeny across all parents
Score_Means <- as.data.frame(rowMeans(Scores))
Scores_std <- apply(Scores, 1, sd, na.rm=TRUE)

# # Calculating Z-scores based on 1 progeny compared to all parents
z_scores <- Scores
for (i in 1:nrow(Scores)){
  for (j in 1:ncol(Scores)){
    z_scores[i,j] <- (Scores[i,j] - Score_Means[i,1]) / Scores_std[i]  
  }
}
x = 301
for (i in 1:17){
  colnames(z_scores)[i] <- x
  if (x == 308 | x == 310 | x == 316){
    x = x+2}
  else
    x = x + 1
}

# # Making the prgoeny names the same as the ones out in the field
z_scores <- merge(x=z_scores, y=file_to_feild_key, by="row.names")
z_scores <- tibble::column_to_rownames(z_scores, "Progeny_feild")
z_scores <- select(z_scores, -c("Row.names"))  

# # Kmeans clustering

x = 301
All_centers <- data.frame(matrix(0, ncol = 17, nrow = 4))
colnames(All_centers) <- colnames(z_scores)
Predicted_Parents <- z_scores
for (i in 1:ncol(z_scores)){
  Predicted_Parent <- subset(z_scores[,i, drop = FALSE])
  set.seed(10)
  clusters <- kmeans(z_scores[,i],4)
  Predicted_Parents[i]=clusters$cluster
  All_centers[i] <- clusters$centers
  if (x == 308 | x == 310 | x == 316){
    x = x+2}
  else
    x = x + 1
}

write.table(Predicted_Parents, file = "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Predicted_Parents.txt")
write.table(All_centers, file = "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/All_centers.txt")














