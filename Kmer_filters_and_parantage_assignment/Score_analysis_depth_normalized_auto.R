library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)
library(car)


# Argument inputs
Args <- commandArgs(trailingOnly=TRUE)
# input variables will need to be a number for the seed and all input files. 
Scores <- read.table (Args[2], sep = ",", header = TRUE, row.names = 1)
half_key_parents <- read.table (Args[3], sep = "\t", header = TRUE)
half_key_progeny <- read.table (Args[4], sep = "\t", header = TRUE)
file_to_feild_key <- read.table(Args[5], sep = ",", header = TRUE)
depth_data <- read.table("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/sample_depth.txt", sep = ",",header = TRUE, row.names = 1 )

# # Organizing data 
key <- cbind(half_key_progeny,half_key_parents)
key <- tibble::column_to_rownames(key, "FileName")
colnames(key)[1] <- "Known_parent"
colnames(file_to_feild_key)[1] <- "Progeny_feild"
file_to_feild_key <- tibble::column_to_rownames(file_to_feild_key, "FileName")

# # Combining depth data with Score data
depth_codes = select(depth_data, c("Customer_Code", "Depth")) 
rownames(depth_codes) <- depth_codes[,1]
depth_codes <- select(depth_codes, -c ("Customer_Code"))

Scores_feild_names <- merge(x=Scores, y=file_to_feild_key, by="row.names")
Scores_feild_names <- merge(x=Scores_feild_names, y=key, by.x ="Row.names", by.y = "row.names")
Scores_feild_names <- select(Scores_feild_names, -c("Row.names"))
Scores_feild_names <- tibble::column_to_rownames(Scores_feild_names, "Progeny_feild")

Scores_depth <- merge(x=Scores_feild_names, y=depth_codes, by="row.names")
Scores_depth <- tibble::column_to_rownames(Scores_depth, "Row.names")

x = 301
for (i in 1:17){
  colnames(Scores_depth)[i] <- x
  if (x == 308 | x == 310 | x == 316){
    x = x+2}
  else
    x = x + 1
}

# # Using depth to normalize the scores
Scores_depth_norm <- Scores_depth[1:17]/Scores_depth$Depth
Known <- Scores_depth[18:19]
Scores_depth_norm <-  merge(x=Scores_depth_norm, y=Known, by="row.names")
Scores_depth_norm <- tibble::column_to_rownames(Scores_depth_norm, "Row.names")
Scores_depth_norm2 <- select(Scores_depth_norm, -c ("Known_parent","Depth"))

# # Kmeans clustering

x = 301
All_centers <- data.frame(matrix(0, ncol = 17, nrow = 4))
colnames(All_centers) <- colnames(Scores_depth_norm2)
Predicted_Parents <- Scores_depth_norm2
for (i in 1:ncol(Scores_depth_norm2)){
  Predicted_Parent <- subset(Scores_depth_norm2[,i, drop = FALSE])
  set.seed(Args[1])
  clusters <- kmeans(Scores_depth_norm2[,i],4)
  Predicted_Parents[i]=clusters$cluster
  All_centers[i] <- clusters$centers
  if (x == 308 | x == 310 | x == 316){
    x = x+2}
  else
    x = x + 1
}
write.table(Predicted_Parents, file = "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Predicted_Parents.txt")
write.table(All_centers, file = "/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/All_centers.txt")




