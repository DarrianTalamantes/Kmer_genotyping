library(reshape2)
library(tidyverse)
library(vegan)
library(Rfast)
library(car)

Scores <- read.table ("/home/drt83172/Documents/Tall_fescue/Usefull_Kmers/R_Files/Score_table.csv", sep = ",", header = TRUE, row.names = 1)
half_key_progeny <- read.table ("/home/drt83172/Documents/Tall_fescue/half_key_parents.txt", sep = "\t", header = TRUE)
half_key_parents <- read.table ("/home/drt83172/Documents/Tall_fescue/half_key_progeny.txt", sep = "\t", header = TRUE)

# # Organizing data 
key <- cbind(half_key_progeny,half_key_parents)



# # Make column of mean values for each progeny across all parents. 
Score_Means <- as.data.frame(rowMeans(Scores))
Scores_std <- apply(Scores, 1, sd, na.rm=TRUE)

# # Calculating Z-scores
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


# # adding key to data
key <- tibble::column_to_rownames(key, "FileName")
z_scores_key <- merge(x=z_scores, y=key, by="row.names")
colnames(z_scores_key)[19] <- "Mother"
z_scores_key <- tibble::column_to_rownames(z_scores_key, "Row.names")


########################## Graphs #################################
# # Density plot that shows known mothers

i = 0
x = 301
for (i in 1:17){
  Parent = colnames(z_scores_key[i])
  KnownParent = subset(z_scores_key, Mother == x)
  print(ggplot(z_scores_key, aes_(x=as.name(Parent))) + geom_histogram(binwidth=.01)  + geom_histogram(data = KnownParent, fill = "green1", binwidth=.01)+ theme_bw())
  if (x == 308 | x == 310 | x == 316){
    x = x+2}
}




