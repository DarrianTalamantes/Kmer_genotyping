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
# # Make standard deviation values for each progeny across all parents
Score_Means <- as.data.frame(rowMeans(Scores))
Scores_std <- apply(Scores, 1, sd, na.rm=TRUE)

# # Make column of mean values for each parent across all progeny. 
# # Make standard deviation values for each parent across all progeny.
Score_Means2 <- as.data.frame(colMeans(Scores))
Scores_std2 <- apply(Scores, 2, sd, na.rm=TRUE)

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

# # Calculating Z-scores based on 1 parent compared to all progeny
z_scores <- Scores
for (i in 1:ncol(Scores)){
  for (j in 1:nrow(Scores)){
    z_scores[j,i] <- (Scores[j,i] - Score_Means2[i,1]) / Scores_std2[i]  
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

z_scores_key_t <- as.data.frame(t(z_scores_key))
z_scores_key_t_m <- t(z_scores_key)

# # Kmeans clustering
x = 301
Predicted_Parents <- z_scores
for (i in 1:ncol(z_scores)){
  Predicted_Parent <- subset(z_scores[,i, drop = FALSE])
  clusters <- kmeans(z_scores[,i],2)
  Predicted_Parents[i]=clusters$cluster
  if (x == 308 | x == 310 | x == 316){
    x = x+2}
  else
    x = x + 1
}

poop <- kmeans(cars[,4],2)
cars <- cars %>% mutate(cluster = poop$cluster)





########################## Graphs #################################
# # Density plots that shows known mothers. each graph is parents by progeny

i = 0
x = 301
for (i in 1:17){
  Parent = colnames(z_scores_key[i])
  print(x)
  KnownParent = subset(z_scores_key, Mother == x)
  print(ggplot(z_scores_key, aes_(x=as.name(Parent))) + geom_histogram(binwidth=.01)  + geom_histogram(data = KnownParent, fill = "green1", binwidth=.01)+ theme_bw())
  if (x == 308 | x == 310 | x == 316){
    x = x+2}
  else {
    x = x + 1}
}
# # Density plot that shows known mothers each graph is progeny by parents
number_o_progeny = nrow(z_scores_key)
i = 1
for (i in 1:number_o_progeny){
  Maternal <- z_scores_key_t_m[18,i]
  Progeny = colnames(z_scores_key_t[i])
  working_scores <- subset(z_scores_key_t, select = c(paste0(Progeny)))
  working_scores <- working_scores[-c(18),, drop = FALSE ]
  for (j in 1:n_parents){
    if (as.integer(Parent_names[j]) == Maternal){
      Parent_name <- Parent_names[j] 
      Progeny_name <- Progeny_names[i]
      small_df <- as.data.frame(working_scores[j,1])
      colnames(small_df) <- c(Progeny_name)
      rownames(small_df) <- c(Parent_name)
     }
    }
  print(ggplot(working_scores, aes_(x=as.name(Progeny))) + geom_histogram(binwidth=.2) + geom_histogram(data =) + 
          geom_histogram(data = small_df, fill = "darkblue", binwidth=.2)+ theme_bw())
}

# # Histogram showing a divide between possible progeny and random stuff
hist(Scores[,1])
sum(Scores[,1]<1400)
342/3000

# # NMDS of parents with progeny as points
for (i in 1:ncol(z_scores)){
  small_df <- z_scores[i]
  nmds = metaMDS(small_df, distance = "bray")
}


cars <- mtcars
poop <- kmeans(cars[,4],2)
cars <- cars %>% mutate(cluster = poop$cluster)


