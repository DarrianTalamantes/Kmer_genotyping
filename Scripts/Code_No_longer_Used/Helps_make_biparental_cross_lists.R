library(tidyverse)
library(vegan)


key <- read.table("/home/drt83172/Documents/Tall_fescue/progeny_key.csv", sep = ",", header = TRUE)
progeny_List <- read.table("/home/drt83172/Documents/Tall_fescue/Plant_Info/Parent314_312.csv", sep = ",", header = TRUE)
progeny <- merge(key,progeny_List, by.x = "Customer_Code", by.y = "X")
progeny <- subset(progeny, select = -c(X0,X1))
progeny
write.csv(progeny,"/home/drt83172/Documents/Tall_fescue/Kmer_analyses/Cross_Lists/Biparental/314x312_Progeny", row.names = FALSE)
