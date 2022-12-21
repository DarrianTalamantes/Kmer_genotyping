# This scrpit connects to match_all_progeny_to_parents.sh
# this just takes the inputs from that script and merges them and outputs the merged file

Args <- commandArgs(trailingOnly=TRUE)
# input variables will need to be a number for the seed and all input files. 
file_to_feild_key <- read.table(Args[1], sep = ",", header = TRUE) #key of file names and feild codes
codes <- read.table(Args[2], sep = ",", header = FALSE) #codes I am looking at currently
save_file <- Args[3] # file  of codes I am using attached to their file names

merged <- merge(file_to_feild_key, codes, by.x = 1, by.y = 1)
merged <- subset(merged,select=c("Customer_Code","FileName"))
write.table(merged, file = save_file, col.names = F, row.names = F)














