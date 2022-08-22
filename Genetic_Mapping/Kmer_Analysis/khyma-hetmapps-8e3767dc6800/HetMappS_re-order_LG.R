
######################################################################################
#
#    HetMappS_re-order_LG <- Removes markers designated by user and re-order markers in selected LG
#
#    Copyright (C) 2014  Paola Barba (plb74@cornell.edu)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
######################################################################################

### Packages required:
library(qtl)



#### define working directory, HetMappS_functions.R should be located there
#setwd(       )

#### Load functions 
source("HetMappS_functions.R")


## This stage of the pipeline reqieres manual input. 
## Basically, look at rf plots, rf summary and identify markers to manually drop
## then, re order the LG
#### this step may take long time



### Variables to define

## load data from data_prep_1
load("Object.data1.R")
data_order <- data1
## chromosome/linkage group to re-order
chr2order <- "NA"



###################### Here introduce markers to drop##############
####
### by name
markers2drop <- "NA"
## by a one column csv file with the markers to drop
markers2drop <- read.csv("to_drop.csv",header = FALSE, colClasses = "character")

####################################################################







########## Generate tables and plots to aid the identification of 
########## spurious markers


## Generate a summary of map distances, mean rf and mean lod for each marker
## data_order objects need to have rf calculated (use est.rf)
data_order.summary <- rfSummary(est.rf(data_order))
write.table(data_order.summary, file = "summary_for_order_LG.csv", sep = ",", col.names = NA, row.names = TRUE)


### create plots and summary data to aid manual drop of markers
## generate a file with both, plot of the rf and genoplot for each chromosome
pdf(file = "plots_for_order_LG.pdf", height = 7.5, width = 10)
for(n in 1:nchr(data_order)){
par(mfrow = c(1,2))
plot.rf(data_order,chr =n)
plotGeno(data_order,chr = n, ind = c(1:25), min.sep = 0)
}
dev.off()

##########################################################################


## find and drop duplicated markers
duplicated <- findDupMarkers(data_order,chr2order,exact.only = FALSE)
data_order <- drop.markers(data_order, unlist(duplicated))

## drop the desired markers
data_order <- drop.markers(data_order, markers2drop)

##########################################################################
## re calculate marker order for selected chr
data_order_out <- orderMarkers(data_order, chr = chr2order, error.prob=0.02, use.ripple = FALSE, sex.sp=FALSE)
save(data_order_out, file = "Object.data_order_out.R")

#########################################################################
## generate a file to evaluate the new order

## Estimate recombination fractions for reordered chr
data_order_subset.rf <- data_order_out[chr = chr2order]
save(data_order_subset.rf , file = "Object.data_order_subset.rf.R") 

## generate a file with both, plot of the rf and genoplot for ordered chr

if(nchr(data_order_subset.rf)==1){
pdf(file = "plots_ordered_LG.pdf", height = 7.5, width = 10)
	plot.rf(data_order_subset.rf,chr =chr2order)
	plotGeno(data_order_subset.rf,chr = chr2order, ind = c(1:25), min.sep = 0)
dev.off()
} else {
pdf(file = "plots_ordered_LG.pdf", height = 7.5, width = 10)
	for(n in 1:nchr(data_order_subset.rf)){
		par(mfrow = c(1,2))
		plot.rf(data_order_subset.rf,chr =chr2order[n])
		plotGeno(data_order_subset.rf,chr = chr2order[n], ind = c(1:25), min.sep = 0)
	}
dev.off()
}












