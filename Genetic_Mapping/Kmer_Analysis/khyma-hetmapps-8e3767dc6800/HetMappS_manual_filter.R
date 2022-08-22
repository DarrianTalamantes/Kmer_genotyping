
######################################################################################
#
#     HetMappS_manual_filter <- Provides tables and figures to help determine spurius markers
#                         Removes markers designated by user.
#			  Determine error rate
#			  Estimate genetic distances 
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
library(snow)


#### define working directory, HetMappS_functions.R should be located there
#setwd(       )

#### Load functions 
source("HetMappS_functions.R")


## This stage of the pipeline requires manual input. 
## Basically, look at rf plots, rf summary and identify markers to manually drop
## then, re order the LG
#### this step may take long time


### Variables to define

# the size of increments in error rate value for calculating genotyping error
k.error.rate.increment <- 0.0025
# the upper limit range for calculating genotyping error. multiple of k.error.rate.increment
k.error.rate.limit <- 0.0125

# number of clusters to use for parallel processing
k.cluster <- 1




load("Object.data3.rf.R")
data.man <- data3.rf



###################### Here introduce markers to drop##############
####
### by creating a vector with their names
markers.manual.filter <- "NA"
### by loading a one-columns csv file with marker's names
markers.manual.filter <- read.csv("drop_man.csv",header = FALSE, colClasses = "character")

	
	
####################################################################







########## Generate tables and plots to aid the identification of 
########## spurious markers


## Generate a summary of map distances, mean rf and mean lod for each marker
## data_order objects need to have rf calculated (use est.rf)
data.manual.summary <- rfSummary(data.man)
write.table(data.manual.summary, file = "summary_for_manual_filter.csv", sep = ",", col.names = NA, row.names = TRUE)


### create plots and summary data to aid manual drop of markers
## generate a file with both, plot of the rf and genoplot for each chromosome
pdf(file = "plots_for_manual_filter.pdf", height = 7.5, width = 10)
for(n in 1:nchr(data.man)){
par(mfrow = c(1,2))
plot.rf(data.man,chr =n)
plotGeno(data.man,chr = n, ind = c(1:25), min.sep = 0)
}
dev.off()

##########################################################################


## drop the desired markers
data.man <- drop.markers(data.man, markers.manual.filter)

## determine the error rate of this cross object
data.man.error.table <- ErrorRate(data.man,err.end = k.error.rate.limit, cluster = k.cluster )
data.man.error <- data.man.error.table[which.max(data.man.error.table[,2]),1]
cat("\n Error rate of cross is ", data.man.error, "\n")
save(data.man.error.table, file = "Table.data.man.error.table.R")



## re-estimate marker distances for all linkage groups
data.map <- NULL
data.map <- est.map(data.man , error.prob = data.man.error, sex.sp = FALSE,  offset = 0, n.cluster = k.cluster)
data.man.out <- replace.map(data.man, data.map)
#save the new cross file
save(data.man.out, file = "Object.data.man.out.R")









