
######################################################################################
#
#    HetMappS_auto_filter <- Removes markers, recalculates genetic distances and recombination fractions
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
## Basically, look at the distributions of the lod score for a desired Ldiff and 
## determine the LOD threshold


### Variables to define
# the lod threshold for dropping markers according to SlideDroponemarker results
lodt <- -15
# Length difference threshold for dropping markers acording to SlideDroponemarker results
k.limLdiff <- 2
# the size of increments in error rate value for calculating genotyping error
k.error.rate.increment <- 0.0025
# the upper limit range for calculating genotyping error. multiple of k.error.rate.increment
k.error.rate.limit <- 0.02
# number of clusters to use for parallel processing
k.cluster <- 1
# error rate from data_prep2
data2.error <- 0.01



load("Object.data2.R")
load("Table.data2.drop1marker.table.R")


# filter cross object based on droponemarker slide table
data2.filter <- FilterByLODSlidingDroponemarker(data2, data2.drop1marker.table, limLOD = lodt, limLdiff = k.limLdiff)
 

## determine the error rate of this cross object
data3.error.table <- ErrorRate(data2.filter,err.end = k.error.rate.limit, cluster = k.cluster,err.increment = k.error.rate.increment)
data3.error <- data3.error.table[which.max(data3.error.table[,2]),1]
cat("Error rate of cross is ", data3.error, "\n")
save(data3.error.table, file = "Table.data3.error.table.R")
### alternatively, the user can use the error rate determined from data prep 2
## data3.error <- data2.error

## re-estimate marker distances for all linkage groups
data.map <- NULL
data.map <- est.map(data2.filter , error.prob = data3.error, sex.sp = FALSE,  offset = 0, n.cluster = k.cluster)
data3 <- replace.map(data2.filter, data.map)
#save the new cross file
save(data3, file = "Object.data3.R")




#### calculate recombination fractions and lod 
data3.rf <- est.rf(data3)
save(data3.rf, file = "Object.data3.rf.R")





## Generate a summary of map distances, mean rf and mean lod for each marker
data3.rf.summary <- rfSummary(data3.rf)
save(data3.rf.summary, file = "Table.data3.rf.summary.R")
#write.table(data3.rf.summary, file = "summary_for_manual_drop.csv", sep = ",", col.names = NA, row.names = TRUE)

## plot the map
#plot.map(data3)





