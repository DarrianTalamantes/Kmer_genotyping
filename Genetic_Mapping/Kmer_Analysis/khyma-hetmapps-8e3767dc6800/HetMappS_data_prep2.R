
######################################################################################
#
#    HetMappS_data_prep2.R <- Determine error rate, flip desired linkage groups, 
#					re-calculates genetic distances and apply SlideDroponemarker
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


#### Load data object
load("Object.data_order_out.R")
data_prep <- data_order_out


### Variables to define

# the size of increments in error rate value for calculating genotyping error
k.error.rate.increment <- 0.0025
# the upper limit range for calculating genotyping error. multiple of k.error.rate.increment 
k.error.rate.limit <- 0.03
# number of clusters to use for parallel processing
k.cluster <- 1
# a vector with linkage groups to invert
flip.vector <- "NA"



## determine the error rate of this cross object
data2.error.table <- ErrorRate(data_prep,err.end = k.error.rate.limit, cluster = k.cluster )
data2.error <- data2.error.table[which.max(data2.error.table[,2]),1]
cat("Error rate of cross is ", data2.error, "\n")
save(data2.error.table, file = "Table.data2.error.table.R")

##### Flip some linkage groups and recalculate distances whole LG set
data2 <- FlipLG(data_prep,flip.vector, error.rate = data2.error, cluster = k.cluster)
save(data2,file = "Object.data2.R")



## apply droponemarker in sliding window
data2.drop1marker.table <- SlideDroponemarker(data2, error.rate = data2.error)
save(data2.drop1marker.table, file = "Table.data2.drop1marker.table.R")

