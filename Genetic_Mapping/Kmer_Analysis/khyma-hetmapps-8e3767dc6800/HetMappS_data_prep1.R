

#############################################################################################
#
#     HetMappS_data_prep1.R <- prepare data for map analysis by:
#					loading output from phased2rqtl
#					dropping markers and individuals with too much missing data
#					dropping individuals with too many cross over 
#
#    Copyright (C) 2015  Paola Barba (plb74@cornell.edu)
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
###############################################################################################

### Packages requiered:
library(qtl)

#### define working directory, phased2rqtl output and HetMappS_functions.R should be located there
#setwd(       )

#### Load functions 
source("HetMappS_functions.R")

### Variables to define

pheno.file <- 
geno.file <-  
k.prop.markers.typed = 0.5
k.prop.individuals.typed = 0.5
k.prop.mean.cross.over = 2





## read the phased2rqtl output
### load the data, create the cross object
data <- read.cross(format = "csvsr", phefile = pheno.file, genfile = geno.file ,genotypes = c("A", "B"))


### filter data
data1 <- PrepCross(data,prop.markers.typed = k.prop.markers.typed, prop.individuals.typed = k.prop.individuals.typed, prop.mean.cross.over = k.prop.mean.cross.over)
save(data1, file = "Object.data1.R")

### calculate recombination fractions
data1.rf <- est.rf(data1)
save(data1.rf, file = "Object.data1.rf.R")

### plot recombination fraction
plot.rf(data1.rf)

pdf(file = "recombination_fraction_per_chr_data_prep_1.pdf", width = 10, height = 7.5)
for(i in 1:nchr(data1.rf)) plot.rf(data1.rf, chr = i)
dev.off()





