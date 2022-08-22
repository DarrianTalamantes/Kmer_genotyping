
#############################################################################################
#
#     HetMappS_functions <- compile functions for curation of HetMappS 
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
###############################################################################################



### Packages requiered for these functions:

	 library(qtl)
	 library(snow)
	
	
####   Functions in this script ####
#   PrepCross
#   FlipLG
#   ErrorRate
#   SlideDroponemarker
#   FilterByLODSlidingDroponemarker
#   rfSummary

	
	
	
	
	
##############################################################################################################
# PrepCross:	Performs standard preparation of cross data
#				for linkage map analysis, as indicated in 
#				according to http://www.rqtl.org/tutorials/geneticmaps.pdf 
#
#  It will drop:
#		- markers with proportion of information below threshold (default = 0.6)
#		- Individuals with proportion of information below threshold  (default = 0.6)
#		- Individuals with too many cross over, defined as proportion
#       over the sample mean (default = 2 )
###############################################################################################################
# Input:
#	cross2prep: a object from read.cross (requiered)
# 
# Options:
#	prop.markers.typed: The minimum proportion of typed markers (default 0.6, from 0 to 1)
#	prop.individuals.typed: The minimum proportion of markers within individuas (default 0.6, from 0 to 1)
#	prop.mean.cross.over: threshold for proportion of cross over (default 2, ie, twice the mean of the sample)  
# 
# Output: 
#	same cross object with markers and individuals removed
#################################################################################################################


PrepCross <- function(cross2prep, prop.markers.typed = 0.6, prop.individuals.typed = 0.6, prop.mean.cross.over = 2 ) {
 	# drop sites with few genotypes
	nt.bymar <- ntyped(cross2prep, "mar")
	todrop <- names(nt.bymar[nt.bymar < nind(cross2prep)*prop.markers.typed])
	cross2prep <- drop.markers(cross2prep, todrop)
	
	# drop individuals with few genotypes
	cross2prep <- subset(cross2prep, ind=(ntyped(cross2prep)>sum(nmar(cross2prep))*prop.individuals.typed))
	
	# drop individuals with too many cross over
	cross2prep <- subset(cross2prep, ind=(countXO(cross2prep) < mean(countXO(cross2prep))*prop.mean.cross.over))
	
	return(cross2prep)
}	
	
	
	
		
	
##################################################################################################
# FlipLG:		Invert the order of specific linkage groups in a genetic 
#				map and re-estimate the genetic distances for the whole cross object
# 
##################################################################################################
# Input:
#	cross.data: a object from read.cross (requiered) 
# 
# Options:
#	LG2flip: a vector contaning the name of the LG flip
#	error.rate: error rate for the genotyping error. Default is 0.01. 
#	cluster: Number of cluster to run paralel calculations. Default is 1. requieres package snow()
#	m.function: c("haldane","kosambi","c-f","morgan"). Default is Haldane
# 
# Output: 
#	same cross object with chosen LG fliped and all genetic distances recalculated
##################################################################################################


FlipLG <- function(cross.data, lg.flip, error.rate = 0.01, cluster = 1, m.function = "haldane"){
	numchr <- nchr(cross.data)
	 if(missing(lg.flip)) { 
		data.map <- NULL
		## estimate marker distances for all linkage groups
		data.map <- est.map(cross.data , error.prob = error.rate, sex.sp = FALSE,  offset = 0, n.cluster = cluster)
		cross.data.new <- replace.map(cross.data, data.map)
		} else {
		
		## loop within lg.flip. Apply swith.order to linkage groups designated for flip
		for (i in 1:numchr){
			if(i %in% lg.flip){
				n <- nmar(cross.data[chr=i])
				order <- rev(seq(1,n,1))
				cross.data <- switch.order(cross.data,i,order, tol=1e-5, sex.sp = FALSE, maxit=2 ) 
				## marker distances will be calculated later. here only 2 iterations 
				## are allowed to speed up the process
			}
		}
		data.map <- NULL
		## estimate marker distances for all linkage groups
		data.map <- est.map(cross.data , error.prob = error.rate, sex.sp = FALSE,  offset = 0, n.cluster = cluster, map.function = m.function)
		cross.data.new <- replace.map(cross.data, data.map)
		## Returns the modified cross object
	}
	return(cross.data.new)
}



######################################################################################
# ErrorRate:	Estimate the log likelihood (lod)of a cross object at a range of error rates
#			The error rate of a given cross object is the value that maximize the lod  
#			according to http://www.rqtl.org/tutorials/geneticmaps.pdf 
#			pag35, accesed April 2014
######################################################################################
# Input:
#	cross.data: the cross.file to analyze
# Options:
#	cluster: Number of cluster to run paralel calculations. Default is 1. requieres package snow()
#	err.end: Last error rate value used for this estimation. Default value 0.02. The more, the slower
#	err.increment: Increments the error rate value used for max. lod estimation. Default value is 0.0025. 
#                The smaller, the slower  
######################################################################################
# Output: 
#	Table with error rates and corresponding lod values. Error rates ranking from 0 to err.end,
#		 with increments of err.increment
#	A pdf file with a plot error rate vs lod. The error rate is the value that maximize lod
######################################################################################




ErrorRate <- function(cross.data, err.end = 0.03, err.increment = 0.0025, cluster = 1){
## Determine the error rates to be tested, starting from 0 up to err.end, with increments of err.increment
loglik <- err <- seq(0,err.end,err.increment)
## loop over error rates to estimate log likelihood
for(i in seq(along=err)) {
        cat(i, "of", length(err), "\n")
        tempmap <- est.map(cross.data, error.prob=err[i],tol=1e-5, maxit = 3000, n.cluster = cluster)
        loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
## lod scores associated to each error rate
lod <- (loglik - max(loglik))/log(10)
## creates output, with error values in first column and lod in second column
error.rate.result <- data.frame(cbind(err,lod))
colnames(error.rate.result) <- c("Error_rate","lod")
## Error rate correspond to the element of err that maximize the lod
## Which can be visually determine in the printed file below, or
pdf(file = "genotyping error rate.pdf",width = 5, height = 5)
plot(err, lod, cex = 1.4, pch = 21, bg = "black", main = "genotyping error rate", xlab="Genotyping error rate", xlim=c(0,err.end),ylab=expression(paste(log[10], " likelihood")))
dev.off()
## from the output table with command
## error.rate[which.max(error.rate[,2]),1]
return(error.rate.result)
}


   

###################################################################################
#  SlideDroponemarker: Droponemarker in a sliding window. It apply the function 
#    droponemarker in subsets of window.size + 4 markers along each LG. The droponemarker 
#    table of each window (trims two markers at each side) are bind to create an unique output 
#    for each linkage group.
#    Values for markers located at the extreme of each linkage group are nor reported
#    (first 2 markers and between 2 and 4 markers at the end) 
# 
####################################################################################
# Input:
#	data.drop: the cross object to analyze
# Options:
#	error.rate: Genotyping error rate of the cross object. default = 0.01
#	window size: default = 5, efective window size
#Output: 
#  dropeonemarker table output, markers at extremes of each LG are not reported
#####################################################################################


SlideDroponemarker <- function(data.drop, error.rate = 0.01, window.size = 5){ 

drop_one_all <- NULL

## loop over each linkage group of the cross object
for (k in 1:nchr(data.drop)){
	test.chr = chrnames(data.drop)[k]
	## estimate recombination fractions
	data.sm  <- est.rf(data.drop[chr = test.chr])
	## determine number of marker in LG. Used to determine the limits of the window loop
	n.markers.sm <- nmar(data.sm)
	## try to replace for nmar, check this if broken
	#n.markers.sm <- summary(data.sm)$n.mar
	
	
	## sliding window. Feeds dropeonemarker with a window.size + 4 snp window.
	##                 Results for 2 snps at extremes are discarted)
	
	## applying the droponemarker function to a sliding window
	marker.id <- colnames(data.sm$rf)
	
	
	## looping over each window. The number of windows is total markers minus 4 
	## (room for last 9 snp window), divided in 5
	
	j = 1
	# contain each window report
	dropping.one <- NULL
	# contains the compiled table for each LG
	drop.one.result <- NULL
	for (i in 1:(as.integer((n.markers.sm-4)/window.size))){
		## creates the window
		data.sm.drop <- pull.markers(data.sm,marker.id[j:(j+window.size + 3)])
		## apply droponemarker to the window
		dropping.one <- droponemarker(data.sm.drop, chr = test.chr, map.function = "haldane", error.prob = error.rate, maxit = 4000, tol = 1e-6, verbose = 1, sex.sp = FALSE)
		## trimming (two markers each side) and binding the results to generate one output per LG
		drop.one.result <- rbind(drop.one.result,dropping.one[3:(window.size + 2),])
		## the window moves window.size positions
		j = window.size*i + 1
	}
	#to process the last block, if it is smaller that window.size, but bigger than 2 snps
	if(n.markers.sm - j > 4) {
				## creates the window
		data.sm.drop <- pull.markers(data.sm,marker.id[j:n.markers.sm])
		## apply droponemarker to the window
		dropping.one <- droponemarker(data.sm.drop, chr = test.chr, map.function = "haldane", error.prob = error.rate, maxit = 4000, tol = 1e-6, verbose = 1, sex.sp = FALSE)
		## trimming (two markers each side) and binding the results to generate one output per LG
		drop.one.result <- rbind(drop.one.result,dropping.one[3:(n.markers.sm - j - 1),])
		## the window moves window.size positions
	}
	## compiling the results from all LG
	drop_one_all <- rbind(drop_one_all,drop.one.result)
	map.data <- pull.map(data.drop,as.table = TRUE)
	drop_one_all[,2] <- map.data[match(row.names(drop_one_all),row.names(map.data)),2]

}
return(drop_one_all)
}




#################################################################################
#  FilterByLODSlidingDroponemarker:  Uses droponemarker or SlideDroponemarker
#			funtions output to filter out markers from a given cross object.
#			Markers with lod and Ldiff values above threshold defined by user are droped
#			lod threshold can be: 
#				uniform: indicate one value for limLOD and no chr.list should be defined
#				specific for each chr: a vector with chromosome names (chr.list) and  
#						a vector with LOD threshold (limLOD) should be proportionate
#################################################################################
# Input:
#	data.filter: the cross object to analyze
#	drop.one.result: output from droponemarker or SlideDroponemarker functions
# Options:
#	limLOD: either a single value or a vector with LOD threshold for each chromosome in chr.list
#	limLdiff: Ldiff threshold, used to drop SNPs with values over this threshold
#	chr.list: a vector with chromosome names
###########
# Output: 
#	the cross object with markers filtered according to given options
###########


FilterByLODSlidingDroponemarker <- function(data.filter,drop.one.result, limLOD = -10, limLdiff = 3, chr.list){
	data.test <- data.filter
	#if no vector is indicated in chr.list, limLOD is used to filter all chromosomes
	if(missing(chr.list)) { 
		#subset the list of markers to drop based on limLOD and limLdiff
		todrop <- drop.one.result[drop.one.result$LOD > limLOD & drop.one.result$Ldiff > limLdiff,]
		#drop markers from the original cross object
		data.test <- drop.markers(data.test,find.marker(data.filter,chr = todrop[,1],pos = todrop[,2]))
	}
	# if selected chromosomes are indicated, run a loop over this list to select markers to drop
	else {
			todrop <- NULL
		#run a loop over chr.list to capture markers to filter 
		for(i in 1:length(chr.list)){
			#subset markers to drop based on pair chr,limLOD and limLdiff
			todrop.test <- drop.one.result[drop.one.result$chr == chr.list[i] & drop.one.result$LOD > limLOD[i] & drop.one.result$Ldiff > limLdiff,]
			todrop <- rbind(todrop,todrop.test)
		}
		# message if chr.list and limLOD have different length
		if(length(chr.list) != length(limLOD)) cat("\n","\n","chr.list and limLOD must have same length","\n","\n")
		#drop markers from the original cross object
		data.test <- drop.markers(data.test,find.marker(data.filter,chr = todrop[,1],pos = todrop[,2]))
	}
	return(data.test)
}



######################################################################################
#   rfSummary:  Creates a summary from est.rf output. Each row is a marker, first 
#				and second columns are chr/linkage group and position, third and 
#				fourth columns are distance to previous and distance to next marker,
#				and fifth and sixth columns are average recombination fraction and 
#				average lod respectively.  
######################################################################################
# Inputs:
#	cross.rf: the cross object to analyze
# Output: 
#  summary table as described above
#######################################################################################



 
rfSummary <- function(data.rf){

	### generate tables with distance gap for each marker, average rf and average lod
	## loop over each linkage group
	
	summary.table <- NULL
	
	for(n in 1:nchr(data.rf)){
		data.rf.chr <- data.rf[chr = chrnames(data.rf)[n]]
		rf.matrix <- lod.matrix <- data.rf.chr$rf
		## create columns with markers distances
		map.table <- pull.map(data.rf.chr,as.table = TRUE)
		diffst <- rep(0,dim(map.table)[1])
		diffend <- rep(0,dim(map.table)[1])
		for(i in 1:(dim(map.table)[1]-1)) diffst[i] = round(map.table[i+1,2] - map.table[i,2],2)
		for(i in 2:(dim(map.table)[1])) diffend[i] = round(map.table[i,2] - map.table[i-1,2],2)
		
		## create column with marker index
		marker.index <- seq(1,dim(map.table)[1],1)
		
		## determine rf average
		rf.matrix[upper.tri(rf.matrix , diag = TRUE)] <- 0
		rf.mean <- t((colSums(rf.matrix)+t(rowSums(rf.matrix)))/(dim(rf.matrix)[1]-1))

		## determine lod average
		lod.matrix[lower.tri(lod.matrix , diag = TRUE)] <- 0
		lod.mean <- t((t(colSums(lod.matrix))+rowSums(lod.matrix))/(dim(rf.matrix)[1]-1))
			
		map.table <- cbind(marker.index, map.table)
		map.table <- cbind(map.table,diffend)
		map.table <- cbind(map.table,diffst)
		map.table <- cbind(map.table,rf.mean)
		map.table <- cbind(map.table,lod.mean)
		
		summary.table <- rbind(summary.table, map.table)
	}
	
	return(summary.table)
	
}




