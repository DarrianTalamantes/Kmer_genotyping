# Objective: This script will use the first bit of phenotype data to do a power analysis on phenotype data
# Author Darrian Talamantes (darrianrtalamantes6@gamil.com)

library(tidyverse)
library(data.table)
# Setting Variables
biomassData="/home/drt83172/Documents/QPCR_Data_Wrangler/Program/output/Data_for_project/Biomass_Data/Biomass_Data.csv"




bmData <- read.table(biomassData, sep = ",", header = TRUE, row.names = 1)

# Checking differences between CP and Copy number.
# Technically copy number should be better to use because the quantities of amplicon can be different for the same CP number
summary(bmData$CP_Ratio)
sd(bmData$CP_Ratio)
summary(bmData$CopyNum)
sd(bmData$CopyNum)

bmData_hist <- pivot_longer(bmData, CP_Ratio:CopyNum)

# This function takes an effect size and spits out an n for said effect size
poweranalysis <- function(effectsize){
t_power <- pwr.t.test(d=effectsize,
                      sig.level = .05,
                      power= .8,
                      type = "two.sample",
                      alternative = "two.sided")
num <- t_power$n
return(num)
}
# Using the function above to learn the correct N for effect sizes from .2 to .8
coD <- c(.2,.3,.4,.5,.6,.7,.8)
powerA <- data.table(coD)
powerA$N <- 0 
powerA[1,2]

for (i in 1:length(coD)){
  D <-(coD[i])
  num <- poweranalysis(D)
  powerA[i,2] <- num
}

# Lets find out where the top 5% is for the cp ratio data 
totSamples <- nrow(bmData) 
target <- ceiling(totSamples * .95)
CpData <- (bmData$CP_Ratio)
CpData <- sort(CpData)
nth(CpData, target)
Std<-sd(CpData)
avg<-mean(CpData)

Std




 ################################# graph ########################################

# This graph shows the difference between the CP value and the copy number ratio
ggplot(bmData_hist, aes(x=value, color=name, fill=name)) +
  geom_histogram(alpha=0.5, position="identity", binwidth=.01) +
  annotate("text", x=.94, y=6, size=3, col="red", label= paste0("CopyNum mean=", round(mean(bmData$CopyNum),2))) +
  annotate("text", x=.9, y=8, size=3, col="blue", label= paste0("CP mean=", round(mean(bmData$CP_Ratio),2))) +
  theme_bw()


# Create a plot using the cohens D's and N's needed to meet them at good levels
ggplot(powerA, aes(x=coD, y=N)) +
  geom_point(shape=16, size=3) +
  xlab("Cohen's D") +
  ggtitle("N needed for specific effect sizes with \n power of .8 and p=.05") + 
  theme(plot.title = element_text(hjust = 0.5))


# makes a verticle line on position of graph
# + geom_vline(xintercept = 3, linetype="dotted", color = "blue", size=1.5)

