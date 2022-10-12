# Objective: This script will use the first bit of phenotype data to do a power analysis on phenotype data
# Author Darrian Talamantes (darrianrtalamantes6@gamil.com)

library(tidyverse)
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

# THis creates a shifted data set 
bmData2 <- bmData
bmData2$CopyNum1 <- bmData$CopyNum+1
bmData2 <- subset(bmData2,select = -c(CP_Ratio))
bmData2_hist <- pivot_longer(bmData2, CopyNum1:CopyNum)

## 


################################# graph ########################################

# This graph shows the difference between the CP value and the copy number ratio
ggplot(bmData_hist, aes(x=value, color=name, fill=name)) +
  geom_histogram(alpha=0.5, position="identity", binwidth=.01) +
  annotate("text", x=.94, y=6, size=3, col="red", label= paste0("CopyNum mean=", round(mean(bmData$CopyNum),2))) +
  annotate("text", x=.9, y=8, size=3, col="blue", label= paste0("CP mean=", round(mean(bmData$CP_Ratio),2))) +
  theme_bw()

ggplot(bmData2_hist, aes(x=value, color=name, fill=name)) +
  geom_histogram(alpha=0.5, position="identity", binwidth=.01) +
  theme_bw()



# makes a verticle line on position of graph
# + geom_vline(xintercept = 3, linetype="dotted", color = "blue", size=1.5)

