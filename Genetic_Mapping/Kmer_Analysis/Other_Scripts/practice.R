library("pwr")
install.packages("pwr")
remove.packages("pwr")

# power is the probability that a test will correctly give you a small p-value
cohen.ES(test = "r", size = "medium")

##### Power calculations for correlation test
# r is correlation coefficient 
# n is observations
r_power <- pwr.r.test(r = .5,
                      sig.level = .05,
                      power = .8,
                      alternative = "two.sided")
r_power
plot(r_power)

# d is the cohens D
# Type is the type of t test (usually two sample or paired)
t_power <- pwr.t.test(d=.3,
                      sig.level = .05,
                      power= .8,
                      type = "two.sample",
                      alternative = "two.sided")
t_power
plot(t_power)

dList <- list()
nList <- list()
for (i in seq(from=.1, to=.5, by=.1)){
  t_power <- pwr.t.test(d=i,
                        sig.level = .05,
                        power= .8,
                        type = "two.sample",
                        alternative = "two.sided")
  
  n <- t_power$n
  d <- t_power$d
  dList <- append(dList,d)
  nList <- append(nList,n)
  
}

# Turns lists into a dataframe
power <- do.call(rbind, Map(data.frame, d=dList, n=nList))



