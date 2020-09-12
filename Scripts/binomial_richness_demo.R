#binomial richness demo

#suppose you have 3 species: A, B, and C
#each has a different probability of presence

probA <- 0.2
probB <- 0.8
probC <- 0.5


#if we treat Richness as the number of species with prob > 0.5, the answer is always the same. R = 1. 

#if we treat Richness as probabilistic, we get a distribution of Richnesses that are possible. If we do 100 draws, we get:

Richness <- NA
for (i in 1:100) {
  Richness[i] <-
    sum (
    rbinom(1,1,probA),
    rbinom(1,1,probB),
    rbinom(1,1,probC)
  )
}

hist(Richness)
mean(Richness)
median(Richness)


#so, instead of doing the total number of species at each parameter set w prob > 50%, we should be using all of those 100 draws at 
#each set of environmental parameters to calculate 100 different richnesses, and taking a mean.