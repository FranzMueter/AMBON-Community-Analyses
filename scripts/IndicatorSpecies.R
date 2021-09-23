#####################################################################

# Compute magnitude of differences in abundance between 2015 and
# 2017 for all taxa in each of eight assemblages, and rank taxa based 
# on estimated differences relative to uncertainty in estimates

# Author: Franz Mueter 
# Last updated: 9/8/2021

#####################################################################


library(tidyverse)
library(mgcv)

source("./scripts/AMBON_Import_all_data.R")


# To use a common framework across assemblages and taxa, I 
# compared log-transformed CPUE between cruises (years) using 
# a simple linear model (equivalent to t-test), despite obvious
# violations due to non-normality & large number of zeros
# I used p-values only as a cutoff and sorted magnitude of
# differences based on t-values

###### For all assemblages with cpue data: (for microbes see below)
## Select assemblage (and change name for output below):
hdr <- epi.hauls
cpue <- epi.cpue

## Select stations sampled in both years to minimize biases
j <- hdr$cruise == "AMBON2015"
stns <- intersect(hdr$station.master[j], hdr$station.master[!j]) 
sub <- hdr$station.master %in% stns
hdr <- hdr[sub,]
cpue <- cpue[sub,]
cpue <- cpue[,apply(cpue,2,sum)>0]
dim(cpue)

out <- matrix(NA, nrow=ncol(cpue), ncol=8)
dimnames(out) <- list(dimnames(cpue)[[2]], c("N", "FO.15", "FO.17", "change", "diff", "se.diff", "t", "p"))
# Outputs are:
# N:              Number of stations sampled in both years (repeated for each taxon)
# FO15, FO17:     Frequency of occurrence in 2015 / 2017
# change:         Change in CPUE (2017 - 2015) relative to 2015
# diff, se.diff:  Difference in log-transformed CPUE between 2015 and 2017 and its standard error
# t, p:           Corresponding t-statistic and p-value (based on normality assumption)


for(i in dimnames(cpue)[[2]]) {
  dat <- cbind(hdr[,c("cruise", "latitude.master", "longitude.master")], cpue = cpue[,i])
#  .form <- as.formula("log(cpue + 0.001) ~ cruise + s(latitude.master, longitude.master, k=15)")
#  fit <- gam(.form, data=dat)
#  fit <- try(gls(log(cpue+0.001) ~ cruise, correl=corAR1(form = ~ latitude.master + longitude.master | cruise), data=dat))
  fit <- lm(log(cpue+0.001) ~ cruise, data=dat)
  mu <- with(dat, tapply(cpue, cruise, mean))
  out[i,1] <- length(stns)
  out[i,2:3] <- with(dat, tapply(cpue, cruise, function(x) sum(x>0)))
  out[i,4] <- round((mu[2] - mu[1]) / mu[1],2)  
  out[i,5:8] <- summary(fit)$coef[2,]
}
out <- out[rev(order(abs(out[,"t"]))),]
# infauna.diff <- out[out[,"p"] < 0.05,]  # Could use different cutoff by group
# fish.diff <- out[out[,"p"] < 0.05,]
# epifauna.diff <- out[out[,"p"] < 0.05,]
# zoop150 <- out[out[,"p"] < 0.05,]
# zoop505 <- out[out[,"p"] < 0.05,]
# birds.diff <- out[out[,"p"] < 0.05,]

inf <- cbind(as.data.frame(infauna.diff),assmbl = "Infauna")
epi <- cbind(as.data.frame(epifauna.diff),assmbl = "Epiauna")
fish <- cbind(as.data.frame(fish.diff),assmbl = "Fish")
zoop150 <- cbind(as.data.frame(zoop150),assmbl = "Zoop150")
zoop505 <- cbind(as.data.frame(zoop505),assmbl = "Zoop505")
bird <- cbind(as.data.frame(birds.diff),assmbl = "Birds")

all <- rbind(zoop150, zoop505, inf, epi, fish, bird)
write.table(all, "results/Changes_2015-2017.csv", sep=",")

# Revised zooplankton changes:
write.table(rbind(zoop150, zoop505), "results/Zoop_changes_2015-2017_revised.csv", sep=",")


###### For microbes (presence / absence data at highly aggregated level, using binomial model)
## Select assemblage (and change name for output below):
hdr <- sfc.microbe.hdr
pa <- sfc.bact.agg

## Select stations sampled in both years to minimize biases
j <- hdr$cruise == "AMBON2015"
stns <- intersect(hdr$station.master[j], hdr$station.master[!j]) 
sub <- hdr$station.master %in% stns
hdr <- hdr[sub,]
pa <- pa[sub,]
pa <- pa[,apply(pa,2,sum)>0]
dim(pa)

out <- matrix(NA, nrow=ncol(pa), ncol=8)
dimnames(out) <- list(dimnames(pa)[[2]], c("N", "FO.15", "FO.17", "change", "diff", "se.diff", "t", "p"))

for(i in dimnames(pa)[[2]]) {
  dat <- cbind(hdr[,c("cruise", "latitude.master", "longitude.master")], pa = pa[,i])
  fit <- glm(pa ~ cruise, data=dat, family=binomial)
  mu <- with(dat, tapply(pa, cruise, mean))
  out[i,1] <- length(stns)
  out[i,2:3] <- with(dat, tapply(pa, cruise, function(x) sum(x>0)))
  out[i,4] <- round((mu[2] - mu[1]) / mu[1],2)  
  out[i,5:8] <- summary(fit)$coef[2,]
}
out <- out[rev(order(abs(out[,"t"]))),]
agg.diff <- out[out[,"p"] < 0.1,]

bact <- cbind(as.data.frame(agg.diff),assmbl = "Bacteria")
#prot <- cbind(as.data.frame(agg.diff),assmbl = "Protists")
  
mic <- rbind(bact,prot)
write.table(mic, "results/Changes_2015-2017_microbes.csv", sep=",")

