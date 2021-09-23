#####################################################################

# Determine number of stations sampled by year and community, 
# diversity metrics (Hill numbers) by year and assemblage
# (standardized to sample 50 stations), and number of taxa in 
# the analysis for two AMBON cruises in 2015 and 2017

# Author: Franz Mueter 
# Last updated: 9/8/2021

#####################################################################

source("./scripts/AMBON_Import_all_data.R")

###########################################################################
# Number of species and number of stations sampled per year by assemblage 
###########################################################################
sfc.microbe.hdr %>% group_by(cruise,station.master) %>% 
  summarize(N=n()) %>% summarize(N=n())
ncol(sfc.bacteria)
ncol(sfc.protist)
sum(apply(sfc.bacteria, 2, function(x) sum(x>0)) < 4)
sum(apply(sfc.protist, 2, function(x) sum(x>0)) < 4)

zoop150.hauls %>% group_by(cruise,station.master) %>% 
  summarize(N=n()) %>% summarize(N=n())
ncol(zoop150.cpue); ncol(zoop505.cpue)
sum(apply(zoop150.cpue, 2, function(x) sum(x>0)) < 4)

zoop505.hauls %>% group_by(cruise,station.master) %>% 
  summarize(N=n()) %>% summarize(N=n())
ncol(zoop150.cpue); ncol(zoop505.cpue)
sum(apply(zoop505.cpue, 2, function(x) sum(x>0)) < 4)

infauna.hauls %>% group_by(cruise,station.master) %>% 
  summarize(N=n()) %>% summarize(N=n())
ncol(infauna.cpue)
sum(apply(infauna.cpue, 2, function(x) sum(x>0)) < 4)

epi.hauls %>% group_by(cruise,station.master) %>% 
  summarize(N=n()) %>% group_by(cruise) %>% summarize(N=n())
ncol(epi.cpue)
sum(apply(epi.cpue, 2, function(x) sum(x>0)) < 4)

fish.hauls %>% group_by(cruise,station.master) %>% 
  summarize(N=n()) %>% summarize(N=n())
ncol(fish.cpue)
sum(apply(fish.cpue, 2, function(x) sum(x>0)) < 4)

bird.hdr %>% group_by(cruise,station.master) %>% 
  summarize(N=n()) %>% summarize(N=n())
ncol(bird.cpue)
sum(apply(bird.cpue, 2, function(x) sum(x>0)) < 4)



###########################################################################
## Summarize Diversity measures (Hill numbers, q) by year and order for
## each assemblage, standardized to sample size of N=50
###########################################################################
require(iNEXT)

sfc.bact <- list('2015' = t(sfc.bacteria[sfc.microbe.hdr$cruise == "AMBON2015",]),
                 '2017' = t(sfc.bacteria[sfc.microbe.hdr$cruise == "AMBON2017",]))
qD.bact <- estimateD(sfc.bact, datatype="incidence_raw", base="size", level=50)

#mid.bact <- list('2015' = t(mid.bacteria[mid.microbe.hdr$cruise == "AMBON2015",]),
#                 '2017' = t(mid.bacteria[mid.microbe.hdr$cruise == "AMBON2017",]))

sfc.prot <- list('2015' = t(sfc.protist[sfc.microbe.hdr$cruise == "AMBON2015",]),
                 '2017' = t(sfc.protist[sfc.microbe.hdr$cruise == "AMBON2017",]))
qD.prot <- estimateD(sfc.prot, datatype="incidence_raw", base="size", level=50)

sfc.meta <- list('2015' = t(sfc.metazoa[sfc.microbe.hdr$cruise == "AMBON2015",]),
                 '2017' = t(sfc.metazoa[sfc.microbe.hdr$cruise == "AMBON2017",]))
qD.meta <- estimateD(sfc.meta, datatype="incidence_raw", base="size", level=50)

zoop150 <- list('2015' = t(zoop150.cpue[zoop150.hauls$cruise == "AMBON2015",]),
                '2017' = t(zoop150.cpue[zoop150.hauls$cruise == "AMBON2017",]))
qD.zoop150 <- estimateD(zoop150, datatype="incidence_raw", base="size", level=50)

zoop505 <- list('2015' = t(zoop505.cpue[zoop505.hauls$cruise == "AMBON2015",]),
                '2017' = t(zoop505.cpue[zoop505.hauls$cruise == "AMBON2017",]))
qD.zoop505 <- estimateD(zoop505, datatype="incidence_raw", base="size", level=50)

inf <- list('2015' = t(infauna.cpue[infauna.hauls$cruise == "AMBON2015",]),
            '2017' = t(infauna.cpue[infauna.hauls$cruise == "AMBON2017",]))
qD.inf <- estimateD(inf, datatype="incidence_raw", base="size", level=50)

epi <- list('2015' = t(epi.cpue[epi.hauls$cruise == "AMBON2015",]),
            '2017' = t(epi.cpue[epi.hauls$cruise == "AMBON2017",]))
qD.epi <- estimateD(epi, datatype="incidence_raw", base="size", level=50)

fish <- list('2015' = t(fish.cpue[fish.hauls$cruise == "AMBON2015",]),
             '2017' = t(fish.cpue[fish.hauls$cruise == "AMBON2017",]))
qD.fish <- estimateD(fish, datatype="incidence_raw", base="size", level=50)

bird <- list('2015' = t(bird.cpue[bird.hdr$cruise == "AMBON2015",]),
             '2017' = t(bird.cpue[bird.hdr$cruise == "AMBON2017",]))
qD.bird <- estimateD(bird, datatype="incidence_raw", base="size", level=50)

qD.all <- rbind(qD.bact, qD.prot, qD.zoop150, qD.zoop505, qD.inf, qD.epi,qD.fish,qD.bird)

# Hill numbers by year and order:
round(matrix(qD.all$qD, ncol=6,byrow=T))


###########################################################################
## Table of taxa with frequency of occurrence by year and overall FO
###########################################################################
## Run previous section first to get species by station matrices

# Select assemblage and choose file name:
dat <- zoop505
name <- "zoop505"
x15 <- apply(dat[['2015']], 1, function(x) sum(x>0))
x17 <- apply(dat[['2017']], 1, function(x) sum(x>0))
x.tot <- apply(cbind(dat[['2015']], dat[['2017']]), 1, function(x) sum(x>0))
# Add total number of samples:
x15 <- c(n1<-ncol(dat[['2015']]), x15)
x17 <- c(n2<-ncol(dat[['2017']]), x17)
x.tot <- c(n1+n2, x.tot)
write.table(data.frame(FO.2015=x15, FO.2017=x17, FO.tot = x.tot), 
            file = paste("./results/FO_", name, ".csv", sep=""), sep=",")

