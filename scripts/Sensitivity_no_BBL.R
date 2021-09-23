#####################################################################

# Prepare data for sensitivity analysis to assess impact of a transect 
# sampled in 2017 only (BBL transect) on estimates of diversity 

# Author: Franz Mueter 
# Last updated: 9/8/2021

#####################################################################

# This script removes all stations along the BBL transect, which was
# only samples in 2017. To assess impacts, re-run diversity analysis
# ('AMBON diversity 2015-2017.R') to estimate species accumulation 
# curves without this transect included:
# (Results are virtually identical)

# Load all data:
source("./scripts/AMBON_Import_all_data.R")

sum(j <- substr(bird.hdr$station.master,1,3) != "BBL")
bird.cpue <- bird.cpue[j,]
bird.hdr <- bird.hdr[j,]

sum(j <- substr(epi.hauls$station.master,1,3) != "BBL")
epi.cpue <- epi.cpue[j,]
epi.hauls <- epi.hauls[j,]

sum(j <- substr(infauna.hauls$station.master,1,3) != "BBL")
infauna.cpue <- infauna.cpue[j,]
infauna.hauls <- infauna.hauls[j,]

sum(j <- substr(fish.hauls$station.master,1,3) != "BBL")
fish.cpue <- fish.cpue[j,]
fish.hauls <- fish.hauls[j,]

sum(j <- substr(zoop150.hauls$station.master,1,3) != "BBL")
zoop150.cpue <- zoop150.cpue[j,]
zoop150.hauls <- zoop150.hauls[j,]

sum(j <- substr(zoop505.hauls$station.master,1,3) != "BBL")
zoop505.cpue <- zoop505.cpue[j,]
zoop505.hauls <- zoop505.hauls[j,]





