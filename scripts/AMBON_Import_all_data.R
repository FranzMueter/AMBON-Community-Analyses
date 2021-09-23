#####################################################################

# Import all data sets for analyses. Usually sourced from other scripts
# Re-run appropriate scripts in 'source' if data for an assemblage is
# updated or for analysis of subsets

# Author: Franz Mueter 
# Last updated: 9/8/2021

#####################################################################

# Load required packages
library(tidyverse)

### load previously saved objects
load("./data_Mar2021/StationData/all_env.RData") # # from 'AMBON env sed chla.R'
load("./data_Mar2021/Fish/fish.RData") # from 'AMBON fish cpue by taxon.R'
load("./data_Mar2021/Epifauna/epifauna.RData") # from 'AMBON epifauna cpue by taxon.R'
load("./data_Mar2021/Macroinfauna/infauna.RData") # from 'AMBON infauna cpue by taxon.R'
load("./data_Mar2021/Zooplankton/zoop150.RData") # from 'AMBON zoop150 cpue by taxon.R' 
load("./data_Mar2021/Zooplankton/zoop505.RData") # from 'AMBON zoop505 cpue by taxon.R' 
load("./data_Mar2021/Microbes/sfc_microbes.RData") # from 'AMBON zoop505 cpue by taxon.R' 

# Seabirds
birds <- read_csv("./data_Mar2021/Seabirds/AMBON_Seabird_Station_Densities_Apportioned.csv") %>% 
  rename(cruise = Year,
         station='Station_Name') %>% 
  mutate(cruise = paste("AMBON", cruise, sep=""))
bird.hdr <- birds[,1:4]
bird.cpue <- as.matrix(birds[,-c(1:4)])
rm(birds)

# Add master station names and latitude / longitude:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
bird.hdr <- left_join(bird.hdr, master, by = "station")
if(any(is.na(bird.hdr$station.master))) stop("Unrecognized station name")  


