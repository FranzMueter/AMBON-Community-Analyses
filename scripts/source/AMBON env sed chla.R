#####################################################################

# Read in complete station list and associated environmental data
# for analysis of AMBON 2015/2017 data

# Author: Franz Mueter 
# Last updated: 5/15/2021

# Results saved as 'all_env.RData'

#####################################################################

require(tidyverse)


# Station list from CTD bottle data (master list)
St15 <- read_csv("./data_Mar2021/CTD/2015/AMBON2015_station_masterlist.csv") %>% 
  mutate(Cruise = "AMBON2015") 
St17 <- read_csv("./data_Mar2021/CTD/2017/AMBON2017_StationList.csv") %>% 
  mutate(Cruise = "AMBON2017") 

Stations <- bind_rows(St15, St17) %>% 
  rename(cruise = Cruise,
         station = 'Station Name',
         alt.station = 'Alternate Station Name',
         latitude=Latitude,
         longitude = Longitude,
         depth = 'Depth (m)') %>% select(-'Station #') %>% 
  relocate(cruise)
rm(St15, St17)
# Eliminate row with '-999' in 2015 (repeat station)
Stations <- filter(Stations, depth != -999)

# Import surface and bottom water properties from fish haul data:
# (may use some different station names)
ts <- read_csv("./data_Mar2021/StationData/AMBON_hauls.csv") %>% 
  rename(cruise = Cruise,
         date = eventDate,
         station = locationID,
         latitude = decimalLatitude,
         longitude = decimalLongitude,
         depth = maximumDepthInMeters) %>% 
  mutate(date = lubridate::date(date),
         cruise = str_replace(cruise, "AMBON", "AMBON20")) %>% 
  select(-samplingProtocol, -Quantitative, -Distance, -Width, -Haul) %>% 
  group_by(cruise, station) %>% summarize(across(everything(), mean)) %>% 
  ungroup()
  

# Import sediment & environmental data from Jackie Grebmeier
# (used some different station names but changed in csv files)
env15 <- read_csv("./data_Mar2021/Macroinfauna/AMBON15_env.csv")
env17 <- read_csv("./data_Mar2021/Macroinfauna/AMBON17_env.csv")

env <- bind_rows(mutate(env15, cruise="AMBON2015"), mutate(env17, cruise="AMBON2017")) %>% 
  select(-DataDate, -GMT, - abundance, -biomass, -Carbon, -richness, -SW.div, -SW.even) %>% relocate(cruise)
rm(env15,env17)

# Change station names to master station names:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
any(duplicated(master$station))  # Should be FALSE

ts <- left_join(ts, master, by = "station")
env <- left_join(env, master, by="station")
Stations <- left_join(Stations, master, by="station")
# Should all be FALSE:
any(is.na(ts$station.master))
any(is.na(env$station.master))
any(is.na(Stations$station.master))

# Get stratification index:
strat <- read_csv("data_Mar2021/CTD/stratification.csv")
strat <- left_join(strat, master, by = "station")

# Get surface and bottom fractions of sea ice melt, runoff and deep Pacific water
f.bot <- read_csv("data_Mar2021/O18/bottom_water_fractions.csv")
f.sfc <- read_csv("data_Mar2021/O18/surface_water_fractions.csv")

# Combine all environmental data
all.env <- select(Stations, cruise, station.master, latitude.master, longitude.master, depth) %>% 
  left_join(select(ts, cruise, station.master, SST, BT, SSS, Bsal), by=c("cruise", "station.master")) %>% 
  left_join(select(env, cruise, station.master, NH4:sed.chla, n.grab:delN15), by=c("cruise", "station.master")) %>% 
  left_join(select(strat, cruise, station.master, strat), by=c("cruise", "station.master")) %>% 
  left_join(select(f.bot, cruise, station.master, melt.bot, runoff.bot), by=c("cruise", "station.master")) %>% 
  left_join(select(f.sfc, cruise, station.master, melt.sfc, runoff.sfc), by=c("cruise", "station.master"))
  

# Some stations sampled more than once (where transects cross). For these stations, I averaged all environmental data:
f <- function(x) mean(x, na.rm=T)
all.env <- all.env %>% group_by(cruise, station.master) %>% summarize(across(everything(), f))

rm(env, master, Stations, ts, strat, f, f.bot, f.sfc)
all.env <- ungroup(all.env)

save(all.env, file="./data_Mar2021/StationData/all_env.RData")


# Alternative simple stratification index based on temperature and salinity
# PC1 of temperature and salinity differences as quick index of stratification:
#all.env <- all.env %>% mutate(strat.t = SST - BT, strat.s = SSS - BT) 

#strat <- select(all.env, cruise, station.master, strat, strat.t, strat.s)
#strat <- na.omit(strat)

# PCA:
#strat.pca <- princomp(strat[,-c(1:3)], cor=T)
#strat <- mutate(strat, strat.ts=strat.pca$scores[,"Comp.1"])

