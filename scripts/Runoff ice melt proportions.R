#####################################################################

# Compute proportions of runoff, sea ice melt and deep Pacific water
# at each station and for surface and bottom waters sampled during 
# AMBON cruises in the eastern Chukchi Sea in 2015 & 2017
# and map results

# Author: Franz Mueter
# last updated: May 13, 2021

#####################################################################
#
# Values of delta O-18 and salinity for end members (deep Pacific water,
# freshwater runoff, sea ice melt) from Lee Cooper. See spreadsheet:
# 'O-18 salinity mixing equation with deep Pacific end members (modified from Lee Cooper).xlsx'

library(tidyverse)

O18 <- read_csv("data_Mar2021/O18/O18.csv")

## Function for mixing equations to compute proportions of
## Pacific water, freshwater and sea ice melt based on delta O-18
## and salinity measurements and assumed delta O-18 and salinity 
## of the respective end members (adjust/update as appropriate)
## (See Mueter et al. (2021). Oceanography paper for details)
mix <- function(delO18, salinity, O18.Pac = -0.17, O18.FW = -21.5,
                O18.ice = -1, sal.Pac = 34.7, sal.FW = 0, sal.ice = 4)
{
  expr1 <- (salinity - sal.Pac) * (O18.FW - O18.Pac)
  expr2 <- (sal.FW - sal.Pac) * (delO18 - O18.Pac)
  expr3 <- (sal.ice - sal.Pac) * (O18.FW - O18.Pac)
  expr4 <- (sal.FW - sal.Pac) * (O18.ice - O18.Pac)
  expr5 <- (salinity - sal.Pac) / (sal.FW - sal.Pac)
  expr6 <- (sal.ice - sal.Pac) / (sal.FW - sal.Pac)
  prop.melt <-  pmax(0, (expr1 - expr2) / (expr3 - expr4))
  prop.runoff <- expr5 - expr6 * prop.melt
  data.frame(sea.ice.melt = prop.melt, runoff = prop.runoff, 
             Pacific = 1 - prop.melt - prop.runoff)
}

O18 <- na.omit(O18)
O18 <- cbind(O18, with(O18, mix(O18, salinity)))

# Change station names to master station names:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
any(duplicated(master$station))  # Should be FALSE

O18 <- left_join(O18, master, by = "station")
any(is.na(O18$station.master))
# Make sure CTD locations are close to 'master' locations:
range(O18$latitude - O18$latitude.master)
range(O18$longitude - O18$longitude.master)

f.bot <- O18 %>% filter(layer == "Bottom") %>% 
  rename(melt.bot = sea.ice.melt, runoff.bot = runoff, Pac.bot = Pacific)
f.sfc <- O18 %>% filter(layer == "Surface") %>% 
  rename(melt.sfc = sea.ice.melt, runoff.sfc = runoff, Pac.sfc = Pacific)

# Aggregate duplicate measurements (stations sampled twice):
table(f.bot$cruise)
dupl <- f.sfc %>% filter(cruise=="AMBON2015") %>% select(station.master) %>% duplicated()
# Duplicated stations:
st.dupl <- f.sfc %>% filter(cruise=="AMBON2015") %>% filter(dupl) %>% pull(station.master)
# These are stations where transects intersect, sampled twice.
# Compare duplicates:
f.sfc %>% filter(cruise=="AMBON2015") %>% filter(station.master %in% st.dupl) %>% arrange(station.master)
# check AMBON2017 and surface values
# O18 values are generally very similar, hence aggregate! 

f.bot <- f.bot %>% select(-station, -layer) %>% 
  group_by(cruise, station.master) %>% 
  summarize(across(everything(), mean)) 

f.sfc <- f.sfc %>% select(-station, -layer) %>% 
  group_by(cruise, station.master) %>% 
  summarize(across(everything(), mean))

write.csv(f.bot, file="data_Mar2021/O18/bottom_water_fractions.csv", row.names=F)
write.csv(f.sfc, file="data_Mar2021/O18/surface_water_fractions.csv", row.names=F)


################################################################################
# Map values to check if they are sensible:
load("./Base_Map/Chukchi_basemap.RData")  # loads 'CSmap', a ggplot object

# Aggregate duplicated measurements by cruise, station, layer:
O18 <- O18 %>% select(-station) %>% 
  group_by(cruise, station.master, layer) %>% 
  summarize(across(everything(), mean)) 

# Map bottom fractions of freshwater runoff and melt water:

O18 <- mutate(O18, layer = ordered(layer, levels = c("Surface", "Bottom")))

png("plots/fraction_runoff.png")
CSmap + geom_point(aes(longitude, latitude, color=runoff), size=4, 
                   data=O18) + facet_grid(vars(layer),vars(cruise)) +  
  scale_color_viridis_c()
dev.off()

png("plots/fraction_meltwater.png")
CSmap + geom_point(aes(longitude, latitude, color=sea.ice.melt), size=4, 
                   data=O18) + facet_grid(vars(layer),vars(cruise)) +  
  scale_color_viridis_c()
dev.off()

png("plots/fraction_Pacific_water.png")
CSmap + geom_point(aes(longitude, latitude, color=Pacific), size=4, 
                   data=O18) + facet_grid(vars(layer),vars(cruise)) +  
  scale_color_viridis_c()
dev.off()


