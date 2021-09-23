######################################################################

# Import / prepare 2015 & 2017 haul and cpue data for large zooplankton
# This includes all taxa retained in 505 micron net, regardless of adult size

# Author: Franz Mueter
# Last update: MARCH 20, 2021

# superseded by 'AMBON zoop505 cpue by taxon.R' for final manuscript version

######################################################################

require(tidyverse)

## Compute CPUE by number and year

# Import 2015 species data:
zoop15 <- read_csv("./data_Mar2021/Zooplankton/AMBON2015505.csv", na=c("#N/A", "n/a","")) %>% 
  select(cruise=Cruise,
         station = Station,
         species = Accepted_Organism_Identification,
         cpue = "Biomass_[mg dw/m3]",
         depth = "Bottom_Depth_[m]",
         net.depth = "Depth_[m]",
         latitude = "Latitude_[decimal _degrees_north]",
         longitude = "Longitude_[decimal_degrees_east]") 

# Import 2017 specimen data:
zoop17 <- read_csv("./data_Mar2021/Zooplankton/AMBON2017505.csv", na=c("n/a", "#N/A", "")) %>% 
  select(cruise=Cruise,
         station = Station,
         species = Accepted_Organism_Identification,
         cpue = "Biomass_[mg dw/m3]",
         depth = "Bottom_Depth_[m]",
         net.depth = "Depth_[m]",
         latitude = "Latitude_[decimal _degrees_north]",
         longitude = "Longitude_[decimal_degrees_east]")

# For initial taxon determinations:
#taxa15 <- zoop15 %>% select(species) %>% group_by(species) %>% summarize(N15=n())
#taxa17 <- zoop17 %>% select(species) %>% group_by(species) %>% summarize(N17=n())
#write.table(full_join(taxa15, taxa17), file="./data_Mar2021/Zooplankton/prelim_taxa.csv", 
#          row.names = F, sep=",")

# Import taxa names and add to data frames:
taxa <- read_csv("./data_Mar2021/Zooplankton/zoop505_taxa.csv")
zoop15 <- left_join(zoop15, taxa)
zoop17 <- left_join(zoop17, taxa)

# Combine both years and drop 'other' and unidentified taxa:
zoop <- rbind(zoop15, zoop17) %>% filter(taxon != "other.unid")

### Get cpue by species for each cruise & station
zoop505.cpue <- zoop %>% group_by(cruise, station, taxon) %>% 
  summarize(N=sum(cpue)) %>% 
  pivot_wider(id_cols = c("cruise", "station", "taxon"), 
              names_from = taxon, values_from = N, values_fill=0)

# Extract and save haul IDs
zoop505.hauls <- select(zoop505.cpue, 1:2)

# Add master station names and latitude / longitude:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
zoop505.hauls <- left_join(zoop505.hauls, master, by = "station")
if(any(is.na(zoop505.hauls$station.master))) stop("Unrecognized station name")  

# Save cpue as matrix for analysis:
rm(taxa,zoop15,zoop17)
zoop505.cpue <- as.matrix(zoop505.cpue[,-c(1:2)])
zoop505.hauls <- ungroup(zoop505.hauls)
save(zoop505.hauls, zoop505.cpue, file="./data_Mar2021/Zooplankton/zoop505.RData")

