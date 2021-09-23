######################################################################

# Import / prepare 2015 & 2017 haul and cpue data or zooplankton 150micron net: 
# Last update: September 6, 2021

# This selects only small zooplankton species that are not caught in the 
# large-mesh net (< 2mm as adults)

# CPUE (biomass) and haul info written to 'zoop150.RData'

######################################################################

require(tidyverse)

## Compute CPUE by number and year

# Import 2015 species data:
zoop15 <- read_csv("./data_Mar2021/Zooplankton/AMBON2015150.csv", na=c("#N/A", "n/a","")) %>% 
  select(cruise=Cruise,
         station = Station,
         species = Accepted_Organism_Identification,
         cpue = "Biomass_[mg dw/m3]",
         depth = "Bottom_Depth_[m]",
         net.depth = "Depth_[m]",
         latitude = "Latitude_[decimal _degrees_north]",
         longitude = "Longitude_[decimal_degrees_east]") 

# Import 2017 specimen data:
zoop17 <- read_csv("./data_Mar2021/Zooplankton/AMBON2017150.csv", na=c("n/a", "#N/A", "")) %>% 
  select(cruise=Cruise,
         station = Station,
         species = Accepted_Organism_Identification,
         cpue = "Biomass_[mg dw/m3]",
         depth = "Bottom_Depth_[m]",
         net.depth = "Depth_[m]",
         latitude = "Latitude_[decimal _degrees_north]",
         longitude = "Longitude_[decimal_degrees_east]")

# For initial taxon determinations:
# taxa15 <- zoop15 %>% select(species) %>% group_by(species) %>% summarize(N15=n())
# taxa17 <- zoop17 %>% select(species) %>% group_by(species) %>% summarize(N17=n())
# write.table(full_join(taxa15, taxa17), file="./data_Mar2021/Zooplankton/prelim_taxa.csv", 
#            row.names = F, sep=",")

# Import taxa names and add to data frames:
taxa <- read_csv("./data_Mar2021/Zooplankton/zoop150_taxa.csv")
# Set all 'large' zooplankton taxa to 'other.unid'
taxa <- taxa %>% mutate(taxon = replace(taxon, large == 1, "other.unid")) 
zoop15 <- left_join(zoop15, select(taxa, species, taxon))
zoop17 <- left_join(zoop17, select(taxa, species, taxon))

# Combine both years and drop 'other' and unidentified taxa:
zoop <- rbind(zoop15, zoop17) %>% filter(taxon != "other.unid")

### Get cpue by species for each cruise & station
zoop150.cpue <- zoop %>% group_by(cruise, station, taxon) %>% 
  summarize(N=sum(cpue)) %>% 
  pivot_wider(id_cols = c("cruise", "station", "taxon"), 
              names_from = taxon, values_from = N, values_fill=0)

# Extract and save haul IDs
zoop150.hauls <- select(zoop150.cpue, 1:2)

# Save cpue as matrix for analysis:
zoop150.cpue <- as.matrix(zoop150.cpue[,-c(1:2)])

# Add master station names and latitude / longitude:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
zoop150.hauls <- left_join(zoop150.hauls, master, by = "station")
if(any(is.na(zoop150.hauls$station.master))) stop("Unrecognized station name")  

rm(taxa, zoop15, zoop17)
zoop150.hauls <- ungroup(zoop150.hauls)
save(zoop150.hauls, zoop150.cpue, file="./data_Mar2021/Zooplankton/zoop150.RData")

