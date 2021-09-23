######################################################################

# Import / prepare 2015 & 2017 haul and catch data for fish assemblage

# Author: Franz Mueter
# Last update: MARCH 19, 2021

# CPUE (by number) and haul info written to 'fish.RData'

######################################################################

require(tidyverse)

## Compute CPUE by number and year

# Import haul data from Excel file:
hauls <- read_csv("./data_Mar2021/StationData/AMBON_hauls.csv") %>% 
  rename(cruise = Cruise,
         station = locationID,
         latitude = decimalLatitude,
         longitude = decimalLongitude,
         depth = maximumDepthInMeters,
         gear = samplingProtocol,
         haul = Haul) %>% 
  filter(gear == "PSBTA", Quantitative == "Quantitative") %>% 
  mutate(area.swept = Distance * Width, 
         date = lubridate::date(eventDate),
         uniqueID = paste(cruise, station, haul, sep="-")) %>% 
  select(-Quantitative, - gear, -Distance, -Width, -eventDate)

# Import specimen data
specimen <- read_csv("./data_Mar2021/Fish/specimen_data.csv",
                     col_types = "ccciciidi") %>% 
  rename(cruise = Cruise,
         gear = samplingProtocol,
         haul = Haul,
         length = TotalLength,
         species = scientificName,
         station = locationID,
         lengthBin = LengthClass10mm,
         count = Count) %>% 
  filter(gear == "PSBTA") %>% 
  mutate(uniqueID = paste(cruise, station, haul, sep="-")) %>% 
  filter(uniqueID %in% hauls$uniqueID)

# Check if station names are identical between the two data sheets:
if(!identical(sort(unique(hauls$station)), sort(unique(specimen$station)))) stop("station names differ")

### Get counts (N) by species for each cruise & station
# This averages over hauls if more than one haul per station!! 
N <- specimen %>% group_by(cruise, station, haul, species) %>% 
  summarize(N=sum(count)) %>% 
  group_by(cruise, station, species) %>% 
  summarize(N = mean(N)) %>% 
  pivot_wider(id_cols = c("cruise", "station", "species"), 
              names_from = species, values_from = N, values_fill=0)

# Get corresponding hauls
temp <- select(N, 1:2) %>%  left_join(hauls) 

# Aggregate duplicated stations with more than one haul, 
# use average area swept from hauls (as well as averages of all
# other variables (although they may be identical)

# Duplicated hauls at:
j <- duplicated(paste(temp$cruise, temp$station))
temp$station[j]  # Duplicated station (2 hauls)

fish.hauls <- temp %>%  group_by(cruise,station) %>% select(-uniqueID) %>% summarize(across(everything(), mean))
if(!identical(fish.hauls$station, N$station)) stop("stations don't match")  # Check
rm(j, temp)

fish.cpue <- N[,-c(1,2)] / fish.hauls$area.swept

# Get groups for aggregating
groups <- read_csv("./data_Mar2021/Fish/Species_groups.csv")
groups <- rename(groups, species=scientificName)
groups
(taxa <- groups$Taxon[match(colnames(fish.cpue), groups$species)])
Match <- match(colnames(fish.cpue), groups$species)
if(any(is.na(Match))) stop("Taxa in catch file must match taxa in Species_groups.csv")
(taxa.fine <-groups$Taxon.fine[match(colnames(fish.cpue), groups$species)])

# Aggregate groups:
# Select grouping level ('taxa' or 'taxa.fine')
fish.cpue <- t(apply(fish.cpue, 1, function(x) tapply(x, taxa.fine, sum)))

# Eliminate juveniles and 'other' from analysis:
fish.cpue <- fish.cpue[, colnames(fish.cpue) != "Juvenile" & colnames(fish.cpue) != "other"]

# Match cruise name to other data:
fish.hauls <- mutate(fish.hauls, cruise = str_replace(cruise, "AMBON", "AMBON20"))

# Add master station names and latitude / longitude:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
fish.hauls <- left_join(fish.hauls, master, by = "station")
if(any(is.na(epi.hauls$station.master))) stop("Unrecognized station name")  

rm(groups, taxa, Match, taxa.fine)
fish.hauls <- ungroup(fish.hauls)
save(fish.hauls, fish.cpue, file="./data_Mar2021/Fish/fish.RData")

