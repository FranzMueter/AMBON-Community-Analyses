######################################################################

# Import / prepare 2015 & 2017 haul and cpue data or epifauna: 
# Last update: MARCH 19, 2021

# CPUE (biomass) and haul info written to 'epifauna.RData'

######################################################################

require(tidyverse)

## Compute CPUE by number and year

# Import taxonomic groups:
taxa <- read_csv("./data_Mar2021/Epifauna/epifauna_taxa.csv")
  
# Import 2015 specimen data:
epi15 <- read_csv("./data_Mar2021/Epifauna/AMBON2015_epifauna_biomass_DWC.csv") %>% 
  select(station = eventID,
         species = scientificName,
         cpue = organismQuantity,
         depth = maximumDepthinMeters,
         latitude = decimalLatitude,
         longitude = decimalLongitude) %>% 
  mutate(cruise = "AMBON 2015") %>% 
  left_join(select(taxa, species=species2015, taxon))

# Import 2017 specimen data:
epi17 <- read_csv("./data_Mar2021/Epifauna/AMBON2017_epifauna_biomass_DWC.csv") %>% 
  select(station = eventID,
         species = ScientificName_accepted,
         cpue = organismQuantity,
         depth = maximumDepthinMeters,
         latitude = decimalLatitude,
         longitude = decimalLongitude) %>% 
  mutate(cruise = "AMBON 2017") %>% 
  left_join(select(taxa, species=species2017, taxon))

# Drop 'other' and unidentified taxa:
epi <- rbind(epi15, epi17) %>% filter(taxon != "other.unid")

### Get cpue by species for each cruise & station
epi.cpue <- epi %>% group_by(cruise, station, taxon) %>% 
  summarize(N=sum(cpue)) %>% 
  pivot_wider(id_cols = c("cruise", "station", "taxon"), 
              names_from = taxon, values_from = N, values_fill=list(N=0)) %>% 
  ungroup()

# Extract haul IDs
epi.hauls <- select(epi.cpue, 1:2) %>% 
  mutate(cruise = str_replace(cruise, "AMBON 2", "AMBON2"))

# Add master station names and latitude / longitude:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
epi.hauls <- left_join(epi.hauls, master, by = "station")
if(any(is.na(epi.hauls$station.master))) stop("Unrecognized station name")  

# Save cpue as matrix for analysis:
epi.cpue <- as.matrix(epi.cpue[,-c(1:2)])
rm(taxa, epi15, epi17)
save(epi.hauls, epi.cpue, file="./data_Mar2021/Epifauna/epifauna.RData")

