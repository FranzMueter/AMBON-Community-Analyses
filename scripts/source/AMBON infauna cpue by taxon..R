######################################################################

# Import / prepare 2015 & 2017 haul and cpue data for macroinfauna: 

# Author: Franz Mueter
# Last update: MARCH 19, 2021

# CPUE (biomass) and haul info written to 'infauna.RData'

######################################################################

require(tidyverse)

## Compute CPUE by number and year

# Import haul / environmental data from Jackie Grebmeier
hauls15 <- read_csv("./data_Mar2021/Macroinfauna/AMBON15_env.csv")
hauls17 <- read_csv("./data_Mar2021/Macroinfauna/AMBON17_env.csv")

hauls <- bind_rows(mutate(hauls15, cruise="AMB15"), mutate(hauls17, cruise="AMB17")) %>% 
                     select(-DataDate, -GMT)

# Import species data:
inf15 <- read_csv("./data_Mar2021/Macroinfauna/AMBON2015_SPECIES_Macroinfaunal_taxa.csv") %>% 
  select(cruise = Cruise,
         station = 'Station Name',
         species = Species,
         cpue = 'Wet Weight',
         depth = Depth,
         latitude = Latitude,
         longitude = Longitude,
         grabs = 'No of Grab')

# remove taxa with "X"
inf15 <- inf15[-grep("(X)", inf15$species, fixed=T),]

inf17 <- read_csv("./data_Mar2021/Macroinfauna/AMBON2017_SPECIES_Macroinfauna_taxa.csv") %>% 
  select(cruise = Cruise,
         station = 'Station Name',
         species = Species,
         cpue = 'Wet Weight',
         depth = Depth,
         latitude = Latitude,
         longitude = Longitude,
         grabs = 'No of Grab')  

# remove taxa with "X"
inf17 <- inf17[-grep("(X)", inf17$species, fixed=T),]

taxa15 <- inf15 %>% select(species) %>% group_by(species) %>% summarize(N15=n())
taxa17 <- inf17 %>% select(species) %>% group_by(species) %>% summarize(N17=n())
write.table(full_join(taxa15, taxa17), file="./data_Mar2021/Macroinfauna/prelim_taxa.csv", 
            row.names = F, sep=",")

# Import taxa names and add to data frames:
taxa <- read_csv("./data_Mar2021/Macroinfauna/infauna_taxa.csv")
inf15 <- left_join(inf15, taxa)
inf17 <- left_join(inf17, taxa)

# Drop 'other' and unidentified taxa:
infauna <- rbind(inf15, inf17) %>% filter(taxon != "other.unid")

### Get cpue by taxon for each cruise & station
infauna.cpue <- infauna %>% group_by(cruise, station, taxon) %>% 
  summarize(N=sum(cpue)) %>% 
  pivot_wider(id_cols = c("cruise", "station", "taxon"), 
              names_from = taxon, values_from = N, values_fill=0)

# Get corresponding hauls
infauna.hauls <- select(infauna.cpue, 1:2) %>%  left_join(hauls) 

# Match cruise name to other data:
infauna.hauls <- mutate(infauna.hauls, cruise = str_replace(cruise, "AMB", "AMBON20"))

rm(taxa,inf15,inf17,hauls15,hauls17)

# Add master station names and latitude / longitude:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
infauna.hauls <- left_join(infauna.hauls, master, by = "station")
if(any(is.na(epi.hauls$station.master))) stop("Unrecognized station name")  

# Save cpue as matrix for analysis:
infauna.cpue <- as.matrix(infauna.cpue[,-c(1:2)])
infauna.hauls <- ungroup(infauna.hauls)
save(infauna.hauls, infauna.cpue, file="./data_Mar2021/Macroinfauna/infauna.RData")

