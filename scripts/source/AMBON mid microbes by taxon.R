######################################################################

# Import / prepare 2015 & 2017 microbe data: 
# Last update: MARCH 20, 2021

# results written to 'mid_microbes.RData'

######################################################################

require(tidyverse)

# Extract presence data for bacteria, protists, and metazoans:
# in midwater (near Chl maximum) only
mic <- read_rds("./data_Mar2021/Microbes/ambon_microbes_presence.rds")
mic <- as_tibble(mic) %>%  rename(cruise=project) %>% 
  select(-side,-depth_original,-filter,-CTD_number,-station_type) %>% 
  filter(depth_type == "mid")
bacteria <- select(mic, contains("16S_prokaryote"))
names(bacteria) <- paste("Sp", 1:ncol(bacteria), sep="")
bacteria <- as.matrix(bacteria)
bacteria <- bacteria[,apply(bacteria,2,sum)>0]

protist <- select(mic, contains("18S_protist"))
names(protist) <- paste("Sp", 1:ncol(protist), sep="")
protist <- as.matrix(protist)
protist <- protist[,apply(protist,2,sum)>0]

metazoa <- select(mic, contains("18S_metazoa"))
names(metazoa) <- paste("Sp", 1:ncol(metazoa), sep="")
metazoa <- as.matrix(metazoa)
metazoa <- metazoa[,apply(metazoa,2,sum)>0]

microbe.hdr <- select(mic, 1:11)

mid.microbe.hdr <- microbe.hdr
mid.bacteria <- bacteria
mid.protist <- protist
mid.metazoa <- metazoa
rm(mic, microbe.hdr, bacteria, protist, metazoa)

# Add master station names and latitude / longitude:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
mid.microbe.hdr <- left_join(mid.microbe.hdr, master, by = "station")
if(any(is.na(mid.microbe.hdr$station.master))) stop("Unrecognized station name")  

# There are duplicated samples at a number of stations within both years. 
# For now, use one of these only, assuming each is from a full sample
mid.microbe.hdr <- mid.microbe.hdr %>% mutate(ID = paste(cruise, station.master))
j <- duplicated(mid.microbe.hdr$ID)
mid.microbe.hdr <- mid.microbe.hdr[!j,]
mid.bacteria <- mid.bacteria[!j,]
mid.protist <- mid.protist[!j,]
mid.metazoa <- mid.metazoa[!j,]

save(mid.microbe.hdr, mid.bacteria, mid.protist, mid.metazoa, file="./data_Mar2021/Microbes/mid_microbes.RData")

