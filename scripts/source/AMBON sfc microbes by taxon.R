######################################################################

# Import / prepare 2015 & 2017 microbe data: 
# Last update: MARCH 20, 2021

# results saved as 'microbes.RData'

######################################################################

require(tidyverse)

# Extract presence data for bacteria (prokaryotes), protists, and metazoans:
# in SURFACE waters only
mic <- read_rds("./data_Mar2021/Microbes/ambon_microbes_presence.rds")
mic <- as_tibble(mic) %>%  rename(cruise=project) %>% 
  select(-side,-depth_original,-filter,-CTD_number,-station_type) %>% 
  filter(depth_type == "surf")
bacteria <- select(mic, contains("16S_prokaryote"))
names(bacteria) <- paste("Sp", 1:ncol(bacteria), sep="")
bacteria <- as.matrix(bacteria)
bacteria <- bacteria[,apply(bacteria,2,sum)>0]

# aggregate by group and safe:
bact.agg <- select(mic, contains("16S_prokaryote"))
.nam <- strsplit(names(bact.agg), ";")
.nam <- matrix(unlist(.nam), ncol=25, byrow = T)
unique(.nam[,10])   # Aggregate at this level?
grps <- .nam[,10]
names(bact.agg) <- paste("Sp", 1:ncol(bact.agg), sep="")
bact.agg <- as.matrix(bact.agg)

f <- function(x) tapply(x, grps, max)
bact.agg <- as.data.frame(t(apply(bact.agg,1,f)))
bact.agg <- bact.agg[,apply(bact.agg,2,sum) > 0]


# Protists, all taxa:
protist <- select(mic, contains("18S_protist"))
names(protist) <- paste("Sp", 1:ncol(protist), sep="")
protist <- as.matrix(protist)
protist <- protist[,apply(protist,2,sum)>0]

# aggregate by group and safe:
prot.agg <- select(mic, contains("18S_protist"))
.nam <- strsplit(names(prot.agg), ";")
.nam <- matrix(unlist(.nam), ncol=10, byrow = T)
unique(.nam[,4])   # Aggregate at this level?
grps <- .nam[,4]
names(prot.agg) <- paste("Sp", 1:ncol(prot.agg), sep="")
prot.agg <- as.matrix(prot.agg)
dim(prot.agg)

f <- function(x) tapply(x, grps, max)
prot.agg <- as.data.frame(t(apply(prot.agg,1,f)))
prot.agg <- prot.agg[,apply(prot.agg,2,sum)>0] 


metazoa <- select(mic, contains("18S_metazoa"))
names(metazoa) <- paste("Sp", 1:ncol(metazoa), sep="")
metazoa <- as.matrix(metazoa)
metazoa <- metazoa[,apply(metazoa,2,sum)>0]

microbe.hdr <- select(mic, 1:11)
rm(mic)

sfc.microbe.hdr <- microbe.hdr
sfc.bacteria <- bacteria
sfc.protist <- protist
sfc.metazoa <- metazoa

rm(microbe.hdr, bacteria, protist, metazoa)

# Add master station names and latitude / longitude:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")
sfc.microbe.hdr <- left_join(sfc.microbe.hdr, master, by = "station")
if(any(is.na(sfc.microbe.hdr$station.master))) stop("Unrecognized station name")  

# There are duplicated samples at a number of stations within both years. 
# For now, use one of these only, assuming each is from a full sample
sfc.microbe.hdr <- sfc.microbe.hdr %>% mutate(ID = paste(cruise, station.master))
j <- duplicated(sfc.microbe.hdr$ID)
sfc.microbe.hdr <- sfc.microbe.hdr[!j,]
sfc.bacteria <- sfc.bacteria[!j,]
sfc.protist <- sfc.protist[!j,]
sfc.metazoa <- sfc.metazoa[!j,]
sfc.bact.agg <- bact.agg[!j,]
sfc.prot.agg <- prot.agg[!j,]

save(sfc.microbe.hdr, sfc.bacteria, sfc.protist, sfc.metazoa, sfc.bact.agg, sfc.prot.agg, file="./data_Mar2021/Microbes/sfc_microbes.RData")

