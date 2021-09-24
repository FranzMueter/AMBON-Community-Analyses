library(tidyverse)
library(readxl)

# Script to import data from multiple assemblages and line up data for a common
# set of stations that were sampled in 2015 
# This imports only 2015 data!

# Import haul data and catch data and compute CPUE by taxon: 
setwd("C:/Dropbox/Research/Projects/AMBON/Analysis")

# Read in numeric CPUE data for PSBTA to estimate species richness,
# mean CPUE for selected species, etc.
source("AMBON_CPUE_by_taxon.R")   # Reads in all data (2015, 2017)
only15 <- PSBTA.hauls$ Cruise == "AMBON 2015"
PSBTA.fish.cpue <- cpue.nr[only15,]
rm(cpue.nr)
PSBTA.hauls <- PSBTA.hauls[only15,]
rm(only15)

# Take average of CPUE at stations with more than one quantiative haul:
if(any(duplicated(PSBTA.hauls$locationID))) {
  PSBTA.fish.cpue <- aggregate(PSBTA.fish.cpue, by = list(Station=PSBTA.hauls$locationID), mean)
  rownames(PSBTA.fish.cpue) <- PSBTA.fish.cpue$Station
  PSBTA.fish.cpue <- PSBTA.fish.cpue[,-1]
  # Drop duplicated haul - note that distance towed may differ between hauls, hence use mean
  # of duplicated hauls at a station:
  d <- tapply(PSBTA.hauls$Distance, list(Station=PSBTA.hauls$locationID), sum)
  PSBTA.hauls <- PSBTA.hauls[!duplicated(PSBTA.hauls$locationID),]
  PSBTA.hauls$Distance <- d
} else rownames(PSBTA.fish.cpue) <- PSBTA.hauls$locationID
PSBTA.fish.cpue[,1:2]

# Check if station names match:
if(!identical(rownames(PSBTA.fish.cpue), PSBTA.hauls$locationID)) stop("Station names don't match")

# Eliminate species not caught in 2015:
j <- apply(PSBTA.fish.cpue,2,sum) != 0
PSBTA.fish.cpue <- PSBTA.fish.cpue[,j]

rm(d,j)

# Use resulting data matrix for 2015 community analysis

###########################################################################################
# Import epibenthic invertebrate data:
setwd("C:/Dropbox/Research/Projects/AMBON/Data/Epifauna")
epi <- read.csv("AMBON2015_epifauna_abundance_database.csv")
str(epi)
table(epi$Station)

# Re-name station 'ML6-10' to station 'ML4-4' (as per Master List):
epi$Station <- as.character(epi$Station)
epi$Station[epi$Station=="ML6-10"] <- "ML4-4"

# Rename long, complex names
epi <- rename(epi, 
              cpue = "Abundance..ind.1000.m.2.",
              species = "Scientific.name")
# Cross-tabulate by station and species:
epifauna.cpue <- xtabs(cpue ~ Station + species, data=epi)
# Change station names to match fish stations:
rownames(epifauna.cpue) <- gsub("-",".",rownames(epifauna.cpue))

# Check if stations (rows) in epifauna data match stations in fish data:
identical(rownames(epifauna.cpue), PSBTA.hauls$locationID)
dimnames(epifauna.cpue)

# Check frequence of occurrence and total number:
tot <- apply(epifauna.cpue,2,sum)
fo <- apply(epifauna.cpue,2,function(x) sum(x>0))
# Write data to examine and combine into groups(?):
write.csv(cbind(tot,fo), "epifauna_taxa.csv")
rm(tot,fo)

# Import file with aggregated taxa for analysis:
# (I grouped these based on my own judgement in May 2019,
#  perhaps run by Katrin for confirmation or changes)
epi.taxa <- read_xlsx("epifauna_taxa.xlsx", "epifauna_taxa")

epifauna.cpue <- t(apply(epifauna.cpue, 1, function(x) tapply(x, list(epi.taxa$Taxa), sum)))

# Check if stations (rows) in epifauna data still match stations in fish data:
if(!identical(rownames(epifauna.cpue), PSBTA.hauls$locationID)) stop("Stations must match")

# Remove 'other/unid'd' group:
epifauna.cpue <- epifauna.cpue[,dimnames(epifauna.cpue)[[2]] != "other.unid"]
dim(epifauna.cpue)

rm(epi)


###########################################################################################
# Import zooplankton data:
# Downloaded data from Workspace in May 2019 and determined groups for aggregating / analysis 
# with Russ's help. Did grouping in Excel and saved file:
setwd("C:/Dropbox/Research/Projects/AMBON/Data/Zooplankton")
zoop.cpue <- read.csv("AMBON2015_zooplankton_150.csv", row.names=1)
dim(zoop.cpue)

# Re-name stations to match up with fish/invert stations:
#  'ML6.10' to 'ML4.4':
row.names(zoop.cpue)[row.names(zoop.cpue) =="ML6.10"] <- "ML4.4"
#  'ML1.3' to  'ML6.2':
row.names(zoop.cpue)[row.names(zoop.cpue) =="ML1.3"] <- "ML6.2"


#  'ML3.4' is same as 'ML6.2', both appear in Zoops data
# I removed ML3.4 because event log suggests no PSBTA taken at that time
zoop.cpue <- zoop.cpue[row.names(zoop.cpue) != "ML3.4",]

#  'ML5.8' is same as  'ML3.8':
# I removed ML5.8 because event log suggests no PSBTA taken at that time
zoop.cpue <- zoop.cpue[row.names(zoop.cpue) != "ML5.8",]

# Finally, I removed ML1.10 and ML4.3 because no quantiatitive PSBTA was taken...
no.psbta <- row.names(zoop.cpue)[is.na(match(row.names(zoop.cpue), PSBTA.hauls$locationID))]
zoop.cpue <- zoop.cpue[!(row.names(zoop.cpue) %in% no.psbta),]

# ... and removed "DBO3.1" and "DBO3.3" from fish and epibentic inverts 
#      because no zoops were collected:
no.zoop <- PSBTA.hauls$locationID[is.na(match(PSBTA.hauls$locationID, row.names(zoop.cpue)))]
epifauna.cpue <- epifauna.cpue[!(row.names(epifauna.cpue) %in% no.zoop),]
PSBTA.fish.cpue <- PSBTA.fish.cpue[!(row.names(PSBTA.fish.cpue) %in% no.zoop),]
PSBTA.hauls <- PSBTA.hauls[!(PSBTA.hauls$locationID %in% no.zoop),]

# Eliminate any epifauna taxa that don't occur at any of the selected stations:
epifauna.cpue <- epifauna.cpue[,apply(epifauna.cpue,2,sum)!=0]

rm(no.psbta, no.zoop)

# Re-order alphabetically:
zoop.cpue <- zoop.cpue[order(row.names(zoop.cpue)),]
if(!identical(row.names(zoop.cpue), row.names(PSBTA.fish.cpue))) stop("Stations must match")
if(!identical(row.names(zoop.cpue), row.names(epifauna.cpue))) stop("Stations must match")

# Eliminate 'other' group:
zoop.cpue <- zoop.cpue[,names(zoop.cpue) !="other"]

# convert to matrix:
zoop.cpue <- as.matrix(zoop.cpue)

###########################################################################################
# Import infaunal data:
# Start with Grab 1 from Jackie for 2015 stations as per e-mail from May 2019
setwd("C:/Dropbox/Research/Projects/AMBON/Data/Infauna")

infauna <- read_xlsx("AMBON2015_Infauna_Abundance_Grab1.xlsx", "Infauna")
str(infauna)

# Replace '-' in station names with '.'
infauna$Station <- gsub("-",".",infauna$Station)

# Cross-tabulate by station and species:
infauna.cpue <- xtabs(cpue ~ Station + Species, data=infauna)

# Re-name stations to match up with fish/invert stations:
#  'ML6.10' to 'ML4.4':
row.names(infauna.cpue)[row.names(infauna.cpue) =="ML6.10"] <- "ML4.4"

# Delete stations without data for other assemblages (PSBTA.haul)
no.dat <- row.names(infauna.cpue)[is.na(match(row.names(infauna.cpue), PSBTA.hauls$locationID))]
infauna.cpue <- infauna.cpue[!(row.names(infauna.cpue) %in% no.dat),]

infauna.cpue <- infauna.cpue[order(row.names(infauna.cpue)),]

# Check if stations (rows) in epifauna data match stations in fish data:
identical(rownames(infauna.cpue), PSBTA.hauls$locationID)
dimnames(infauna.cpue)

# Check frequence of occurrence and total number:
tot <- apply(infauna.cpue,2,sum)
fo <- apply(infauna.cpue,2,function(x) sum(x>0))
# Write data to examine and combine into groups(?):
write.csv(cbind(tot,fo), "infauna_taxa.csv")
rm(tot,fo)


# Import file with aggregated taxa for analysis:
# (I grouped these based on my own judgement in May 2019,
#  perhaps run by Jackie for confirmation or changes)
in.taxa <- read_xlsx("infauna_taxa.xlsx", "infauna_taxa")

infauna.cpue <- t(apply(infauna.cpue, 1, function(x) tapply(x, list(in.taxa$Species), sum)))

# Remove 'other/unid'd' group:
infauna.cpue <- infauna.cpue[,dimnames(infauna.cpue)[[2]] != "other.unid"]
dim(infauna.cpue)

# Check if stations (rows) in epifauna data still match stations in fish data:
if(!identical(rownames(infauna.cpue), PSBTA.hauls$locationID)) stop("Stations must match")

# Eliminate species not caught in 2015:
j <- apply(infauna.cpue,2,sum) != 0
infauna.cpue <- infauna.cpue[,j]
rm(j, no.dat, infauna)


###########################################################################################
# Import seabird data:
setwd("C:/Dropbox/Research/Projects/AMBON/Data/SeaBirds")

birds <- read_csv("AMBON_Seabird_Station_Densities_Apportioned.csv")
str(birds)
birds <- base::subset(birds, Year == "2015")

# Re-name stations to match up with fish/invert stations:
#  'ML1.3' to 'ML6.2'  and  'ML3.4' to 'ML6.7':
birds$Station_Name[birds$Station_Name =="ML1.3"] <- "ML6.2"
birds$Station_Name[birds$Station_Name =="ML3.4"] <- "ML6.7"

birds <- left_join(PSBTA.hauls[,3:5], birds, by=c("locationID" = "Station_Name"))
dim(birds)

birds.hdr <- birds[,1:6]
rownames(birds)
birds.dens <- as.matrix(birds[,7:31])
birds.dens <- birds.dens[,apply(birds.dens, 2, sum) != 0]
dimnames(birds.dens)[[1]] <- birds.hdr$locationID

############################################################################################
# Import 2015 sediment data:
setwd("C:/Dropbox/Research/Projects/AMBON/Data/Sediment")
sed <- read_csv("AMBON2015_Surface_Sediment Parameters.csv")
sed

# Stations that are not matched in PSBTA.hauls:
j <- match(sed$'Station Name', PSBTA.hauls$locationID)
sed$'Station Name'[is.na(j)]

# Select sediment data corresponding to data in PSBTA hauls:
sediment <- left_join(as_tibble(PSBTA.hauls[,3:5]), sed, by=c("locationID" = "Station Name"))
str(sediment)
identical(sediment$locationID, PSBTA.hauls$locationID)

if(any(sediment$Latitude.x-sediment$Latitude.y>0.01 & sediment$Longitude.x-sediment$Longitude.y>0.015)) {
  stop("Sediment stations too far from fish haul?")
}

rm(sed, j)

