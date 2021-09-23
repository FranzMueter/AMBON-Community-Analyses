# Basemap for Northeast Chukchi Sea with bathymetry

require(tidyverse)

# Load coastline and contour lines:
load("./Base_Map/akmap.RData")
load("./Base_Map/Chukchi_contours.RData")
# Contours created by 'Research/Data/Bathymetry/GetContourLines.R'

load("./Base_Map/Chukchi_basemap.RData")

## Map with contours:
Z <- contour.lines
CSmap_bathy <- CSmap +
  geom_path(data=Z$z500, aes(x,y), linetype=1, col="grey60") +
  geom_path(data=Z$z200, aes(x,y), linetype=2, col="grey60") +
  geom_path(data=Z$z100, aes(x,y), linetype=4, col="grey60") +
  geom_path(data=Z$z40, aes(x,y), linetype=3, col="grey60") 

save(CSmap_bathy, file = "Chukchi_basemap_bathy.RData")

