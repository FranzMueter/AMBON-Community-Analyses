#####################################################################

# Maps of surface and bottom temperatures in the eastern Chukchi Sea 
# in summer 2015 and summer 2017 with differences between years

# Author: Franz Mueter 
# Last updated: 9/15/2021

#####################################################################

library(tidyverse)

load("./Base_Map/Chukchi_basemap_bathy.RData")  # loads 'CSmap_bathy', a ggplot object

# Import temperature and salinity data
ts <- read_csv("data_Mar2021/CTD/AMBON_TS.csv")

# Change station names to master station names:
master <- read_csv("./data_Mar2021/StationData/MasterStationNames.csv")

ts <- left_join(ts, master, by = "station")
ts <- na.omit(ts) # Drop a snd CTD cast at one station 

my.theme1 <- function() theme(legend.position = c(0.93,0.3), 
                              legend.key.size = unit(0.36,'cm'),
                              strip.text = element_blank())

### SST
pl.SST <- CSmap_bathy + geom_point(aes(longitude, latitude, color=SST), size=3, 
                                    data=ts) + facet_wrap(~cruise) +
  scale_color_viridis_c(option="magma") + my.theme1()

pl.BT <- CSmap_bathy + geom_point(aes(longitude, latitude, color=BT), size=3, 
                                   data=ts) + facet_wrap(~cruise) +
  scale_color_viridis_c(option="magma") +  my.theme1()

pl.SSS <- CSmap_bathy + geom_point(aes(longitude, latitude, color=SSS), size=3, 
                                   data=ts) + facet_wrap(~cruise) +
  scale_color_viridis_c(option="viridis") + my.theme1()

pl.BS <- CSmap_bathy + geom_point(aes(longitude, latitude, color=BS), size=3,
                                  data=ts) + facet_wrap(~cruise) +
  scale_color_viridis_c(option="viridis") + my.theme1()


## Differences between years for stations sampled both years:
j <- ts$cruise == "AMBON2015"
st <- intersect(ts$station.master[j],ts$station.master[!j])
dat <- ts %>% select(cruise, station.master, latitude, longitude, SST, SSS, BT, BS) %>%
  filter(station.master %in% st) 

table(dat$cruise, dat$station.master)

dat <- split(dat, dat$cruise)
# Select first sample only:
dat$AMBON2015 <- dat$AMBON2015[!duplicated(dat$AMBON2015$station.master),]
dat$AMBON2017 <- dat$AMBON2017[!duplicated(dat$AMBON2017$station.master),]

del <- left_join(dat$AMBON2015, 
                     select(dat$AMBON2017,-longitude,-latitude), 
                     by = "station.master") %>% 
  mutate(del.SST = SST.y - SST.x, 
         del.SSS = SSS.y - SSS.x,
         del.BT = BT.y - BT.x, 
         del.BS = BS.y - BS.x) %>% 
  select(station.master, longitude, latitude, del.SST, del.SSS, del.BT, del.BS)

my.theme2 <- function () theme(legend.position = c(0.83,0.3), 
                               legend.key.size = unit(0.36,'cm'),
                               strip.text = element_blank(), 
                               axis.ticks = element_blank())

pl.SST.diff <- CSmap_bathy + 
  geom_point(aes(longitude, latitude, color=del.SST), data=del, size=3) +
  geom_point(aes(longitude, latitude), data=del, size=3, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(color="d SST") + my.theme2()

pl.BT.diff <- CSmap_bathy + 
  geom_point(aes(longitude, latitude, color=del.BT), data=del, size=3) +
  geom_point(aes(longitude, latitude), data=del, size=3, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(color="d BT") + my.theme2()

pl.SSS.diff <- CSmap_bathy + 
  geom_point(aes(longitude, latitude, color=del.SSS), data=del, size=3) +
  geom_point(aes(longitude, latitude), data=del, size=3, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(color="d SSS") + my.theme2()

pl.BS.diff <- CSmap_bathy + 
  geom_point(aes(longitude, latitude, color=del.BS), data=del, size=3) +
  geom_point(aes(longitude, latitude), data=del, size=3, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(color="d BS") + my.theme2()

library(gridExtra)
library(grid)
png(file="./plots/Fig2_TS_maps.png", width=7, height=9.5, units="in", res=300)
lay = rbind(c(1,1,2),
            c(3,3,4),
            c(5,5,6),
            c(7,7,8))
grid.arrange(pl.SST, pl.SST.diff, 
             pl.BT, pl.BT.diff,
             pl.SSS, pl.SSS.diff,
             pl.BS,pl.BS.diff, 
             layout_matrix = lay,
             top = textGrob("               2015                2017                        (2017-2015)", 
                                gp=gpar(fontsize=18)))
dev.off()


