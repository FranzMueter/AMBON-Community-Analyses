#####################################################################

# Spatial maps of diversity by year and assemblage for eight 
# assemblages sampled in the Chukchi Sea in 2015 and 2017

# Author: Franz Mueter 
# Last updated: 9/15/2021

#####################################################################

source("./scripts/AMBON_Import_all_data.R")

# Load maps & labels, etc:
loc <- read.csv("./Base_Map/map_labels.csv")
load("./Base_Map/Chukchi_basemap_bathy.RData")  # loads 'CSmap_bathy', a ggplot object
load("./Base_Map/Chukchi_basemap.RData")  # loads 'CSmap', a ggplot object

# Overview basemap with bathymetry
basemap <- CSmap_bathy + 
  geom_point(data=loc, aes(longitude, latitude), size=3, col="black") +
  geom_point(data=loc, aes(longitude, latitude), size=1.6, col="white") +
  geom_text(data=loc, aes(longitude, latitude, label = label), 
            size=3, hjust=0, nudge_x=c(0.4,0.2,0.2,0.2,-0.1), 
            nudge_y=c(0.03,-0.1,0,0,-0.1), color="white")
basemap

# Compute richness
sfc.microbe.hdr$R.bacteria <- apply(sfc.bacteria,1,function(x) sum(x>0))
sfc.microbe.hdr$R.protists <- apply(sfc.protist,1,function(x) sum(x>0))
sfc.microbe.hdr$R.metazoa <- apply(sfc.metazoa,1,function(x) sum(x>0))
zoop150.hauls$R <- apply(zoop150.cpue, 1, function(x) sum(x>0))
zoop505.hauls$R <- apply(zoop505.cpue, 1, function(x) sum(x>0))
infauna.hauls$R <- apply(infauna.cpue, 1, function(x) sum(x>0))
epi.hauls$R <- apply(epi.cpue, 1, function(x) sum(x>0))
fish.hauls$R <- apply(fish.cpue, 1, function(x) sum(x>0))
bird.hdr$R <- apply(bird.cpue, 1, function(x) sum(x>0))

################################################################################
# Adjusting for effort
################################################################################
# Bird richness varies with effort (transect length), hence adjust for effort:
library(mgcv)
# fit1 <- gam(R ~ s(tLength, k=4, by=factor(cruise)), data=bird.hdr, family=poisson())
# fit2 <- gam(R ~ factor(cruise) + s(tLength, k=4), data=bird.hdr, family=poisson())
fit3 <- gam(R ~ s(tLength, k=4), data=bird.hdr, family=poisson())
#AIC(fit1, fit2, fit3) # fit 3 has lowest AIC!
# summary(fit3)
# plot(fit3)
# gam.check(fit3)
mu <- predict(fit3, newdata=data.frame(tLength = mean(bird.hdr$tLength)), type="response")
bird.hdr$R <- round(as.vector(mu) + resid(fit3, type="response"))

# Fish richness does not vary with effort (area swept), hence no adjustment needed:
# fit1 <- gam(R ~ s(area.swept, k=4, by=factor(cruise)), data=fish.hauls, family=poisson())
# summary(fit1)
# fit2 <- gam(R ~ factor(cruise) + s(area.swept, k=4), data=fish.hauls, family=poisson())
# summary(fit2)

# Effort adjustment for epifauna?
epi.hauls <- left_join(epi.hauls, select(fish.hauls, cruise, station, area.swept))
fit1 <- gam(R ~ s(area.swept, k=4, by=factor(cruise)), data=na.omit(epi.hauls), family=poisson())
# fit2 <- gam(R ~ factor(cruise) + s(area.swept, k=4), data=na.omit(epi.hauls), family=poisson())
# fit3 <- gam(R ~ s(area.swept, k=4), data=na.omit(epi.hauls), family=poisson())
# fit4 <- gam(R ~ 1, data=na.omit(epi.hauls), family=poisson())
# AIC(fit1, fit2, fit3, fit4) # fit 1 has lowest AIC by far!
summary(fit1)
plot(fit1)
gam.check(fit1)
# Richness increases significantly with effort in 2015, but by less than 1
# over the range of area.swept values observed in 2015. There is 
# no significant effect for 2017. Because the effect is small
# relative to the average richness (34.2 in 2015; 27.7 in 2017)
# and was not evident in 2017, I did not make any adjustments.

# Effort for infauna (number of grabs)
# fit1 <- gam(R ~ s(n.grab, k=4, by=factor(cruise)), data=infauna.hauls, family=poisson())
# fit2 <- gam(R ~ factor(cruise) + s(n.grab, k=4), data=infauna.hauls, family=poisson())
fit3 <- gam(R ~ s(n.grab, k=4), data=infauna.hauls, family=poisson())
# fit4 <- gam(R ~ 1, data=infauna.hauls, family=poisson())
# AIC(fit1, fit2, fit3, fit4) # fit 1 has lowest AIC by far!
summary(fit3)
plot(fit3)
gam.check(fit3)
# Species richness is marginally smaller on average when only one
# grab was analyzed, and the effect is negligible (< 0.4 species)
# relative to the mean number of species per haul (54)

# Effort for other groups is standardized (microbes) or reasonably 
# constant across stations (zooplankton)


######################################################################################
### Maps by group and year

ps <- 2 # Point size for color circles

# Theme for diversity maps:
my.theme1 <- function()  theme(legend.position = c(0.93,0.25), 
                               legend.key.size = unit(0.28,'cm'),
                               strip.text = element_blank(),
                               plot.margin=unit(c(0,0.16,0,0), "cm"))

# Theme for difference plots:
my.theme <- function()   theme(legend.position = c(0.82, 0.25), 
                               legend.key.size = unit(0.28,'cm'),
                               strip.text = element_blank(),
                               axis.text.y = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title.y.right = element_text(size=14),
                               plot.margin=unit(c(0,0,0,0.05), "cm"))


## Bacteria by year (2015 on left):
pl.bact <- CSmap_bathy + geom_point(aes(longitude.master, latitude.master, color=R.bacteria), size=ps, 
                                data=sfc.microbe.hdr) + facet_wrap(~cruise) +
  scale_color_viridis_c() +
  labs(color=NULL) +
  my.theme1()

## Protists by year
pl.prot <- CSmap_bathy + geom_point(aes(longitude.master, latitude.master, color=R.protists), size=ps, 
                                    data=sfc.microbe.hdr) + facet_wrap(~cruise) +
  scale_color_viridis_c() +
  labs(color=NULL) +
  my.theme1()

## Differences between years for stations sampled both years:
j <- sfc.microbe.hdr$cruise == "AMBON2015"
st <- intersect(sfc.microbe.hdr$station.master[j],sfc.microbe.hdr$station.master[!j])
dat <- sfc.microbe.hdr %>% 
  select(cruise, station.master, longitude.master, latitude.master, R.bacteria, R.protists) %>%
  filter(station.master %in% st) 

table(dat$cruise, dat$station.master)

dat <- split(dat, dat$cruise)
# Select first sample only:
dat$AMBON2015 <- dat$AMBON2015[!duplicated(dat$AMBON2015$station.master),]
dat$AMBON2017 <- dat$AMBON2017[!duplicated(dat$AMBON2017$station.master),]

# compute differences
del.mic <- left_join(dat$AMBON2015, 
                     select(dat$AMBON2017,-longitude.master,-latitude.master), 
                     by = "station.master") %>% 
  mutate(del.bacteria = R.bacteria.y - R.bacteria.x, 
         del.protists = R.protists.y - R.protists.x) %>% 
  select(station.master, longitude.master, latitude.master, del.bacteria, del.protists)

# Difference map, bacteria
pl.bact.diff <- CSmap_bathy + 
  geom_point(aes(longitude.master, latitude.master, color=del.bacteria), data=del.mic, size=ps) +
  geom_point(aes(longitude.master, latitude.master), data=del.mic, size=ps, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue",
                        mid = "white",
                        high = "red") +
  labs(color=NULL) +
  ylab("Bacteria") +
  scale_y_continuous(position = "right") +
  my.theme()

# Difference map, protists
pl.prot.diff <- CSmap_bathy + 
  geom_point(aes(longitude.master, latitude.master, color=del.protists), data=del.mic, size=ps) +
  geom_point(aes(longitude.master, latitude.master), data=del.mic, size=ps, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue",
                        mid = "white",
                        high = "red") +
  labs(color=NULL) +
  ylab("Protists") +
  scale_y_continuous(position = "right") +
  my.theme()


##############################################################################
## Zooplankton (150u) by year:
# Maps by year
pl.zoop150 <- CSmap_bathy + geom_point(aes(longitude.master, latitude.master, color=R), size=ps, 
                                       data=zoop150.hauls) + facet_wrap(~cruise) +
  scale_color_viridis_c() + labs(color=NULL) + my.theme1()

## Differences between years for stations sampled both years:
j <- zoop150.hauls$cruise == "AMBON2015"
st <- intersect(zoop150.hauls$station.master[j],zoop150.hauls$station.master[!j])
dat <- zoop150.hauls %>% 
  select(cruise, station.master, longitude.master, latitude.master, R) %>%
  filter(station.master %in% st) 

# Check for duplicate stations in a year
table(dat$cruise, dat$station.master)

# Some intersecting stations were sampled twice  
# take average richness:
dat <- dat %>% group_by(cruise,station.master) %>% 
  summarize(across(everything(), mean)) %>%  ungroup()
    
dat <- split(dat, dat$cruise)

# compute differences
del <- left_join(dat$AMBON2015,
                 select(dat$AMBON2017,-longitude.master,-latitude.master), 
                 by = "station.master") %>% 
  mutate(del.R = R.y - R.x) %>% 
  select(station.master, longitude.master, latitude.master, del.R)

# Difference map, small zooplankton
pl.zoop150.diff <- CSmap_bathy + 
  geom_point(aes(longitude.master, latitude.master, color=del.R), data=del, size=ps) +
  geom_point(aes(longitude.master, latitude.master), data=del, size=ps, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue",
                        mid = "white",
                        high = "red") +
  labs(color=NULL) + ylab("Small Zoopl.") + scale_y_continuous(position = "right") + my.theme()
  

##############################################################################
## Zooplankton (505u) by year:
# Maps by year
pl.zoop505 <- CSmap_bathy + geom_point(aes(longitude.master, latitude.master, color=R), size=ps, 
                                       data=zoop505.hauls) + facet_wrap(~cruise) +
  scale_color_viridis_c() +
  labs(color=NULL) +
  my.theme1()


## Differences between years for stations sampled both years:
j <- zoop505.hauls$cruise == "AMBON2015"
st <- intersect(zoop505.hauls$station.master[j],zoop505.hauls$station.master[!j])
dat <- zoop505.hauls %>% 
  select(cruise, station.master, longitude.master, latitude.master, R) %>%
  filter(station.master %in% st) 

# Check for duplicate stations in a year
any(table(dat$cruise, dat$station.master) > 1)

# Some intersecting stations were sampled twice  
# take average richness:
dat <- dat %>% group_by(cruise,station.master) %>% 
  summarize(across(everything(), mean)) %>%  ungroup()

dat <- split(dat, dat$cruise)

# compute differences
del <- left_join(dat$AMBON2015,
                 select(dat$AMBON2017,-longitude.master,-latitude.master), 
                 by = "station.master") %>% 
  mutate(del.R = R.y - R.x) %>% 
  select(station.master, longitude.master, latitude.master, del.R)

# Difference map, large zooplankton
pl.zoop505.diff <- CSmap_bathy + 
  geom_point(aes(longitude.master, latitude.master, color=del.R), data=del, size=ps) +
  geom_point(aes(longitude.master, latitude.master), data=del, size=ps, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue",
                        mid = "white",
                        high = "red") +
  labs(color=NULL) + ylab("Large Zoopl.") + scale_y_continuous(position = "right") + my.theme()


##############################################################################
## Macroinfauna by year:
# Maps by year
pl.infauna <- CSmap_bathy + geom_point(aes(longitude.master, latitude.master, color=R), size=ps, 
                                       data=infauna.hauls) + facet_wrap(~cruise) +
  scale_color_viridis_c() +
  labs(color=NULL) +
  my.theme1()

## Differences between years for stations sampled both years:
j <- infauna.hauls$cruise == "AMBON2015"
st <- intersect(infauna.hauls$station.master[j],infauna.hauls$station.master[!j])
dat <- infauna.hauls %>% 
  select(cruise, station.master, longitude.master, latitude.master, R) %>%
  filter(station.master %in% st) 

# Check for duplicate stations in a year
any(table(dat$cruise, dat$station.master) > 1)

# IF intersecting stations are sampled twice  
# take average richness:
dat <- dat %>% group_by(cruise,station.master) %>% 
  summarize(across(everything(), mean)) %>%  ungroup()

dat <- split(dat, dat$cruise)

# compute differences
del <- left_join(dat$AMBON2015,
                 select(dat$AMBON2017,-longitude.master,-latitude.master), 
                 by = "station.master") %>% 
  mutate(del.R = R.y - R.x) %>% 
  select(station.master, longitude.master, latitude.master, del.R)

# Difference map, macroinfauna
pl.infauna.diff <- CSmap_bathy + 
  geom_point(aes(longitude.master, latitude.master, color=del.R), data=del, size=ps) +
  geom_point(aes(longitude.master, latitude.master), data=del, size=ps, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue",
                        mid = "white",
                        high = "red") +
  labs(color=NULL) + ylab("Infauna") + scale_y_continuous(position = "right") + my.theme()


##############################################################################
## Epifauna by year:
# Map by year
pl.epi <- CSmap_bathy + geom_point(aes(longitude.master, latitude.master, color=R), size=ps, 
                                        data=epi.hauls) + facet_wrap(~cruise) +
  scale_color_viridis_c() +
  labs(color=NULL) +
  my.theme1()

## Differences between years for stations sampled both years:
j <- epi.hauls$cruise == "AMBON2015"
st <- intersect(epi.hauls$station.master[j],epi.hauls$station.master[!j])
dat <- epi.hauls %>% 
  select(cruise, station.master, longitude.master, latitude.master, R) %>%
  filter(station.master %in% st) 

# Check for duplicate stations in a year
any(table(dat$cruise, dat$station.master) > 1)

# IF intersecting stations are sampled twice  
# take average richness:
dat <- dat %>% group_by(cruise,station.master) %>% 
  summarize(across(everything(), mean)) %>%  ungroup()

dat <- split(dat, dat$cruise)

# compute differences
del <- left_join(dat$AMBON2015,
                 select(dat$AMBON2017,-longitude.master,-latitude.master), 
                 by = "station.master") %>% 
  mutate(del.R = R.y - R.x) %>% 
  select(station.master, longitude.master, latitude.master, del.R)

# Difference map, benthic epifauna
pl.epi.diff <- CSmap_bathy + 
  geom_point(aes(longitude.master, latitude.master, color=del.R), data=del, size=ps) +
  geom_point(aes(longitude.master, latitude.master), data=del, size=ps, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue",
                        mid = "white",
                        high = "red") +
  labs(color=NULL) + ylab("Epifauna") + scale_y_continuous(position = "right") + my.theme()


##############################################################################
## Demersal fish by year:
# Map by year
pl.fish <- CSmap_bathy + geom_point(aes(longitude.master, latitude.master, color=R), size=ps, 
                                   data=fish.hauls) + facet_wrap(~cruise) +
  scale_color_viridis_c() + labs(color=NULL) + my.theme1()

## Differences between years for stations sampled both years:
j <- fish.hauls$cruise == "AMBON2015"
st <- intersect(fish.hauls$station.master[j],fish.hauls$station.master[!j])
dat <- fish.hauls %>% 
  select(cruise, station.master, longitude.master, latitude.master, R) %>%
  filter(station.master %in% st) 

# Check for duplicate stations in a year
any(table(dat$cruise, dat$station.master) > 1)

# IF intersecting stations are sampled twice  
# take average richness:
dat <- dat %>% group_by(cruise,station.master) %>% 
  summarize(across(everything(), mean)) %>%  ungroup()

dat <- split(dat, dat$cruise)

# compute differences
del <- left_join(dat$AMBON2015,
                 select(dat$AMBON2017,-longitude.master,-latitude.master), 
                 by = "station.master") %>% 
  mutate(del.R = R.y - R.x) %>% 
  select(station.master, longitude.master, latitude.master, del.R)

# Difference map, demersal fish
pl.fish.diff <- CSmap_bathy + 
  geom_point(aes(longitude.master, latitude.master, color=del.R), data=del, size=ps) +
  geom_point(aes(longitude.master, latitude.master), data=del, size=ps, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(color=NULL) + ylab("Fish") + scale_y_continuous(position = "right") + my.theme()


##############################################################################
## Seabirds by year:
# Map by year
pl.bird <- CSmap_bathy + geom_point(aes(longitude.master, latitude.master, color=R), size=ps, 
                                    data=bird.hdr) + facet_wrap(~cruise) +
  scale_color_viridis_c() + labs(color=NULL) + my.theme1()

## Differences between years for stations sampled both years:
j <- bird.hdr$cruise == "AMBON2015"
st <- intersect(bird.hdr$station.master[j], bird.hdr$station.master[!j])
dat <- bird.hdr %>% 
  select(cruise, station.master, longitude.master, latitude.master, R) %>%
  filter(station.master %in% st) 

# Check for duplicate stations in a year
any(table(dat$cruise, dat$station.master) > 1)

# IF intersecting stations are sampled twice  
# take average richness:
dat <- dat %>% group_by(cruise,station.master) %>% 
  summarize(across(everything(), mean)) %>%  ungroup()

dat <- split(dat, dat$cruise)

# compute differences
del <- left_join(dat$AMBON2015,
                 select(dat$AMBON2017,-longitude.master,-latitude.master), 
                 by = "station.master") %>% 
  mutate(del.R = R.y - R.x) %>% 
  select(station.master, longitude.master, latitude.master, del.R)

# Difference map, seabirds
pl.bird.diff <- CSmap_bathy + 
  geom_point(aes(longitude.master, latitude.master, color=del.R), data=del, size=ps) +
  geom_point(aes(longitude.master, latitude.master), data=del, size=ps, pch=1, color=grey(0.6)) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(color=NULL) + ylab("Seabirds") + scale_y_continuous(position = "right") + my.theme()


# Plot maps by assemblage and year
require(gridExtra)
png(file="./plots/Fig4_DiversityMaps.png", width=5, height=12, units="in", res=300)
lay = rbind(c(1,1,2),
            c(3,3,4),
            c(5,5,6),
            c(7,7,8),
            c(9,9,10),
            c(11,11,12))
grid.arrange(pl.bact, pl.bact.diff, 
             pl.zoop150, pl.zoop150.diff,
             pl.infauna, pl.infauna.diff,
             pl.epi,pl.epi.diff, 
             pl.fish,pl.fish.diff,
             pl.bird,pl.bird.diff, layout_matrix = lay)
dev.off()

