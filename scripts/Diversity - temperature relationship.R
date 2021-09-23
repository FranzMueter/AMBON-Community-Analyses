#####################################################################

# Quantify and graphically assess:
# Temperature - diversity relationships for each assemblage by year / overall
# Salinity - diversity relationships for each assemblage by year / overall

# Author: Franz Mueter 
# Last updated: 9/17/2021

#####################################################################


# assess significance of relationships individually and plot
# linear relationships by assemblage and year
library(tidyverse)
library(corrplot)

source("./scripts/AMBON_Import_all_data.R")

# Compute species richness at each station:
sfc.microbe.hdr$R.bacteria <- apply(sfc.bacteria,1,function(x) sum(x>0))
sfc.microbe.hdr$R.protists <- apply(sfc.protist,1,function(x) sum(x>0))
sfc.microbe.hdr$R.metazoa <- apply(sfc.metazoa,1,function(x) sum(x>0))
zoop150.hauls$R <- apply(zoop150.cpue, 1, function(x) sum(x>0))
zoop505.hauls$R <- apply(zoop505.cpue, 1, function(x) sum(x>0))
infauna.hauls$R <- apply(infauna.cpue, 1, function(x) sum(x>0))
epi.hauls$R <- apply(epi.cpue, 1, function(x) sum(x>0))
fish.hauls$R <- apply(fish.cpue, 1, function(x) sum(x>0))
bird.hdr$R <- apply(bird.cpue, 1, function(x) sum(x>0))

# Add environmental data for each assemblage:
env <- select(sfc.microbe.hdr, cruise, station.master, R.bacteria, R.protists) %>% 
  left_join(all.env) %>% 
  select(cruise, R.bacteria, R.protists, latitude = latitude.master, Temperature = SST, BT, Salinity=SSS, Bsal, int.chla, sand, strat:runoff.sfc)
corrplot(cor(na.omit(env[,-1])))
bact <- na.omit(select(env, cruise, R=R.bacteria, Temperature, Salinity, latitude))
prot <- na.omit(select(env, cruise, R=R.protists, Temperature, Salinity, latitude))
bact <- mutate(bact, group = "Bacteria")
prot <- mutate(prot, group = "Protists")

# Bacteria, temperature:
fit1 <- lm(R~cruise/Temperature, data = bact)
fit2 <- lm(R~Temperature + cruise, data = bact)
AIC(fit1,fit2)
anova(fit1)
summary(fit1)
# Non-ignificant difference in slopes
# However, variances are quite different between 2015 and 2017 with much higher
# variance in richness across stations in 2015. Therefore, I fit the models 
# separately by year to test for signficant relationships with temperature:
summary(lm(R~Temperature, data = bact, subset = cruise=="AMBON2015"))
## No significant effect in 2015
summary(lm(R~Temperature, data = bact, subset = cruise=="AMBON2017"))
## Highly significant, positive temperature effect in 2017

# Bacteria, salinity:
fit1 <- lm(R~cruise/Salinity, data = bact)
fit2 <- lm(R~Salinity + cruise, data = bact)
fit3 <- lm(R~cruise, data = bact)
AIC(fit1,fit2,fit3)
anova(fit1)
anova(fit2)
summary(fit1)
# No apparent salinity effect

# Protists, temperature: 
fit1 <- lm(R~cruise/Temperature, data = prot)
fit2 <- lm(R~Temperature + cruise, data = prot)
AIC(fit1,fit2)
anova(fit1)
summary(fit1)
summary(fit2)
# No significant interaction
# Significant slope

# Protists, salinity:
fit1 <- lm(R~cruise/Salinity, data = prot)
fit2 <- lm(R~Salinity + cruise, data = prot)
fit3 <- lm(R~cruise, data = prot)
AIC(fit1,fit2,fit3)
anova(fit1)
anova(fit2)
summary(fit2)
# No significant salinity effect

# Small Zooplankton, temperature
env <- select(zoop150.hauls, cruise, station.master, R) %>% 
  left_join(all.env) %>% 
  select(cruise, R, latitude = latitude.master, Temperature = SST, BT, Salinity=SSS, Bsal, int.chla, strat:runoff.sfc, sand)
corrplot(cor(na.omit(env[,-1])))
z150 <- na.omit(select(env, cruise, R, Temperature, Salinity, latitude)) %>% 
  mutate(group = "Sm. Zoopl.")
fit1 <- lm(R~cruise/Temperature, data = z150)
fit2 <- lm(R~Temperature + cruise, data = z150)
AIC(fit1,fit2)
anova(fit1)
summary(fit1)
summary(fit2)
# Significant interaction
# significant slope in 2017

# Small Zooplankton, salinity
fit1 <- lm(R~cruise/Salinity, data = z150)
fit2 <- lm(R~Salinity + cruise, data = z150)
fit3 <- lm(R~cruise, data = z150)
AIC(fit1,fit2,fit3)
anova(fit1)
anova(fit2)
summary(fit2)
# No significant salinity effect

# Large Zooplankton, temperature
env <- select(zoop505.hauls, cruise, station.master, R) %>% 
  left_join(all.env) %>% 
  select(cruise, R, latitude = latitude.master, Temperature = SST, BT, Salinity=SSS, Bsal, int.chla, strat:runoff.sfc, sand)
corrplot(cor(na.omit(env[,-1])))
z505 <- na.omit(select(env, cruise, R, Temperature, Salinity, latitude)) %>% 
  mutate(group = "Lg. Zoopl.")
fit1 <- lm(R~cruise/Temperature, data = z505)
fit2 <- lm(R~Temperature + cruise, data = z505)
AIC(fit1,fit2)
anova(fit1)
summary(fit1)
summary(fit2)
# Significant interaction
# Significant slope in 2017

# Large Zooplankton, salinity
fit1 <- lm(R~cruise/Salinity, data = z505)
fit2 <- lm(R~Salinity + cruise, data = z505)
fit3 <- lm(R~cruise, data = z505)
AIC(fit1,fit2,fit3)
anova(fit1)
anova(fit2)
summary(fit2)
# No significant salinity effect

# Infauna, temperature
env <- select(infauna.hauls, cruise, station.master, R) %>% 
  left_join(all.env) %>% 
  select(cruise, R, latitude = latitude.master, Temperature = BT, SST, SSS, Salinity=Bsal, int.chla, strat:runoff.sfc, sand)
corrplot(cor(na.omit(env[,-1])))
inf <- na.omit(select(env, cruise, R, Temperature, Salinity, latitude)) %>% 
  mutate(group = "Macroinfauna")
fit1 <- lm(R~cruise*Temperature, data = inf)
fit2 <- lm(R~Temperature + cruise, data = inf)
AIC(fit1,fit2)
anova(fit1)
summary(fit1)
summary(fit2)
# Non-significant interaction
# Significant negative slope in both years

# Infauna, salinity
fit1 <- lm(R~cruise*Salinity, data = inf)
fit2 <- lm(R~Salinity + cruise, data = inf)
fit3 <- lm(R~cruise, data = inf)
AIC(fit1,fit2,fit3)
anova(fit1)
anova(fit2)
summary(fit2)
# No significant interaction
# Significant, positive salinity effect
fit4 <- lm(R ~ Salinity + Temperature + cruise, data = inf)
summary(fit4)
fit5 <- lm(R ~ Salinity * Temperature + cruise, data = inf)
summary(fit5)
AIC(fit4, fit5)
# Salinity-temperature interaction not significant
# Salinity not significant in joint model!


# Epifauna, temperature
env <- select(epi.hauls, cruise, station.master, R) %>% 
  left_join(all.env) %>% 
  select(cruise, R, latitude = latitude.master, Temperature = BT, SST, SSS, Salinity = Bsal, int.chla, strat:runoff.sfc, sand)
corrplot(cor(na.omit(env[,-1])))
epi <- na.omit(select(env, cruise, R, Temperature, Salinity, latitude)) %>% 
  mutate(group = "Epibenthos")
fit1 <- lm(R~cruise*Temperature, data = epi)
fit2 <- lm(R~Temperature + cruise, data = epi)
AIC(fit1,fit2)
anova(fit1)
summary(fit1)
summary(fit2)
# No significant interaction
# Highly significant negative slope with temperature

# epifauna, salinity
fit1 <- lm(R~cruise*Salinity, data = epi)
fit2 <- lm(R~Salinity + cruise, data = epi)
fit3 <- lm(R~cruise, data = epi)
AIC(fit1,fit2,fit3)
anova(fit1)
anova(fit2)
summary(fit2)
# No significant interaction
# Significant, positive salinity effect
fit4 <- lm(R ~ Salinity + Temperature + cruise, data = epi)
summary(fit4)
fit5 <- lm(R ~ Salinity * Temperature + cruise, data = epi)
summary(fit5)
AIC(fit4, fit5)
# Salinity-temperature interaction not significant
# Salinity not significant in joint model!

# Fish, temperature
env <- select(fish.hauls, cruise, station.master, R) %>% 
  left_join(all.env) %>% 
  select(cruise, R, latitude = latitude.master, Temperature = BT, SST, SSS, Salinity = Bsal, int.chla, strat:runoff.sfc, sand)
corrplot(cor(na.omit(env[,-1])))
fish <- na.omit(select(env, cruise, R, Temperature, Salinity, latitude)) %>% 
  mutate(group = "Demersal fish")
fit1 <- lm(R~cruise*Temperature, data = fish)
fit2 <- lm(R~Temperature + cruise, data = fish)
AIC(fit1,fit2)
anova(fit1)
summary(fit1)
summary(fit2)
# No significant interaction
# Significant, positive slope with temperature

# Fish, salinity
fit1 <- lm(R~cruise*Salinity, data = fish)
fit2 <- lm(R~Salinity + cruise, data = fish)
fit3 <- lm(R~cruise, data = fish)
AIC(fit1,fit2,fit3)
anova(fit1)
anova(fit2)
summary(fit2)
# No significant interaction
# Significant, negative salinity effect
fit4 <- lm(R ~ Salinity + Temperature + cruise, data = fish)
summary(fit4)
fit5 <- lm(R ~ Salinity * Temperature + cruise, data = fish)
summary(fit5)
AIC(fit4, fit5)
# Salinity not significant in joint model


# Birds, temperature
env <- select(bird.hdr, cruise, station.master, R) %>% 
  left_join(all.env) %>% 
  select(cruise, R, latitude = latitude.master, Temperature = SST, BT, Salinity = SSS, Bsal, int.chla, strat:runoff.sfc, sand)
corrplot(cor(na.omit(env[,-1])))
bird <- na.omit(select(env, cruise, R, Temperature, Salinity, latitude)) %>% 
  mutate(group = "Seabirds")
fit1 <- lm(R~cruise*Temperature, data = bird)
fit2 <- lm(R~Temperature + cruise, data = bird)
AIC(fit1,fit2)
anova(fit1)
summary(fit1)
summary(fit2)
# No significant interaction (but interaction model has slightly lower AIC)
# Significant, positive slope with temperature

# Birds, salinity
fit1 <- lm(R~cruise*Salinity, data = bird)
fit2 <- lm(R~Salinity + cruise, data = bird)
fit3 <- lm(R~cruise, data = bird)
AIC(fit1,fit2,fit3)
anova(fit1)
anova(fit2)
summary(fit2)
# No significant interaction
# Significant, positive salinity effect
fit4 <- lm(R ~ Salinity + Temperature + cruise, data = bird)
summary(fit4)
fit5 <- lm(R ~ Salinity * Temperature + cruise, data = bird)
summary(fit5)
AIC(fit4, fit5)
# No significant interaction between temp & salinity
# Both Salinity & Temperature significant in joint model
# (positive effects in both cases)
# Higher bird diversity associated with high salinities 
# and higher temperatures
library(visreg)
visreg(fit4)

all.R <- bind_rows(bact, prot, z150, z505, inf, epi, fish, bird)
all.R <- mutate(all.R, group = ordered(group, levels = unique(group)))

# High resolution  / size for manuscript:
png("./plots/Richness_temperature.png", width=7.5, height=4, units = "in", res=300)
ggplot(all.R, aes(Temperature, R, color = cruise)) + 
  geom_point(alpha=0.5) + 
  geom_smooth(method="lm") +
  facet_wrap(~group, scales="free_y", nrow=2) +
  theme(legend.position = "none", 
        strip.text = element_text(size=15),
        axis.title = element_text(size=15)) +
  scale_color_manual(values = c("AMBON2015" = "skyblue3", "AMBON2017" = "orangered")) +
  scale_fill_manual(values = c("AMBON2015" = "skyblue3", "AMBON2017" = "orangered")) +
  ylab("Species richness")
dev.off()

# Low resolution for salinity plot:
png("./plots/Richness_salinity.png", width=750, height=400)
ggplot(all.R, aes(Salinity, R, color = cruise)) + 
  geom_point(alpha=0.5) + 
  geom_smooth(method="lm") +
  facet_wrap(~group, scales="free_y", nrow=2) +
  theme(legend.position = "none", 
        strip.text = element_text(size=16),
        axis.title = element_text(size=16)) +
  scale_color_manual(values = c("AMBON2015" = "skyblue3", "AMBON2017" = "orangered")) +
  scale_fill_manual(values = c("AMBON2015" = "skyblue3", "AMBON2017" = "orangered")) +
  ylab("Species richness")
dev.off()

