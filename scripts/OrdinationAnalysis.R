#####################################################################
# Script to conduct NMDS ordination for eight assemblages in the 
# Northeast Chukchi Sea. The relationship between each ordination
# and a suite of explanatory (environmental) variables is tested 
# and quantified as follows:

# To reduce the chance of identifying spurious relationships, we 
# first selected the set of environmental variables that has the 
# strongest Mantel correlation with the biological data matrix for 
# each assemblage (using envfit). We then tested the significance 
# of these variables using PERMANOVA and quantified the marginal 
# contribution of each of the significant variable to overall 
# assemblage composition (marginal R^2).

# Ordination results and the output from'envfit' for each assemblage 
# are saved for plotting (see 'Ordination Env plots.R) 

# Author: Franz Mueter 
# Last updated: 9/18/2021

#####################################################################

# Load required packages
library(tidyverse)
library(mgcv)
library(vegan)
library(maps)
library(mapdata)

# Process & import data for each assemblage plus environmental data
source("./scripts/AMBON_Import_all_data.R")

############################################################################
### Select sediment parameters for analysis and reduce dimensionality (PCA):
sed <- select(all.env, cruise, station.master, phi.0:phi.5, TOC:delN15)
summary(sed)  # 3 extra missing values in 13C, hence eliminate 
# remove sand (=sum of phi1 - phi4)
sed <- select(sed, -delC13, -sand)
sed <- na.omit(sed)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(sed[,-c(1:2),], diag.panel=panel.hist)

# Transform coarse fractions, which tend to be highly right-skewed
sed <- mutate(sed, phi.0 = phi.0^0.25, phi.1 = phi.1^0.25, phi.2 = phi.2^0.25, 
              phi.3 = phi.3^0.25) 
pairs(sed[,-c(1,2)], diag.panel=panel.hist)

sed.pca <- princomp(sed[,-c(1,2)], cor=T)
summary(sed.pca)
plot(sed.pca)
biplot(sed.pca)
round(loadings(sed.pca)[,1:3],3)
# First three PCs are easily interpreted:
# PC1 = gradient coarse sand (high phi1 - phi3) and no mud/silt with 
#       low TOC/TON to high gravel content and high TOC/TON 
# PC2 = gradient from very coarse sand/gravel (highest pos. loading on phi.0) with high TOC / TON 
#       and C/N ratio to very fine sand (highest neg. loading on phi.4) substrate 
#       with lower organic content
# PC3 =  gradient from high fraction of coarse sand/gravel (phi.0 and phi.1) 
#       with low C/N ratio and high del N15 to little coarse material 
#       with high C/N ratio and low delN15

# Add first three PCs to other environmental data:
sed <- mutate(sed, sedPC1=sed.pca$scores[,"Comp.1"], sedPC2=sed.pca$scores[,"Comp.2"],
              sedPC3=sed.pca$scores[,"Comp.3"])
all.env <- left_join(all.env, select(sed, cruise, station.master, sedPC1:sedPC3))

rm(sed, sed.pca)


####################################################################################
##### Community analysis:

#-----------------------------------------------------------------------------------
####  Bacteria:

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(sfc.bacteria,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
sfc.bacteria <- sfc.bacteria[,j>3]

# Remove three stations with missing environmental data first for simplicity:
# (from protists as well!!)
j <- sfc.microbe.hdr$station.master == "CL2.5" | 
  (sfc.microbe.hdr$station.master == "ML6.8" & sfc.microbe.hdr$cruise == "AMBON2015") |
  (sfc.microbe.hdr$station.master == "ML5.7" & sfc.microbe.hdr$cruise == "AMBON2015")
sfc.bacteria <- sfc.bacteria[!j,]
sfc.protist <- sfc.protist[!j,]
sfc.microbe.hdr <- sfc.microbe.hdr[!j,]

# NMDS ordination, using Bray-Curtis for presence/absence data:
bact.mds <- metaMDS(sfc.bacteria, k=2, binary=T, trymax=100)
plot(bact.mds, display="site")

# Identify combination of variables most strongly correlated with community:
bact.env <- select(sfc.microbe.hdr, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, strat, melt.sfc, runoff.sfc)
summary(bact.env)

bact.bc <- vegdist(sfc.bacteria, binary=T)  # Bray-Curtis distances
bact.be <- bioenv(bact.bc, select(bact.env, SST:runoff.sfc), upto=8)
summary(bact.be)

# Select best variables:
(bact.best <- bact.be$names[bact.be$models[[bact.be$whichbest]]$best])

sub <- bact.env[,bact.best]
(bact.ef1 <- envfit(bact.mds, sub, choices=c(1,2)))
plot(bact.mds, display="site", cex=1.2)
plot(bact.ef1)
title("Bacterial community ordination (NMDS 1 vs 2)")

# Sort from highest to lowest correlation for model fitting:
bact.best <- rev(names(sort(bact.ef1$vectors$r)))

# Partition variance along environmental gradients:
form <- as.formula(paste("bact.bc", paste(c(bact.best, "cruise"), collapse=" + "), sep=" ~ "))
adonis2(form, data=bact.env, by="margin", permutations = 4999)
#adonis(form, data=env.dat, perm=4999)

# Remove phosphate:
bact.best <- c("SST", "SSS", "melt.sfc")
form <- as.formula(paste("bact.bc", paste(c(bact.best, "cruise"), collapse=" + "), sep=" ~ "))
adonis2(form, data=bact.env, by="margin", permutations = 4999)

# Re-select variabls for plotting
sub <- bact.env[,bact.best]
(bact.ef1 <- envfit(bact.mds, sub, choices=c(1,2)))


#-----------------------------------------------------------------------------------
####  Protists:

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(sfc.protist,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
sfc.protist <- sfc.protist[,j>3]

# Removed three stations with missing environmental data first for simplicity:
# (from bacteria and protists - see under bacteria above!!)

# NMDS ordination, using Jaccard for presence/absence data:
prot.mds <- metaMDS(sfc.protist, k=2, binary=T, trymax=100)
plot(prot.mds, display="site")

# Identify combination of variables most strongly correlated with community:
prot.env <- select(sfc.microbe.hdr, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, strat:runoff.sfc)

prot.bc <- vegdist(sfc.protist, binary=T)  # Bray-Curtis distances
prot.be <- bioenv(prot.bc, select(prot.env, SST:runoff.sfc), upto=6)
summary(prot.be)

# Select best variables:
(prot.best <- prot.be$names[prot.be$models[[prot.be$whichbest]]$best])

sub <- prot.env[,prot.best]
(prot.ef1 <- envfit(prot.mds, sub, choices=c(1,2)))
plot(prot.mds, display="site", cex=1.2)
plot(prot.ef1)
title("Protist community ordination (NMDS 1 vs 2)")

# Sort from highest to lowest correlation for model fitting:
prot.best <- rev(names(sort(prot.ef1$vectors$r)))

# Partition variance along environmental gradients:
form <- as.formula(paste("prot.bc", paste(c(prot.best, "cruise"), collapse=" + "), sep=" ~ "))
adonis2(form, data=prot.env, by="margin", permutations = 4999)
#adonis(form, data=env.dat, perm=4999)

# Remove nitrite/ate:
prot.best <- c("SST", "phosphate","silicate")
form <- as.formula(paste("prot.bc", paste(c(prot.best, "cruise"), collapse=" + "), sep=" ~ "))
adonis2(form, data=prot.env, by="margin", permutations = 4999)

# Re-select best variables for overlaying vectors on ordination plot:
sub <- prot.env[,prot.best]
(prot.ef1 <- envfit(prot.mds, sub, choices=c(1,2)))



#-----------------------------------------------------------------------------------
####  Infauna:

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(infauna.cpue,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
infauna.cpue <- infauna.cpue[,j>3]

# Select environmental variables and eliminate cases with missing values:
inf.env <- select(infauna.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, sed.chla, strat:sedPC3)
(j <- attributes(na.omit(inf.env))$na.action)

inf.env <- inf.env[-j,]
infauna.red <- infauna.cpue[-j,]

# Fourth-root transform, then standardize to maximum:
inf.tf.std <- apply(infauna.red^0.25, 2, function(x) x/max(x)) 

# NMDS ordination:
inf.mds <- metaMDS(inf.tf.std, k=2, trymax=100)
plot(inf.mds, display="site")

# Identify combination of variables most strongly correlated with community:
inf.env <- select(inf.env, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, sed.chla, strat:sedPC3)
inf.be <- bioenv(inf.tf.std, select(inf.env, SST:sedPC3), upto=6)
summary(inf.be)

# Select best variables:
(inf.best <- inf.be$names[inf.be$models[[inf.be$whichbest]]$best])

# Re-select variables with no missing data for the "best" subset
inf.env <- select(infauna.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise, all_of(inf.best))
(j <- attributes(na.omit(inf.env))$na.action)
# None missing

# NMDS with matching data:
inf.tf.std <- apply(infauna.cpue^0.25, 2, function(x) x/max(x)) 
inf.mds <- metaMDS(inf.tf.std, k=2, trymax=100, autotransform = F)

sub <- inf.env[,inf.best]
(inf.ef1 <- envfit(inf.mds, sub, choices=c(1,2)))
plot(inf.mds, display="site", cex=1.2)
plot(inf.ef1)
title("Infaunal community ordination (NMDS 1 vs 2)")

# Sort from highest to lowest correlation for model fitting:
inf.best <- rev(names(sort(inf.ef1$vectors$r)))

# Partition variance along environmental gradients:
inf.bc <- vegdist(inf.tf.std)  # Bray-Curtis distances
form <- as.formula(paste("inf.bc", paste(c(inf.best, "cruise"), collapse=" + "), sep=" ~ "))
adonis2(form, data=inf.env, by="margin", permutations = 4999)


#-----------------------------------------------------------------------------------
#### 2. Epifauna:

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(epi.cpue,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
epi.cpue <- epi.cpue[,j>3]

# Select environmental variables and eliminate cases with missing values:
epi.env <- select(epi.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, sed.chla, 
         strat:runoff.bot, sedPC1:sedPC3)
(j <- attributes(na.omit(epi.env))$na.action)

epi.env <- epi.env[-j,]
epi.red <- epi.cpue[-j,]

# Fourth-root transform, then standardize to maximum:
epi.tf.std <- apply(epi.red^0.25, 2, function(x) x/max(x)) 
hist(epi.tf.std)

# Identify combination of variables most strongly correlated with community:
epi.env <- select(epi.env, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, sed.chla, 
         strat:runoff.bot, sedPC1:sedPC3)
epi.be <- bioenv(epi.tf.std, select(epi.env, SST:sedPC3), upto=6)
summary(epi.be)

# Select best variables:
(epi.best <- epi.be$names[epi.be$models[[epi.be$whichbest]]$best])


# Re-select variables with no missing data for the "best" subset
epi.env <- select(epi.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise, all_of(epi.best))
(j <- attributes(na.omit(epi.env))$na.action)
# None missing

# NMDS with matching data:
# NMDS ordination:
epi.tf.std <- apply(epi.cpue^0.25, 2, function(x) x/max(x)) 
epi.mds <- metaMDS(epi.tf.std, k=2, trymax=100, autotransform = F)
plot(epi.mds, display="site")

sub <- epi.env[,epi.best]

(epi.ef <- envfit(epi.mds, sub))
plot(epi.mds, display="site", cex=1.2)
plot(epi.ef)
title("Epifaunal community ordination")

# Sort from highest to lowest correlation for model fitting:
(epi.best <- rev(names(sort(epi.ef$vectors$r))))

# Partition variance along environmental gradients:
epi.bc <- vegdist(epi.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("epi.bc", paste(c(epi.best, "cruise"), collapse=" + "), sep=" ~ ")))
adonis2(form, data=epi.env, by="margin", permutations = 4999)


#-----------------------------------------------------------------------------------
#### 3. Zooplankton (150 micron)

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(zoop150.cpue,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
zoop150.cpue <- zoop150.cpue[,j>3]

# Select environmental variables and eliminate cases with missing values:
# Modify interactively once key variables are determined
zoop150.env <- select(zoop150.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, strat:runoff.sfc)
(j <- attributes(na.omit(zoop150.env))$na.action)

zoop150.env <- zoop150.env[-j,]
zoop150.red <- zoop150.cpue[-j,]

# Fourth-root transform, then standardize to maximum:
zoop150.tf.std <- apply(zoop150.red^0.25, 2, function(x) x/max(x)) 
hist(zoop150.tf.std)


# Identify combination of variables most strongly correlated with community:
zoop150.env <- select(zoop150.env, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, strat:runoff.sfc)
zoop150.be <- bioenv(zoop150.tf.std, select(zoop150.env, SST:runoff.sfc), upto=6)
summary(zoop150.be)

# Select best variables:
(zoop150.best <- zoop150.be$names[zoop150.be$models[[zoop150.be$whichbest]]$best])
sub <- zoop150.env[,zoop150.best]

# Re-select variables with no missing data for the "best" subset
zoop150.env <- select(zoop150.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise, all_of(zoop150.best))
(j <- attributes(na.omit(zoop150.env))$na.action)

# NMDS with matching data:
zoop150.tf.std <- apply(zoop150.cpue^0.25, 2, function(x) x/max(x)) 
zoop150.mds <- metaMDS(zoop150.tf.std, k=3, trymax=150, autotransform = F)

(zoop150.ef <- envfit(zoop150.mds, select(zoop150.env,-cruise)))
plot(zoop150.mds, display="site", cex=1.2)
plot(zoop150.ef)
title("Zooplankton community ordination")

# Sort from highest to lowest correlation for model fitting:
(zoop150.best <- rev(names(sort(zoop150.ef$vectors$r))))

# Partition variance along environmental gradients:
zoop150.bc <- vegdist(zoop150.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("zoop150.bc", paste(c(zoop150.best, "cruise"), collapse=" + "), sep=" ~ ")))
adonis2(form, data=zoop150.env, by="margin", permutations = 4999)


#-----------------------------------------------------------------------------------
#### Zooplankton (505 micron)

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(zoop505.cpue,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
zoop505.cpue <- zoop505.cpue[,j>3]

# Select environmental variables and eliminate cases with missing values:
# Modify interactively once key variables are determined
zoop505.env <- select(zoop505.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, strat:runoff.sfc) %>% 
  mutate(Pac.sfc = 1-runoff.sfc-melt.sfc) %>% 
  select(-melt.bot, -runoff.bot, -melt.sfc, -runoff.sfc)
(j <- attributes(na.omit(zoop505.env))$na.action)

# Initial explorations suggested that surface runoff fraction may be important,
# if not significant at 95%. I replaced it with fraction of Pacific water,
# which did not emerge as an important variable, hence proceeded without

zoop505.env <- zoop505.env[-j,]
zoop505.red <- zoop505.cpue[-j,]

# Fourth-root transform, then standardize to maximum:
zoop505.tf.std <- apply(zoop505.red^0.25, 2, function(x) x/max(x)) 
hist(zoop505.tf.std)

# Identify combination of variables most strongly correlated with community:
zoop505.be <- bioenv(zoop505.tf.std, select(zoop505.env, SST:Pac.sfc), upto=8)
summary(zoop505.be)

# Select best variables:
(zoop505.best <- zoop505.be$names[zoop505.be$models[[zoop505.be$whichbest]]$best])

# Re-select variables with no missing data for the "best" subset
zoop505.env <- select(zoop505.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise, all_of(zoop505.best))
(j <- attributes(na.omit(zoop505.env))$na.action)

zoop505.env <- zoop505.env[-j,]
zoop505.red <- zoop505.cpue[-j,]

# NMDS with matching data:
zoop505.tf.std <- apply(zoop505.red^0.25, 2, function(x) x/max(x)) 
zoop505.mds <- metaMDS(zoop505.tf.std, k=2, trymax=100, autotransform = F)

(zoop505.ef <- envfit(zoop505.mds, zoop505.env[,-1]))
plot(zoop505.mds, display="site", cex=1.2)
plot(zoop505.ef)
title("Zooplankton community ordination")

# Sort from highest to lowest correlation for model fitting:
(zoop505.best <- rev(names(sort(zoop505.ef$vectors$r))))


# Partition variance along environmental gradients:
zoop505.bc <- vegdist(zoop505.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("zoop505.bc", paste(c(zoop505.best, "cruise"), collapse=" + "), sep=" ~ ")))
adonis2(form, data=zoop505.env, by="margin", permutations = 4999)


#-----------------------------------------------------------------------------------
#### Fish

# Eliminate rare species (species that occur at fewer than 3 stations):
j <- apply(fish.cpue,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
fish.cpue <- fish.cpue[,j>3]

# Select environmental variables and eliminate cases with missing values:
# Modify interactively once key variables are determined
fish.env <- select(fish.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, sed.chla, 
         strat:runoff.bot, sedPC1:sedPC3)
(j <- attributes(na.omit(fish.env))$na.action)

fish.env <- fish.env[-j,]
fish.red <- fish.cpue[-j,]

# Fourth-root transform, then standardize to maximum:
fish.tf.std <- apply(fish.red^0.25, 2, function(x) x/max(x)) 
hist(fish.tf.std)

# Identify combination of variables most strongly correlated with community:
fish.be <- bioenv(fish.tf.std, select(fish.env, SST:sedPC3), index=, upto = 6)
summary(fish.be)
# Marginal increase in correlation when adding more variables, hence
# for parsimony, select two most important variables

# Select best variables:
#(fish.best <- fish.be$names[fish.be$models[[fish.be$whichbest]]$best])
fish.best <- c("BT", "sedPC1")

# Re-select variables with no missing data for the "best" subset
fish.env <- select(fish.hauls, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise, all_of(fish.best))
(j <- attributes(na.omit(fish.env))$na.action)


# NMDS with matching (full) data:
fish.tf.std <- apply(fish.cpue^0.25, 2, function(x) x/max(x)) 
fish.mds <- metaMDS(fish.tf.std, k=3, trymax=100, autotransform = F)

(fish.ef1 <- envfit(fish.mds, select(fish.env, fish.best), choices=c(1,2)))
plot(fish.mds, display="site", cex=1.2, choices=c(1,2))
# plot(fish.mds, display="site", cex=1.2, choices=c(2,3))
plot(fish.ef1)
title("Fish community ordination (NMDS 1 vs 2)")

# Interannual differences are evident in NMDS 3 only:
boxplot(split(fish.mds$points[,3], fish.env$cruise))
# For final figure, could display sites 3 vs 1 if desired by  dropping 2nd dimension:
# fish.mds$points <- fish.mds$points[,-2]

# Partition variance along environmental gradients:
fish.bc <- vegdist(fish.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("fish.bc", paste(c(fish.best, "cruise"), collapse=" + "), sep=" ~ ")))
adonis2(form, data=fish.env, by="margin", permutations = 4999)


#-----------------------------------------------------------------------------------
#### Seabirds

# Eliminate rare species (species that occur at fewer than 3 stations):
j <- apply(bird.cpue,2, function(x) sum(x > 0))
sum(j<3)  # Number of species eliminated:
bird.cpue <- bird.cpue[,j>2]

# Add dummy to deal with station that had no birds:
bird.cpue <- cbind(bird.cpue,dummy=1)

# Select environmental variables and eliminate cases with missing values:
# Modify interactively once key variables are determined
bird.env <- select(bird.hdr, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, strat:sedPC1)
(j <- attributes(na.omit(bird.env))$na.action)

bird.env <- bird.env[-j,]
bird.red <- bird.cpue[-j,]

# Fourth-root transform, then standardize to maximum:
birds.tf.std <- apply(bird.red^0.25, 2, function(x) x/max(x)) 
hist(birds.tf.std[,-20])

# Identify combination of variables most strongly correlated with community:
bird.env <- select(bird.env, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise:depth, SST:silicate, int.chla, strat:sedPC1)
birds.be <- bioenv(birds.tf.std, select(bird.env, SST:sedPC1), upto=6)
summary(birds.be)

# Select best variables:
(birds.best <- birds.be$names[birds.be$models[[birds.be$whichbest]]$best])

# Re-select variables with no missing data for the "best" subset
bird.env <- select(bird.hdr, cruise, station.master) %>% 
  left_join(all.env) %>% 
  select(cruise, all_of(birds.best))
(j <- attributes(na.omit(bird.env))$na.action)

bird.env <- bird.env[-j,]
bird.red <- bird.cpue[-j,]

# NMDS with matching data:
birds.tf.std <- apply(bird.red^0.25, 2, function(x) x/max(x)) 
birds.mds <- metaMDS(birds.tf.std, k=2, trymax=100, autotransform = F)

sub <- bird.env[,birds.best]
(birds.ef1 <- envfit(birds.mds, sub, choices=c(1,2), na.rm=T))
plot(birds.mds, display="site", cex=1.2)
plot(birds.ef1)
title("Bird community ordination (NMDS 1 vs 2)")

# Sort from highest to lowest correlation for model fitting:
(birds.best <- rev(names(sort(birds.ef1$vectors$r))))

# Partition variance along environmental gradients:
birds.bc <- vegdist(birds.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("birds.bc", paste(c(birds.best, "cruise"), collapse=" + "), sep=" ~ ")))
adonis2(form, data=bird.env, by="margin", permutations = 4999)


# Save output for later and for plotting results:
save.image("results/CommunityAnalysis results.RData")

# See 'Ordinations Env plots.R' for graphical summary of results
