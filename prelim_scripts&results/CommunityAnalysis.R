library(mgcv)
library(vegan)
library(maps)
library(mapdata)

setwd("C:/Dropbox/Research/Projects/AMBON/Analysis/Community analyses")
source("Prep_CPUE_data.R")

# Select sediment parameters for analysis:
sed <- sediment[,c(9:15,17:22)]
summary(sed)  # 2 missing values in 13C, hence eliminate 
sed <- select(sed, -`d13C (per mil vs VPDB)`)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(sed, diag.panel=panel.hist)
# Transform coarse fractions, which tend to be highly right-skewed
sed[,1:4] <- sed[,1:4]^0.25
pairs(sed, diag.panel=panel.hist)

sed.pca <- princomp(sed, cor=T)
summary(sed.pca)
plot(sed.pca)
biplot(sed.pca)
loadings(sed.pca)
# First two PCs are easily interpreted:
# PC1 = gradient from fine silt/clay with high TOC / TON and d15N to more sandy substrate
#       with lower organic content
# PC2 = gradient from finer sand to high gravel content and high C/N ratio

# Combine first two PCs with other environmental variables in a single data frame:

env.dat <- select(PSBTA.hauls, c("Depth", "SST", "BT", "SSS", "Bsal")) %>% 
  mutate(sedPC1=sed.pca$scores[,"Comp.1"], sedPC2=sed.pca$scores[,"Comp.2"]) %>% 
  scale()
pairs(env.dat, diag.panel=panel.hist)
row.names(env.dat) <- PSBTA.hauls$locationID
env.dat <- as.data.frame(env.dat)

rm(sed)
####################################################################################
##### Community analysis:

#-----------------------------------------------------------------------------------
#### 1. Infauna:

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(infauna.cpue,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
infauna.cpue <- infauna.cpue[,j>3]

# Examine values to determine transformation/standardization:
hist(infauna.cpue)
hist(infauna.cpue^0.25)  # Still highly right skewed, dominated by zeros

# Square-root transform, then standardize to maximum:
inf.tf.std <- apply(infauna.cpue^0.25, 2, function(x) x/max(x)) 
hist(inf.tf.std)

# NMDS ordination:
inf.mds <- metaMDS(inf.tf.std, k=3, trymax=100)
plot(inf.mds, display="site")

# Identify combination of variables most strongly correlated with community:
inf.be <- bioenv(inf.tf.std, env.dat)
summary(inf.be)

# Select best variables:
(inf.best <- names(env.dat)[inf.be$models[[inf.be$whichbest]]$best])

sub <- env.dat[,inf.best]
(inf.ef1 <- envfit(inf.mds, sub, choices=c(1,2)))
plot(inf.mds, display="site", cex=1.2)
plot(inf.ef1)
title("Infaunal community ordination (NMDS 1 vs 2)")

(inf.ef2 <- envfit(inf.mds, sub, choices=c(2,3)))
plot(inf.mds, display="site", cex=1.2, choices=c(2,3))
plot(inf.ef2)
title("Infaunal community ordination (NMDS 2 vs 3)")

(inf.ef3 <- envfit(inf.mds, sub, choices=c(1,3)))
plot(inf.mds, display="site", cex=1.2, choices=c(1,3))
plot(inf.ef3)
title("Infaunal community ordination (NMDS 1 vs 3)")

# Fit in 3-D:
inf.ef.3D <- envfit(inf.mds, sub, choices=1:3)

# Sort from highest to lowest correlation for model fitting:
inf.best <- rev(names(sort(inf.ef.3D$vectors$r)))

# Partition variance along environmental gradients:
inf.bc <- vegdist(inf.tf.std)  # Bray-Curtis distances
form <- as.formula(paste("inf.bc", paste(inf.best, collapse=" + "), sep=" ~ "))
#adonis2(form, data=env.dat, by="margin", permutations = 4999)
#adonis(form, data=env.dat, perm=4999)

# 5 environmental gradients account for approximately 22.5% of the 
# variability in the community

#-----------------------------------------------------------------------------------
#### 2. Epifauna:

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(epifauna.cpue,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
epifauna.cpue <- epifauna.cpue[,j>3]

# Examine values to determine transformation/standardization:
hist(epifauna.cpue)
hist(epifauna.cpue^0.25)  # Still highly right skewed, dominated by zeros

# Square-root transform, then standardize to maximum:
epi.tf.std <- apply(epifauna.cpue^0.25, 2, function(x) x/max(x)) 
hist(epi.tf.std)

# NMDS ordination:
epi.mds <- metaMDS(epi.tf.std, k=2, trymax=100, autotransform = F)
plot(epi.mds, display="site")

# Identify combination of variables most strongly correlated with community:
epi.be <- bioenv(epi.tf.std, env.dat)
summary(epi.be)

# Select best variables:
(epi.best <- names(env.dat)[epi.be$models[[epi.be$whichbest]]$best])
sub <- env.dat[,epi.best]

(epi.ef <- envfit(epi.mds, sub))
plot(epi.mds, display="site", cex=1.2)
plot(epi.ef)
title("Epifaunal community ordination")

# Sort from highest to lowest correlation for model fitting:
(epi.best <- rev(names(sort(epi.ef$vectors$r))))

# Partition variance along environmental gradients:
epi.bc <- vegdist(epi.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("epi.bc", paste(epi.best, collapse=" + "), sep=" ~ ")))
#adonis2(form, data=env.dat, by="margin", permutations = 4999)
#adonis(form, data=env.dat, perm=4999)

# Four environmental variables account for approximately 29% of the
# variability in the epifaunal community
# Primary gradients are the depth/temperature gradient and sediment gradient

#-----------------------------------------------------------------------------------
#### 3. Zooplankton

# Eliminate rare species (species that occur at fewer than 4 stations):
j <- apply(zoop.cpue,2, function(x) sum(x > 0))
sum(j<4)  # Number of species eliminated:
zoop.cpue <- zoop.cpue[,j>3]

dim(zoop.cpue)

# Examine values to determine transformation/standardization:
hist(zoop.cpue)
hist(zoop.cpue^0.25)  # Still highly right skewed, dominated by zeros

# Square-root transform, then standardize to maximum:
zoop.tf.std <- apply(zoop.cpue^0.25, 2, function(x) x/max(x)) 
hist(zoop.tf.std)

# NMDS ordination:
zoop.mds <- metaMDS(zoop.tf.std, k=2, trymax=100, autotransform = F)
plot(zoop.mds, display="site")

# Identify combination of variables most strongly correlated with community:
zoop.be <- bioenv(zoop.tf.std, env.dat)
summary(zoop.be)

# Select best variables:
(zoop.best <- names(env.dat)[zoop.be$models[[zoop.be$whichbest]]$best])
sub <- env.dat[,zoop.best]

(zoop.ef <- envfit(zoop.mds, sub))
plot(zoop.mds, display="site", cex=1.2)
plot(zoop.ef)
title("Zooplankton community ordination")

# Sort from highest to lowest correlation for model fitting:
(zoop.best <- rev(names(sort(zoop.ef$vectors$r))))

# Partition variance along environmental gradients:
zoop.bc <- vegdist(zoop.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("zoop.bc", paste(zoop.best, collapse=" + "), sep=" ~ ")))
#adonis2(form, data=env.dat, by="margin", permutations = 4999)
#adonis(form, data=env.dat, perm=4999)

# Three environmental variables account for approximately 25.5% of the
# variability in the zooplankton community
# Primary gradients structuring the zooplankton community are bottom temperature 
# and bottom salinity (rather than surface layer characteristics!)

#-----------------------------------------------------------------------------------
#### 4. Fish
fish.cpue <- as.matrix(PSBTA.fish.cpue)

# Eliminate rare species (species that occur at fewer than 3 stations):
j <- apply(fish.cpue,2, function(x) sum(x > 0))
sum(j<3)  # Number of species eliminated:
fish.cpue <- fish.cpue[,j>2]

dim(fish.cpue)

# Examine values to determine transformation/standardization:
hist(fish.cpue)
hist(fish.cpue^0.25)  # Still highly right skewed, dominated by zeros

# Square-root transform, then standardize to maximum:
fish.tf.std <- apply(fish.cpue^0.25, 2, function(x) x/max(x)) 
hist(fish.tf.std)

# NMDS ordination:
fish.mds <- metaMDS(fish.tf.std, k=3, trymax=100, autotransform = F)
plot(fish.mds, display="site")

# Identify combination of variables most strongly correlated with community:
fish.be <- bioenv(fish.tf.std, env.dat)
summary(fish.be)

# Select best variables:
(fish.best <- names(env.dat)[fish.be$models[[fish.be$whichbest]]$best])
sub <- env.dat[,fish.best]

(fish.ef1 <- envfit(fish.mds, sub, choices=c(1,2)))
plot(fish.mds, display="site", cex=1.2)
plot(fish.ef1)
title("Fish community ordination (NMDS 1 vs 2)")

(fish.ef2 <- envfit(fish.mds, sub, choices=c(2,3)))
plot(fish.mds, display="site", cex=1.2, choices=c(2,3))
plot(fish.ef2)
title("Fish community ordination (NMDS 2 vs 3)")

(fish.ef3 <- envfit(fish.mds, sub, choices=c(1,3)))
plot(fish.mds, display="site", cex=1.2, choices=c(1,3))
plot(fish.ef3)
title("Fish community ordination (NMDS 1 vs 3)")

# Fit in 3-D:
fish.ef.3D <- envfit(fish.mds, sub, choices=1:3)


# Sort from highest to lowest correlation for model fitting:
(fish.best <- rev(names(sort(fish.ef.3D$vectors$r))))

# Partition variance along environmental gradients:
fish.bc <- vegdist(fish.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("fish.bc", paste(fish.best, collapse=" + "), sep=" ~ ")))
#adonis2(form, data=env.dat, by="margin", permutations = 4999)
#adonis(form, data=env.dat, perm=4999)

# Three environmental variables account for approximately 23% of the
# variability in the fish community
# Primary gradients structuring the fish community are sediment size 
# and salinity, with surface salinities having a stronger relationship
# than bottom salinities


#-----------------------------------------------------------------------------------
#### 4. Seabirds

# Eliminate rare species (species that occur at fewer than 3 stations):
j <- apply(birds.dens,2, function(x) sum(x > 0))
sum(j<3)  # Number of species eliminated:
birds.dens <- birds.dens[,j>2]

# Add dummy to deal with station that had no birds:
birds.dens <- cbind(birds.dens,1)

# Examine values to determine transformation/standardization:
hist(birds.dens)
hist(birds.dens^0.25)  # Still highly right skewed, dominated by zeros

# Square-root transform, then standardize to maximum:
birds.tf.std <- apply(birds.dens^0.25, 2, function(x) x/max(x)) 
hist(birds.tf.std)

# NMDS ordination:
birds.mds <- metaMDS(birds.tf.std, k=3, trymax=100, autotransform = F)
plot(birds.mds, display="site")

# Identify combination of variables most strongly correlated with community:
birds.be <- bioenv(birds.tf.std, env.dat)
summary(birds.be)

# Select best variables:
(birds.best <- names(env.dat)[birds.be$models[[birds.be$whichbest]]$best])
sub <- env.dat[,birds.best]

(birds.ef1 <- envfit(birds.mds, sub, choices=c(1,2)))
plot(birds.mds, display="site", cex=1.2)
plot(birds.ef1)
title("Bird community ordination (NMDS 1 vs 2)")

(birds.ef2 <- envfit(birds.mds, sub, choices=c(2,3)))
plot(birds.mds, display="site", cex=1.2, choices=c(2,3))
plot(birds.ef2)
title("Bird community ordination (NMDS 2 vs 3)")

(birds.ef3 <- envfit(birds.mds, sub, choices=c(1,3)))
plot(birds.mds, display="site", cex=1.2, choices=c(1,3))
plot(birds.ef3)
title("Bird community ordination (NMDS 1 vs 3)")

# Fit in 3-D:
birds.ef.3D <- envfit(birds.mds, sub, choices=1:3)


# Sort from highest to lowest correlation for model fitting:
(birds.best <- rev(names(sort(birds.ef.3D$vectors$r))))

# Partition variance along environmental gradients:
birds.bc <- vegdist(birds.tf.std)  # Bray-Curtis distances
# Model formula:
(form <- as.formula(paste("birds.bc", paste(birds.best, collapse=" + "), sep=" ~ ")))
adonis2(form, data=env.dat, by="margin", permutations = 4999)
adonis(form, data=env.dat, perm=4999)

# Three environmental variables account for approximately 30% of the
# variability in the bird community

# Primary gradients associated with the bird community are bottom temperature 
# and bottom salinity

rm(form, sub)

# save.image("CommunityAnalysis results.RData")

# Plot all ordination results (biplots) in a nice format
if(F) {
  library(plotrix)
  plot.ord <- function(mds, ef, best, TITLE, vec.lab=T) {
    col.fun <- colorRampPalette(c("lightblue", "darkblue"))
    # COL <- col.fun(65)[rank(env.dat[,best[1]])]
    COL <- col.fun(65)[rank(PSBTA.hauls$Latitude)]
    plot(mds$points[,1:2], cex=1.6, pch=16, col=COL, axes=F, asp=1); box()
    title(TITLE, line=0.5, cex.main=2)
    if(vec.lab) {
      plot(ef, col="red") 
      } else {
        X <- ef$vectors$arrows[,1] * sqrt(ef$vectors$r)
        Y <- ef$vectors$arrows[,2] * sqrt(ef$vectors$r)
        arrows(0, 0, X, Y, lwd=2, length=0.12, col=2)
      }
  }

  par(mfrow=c(2,3), bg="grey90", mar=c(0,0,2,0.5))

  # Infauna:
  plot.ord(mds = inf.mds, ef = inf.ef1, best = inf.best, 
       TITLE = "Infauna", vec.lab=T)

  # Epifauna
  plot.ord(mds = epi.mds, ef = epi.ef, best = epi.best, 
           TITLE = "Epifauna", vec.lab=T)

  # Demersal Fish
  plot.ord(mds = fish.mds, ef = fish.ef1, best = fish.best, 
           TITLE = "Demersal fish", vec.lab=T)

  # Zoopankton
  # rotate for plotting to align with seabirds...
  zoop.mds.rot <- zoop.mds
  zoop.mds.rot$points <- -zoop.mds.rot$points
  zoop.ef.rot <- zoop.ef
  zoop.ef.rot$vectors$arrows <- -zoop.ef$vectors$arrows 
  plot.ord(mds = zoop.mds.rot, ef = zoop.ef.rot, best = zoop.best, 
           TITLE = "Zooplankton", vec.lab=T)

  # Seabirds
  plot.ord(mds = birds.mds, ef = birds.ef1, best = birds.best, 
           TITLE = "Seabirds", vec.lab=T)



  # Add legend
  plot.new()
  # Lenged IF plotting latitudinal gradient:
  col.fun <- colorRampPalette(c("lightblue", "darkblue"))
  COL <- col.fun(65)
  color.legend(0.2, 0.85, 0.8, 0.95, c("South", "North"), COL, 
               align="rb", col = COL[c(15,65)])

}  # Close 'if' statement (to turn off/on plotting)

