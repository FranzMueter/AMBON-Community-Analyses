#######################################################################################
# Compute correlations among communities sampled (currently includes infauna,
# epifauna, demersal fish, zooplankton, and seabirds - NOT including microbes)

#### NOTE: This script file is superseded by 'Cross-community indicator species.Rmd'######


library(vegan)
library(mgcv)

# If output from 'CommunityAnalysis.R' was saved, open up the workspace, otherwise, 
# source 'CommunityAnalysis.R':

load("CommunityAnalysis results.RData")
# OR:
# source("CommunityAnalysis.R")


# Mantel correlations among communities (after removal of 'rare' species):
mantel.mat <- matrix(NA, 5,5)
dimnames(mantel.mat) <- list(c("infauna", "epifauna", "fish", "zoop", "birds"), 
                             c("infauna", "epifauna", "fish", "zoop", "birds"))
mantel.mat

mnt <- mantel(birds.bc, zoop.bc, method="spearman", perm=1999)
mantel.mat["birds","zoop"] <- mnt$stat
mantel.mat["zoop","birds"] <- mnt$signif

mnt <- mantel(birds.bc, fish.bc, method="spearman", perm=1999)
mantel.mat["birds","fish"] <- mnt$stat
mantel.mat["fish","birds"] <- mnt$signif

mnt <- mantel(birds.bc, epi.bc, method="spearman", perm=1999)
mantel.mat["birds","epifauna"] <- mnt$stat
mantel.mat["epifauna","birds"] <- mnt$signif

mnt <- mantel(birds.bc, inf.bc, method="spearman", perm=1999)
mantel.mat["birds","infauna"] <- mnt$stat
mantel.mat["infauna","birds"] <- mnt$signif

mnt <- mantel(zoop.bc, fish.bc, method="spearman", perm=1999)
mantel.mat["zoop","fish"] <- mnt$stat
mantel.mat["fish","zoop"] <- mnt$signif

mnt <- mantel(fish.bc, epi.bc, method="spearman", perm=1999)
mantel.mat["fish","epifauna"] <- mnt$stat
mantel.mat["epifauna","fish"] <- mnt$signif

mnt <- mantel(fish.bc, inf.bc, method="spearman", perm=1999)
mantel.mat["fish","infauna"] <- mnt$stat
mantel.mat["infauna","fish"] <- mnt$signif

mnt <- mantel(zoop.bc, epi.bc, method="spearman", perm=1999)
mantel.mat["zoop","epifauna"] <- mnt$stat
mantel.mat["epifauna","zoop"] <- mnt$signif

mnt <- mantel(epi.bc, inf.bc, method="spearman", perm=1999)
mantel.mat["epifauna", "infauna"] <- mnt$stat
mantel.mat["infauna", "epifauna"] <- mnt$signif

mnt <- mantel(zoop.bc, inf.bc, method="spearman", perm=1999)
mantel.mat["zoop","infauna"] <- mnt$stat
mantel.mat["infauna","zoop"] <- mnt$signif

mantel.mat

