#####################################################################

# NMDS ordinations by assemblage with environmental vectors
# (Figure 6 for cross-community manuscript)

# Author: Franz Mueter 
# Last updated: 9/17/2021

#####################################################################

# Load previous results or, if needed, update analysis using 
# the script 'OrdinationAnalysis.R'
load("results/CommunityAnalysis results.RData")

library(plotrix)

# Plot all ordination results (biplots) in consistent format
plot.ord <- function(mds, ef, env, TITLE, vec.lab=F, LEN=1) {
  COL <- c("skyblue3", "orangered")[as.numeric(factor(env$cruise))]
  plot(mds$points[,1:2], cex=1.8, pch=16, col=COL, axes=F, asp=1); box()
  title(TITLE, line=0.5, cex.main=1.8)
  if(vec.lab) {
    plot(ef, col="black", cex=1.5) 
  } else {
    X <- 1.2 * LEN * (ef$vectors$arrows[,1] * sqrt(ef$vectors$r))
    Y <- 1.2 * LEN * (ef$vectors$arrows[,2] * sqrt(ef$vectors$r))
    arrows(0, 0, X, Y, lwd=1.6, length=0.12, col="black")
  }
}

png(file="./plots/Ordinations_Fig6_no_labels.png",width=8,height=4,units = "in",res=300)
#png(file="./plots/Ordinations_Fig6.png", width=600, height=300)

par(mfrow=c(2,4), bg="white", mar=c(0,0,2,0.5))

# Bacteria:
plot.ord(mds = bact.mds, ef = bact.ef1, env = bact.env, 
         TITLE = "Bacteria", LEN=1.1)

# Protists:
plot.ord(mds = prot.mds, ef = prot.ef1, env = prot.env, 
         TITLE = "Protists", LEN=1.2)

# Zoopankton (150 micron)
plot.ord(mds = zoop150.mds, ef = zoop150.ef, env = zoop150.env,
         TITLE = "Small Zoopl.", LEN=0.9)

# Zoopankton (505 micron)
# rotate for plotting to align with seabirds...
zoop505.mds.rot <- zoop505.mds
zoop505.mds.rot$points[,2] <- -zoop505.mds.rot$points[,2]
zoop505.ef.rot <- zoop505.ef
zoop505.ef.rot$vectors$arrows[,2] <- -zoop505.ef.rot$vectors$arrows[,2] 
plot.ord(mds = zoop505.mds.rot, ef = zoop505.ef.rot, env = zoop505.env,
         TITLE = "Large Zoopl.", LEN=1.2)

# Infauna:
plot.ord(mds = inf.mds, ef = inf.ef1, env = inf.env, 
         TITLE = "Macroinfauna", LEN=1.1)

# Epifauna
plot.ord(mds = epi.mds, ef = epi.ef, env = epi.env, 
         TITLE = "Epifauna", LEN=1.2)

# Demersal Fish
plot.ord(mds = fish.mds, ef = fish.ef1, env = fish.env, 
         TITLE = "Demersal fish", LEN=1.2)

# Seabirds
plot.ord(mds = birds.mds, ef = birds.ef1, env = bird.env,
         TITLE = "Seabirds", LEN=1.2)

dev.off()


