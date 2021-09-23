#####################################################################

# Compute pairwise Mantel correlations between eight assemblages sampled
# in the Chukchi Sea in 2015 and 2017 by year

# Author: Franz Mueter 
# Last updated: 9/8/2021

#####################################################################

# Load required packages
library(vegan)

# Run 'OrdinationAnalysis.R' or load stored results:
load("results/CommunityAnalysis results.RData")

# Mantel correlations among communities (after removal of 'rare' species):
mantel15.mat <- mantel17.mat <- matrix(NA, 8,8)
assemblages <- c("bacteria", "protists", "zoop150", "zoop505", "infauna", "epifauna", "fish", "birds")
dimnames(mantel15.mat) <- dimnames(mantel17.mat) <- list(assemblages,assemblages)
mantel17.mat

# Compile list of header info and CPUEs for analysis:

bind_cols(select(sfc.microbe.hdr,cruise,station.master), as.data.frame(sfc.bacteria))

st.cpue.all <- list(bacteria = bind_cols(select(sfc.microbe.hdr,cruise,station.master), as.data.frame(sfc.bacteria)),
                    protists = bind_cols(select(sfc.microbe.hdr,cruise,station.master), as.data.frame(sfc.protist)),
                    zoop150 = bind_cols(select(zoop150.hauls,cruise,station.master), as.data.frame(zoop150.cpue)),
                    zoop505 = bind_cols(select(zoop505.hauls,cruise,station.master), as.data.frame(zoop505.cpue)),
                    infauna = bind_cols(select(infauna.hauls,cruise,station.master), as.data.frame(infauna.cpue)),
                    epifauna = bind_cols(select(epi.hauls,cruise,station.master), as.data.frame(epi.cpue)),
                    fish = bind_cols(select(fish.hauls,cruise,station.master), as.data.frame(fish.cpue)),
                    birds = bind_cols(select(bird.hdr,cruise,station.master), as.data.frame(bird.cpue)))

for(i in 2:8) {
  for(j in 1:(i-1)) {
    n.i <- ncol(st.cpue.all[[i]]) - 2
    n.j <- ncol(st.cpue.all[[j]]) - 2
    intrsect <- left_join(st.cpue.all[[i]], st.cpue.all[[j]], by=c("cruise","station.master"))
    # 2015
    dat <- intrsect %>% filter(cruise=="AMBON2015") %>% select(-cruise,-station.master)
    dat <- na.omit(dat)
    x.i <- dat[,1:n.i]
    x.j <- dat[,(n.i+1):(n.i+n.j)]
    bc.i <- vegdist(x.i)
    bc.j <- vegdist(x.j)
    mnt <- mantel(bc.i, bc.j, method="spearman", perm=1999)
    mantel15.mat[i,j] <- mnt$stat
    mantel15.mat[j,i] <- mnt$signif
    # 2017
    dat <- intrsect %>% filter(cruise=="AMBON2017") %>% select(-cruise,-station.master)
    dat <- na.omit(dat)
    x.i <- dat[,1:n.i]
    x.j <- dat[,(n.i+1):(n.i+n.j)]
    bc.i <- vegdist(x.i)
    bc.j <- vegdist(x.j)
    mnt <- mantel(bc.i, bc.j, method="spearman", perm=1999)
    mantel17.mat[i,j] <- mnt$stat
    mantel17.mat[j,i] <- mnt$signif
  }
}


round(mantel15.mat,3)
round(mantel17.mat,3)



