#####################################################################

# Compute three metrics of species diversity (Hill numbers) for eight
# assemblages in the eastern Chukchi Sea between 2015 and 2017
# and plot species accumulation curves for each metric by assemblage
# and sampling year.

# Author: Franz Mueter 
# Last updated: 9/8/2021

#####################################################################

# Load required packages
require(tidyverse)
library(iNEXT)
library(grid)
library(gridExtra)

# Read in all data:
source("./scripts/AMBON_Import_all_data.R")

# For sensitivity analyses, run the following script to eliminate 
# stations along the BB line to assess the impact on results:
# source(./Scripts/Sensitivity_no_BBL.R)


##############################################################
## Bacterial diversity, surface waters
##############################################################

sfc.bact <- list('2015' = t(sfc.bacteria[sfc.microbe.hdr$cruise == "AMBON2015",]),
                '2017' = t(sfc.bacteria[sfc.microbe.hdr$cruise == "AMBON2017",]))

sfc.bact.div <- iNEXT(sfc.bact, q=c(0,1,2), datatype="incidence_raw", endpoint=82)

theme_get()$plot.margin
my.theme <-  function() theme(legend.position = "none", 
                   text=element_text(size=18),
                   axis.title.x = element_blank(),
                   strip.text = element_blank(),
                   plot.margin = unit(c(5.5, 5.5, 5.5, 10), "pt"))
  
  
p.sfc.bact <- ggiNEXT(sfc.bact.div, type=1, se=TRUE, facet.var = "order",
                     color.var="site", grey=FALSE)  + ylab("Sfc. bacteria") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))



##############################################################
## Bacterial diversity, mid water (chl max?) - NOT USED!
##############################################################

if(F) {
mid.bact <- list('2015' = t(mid.bacteria[mid.microbe.hdr$cruise == "AMBON2015",]),
                 '2017' = t(mid.bacteria[mid.microbe.hdr$cruise == "AMBON2017",]))
mid.bact[['2015']][mid.bact[['2015']]>0] <- 1
mid.bact[['2017']][mid.bact[['2017']]>0] <- 1

mid.bact.div <- iNEXT(mid.bact, q=c(0,1,2), datatype="incidence_raw", endpoint=82)

p.mid.bact <- ggiNEXT(mid.bact.div, type=1, se=TRUE, facet.var = "order",
                      color.var="site", grey=FALSE)  + ylab("Midwater bacteria") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))

}

##############################################################
## Protist diversity, surface waters
##############################################################

sfc.prot <- list('2015' = t(sfc.protist[sfc.microbe.hdr$cruise == "AMBON2015",]),
                 '2017' = t(sfc.protist[sfc.microbe.hdr$cruise == "AMBON2017",]))
sfc.prot[['2015']][sfc.prot[['2015']]>0] <- 1
sfc.prot[['2017']][sfc.prot[['2017']]>0] <- 1

sfc.prot.div <- iNEXT(sfc.prot, q=c(0,1,2), datatype="incidence_raw", endpoint=82)

p.sfc.prot <- ggiNEXT(sfc.prot.div, type=1, se=TRUE, facet.var = "order",
                      color.var="site", grey=FALSE)  + ylab("Protists") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))

##############################################################
## Metazoan diversity (from eDNA), surface waters
##############################################################

sfc.meta <- list('2015' = t(sfc.metazoa[sfc.microbe.hdr$cruise == "AMBON2015",]),
                 '2017' = t(sfc.metazoa[sfc.microbe.hdr$cruise == "AMBON2017",]))
sfc.meta[['2015']][sfc.meta[['2015']]>0] <- 1
sfc.meta[['2017']][sfc.meta[['2017']]>0] <- 1

sfc.meta.div <- iNEXT(sfc.meta, q=c(0,1,2), datatype="incidence_raw", endpoint=82)

p.sfc.meta <- ggiNEXT(sfc.meta.div, type=1, se=TRUE, facet.var = "order",
                      color.var="site", grey=FALSE)  + ylab("Metazoa") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))

##############################################################
## Zooplankton (150 micron) diversity
##############################################################

zoop150 <- list('2015' = t(zoop150.cpue[zoop150.hauls$cruise == "AMBON2015",]),
            '2017' = t(zoop150.cpue[zoop150.hauls$cruise == "AMBON2017",]))
zoop150[['2015']][zoop150[['2015']]>0] <- 1
zoop150[['2017']][zoop150[['2017']]>0] <- 1

zoop150.div <- iNEXT(zoop150, q=c(0,1,2), datatype="incidence_raw", endpoint=82)

p.zoop150 <- ggiNEXT(zoop150.div, type=1, se=TRUE, facet.var = "order",
                 color.var="site", grey=FALSE)  + ylab("Sm. Zoopl.") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))


##############################################################
## Zooplankton (505 micron) diversity
##############################################################

zoop505 <- list('2015' = t(zoop505.cpue[zoop505.hauls$cruise == "AMBON2015",]),
                '2017' = t(zoop505.cpue[zoop505.hauls$cruise == "AMBON2017",]))
zoop505[['2015']][zoop505[['2015']]>0] <- 1
zoop505[['2017']][zoop505[['2017']]>0] <- 1

zoop505.div <- iNEXT(zoop505, q=c(0,1,2), datatype="incidence_raw", endpoint=82)

p.zoop505 <- ggiNEXT(zoop505.div, type=1, se=TRUE, facet.var = "order",
                     color.var="site", grey=FALSE)  +  ylab("Lg. Zoopl.") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))

##############################################################
## Macroinfaunal diversity
##############################################################

inf <- list('2015' = t(infauna.cpue[infauna.hauls$cruise == "AMBON2015",]),
            '2017' = t(infauna.cpue[infauna.hauls$cruise == "AMBON2017",]))
inf[['2015']][inf[['2015']]>0] <- 1
inf[['2017']][inf[['2017']]>0] <- 1

inf.div <- iNEXT(inf, q=c(0,1,2), datatype="incidence_raw", endpoint=82)

p.inf <- ggiNEXT(inf.div, type=1, se=TRUE, facet.var = "order",
                 color.var="site", grey=FALSE)  + ylab("Macroinfauna") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))


##############################################################
## Epifaunal diversity
##############################################################

epi <- list('2015' = t(epi.cpue[epi.hauls$cruise == "AMBON2015",]),
             '2017' = t(epi.cpue[epi.hauls$cruise == "AMBON2017",]))
epi[['2015']][epi[['2015']]>0] <- 1
epi[['2017']][epi[['2017']]>0] <- 1

epi.div <- iNEXT(epi, q=c(0,1,2), datatype="incidence_raw", endpoint=82)

p.epi <- ggiNEXT(epi.div, type=1, se=TRUE, facet.var = "order",
                  color.var="site", grey=FALSE)  + ylab("Epibenthos") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))


## Fish diversity

##############################################################
## Fish diversity
##############################################################

fish <- list('2015' = t(fish.cpue[fish.hauls$cruise == "AMBON2015",]),
     '2017' = t(fish.cpue[fish.hauls$cruise == "AMBON2017",]))
fish[['2015']][fish[['2015']]>0] <- 1
fish[['2017']][fish[['2017']]>0] <- 1

fish.div <- iNEXT(fish, q=c(0,1,2), datatype="incidence_raw", endpoint=80)

p.fish <- ggiNEXT(fish.div, type=1, se=TRUE, facet.var = "order",
        color.var="site", grey=FALSE)  + ylab("Demersal fish") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))
  
##############################################################
## Seabird diversity
##############################################################

bird <- list('2015' = t(bird.cpue[bird.hdr$cruise == "AMBON2015",]),
             '2017' = t(bird.cpue[bird.hdr$cruise == "AMBON2017",]))
bird[['2015']][bird[['2015']]>0] <- 1
bird[['2017']][bird[['2017']]>0] <- 1

bird.div <- iNEXT(bird, q=c(0,1,2), datatype="incidence_raw", endpoint=80)

p.bird <- ggiNEXT(bird.div, type=1, se=TRUE, facet.var = "order",
                  color.var="site", grey=FALSE)  + ylab("Seabirds") +
  theme_bw() +
  my.theme() +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red"))



##############################################################
### Combine plots:
# (5 inches wide, 2 inches tall per community?)
png(file="./plots/Fig3_Diversity.png", width=9.5, height=8, units = "in", res=300)
# png(file="./plots/Fig3_Diversity_noBBL.png", width=950, height=800)
grid.arrange(p.sfc.bact, p.sfc.prot, p.zoop150, p.zoop505, p.inf, p.epi, p.fish, p.bird, 
             ncol=2, bottom = textGrob("Number of stations sampled",
                                       gp=gpar(fontsize=18)))
grid.text("Richness", 0.13, 0.97, gp=gpar(fontsize=12))
grid.text("Simpson", 0.27, 0.97, gp=gpar(fontsize=12))
grid.text("Shannon", 0.41, 0.97, gp=gpar(fontsize=12))
grid.text("Richness", 0.63, 0.97, gp=gpar(fontsize=12))
grid.text("Simpson", 0.77, 0.97, gp=gpar(fontsize=12))
grid.text("Shannon", 0.91, 0.97, gp=gpar(fontsize=12))
dev.off()
