#####################################################################

# Summarize environmental variables measured during two AMBON cruises
# in the Chukchi Sea by year:

# Author: Franz Mueter 
# Last updated: 5/23/2021

#####################################################################

library(tidyverse)

load("./data_Mar2021/StationData/all_env.RData") # # from 'AMBON env sed chla.R'

stns <- all.env %>% 
  filter(cruise == "AMBON2015", !is.na(phi.1)) %>% 
  pull(station.master)
  
common.env <- all.env %>% group_by(cruise) %>%
  filter(station.master %in% stns) %>% 
  select(latitude.master:sed.chla, phi.0:phi.4, phi.5, TOC:runoff.sfc)

f <- function(x) c(mean(x,na.rm=T), sd(x, na.rm=T), min(x,na.rm=T), max(x,na.rm=T), sum(!is.na(x)))
env.summ <- common.env %>% 
  summarize(across(everything(), f))

pv <- rep(NA, ncol(common.env))
library(mgcv)
for(i in 4:ncol(common.env)) {
  y <- names(common.env)[i]
  .form <- as.formula(paste(y, " ~ ", "cruise + s(latitude.master, longitude.master)"))
  fit <- gam(.form, data=common.env)
  pv[i] <- summary(fit)$p.pv[2]
}

out <- t(env.summ[,-1])
dimnames(out)[[2]] <- pull(env.summ,1)
out <- cbind(out, p.val.means = pv[-1])

write.table(out, "results/Env_vars_summary.csv", sep=",")
