# Wrapper script for the Etelis coruscans analyses

# ---- libraries ----
## global
library(tidyverse)      # the beautiful beautiful tidyverse
library(glue)           # use braces for PXA axis

## plot
library(ggplot2)        # plots
# library(ggh4x)          # nested facets
# library(patchwork)      # multiple plot by hand
# library(gridExtra)      # multiple plot from list
library(viridis)        # palette
# library(ggrepel)        # text labels on plots

## genetics
library(dartR)          # amazing wrapper for pop genetics analysis
library(pcadapt)
library(qvalue)         # used for pcadapt
# library(adegenet)
# library(hierfstat)
# library(mmod)
# library(LEA)            # SNMF



# ---- data ----
# data_samples <- read_csv("data/")



# ---- functions ----

## Melt distance matrix
melt.dist <- function(distmat, metric) {
  if(class(distmat)[1] == "dist") {distmat <- as.matrix(distmat)}
  distmat[upper.tri(distmat, diag = T)] <- NA
  distmat <- 
    as.data.frame(distmat) %>% 
    rownames_to_column(paste0(level, "1")) %>% 
    pivot_longer(cols = -paste0(level, "1"), 
                 names_to = paste0(level, "2"), 
                 values_to = metric) %>% 
    na.omit()
  
  return(distmat)
}


# species2[!(species2 %in% species1)]




# ---- arborescence ----
dir.create("data/", showWarnings = F)
dir.create("intermediate/", showWarnings = F)
dir.create("results/", showWarnings = F)




