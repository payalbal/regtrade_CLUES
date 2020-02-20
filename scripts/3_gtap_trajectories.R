# This script prepares demsnd trajectories as required for CLUES input. 

rm(list = ls())
setwd(file.path("~", "Dropbox", "trade-biodiversity discovery", "discovery_paper_1"))

#load function script
library("raster")
source(file.path("R", "dp1_functions.R"))

#Load FAO Harvested data
harvested <- read.csv(file.path("raw_data", "FAO_SEA_harvestedarea2015.csv"), header = T)

landuse <- raster(file.path("processed_data", "layers_clues", "cov_all.0.asc")) # i.e. "landuse_vnm_luc_000_00" == "cov_all.0.asc"
scens <- c("tpp", "bau")

#crop, grass/shrub, Forest, artificial, unused
clues_path <- file.path("processed_data", "layers_clues")
timesteps <- c(1, 16, 26, 36, 46, 56)
demand_2005 <- table(landuse[]) #5 classes: 0 crop, 1 grass/shrubland, 2 forest, 3 urban/artificial, 4 barren etc

for (i in 1:2){     #counter for scens i.e. two scenarios: tpp and bau
  #first row is ts 0
  yrs <- length()
  demand_2050 <- matrix(NA, nrow = yrs, ncol = length(demand_2005))
  
  #URBAN
  #Australien 2017 - 2050, then extrapolate to 2080
  urb_change <- (55739 - 33121.4)/33121.4
  an_urb_change <- urb_change/yrs
  demand_2050[,4] <- round(demand_2005[4] * (1 + c(1:(yrs)) * c(0, rep(an_urb_change, yrs-1)))) # assumes linear growth
  
  #CROPLAND
  #calculate crop areas and weight by their initial area occupation
  harv <- harvested
  
  #sum the area (in ha) produced in available GTAP classes (2005 estamtes)
  total_bysector <- aggregate(harv$Value, by = list(harv$GTAP.sector), FUN = function(x) sum(na.omit(x)))
  
  #relative area harvested of agricultural commodities as share of total area harvested in one of the GTAP sectors.
  total_bysector[,2] <- total_bysector[,2]/sum(total_bysector[,2])
  
  dat <- read.csv(file.path("raw_data", paste0("qo_", scens[i], "_2018.csv")))
  
  #subset prod values to available GTAP classes
  dat[,1] <- trimws(as.character(dat[,1]))
  s <- match(total_bysector[,1], as.character(dat[,1]))
  gtap_classes <- dat[,1]
  P1 <- dat
  
  P1 <- P1[,-1]
  
  P1 <- t(apply(P1, MARGIN = 1, FUN = function(x) approx_fun2(x,colnames(P1))))
  P1 <- P1[, 1:length(2015:2070)]
  
  #calculate mean sector output, weighted by fao estimates of harvests
  demand_2050[,1] <- round(demand_2005[1] * colMeans(total_bysector[,2] * (1+P1[s,]/100) * nrow(total_bysector)))
  
  #PASTURE
  s <- which(gtap_classes == "ctl")
  demand_2050[,2] <- round(demand_2005[2] * (1+P1[s,]/100))
  
  #FOREST
  #change proportional to output and prod change in forestry sector, also assuming that changes in crop and pasture will be at the expense of prim hab or lead to the establishment of prim hab
  
  s <- which(gtap_classes == "frs")
  demand_2050[,3] <- round(demand_2005[3] * (1 - P1[s,]/100))
  
  #GRSDD/SHRUBLAND
  #the "sink" for everything else
  
  demand_2050[,5] <- round(demand_2005[5])
  demand_2050[,2] <- round(demand_2050[,2] + (sum(demand_2005, na.rm = T) - rowSums(demand_2050, na.rm = T)))
  
  #i.e: land use in all land use classes is driven, apart from secondary habitat which is a result of land abandonment in other classes (i.e. succession, reduction in primary forest), so the net difference between 1 and the driven land use classes is the projected change in secondary habitat.
  demand <- demand_2050[timesteps,]
  demand <- demand * 100
  colnames(demand) <- NULL
  print(rowSums(demand))
  write.table(demand, file.path("processed_data", "layers_clues", paste0("demand_vnm_", scens[i])), sep = " ", row.names = F, col.names = F)
}
