## Set working environment 
install.packages("pacman")
library(pacman)
p_load("data.table","sp","raster","ochRe")

layer_path <- "./figures_rcp85" ## specify file path to the figures folder "./folder-path" 
  ## "./" reresents current folder. 
  ## Do not add a "/" at the end of the file path. 



## --------------------------------------------------------------------------------
## STACKED SDMs
## --------------------------------------------------------------------------------

## Load and mask rasters
mask_vnm <- raster(file.path(layer_path, "mask_vnm.tif"))
sum_000 <- mask(raster(file.path(layer_path, "sum_areaweighted_000.tif")), mask_vnm)  # current
sum_bau <- mask(raster(file.path(layer_path, "sum_areaweighted_bau.tif")), mask_vnm)  # bau sceanrio in 2070
sum_tpp <- mask(raster(file.path(layer_path, "sum_areaweighted_tpp.tif")), mask_vnm)  # tpp sceanrio in 2070


## Plot (note difference in range of values in llegend)
plot(sum_000, col = topo.colors(10))  # current
plot(sum_bau, col = topo.colors(10))  # bau sceanrio in 2070
plot(sum_tpp, col = topo.colors(10))  # tpp sceanrio in 2070

## Plot difference
plot(overlay(sum_000, sum_bau, fun=function(x,y) as.logical(round(x,3) == round(y,3))), 
     col= c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE) # current - bau
plot(overlay(sum_000, sum_tpp, fun=function(x,y) as.logical(round(x,3) == round(y,3))), 
     col=c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE)  # current - tpp



## --------------------------------------------------------------------------------
## SPECIES_SPECIFIC SDMs
## --------------------------------------------------------------------------------

## Load and plot rasters for 3 species
## "Aquila_nipalensis" - Endangered and decreasing acc to IUCN listing
## "Carpococcyx_renauldi" - Vulnerable and decreasing acc to IUCN listing
## "Terpsiphone_atrocaudata" - Near threatened acc to IUCN listing

load(file.path(layer_path, "selected.results_til.RData"))

sp.names.full <- c("Aquila_nipalensis", "Carpococcyx_renauldi", "Terpsiphone_atrocaudata")
sp.names.short <- c("aquila", "carpococcyx", "terpsiphone")
par(mfrow=c(1,3), cex = 0.8)
for (i in 1:3) {
  assign(paste0(sp.names.short[i], "_000"), selected.results[[i]][[4]][[1]])  # current
  assign(paste0(sp.names.short[i], "_bau"), selected.results[[i]][[4]][[2]])  # bau sceanrio in 2070
  assign(paste0(sp.names.short[i], "_tpp"), selected.results[[i]][[4]][[3]])  # tpp sceanrio in 2070
  
  plot(get(paste0(sp.names.short[i], "_000")), main = paste0(sp.names.full[i], "_000"),
       cex = 0.8, legend = FALSE)   # current
  plot(get(paste0(sp.names.short[i], "_bau")), main = paste0(sp.names.full[i], "_bau"),
       cex = 0.8, legend = FALSE)   # bau sceanrio in 2070
  plot(get(paste0(sp.names.short[i], "_tpp")), main = paste0(sp.names.full[i], "_tpp"),
       cex = 0.8)   # tpp sceanrio in 2070
}



## --------------------------------------------------------------------------------
## LAND USE
## --------------------------------------------------------------------------------

## Load rasters
lu_000 <- raster(file.path(layer_path, "landuse_vnm_luc_000_00.tif"))   # current
lu_bau <- raster(file.path(layer_path, "landuse_vnm_luc_bau_56.tif"))   # bau sceanrio in 2070
lu_tpp <- raster(file.path(layer_path, "landuse_vnm_luc_tpp_56.tif"))   # tpp sceanrio in 2070

## Load colour scheme
class.palette <- colorRampPalette(c(ochre_palettes[["lorikeet"]][3], ochre_palettes[["parliament"]][c(2,3)], "red", ochre_palettes[["parliament"]][6]))

## Plot
par(mfrow=c(1,1))
par(xpd=FALSE)
plot(lu_000, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE)  # current
plot(lu_bau, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE)  # bau sceanrio in 2070
plot(lu_tpp, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE)  # tpp sceanrio in 2070
plot(NULL)
legend("center", legend = c("crop", "grass/shrub", "forest", "urban", "bare"), fill = class.palette(5), bty = "n", cex = 3)

## Plot difference in land use
plot(overlay(lu_000, lu_bau,fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE) # current - bau
plot(overlay(lu_000, lu_tpp,fun=function(x,y) as.logical(x==y)), 
     col=c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE)  # current - tpp
plot(NULL)
legend("topright", legend = c("change", "no change"), fill=c("salmon", "steelblue"), bty = "n", cex = 3)

