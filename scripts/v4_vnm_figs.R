## Set working environment 

pacman::p_load("data.table","sp","raster","ochRe")
layer_path <- "./vnm_rcp85" ## specify file path to the figures folder "./folder-path" 
  ## "./" reresents current folder. 
  ## Do not add a "/" at the end of the file path. 



## --------------------------------------------------------------------------------
## STACKED SDMs
## --------------------------------------------------------------------------------

## Load and mask rasters
mask_vnm <- raster("./data/processed/masks/mask_vnm.tif")
sum_000 <- mask(raster(file.path(layer_path, "sum_areaweighted_000.tif")), mask_vnm)  # current
sum_bau <- mask(raster(file.path(layer_path, "sum_areaweighted_bau.tif")), mask_vnm)  # bau sceanrio in 2070
sum_tpp <- mask(raster(file.path(layer_path, "sum_areaweighted_tpp.tif")), mask_vnm)  # tpp sceanrio in 2070


## Plot (note difference in range of values in legend)
plot(sum_000, col = topo.colors(10), axes=F, box=F)  # current
plot(sum_bau, col = topo.colors(10), axes=F, box=F)  # bau sceanrio in 2070
plot(sum_tpp, col = topo.colors(10), axes=F, box=F)  # tpp sceanrio in 2070

  ## Richness_000 < Richness_bau or Richess_tpp!!??


s <- stack(sum_000, sum_bau, sum_tpp) 
sp <- as(s, 'SpatialGridDataFrame')
spplot(sp, names.attr=c('current', 'bau', 'tpp'))

par(mfrow = c(1,3))
sum_max <- max(c(values(sum_000), values(sum_bau),values(sum_tpp)), na.rm = TRUE)
sum_min <- min(c(values(sum_000), values(sum_bau),values(sum_tpp)), na.rm = TRUE)
plot(sum_000, legend = F, zlim=c(sum_min,sum_max), axes=F, box=F)  # current
plot(sum_bau, legend = F, zlim=c(sum_min,sum_max), axes=F, box=F)  # bau sceanrio in 2070
plot(sum_tpp, legend = F, zlim=c(sum_min,sum_max), axes=F, box=F)  # tpp sceanrio in 2070
dev.off()
plot(NULL)
plot(sum_000, legend.only = T, zlim=c(sum_min,sum_max), axes=F, box=F)
plot(NULL)
plot(sum_bau, legend.only = T, zlim=c(sum_min,sum_max), axes=F, box=F)

## Plot difference
plot(overlay(sum_000, sum_bau, fun=function(x,y) as.logical(round(x,3) == round(y,3))), 
     col= c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE) # current - bau
plot(overlay(sum_000, sum_tpp, fun=function(x,y) as.logical(round(x,3) == round(y,3))), 
     col=c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE)  # current - tpp



## --------------------------------------------------------------------------------
## SPECIES_SPECIFIC SDMs
## --------------------------------------------------------------------------------
## To identify grassland versus forest species see: https://figshare.com/articles/Data_Paper_Data_Paper/3559887

## Load and plot rasters for 3 species
## "Aquila_nipalensis" - Endangered and decreasing acc to IUCN listing
## "Carpococcyx_renauldi" - Vulnerable and decreasing acc to IUCN listing
## "Terpsiphone_atrocaudata" - Near threatened acc to IUCN listing

load("./output/selected.results_til.RData")

sp.names.full <- c("Aquila_nipalensis", "Carpococcyx_renauldi", "Terpsiphone_atrocaudata")
## 
sp.names.short <- c("aquila", "carpococcyx", "terpsiphone")

for (i in 1:3) {
  assign(paste0(sp.names.short[i], "_000"), selected.results[[i]][[4]][[1]])  # current
  assign(paste0(sp.names.short[i], "_bau"), selected.results[[i]][[4]][[2]])  # bau sceanrio in 2070
  assign(paste0(sp.names.short[i], "_tpp"), selected.results[[i]][[4]][[3]])  # tpp sceanrio in 2070
  
  par(mfrow=c(1,3), oma = c(0, 0, 2, 0))
  plot(get(paste0(sp.names.short[i], "_000")), legend = F, axes=F, box=F)   # current
  plot(get(paste0(sp.names.short[i], "_bau")), legend = F, axes=F, box=F)   # bau sceanrio in 2070
  plot(get(paste0(sp.names.short[i], "_tpp")), axes=F, box=F)   # tpp sceanrio in 2070
  mtext(sp.names.full[i], outer = TRUE, cex=1)
}

## Grassland species: Aquila_nipalensis
sum_max <- max(c(values(aquila_000), values(aquila_bau), 
                 values(aquila_tpp)), na.rm = TRUE)
sum_min <- min(c(values(aquila_000), values(aquila_bau), 
                 values(aquila_tpp)), na.rm = TRUE)
par(mfrow=c(1,1), oma = c(0, 0, 2, 0))
plot(aquila_000, zlim=c(sum_min,sum_max), legend = F, axes=F, box=F)   # current
plot(aquila_bau, zlim=c(sum_min,sum_max), legend = F, axes=F, box=F)   # bau sceanrio in 2070
plot(aquila_tpp, legend = F, axes=F, box=F)   # tpp sceanrio in 2070
dev.off()
plot(NULL)
plot(aquila_tpp, legend.only = T, zlim=c(sum_min,sum_max), axes=F, box=F)   # tpp sceanrio in 2070

## Zoom
e <- drawExtent()
plot(carpococcyx_tpp, ext = e, legend = F, axes=F, box=F)
plot(lu_bau, ext= e, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE, frame.plot=F)  # bau sceanrio in 2070

## Forest species: Carpococcyx_renauldi
sum_max <- max(c(values(carpococcyx_000), values(carpococcyx_bau), 
                 values(carpococcyx_tpp)), na.rm = TRUE)
sum_min <- min(c(values(carpococcyx_000), values(carpococcyx_bau), 
                 values(carpococcyx_tpp)), na.rm = TRUE)
par(mfrow=c(1,1), oma = c(0, 0, 2, 0))
plot(carpococcyx_000, legend = F, axes=F, box=F)   # current
plot(carpococcyx_bau, legend = F, axes=F, box=F)   # bau sceanrio in 2070
plot(carpococcyx_tpp, legend = F, axes=F, box=F)   # tpp sceanrio in 2070
dev.off()
plot(NULL)
plot(carpococcyx_bau, legend.only = T, axes=F, box=F)   # tpp sceanrio in 2070



## ---- XXXX ------
library(ochRe)
class.palette <- colorRampPalette(ochre_palettes[["winmar"]])
class.palette2 <- colorRampPalette(ochre_palettes[["dead_reef"]])
plot(aquila_000, col = class.palette(50), legend = T, axes=F, box=F)   # current
plot(carpococcyx_000, col = class.palette2(50), legend = T, axes=F, box=F)   # current
dev.off()
plot(NULL)
plot(aquila_000, legend.only = T, zlim=c(0,0.85), legend = T, axes=F, box=F)   # current




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
plot(lu_000, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE, frame.plot=F)  # current
plot(lu_bau, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE, frame.plot=F)  # bau sceanrio in 2070
plot(lu_tpp, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE, frame.plot=F)  # tpp sceanrio in 2070
plot(NULL)
legend("center", legend = c("crop", "grass/shrub", "forest", "urban", "bare"), fill = class.palette(5), bty = "n", cex = 3)

# rev.palette <- c("#304830", "#F0A800", "#F0A800", "#FF0000", "#F0D8D8")
rev.palette <- c("#F0A800", "#304830", "#F0A800", "#FF0000", "#90A8A8")
plot(lu_000, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE, frame.plot=F)
plot(lu_000, col = rev.palette, legend = NULL, box=FALSE, axes=FALSE, frame.plot=F)
plot(lu_tpp, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE, frame.plot=F)
plot(lu_tpp, col = rev.palette, legend = NULL, box=FALSE, axes=FALSE, frame.plot=F)


  ## Note: Forest_000 < Forest_tpp < Forest_bau..!!??


## Plot difference in land use
plot(overlay(lu_000, lu_bau,fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE) # current - bau
plot(overlay(lu_000, lu_tpp,fun=function(x,y) as.logical(x==y)), 
     col=c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE)  # current - tpp
plot(NULL)
legend("topright", legend = c("change", "no change"), fill=c("salmon", "steelblue"), bty = "n", cex = 3)

