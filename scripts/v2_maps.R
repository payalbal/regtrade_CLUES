library(sp)
library(raster)
library(mime)
library(base64enc)
library(mapview)
library(htmltools)
library(mapedit)
library(lattice)
library(RColorBrewer)
library(latticeExtra)
library(rasterVis)
library(viridisLite)


## SPECIES DISCTRIBUTIONS
## ----------------------------

## Load summed rasters and mask to remove NAs
mask_vnm <- raster(file.path("processed_data", "masks", paste0("mask_vnm.tif")))

sum1_000 <- mask(raster("./output/sum_000.tif"), mask_vnm)
sum1_bau <- mask(raster("./output/sum_bau.tif"), mask_vnm)
sum1_tpp <- mask(raster("./output/sum_tpp.tif"), mask_vnm)

sum2_000 <- mask(raster("./output/sum_areaweighted_000.tif"), mask_vnm)
sum2_bau <- mask(raster("./output/sum_areaweighted_bau.tif"), mask_vnm)
sum2_tpp <- mask(raster("./output/sum_areaweighted_tpp.tif"), mask_vnm)

# 1. Summed relative likelihood maps with consistent legend
r.range <- c(min(minValue(stack(sum1_000, sum1_bau, sum1_tpp))),max(maxValue(stack(sum1_000, sum1_bau, sum1_tpp))))
par(mfrow=c(1,3))
plot(sum1_000, col = topo.colors(60), breaks=seq(r.range[1],r.range[2], by=5), legend = FALSE, axes=F, box=F)
plot(sum1_bau, col = topo.colors(60), breaks=seq(r.range[1],r.range[2], by=5), legend = FALSE, axes=F, box=F)
plot(sum1_tpp, col = topo.colors(60), breaks=seq(r.range[1],r.range[2], by=5), legend = FALSE, axes=F, box=F)
## add legend
plot(sum1_000, legend.only=TRUE, col=topo.colors(60), legend.width=1.5, legend.shrink=0.5, axis.args=list(at=seq(floor(r.range[1]), ceiling(r.range[2]), ceiling(abs(r.range[1]-r.range[2])/4)), labels=seq(floor(r.range[1]), ceiling(r.range[2]), ceiling(abs(r.range[1]-r.range[2])/4)), cex.axis=1), legend.args=list(text='Richness\n indicator', side=3, font=2, line=1, cex=1.2))

## Find differences between scenarios
plot(overlay(sum1_000, sum1_tpp,fun=function(x,y) as.logical(x==y), col=class.palette(2)))



## 2. Area corrected richness maps - sum of relative contribution (importance) of a pixel to a species' distribution
r.range <- c(min(minValue(stack(sum2_000, sum2_bau, sum2_tpp))),max(maxValue(stack(sum2_000, sum2_bau, sum2_tpp))))
plot(sum2_000, col = viridis(10), breaks=seq(minValue(sum2_000), maxValue(sum2_000), length.out = 10), axes=F, box=F)
plot(sum2_bau, col = topo.colors(10))
plot(sum2_tpp, col = topo.colors(10))
mapview(sum2_000, legend=F)

## scale to higher resolution
temp <- aggregate(sum2_000, fact=6, fun=mean) # check res(temp): 6 gives 0.05, 0.05 resolution
plot(temp, col = viridis(10), axes=F, box=F)

levelplot(temp^2, contour=TRUE)


library(colorspace)
myTheme <- rasterTheme(region=sequential_hcl(10, power=2.2))
levelplot(temp, contour = TRUE)


## Find differences between scenarios
plot(overlay(sum2_bau, sum2_tpp,fun=function(x,y) as.logical(x==y), col=topo.colors(2)))
  ## This doesn't tell us much because values may be different in the 4-5 decimal place..
  ## TO DO: round(values(x),2)==round(values(y),2)

## LANDUSE CHANGE
layer_path <- "./data/processed/layers_sdm/"
library(ochRe)
class.palette <- colorRampPalette(c(ochre_palettes[["lorikeet"]][3], ochre_palettes[["parliament"]][c(2,3)], "red", ochre_palettes[["parliament"]][6]))
par(mfrow=c(1,1))

lu_000 <- raster(paste0(layer_path, "landuse_vnm_luc_000_00.tif"))
lu_bau <- raster(paste0(layer_path, "landuse_vnm_luc_bau_56.tif"))
lu_tpp <- raster(paste0(layer_path, "landuse_vnm_luc_tpp_56.tif"))

par(xpd=FALSE)
plot(lu_000, col = class.palette(5), legend = NULL, box=FALSE, axes=FALSE)
par(xpd=TRUE)
legend(101,17, legend = c("crop", "grass/shrub", "forest", "urban", "bare"), fill = class.palette(5), bty = "n", cex = 1.3)

plot(lu_bau, col = class.palette(5), legend = NULL)
plot(lu_tpp, col = class.palette(5), legend = NULL)

mapview(lu_000)
mapview(lu_bau)

## Difference between scenarios
plot(overlay(lu_000, lu_bau,fun=function(x,y) as.logical(x==y)), col= c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE)
plot(NULL)
legend("topright", legend = c("change", "no change"), fill=c("salmon", "steelblue"), bty = "n", cex = 3)
plot(overlay(lu_000, lu_tpp,fun=function(x,y) as.logical(x==y)), col=c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE)
plot(overlay(lu_bau, lu_tpp,fun=function(x,y) as.logical(x==y)), col=c("salmon", "steelblue"), legend = NULL, box=FALSE, axes=FALSE)
  ## overlay can also be used to find difference netween values, but values are categorical here so not appropriate


mapview(stack(sum1_000, sum1_bau,sum1_tpp), maxpixels =  1563741)


