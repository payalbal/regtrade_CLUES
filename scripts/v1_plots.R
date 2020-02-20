
## PAPER 1 PLOTS
## -----------------------------------------------------------------------------
## Code developed by Simon Kapitza
## Modified by Payal Bal
## Script generates plots (for bird species) for Vietnam under RCP 85 for three trade-scenarios: 
## 1. current (2014), 2. buisness-as-usual - bau (2070), 3. tpp (2070)
## Results generated and saved in dp1_sdms.R
## -----------------------------------------------------------------------------


rm(list = ls())
pacman::p_load(stringi, caret, nlme, plotrix, viridisLite, RColorBrewer, raster, scales, data.table, viridis)


## SET WORKING DIRECTORY FOR BOAB
setwd("~/discovery_paper_1")
cols <- viridis(3, begin = 0.1, end =0.9, direction = -1)
country_abbr <- "til"


## CHANGE IN AREA OCCUPIED BY SPECIES 
## -------------------------------------
## Caluclate % change in area occupiued by each species from 2014-1070 
##    as per SDM predictions under present, bau and tpp scenarios.
## CLARIFY: % change is in number of cells/area in km2 (because each cell is 1km2)

## Load SDM results
load(paste0("./output/sparea_", country_abbr, ".RData"))

## Calculate mean and se for change under each scenario (pr, bau, tpp)
spar.change <- (spar[,3:5]-spar[,3])/spar[,3] * 100
spar.change <- cbind(spar[,2],spar.change)
spar.mean <- numeric()
spar.se <- numeric()
for (i in 3:5){
spar.mean <- c(spar.mean, mean(spar.change[!is.infinite(spar.change[,i]),i], na.rm = T))
spar.se <- c(spar.se, std.error(spar.change[!is.infinite(spar.change[,i]),i], na.rm = T))
}

## Plot
page_width <- 17.83 #(15.7)

vnm_x <- c(1:3) + 0.125
par(mar =  c(4.1, 4.1, 0.1, 0.1))
plot(vnm_x, spar.mean, ylim = c(-40, 30), type = "n", xaxt = "n", ylab = "% change (2014 - 2070)", xlab = "Scenario", cex.lab = 1.5)
axis(1, at =vnm_x, labels = c("pr", "bau", "tpp"), cex.axis = 1.5)
arrows(vnm_x, spar.mean - 1 * spar.se ,vnm_x, spar.mean + 1.96 * spar.se, length=0.05, angle=90, code=3)
points(vnm_x, spar.mean, pch = 19, col =cols)
abline(h=0, lty = 2)

# png(file=paste0(path_plots, "mean_change.png"), units="cm", width=page_width * 1.4, height=page_width * 1.4 /1.5, res=600, pointsize = 18)
# legend("topright", legend = c("", "RCP 4.5", "RCP 8.5", "Australia", "Vietnam"), fill = c(cols, NA, NA), pch = c(NA,NA,NA, 17, 19), bty = "n", border = c(NA, NA, NA, NA, NA), cex = 1)
# dev.off()


## CHANGE IN BIOCLIM VARIABLES
## --------------------------------------
## Calculates % change in bioclim variables from 2014-2070 as per data downloaded 
##    from WorlClim website: http://www.worldclim.org/bioclim

layer_path <- "processed_data/layers_sdm"
files <- list.files(layer_path , full.names = T)
country_abbr <- "vnm"

## Load 2014 bioclim layers and calculate mean
layers_pres <- stack(files[grepl(paste0("(?=.*",country_abbr,")(?=.*bio)(?=.*0_0)"), files, perl = TRUE)])
layers_pres_df <- as.data.frame(layers_pres)
means <- colMeans(layers_pres_df, na.rm = T)
  
## Load 2070 bioclim layers and calculate mean
layers_fut <- stack(files[grepl(paste0("(?=.*",country_abbr,")(?=.*bio)(?=.*56)"), files, perl = TRUE)])
layers_fut_df <- as.data.frame(layers_fut)
means <- rbind(means, colMeans(layers_fut_df, na.rm = T))

## Calculate mean change
means <- t(means)
diff_means <- (means - means[,1]) /means[,1] * 100

## Plot
par(mar =  c(4.1, 4.1, 0.1, 0.1))
matplot(diff_means[,-1], ylim = c(floor(min(diff_means)), ceiling(max(diff_means))), pch = c(19), col = cols , xlab = "bioclim variables", ylab = "% change (2014 - 2070)", xaxt = "n")
abline(h = 0, lty = 2)
axis(side=1, at=1:19, labels=paste0("bio", c(1:19)),cex.axis=1, las = 2, cex.axis = 0.9)
grid()
  
  # legend("bottomleft", legend = c("RCP 2.6", "RCP 4.5", "RCP 8.5"), fill = cols, cex = 1)
  # png(file=paste0(path_plots, "/cc_", country_abbr, ".png"), units="cm", width=page_width/2, 
  #     height=page_width/2, res=600, pointsize = 10)
  # dev.off()


## CHANGE IN COMMODITY DEMANDS
## ----------------------------------------
## Calculates % change in commodity demands from 2014-2070 as per ...??
## CLARIFY: which layers are being loaded here: CGE outputs
##    provided by Tom and Ha ... demand

# rcps <- c("1oC", "2oC", "4oC")
# scens <- c("26", "45", "85")  # corresponding rcp forcings
# cn <- c("aus", "vn")

rcps <- c("85")   # previously as 'scens'
rcpsT <- c("4oC") # previously as 'rcps'
country <- "vn"   # previously as 'cn'
gtap_change <- list()


  P3 <- list()
  for (j in 1:length(rcps)){
    ## sector output 
    prod <- read.csv(paste0("./data/raw/prod and output/ao_", country, rcpsT[j], ".csv"))  
    ## sector productivity 
    dat <- read.csv(paste0("./data/raw/prod and output/qo_", country,"_", rcpsT[j], ".csv"), sep = ",")
    
    ## subset prod values to available GTAP classes
    prod[,2] <- trimws(as.character(prod[,2]))
    dat[,2] <- trimws(as.character(dat[,2]))
    s <- match(prod[,2], dat[,2])
    P1 <- dat[s[-c(which(prod[,2] == ""))], c(2,3, which(dat[1,] == 2071))]
    s <- match(P1[,1], prod[,2])
    P2 <- prod[s, c(2,3, which(dat[1,] == 2071))]
    
    #now we have the productivity (P2) and change (P1) in sector output, linearly interpolated to fill missing time steps - this is assumed to be proportional to land use change
    gtap_classes <- P1[,1]
    gtap_fullnames <- P1[,2]
    P1 <- P1[,-c(1,2)]
    P2 <- P2[,-c(1,2)]
    
    #calculate mean sector output, weighted by fao estimates of harvests
    P3[[j]] <- P1/100 + P2/100 + P1 * P2/10000
  }
  gtap_change <- cbind(gtap_classes, gtap_fullnames, as.data.frame(do.call("cbind", P3) * 100))

commodities <- c("pdr", "wht", "gro", "v_f", "osd", "c_b", "pfb", "ocr", "ctl", "frs")

classes <- which(gtap_classes%in%commodities)
comm_names <- trim(as.character(gtap_fullnames[classes]))
x1 <- c(1:length(classes) - 0.14)
x2 <- c(1:length(classes) + 0.14)

# png(file=paste0(plots, "/gtap.png"), units="cm", width=page_width, height=page_width - 2, res=600, pointsize = 18)

par(mar =  c(9, 4.1, 0.1, 0.1))
matplot(x1, gtap_change[classes, 3],  pch = 16, cex = 2, col = 'black', ylim = c(-70, 20), xlim = c(0.5, length(classes) + 0.5), ylab = "% Change (present - 2070)", xlab = "", xaxt='n')
axis(1, at = 1:length(classes), labels = F)
text(x=1:length(classes)-0.3, y=par()$usr[3]-0.09*(par()$usr[4]-par()$usr[3]),
     labels=comm_names, srt=45, adj=0.8, xpd=TRUE, cex = 1)
abline(h = 0, lty = 2)

# legend("bottomright", legend = c("RCP 2.6", "RCP 4.5", "RCP 8.5", "Australia", "Vietnam"), fill = c(cols, NA, NA), pch = c(NA,NA,NA, 17, 19), bty = "n", border = c(NA, NA, NA, NA, NA), cex = 0.9)
# dev.off()


## CHANGE IN LANDUSE
## ----------------------------
## Calculates change in area for each of the 5 landuse classes basesd on...
## CLARIFY...

lu_classes <- c("crop", "grass/shrub", "forest", "urban", "other")
library(viridisLite)
cols <- viridis(2, begin = 0.1, end =0.9, direction = -1)

## Load data
bau_landuse <- read.table("./data/processed/layers_clues/demand_vnm_bau")
names(bau_landuse) <- lu_classes
rownames(bau_landuse) <- c("2014","2030","2040", "2050", "2060", "2070")
tpp_landuse <- read.table("./data/processed/layers_clues/demand_vnm_tpp")
names(tpp_landuse) <- lu_classes
rownames(tpp_landuse) <- c("2014","2030","2040", "2050", "2060", "2070")

## Calculate percent change in lu class areas for both sceanrios
lu_change <- rbind((((bau_landuse[6,]-bau_landuse[1,])/bau_landuse[1,])*100), (((tpp_landuse[6,]-tpp_landuse[1,])/tpp_landuse[1,])*100))
rownames(lu_change) <- c("bau", "tpp")

## Plot: % change
par(mar =  c(5.1, 4, 0.1, 0.1))
matplot(1:5, xlim = c(1, 5), t(lu_change), pch = 17, col = c("grey40", "gold2"), ylab = "", xlab = "", xaxt = "n", cex=2)
grid()
# matpoints(1:5, xlim = c(1, 5), t(lu_change), pch = 4, col = c("grey40", "gold2"), cex=3)
axis(1, at = 1:5, labels = lu_classes, cex = 4)
title(ylab = "% change in area", 2, line = 2, cex.main = 4)
legend("topleft", legend = c("BAU", "TPP"), col = c("grey40", "gold2"), pch = c(17,17), bty = "n", border = c(NA, NA), cex = 2)
abline(h=0, lty = 2, col="grey50")
# text(4,6,"current", col="grey50", cex = 2)

# text(x=1:5, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
#      labels=lu_classes, srt=45, adj=1, xpd=TRUE, cex = 0.8)
# png(file=paste0(path_plots, "/land_use_change_vnm.png"), units="cm", width=page_width/1.9, height=page_width/1.9, res=600, pointsize = 18)
# legend("topright", legend = c("Present", "RCP 2.6", "RCP 4.5", "RCP 8.5"), fill = c(cols, NA), pch = c(NA,NA,NA, NA), bty = "n", border = c(NA, NA, NA, NA), cex = 0.8)
# dev.off()

