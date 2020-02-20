## CREATE BIAS LAYER TO SAMPLE BACKGROUND POINTS FOR THE SDMs 


require("raster")
require("dismo")
require('dismo')
require("rJava")

#in case of .Jinit error: #https://support.apple.com/kb/DL1572?viewlocale=en_US&locale=en_US
setwd(file.path("~", "Dropbox", "trade-biodiversity discovery", "discovery_paper_1"))

bias_preds <- c("distbuiltup", "distroads", "protecteda", "popdens", "roughness")

#get files
file_list <- list.files("./processed_data/layers_sdm/", full.names = T, pattern = "til")
layers <- stack(file_list[grepl(paste0(c(bias_preds), collapse = "|"), file_list)])

names(layers) <- gsub('.{15}$', '', names(layers))


birds <- read.csv(file.path("processed_data", "aves_processed_til.csv"))
coords <- birds[,c(1,2)]
subs <- sample(1:nrow(coords), size = 10000)
spdf <- SpatialPointsDataFrame(coords = coords[subs,], 
                               data = birds[subs,],proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
plot(spdf)

#Only need spp names & lat/long cols

mod <- dismo::maxent(x= layers, p= spdf, nbg = 30000, factors = "protecteda", args= c("randomtestpoints=25","-J", "-P", "-p", "-h", "threshold=FALSE"))

pred <- dismo::predict(mod, layers)
writeRaster(pred, paste0(file.path("processed_data", "layers_bias", "bias_til.tif"), format = "GTiff", overwrite = T))
