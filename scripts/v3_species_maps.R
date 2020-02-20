library("data.table")
library(sp)
library(raster)

load("./output/sparea_til.RData")

## Calculate (proportional) change in area to find species with the most notable differences from pr scenario
spar$change_bau <- (spar[,pr]-spar[,bau])/spar[,pr]
spar$change_tpp <- (spar[,pr]-spar[,tpp])/spar[,pr]

## for TPP
## Number of species showing > 0.9 proprotional change in habitat area compared to present scenario
sum(spar[,(change_tpp > 0.9 & change_tpp < 1)], na.rm=T)
## Index of species with > 0.9 prop change
spar[which(spar[,(change_tpp > 0.9 & change_tpp < 1)])]$spID
## Extract data for species with > 0.9 prop change
tpp.sp <- spar[which(spar[,(change_tpp > 0.9 & change_tpp < 1)]),]

## look up species names at https://www.iucnredlist.org/ to find IUCN listed species
tpp.endangered <- c("Aquila nipalensis", "Leptoptilos dubius")
tpp.vulnerable <- c("Carpococcyx renauldi", "Sitta formosa")
tpp.nearthreatened <- c("Terpsiphone atrocaudata", "Vanellus vanellus")
sp.ind <- c(spar[which(spar[,(species %in% tpp.endangered)])]$spID, spar[which(spar[,(species %in% tpp.vulnerable)])]$spID, spar[which(spar[,(species %in% tpp.nearthreatened)])]$spID)


## for BAU
sum(spar[,(change_bau > 0.9 & change_bau < 1)], na.rm=T)
spar[which(spar[,(change_bau > 0.9 & change_bau < 1)])]$spID
bau.sp <- spar[which(spar[,(change_bau > 0.9 & change_bau < 1)]),]
## Find additional species to tpp.species
bau.sp[which(!(bau.sp$species %in% tpp.sp$species))]$species #is least concern, so ignore


## Extract info from results for selected species - on boab
# load("./output/results_til.RData")
# selected.results <- results[sp.ind]
# save(selected.results, file = paste0("output/", "selected.results_", country_abbr, ".RData"))
load("./output/selected.results_til.RData")

sp.names.full <- c(tpp.endangered, tpp.vulnerable, tpp.nearthreatened)
sp.names.full <- sub(" ","_", sp.names.full)
sp.names.short <- sub("_.*","", sp.names.full) # for second part: sub(".*_","", sp.names.full)
sp.names.short <- sapply(sp.names.short, tolower)

for (i in 1:6) {
  assign(paste0(sp.names.short[i], "_000"), selected.results[[i]][[4]][[1]])
  assign(paste0(sp.names.short[i], "_bau"), selected.results[[i]][[4]][[2]])
  assign(paste0(sp.names.short[i], "_tpp"), selected.results[[i]][[4]][[3]])
  
  par(mfrow=c(1,3))
  plot(get(paste0(sp.names.short[i], "_000")), main = paste0(sp.names.full[i], "_000"))
  plot(get(paste0(sp.names.short[i], "_bau")), main = paste0(sp.names.full[i], "_bau"))
  plot(get(paste0(sp.names.short[i], "_tpp")), main = paste0(sp.names.full[i], "_tpp"))
  
  # writeRaster(get(paste0(sp.names.short[i], "_000")), paste0("figures/", sp.names.full[i], "_000", ".tiff"))
  # writeRaster(get(paste0(sp.names.short[i], "_bau")), paste0("figures/", sp.names.full[i], "_bau", ".tiff"))
  # writeRaster(get(paste0(sp.names.short[i], "_tpp")), paste0("figures/", sp.names.full[i], "_tpp", ".tiff"))
}




