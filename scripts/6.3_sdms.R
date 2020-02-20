##

system("ps")
system("pkill -f R")  # kills all R proccesses

# ## UNCOMMENT FOR REINSTALLING MAXENT FOR BOAB ---------------------------------------- ##
# remove.packages("dismo")
# remove.packages("rJava")
# quit(save = "no", status = 0, runLast = TRUE)
# install.packages("dismo")
# install.packages("rJava")
# quit(save = "no", status = 0, runLast = TRUE)
# 
# ## IMPORTANT STEPS:
# ## 1. download earlier version of maxent jar (3.3.3k) from here: https://github.com/mrmaxent/Maxent/tree/master/ArchivedReleases
# ## 2. copy maxent.jar file to the specified location, here in home folder
# ## 3. Now run the command to file.copy below:
# file.copy("~/maxent.jar", "~/.r-dir/R/library/dismo/java/maxent.jar") 
# file.exists("~/.r-dir/R/library/dismo/java/maxent.jar")
# 
# ## FOR FIRST TIME BOAB USERS/OR IF {rgdal} DOESN"T LOAD DUE TO GDAL ISSUES
# ## There might be problems loading the rgdal package required for the raster package
# ## This is because of its dependency on GDAL (not an R package)
# ## In this case, get in touch with Nick/Jian/Casey to help you out
# 
# ## ----------------------------------------------------------------------------------- ##

## NOTE: If rJava fails to load see: https://github.com/MTFA/CohortEx/wiki/Run-rJava-with-RStudio-under-OSX-10.10,-10.11-(El-Capitan)-or-10.12-(Sierra)

rm(list = ls())
gc()

library("pacman")
p_load("data.table", "raster", "doParallel", "dismo", "rJava","iterators", "parallel", "foreach")
maxent()


## SET WORKING DIRECTORY FOR BOAB
setwd("~/discovery_paper_1")

start_time <- Sys.time()

## LOAD DATA
## 1.a) Data loading parameters
country_abbr <- "til" #til here, is changed automatically in code for preds
layer_path <- "./data/processed/layers_sdm/"
bias_preds <- c("distbuiltup", "distroads", "protecteda", "popdens", "roughness")
not_used <- c("srtm", "popdens")

## 1.b) Bias data
## i) Bias layer
dens_rast <- raster(list.files("./data/processed/layers_bias/", pattern = country_abbr, full.names = T))

## ii) Bias points
load(paste0("./output/bgp_", country_abbr, "_bias.RData"))
inds <- inds_all

## 1.c) Observation data
obs <- read.csv(paste0("./data/processed/aves_processed_",country_abbr,".csv"))
species <- sort(as.character(unique(obs$species)))
str(obs)

## 1.d) Covariates
load(paste0("./output/preds_", country_abbr, ".RData"))
files <- list.files(layer_path, full.names = T)
sta <- stack(files[grepl(paste0("(?=.*",country_abbr,")(?=.*sta)(?=.*000_00)"), files, perl = TRUE)])
dyn <- stack(files[grepl(paste0("(?=.*",country_abbr,")(?=.*dyn)(?=.*000_00)"), files, perl = TRUE)])
luc <- stack(files[grepl(paste0("(?=.*",country_abbr,")(?=.*landuse)(?=.*000_00)"), files, perl = TRUE)])
covariates <- stack(sta, dyn, luc)
cov_names <- names(covariates) <- gsub('.{15}$', '', names(covariates))
covariates <- covariates[[-which(grepl(paste(c(bias_preds, not_used), collapse = "|"), names(covariates)))]]
?? ... covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]


## MODEL & PREDICT
## 2.a) Model parameters
factors <- "landuse"
scens <- c("000", "bau", "tpp")
output <- paste0("./output/")
logfile <- paste0(output, country_abbr, "_log.txt")
writeLines(c(""), logfile)
l <- length(species)
# spec <- round(runif(l, 1, length(species)))                         #subset step
sp.stack <- list()  # list to store 3 rasters (one per scenario) for each species


## 2.b) Load Cluster for Modelling
cl <- makeCluster(30) # or detectCores() to use the default number of cores available
registerDoParallel(cl)
results <- foreach(i = 1:l, .packages = c("dismo", "raster")) %dopar% {
  # results <- foreach(i = spec, .packages = c("dismo", "raster")) %dopar% {       #subset step
  cov <- covariates
  cat(paste("Starting model",i,"\n"), file = logfile, append = T)
  subs <- obs[obs$species == species[i],]
  bgp <- SpatialPoints(xyFromCell(dens_rast, cell = inds))
  
  ## 2.c) Determine covariate importance
  ## i) Initial model run
  mod <- tryCatch(maxent(cov,
                         p= subs[,c(1,2)], 
                         a= bgp, 
                         factors = factors, 
                         args= c("randomtestpoints=25",
                                 "-J", "-p", "-h", 
                                 "threshold=FALSE")),
                  error = function(e) NA)
  
  
  ## ii) Determine covariates with perm importance < 1 and remove from covariate set
  if(class(mod) != "MaxEnt"){
    return(NA)
  } else {
    
    res <- data.frame("names" = as.character(rownames(mod@results)), "results" = as.numeric(mod@results))
    perm <- res[grep("permutation.importance", res$names),]
    kickout <- perm$names[which(perm$results < 1)]
    kickout <- gsub('.permutation.importance', '', kickout)
    
    if(length(kickout) > 0){
      covariates.new <- cov[[-which(grepl(paste0(kickout, collapse = "|"), names(cov)))]]
      names.new <- names(covariates.new)
    }else{
      covariates.new <- cov
      names.new <- names(covariates.new)
    }
    
    factors.new <- factors[factors%in%names(covariates.new)]
    
    ## 2.d) Cross-validated model
    mod.new <- tryCatch(dismo::maxent(covariates.new, 
                                      p= subs[,c(1,2)], 
                                      a= bgp, 
                                      factors = factors.new, 
                                      args= c("replicatetype=crossvalidate", 
                                              "replicates=4", 
                                              "-J","-P", "-p", "-h", 
                                              "threshold=FALSE")), 
                        error = function(e) NA)
    
    modresults <- mod.new@results
    rm(mod, mod.new, cov)
    gc()
    
    ## 2.e) Final model
    mod.final <- tryCatch(dismo::maxent(covariates.new, 
                                        p= subs[,c(1,2)], 
                                        a= bgp, 
                                        factors = factors.new,
                                        args= c("randomtestpoints=25",
                                                "-J","-P", "-p", "-h", 
                                                "threshold=FALSE")), 
                          error = function(e) NA)
    
    if(class(mod.final) != "MaxEnt"){
      return(NA)
    } else {
      
      ## 2.f) Predictions
      cat(paste("Starting prediction",i,"\n"), file = logfile, append = T)
      
      thresh_maxSSS <- mod.final@results[names(mod.final@results[,1])%in%"Maximum.test.sensitivity.plus.specificity.logistic.threshold"]
      
      if (country_abbr == "til"){
        country_abbr_preds <- "vnm"
        sta <- stack(files[grepl(paste0("(?=.*vnm)(?=.*sta)(?=.*000_00)"), files, perl = TRUE)])
        area_maxSSS <- numeric()
      }
      
      for (j in 1:length(scens)){
        if(scens[j] == "000" ) {
          dyn <- stack(files[grepl(paste0("(?=.*",country_abbr_preds,")(?=.*dyn)(?=.*000_00)"), files, perl = TRUE)])
          luc <- stack(files[grepl(paste0("(?=.*",country_abbr_preds,")(?=.*landuse)(?=.*000_00)"), files, perl = TRUE)])
          
        } else if(scens[j] == "bau") {
          dyn <- stack(files[grepl(paste0("(?=.*",country_abbr_preds,")(?=.*dyn)(?=.*000_56)"), files, perl = TRUE)])
          luc <- stack(files[grepl(paste0("(?=.*",country_abbr_preds,")(?=.*landuse)(?=.*bau)"), files, perl = TRUE)])
          
        } else if(scens[j] == "tpp") {
          dyn <- stack(files[grepl(paste0("(?=.*",country_abbr_preds,")(?=.*dyn)(?=.*000_56)"), files, perl = TRUE)])
          luc <- stack(files[grepl(paste0("(?=.*",country_abbr_preds,")(?=.*landuse)(?=.*tpp)"), files, perl = TRUE)])
        }
        
        covariates.new <- stack(sta, dyn, luc)
        
        names(covariates.new) <- cov_names
        covariates.new <- covariates.new[[which(names(covariates.new)%in%names.new)]]
        mapout <- dismo::predict(mod.final, covariates.new)
        area_maxSSS <- c(area_maxSSS, length(which(mapout[] >= thresh_maxSSS)))
        
        if(scens[j] == "000" ) {sp.stack[[1]] <- mapout}
        else if(scens[j] == "bau") {sp.stack[[2]] <- mapout}
        else if(scens[j] == "tpp") {sp.stack[[3]] <- mapout}
      }
    }
  }
  
  rm(dyn, luc, covariates.new, mod.final)
  gc()
  list(modresults, area_maxSSS, i, sp.stack)
}

save(results, file = paste0(output, "results_", country_abbr, ".RData"))


## EXTRACT AND SAVE OUTPUT DATA --- RUNNING THIS SECTION WILL REWRITE SAVED FILES IN OUTPUT FOLDER ---- ##
# country_abbr <- "til" #til here, is changed automatically in code for preds
# load(paste0("output/", "results_", country_abbr, ".RData"))

## Species areas
spar <- sapply(results, "[", 2)
spar <- as.data.table(do.call(rbind,spar))
spar <- cbind(1:761, species, spar)
names(spar) <- c("spID", "species", "pr", "bau", "tpp")
# na.omit(spar, invert=TRUE will return the rows with NAs)
# can also use: spar <- spar[complete.cases(spar*0)]
## CLARIFY: why are there NA values?
na.ind <- unique(which(is.na(spar), arr.ind=TRUE)[,1]) # gives indices of rows with NAs (i.e. for 76--724 = 37 species)
# 'unique is used in the above expression as otherwise row indices are repeated for every column with NAs. 
spar <- na.omit(spar)
spar <- as.data.frame(spar)
save(spar, file = paste0(output, "sparea_", country_abbr, ".RData"))


## Summed rasters by sceanrio
source("R/dp1_functions.R")

temp <-  sapply(results, "[", 4)
temp <- temp[-na.ind] # remove rasters for species with NA values

# rm(results)
gc()

vnm.stack000 <- sapply(temp, "[", 1)
vnm.stack000 <- stack(vnm.stack000)  # list or stack have the same size
sum_000 <- sum(vnm.stack000, na.rm = T)
sum_areaweighted_000 <- sum(area_corrted_richnes(vnm.stack000), na.rm = T) ## see function for area_corrected_richness

vnm.stackbau <- sapply(temp, "[", 2)
vnm.stackbau <- stack(vnm.stackbau)
sum_bau <- sum(vnm.stackbau, na.rm = T)
sum_areaweighted_bau <- sum(area_corrted_richnes(vnm.stackbau), na.rm = T)

vnm.stacktpp <- sapply(temp, "[", 3)
vnm.stacktpp <- stack(vnm.stacktpp)
sum_tpp <- sum(vnm.stacktpp, na.rm = T)
sum_areaweighted_tpp <- sum(area_corrted_richnes(vnm.stacktpp), na.rm = T)

save(vnm.stack000, file = paste0("output/", "stack000_", country_abbr, ".RData"))
save(vnm.stackbau, file = paste0("output/", "stackbau_", country_abbr, ".RData"))
save(vnm.stacktpp, file = paste0("output/", "stacktpp_", country_abbr, ".RData"))

summed.out <- stack(sum_000, sum_areaweighted_000, sum_bau, sum_areaweighted_bau, sum_tpp, sum_areaweighted_tpp)
vec <- c("sum_000", "sum_areaweighted_000", "sum_bau", "sum_areaweighted_bau", "sum_tpp", "sum_areaweighted_tpp")
writeRaster(summed.out, filename = paste0("output/", vec, ".tif"), format="GTiff", overwrite = TRUE, bylayer = TRUE)

rm(temp)
gc()


# ## Raster stacks as matrices by scenario
# ## ref: https://stackoverflow.com/questions/19833784/how-to-extract-values-from-rasterstack-with-xy-coordinates 
# vnm.mat000 <- rasterToPoints(vnm.stack000)  ## mapout@data@values for each species as a matrix,
# ## nrwos = ncelss, ncols = nstacks (i.e. nspecies) + 1 lat col + 1 long col
# colnames(vnm.mat000)[3:dim(vnm.mat000)[2]] <- rev(species[-na.ind])
# # colnames(vnm.mat000)[3:dim(vnm.mat000)[2]] <- rev(c(species[spec]))             #subset step
# 
# vnm.matbau <- rasterToPoints(vnm.stackbau)
# colnames(vnm.matbau)[3:dim(vnm.matbau)[2]] <- rev(species[-na.ind])
# # colnames(vnm.matbau)[3:dim(vnm.matbau)[2]] <- rev(c(species[spec]))             #subset step
# 
# vnm.mattpp <- rasterToPoints(vnm.stacktpp)
# colnames(vnm.mattpp)[3:dim(vnm.mattpp)[2]] <- rev(species[-na.ind])
# # colnames(vnm.mattpp)[3:dim(vnm.mattpp)[2]] <- rev(c(species[spec]))             #subset step
# 
# save(vnm.mat000, file = paste0("output/", "mat000_", country_abbr, ".RData"))
# save(vnm.matbau, file = paste0("output/", "matbau_", country_abbr, ".RData"))
# save(vnm.mattpp, file = paste0("output/", "mattpp_", country_abbr, ".RData"))

# rm(vnm.stack000, vnm.stackbau, vnm.stacktpp, vnm.mat000, vnm.matbau, vnm.mattpp)
gc()

end_time <- Sys.time()
end_time - start_time
