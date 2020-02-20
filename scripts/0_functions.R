## Functions used in discovery_paper_1 analyses

## Sync NA
## Author: Simon Kapitza
syncNA1 <- function (x) 
{
  val <- getValues(x)
  NA.pos <- unique(which(is.na(val), arr.ind = T)[, 1])
  val[NA.pos, ] <- NA
  x <- setValues(x, val)
  return(x)
}

## Reduce predictor set automatically - by Simon Kapitza
reduce_predset <- function(cors = matrix(), thresh = 0.7) {
  while (min(abs(cors[abs(cors) >= thresh])) != 1){
    values <- cors[which(abs(cors) > thresh)]
    corellated <- which(abs(cors) > thresh)
    values[values ==1] <- NA
    corellated[which(values== max(values, na.rm = T))]
    rows_highest_cor <- which(cors == max(values, na.rm = T), arr.ind = T)[,1]
    cors_cur <- abs(cors[rows_highest_cor,])
    '%ni%' <- Negate('%in%')
    m1 <- max(cors_cur[1,][cors_cur[1,]%ni%c(max(values, na.rm = T),1)])
    m2 <- max(cors_cur[2,][cors_cur[2,]%ni%c(max(values, na.rm = T),1)])
    out <- ifelse(m1 > m2, 1, 2)
    cors <- cors[-which(colnames(cors) == names(rows_highest_cor)[out]), -which(colnames(cors) == names(rows_highest_cor)[out])]
    nrow(cors)
  }
  return(cors)
}


## Interpolation and output of results
## Author: Simon Kapitza
approx_fun <- function(x, n){
  vec <- c(seq(0,28, 2),seq(30,54, 4))
  x <- x
  out <- approx(vec, x, method = "linear", n = n)$y[1:n]
  return(out)
}

# approx_fun2 <- function(x, names){
#   vec <- as.numeric(sub("X", "", names)) -2015
#   out <- approx(vec, x, method = "linear", n = max(vec))$y
#   return(out)
# }


## https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
## existing_agr: matrix with 4 layers: counts, means, ssq, var; nv: layer with nbew preds
## continue here.
update_stats <- function(existing_aggr = matrix(), new_value = numeric()){
  v <- existing_aggr
  nv <- new_value
  v[,1] <- v[,1] + 1
  d <- nv - v[,2]
  v[,2] <- v[,2] + d / v[,1]
  d2 <- nv - v[,2]
  v[,3] <- v[,3] + d * d2 #sum of squares
  v[,4] <- v[,3] / (v[,1] - 1)
  return(v)
}


## Calculate area-corrected richness for each raster in the stack
## Author: Payal Bal
## NOTES Area-corrected richness is the ratio of value in a pixel vs the sum of values across the area for a species. It gives a relative value that can be compared across species (unlike the relative likelihoods). So we can sum this for species. This is still an index of area corrected richness and not a measure of the richness directly. 
area_corrted_richnes <- function (inputstack){
  output <- stack()
  if(class(inputstack) != "RasterStack"){
    return(NA)
  } else {
    for (i in 1:dim(inputstack)[3]){
      output <- stack(output, inputstack[[i]]/cellStats(inputstack[[i]], stat = "sum", na.rm = T))
      ## values come out the same as inputstack[[i]]@data@values/sum(inputstack[[i]]@data@values, na.rm = T)
      ## check by subsetting i=1 and comparing min and max with that of resultant raster
    }
    return(output)
  }
}

## Save Raster
## Author: Simon Kapitza
writeToDisk <- function(covariates, folder){
  dir.create(folder, recursive = T)
  writeRaster(covariates, filename = file.path(folder, paste0(names(covariates), ".tif")), bylayer = T, driver = "GTiff", overwrite = T)
}


## Download WorldClim data
## Author: Simon Kapitza
get_wctiles <- function(tiles, var, path, ...){
  
  if(missing(path)){
    path <- getwd()
  }
  
  rs <- raster(nrows = 5, ncols = 12, xmn = -180, xmx = 180, ymn = -60, ymx = 90)
  rs[] <- 1:length(rs)
  tiles_names <- c(paste0(0, 0:11), paste0(1, 0:11), paste0(2, 0:11), paste0(3, 0:11), paste0(4, 0:11))
  points <- xyFromCell(rs, which(tiles_names%in%tiles))
  message("Loading required worldclim tiles...")
  biotiles <- list()
  for (i in 1:nrow(points)){
    biotiles[[i]] <- getData(name = "worldclim", var = var, path = path, res = 0.5, lon = points[i,1], lat = points[i,2], ...)
  }
  message(paste0("Data successfully downloaded to ", path))
  biotiles
}


## Merge downloaded WorldClim data by tiles
## Author: Simon Kapitza
merge_wctiles <- function(biotiles){
  
  if(!is.list(biotiles)) stop("Please provide list with stacks for each tile")
  
  if (length(biotiles) == 1){
    out <- biotiles[[1]]
    message("Only 1 tile detected, merging not necessary")
  }
  
  if(length(biotiles) > 1){
    bionames <- names(biotiles[[1]])
    out <- list()
    message("Merging tiles...")
    for (i in 1:nlayers(biotiles[[1]])){
      b <- biotiles[[1]][[i]]
      message(paste0(substr(bionames[i], 1, nchar(bionames[i])-3)," | " , i, "/", length(bionames)))
      
      for(j in 2:length(biotiles)){
        b <- raster::merge(b, biotiles[[j]][[i]])
      }
      
      out[[i]] <- b
    }
    out <- stack(out)
    names(out) <- bionames
  }
  out
}


## Calculate distance to feature and create a raster based on mask
## Author: Roozbeh Valavi (Feb 2019)
rasterDistance <- function(feature, rastermask){
  require(raster)
  require(sf)
  require(progress)
  p <- st_as_sf(rasterToPoints(rastermask, spatial = TRUE))
  p$indx <- st_nearest_feature(st_geometry(p), st_geometry(feature))
  pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
                                   total=nrow(p), clear=FALSE, width=75) # add progress bar
  for(i in 1:nrow(p)){
    p$dist[i] <- st_distance(st_geometry(p[i,]), st_geometry(feature[p$indx[i],]))
    pb$tick() # update progress bar
  }
  output <- raster::rasterize(p, rastermask, field = "dist")
  return(output)
}




