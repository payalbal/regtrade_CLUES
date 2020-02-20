## Data processing script
rm(list = ls())
gc()

x <- c('data.table', 'rgdal', 'rgeos', 'sp', 'raster')
lapply(x, require, character.only = TRUE)

layer_path <- file.path("./data", "processed", "layers_sdm")

## Start pre-processing loop
country_abbr <- c("til", "vnm") {
  
  ## Load mask
  ## masks were pre-processed by SK in QGIS using a high resolution raster 
  reg.mask <- raster(file.path("data/processed", "masks", paste0("mask_", country_abbr, ".tif"))) 
  
  ## BIODIVERSITY DATA
  ## ------------------------------
  ## Load downloaded GBIF occurrence data
  if (country_abbr == "til") {
    gbif.raw <- fread(file.path("data/raw", "gbif", "0009081-181003121212138.csv"), header = T, na.strings=c("NA", "", " "))
  } else {
    gbif.raw <- fread(file.path("data/raw", "gbif", "0009076-181003121212138.csv"), header = T, na.strings=c("NA", "", " "))
  }
  
  
  ## Filter GBIF data and save as reduced data.table
  backbone <-fread("data/raw/gbif/gbif_backbone_taxonomy.tsv") # GBIF backbone taxonomy
  source("./scripts/0_filtergbif_data.R")
  filter_gbif_data (gbif.downloaded.data = gbif.raw, gbif.nub.taxonomy = backbone, 
                    output_folder = "data/processed", output_name = paste0("filtered_gbif_", country_abbr), 
                    domain.mask = reg.mask, 
                    start.year = 1950, end.year = 2018,
                    spatial.uncertainty.m = 1000, 
                    verbose = TRUE)
  
  rm(gbif.raw)
  gc()
  
  ## Bin species according to number of occurrencs
  ## See Merow et al. 2013
  ## See https://onlinelibrary-wiley-com.ezp.lib.unimelb.edu.au/doi/full/10.1111/j.0906-7590.2006.04700.x
  
  
  
  ## STATIC PREDICTOR DATA
  ## ------------------------------
  ## Topographic variables - Elevation, slope and roughness
  ## source: SRTM Global 30 arc secs NASA;	https://webmap.ornl.gov/ogc/wcsdown.jsp?dg_id=10008_1
  srtm <- raster(file.path("data/raw", "srtm_aus_vnm_tile29.tif"))
  srtm <- crop(srtm, reg.mask)
  names(srtm) <- "srtm"
  extent(srtm) <- extent(reg.mask)
  elevation <- mask(srtm, reg.mask)
  slope <- terrain(srtm, opt = "slope")
  roughness <- terrain(srtm, opt = "roughness")
  aspect <- terrain(srtm, opt = "aspect")
  terrain <- stack(elevation, slope, roughness, aspect)
  
  for (i in 1:3){
    r <- terrain[[i]]
    r <- projectRaster(r, reg.mask)
    r <- mask(r, reg.mask)
    print(cellStats(r, mean))
    writeRaster(r, paste0(layer_path, "/", paste(c(names(terrain)[i], country_abbr, "sta","000", "00"), collapse = "_"), ".tif"), format = "GTiff", overwrite = T)
  }
  
  rm(r, srtm, elevation, slope, roughness, aspect)
  gc()
  
  ## Soil varaibales
  ## source: Global Soil Data Task Group;	https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
  names <- c( "soilbulkdens", "soilavailablewataer", "soilcarbondens",  "soilnitro")
  soil <- stack(list.files(file.path("data/raw", "soil", "data"), pattern = "*.dat", full.names = T))
  
  for (i in 1:4){
    r <- raster(soil[[i]])
    crs(r) <- crs(reg.mask)
    r <- projectRaster(r, reg.mask)
    r <- mask(r, reg.mask)
    writeRaster(r, paste0(layer_path, "/", paste(c(names[i], country_abbr, "sta", "000", "00"), collapse = "_"), ".tif"), overwrite = T)
  }
  
  ## 'Distance to' variables
  ## sources: 
  ## Protected areas: https://www.protectedplanet.net/c/about
  ## Road network: Socioecnomic data and application centre; http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download#openModal
  ## Built-up: 
  ## Drainage, Sea, Lakes: NOAA - Global Self-consistent Hierarchical High-resolution Geography;	http://www.soest.hawaii.edu/wessel/gshhg/
  
  pa <- shapefile("data/raw/protectedarea/WDPA_Nov2018-shapefile-polygons.shp")
  pa <- crop(pa, reg.mask) # to clip vector to ratser extent
  r <- raster(extent(reg.mask))
  crs(r) <- crs(reg.mask)
  r <- projectRaster(r, reg.mask, method = "ngb")
  dd = gDistance(pa, as(r,"SpatialPoints"), byid=TRUE)
  r[] = apply(dd,1,min)
  writeRaster(r, paste0(layer_path, "_new", "/", paste(c("distpas", country_abbr, "sta", "000", "00"), collapse = "_"), ".tif"), format = "GTiff", overwrite = T)
  # pas <- raster("data/processed/layers_sdm_new/distpas_til_sta_000_00.tif") ## does not work, use the code below instead
  dir <- paste0(layer_path, "_new", "/")
  files<-list.files(dir, pattern = ".tif")
  pas <- lapply(paste0(dir, files), raster)
  plot (pas[[1]])
  
  
  
  
  distances <- stack(list.files(file.path("data/raw", "qgis_preprocessed"),  full.names = T, pattern = country_abbr))
  distances[is.na(reg.mask[])] <- NA
  ## OR? distances <- mask(distances, reg.mask)??
  names <- c("distbuiltup", "distcoastline", "distlakes", "distrivers", "distroads", "distpas")
  ## ?? pa[] <- (pa[] * -1) + 1
  for(i in 1:nlayers(stack)){
    writeRaster(stack[[i]], paste0(layer_path, "/",paste(c(names[i], country_abbr, "sta", "000", "00"), collapse = "_"), ".tif"), overwrite = T)
  }
  ## --------------------------------------------------------------------------------------------------------------------------------------------
  ## steps for pre-processing:
  ## combine shapefiles (different levels) in qgis (Levels 2-9 for rivers)
  ## rasterize in qgis (automatically subsetted from global data set), because r is too slow for that
  ## proximity in qqis (r too slow)
  ## then load here to mask NAs, align with other layers and write
  
  ## this is so qgis can read the attributes for some of the croppsed shapefiles
  ## have done this already, only needs done once, so code just for documentation purposes
  
  # name <- "WDBII_rivers_global_L2-L9"
  # country_abbr <- "vnm"
  # shp <- readOGR("/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name)
  # shp@data[,2] <- as.numeric(shp@data[,2])
  # shp@data[,1] <- as.numeric(shp@data[,1])
  # writeOGR(shp, dsn = "/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name, driver = "ESRI Shapefile", overwrite = T)
  ## load qgis processed files, and set NA
  # lakes <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/lakes_", country_abbr, ".tif"))
  # coast <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/coast_", country_abbr, ".tif"))
  # rivers <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/rivers_", country_abbr, ".tif"))
  # PA <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/PA_", country_abbr, ".tif"))
  # roads <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/PA_", country_abbr, ".tif"))
  # builtup <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/PA_", country_abbr, ".tif"))
  # names <- c("distrivers", "distlakes", "distcoastline", "PA", "distbuiltup", "distroads")
  ## --------------------------------------------------------------------------------------------------------------------------------------------
  
  
  ## Population density
  ## source: Global Rural-Urban Mapping Project, Version 1; http://sedac.ciesin.columbia.edu/data/set/grump-v1-population-density
  r <- raster(file.path("data/raw", "Pop_density", "gluds00ag.bil"))
  pop_dens <- crop(r, reg.mask, snap = "near")
  extent(pop_dens) <- extent(reg.mask)
  pop_dens <- mask(pop_dens, reg.mask)
  writeRaster(pop_dens, paste0(layer_path, "/", paste(c("popdens", country_abbr, "sta", "000", "00"), collapse = "_"), ".tif"), overwrite = T)
  rm(r).
  gc()
  
  ## Landuse 
  ## source: Global Land Cover Characterization (GLCC); https://www.usgs.gov/media/images/global-land-cover-characteristics-data-base-version-20 
  world_lu <- raster(file.path("data/raw", "LCType.tif"))
  lu <- crop(world_lu, reg.mask)
  lu <- projectRaster(lu, reg.mask, method = "ngb")
  lu <- mask(lu, reg.mask)
  rm(world_lu)
  gc()
  
  ## Reclassify land use classes
  ## Original land use classes in GLCC data: 
  ##    # 0	Water
  ##    # 1	Evergreen Needle leaf Forest
  ##    # 2	Evergreen Broadleaf Forest
  ##    # 3	Deciduous Needle leaf Forest
  ##    # 4	Deciduous Broadleaf Forest
  ##    # 5	Mixed Forests
  ##    # 6	Closed Shrublands
  ##    # 7	Open Shrublands
  ##    # 8	Woody Savannas
  ##    # 9	Savannas
  ##    # 10	Grasslands
  ##    # 11	Permanent Wetland
  ##    # 12	Croplands
  ##    # 13	Urban and Built-Up
  ##    # 14	Cropland/Natural Vegetation Mosaic
  ##    # 15	Snow and Ice
  ##    # 16	Barren or Sparsely Vegetated
  
  lu[lu[]%in%c(12,14)] <- 100       # 0 crop
  lu[lu[]%in%c(6,7,8,9,10)] <- 200  # 1 grass/shrub
  lu[lu[]%in%c(1,2,3,4,5)] <- 300   # 2 forest
  lu[lu[]%in%c(13)] <- 400          # 3 urban/artificial
  lu[lu[]%in%c(11,15,16,0)] <- 500  # 4 barren
  lu <- lu/100 - 1
  table(lu[])
  names(lu) <- "landuse"
  writeRaster(lu, paste0(layer_path, "/", paste(c("landuse", country_abbr, "luc", "000", "00"), collapse = "_"), ".tif"), format = "GTiff", overwrite = T)
  
  
  ## DYNAMIC PREDICTOR DATA
  ## ------------------------------
  ## Bioclim data - present 
  ## source: downloaded for tile29; http://www.worldclim.org/cmip5_30s
  
  ## SK's code to downlaod and merge worldclim data by tile
  # source("./R/00_functions.R")
  # tiles <- "29"
  # temp_folder <- file.path(".", "temp")
  # dir.create(temp_folder)
  # bioclim <- get_wctiles(tile = tiles, var = "bio", path = temp_folder)
  # bioclim <- merge_wctiles(bioclim)
  # names(bioclim) <- paste0("bio", c(1:19))
  # bioclim <- crop(bioclim, reg.mask)
  # bioclim <- mask(bioclim, reg.mask)
  # unlink(temp_folder, recursive = T)
  # covariates <- stack(terrain, soil, distances, pa, popdens, bioclim, lu)
  # saveRDS(covariates, file = paste0(layer_path, "covariates_", country_abbr, ".rds"))
  
  
  bionames <- gsub('.{7}$', '', list.files(file.path("data/raw", "bioclim", "present_tile29"), pattern = ".bil"))
  tile29_bio <- stack(list.files(file.path("data/raw/", "bioclim", "present_tile29"), pattern = ".bil", full.names = T))
  crs(tile29_bio) <- crs(reg.mask)
  files_sub <- crop(tile29_bio, reg.mask)
  files_sub <- mask(files_sub, reg.mask)
  for(i in 1:nlayers(files_sub)){
    writeRaster(files_sub[[i]], paste0(layer_path, "/", paste(c(bionames[i], country_abbr, "dyn", "000", "00"), collapse = "_"), ".tif"), overwrite = T)
  }
  
  rcp <- c(26, 85)
  ## Bioclim data - 2070, RCP 2.6 
  ## source: downloaded for tile29; http://www.worldclim.org/cmip5_30s
  files <- stack(list.files(file.path("data/raw", "bioclim", "he26bi70"), pattern = "tif", full.names = T, recursive = T))
  files <- crop(files, reg.mask)
  files <- mask(files, reg.mask)
  for(i in 1:nlayers(files)){
    writeRaster(files[[i]], paste0(layer_path, "/", paste(c(bionames[i], "26",country_abbr, "dyn", "000", "56"), collapse = "_"), ".tif"), overwrite = T, format = "GTiff")
  }
  
  ## Bioclim data - 2070, RCP 8.5 
  ## source: downloaded for tile29; http://www.worldclim.org/cmip5_30s
  files <- stack(list.files(file.path("data/raw", "bioclim", "he85bi70"), pattern = "tif", full.names = T, recursive = T))
  files <- crop(files, reg.mask)
  files <- mask(files, reg.mask)
  for(i in 1:nlayers(files)){
    writeRaster(files[[i]], paste0(layer_path, "/", paste(c(bionames[i], "85", country_abbr, "dyn", "000", "56"), collapse = "_"), ".tif"), overwrite = T, format = "GTiff")
  }
  
  
  
  ## SYNC NAs FOR PREDICTOR LAYERS
  ## ------------------------------
  nas <- list()
  files <- list.files(layer_path, full.names = T, pattern = country_abbr)
  pred.stack <- stack(files[grepl(paste0("(?=.*",country_abbr,")"), files, perl = TRUE)])
  sums <- calc(pred.stack, sum)
  sums_nona <- which(is.na(sums[]))
  
  for(i in 1:nlayers(pred.stack)){
    r <- pred.stack[[i]]
    r[sums_nona] <- NA
    writeRaster(r, paste0(layer_path, "/", names(pred.stack)[[i]]), format = "GTiff", overwrite = T)
    print(paste0(i, "/", nlayers(pred.stack)))
  }
  rm(r, pred.stack, sums, sums_nona)
  gc()
  
  
  
  ## CREATE PREDICTOR LAYERS FOR ALL TIME STEPS VIA LINEAR INTERPOLATION
  #-------------------------------------------------------------------------
  
  ## Load present climate data
  pres_clim <- stack(list.files(layer_path, pattern = paste(c( country_abbr,"dyn", "000", "00"), collapse = "_"), full.names = T))
  pres_clim_df <- as.data.frame(pres_clim)
  time_steps <- c(16, 26, 36, 46) # from 2015 onwards: 2030, 2040, 2050, 2060
  ras_template <- pres_clim[[1]]
  
  ## Loop through scenarios and produce linearly interpolated covariates for 10 yr time steps
  for (i in 1:length(rcp)){
    fut_clim <- stack(list.files(layer_path, pattern = paste(c(rcp[i], country_abbr, "dyn", "000", "56"), collapse = "_"), full.names = T))
    fut_clim_df <- as.data.frame(fut_clim)
    incr_df <- (fut_clim_df - pres_clim_df)/ (2070 - 2015)
    for (j in 1:length(time_steps)){
      ts_clim_df <- pres_clim_df + time_steps[j] * incr_df
      for (k in 1:nlayers(pres_clim)){
        ras_template[] <-  ts_clim_df[,k]
        writeRaster(ras_template, paste0(layer_path,"/", paste(c(bionames[k], "_", rcp[i], country_abbr, "dyn", "000", time_steps[j]), collapse = "_")), format = "GTiff", overwrite = T)
      }
    } # end time_steps
  } # end rcp
} # end country_abbr

rm(backbone)
gc()
