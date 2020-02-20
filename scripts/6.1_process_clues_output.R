require(raster)
m <- raster("~/Dropbox/discovery_paper_1/processed_data/landuse_vnm_bau_000_56.asc")
n <- raster("~/Dropbox/discovery_paper_1/processed_data/landuse_vnm_tpp_000_56.asc")
crs(m) <- "+proj=utm +zone=49 +datum=WGS84 +no_defs +ellps=WGS84 +units=m"
crs(n) <- "+proj=utm +zone=49 +datum=WGS84 +no_defs +ellps=WGS84 +units=m"
mask <- raster("~/Dropbox/discovery_paper_1/processed_data/masks/mask_vnm.tif")

m2 <- projectRaster(m, mask, method = "ngb")
n2 <- projectRaster(n, mask, method = "ngb")

table(n2[])

writeRaster(m2, "~/Dropbox/discovery_paper_1/processed_data/layers_sdm/landuse_vnm_luc_bau_56.tif", format = "GTiff")
writeRaster(n2, "~/Dropbox/discovery_paper_1/processed_data/layers_sdm/landuse_vnm_luc_tpp_56.tif", format = "GTiff")

str(results[[1]])

#first object: cross-validated model results
#second object:
results[[5]][[3]]


