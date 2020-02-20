#-------------------------------------------#
#---I. PROCESSING LAYERS FOR USE IN CLUES####
#-------------------------------------------#
rm(list = ls())
gc()
setwd(file.path("~", "Dropbox", "trade-biodiversity discovery", "discovery_paper_1"))



library("raster")
library("SDMTools")
getwd()

#load functions
source(file.path("R", "dp1_functions.R"))

#1) LOAD SDM LAYERS

layer_path <- file.path("processed_data", "layers_sdm")
mask_path <- file.path("processed_data", "masks")
country_abbr <- "vnm" #alwyas vnm, because we run clues only for vnm)

clues_path <- file.path("processed_data", "layers_clues")
crs_clues <- "+proj=utm +zone=49 +datum=WGS84 +no_defs +ellps=WGS84 +units=m" #Target UTM CRS

files <- list.files(layer_path, full.names = T) #all files

#filter the ones we need
bio_stack <- stack(files[grepl(paste0("(?=.*vnm)(?=.*bio)(?=.*000_00)"), files, perl = TRUE)]) # Stacks all files with vnm, bio and 000_00 in the file name
sta_stack <- stack(files[grepl(paste0("(?=.*vnm)(?=.*sta)(?=.*000_00)"), files, perl = TRUE)])
sta_stack <- sta_stack[[-which(grepl("popdens|protecteda", names(sta_stack)))]] # Removes files with poplins or protected in their name

mask <- raster(file.path(mask_path, "mask_vnm.tif"))
mask <- mask-1

lu <- raster(files[grepl("landuse_vnm_luc_000_00", files)])
pa <- raster(files[grepl(paste0("protecteda_vnm"), files)])
pa[pa ==1] <- -9998
cov_stack <- stack(bio_stack, sta_stack)
crs(cov_stack) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#2) CORELLATION ANALYSIS
cov_df <- getValues(cov_stack)
rows <- nrow(cov_df)
cov_df <- na.omit(cov_df)
subs <- sample(1:nrow(cov_df), 10000)
cors <- cor(cov_df[subs,])
preds <- rownames(reduce_predset(cors))

save(preds, file = file.path(clues_path, paste0("preds_clues_", country_abbr, ".RData")))

#3) REPROJECTING RASTERS
maskpalu <- stack(mask, pa, lu)
maskpalu <- projectRaster(maskpalu, crs = crs_clues, res = 1000, method = "ngb")
sta_stack <- sta_stack[[which(names(sta_stack)%in%preds)]]  # removes from stack the predictor variables that were eliminated in the correlation analysis step
bio_stack <- bio_stack[[which(names(bio_stack)%in%preds)]]
cov_stack <- stack(bio_stack, sta_stack)
crs(cov_stack) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "


#reproject by raster, because of RAM limits
cov_stack_repr <- stack()
for (i in 1:nlayers(cov_stack)){
  s <- projectRaster(cov_stack[[i]], maskpalu)
  cov_stack_repr <- stack(cov_stack_repr, s)
  print(i)
}

cov_stack <- stack(cov_stack_repr[[1:nlayers(bio_stack)]], cov_stack_repr, maskpalu)
names(cov_stack)
cov_stack <- syncNA1(cov_stack)

#determine which cells are NA
# This step omits NAs from data frame and finds the indices for rows in the data frame with the NAs
cov_df <- as.data.frame(cov_stack)
all_inds <- 1:nrow(cov_df)
cov_df <- na.omit(cov_df)
inds <- as.numeric(as.character(rownames(cov_df)))
na_inds <- all_inds[!all_inds%in%inds]

#4)WRITE CURRENT TIME STEP
names(cov_stack)

#WRITE sc1 .fil files 
names <- c(paste0("Sc1gr",0:(nlayers(bio_stack)-1), ".0.asc"),
           paste0("Sc1gr",0:(nlayers(cov_stack)-(4+nlayers(bio_stack))),".fil.asc"),
           paste0("region", country_abbr, ".fil.asc"),
           paste0("regionpa", country_abbr, ".fil.asc"),
           "cov_all.0.asc"
)

var_names <- data.frame(names(cov_stack), names) # original names and corresponding CLUES names for selected (from correlation step) precitor variables 
write.table(var_names, file.path(clues_path, paste0("preds_clues_", country_abbr)))

#3a) BioClim variable
rm(bio_stack, sta_stack, cov_df, cov_stack_repr, pa)

gc()

for (i in 1:nlayers(cov_stack)){
  r <- cov_stack[[i]]
  r <- asc.from.raster(r)
  write.asc(r, file = file.path(clues_path, names[i]))
  #writeRaster(r, file = paste0(clues_path, names[i]), format = "ascii", overwrite = T)
}


#4) DYNAMIC VARIABLES
var_names <- read.table(file.path(clues_path, paste0("preds_clues_", country_abbr)))

i <- 1
j <- 1
var <- gsub('.{16}$', '', var_names[,1])
time_steps <- c(16, 26, 36, 46, 56)
mask <- raster(file.path(clues_path, "cov_all.0.asc"))
crs(mask) <- crs_clues

for (j in 1:length(time_steps)){
  fut_stack <- stack(files[grepl(paste0("(?=.*vnm)(?=.*bio)(?=.*_", time_steps[j], ")"), files, perl = TRUE)])
  fut_stack <- fut_stack[[which(grepl(paste(var, collapse = "|"), names(fut_stack)))]]
  crs(fut_stack) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  fut_stack <- projectRaster(fut_stack , mask)
  fut_stack <- syncNA1(stack(fut_stack, mask))[[1:nlayers(fut_stack)]]
  for (k in 1:nlayers(fut_stack)){
    r <- fut_stack[[k]]
    r <- asc.from.raster(r)
    print(paste("writing:", k-1, j))
    write.asc(r, file.path(clues_path, paste0("Sc1gr",k-1,".", j, ".asc")))
  }
}

#REGRESSION
# This step does a simple regression analysis (lu class ~ variables) and produces the alloc.reg1 file 
# with significant regression parameters as input for clues 
lu <- raster(file.path(clues_path, "cov_all.0.asc"))
landuse_val <- lu[]

landuse_ma <- matrix(0, nrow = length(landuse_val), ncol = 5)
for (i in 1:ncol(landuse_ma)){
  landuse_ma[which(landuse_val == i-1), i] <- 1
  print(i)
}

reg_stack <- stack(list.files(clues_path, recursive = F, pattern = ".fil", full.names = T))
reg_stack <- reg_stack[[-which(grepl("region", names(reg_stack)))]]
driver_id <- as.character(0:(nlayers(reg_stack)-1))
driver_id <- as.numeric(sort(driver_id))

reg_df <- as.data.frame(reg_stack)
reg_df <- cbind(landuse_ma,reg_df)
names(reg_df)[1:5] <- c(paste0("lu", 0:4))
reg_df <- na.omit(reg_df)
subs <- sample(1:nrow(reg_df), 100000)

coeffs.df <- list()
pvalues <- list()

for (i in 1:5){
  model <- glm(as.formula(paste(colnames(reg_df)[i], "~", paste(names(reg_df[,-c(1:5)]), collapse="+"))), data = reg_df[subs,], family = "binomial")
  coeffs.df[[i]] <- round(coefficients(model),8)
  pvalues[[i]] <- round(summary(model)$coefficients[,4],8)
}

lu_id <- 0:4
names <- names(pvalues[[5]])[-1]
sink(paste0(clues_path, "/alloc1.reg"))
for (i in 1:length(lu_id)){
  intercept <- coeffs.df[[i]][1]
  coeffs <- coeffs.df[[i]][-1][which(pvalues[[i]][-1] < 0.05)]
  id <- driver_id[which(names%in%names(coeffs))]
  cat(lu_id[i],"\r","\n","\t",intercept,"\r","\n", length(coeffs), "\r", "\n", sep = "")
  for (j in sort(id, index.return = T)$ix){
    cat("\t",coeffs[j], " ", id[j],"\n","\r", "\n", sep = "")
  }
}
sink()
