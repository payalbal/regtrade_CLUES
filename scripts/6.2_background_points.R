require("raster")


######################################
#FINDING CORRELATIONS BETWEEN PARS####
######################################
source("~/discovery_paper_1/R/dp1_functions.R")
country_abbr <- "til"
bias_preds <- c("distbuiltup", "distroads", "protecteda", "popdens", "roughness")
not_used <- c("srtm", "popdens")
layer_path <- "~/discovery_paper_1/data/layers_sdm"
files <- list.files(layer_path, full.names = T)
covariates <- stack(files[grepl(paste0("(?=.*",country_abbr,")(?=.*000_00)"), files, perl = TRUE)])
covariates <- covariates[[-which(grepl(paste(c(bias_preds, not_used), collapse = "|"), names(covariates)))]]
names(covariates) <- gsub('.{15}$', '', names(covariates))
cov_df <- getValues(covariates)
cov_df <- na.omit(cov_df)
subs <- sample(1:nrow(cov_df), size = 10000)
preds <- rownames(reduce_predset(cor(cov_df[subs,])))

save(preds, file = paste0("~/discovery_paper_1/output/preds_", country_abbr, ".RData"))

#Sample background poiints
dens_rast <- raster(list.files("~/discovery_paper_1/data/layers_bias", pattern = country_abbr, full.names = T))
inds <- which(!is.na(dens_rast[]))
probs <- round(dens_rast[inds],3)
inds_all <- sample(inds, p = probs, size = 10000)
save(inds_all, file = paste0("~/discovery_paper_1/output//bgp_", country_abbr, "_bias.RData"))
plot(probs)

