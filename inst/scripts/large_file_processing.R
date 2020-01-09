# Creating the largest satellite data covariate files (not stored on Github)

# Might take an hour or so to process
outfolder <- here::here("inst","extdata", "processed")

band1_5000 <- CovariateCreationFunction(noPixEachDim = 5000, covariates.type="band1")
saveRDS(object = band1_5000, file = file.path(outfolder,"cov_grid_band1_5000.rds"))

band2_5000 <- CovariateCreationFunction(noPixEachDim = 5000, covariates.type="band2")
saveRDS(object = band2_5000, file = file.path(outfolder,"cov_grid_band2_5000.rds"))
