

#' Function for creating a dense grid with all covariate values
#'
#' Transforms a tif file with lonlat values to a dense grid with the same coordinate system as the seal data (created with the lb2xykm function),
#' with mean values imputed where satellite image value is unknown. The output is a spatstat::im variable.
#'
#' @param xkmRange Vector of size 2, corresponding to the x-range of the mesh gotton with the standard input variables
#' @param ykmRange Vector of size 2, corresponding to the y-range of the mesh gotton with the standard input variables
#' @param covariates.type String indicating the type of covariate to create (either "band1" or "band2")
#' @param chunksPerDim Numeric. The number of chunkgs to split the data into when transforming the covariates
#' @param noPixEachDim Numeric. The number of pixels to use per dimension.
#' @param expandBy Numeric, indicating how many pixels the chunks should overlap with.
#' @return Returns an image objects (which should be saved to file)
#' @import rgdal
#' @import akima
#' @import fields
#' @import spatstat
#' @export

CovariateCreationFunction <- function(satelliteFile = system.file("extdata", "original","reflectance_0.0025deg_grid_modis_20120328_1310.tif", package = "SealCoxProcess"),
                                      sealFile = system.file("extdata", "original","WestIce2012.csv", package = "SealCoxProcess"),
                                      xkmRange = c(-65.33908,  67.99593), #
                                      ykmRange = c(-103.87783, 85.35911), #
                                      covariates.type, # "band1", "band2",
                                      chunksPerDim = 5,
                                      noPixEachDim = 1000,
                                      expandBy = 5){
  ## Load the covariates
  covariates <- rgdal::readGDAL(satelliteFile)

  ## Extracts the covariates
  covFullCoordLatLon <- sp::coordinates(covariates)


  ## Defines the reference points for the transformation (the same used for transforming the original seals data)
  seals0=read.csv(file=sealFile)
  lon0=mean(seals0$lon)
  lat0=mean(seals0$lat)

  covFullCoordxykm <- lb2xykm(covFullCoordLatLon[,1],covFullCoordLatLon[,2],lon0,lat0)

  ## Defines the gridpoins we are mapping covariates to initially
  #Preset#xkmRange <- range(mesh$loc[,1])
  #Preset#ykmRange <- range(mesh$loc[,2])
  xkmVal <- seq(xkmRange[1],xkmRange[2],length.out=noPixEachDim)
  ykmVal <- seq(ykmRange[1],ykmRange[2],length.out=noPixEachDim)
  gridPoints <- list(x=xkmVal, y=ykmVal)

  ## Creates an "image" for the covariate type to be used (this is perhaps a bit tedious, could be handled better, but works)
  if (covariates.type=="band1"){
    gridcov0 <- fields::as.image(covariates$band1, x = covFullCoordxykm, grid = gridPoints)
  }
  if (covariates.type=="band2"){
    gridcov0 <- fields::as.image(covariates$band2, x = covFullCoordxykm, grid = gridPoints)
  }

  ## Need to fill in the interior NA values. Need to do this in chunks and combine them afterwords.
  valPerChunkEachDim <- noPixEachDim/chunksPerDim
  gridcovNew <- gridcov0
  k=0
  for (i in 1:chunksPerDim){
    for (j in 1:chunksPerDim){
      xDim <- 1:valPerChunkEachDim+(i-1)*valPerChunkEachDim
      yDim <- 1:valPerChunkEachDim+(j-1)*valPerChunkEachDim

      ## Expands the dimension of the x and y to not get internal boundary effects when
      xDimExpanded <- (max(min(xDim)-expandBy,1)):(min(max(xDim)+expandBy,noPixEachDim))
      yDimExpanded <- (max(min(yDim)-expandBy,1)):(min(max(yDim)+expandBy,noPixEachDim))

      ## Maps to these
      xo <- gridcov0$x[xDim]
      yo <- gridcov0$y[yDim]

      ## Maps from these
      xx <- gridcov0$x[xDimExpanded]
      yy <- gridcov0$y[yDimExpanded]
      zz <- gridcov0$z[xDimExpanded,yDimExpanded]

      nas <- !is.na(zz)
      gridcovNew$z[xDim,yDim] <- akima::interp(x=xx[row(zz)[nas]],y=yy[col(zz)[nas]],z=zz[nas],xo=xo,yo=yo,linear=TRUE,extrap=FALSE)$z  # Extracpolation does not work...
      k=k+1
      print(paste((k/chunksPerDim^2)*100," % covariate extraction complete",sep=""))
    }
  }

  ## Finally we fill inn the NA-values which are not covered by the original covariates, but included in the mesh domain
  # We just fill in the mean value as these will not be used in the modelling anyway (outside weight=0-domain etc.)
  gridcovNew$z[  is.na(gridcovNew$z) ] <- mean(gridcovNew$z,na.rm=TRUE)

  ## Creates an im of gridcocNew
  gridcovNew.im <- spatstat::as.im(gridcovNew)

  return(gridcovNew.im)
}


