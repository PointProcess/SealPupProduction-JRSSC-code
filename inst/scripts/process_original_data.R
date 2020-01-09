
library(SealCoxProcess)

library(sp)
library(rgdal)
library(spatstat)
library(fields)
library(akima)

processed_folder <- here::here("inst","extdata", "processed")
original_folder <- here::here("inst","extdata", "original")


#### Processing orignal data ####

# Convert original seals data to km data
infile <- "WestIce2012.csv"
outfile <- "OriginalSealsKm.rds"
seals=read.csv(file=file.path(original_folder,infile))


## Convert coordinates to km (from Nm and m)
seals$xkm <- seals$x*1852/1000 # from Nm
seals$ykm <- seals$y*1852/1000 # from Nm
seals$lengthkm <- seals$length/1000 # from m
seals$widthkm <- seals$width/1000 # from m

saveRDS(seals,file=outfile)

###### Creating transect data counts on same x-y-km-format #####

infile <- "OigardTransectCounts.rds"
outfile <- "OigardTablesTransformed.rds"

tab1Df <- readRDS(file.path(processed_folder,infile))

meanHeightPic=mean(seals$width)
meanWidthPic=mean(seals$length)
scalefactor <- (3*1852/1000)/meanHeightPic

lon0=mean(seals$lon)
lat0=mean(seals$lat)

startCoord <- lb2xykm(tab1Df$lonStart,tab1Df$lat,lon0,lat0)
endCoord <- lb2xykm(tab1Df$lonEnd,tab1Df$lat,lon0,lat0)
transectLength <- rowSums((startCoord-endCoord)^2)
transectArea <- meanHeightPic*transectLength
noHarpsAreaCovered <- tab1Df$noHarps*scalefactor
noHoodedAreaCovered <- tab1Df$noHooded*scalefactor

retDf <- data.frame(x.start = pmin(startCoord$x,endCoord$x), x.end = pmax(startCoord$x,endCoord$x),
                    y.start = pmin(startCoord$y,endCoord$y), y.end = pmax(startCoord$y,endCoord$y),
                    noHarps=tab1Df$noHarps,noHooded=tab1Df$noHooded,noHarpsAreaCovered,noHoodedAreaCovered)
retDf$y.start.area = retDf$y.start-1.5*1.852
retDf$y.end.area = retDf$y.end+1.5*1.852
retDf$noPhotos <- tab1Df$noPhotos

# Saving these
saveRDS(retDf,file=file.path(processed_folder,outfile))

###### Creating the covariates ##########3

outfolder <- here::here("inst","extdata", "processed")

band1_2000 <- CovariateCreationFunction(noPixEachDim = 2000,covariates.type="band1")
saveRDS(object = band1_2000, file = file.path(outfolder,"cov_grid_band1.rds"))

band2_2000 <- CovariateCreationFunction(noPixEachDim = 2000,covariates.type="band2")
saveRDS(object = band2_2000, file = file.path(outfolder,"cov_grid_band2.rds"))


