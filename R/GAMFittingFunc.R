
#' Do not both to write help function
#'
#' @import mgcv
#' @export

GAMFittingFunc <- function(resultsBaseFolder = "./Results/GAM",
                           sealPhotoDataFile = system.file("extdata", "processed","OriginalSealsKm.rds", package = "SealCoxProcess"),
                           sealTransectDataFile = system.file("extdata", "processed","OigardTablesTransformed.rds", package = "SealCoxProcess"),
                           satelliteDataFolder = system.file("extdata", "processed", package = "SealCoxProcess"),
                           sealType = "hooded",
                           covariate.fitting = "linear",
                           covariates.type = "band1",
                           additional_comment = "",
                           grid.pixelsize = 0.2,
                           leaveOutTransect = FALSE,
                           transectAsCountDomain = FALSE,
                           noSamp=1000,
                           subSampPerSamp = 5,
                           results.CI.level = 0.95,
                           parallelize.numCores = 10,
                           parallelize.noSplits = 2*parallelize.numCores,
                           delete.temp = TRUE,
                           save.data = TRUE,
                           fam="negbin", # "negbin" or "poisson"
                           CVFold = list(0,0,0)) {

  #### Initial definitions ####

  runType <- "GAM" # Setting the type of function

  # Starting time measuring
  time0 <- proc.time()

  # Computing use.covariates variable
  if (covariate.fitting %in% c("linear", "quadratic", "linAndLog", "nonlinear","linearAndSpatial")){
    use.covariates <- T
  } else {
    use.covariates <- F
  }

  inputVar <- formals(GAMFittingFunc)


  # Creating a vector with the input variables
  #inputVar <- formals(GAMFittingFunc)


  #### Initial folder and seed setup ####

  seed <- sum(as.numeric(proc.time()))*1000 # Just to get a random seed (I get the same every time I open R in Thinlinc)
  set.seed(seed) # Only used for creating the foldername

  folder_name0 <- paste("runType=",runType," seal=",sealType," cov=",covariate.fitting," comment=", additional_comment,sep="")
  folder_name <- paste(folder_name0,"_",sample(1:10^3,1),sep="")
  savingFolder <- paste(resultsBaseFolder,folder_name,sep="/") # Where to save all results in the end

  comment=folder_name0 #

  if (!dir.exists(savingFolder)){
    dir.create(savingFolder,recursive = TRUE)
  }

  tempFolder <- paste(savingFolder,"temp",sep="/") # Where to save temporary files
  if (!dir.exists(tempFolder)){
    dir.create(tempFolder,recursive = TRUE)
  }


  #### Loading and preparing the seals data ####

  ## Loading the seal photo data
  if (tools::file_ext(sealPhotoDataFile)=="rds"){
    seals <- readRDS(file=sealPhotoDataFile) # RDS-file
  } else {
    load(file=sealPhotoDataFile) # RData, with a seal object.
  }

  # Creating a list where all important data to be used later are stored explanatory names
  dataList <- list()
  dataList$org <- list()
  dataList$org$noPhoto <- dim(seals)[1]   # Total number of photo taken (with or without seals)
  dataList$org$coordPhotoX <- seals$xkm   # X coordinate of center of each photo, given in kilometers
  dataList$org$coordPhotoY <- seals$ykm   # Y coordinate of center of each photo,  given in kilometers
  dataList$org$photoWidth <- seals$lengthkm # The width (X-direction) of each of the photos, given in kilometers
  dataList$org$photoHeight <- seals$widthkm # The height (Y-direction) of each of the photos, given in kilometers

  ## Number of observed seals of correct type for each photo
  if (sealType=="hooded") dataList$org$noObsPerPhoto <- seals$hooded
  if (sealType=="harps")  dataList$org$noObsPerPhoto <- seals$harps

  ## Loading the transect seal data
  transectData <- readRDS(sealTransectDataFile)

  dataList$transect <- list()
  dataList$transect$noTransects <- nrow(transectData)
  dataList$transect$transectStartCoordX <- transectData$x.start
  dataList$transect$transectStartCoordY <- transectData$y.start
  dataList$transect$transectEndCoordX <- transectData$x.end
  dataList$transect$transectEndCoordY <- transectData$y.end

  ## Loading the satellite covariate data, if applicable

  if (covariates.type%in%c("band1","band2")){
    if (grid.pixelsize<0.1){
      covGrid <- readRDS(file.path(satelliteDataFolder,paste("cov_grid_",covariates.type,"_5000.rds",sep=""))) # A higher resolution covariate image is used to allow smaller pixelsizes
    } else {
      covGrid <- readRDS(file.path(satelliteDataFolder,paste("cov_grid_",covariates.type,".rds",sep="")))
    }
  }

  print("Finsihed loading and preparing seal data")

  #### Creat polygon defining the counting domain for the seals (based on the area spanned by the transects) ####

  # Assigning each photo to a single transect
  transectYmid <- (dataList$transect$transectStartCoordY + dataList$transect$transectEndCoordY)/2

  photoinTransectVec <- rep(NA,dataList$org$noPhoto)
  for (i in 1:dataList$org$noPhoto){
    photoinTransectVec[i] <-which.min(abs(dataList$org$coordPhotoY[i]-transectYmid))
  }

  photoOrderInTransList <- list()
  whichphotoInTransList <- list()
  for (i in 1:dataList$transect$noTransects){
    whichphotoInTransList[[i]] <- which(photoinTransectVec == i)
    photoOrderInTransList[[i]] <- order(dataList$org$coordPhotoX[whichphotoInTransList[[i]]])
  }

  includeThesePhotos <- !(photoinTransectVec%in% leaveOutTransect)

  if(CVFold[[2]]>0){
    n <- nrow(seals)

    if (CVFold[[3]]<1){
      set.seed(123)
      randomOrder <- sample(1:n,n)
      splittedRandomOrder=split(randomOrder,ceiling((1:n)/n*CVFold[[2]]))

      noNeighborPhotos <- -CVFold[[3]]

    } else {
      randomOrder <- rep(NA,n)
      prev <- 0
      for (i in 1:max(photoinTransectVec)){
        these <- prev + 1:sum(photoinTransectVec==i)
        randomOrder[these] <- which(photoinTransectVec==i)
        prev <- max(these)
      }

      splittedRandomOrder=split(randomOrder,sort(photoinTransectVec))
      noNeighborPhotos <- 0
    }

    predThesePhotos <- NULL
    for (i in 1:length(CVFold[[1]])){
      predThesePhotos <- c(predThesePhotos,splittedRandomOrder[[CVFold[[1]][i]]])
    }
    predThesePhotosLogical <- (1:n %in% predThesePhotos)

    photosNotModelled <- NULL

    for (i in 1:length(predThesePhotos)){
      thisTrans <- photoinTransectVec[predThesePhotos[i]]

      thisPhotoInTrans <- which(whichphotoInTransList[[thisTrans]]==predThesePhotos[i])
      thisPhotoInTransOrder <- photoOrderInTransList[[thisTrans]][thisPhotoInTrans]

      neighborPhotos <- unique((thisPhotoInTransOrder-noNeighborPhotos):(thisPhotoInTransOrder+noNeighborPhotos))
      whichneighborPhotos <- which(photoOrderInTransList[[thisTrans]] %in% neighborPhotos)
      newPhotosNotModelled <- whichphotoInTransList[[thisTrans]][whichneighborPhotos]
      photosNotModelled <- c(photosNotModelled,newPhotosNotModelled)
    }

    includeThesePhotos <- !(1:n %in% photosNotModelled)

  }




  ### The data used in the actual modelling
  dataList$mod <- list()
  dataList$mod$coordPhotoX <- dataList$org$coordPhotoX[includeThesePhotos]
  dataList$mod$coordPhotoY <- dataList$org$coordPhotoY[includeThesePhotos]
  dataList$mod$photoWidth <- dataList$org$photoWidth[includeThesePhotos]
  dataList$mod$photoHeight <- dataList$org$photoHeight[includeThesePhotos]
  dataList$mod$noPhoto <- length(dataList$mod$coordPhotoX)
  dataList$mod$noObsPerPhoto <- dataList$org$noObsPerPhoto[includeThesePhotos]


  photos <- list()
  photos$pic.x.left <- dataList$org$coordPhotoX-0.5*dataList$org$photoWidth
  photos$pic.x.right <- dataList$org$coordPhotoX+0.5*dataList$org$photoWidth
  photos$pic.y.bottom <- dataList$org$coordPhotoY-0.5*dataList$org$photoHeight
  photos$pic.y.top <- dataList$org$coordPhotoY+0.5*dataList$org$photoHeight
  orgPhotos <- as.data.frame(photos)
  modPhotos <- orgPhotos[includeThesePhotos,]

  modObservationDomain <- as.data.frame(MergePolygons(modPhotos))

  if (sum(leaveOutTransect)>0 | CVFold[[2]]>0){
    predPhotos <- orgPhotos[predThesePhotosLogical,]
    predDomain <- as.data.frame(MergePolygons(predPhotos))
    if (transectAsCountDomain){
      theseTransInCountDomain <- leaveOutTransect
    } else {
      theseTransInCountDomain <- 1:dataList$transect$noTransects
    }
  } else {
    predDomain <- predPhotos <- NULL
    theseTransInCountDomain <- 1:dataList$transect$noTransects
  }


  countingDomain <- CreateCountDomainPolygon(transectStartCoordX = dataList$transect$transectStartCoordX,
                                             transectStartCoordY = dataList$transect$transectStartCoordY,
                                             transectEndCoordX = dataList$transect$transectEndCoordX,
                                             transectEndCoordY = dataList$transect$transectEndCoordY,
                                             coordPhotoX = dataList$org$coordPhotoX,
                                             coordPhotoY = dataList$org$coordPhotoY,
                                             photoWidth = dataList$org$photoWidth,
                                             photoHeight = dataList$org$photoHeight,
                                             transectYSpan = 1.5*1.852,
                                             theseTransInCountDomain=theseTransInCountDomain)

  areaCountDomain <- rgeos::gArea(coo2sp(countingDomain)) # Should be close to 4569 which Tor-Arne uses in his papers

  print("Finished computation of counting domain")


  logAreakm <- log(dataList$mod$photoHeight*dataList$mod$photoWidth)

  if (use.covariates){
    nearestPixelObsPoints <- spatstat::nearest.pixel(dataList$mod$coordPhotoX, dataList$mod$coordPhotoY,covGrid)
    covAtObsPoints <- covGrid[Reduce('cbind', nearestPixelObsPoints)]
    covAtObsPoints2 <- covAtObsPoints^2
    covAtObsPointsLog <- log(covAtObsPoints)

    fitData <- data.frame(count = dataList$mod$noObsPerPhoto, logAreakm = logAreakm, x = dataList$mod$coordPhotoX, y = dataList$mod$coordPhotoY,
                          covariate = covAtObsPoints, covariate2 = covAtObsPoints2, covariateLog = covAtObsPointsLog,
                          spatialX = dataList$mod$coordPhotoX, spatialY = dataList$mod$coordPhotoY, spatialXY = sqrt(dataList$mod$coordPhotoX^2+dataList$mod$coordPhotoY^2))
    } else {
    fitData <- data.frame(count = dataList$mod$noObsPerPhoto, logAreakm = logAreakm, x = dataList$mod$coordPhotoX, y = dataList$mod$coordPhotoY)
  }



  #### Fitting the GAM model ####

  # Fitting the GAM model
  # negbin familty searching for parameter between (0.1 and 10),
  # Inflate the model degrees of freedom with 1.4 (used in the GCV-method), such smoother fields are favoured
  if (fam=="negbin"){
    if (!use.covariates){
      GAMfit <- gam(count ~ s(x,y)+offset(logAreakm),data = fitData,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
    } else {
      if (covariate.fitting=="linear"){
        GAMfit <- gam(count ~ s(x,y)+covariate+offset(logAreakm),data = fitData,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
      }
      if (covariate.fitting=="quadratic"){
        GAMfit <- gam(count ~ s(x,y)+covariate+covariate2+offset(logAreakm),data = fitData,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
      }
      if (covariate.fitting=="linAndLog"){
        GAMfit <- gam(count ~ s(x,y)+covariateLog+offset(logAreakm),data = fitData,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
      }
      if (covariate.fitting=="nonlinear"){
        GAMfit <- gam(count ~ s(x,y)+s(covariate)+offset(logAreakm),data = fitData,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
      }
      if (covariate.fitting=="linearAndSpatial"){
        #GAMfit <- gam(count ~ s(x,y)+covariate+offset(logAreakm) + spatialX + spatialY + spatialXY,data = fitData,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.42) # For a fit that did not converge with gamma = 1.4
        GAMfit <- gam(count ~ s(x,y)+covariate+offset(logAreakm) + spatialX + spatialY + spatialXY,data = fitData,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
      }

    }
  } else {
    if (!use.covariates){
      GAMfit <- gam(count ~ s(x,y)+offset(logAreakm),data = fitData,family=poisson,optimizer = "perf",gamma = 1.4)
    } else {
      if (covariate.fitting=="linear"){
        GAMfit <- gam(count ~ s(x,y)+covariate+offset(logAreakm),data = fitData,family=poisson,optimizer = "perf",gamma = 1.4)
      }
      if (covariate.fitting=="quadratic"){
        GAMfit <- gam(count ~ s(x,y)+covariate+covariate2+offset(logAreakm),data = fitData,family=poisson,optimizer = "perf",gamma = 1.4)
      }
      if (covariate.fitting=="linAndLog"){
        GAMfit <- gam(count ~ s(x,y)+covariateLog+offset(logAreakm),data = fitData,family=poisson,optimizer = "perf",gamma = 1.4)
      }
      if (covariate.fitting=="nonlinear"){
        GAMfit <- gam(count ~ s(x,y)+s(covariate)+offset(logAreakm),data = fitData,family=poisson,optimizer = "perf",gamma = 1.4)
      }
      if (covariate.fitting=="linearAndSpatial"){
        GAMfit <- gam(count ~ s(x,y)+covariate+offset(logAreakm) + spatialX + spatialY + spatialXY,data = fitData,family=poisson,optimizer = "perf",gamma = 1.4)
      }
     }
  }


  #### Defining a grid ####

  if (grid.pixelsize==0){
    grid.pixelsize = mean(sqrt(dataList$mod$photoWidth*dataList$mod$photoHeight)) # Uses a pixelsize corresponding to the average area of a photo
  }


  gridList <- BasicGridCreation(countingDomain = countingDomain,
                                grid.pixelsize = grid.pixelsize,
                                areaCountDomain = areaCountDomain)


  imputedList <- GAMImpute(gridvalX = gridList$gridvalX,
                           gridvalY = gridList$gridvalY,
                           use.covariates = use.covariates,
                           covGrid = covGrid,
                           GAMfit = GAMfit,
                           logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain)

  resultList <- BasicResultExtractionGAM(GAMfit = GAMfit,
                                         imputedList = imputedList,
                                         use.covariates = use.covariates,
                                         covariate.fitting = covariate.fitting,
                                         logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain,
                                         gridList = gridList)


  if (sum(leaveOutTransect)>0){
    predPhotoGridPointList <- GridPointsInPredPhotosGAM(gridList = gridList,
                                                        predPhotos = predPhotos,
                                                        grid.pixelsize = grid.pixelsize)

    allPhotosinGrid <- 0
    for (i in 1:length(predPhotoGridPointList)){
      allPhotosinGrid <- allPhotosinGrid + as.numeric(predPhotoGridPointList[[i]]$logical)
    }
    allPhotosinGrid <- as.logical(allPhotosinGrid)

  } else {
    predPhotoGridPointList <- NULL
    allPhotosinGrid <- rep(TRUE,gridList$nxy[1]*gridList$nxy[2])
  }

  if (CVFold[[2]]>0){

    predPhotoCenterCord <- cbind(dataList$org$coordPhotoX,dataList$org$coordPhotoY)[predThesePhotosLogical,]

    imputedPredPhotoList <- GAMImputePred(allX = predPhotoCenterCord[,1],
                                          allY = predPhotoCenterCord[,2],
                                          use.covariates = use.covariates,
                                          covGrid = covGrid,
                                          GAMfit = GAMfit)


    Rbeta <- MASS::mvrnorm(n = noSamp, coef(GAMfit), vcov(GAMfit))
    Xp <- predict(GAMfit, newdata = imputedPredPhotoList$imputeData, type = "lpmatrix")
    sampInPhotoGridPoints <- Xp %*% t(Rbeta)

    areaPredPhotos <- abs(predPhotos[,1]-predPhotos[,2])*abs(predPhotos[,3]-predPhotos[,4])


    if (fam=="negbin"){

      areaPerSubSampPredPhoto <- samplePostPredDistGAMPredPhoto(arrayGrid = sampInPhotoGridPoints,
                                                    est.theta = GAMfit$family$Theta,
                                                    subSampPerSamp = subSampPerSamp,
                                                    areaPredPhotos = areaPredPhotos)
    } else {
      areaPerSubSampPredPhoto <- samplePostPredDistGAMPoissonPredPhoto(arrayGrid = sampInPhotoGridPoints,
                                                           subSampPerSamp = subSampPerSamp,
                                                           areaPredPhotos = areaPredPhotos)
    }

    areaPerSubSampPhotoMatrixMEANPredPhoto <- exp(sampInPhotoGridPoints)*areaPredPhotos


    ## Creating also the subsampled area intensity for the union of the photos
    areaPerSubSampTransectVec <- rowSums(areaPerSubSampPredPhoto)


    posthistTransect <- hist(areaPerSubSampTransectVec,breaks=50,plot=F)
    tab.areaPerSubSampTransectVec <- table(areaPerSubSampTransectVec)
    evalFullTransect <- as.numeric(names(tab.areaPerSubSampTransectVec))
    posteriorDistTransect <- as.numeric(tab.areaPerSubSampTransectVec)/sum(as.numeric(tab.areaPerSubSampTransectVec))

    posthistPhotoList <- list()
    evalFullPhotoList <- list()
    posteriorDistPhotoList <- list()
    for (i in 1:ncol(areaPerSubSampPredPhoto)){
      posthistPhotoList[[i]] <- hist(areaPerSubSampPredPhoto[,i],breaks=50,plot=F)
      tab.areaPerSubSampPhotoVec <- table(areaPerSubSampPredPhoto[,i])
      evalFullPhotoList[[i]] <- as.numeric(names(tab.areaPerSubSampPhotoVec))
      posteriorDistPhotoList[[i]] <- as.numeric(tab.areaPerSubSampPhotoVec)/sum(as.numeric(tab.areaPerSubSampPhotoVec))
    }


    postPredDistListPredPhoto <- list()
    postPredDistListPredPhoto$posteriorevalFullPhotoList <- evalFullPhotoList
    postPredDistListPredPhoto$posteriorDistPhotoList <- posteriorDistPhotoList
    postPredDistListPredPhoto$posteriorevalFullTransect <- evalFullTransect
    postPredDistListPredPhoto$posteriorDistTransect <- posteriorDistTransect
    postPredDistListPredPhoto$posthistTransect <- posthistTransect
    postPredDistListPredPhoto$posthistPhotoList <- posthistPhotoList
    postPredDistListPredPhoto$areaPerSubSampPhotoMatrix <- areaPerSubSampPhotoMatrixMEANPredPhoto # Note that this is not the same as areaPerSubSampPhotoMatrix used internally in this function... (but the same as used for the INLA approach)

  }




  #### Doing the sampling ####

  samp <- postPredDistGAMFunc(GAMfit=GAMfit,
                              noSamp = noSamp,
                              imputedList = imputedList,
                              allPhotosinGrid = allPhotosinGrid)

  if (sum(leaveOutTransect)>0){

    postPredDistList <- PhotoPostPredDistGAM(samp = samp,
                                             parallelize.noSplits = parallelize.noSplits,
                                             parallelize.numCores = parallelize.numCores,
                                             tempFolder = tempFolder,
                                             gridList = gridList,
                                             predPhotoGridPointList = predPhotoGridPointList,
                                             est.theta = GAMfit$family$Theta,
                                             allPhotosinGrid = allPhotosinGrid,
                                             subSampPerSamp = subSampPerSamp,
                                             fam = fam)
  } else {

    postPredDistList <- FullPostPredDistGAM(samp = samp,
                                            parallelize.noSplits = parallelize.noSplits,
                                            parallelize.numCores = parallelize.numCores,
                                            tempFolder = tempFolder,
                                            gridList = gridList,
                                            predPhotoGridPointList = predPhotoGridPointList,
                                            est.theta = GAMfit$family$Theta,
                                            subSampPerSamp = subSampPerSamp,
                                            fam = fam)
  }

  if (CVFold[[2]]>0){

    photoOrder <- 1:nrow(predPhotos)#order(dataList$org$coordPhotoX[!includeThesePhotos])

    ComparePhotoCountsAndPredGAM(posteriorevalFullPhoto = postPredDistListPredPhoto$posteriorevalFullPhotoList,
                                 posteriorDistPhoto = postPredDistListPredPhoto$posteriorDistPhotoList,
                                 posteriorevalFullTransect = postPredDistListPredPhoto$posteriorevalFullTransect,
                                 posteriorDistTransect = postPredDistListPredPhoto$posteriorDistTransect,
                                 areaPerSubSampPhotoMatrix = t(postPredDistListPredPhoto$areaPerSubSampPhotoMatrix),
                                 photoCounts = dataList$org$noObsPerPhoto[predThesePhotosLogical],
                                 photoOrder = photoOrder,
                                 leaveOutTransect = CVFold[[1]],
                                 inputVar = inputVar,
                                 savingFolder = savingFolder)
  }



  if (sum(leaveOutTransect)>0){

    photoOrder <- order(dataList$org$coordPhotoX[predThesePhotosLogical])

    ComparePhotoCountsAndPredGAM(posteriorevalFullPhoto = postPredDistList$posteriorevalFullPhotoList,
                              posteriorDistPhoto = postPredDistList$posteriorDistPhotoList,
                              posteriorevalFullTransect = postPredDistList$posteriorevalFullTransect,
                              posteriorDistTransect = postPredDistList$posteriorDistTransect,
                              areaPerSubSampPhotoMatrix = postPredDistList$areaPerSubSampPhotoMatrix,
                              photoCounts = dataList$org$noObsPerPhoto[predThesePhotosLogical],
                              photoOrder = photoOrder,
                              leaveOutTransect = leaveOutTransect,
                              inputVar = inputVar,
                              savingFolder = savingFolder)
  }


  if (sum(leaveOutTransect)>0){
    postResList <- SummaryStat(evalPoints = postPredDistList$posteriorevalFullTransect,
                               dist = postPredDistList$posteriorDistTransect,
                               results.CI.level = results.CI.level,
                               posterior = TRUE)
  } else {
    postResList <- SummaryStat(evalPoints = postPredDistList$posteriorEvalPoints,
                               dist = postPredDistList$posteriorDist,
                               results.CI.level = results.CI.level,
                               posterior = TRUE)
  }


  finalResList <- c(resultList,postPredDistList,postResList)

  timeUsed <- proc.time()-time0





  SummaryPlotFuncGAM(covariatesplot = use.covariates,
                     summaryplot = TRUE,
                     savingFolder = savingFolder,
                     sealPhotoDataFile = sealPhotoDataFile,
                     sealTransectDataFile = sealTransectDataFile,
                     dataList = dataList,
                     orgPhotos = orgPhotos,
                     modPhotos = modPhotos,
                     results.CI.level = results.CI.level,
                     gridList = gridList,
                     finalResList = finalResList,
                     countingDomain = countingDomain,
                     logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain,
                     covNewGridval = imputedList$imputeData$covariate,
                     GAMfit = GAMfit,
                     sealType = sealType,
                     use.covariates = use.covariates,
                     covariates.type = covariates.type,
                     covariate.fitting = covariate.fitting,
                     grid.pixelsize =grid.pixelsize,
                     parallelize.noSplits  =parallelize.noSplits,
                     parallelize.numCores = parallelize.numCores,
                     noSamp = noSamp,
                     subSampPerSamp = subSampPerSamp,
                     time = timeUsed,
                     comment = comment,
                     leaveOutTransect = leaveOutTransect,
                     fam = fam)

  #### Saves RData to file ####

  if (save.data){
    save.list <- c("dataList","finalResList","savingFolder","inputVar","gridList","imputedList","timeUsed","countingDomain","areaCountDomain","GAMfit")
    save(list=save.list,file=file.path(savingFolder,"output.RData"))
  }

  if (delete.temp){
    unlink(tempFolder,recursive=TRUE)
  }


  #### End of function ####

  cat(paste("Finished running the GAMfitting procedure function. Output found in \n\n",savingFolder,sep=""))
}




