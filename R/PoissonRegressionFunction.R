
#' Seal counting using INLA and Poisson regression formulation
#'
#' Approximates a Log-Gaussian Cox process with a regression formulation and fits seal counts from transects to it.
#' Samples from the the fitted model to compute the posterior predictive distribution of the total number of seals within an
#' area which the transects are spanning
#'
#' @param resultsBaseFolder String, indicating the path to the base folder where all results are to be stored
#' @param sealPhotoDataFile String, indicating the file where seal photo data are stored
#' @param sealTransectDataFile String, indicating the file where seal transect data are stored
#' @param satelliteDataFolder String, indicating the path to the folder where the satellite data are stored (in files on the form "cov_grid_[covariates.type].rds")
#' @param sealType String, indicating which seal type to model "hooded" (default as it is quicker), or "harps"
#' @param spatial Logical indicating whether a spatial spde model should be used
#' @param covariate.fitting String, indicating how to model covariates. "no", "linear", "quadratic" (default), "linAndLog" or "nonlinear"
#' @param covariates.type String equal to "band1" or "band2" indicating which of the bands from the satellite data should be used.
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param noSamp Numeric, indicating the number of samples from posterior model to use when computing the posterior predictive distribution. 5000 (default) typically sufficient.
#' @param standardMesh Logical indicating whether a standard mesh only fixing mesh points at the photo centers should be created. Defaults to FALSE indicating that the special mesh with Voronoi tesselations corresponding to the photos should be used.
#' @param convHullVar.convex Numeric, corrresponding to the convex parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.concave Numeric, corrresponding to the concave parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.resolution Numeric vector of length 2, corrresponding to the resolution parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param meshVar.max.edge Numeric, corrresponding to the max.edge parameter of inla.mesh.2d. Smaller values gives smaller triangles outside triangles. See ?inla.mesh.2d for details
#' @param meshVar.offset Numeric, corrresponding to the offset parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param meshVar.cutoff Numeric vector, where the indicies corrresponds to the cutoff parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param y.cutoff.boundary Numeric vector deciding the y-values which should be the coundary for the  different cutoff-values. NULL means the first is applied everywhere
#' @param Matern.alpha Numeric, corresponding to the alpha parameter in the Matern covariance function (2 is the default)
#' @param grid.pixelsize Numeric, denoting the size (in km, in both x- and y-direction) of each pixel of the grid being used
#' @param INLA.theta.startval List containing the start values for the theta parameter in the INLA run. (NULL indicates that automatic start values should be used, and is the default)
#' @param INLA.verbose Logical, indicating whether verbose printing should be used when running the actual inla function
#' @param parallelize.numCores Numeric, corresponding to the number of cores any parallelization should be run at
#' @param parallelize.noSplits Numeric, deciding how many sublists samp should be splitted into. Should be a multiple of parallelize.numCores for the highest efficiency. The larger number the less memory is used (and longer time).
#' @param poisson.maxEvals Numeric, corresponding to maximum number of points the Poisson distribution should be evaluated at (a much smaller number is typically used)
#' @param results.CI.level Numeric, denoting the confidence/credibility degree to use in the final credibility interval for the total number of counts
#' @param additional_comment String, where any comments related to the actual run can be given
#' @param save.data Logical, indicating whether input variables, data and results should be saved (TRUE is default)
#' @param delete.temp Logical, indicating whether temporary stored samples and associated eta grids should be deleted
#' @param testing Logical, indicating whether the testing parts of this function should be used
#' @param leaveOutTransect The transect number for the transect to remove from the modelling fitting and then fit (defaults to FALSE, i.e. none)
#' @param transectAsCountDomain Logical indicating whether the transects should be used as counting domain when they are left out of the modeling (does not apply if leaveOutTransect=F)
#' @param CVFold List of size 3. The is a vector indicating the fold number(s) to remove; the second a single number with the total number of folds; the third is a single number indicates whether it should be a fixed (1) "sample" based on the transect number,
#' or random (0) sample. Defaults to c(0,0,0) which means no cross validation is performed
#' @param INLA.constr Logical deciding whether the SPDE model should be fitted with the constaint that it integrates to zero or not.
#' @return Nothing really, but saves results in the resultsBaseFolder
#' @keywords inla cox-process Poisson-regression
#' @import INLA
#' @import spatstat
#' @import fields
#' @import sp
#' @import parallel
#' @import rgeos
#' @import doParallel
#' @import Hmisc
#' @import tools
#' @export


INLAPPSealsPoissonReg <- function(resultsBaseFolder = "./Results/PoissonReg",
                                  sealPhotoDataFile = system.file("extdata", "processed","OriginalSealsKm.rds", package = "SealCoxProcess"),
                                  sealTransectDataFile = system.file("extdata", "processed","OigardTablesTransformed.rds", package = "SealCoxProcess"),
                                  satelliteDataFolder = system.file("extdata", "processed", package = "SealCoxProcess"),
                                  sealType = "hooded",
                                  spatial = TRUE,
                                  covariate.fitting = "linear",
                                  covariates.type = "band1",
                                  additional.iid.term = FALSE,
                                  noSamp = 5000,
                                  standardMesh = FALSE,
                                  convHullVar.convex = -0.15,
                                  convHullVar.concave = convHullVar.convex,
                                  convHullVar.resolution = c(120,120),
                                  meshVar.max.edge = c(2,10),
                                  meshVar.offset = 6,
                                  meshVar.cutoff = c(0.195,0.23),
                                  y.cutoff.boundary=-60,
                                  Matern.alpha = 2,
                                  grid.pixelsize = 0.2,
                                  INLA.theta.startval = NULL,
                                  INLA.verbose = FALSE,
                                  parallelize.numCores = 10,
                                  parallelize.noSplits = parallelize.numCores,
                                  poisson.maxEvals = 5*10^5,
                                  results.CI.level = 0.95,
                                  additional_comment = "",
                                  save.data = TRUE,
                                  delete.temp = TRUE,
                                  testing = FALSE,
                                  leaveOutTransect = FALSE,
                                  transectAsCountDomain = TRUE,
                                  CVFold = list(0,0,0), # Should use transectAsCountDomain=FALSE when this is used!
                                  INLA.constr = TRUE
){

  #### Initial definitions ####
  thisEnvir <- parent.frame()

  runType <- "PoissonReg" # Setting the type of function

  # Starting time measuring
  time0 <- proc.time()

  # Computing use.covariates variable
  if (covariate.fitting %in% c("linear", "quadratic", "linAndLog", "nonlinear","linearAndSpatial")){
    use.covariates <- T
  } else {
    use.covariates <- F
  }

  # Creating a vector with the input variables
  inputVar <- formals(INLAPPSealsPoissonReg)

  # A list with the input variables to be used directly
  #inputList <- list()
  #for (i in 1:length(inputVar)){
  #  eval(parse(text=paste("inputList$",inputVar[i]," <- ",eval(parse(text=inputVar[i])),sep="")))
  #}


  #### Initial folder and seed setup ####

  seed <- sum(as.numeric(proc.time()))*1000 # Just to get a random seed (I get the same every time I open R in Thinlinc)
  set.seed(seed) # Only used for creating the foldername

  folder_name0 <- paste("runType=",runType," seal=",sealType," cov=",covariate.fitting," extra.iid=", additional.iid.term ," samp=",noSamp, " comment=", additional_comment,sep="")
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


  #### Create polygon defining the counting domain for the seals (based on the area spanned by the transects) ####




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

  countingDomainExpanded <- CreateCountDomainPolygon(transectStartCoordX = dataList$transect$transectStartCoordX,
                                                     transectStartCoordY = dataList$transect$transectStartCoordY,
                                                     transectEndCoordX = dataList$transect$transectEndCoordX,
                                                     transectEndCoordY = dataList$transect$transectEndCoordY,
                                                     coordPhotoX = dataList$org$coordPhotoX,
                                                     coordPhotoY = dataList$org$coordPhotoY,
                                                     photoWidth = dataList$org$photoWidth,
                                                     photoHeight = dataList$org$photoHeight,
                                                     transectYSpan = 2*1.852,
                                                     transectXSpan = 0.5*1.852,
                                                     theseTransInCountDomain=theseTransInCountDomain)

  areaCountDomain <- rgeos::gArea(coo2sp(countingDomain)) # Should be close to 4569 which Tor-Arne uses in his papers

  print("Finished computation of counting domain")


  #### Cunstructing the mesh ####
  if (!testing){
    if (standardMesh){
      rectangleCentersX = dataList$org$coordPhotoX
      rectangleCentersY = dataList$org$coordPhotoY
      rectangleWidth = dataList$org$photoWidth
      rectangleHeight = dataList$org$photoHeight

      domain.outer  <- inla.nonconvex.hull(points=cbind(rectangleCentersX,rectangleCentersY),
                                           convex=convHullVar.convex+0.03,
                                           concave=convHullVar.concave+0.03,
                                           resolution=convHullVar.resolution)

      domain.outer.final  <- inla.nonconvex.hull(points=cbind(rectangleCentersX,rectangleCentersY),
                                                 convex=convHullVar.convex,
                                                 concave=convHullVar.concave,
                                                 resolution=convHullVar.resolution)
      domain.inner <- INLA::inla.mesh.segment(loc = as.matrix(countingDomainExpanded))

      obligMeshLoc <- as.data.frame(cbind(x=rectangleCentersX,y=rectangleCentersY))

      mesh <- inla.mesh.2d(loc=obligMeshLoc,
                           boundary = list(domain.inner,domain.outer),
                           max.edge=meshVar.max.edge, # Try c(1.5,10)
                           offset=meshVar.offset,
                           cutoff=meshVar.cutoff[1]) # Try 0.17


      predPhotoMeshPointList <- MeshPointsInPredPhotos(mesh = mesh,
                                                       predPhotos = orgPhotos)
      noMeshPoints = rep(NA,length(predPhotoMeshPointList))
      for (i in 1:length(predPhotoMeshPointList)){
        noMeshPoints[i] <- predPhotoMeshPointList[[i]]$noMeshPoints
      }
      print(paste("ALL prediction photos matching exactly 1 mesh point? ",all.equal(noMeshPoints,rep(1,length(predPhotoMeshPointList))),sep=""))

      thisMeshPoint <- rep(NA,length(predPhotoMeshPointList))
      for (i in 1:length(predPhotoMeshPointList)){
        thisMeshPoint[i] <- which(predPhotoMeshPointList[[i]]$logical)
      }


      # par(mfrow=c(1,2))
      # plot(mesh1,xlim=-15+c(-5,5),ylim=-15+c(-5,5))
      # rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=2)
      # points(mesh0$loc[2013,1],mesh0$loc[2013,2])
      #
      # plot(mesh,xlim=-15+c(-5,5),ylim=-15+c(-5,5))
      #
      # rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=2)
      #
      #
      #
      # par(mfrow=c(1,2))
      # plot(mesh0,xlim=c(-40,-28),ylim=c(-100,-90)+20)
      # rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=2)
      #
      # plot(mesh,xlim=c(-40,-28),ylim=c(-100,-90)+20)
      #
      # rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=2)


    } else {
      mesh <- MeshCreationMatchingRectangles(rectangleCentersX = dataList$org$coordPhotoX,
                                             rectangleCentersY = dataList$org$coordPhotoY,
                                             rectangleWidth = dataList$org$photoWidth,
                                             rectangleHeight = dataList$org$photoHeight,
                                             convHullVar.convex = convHullVar.convex,
                                             convHullVar.concave = convHullVar.concave,
                                             convHullVar.resolution = convHullVar.resolution,
                                             meshVar.max.edge = meshVar.max.edge,
                                             meshVar.offset = meshVar.offset,
                                             meshVar.cutoff = meshVar.cutoff,
                                             y.cutoff.boundary = y.cutoff.boundary,
                                             countingDomainExpanded = countingDomainExpanded)
    }
  }
  if (testing){




    mesh <- MeshCreationMatchingRectangles(rectangleCentersX = dataList$mod$coordPhotoX,
                                           rectangleCentersY = dataList$mod$coordPhotoY,
                                           rectangleWidth = dataList$mod$photoWidth,
                                           rectangleHeight = dataList$mod$photoHeight,
                                           convHullVar.convex = convHullVar.convex,
                                           convHullVar.concave = convHullVar.concave,
                                           convHullVar.resolution = convHullVar.resolution,
                                           meshVar.max.edge = meshVar.max.edge,
                                           meshVar.offset = meshVar.offset,
                                           meshVar.cutoff = meshVar.cutoff,
                                           y.cutoff.boundary = y.cutoff.boundary,
                                           countingDomainExpanded = countingDomainExpanded)
  }


  if (covariate.fitting=="nonlinear"){
    noMeshpoints.nonlinear = 30
    degree.nonlinear = 2

    rangeCovGrid <- range(c(covGrid$v))
    covMesh <- inla.mesh.1d(loc = seq(rangeCovGrid[1],rangeCovGrid[2],length.out=noMeshpoints.nonlinear),degree = degree.nonlinear)
  } else {
    covMesh <- NULL
  }

  print("Finished constructing the mesh")




  #### Creating the Voronoi tessellations and computing the area each of them cover within the  and specfying the weights e.pp for each of the polygons (based on its size) ####

  if (standardMesh){
   # meshweight <- rep(0,mesh$n)
  #  for (i in 1:length(thisMeshPoint)){
  #    meshweight[i] <- dataList$org$photoWidth[i]*dataList$org$photoHeight[i]
  #  }

    modPhotosAreTheseMeshPoints <- thisMeshPoint[includeThesePhotos]

    unobservedMeshPoints <- which(!(1:mesh$n %in% modPhotosAreTheseMeshPoints))

    dataList$mod$NAMeshLoc <- unobservedMeshPoints # Unobserved mesh locations
    dataList$mod$y.pp <- rep(NA,mesh$n)
    dataList$mod$e.pp <- rep(NA,mesh$n)
    dataList$mod$y.pp[modPhotosAreTheseMeshPoints] <- dataList$mod$noObsPerPhoto
    dataList$mod$e.pp[modPhotosAreTheseMeshPoints] <- dataList$mod$photoWidth*dataList$mod$photoHeight

    dataList$mod$y.pp <- dataList$mod$noObsPerPhoto
    dataList$mod$e.pp <- dataList$mod$photoWidth*dataList$mod$photoHeight
    dataList$mod$obsLoc <- mesh$loc[modPhotosAreTheseMeshPoints,]

    # Just need to define these
    voronoiTess <- list()
    weightAtMeshLoc <- NULL


  } else {
    voronoiTess <- CreateVoronoiTessellation(locationsCoordX = mesh$loc[,1],
                                             locationsCoordY = mesh$loc[,2],
                                             observationDomain = modObservationDomain,
                                             parallelize.numCores = parallelize.numCores) ### Might consider replacing countingDomian with a somewhat larger modelling domain here instead

    weightAtMeshLoc <- voronoiTess$tileSize

    print("Finished computation of the Voronoi tesselation")

    ## Checking which tile the photo centers belong to and assign their values to them.

    yAtMeshLoc <- rep(0,mesh$n)
    eAtMeshLoc <- rep(0,mesh$n)

    for (i in 1:mesh$n){
      insidePhotos <- which(as.logical(point.in.polygon(dataList$mod$coordPhotoX,dataList$mod$coordPhotoY,
                                                        voronoiTess$tiles[[i]]$x,voronoiTess$tiles[[i]]$y)))

      if (length(insidePhotos)==0){
        yAtMeshLoc[i] <- NA
        eAtMeshLoc[i] <- NA
      } else {
        yAtMeshLoc[i] <- sum(dataList$mod$noObsPerPhoto[insidePhotos])
        eAtMeshLoc[i] <- weightAtMeshLoc[i]
      }
    }

    ## Exclude all Voronoi tesselations where there are no photo centers inside, as these are actually unobserved
    dataList$mod$NAMeshLoc <- which(is.na(yAtMeshLoc)) # Unobserved mesh locations
    dataList$mod$y.pp <- yAtMeshLoc[-dataList$mod$NAMeshLoc]
    dataList$mod$e.pp <- eAtMeshLoc[-dataList$mod$NAMeshLoc]
    dataList$mod$obsLoc <- mesh$loc[-dataList$mod$NAMeshLoc,]

    # ### #Testing
    # ## Exclude all Voronoi tesselations where there are no photo centers inside, as these are actually unobserved
    # dataList$mod2=list()
    # dataList$mod2$NAMeshLoc <- which(is.na(yAtMeshLoc)) # Unobserved mesh locations
    # dataList$mod2$y.pp <- yAtMeshLoc[-dataList$mod$NAMeshLoc]
    # dataList$mod2$e.pp <- eAtMeshLoc[-dataList$mod$NAMeshLoc]
    # dataList$mod2$obsLoc <- meshold$loc[-dataList$mod$NAMeshLoc,]
    #
    # plot(dataList$mod2$y.pp,dataList$mod2$e.pp)
    # points(dataList$mod$y.pp,dataList$mod$e.pp,col=2)
    #
    # plot(dataList$mod2$y.pp,dataList$mod$y.pp)
    #
    # plot(dataList$mod$obsLoc)
    #
    #     plot(dataList$mod$obsLoc)
    # points(dataList$mod$obsLoc[dataList$mod$e.pp<0.07,1:2],col=2)
    # points(dataList$mod$obsLoc[dataList$mod$y.pp>2,1:2],col=3)
    #
    #
    # plot(dataList$mod2$obsLoc)
    #
    #
    # #####

    }


  #### Preparing and executing the model fitting ####


  inlaPrepList <- PrepareINLAFunc(mesh = mesh,
                                  covMesh = covMesh,
                                  obsLoc = dataList$mod$obsLoc,
                                  y.pp = dataList$mod$y.pp,
                                  e.pp = dataList$mod$e.pp,
                                  covGrid = covGrid,
                                  spatial = spatial,
                                  covariate.fitting = covariate.fitting,
                                  additional.iid.term = additional.iid.term,
                                  Matern.alpha = Matern.alpha,
                                  covariates.type = covariates.type,
                                  INLA.theta.startval = INLA.theta.startval,
                                  INLA.constr = INLA.constr)

  if(spatial){
    spde <- inlaPrepList$spde
  } else {
    spde <- NULL
  }
  if(covariate.fitting=="nonlinear"){
    covSpde <- inlaPrepList$covSpde
  } else {
    covSpde <- NULL
  }


  print("Finished building INLA stack and preparing the INLA run")

  #### Running the INLA function ####

  ss <- proc.time()

 # if(!testing){
    pp.res <- inla(inlaPrepList$formula,
                   family="poisson", data=inla.stack.data(inlaPrepList$stk.pp),
                   control.predictor=list(A=inla.stack.A(inlaPrepList$stk.pp)),
                   E=inla.stack.data(inlaPrepList$stk.pp)$e,
                   verbose=INLA.verbose,
                   control.compute=inlaPrepList$control.compute.list,
                   #control.inla=list(int.strategy='eb'), # just for speed-up. Gives close to identical results, BUT ASSUMES THETA FIXED WHEN SAMPLING
                   control.mode = inlaPrepList$control.mode.list)
#  }

  # if (testing){
  #
  #   control.fixed = list(mean.intercept=-1.6411,prec.intercept = 1/0.1915^2,
  #                        mean = 9.6770, prec = 1/0.9942^2)
  #
  #   pp.res <- inla(inlaPrepList$formula,
  #                  family="poisson", data=inla.stack.data(inlaPrepList$stk.pp),
  #                  control.predictor=list(A=inla.stack.A(inlaPrepList$stk.pp)),
  #                  E=inla.stack.data(inlaPrepList$stk.pp)$e,
  #                  verbose=INLA.verbose,
  #                  control.compute=inlaPrepList$control.compute.list,
  #                  #control.inla=list(int.strategy='eb'), # just for speed-up. Gives close to identical results, BUT ASSUMES THETA FIXED WHEN SAMPLING
  #                  control.mode = inlaPrepList$control.mode.list,
  #                  control.fixed = control.fixed)
  # }


#    control.mode.list <- list(theta=pp.res$mode$theta,
#                              x=pp.res$mode$x,
#                              restart = TRUE,
#                              fixed = TRUE)

#
#
#     pp.res.theta.fixed <- inla(inlaPrepList$formula,
#                         family="poisson", data=inla.stack.data(inlaPrepList$stk.pp),
#                         control.predictor=list(A=inla.stack.A(inlaPrepList$stk.pp)),
#                         E=inla.stack.data(inlaPrepList$stk.pp)$e,
#                         verbose=INLA.verbose,
#                         control.compute=inlaPrepList$control.compute.list,
#                         control.mode = inlaPrepList$control.mode.list,
#                         control.fixed = control.fixed)
#                         #control.inla=list(int.strategy='eb'))
#
#
#     # pp.res.new2 <- inla(inlaPrepList$formula,
    #                     family="poisson", data=inla.stack.data(inlaPrepList$stk.pp),
    #                     control.predictor=list(A=inla.stack.A(inlaPrepList$stk.pp)),
    #                     E=inla.stack.data(inlaPrepList$stk.pp)$e,
    #                     verbose=INLA.verbose,
    #                     control.compute=inlaPrepList$control.compute.list,
    #                     control.mode = control.mode.list,
    #                     control.inla=list(int.strategy='eb'))
    #
    #
    # control.mode.list <- list(theta=pp.res$mode$theta,
    #                           x=pp.res$mode$x,
    #                           restart = TRUE,
    #                           fixed = FALSE)
    #
    # pp.res.new2 <- inla(inlaPrepList$formula,
    #                     family="poisson", data=inla.stack.data(inlaPrepList$stk.pp),
    #                     control.predictor=list(A=inla.stack.A(inlaPrepList$stk.pp)),
    #                     E=inla.stack.data(inlaPrepList$stk.pp)$e,
    #                     verbose=INLA.verbose,
    #                     control.compute=inlaPrepList$control.compute.list,
    #                     control.mode = control.mode.list,
    #                     control.inla=list(int.strategy='eb'))
    #
    #
    # control.mode.list <- list(theta=pp.res.new2$mode$theta,
    #                           x=pp.res.new2$mode$x,
    #                           restart = TRUE,
    #                           fixed = FALSE)
    #
    # pp.res.new3 <- inla(inlaPrepList$formula,
    #                     family="poisson", data=inla.stack.data(inlaPrepList$stk.pp),
    #                     control.predictor=list(A=inla.stack.A(inlaPrepList$stk.pp)),
    #                     E=inla.stack.data(inlaPrepList$stk.pp)$e,
    #                     verbose=INLA.verbose,
    #                     control.compute=inlaPrepList$control.compute.list,
    #                     control.mode = control.mode.list,
    #                     control.inla=list(int.strategy='eb'))
    #
    # control.mode.list <- list(theta=pp.res.new3$mode$theta+1,
    #                           x=pp.res.new3$mode$x,
    #                           restart = TRUE,
    #                           fixed = FALSE)
    #
    # pp.res.new4 <- inla(inlaPrepList$formula,
    #                     family="poisson", data=inla.stack.data(inlaPrepList$stk.pp),
    #                     control.predictor=list(A=inla.stack.A(inlaPrepList$stk.pp)),
    #                     E=inla.stack.data(inlaPrepList$stk.pp)$e,
    #                     verbose=INLA.verbose,
    #                     control.compute=inlaPrepList$control.compute.list,
    #                     control.mode = control.mode.list,
    #                     control.inla=list(int.strategy='eb'))
    #
    #
    #
    #
    # control.fixed = list(mean.intercept=-1.579697,prec.intercept = 1/0.04683453,
    #                      mean = 9.609047, prec = 1/1.449894)
    #
    # totAreaObs <- sum(abs(orgPhotos$pic.x.left-orgPhotos$pic.x.right)*abs(orgPhotos$pic.y.top-orgPhotos$pic.y.bottom))
    # totSealCount <- seals$hooded
    # noSealsPerKmsq <- sum(totSealCount)/sum(totAreaObs)
    # meanFieldVal <- log(noSealsPerKmsq)
    #
    # fullcountingDomain <- CreateCountDomainPolygon(transectStartCoordX = dataList$transect$transectStartCoordX,
    #                                            transectStartCoordY = dataList$transect$transectStartCoordY,
    #                                            transectEndCoordX = dataList$transect$transectEndCoordX,
    #                                            transectEndCoordY = dataList$transect$transectEndCoordY,
    #                                            coordPhotoX = dataList$org$coordPhotoX,
    #                                            coordPhotoY = dataList$org$coordPhotoY,
    #                                            photoWidth = dataList$org$photoWidth,
    #                                            photoHeight = dataList$org$photoHeight,
    #                                            transectYSpan = 1.5*1.852,
    #                                            theseTransInCountDomain=1:27)
    #
    # areaFullCountDomain <- rgeos::gArea(coo2sp(fullcountingDomain)) # Should be close to 4569 which Tor-Arne uses in his papers
    # homogenExpectedCount <- noSealsPerKmsq*areaFullCountDomain
    #
    #
    # pp.res.new2 <- inla(inlaPrepList$formula,
    #                     family="poisson", data=inla.stack.data(inlaPrepList$stk.pp),
    #                     control.predictor=list(A=inla.stack.A(inlaPrepList$stk.pp)),
    #                     E=inla.stack.data(inlaPrepList$stk.pp)$e,
    #                     verbose=INLA.verbose,
    #                     control.compute=inlaPrepList$control.compute.list,
    #                     control.mode = inlaPrepList$control.mode.list,
    #                     control.inla=list(int.strategy='eb'),
    #                     control.fixed = control.fixed)

  #pp.res = pp.res.theta.fixed
 # }

  runINLATime <- proc.time()-ss

  print(paste("Finished running the INLA-function in ",round(runINLATime[3]/60)," minutes",sep=""))

  #### Extracting results from the INLA model, plots and estimates the range parameter ####

  if (transectAsCountDomain){
    gridSpan <-  "countingDomain"
  } else {
    gridSpan <- "mesh"
  }


  gridList <- GridCreation(mesh = mesh,
                           countingDomain = countingDomain,
                           areaCountDomain = areaCountDomain,
                           grid.pixelsize = grid.pixelsize,
                           gridSpan = gridSpan)


  if (sum(leaveOutTransect)>0){
    predPhotoGridPointList <- GridPointsInPredPhotos(gridList = gridList,
                                                     predPhotos = predPhotos)
  } else {
    predPhotoGridPointList <- NULL
  }

  if (CVFold[[2]]>0){
    predPhotoMeshPointList <- MeshPointsInPredPhotos(mesh = mesh,
                                                     predPhotos = predPhotos)
    noMeshPoints = rep(NA,length(predPhotoMeshPointList))
    for (i in 1:length(predPhotoMeshPointList)){
      noMeshPoints[i] <- predPhotoMeshPointList[[i]]$noMeshPoints
    }
    print(paste("ALL prediction photos matching exactly 1 mesh point? ",all.equal(noMeshPoints,rep(1,length(predPhotoMeshPointList))),sep=""))

    thisMeshPoint <- rep(NA,length(predPhotoMeshPointList))
    for (i in 1:length(predPhotoMeshPointList)){
      thisMeshPoint[i] <- which(predPhotoMeshPointList[[i]]$logical)
    }

    areaPredPhotos <- abs(predPhotos[,1]-predPhotos[,2])*abs(predPhotos[,3]-predPhotos[,4])
    nearestPixelPredLoc <- spatstat::nearest.pixel(mesh$loc[thisMeshPoint,1], mesh$loc[thisMeshPoint,2],covGrid)
    covAtPredLoc <- covGrid[Reduce('cbind', nearestPixelPredLoc)]
    covAtPredLoc2 <- covAtPredLoc^2
    covAtPredLocLog <- log(covAtPredLoc)

    }



  covGridList <- covAtNewGrid(use.covariates = use.covariates,
                              covGrid = covGrid,
                              gridvalX = gridList$gridvalX,
                              gridvalY = gridList$gridvalY,
                              modelledGridVal = gridList$modelledGridVal,
                              logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain)



  if (covariate.fitting=="nonlinear"){
    extraNonlinear.covGridList <- covGridCreationNonlinear(covMesh = covMesh,
                                                           covariateValues = covGridList$covariateValues)
  }  else {
    extraNonlinear.covGridList <- NULL
  }


  resultList <- BasicResultExtraction(pp.res = pp.res,
                                      inlaPrepList = inlaPrepList,
                                      projgrid = gridList$projgrid,
                                      use.covariates = use.covariates,
                                      covariate.fitting = covariate.fitting,
                                      additional.iid.term = additional.iid.term,
                                      covariateValues = covGridList$covariateValues,
                                      logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain,
                                      nxy = gridList$nxy,
                                      extraNonlinear.covGridList = extraNonlinear.covGridList)


  print("Finished initial result extraction")



  #image.plot(x=projgrid$x, y=projgrid$y, z=resultList$mean.field,main="Mean of latent field (log-scale)",nlevel=200)
  #lines(countingDomain,lwd=2)
  #aa=resultList$mean.field
  #aa[!logicalGridPointsInsideCountingDomain] = NA
  #image.plot(x=projgrid$x, y=projgrid$y, z=aa,main="Mean of latent field (log-scale)",nlevel=200)

  #### Samples from the the posterior field ####

  set.seed(123)

  print("Starting INLA posterior sample")

  samp <- inla.posterior.sample(n=noSamp,result=pp.res) # This typically takes a while

  print("Finished running INLA posterior sample")

  #### Handling the posterior  ####


  if (CVFold[[2]]>0){ # Does NOT delete the samp object
    if (use.covariates){
      postPredDistListPredPhotos <- ComputePostPredDistforPredPhotos(samp = samp,
                                                                     thisMeshPoint = thisMeshPoint,
                                                                     covAtPredLoc = covAtPredLoc,
                                                                     areaPredPhotos = areaPredPhotos,
                                                                     poisson.maxEvals = poisson.maxEvals,
                                                                     parallelize.numCores = parallelize.numCores,
                                                                     parallelize.noSplits = parallelize.noSplits,
                                                                     covariate.fitting = covariate.fitting,
                                                                     mesh = mesh)
    }
  }


  # NOTE: This function deletes the samp object!

  if (sum(leaveOutTransect)>0){
    postPredDistList <- PhotoPostPredDist(samp = samp,
                                          spatial = spatial,
                                          parallelize.noSplits = parallelize.noSplits,
                                          parallelize.numCores = parallelize.numCores,
                                          tempFolder = tempFolder,
                                          use.covariates = use.covariates,
                                          additional.iid.term = additional.iid.term,
                                          covariate.fitting = covariate.fitting,
                                          gridList = gridList,
                                          covGridList = covGridList,
                                          nxy = gridList$nxy,
                                          noMeshPoints = mesh$n,
                                          poisson.maxEvals = poisson.maxEvals,
                                          extraNonlinear.covGridList = extraNonlinear.covGridList,
                                          predPhotoGridPointList = predPhotoGridPointList)
  } else {
    postPredDistList <-  ComputePostPredDist(samp = samp,
                                             spatial = spatial,
                                             parallelize.noSplits = parallelize.noSplits,
                                             parallelize.numCores = parallelize.numCores,
                                             tempFolder = tempFolder,
                                             use.covariates = use.covariates,
                                             additional.iid.term = additional.iid.term,
                                             covariate.fitting = covariate.fitting,
                                             gridList = gridList,
                                             covGridList = covGridList,
                                             nxy = gridList$nxy,
                                             noMeshPoints = mesh$n,
                                             poisson.maxEvals = poisson.maxEvals,
                                             extraNonlinear.covGridList = extraNonlinear.covGridList)
  }

  print("Finished running last parallel session: All Poisson distributions computed")

  #### Computes basic summary statistics ####

  if (CVFold[[2]]>0){

    photoOrder <- 1:nrow(predPhotos) # order(dataList$org$coordPhotoX[!includeThesePhotos])

    ComparePhotoCountsAndPred(posteriorevalFullPhoto = postPredDistListPredPhotos$evalFullPredPhoto,
                              posteriorDistPhoto = postPredDistListPredPhotos$posteriorDistPredPhoto,
                              posteriorevalFullTransect = postPredDistListPredPhotos$evalFullAllPredPhoto,
                              posteriorDistTransect = postPredDistListPredPhotos$posteriorDistAllPredPhoto,
                              areaPerSubSampPhotoMatrix = t(postPredDistListPredPhotos$areaPerPredPhoto),
                              photoCounts = dataList$org$noObsPerPhoto[predThesePhotosLogical],
                              photoOrder = photoOrder,
                              leaveOutTransect = CVFold[[1]], # For CV-validation this is rather the CV-fold number and has nothings to do with the transect number.
                              inputVar = inputVar,
                              savingFolder = savingFolder)
  }


    if (sum(leaveOutTransect)>0){

    photoOrder <- order(dataList$org$coordPhotoX[predThesePhotosLogical])

    ComparePhotoCountsAndPred(posteriorevalFullPhoto = postPredDistList$posteriorevalFullPhoto,
                              posteriorDistPhoto = postPredDistList$posteriorDistPhoto,
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

  #### Write plots to pdf ####

  SummaryPlotFunc(meshplot = TRUE,
                  boundariesplot = !standardMesh,
                  covariatesplot = use.covariates,
                  summaryplot = TRUE,
                  savingFolder = savingFolder,
                  sealPhotoDataFile = sealPhotoDataFile,
                  sealTransectDataFile = sealTransectDataFile,
                  dataList = dataList,
                  orgPhotos = orgPhotos,
                  modPhotos = modPhotos,
                  gridList = gridList,
                  finalResList = finalResList,
                  mesh = mesh,
                  covMesh = covMesh,
                  tilesList = voronoiTess$tiles,
                  weightAtMeshLoc = weightAtMeshLoc,
                  countingDomain = countingDomain,
                  logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain,
                  covNewGridval = covGridList$covariateValues,
                  pp.res = pp.res,
                  sealType = sealType,
                  use.covariates = use.covariates,
                  covariates.type = covariates.type,
                  covariate.fitting = covariate.fitting,
                  spatial = spatial,
                  additional.iid.term = additional.iid.term,
                  convHullVar.convex = convHullVar.convex,
                  convHullVar.concave = convHullVar.concave,
                  convHullVar.resolution = convHullVar.resolution,
                  meshVar.max.edge = meshVar.max.edge,
                  meshVar.offset = meshVar.offset,
                  meshVar.cutoff = meshVar.cutoff,
                  Matern.alpha = Matern.alpha,
                  grid.pixelsize = grid.pixelsize,
                  INLA.theta.startval = INLA.theta.startval,
                  parallelize.noSplits = parallelize.noSplits,
                  parallelize.numCores = parallelize.numCores,
                  poisson.maxEvals = poisson.maxEvals,
                  noSamp = noSamp,
                  time = timeUsed,
                  testing = testing,
                  comment = comment,
                  leaveOutTransect = leaveOutTransect)

  #### Saves RData to file ####

  if (save.data){
    save.list <- c("dataList","mesh","finalResList","savingFolder","inputVar","gridList","covGridList","timeUsed","countingDomain","areaCountDomain","spde","pp.res")
    save(list=save.list,file=file.path(savingFolder,"output.RData"))
  }

  if (delete.temp){
    unlink(tempFolder,recursive=TRUE)
  }



  #### End of function ####

  cat(paste("Finished running the INLAPPSealsPoissonReg function. Output found in \n\n",savingFolder,sep=""))

}







