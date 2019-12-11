
#' Seal counting using INLA and continuous field Poisson process formulation (work in progress)
#'
#' Approximates a Log-Gaussian Cox process by a continuous field approximation (Simpson et. al (2016)), and fits seal counts from transects to it.
#' Samples from the the fitted model to compute the posterior predictive distribution of the total number of seals within an
#' area which the transects are spanning
#'
#' @param resultsBaseFolder String, indicating the path to the base folder where all results are to be stored
#' @param sealPhotoDataFile String, indicating the file where seal photo data are stored
#' @param sealTransectDataFile String, indicating the file where seal transect data are stored
#' @param satelliteDataFolder String, indicating the path to the folder where the satellite data are stored (in files on the form "cov_grid_[covariates.type].rds")
#' @param sealType String, indicating which seal type to model "hooded" (default as it is quicker), or "harps"
#' @param spatial Logical indicating whether a spatial spde model should be used
#' @param covariate.fitting String, indicating how to model covariates. "linear", quadratic (default) or "linAndLog", or FALSE for no covariates
#' @param covariates.type String equal to "band1" or "band2" indicating which of the bands from the satellite data should be used.
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param observationMethod String, indicating which method to use when defining the locations. One of "asIs", "repeatCenter" or "randomWithinPhoto".
#'        "asIs" takes the total number of counts in the center of each photo.
#'        "repeatedCenter" repates the center location the same number of times as there are observations
#'        "randomWithinPhoto" samples the location uniformly within each photo.
#' @param intPointMethod String, indicating which method to use when defining the integration points. One of "Fabian" or "mesh".
#'        "Fabian" creates a new mesh for the purpose of selecting the integration points, located at the border of observation domain.
#'        "meshPoints" uses the original mesh nodes as integration points.
#' @param intPointWeightMethod String, indicating which method to use when defining the weights for the integration points. One of "Fabian" or "Voronoi".
#'        "Fabian" gives 1/3 of the area of each triangle as weight to the corner points (using the new mesh points giving triangles either completely insside or outside of the observation domain)
#'        "Voronoi" uses the area of the Voronoi tiles (inside the observation domain) at each mesh node (with either integration point method) as the weight for each integration point.
#' @param noSamp Numeric, indicating the number of samples from posterior model to use when computing the posterior predictive distribution. 5000 (default) typically sufficient.
#' @param convHullVar.convex Numeric, corrresponding to the convex parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.concave Numeric, corrresponding to the concave parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.resolution Numeric vector of length 2, corrresponding to the resolution parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param meshVar.max.edge Numeric, corrresponding to the max.edge parameter of inla.mesh.2d. Smaller values gives smaller triangles outside triangles. See ?inla.mesh.2d for details
#' @param meshVar.offset Numeric, corrresponding to the offset parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param meshVar.cutoff Numeric, corrresponding to the cutoff parameter of inla.mesh.2d. See ?inla.mesh.2d for details
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


INLAPPSealsContLikApprox <- function(resultsBaseFolder = "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/PoissonContLikApprox",
                                     sealPhotoDataFile = "/nr/project/stat/PointProcess/Data/Seals/OriginalSealsKm.rds",
                                     sealTransectDataFile = "/nr/project/stat/PointProcess/Data/Seals/OigardTablesTransformed.rds",
                                     satelliteDataFolder = "/nr/project/stat/PointProcess/Data/Seals/Satellite/",
                                     sealType = "hooded",
                                     spatial = TRUE,
                                     covariate.fitting = "quadratic",
                                     covariates.type = "band1",
                                     additional.iid.term = FALSE,
                                     observationMethod = "randomWithinPhoto",
                                     intPointMethod = "meshPoints",
                                     intPointWeightMethod = "Voronoi",
                                     noSamp = 5000,
                                     convHullVar.convex = -0.15,
                                     convHullVar.concave = convHullVar.convex,
                                     convHullVar.resolution = c(120,120),
                                     meshVar.max.edge = c(2,10),
                                     meshVar.offset = 6,
                                     meshVar.cutoff = 0.2,
                                     Matern.alpha = 2,
                                     grid.pixelsize = 0.2,
                                     INLA.theta.startval = NULL,
                                     INLA.verbose = FALSE,
                                     parallelize.numCores = 2,
                                     parallelize.noSplits = parallelize.numCores,
                                     poisson.maxEvals = 5*10^5,
                                     results.CI.level = 0.95,
                                     additional_comment = "",
                                     save.data = TRUE,
                                     delete.temp = TRUE,
                                     testing = FALSE){

#### Initial definitions ####

  runType <- "ContLikApprox" # Setting the type of function

  # Starting time measuring
  time0 <- proc.time()

  # Computing use.covariates variable
  if (covariate.fitting %in% c("linear", "quadratic", "linAndLog")){
    use.covariates <- T
  } else {
    use.covariates <- F
  }

  # Creating a vector with the input variables
  inputVar <- names(formals(INLAPPSealsContLikApprox))

  # A list with the input variables to be used directly
  #inputList <- list()
  #for (i in 1:length(inputVar)){
  #  eval(parse(text=paste("inputList$",inputVar[i]," <- ",eval(parse(text=inputVar[i])),sep="")))
  #}

#### Initial folder and seed setup ####

  seed <- sum(as.numeric(proc.time()),na.rm = TRUE)*1000 # Just to get a random seed (I get the same every time I open R in Thinlinc)
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
  dataList$noPhoto <- dim(seals)[1]   # Total number of photo taken (with or without seals)
  dataList$coordPhotoX <- seals$xkm   # X coordinate of center of each photo, given in kilometers
  dataList$coordPhotoY <- seals$ykm   # Y coordinate of center of each photo,  given in kilometers
  dataList$photoWidth <- seals$lengthkm # The width (X-direction) of each of the photos, given in kilometers
  dataList$photoHeight <- seals$widthkm # The height (Y-direction) of each of the photos, given in kilometers

  ## Number of observed seals of correct type for each photo
  if (sealType=="hooded") dataList$noObsPerPhoto <- seals$hooded
  if (sealType=="harps")  dataList$noObsPerPhoto <- seals$harps

  ## Loading the transect seal data
  transectData <- readRDS(sealTransectDataFile)

  dataList$noTransects <- nrow(transectData)
  dataList$transectStartCoordX <- transectData$x.start
  dataList$transectStartCoordY <- transectData$y.start
  dataList$transectEndCoordX <- transectData$x.end
  dataList$transectEndCoordY <- transectData$y.end

  ## Loading the satellite covariate data, if applicable

  if (covariates.type%in%c("band1","band2")){
    covGrid <- readRDS(file.path(satelliteDataFolder,paste("cov_grid_",covariates.type,".rds",sep="")))
  }



  print("Finished loading and preparing seal data")

#### Creat polygon defining the counting domain for the seals (based on the area spanned by the transects) ####

  # The locations of each corner of the photo
  photos <- list()
  photos$pic.x.left <- dataList$coordPhotoX-0.5*dataList$photoWidth
  photos$pic.x.right <- dataList$coordPhotoX+0.5*dataList$photoWidth
  photos$pic.y.bottom <- dataList$coordPhotoY-0.5*dataList$photoHeight
  photos$pic.y.top <- dataList$coordPhotoY+0.5*dataList$photoHeight
  photos <- as.data.frame(photos)

  observationDomain <- as.data.frame(MergePolygons(photos))


  countingDomain <- CreateCountDomainPolygon(transectStartCoordX = dataList$transectStartCoordX,
                                             transectStartCoordY = dataList$transectStartCoordY,
                                             transectEndCoordX = dataList$transectEndCoordX,
                                             transectEndCoordY = dataList$transectEndCoordY,
                                             coordPhotoX = dataList$coordPhotoX,
                                             coordPhotoY = dataList$coordPhotoY,
                                             photoWidth = dataList$photoWidth,
                                             photoHeight = dataList$photoHeight,
                                             transectYSpan = 1.5*1.852)

  countingDomainExpanded <- CreateCountDomainPolygon(transectStartCoordX = dataList$transectStartCoordX,
                                                     transectStartCoordY = dataList$transectStartCoordY,
                                                     transectEndCoordX = dataList$transectEndCoordX,
                                                     transectEndCoordY = dataList$transectEndCoordY,
                                                     coordPhotoX = dataList$coordPhotoX,
                                                     coordPhotoY = dataList$coordPhotoY,
                                                     photoWidth = dataList$photoWidth,
                                                     photoHeight = dataList$photoHeight,
                                                     transectYSpan = 2*1.852,
                                                     transectXSpan = 0.5*1.852)


  areaCountDomain <- rgeos::gArea(coo2sp(countingDomain)) # Should be close to 4569 which Tor-Arne uses in his papers


  print("Finished computation of counting domain")

#### Cunstructing the mesh ####
  mesh <- IndepMeshCreation(rectangleCentersX = dataList$coordPhotoX,
                            rectangleCentersY = dataList$coordPhotoY,
                            countingDomainExpanded = countingDomainExpanded,
                            convHullVar.convex = convHullVar.convex,
                            convHullVar.concave = convHullVar.convex,
                            convHullVar.resolution = convHullVar.resolution,
                            meshVar.max.edge = meshVar.max.edge,
                            meshVar.offset = meshVar.offset,
                            meshVar.cutoff = meshVar.cutoff)

  print("Finished constructing the mesh")

#### Defining integration points ####

  # For now using a single connected countingDomain

  if(intPointMethod == "Fabian"){
    int.Fabian <- int.polygon.MJ(mesh,loc=cbind(x=observationDomain$X,y=observationDomain$Y),group=observationDomain$SID,parallelize.numCores=parallelize.numCores)
    intPoints <- data.frame(x=int.Fabian$int.polygon$x,y=int.Fabian$int.polygon$y)

    if (intPointWeightMethod == "Fabian"){
      intPoints$weight <- int.Fabian$int.polygon$weight
    }
  }

  if(intPointMethod == "meshPoints"){
    intPoints <- data.frame(x=mesh$loc[,1],y=mesh$loc[,2])

    if (intPointWeightMethod == "Fabian"){
      print("weightMethod=Fabian cannot be combined with intPointMethod=meshPoints. Using weightMethod=Voronoi instead.")
      weightMethod = "Voronoi"
    }
  }

  if (intPointWeightMethod == "Voronoi"){
    int.Voronoi <- CreateVoronoiTessellation(locationsCoordX=intPoints$x,
                                             locationsCoordY=intPoints$y,
                                             observationDomain=observationDomain)
    intPoints$weight <- int.Voronoi$tileSize
  }
  intPoints <- intPoints[intPoints$weight > 0,]

  print("Finished defining the integration points and their weights")


  # plot(mesh,xlim=c(-2,2),ylim=c(-2,2))
  #  symbols(x=intPoints$x,y=intPoints$y,circles = sqrt(intPoints$weight/pi),inches=F,add=T,lwd=2,fg=3)
  #  symbols(x=intPoints0$x,y=intPoints0$y,circles = sqrt(intPoints0$weight/pi),inches=F,add=T,fg=2,lwd=2)

#### Defining the observation ####

  if (observationMethod == "asIs"){ # Should give same likelihood as "asIs"
    thesePositive <- dataList$noObsPerPhoto>0
    obsPoints <- data.frame(x = dataList$coordPhotoX[thesePositive], y = dataList$coordPhotoY[thesePositive], count = dataList$noObsPerPhoto[thesePositive])
  }
  if (observationMethod == "repeatCenter"){ # Should give same likelihood as "asIs"
    obsPoints <- NULL
    for (i in 1:dataList$noPhoto){
      reps <- dataList$noObsPerPhoto[i]
      if (reps>0){
        obsPoints <- rbind(obsPoints,cbind(rep(dataList$coordPhotoX[i],reps),rep(dataList$coordPhotoY[i],reps), rep(1,reps) ))
      }
    }
    obsPoints <- as.data.frame(obsPoints)
    names(obsPoints) <- c("x","y","count")
  }

  if (observationMethod == "randomWithinPhoto"){
    obsPoints <- NULL
    for (i in 1:dataList$noPhoto){
      reps <- dataList$noObsPerPhoto[i]
      if (reps>0){
        obsPoints <- rbind(obsPoints,cbind(runif(reps,
                                                 min = dataList$coordPhotoX[i]-0.5*dataList$photoWidth[i],
                                                 max = dataList$coordPhotoX[i]+0.5*dataList$photoWidth[i]),
                                           runif(reps,
                                                 min = dataList$coordPhotoY[i]-0.5*dataList$photoHeight[i],
                                                 max = dataList$coordPhotoY[i]+0.5*dataList$photoHeight[i]),
                                           rep(1,reps)))
      }
    }
    obsPoints <- as.data.frame(obsPoints)
    names(obsPoints) <- c("x","y","count")
  }

  print ("Finished defining the observations")

#### Defining the fake data ####

  inlaPrepList <-   PrepareINLAFuncContLikApprox(mesh = mesh,
                                                 intPoints = intPoints,
                                                 obsPoints = obsPoints,
                                                 covGrid = covGrid,
                                                 spatial = spatial,
                                                 use.covariates = use.covariates,
                                                 covariate.fitting = covariate.fitting,
                                                 additional.iid.term = additional.iid.term,
                                                 Matern.alpha = Matern.alpha,
                                                 covariates.type = covariates.type ,
                                                 INLA.theta.startval = INLA.theta.startval)

  dataList$y.pp <- inlaPrepList$stk.pp$data$data$y
  dataList$e.pp <- inlaPrepList$stk.pp$data$data$e
  dataList$obsPoints <- obsPoints
  dataList$intPoints <- intPoints

  if(spatial){
    spde <- inlaPrepList$spde
  } else {
    spde <- NULL
  }

  print("Finished building INLA stack and preparing the INLA run")

#### Running the INLA function ####

  ss <- proc.time()

  pp.res <- INLA::inla(inlaPrepList$formula,
                       family="poisson", data=INLA::inla.stack.data(inlaPrepList$stk.pp),
                       control.predictor=list(A=INLA::inla.stack.A(inlaPrepList$stk.pp)),
                       E=INLA::inla.stack.data(inlaPrepList$stk.pp)$e,verbose=INLA.verbose,
                       control.compute=inlaPrepList$control.compute.list,
                       control.mode = inlaPrepList$control.mode.list)

  runINLATime <- proc.time()-ss

  print(paste("Finished running the INLA-function in ",round(runINLATime[3]/60)," minutes",sep=""))

#### Extracting results from the INLA model, plots and estimates the range parameter ####

  gridList <- GridCreation(mesh = mesh,
                           countingDomain = countingDomain,
                           areaCountDomain = areaCountDomain,
                           grid.pixelsize = grid.pixelsize)


  covGridList <- covAtNewGrid(use.covariates = use.covariates,
                              covGrid = covGrid,
                              gridvalX = gridList$gridvalX,
                              gridvalY = gridList$gridvalY,
                              modelledGridVal = gridList$modelledGridVal,
                              logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain)

  resultList <- BasicResultExtraction(pp.res = pp.res,
                                      inlaPrepList = inlaPrepList,
                                      projgrid = gridList$projgrid,
                                      use.covariates = use.covariates,
                                      covariate.fitting = covariate.fitting,
                                      additional.iid.term = additional.iid.term,
                                      covariateValues = covGridList$covariateValues,
                                      logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain,
                                      nxy = gridList$nxy)

  print("Finished initial result extraction")

  #fields::image.plot(z=resultList$mean.field,main="Mean of latent field (log-scale)",nlevel=200)
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

  # NOTE: This function deletes the samp object!
  postPredDistList <- ComputePostPredDist(samp = samp,
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
                                          poisson.maxEvals = poisson.maxEvals)

  print("Finished running last parallel session: All Poisson distributions computed")

#### Computes basic summary statistics ####

  postResList <- SummaryStat(evalPoints = postPredDistList$posteriorEvalPoints,
                             dist = postPredDistList$posteriorDist,
                             results.CI.level = results.CI.level,
                             posterior = TRUE)


  finalResList <- c(resultList,postPredDistList,postResList)

  timeUsed <- proc.time()-time0

#### Write plots to pdf ####

  SummaryPlotFuncContLikApprox(meshplot = TRUE,
                               boundariesplot = TRUE,
                               covariatesplot = TRUE,
                               summaryplot = TRUE,
                               savingFolder = savingFolder,
                               sealPhotoDataFile = sealPhotoDataFile,
                               sealTransectDataFile = sealTransectDataFile,
                               results.CI.level = results.CI.level,
                               observationDomain = observationDomain,
                               dataList = dataList,
                               gridList = gridList,
                               finalResList = finalResList,
                               mesh = mesh,
                               countingDomain = countingDomain,
                               logicalGridPointsInsideCountingDomain = gridList$logicalGridPointsInsideCountingDomain,
                               covNewGridval = covGridList$covariateValues,
                               fixed.effects.vec = pp.res$summary.fixed$mean,
                               sealType = sealType,
                               use.covariates = use.covariates,
                               covariates.type = covariates.type,
                               covariate.fitting = covariate.fitting,
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
                               comment = comment)

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








#
#
#
#
#
#
#
#
#
#
#
#
# # The locations of each corner of the photo
# photos <- list()
# photos$pic.x.left <- dataList$coordPhotoX-0.5*dataList$photoWidth
# photos$pic.x.right <- dataList$coordPhotoX+0.5*dataList$photoWidth
# photos$pic.y.bottom <- dataList$coordPhotoY-0.5*dataList$photoHeight
# photos$pic.y.top <- dataList$coordPhotoY+0.5*dataList$photoHeight
# photos <- as.data.frame(photos)
#
# photoPolygons <- MergePolygons(photos)
#
# photosCovering <- photos
# #  photosCovering[,1] <- photosCovering[,1] - 0.5*1.852
# #  photosCovering[,2] <- photosCovering[,2] + 0.5*1.852
# photosCovering[,3] <- photosCovering[,3] - 1.5*1.852
# photosCovering[,4] <- photosCovering[,4] + 1.5*1.852
#
# countingDomainExact <- MergePolygons(photosCovering)
#
# #  PBSmapping::plotPolys(countingDomainExact,xlim=c(-2,2),ylim=c(-5,5)+20,bg=2,col=3)
# #  PBSmapping::plotPolys(countingDomainExact,xlim=c(-1,1),ylim=c(-50,50),bg=2)
#
#
# plot(mesh,xlim=c(-7,-4),ylim=c(7,11))
# points(intpoints$int.polygon$x,intpoints$int.polygon$y,pch="+",cex=2)
# #  symbols(x=intpoints$int.polygon$x,y=intpoints$int.polygon$y,circles = sqrt(intpoints$int.polygon$weight/pi),add=T,inches=F)
# plot(intpoints$imesh,add=T,edge.color=2,lwd=2)
# plot(mesh,add=T,lwd=2,edge.color=1)
#
# plot(mesh,xlim=c(-5.5,-4.5),ylim=c(8.5,9.5))
# symbols(x=intpoints$int.polygon$x,y=intpoints$int.polygon$y,circles = sqrt(intpoints$int.polygon$weight/pi),add=T,inches=F)
# plot(intpoints$imesh,add=T,edge.color=2,lwd=2)
# points(intpoints$int.polygon$x,intpoints$int.polygon$y,pch="+",cex=3)
#
# plot(mesh,add=T,lwd=2,edge.color=1)
#
#
#
# #### Creating the Voronoi tessellations and computing the area each of them cover within the and specfying the weights e.pp for each of the polygons (based on its size) ####
#
#
#
#
#
# voronoiTess <- CreateVoronoiTessellation(locationsCoordX = mesh$loc[,1],
#                                          locationsCoordY = mesh$loc[,2],
#                                          domainCoordX = countingDomain$x,
#                                          domainCoordY = countingDomain$y) ### Might consider replacing countingDomian with a somewhat larger modelling domain here instead
#
# weightAtMeshLoc <- voronoiTess$tileSize
#
# print("Finished computation of the Voronoi tesselation")
#
# ## Checking which tile the photo centers belong to and assign their values to them.
#
#
# yAtMeshLoc <- rep(0,mesh$n)
# for (i in 1:mesh$n){
#   insidePhotos <- which(as.logical(point.in.polygon(dataList$coordPhotoX,dataList$coordPhotoY,voronoiTess$tiles[[i]]$x,voronoiTess$tiles[[i]]$y)))
#
#   if (length(insidePhotos)==0){
#     yAtMeshLoc[i] <- NA
#   } else {
#     yAtMeshLoc[i] <- sum(dataList$noObsPerPhoto[insidePhotos])
#   }
# }
#
# ## Exclude all Voronoi tesselations where there are no photo centers inside, as these are actually unobserved
# dataList$NAMeshLoc <- which(is.na(yAtMeshLoc)) # Unobserved mesh locations
# dataList$y.pp <- yAtMeshLoc[-dataList$NAMeshLoc]
# dataList$e.pp <- weightAtMeshLoc[-dataList$NAMeshLoc]
# dataList$obsLoc <- mesh$loc[-dataList$NAMeshLoc,]
#
#
#
# # ## The locations of each corner of the photo
# # photos <- list()
# # photos$pic.x.left <- dataList$coordPhotoX-0.5*rectangleWidth
# # photos$pic.x.right <- dataList$coordPhotoX+0.5*rectangleWidth
# # photos$pic.y.bottom <- dataList$coordPhotoY-0.5*rectangleHeight
# # photos$pic.y.top <- dataList$coordPhotoY+0.5*rectangleHeight
# # photos <- as.data.frame(photos)
# # colcol = rep("white",mesh$n)
# # colcol[yAtMeshLoc==0] <- "blue"
# # colcol[yAtMeshLoc>0] <- "red"
#
# # plot(1,1,type='n',xlim=c(-10,10),ylim=c(-10,10))
# # plot(voronoiTess$tiles,col=2,xlim=c(0,20),ylim=c(0,20),add=TRUE,fillcol=colcol,showpoints=FALSE)
# # for (i in 1:dataList$noPhoto){
# #   lines(photos[i,c(1,2,2,1,1)],photos[i,c(3,3,4,4,3)],col=3)
# # }
#
# #### Preparing and executing the model fitting ####
#
# inlaPrepList <- PrepareINLAFunc(mesh = mesh,
#                                 obsLoc = dataList$obsLoc,
#                                 y.pp = dataList$y.pp,
#                                 e.pp = dataList$e.pp,
#                                 covGrid = covGrid,
#                                 use.covariates = use.covariates,
#                                 covariate.fitting = covariate.fitting,
#                                 additional.iid.term = additional.iid.term,
#                                 Matern.alpha = Matern.alpha,
#                                 covariates.type = covariates.type,
#                                 INLA.theta.startval = INLA.theta.startval)

