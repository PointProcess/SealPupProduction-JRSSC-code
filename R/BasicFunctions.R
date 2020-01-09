#' String extraction function
#'
#' Function extracting string between two specific characters, minor customization of this one
#' http://www.r-bloggers.com/how-to-extract-a-string-between-2-characters-in-r-and-sas/
#'
#' @param mystring Character vector to extract from.
#' @param initial.character Character determining the starting point of extractions
#' @param final.character Character determining the end point of extractions
#' @return snippet
#' @export


getstr = function(mystring, initial.character="_", final.character="_")
{
  # check that all 3 inputs are character variables
  if (!is.character(mystring))
  {
    stop('The parent string must be a character variable.')
  }

  if (!is.character(initial.character))
  {
    stop('The initial character must be a character variable.')
  }


  if (!is.character(final.character))
  {
    stop('The final character must be a character variable.')
  }

  add=0
  if(initial.character==final.character){add=1}

  # pre-allocate a vector to store the extracted strings
  snippet = rep(0, length(mystring))

  for (i in 1:length(mystring))
  {
    # extract the initial position
    initial.position = gregexpr(initial.character, mystring[i])[[1]][1] + 1

    # extract the final position
    final.position = gregexpr(final.character, mystring[i])[[1]][1+add] - 1

    # extract the substring between the initial and final positions, inclusively
    snippet[i] = substr(mystring[i], initial.position, final.position)
  }
  return(snippet)
}




#' Prepare for running INLA function with continuous likelihood approximation
#'
#' Gathers and stacks the input to the INLA function
#'
#' @param mesh Mesh object, being the output from inla.mesh.2d
#' @param intPoints Data frame with the integration point locations and their associated weights
#' @param obsPoints Data frame with the point pattern locations (locations with at least one seal observed)
#' @param covGrid, im object representing the covariate values on a dense grid where the counts live.
#' @param spatial Logical indicating whether a spatial spde model should be used
#' @param use.covariates Logical, indicating whether covariates should be used in the fitting (default is true)
#' @param covariate.fitting String, indicating how to model covariates. "linear", quadratic (default) or "linAndLog", or FALSE for no covariates
#' @param covariates.type String equal to "band1" or "band2" indicating which of the bands from the satellite data should be used.
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param Matern.alpha Numeric, corresponding to the alpha parameter in the Matern covariance function (2 is the default)
#' @param INLA.theta.startval List containing the start values for the theta parameter in the INLA run. (NULL indicates that automatic start values should be used, and is the default)
#' @return List with all variables necessary to run the INLA function, in addition to the spde object
#' @keywords inla
#' @import spatstat
#' @import INLA
#' @export

PrepareINLAFuncContLikApprox <- function(mesh,
                                         intPoints,
                                         obsPoints,
                                         covGrid,
                                         spatial = TRUE,
                                         use.covariates = TRUE,
                                         covariate.fitting = "quadratic",
                                         additional.iid.term = FALSE,
                                         Matern.alpha = 2,
                                         covariates.type = "band1",
                                         INLA.theta.startval = NULL) {

  y.pp <- c(rep(0,nrow(intPoints)),obsPoints$count)
  e.pp <- c(intPoints$weight,rep(0,nrow(obsPoints)))

  n.fakeData <- length(y.pp)

  if (spatial){
    A.ppInt <- INLA::inla.spde.make.A(mesh=mesh, loc = as.matrix(intPoints[,1:2]))
    A.ppObs <- INLA::inla.spde.make.A(mesh=mesh, loc = as.matrix(obsPoints[,1:2]))

    A.ppComb <- Matrix::rBind(A.ppInt,A.ppObs)

    spde <- INLA::inla.spde2.matern(mesh=mesh, alpha=Matern.alpha) # Constructing the Matern SPDE object
  } else {
    spde <- NULL
  }


  if (use.covariates) {

    nearestPixelIntPoints <- spatstat::nearest.pixel(intPoints[,1], intPoints[,2],covGrid)
    nearestPixelObsPoints <- spatstat::nearest.pixel(obsPoints[,1], obsPoints[,2],covGrid)

    covAtIntPoints <- covGrid[Reduce('cbind', nearestPixelIntPoints)]
    covAtObsPoints <- covGrid[Reduce('cbind', nearestPixelObsPoints)]

    covAtIntPoints2 <- covAtIntPoints^2
    covAtObsPoints2 <- covAtObsPoints^2

    covAtIntPointsLog <- log(covAtIntPoints)
    covAtObsPointsLog <- log(covAtObsPoints)
  }


  ## Setting up the inla stack:

# First the direct variables

  direct.A.list <- 1

  if (!use.covariates){
    direct.effects.list <- list(intercept=1)
    direct.formula <- "0 + intercept"
  }

  if (use.covariates){
    if (covariate.fitting=="linear"){
      direct.effects.list <- list(intercept=1, covariate = c(covAtIntPoints,covAtObsPoints))
      direct.formula <- "0 + intercept + covariate"
    }
    if (covariate.fitting=="quadratic"){
      direct.effects.list <- list(intercept=1, covariate = c(covAtIntPoints,covAtObsPoints),
                                covariate2 = c(covAtIntPoints2,covAtObsPoints2))
      direct.formula <- "0 + intercept + covariate + covariate2"
    }
    if (covariate.fitting=="linAndLog"){
      direct.effects.list <- list(intercept=1, covariate = c(covAtIntPoints,covAtObsPoints),
                                covariateLog = c(covAtIntPointsLog,covAtObsPointsLog))
      direct.formula <- "0 + intercept + covariate + covariateLog"
    }
  }

# Spatial variables

  if (spatial){
    spatial.formula <- "f(rf,model=spde)"
    spatial.A.list <- A.ppComb
    spatial.effects.list <- list(rf=1:mesh$n)
  } else {
    spatial.formula = spatial.A.list = spatial.effects.list = NULL
  }

# Additional iid term

  if (additional.iid.term){
    iid.formula <- "f(iid,model='iid')"
    iid.A.list <- A.ppObs
    iid.effects.list <- list(iid=1:nrow(obsPoints))
  } else {
    iid.formula = iid.A.list = iid.effects.list = NULL
  }


  A.list <- list(direct.A.list,spatial.A.list,iid.A.list)
  A.list <- A.list[!sapply(A.list,is.null)]

  effects.list <- list(direct.effects.list,spatial.effects.list,iid.effects.list)
  effects.list <- effects.list[!sapply(effects.list,is.null)]

  formula.list <- list(direct.formula,spatial.formula,iid.formula)
  formula.list <- formula.list[!sapply(formula.list,is.null)]


  stk.pp <- INLA::inla.stack(data=list(y=y.pp, e=e.pp),
                             A=A.list,
                             tag='pp',
                             effects=effects.list)


  formula = as.formula(paste('y ~ ',paste(unlist(formula.list),collapse=' + '),sep=''))
  control.compute.list <- list(config = TRUE, dic = TRUE, waic = TRUE, mlik = TRUE) # control.compute=list(config = TRUE) is needed for posterior sampling of the posterior random field


  if (all(is.null(INLA.theta.startval))){
    control.mode.list <- INLA::inla.set.control.mode.default()
  } else {
    control.mode.list <- list(theta = INLA.theta.startval,restart = TRUE) # Use the start values for the (internally specified) theta parameters here
  }

  retList <- list()
  retList$stk.pp <- stk.pp
  retList$formula <- formula
  retList$control.compute.list <- control.compute.list
  retList$control.mode.list <- control.mode.list
  retList$spde <- spde
  return(retList)
}






#' Split lines at mesh edges
#'
#' By Fabian Bachl
#'
#' @aliases split.lines
#' @export
#' @param mesh An inla.mesh object
#' @param sp Start points of lines
#' @param ep End points of lines
#' @param filter.zero.length Filter out segments with zero length? (Bool)
#' @param ... argments to int.quadrature
#' @return List of start and end points resulting from splitting the given lines
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

split.lines = function(mesh, sp, ep, filter.zero.length = TRUE) {

  # locations for splitting
  loc = as.matrix(rbind(sp,ep))
  idx = 1:dim(sp)[1]

  # Filter out segments not on the mesh
  t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
  t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
  # if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
  sp = sp[!((t1==0) | (t2==0)),]
  ep = ep[!((t1==0) | (t2==0)),]
  idx = idx[!((t1==0) | (t2==0))]
  loc = as.matrix(rbind(sp,ep))

  # Split them segments into parts
  if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
  np = dim(sp)[1]
  sp.idx = t(rbind(1:np,np+1:np))
  splt = INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
  #plot(data$mesh)
  #points(loc)
  #points(splt$split.loc,col="blue)

  sp = splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]] # Start point of new segments
  ep = splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]] # End points of new segments
  idx = idx[splt$split.idx[,1]]
  origin = splt$split.origin

  # Filter out zero length segments
  if ( filter.zero.length ) {
    sl = apply((ep-sp)^2,MARGIN=1,sum)
    sp = sp[!(sl==0),]
    ep = ep[!(sl==0),]
    origin = origin[!(sl==0)]
    idx = idx[!(sl==0)]
  }

  return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))

}


#' Integration points for polygons inside an inla.mesh usign parallelization
#'
#' Parallelized function for finding integration points for a given mesh and locations defining the polygons where observations are being made.
#' Edited version of Fabian Bachls code.
#'
#' @aliases int.polygon
#' @export
#' @param mesh An inla.mesh object
#' @param loc Locations defining the polygons
#' @param group If loc defines multiple polygons then this is the ID of the group for each location in loc

int.polygon.MJ = function(mesh, loc, group = NULL,parallelize.numCores=2){

  if ( is.null(group) ) { group = rep(1, nrow(loc)) }
  ipsl = list()
  print(paste0("Number of polygons to integrate over: ", length(unique(group)) ))


  doMC::registerDoMC(parallelize.numCores)

  export.var=c("mesh","loc","group") # Not functions here
  non.export.var=ls()[!(ls()%in%export.var)]

  uu=proc.time()

  ipsl <- foreach::foreach(i=unique(group),.noexport=non.export.var,.packages="INLA",.verbose=FALSE,.inorder=TRUE) %dopar% {
    g = i
    gloc = loc[group==g, ]
    # Check where the polygon intersects with the mesh edges
    sp = gloc[1:nrow(gloc),]
    ep = rbind(gloc[2:nrow(gloc),], gloc[1,])
    sloc = split.lines(mesh, sp, ep, filter.zero.length = FALSE)$split.loc[,1:2]
    if (!is.null(sloc)){ colnames(sloc) = colnames(loc) }
    #    plot(mesh,xlim=range(sloc[,1])+c(-0.5,0.5),ylim=range(sloc[,2])+c(-0.5,0.5))
    #plot(mesh,xlim=c(-15,-5),ylim=c(-5,5))
    #lines(loc)
    #points(gloc,pch=3,col=2,cex=2)

    #points(sloc)

    #points(sp,col=2,cex=1) ; points(ep,col=3,cex=2) ; points(sloc,col=4,cex=3)
  #  axis(side=1);axis(side=2)
    bloc = sloc[,1:2]    ### Martin uses only sloc here, and does not care about gloc as it is already included...
    #points(bloc,col=5,pch=2,cex=1)
    bnd = INLA::inla.mesh.segment(loc = bloc)
    imesh = INLA::inla.mesh.create(boundary = bnd, loc = mesh$loc[,1:2])
    #   plot(mesh);
    #  plot(imesh, add = TRUE,col=5)
    ips = data.frame(imesh$loc[,1:2])
    colnames(ips) = colnames(gloc)
    ips$weight = diag(as.matrix(INLA::inla.mesh.fem(imesh)$c0))
    ips$group = g
    print(g)
    ips
  }
  print(paste("Finished computing integration points in ",round((proc.time()-uu)[3]/60), " minutes."))

  retList <- list()
  retList$int.polygon <- do.call(rbind,ipsl)
#  retList$imesh <- imesh
  return(retList)
}



#' Integration points for polygons inside an inla.mesh
#'
#' Function for finding integration points for a given mesh and locations defining the polygons where observations are being made.
#' Created by Fabian Bachl.
#'
#' @aliases int.polygon
#' @export
#' @param mesh An inla.mesh object
#' @param loc Locations defining the polygons
#' @param group If loc defines multiple polygons then this is the ID of the group for each location in loc
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.polygon = function(mesh, loc, group = NULL){

  if ( is.null(group) ) { group = rep(1, nrow(loc)) }
  ipsl = list()
  k = 1
  print(paste0("Number of polygons to integrate over: ", length(unique(group)) ))
  for ( g in unique(group) ) {
    gloc = loc[group==g, ]
    # Check where the polygon intersects with the mesh edges
    sp = gloc[1:nrow(gloc),]
    ep = rbind(gloc[2:nrow(gloc),], gloc[1,])
    sloc = split.lines(mesh, sp, ep, filter.zero.length = FALSE)$split.loc[,1:2]
    if (!is.null(sloc)){ colnames(sloc) = colnames(loc) }
#    plot(mesh,xlim=range(sloc[,1])+c(-0.5,0.5),ylim=range(sloc[,2])+c(-0.5,0.5))
    #plot(mesh,xlim=c(-15,-5),ylim=c(-5,5))
    #lines(loc)
    #points(gloc,pch=3,col=2,cex=2)

    #points(sloc)

    #points(sp,col=2,cex=1) ; points(ep,col=3,cex=2) ; points(sloc,col=4,cex=3)
    #axis(side=1);axis(side=2)
    bloc = rbind(gloc, sloc[,1:2])
   # points(bloc,col=5,pch=2,cex=1)
    bnd = INLA::inla.mesh.segment(loc = bloc)
    imesh = INLA::inla.mesh.create(boundary = bnd, loc = mesh$loc[,1:2])
  #   plot(mesh);
  #  plot(imesh, add = TRUE,col=5)
    ips = data.frame(imesh$loc[,1:2])
    colnames(ips) = colnames(gloc)
    ips$weight = diag(as.matrix(INLA::inla.mesh.fem(imesh)$c0))
    ips$group = g
    ipsl = c(ipsl, list(ips))
    print(k)
    k = k + 1
  }
  retList <- list()
  retList$int.polygon <- do.call(rbind,ipsl)
  retList$imesh <- imesh
  return(retList)
}



#' Basic data independent mesh creation
#'
#' Creates a 2D mesh with the inla.mesh.2d function with an inner and an outer domain with different triangle density.
#' The outer is defined as a larger nonconvex hull of the rectangleCenters while the inner is defined based on a slightly expanded counting domain
#'
#' @param rectangleCentersX Vector with X-coordinate of the centerpoint of each of the obligatory rectangles
#' @param rectangleCentersY Vector with X-coordinate of the centerpoint of each of the obligatory rectangles
#' @param rectangleWidth Vector with the width of each of the obligatory rectangles
#' @param rectangleHeight Vector with the height of each of the obligatory rectanglges
#' @param convHullVar.convex convex parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.concave concave parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.resolution resolution parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param meshVar.max.edge max.edge parameter of inla.mesh.2d. Smaller values gives smaller triangles outside triangles. See ?inla.mesh.2d for details
#' @param meshVar.offset offset parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param meshVar.cutoff cutoff parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @return An INLA mesh object
#' @keywords mesh
#' @import INLA
#' @export


IndepMeshCreation <- function(rectangleCentersX,
                              rectangleCentersY,
                              countingDomainExpanded,
                              convHullVar.convex = -0.15,
                              convHullVar.concave = convHullVar.convex,
                              convHullVar.resolution = c(120,120),
                              meshVar.max.edge = c(1,10),
                              meshVar.offset = c(6,6),
                              meshVar.cutoff = 0.3){


  #### Creating the mesh based on the positions in obligMeshLoc, with cutoff to remove close duplicates ####

  # The function makes larger trinalges outside these observations automatically.
  domain.outer  <- INLA::inla.nonconvex.hull(points=cbind(rectangleCentersX,rectangleCentersY),
                                       convex=convHullVar.convex,
                                       concave=convHullVar.concave,
                                       resolution=convHullVar.resolution)

  domain.inner <- INLA::inla.mesh.segment(loc = as.matrix(countingDomainExpanded))


  mesh <- INLA::inla.mesh.2d(boundary = list(domain.inner,domain.outer),
                       max.edge=meshVar.max.edge,
                       offset=meshVar.offset,
                       cutoff=meshVar.cutoff)

  return(mesh)
}








#' Summary plot
#'
#' Produces several different summary plots showing the results from INLA approach
#'
#' @param meshplot Logical, indicating whether a plot with the mesh should be produced
#' @param boundariesplot Logical, indicating whether plots with the used boundaries and domains should be produced
#' @param covariatesplot Logical, indicating whether various plots showing the covariates and their effects should be produced
#' @param summaryplot Logical, indicating whether a plot showing the final posterior distribution and mean and pointwise sd of the latent field should be produced
#' @param savingFolder String, indicating the complete path to where the plots are to be stored
#' @param sealPhotoDataFile String, indicating the file where seal photo data are stored
#' @param sealTransectDataFile String, indicating the file where seal transect data are stored
#' @param gridList List containing information about the grid, being the output of the function GridCreation
#' @param dataList List with various data variables
#' @param orgPhotos data frame with the coordinates of all original photos (both those used in the modelling and those not used)
#' @param modPhotos data frame with the coordinates of the photos used in the modelling
#' @param finalResList List gathering all final results
#' @param mesh Mesh object, being the output from inla.mesh.2d
#' @param covMesh Mesh object representing the mesh for the covariate when applicable, being the output from inla.mesh.1d
#' @param tilesList List containing all the voronoi tesselation tiles
#' @param weightAtMeshLoc Numeric vector with the weight assigned to each tile in the tileList
#' @param countingDomain data.frame containing x- and y-coordinates of the counting domain polygon
#' @param logicalGridPointsInsideCountingDomain Logical vector, indicating which of the grid elements are within the counting domain
#' @param covNewGridval Matrix with all covariate values at grid points
#' @param pp.res Result object from INLA run
#' @param sealType String, indicating which seal type to model "hooded" (default as it is quicker), or "harps"
#' @param use.covariates Logical, indicating whether covariates are used or not (see description!)
#' @param covariate.fitting String, indicating how to model covariates. "linear", quadratic (default) or "linAndLog", or FALSE for no covariates
#' @param spatial Logical indicating whether this model includes a spatial term or not.
#' @param covariates.type String equal to "band1" or "band2" indicating which of the bands from the satellite data should be used.
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param convHullVar.convex Numeric, corrresponding to the convex parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.concave Numeric, corrresponding to the concave parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.resolution Numeric vector of length 2, corrresponding to the resolution parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param meshVar.max.edge Numeric, corrresponding to the max.edge parameter of inla.mesh.2d. Smaller values gives smaller triangles outside triangles. See ?inla.mesh.2d for details
#' @param meshVar.offset Numeric, corrresponding to the offset parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param meshVar.cutoff Numeric, corrresponding to the cutoff parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param Matern.alpha Numeric, corresponding to the alpha parameter in the Matern covariance function (2 is the default)
#' @param grid.pixelsize Numeric, denoting the size (in km, in both x- and y-direction) of each pixel of the grid being used
#' @param INLA.theta.startval List containing the start values for the theta parameter in the INLA run. (NULL indicates that automatic start values should be used, and is the default)
#' @param parallelize.numCores Numeric, corresponding to the number of cores any parallelization should be run at
#' @param parallelize.noSplits Numeric, deciding how many sublists samp should be splitted into. Should be a multiple of parallelize.numCores for the highest efficiency. The larger number the less memory is used (and longer time).
#' @param poisson.maxEvals Numeric, corresponding to maximum number of points the Poisson distribution should be evaluated at (a much smaller number is typically used)
#' @param noSamp Numeric, indicating the number of samples from posterior model to use when computing the posterior predictive distribution. 5000 (default) typically sufficient.
#' @param time The number of seconds the function has been running
#' @param testing Logical, indicating whether the testing parts of this function should be used
#' @param comment String to add as a comment to the summary plot
#' @return Produces a set of plots and writes them as pdf to file. Nothing else is returned
#' @keywords plot
#' @export


SummaryPlotFunc <- function(meshplot = TRUE,
                            boundariesplot = TRUE,
                            covariatesplot = TRUE,
                            summaryplot = TRUE,
                            savingFolder,
                            sealPhotoDataFile,
                            sealTransectDataFile,
                            dataList,
                            orgPhotos,
                            modPhotos,
                            results.CI.level = 0.95,
                            gridList,
                            finalResList,
                            mesh,
                            covMesh,
                            tilesList,
                            weightAtMeshLoc,
                            countingDomain,
                            logicalGridPointsInsideCountingDomain,
                            covNewGridval,
                            pp.res,
                            sealType,
                            use.covariates,
                            covariates.type,
                            covariate.fitting,
                            spatial,
                            additional.iid.term,
                            convHullVar.convex,
                            convHullVar.concave,
                            convHullVar.resolution,
                            meshVar.max.edge,
                            meshVar.offset,
                            meshVar.cutoff,
                            Matern.alpha,
                            grid.pixelsize,
                            INLA.theta.startval,
                            parallelize.noSplits,
                            parallelize.numCores,
                            poisson.maxEvals,
                            noSamp,
                            time,
                            testing,
                            comment,
                            leaveOutTransect,
                            Data_2018 = F){

  ## Creating some helping variables
  obsMeshLoc <- mesh$loc[-dataList$mod$NAMeshLoc,]
  fixed.effects.mesh <- finalResList$intercept.mean + finalResList$covariate.mean*covNewGridval + finalResList$covariate2.mean*covNewGridval^2 + finalResList$covariateLog.mean*log(covNewGridval) + finalResList$nonlinearCovgrid.mean
  fixed.effects.domain <- fixed.effects.mesh
  fixed.effects.domain[!logicalGridPointsInsideCountingDomain] = NA

  covNewGridvalDomain <- covNewGridval
  covNewGridvalDomain[!logicalGridPointsInsideCountingDomain] = NA



  ## Producing a mesh plot

  if (meshplot){
    ## Plotting the mesh and defined boundaries
    pdf(file=file.path(savingFolder,"mesh.pdf"),width=6,height=6)
    par(mar=c(0,0,0,0))
    plot(mesh, asp=1, main='')
    points(dataList$mod$coordPhotoX[dataList$mod$noObsPerPhoto==0],dataList$mod$coordPhotoY[dataList$mod$noObsPerPhoto==0],col="red",cex=0.5,pch=16)
    points(dataList$mod$coordPhotoX[dataList$mod$noObsPerPhoto>0],dataList$mod$coordPhotoY[dataList$mod$noObsPerPhoto>0],col="green",cex=0.5,pch=16)
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("y=0","y>0"),col=c("red","green"),pch=16)
    legend("topleft",c("Model boundary","counting domain"),col=c("blue",6),lty=c(1),lwd=3)
    dev.off()
  }

  ## Producing 2 x 2 boundaries plot

  if (boundariesplot){

    # Fixing the weigh per photo such that it corresponds to what is used when fitting the model
    weightAtMeshLocUsed <- weightAtMeshLoc
    weightAtMeshLocUsed[dataList$mod$NAMeshLoc] <- 0

    pdf(file=file.path(savingFolder,"boundaries.pdf"),width=12,height=12)
    par(mfrow=c(2,2))

    ## Plot of the tiles and give their corresponding obesrvation value  -- complete area
    plot(mesh$loc,type='n',main="Observations")
    for (i in 1:mesh$n){
      lines(c(tilesList[[i]]$x, tilesList[[i]]$x[1]), c(tilesList[[i]]$y, tilesList[[i]]$y[1]),col=1,lwd=0.5)
    }

    points(obsMeshLoc[dataList$mod$y.pp==0,],col="red",cex=0.5,pch=16)
    points(obsMeshLoc[dataList$mod$y.pp>0,],col="green",cex=0.5,pch=16)
    legend("bottomright",c("y=0","y>0"),col=c("red","green"),pch=16,bg="white")

    ## Plot of the tiles and give their corresponding weight value ( e ) -- complete area
    plot(mesh$loc,type='n',main="Weights/offset ( e )")
    for (i in 1:mesh$n){
      lines(c(tilesList[[i]]$x, tilesList[[i]]$x[1]), c(tilesList[[i]]$y, tilesList[[i]]$y[1]),col=1,lwd=0.5)
    }
    points(mesh$loc[weightAtMeshLocUsed==0,],col="blue",cex=0.75,pch=16)
    points(mesh$loc[weightAtMeshLocUsed>0,],col="purple",cex=0.75,pch=16)


    legend("bottomright",c("e=0","e>0"),col=c("blue","purple"),pch=16,bg="white")

    ## Plot of the tiles and give their corresponding obesrvation value  -- zooms in on speciic area
    colcol=rep("white",length(tilesList))
    colcol[-dataList$mod$NAMeshLoc][dataList$mod$y.pp==0]="red"
    colcol[-dataList$mod$NAMeshLoc][dataList$mod$y.pp>0]="green"

    par(mar=c(2,2,1,1), mgp=2:0)
    if(Data_2018){
      plot(1,1,xlim=c(0,20), ylim=c(-20,0),type='n',main="Observations")
    } else {
      plot(1,1,xlim=c(25,57), ylim=c(-20,5)-3,type='n',main="Observations")
    }
    plot(tilesList,fillcol=colcol,showpoints=FALSE,add=TRUE)
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("y=NA","y=0","y>0"),col=c("white","red","green"),pch=15,bg="white")
    legend("bottomleft",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")

    colcol=rep("red",length(tilesList))
    colcol[weightAtMeshLocUsed==0]="blue"
    colcol[weightAtMeshLocUsed>0]="purple"

    par(mar=c(2,2,1,1), mgp=2:0)
    if(Data_2018){
      plot(1,1,xlim=c(0,20), ylim=c(-20,0),type='n',main="Weights/offset ( e )")
    } else {
      plot(1,1,xlim=c(25,57), ylim=c(-23,2),type='n',main="Weights/offset ( e )")
    }
    plot(tilesList,fillcol=colcol,showpoints=FALSE,add=TRUE)
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("e=0","e>0"),col=c("blue","purple"),pch=15,bg="white")
    legend("bottomleft",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")

    dev.off()

  }

  ## Producing 4 1x2 covariate plots

  if (covariatesplot){
    par(mar=c(2,2,2,2), mgp=2:0)

    ## For covariate_effect_plot.pdf
    rangeCov <- range(covNewGridval,na.rm=T)
    covVal <- seq(rangeCov[1],rangeCov[2],length.out = 500)
    linearEffect <- finalResList$intercept.mean + finalResList$covariate.mean*covVal + finalResList$covariate2.mean*covVal^2 + finalResList$covariateLog.mean*log(covVal^2)

    if (covariate.fitting=="nonlinear"){
      covProj <-  inla.mesh.projector(covMesh,loc = covVal)
      nonlinearEffect <-  inla.mesh.project(covProj,pp.res$summary.random$nonlinear$mean)
    } else {
      nonlinearEffect <- 0
    }
    ##

    pdf(file=file.path(savingFolder,"covariate_effect_plot.pdf"),width=7,height=4)
    plot(covVal,linearEffect+nonlinearEffect, type='l',main="Covariate effect",xlab="Satellite image value",ylab="Log-intensity effect")
    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_mesh_with_obs.pdf"),width=14,height=6)
    par(mfrow=c(1,2),oma=c(0,0,0,1.5))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=covNewGridval,main="Covariate values",nlevel=200,col=topo.colors(200))
    lines(countingDomain,col=6,lwd=3)

    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)
    legend("bottomright",c("Observed","y=0","y>0"),col=c("white","black","black"),pch=c(0,0,15),bg="white")
    legend("bottomleft",c("Unmodelled","y=0","y>0"),col=c("white","red","red"),pch=c(0,0,15),bg="white")
    legend("topleft",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")


    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=fixed.effects.mesh,main="Mean of fixed effects (log-scale)",nlevel=200,col=topo.colors(200),ylab="")
    lines(countingDomain,col=6,lwd=3)
    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_mesh_without_obs.pdf"),width=14,height=6)
    par(mfrow=c(1,2),oma=c(0,0,0,1.5))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=covNewGridval,main="Covariate values",nlevel=200)
    lines(countingDomain,col=6,lwd=3)

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=fixed.effects.mesh,main="Mean of fixed effects (log-scale)",nlevel=200,ylab="")
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_count_domain_with_obs.pdf"),width=14,height=6)
    par(mfrow=c(1,2),oma=c(0,0,0,1.5))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=covNewGridvalDomain,main="Covariate values",nlevel=200,col=topo.colors(200))
    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)
    legend("bottomright",c("Observed","y=0","y>0"),col=c("white","black","black"),pch=c(0,0,15),bg="white")
    legend("bottomleft",c("Unmodelled","y=0","y>0"),col=c("white","red","red"),pch=c(0,0,15),bg="white")


    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=fixed.effects.domain,main="Mean of fixed effects (log-scale)",nlevel=200,col=topo.colors(200),ylab="")
    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_count_domain_without_obs.pdf"),width=14,height=6)
    par(mfrow=c(1,2),oma=c(0,0,0,1.5))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=covNewGridvalDomain,main="Covariate values",nlevel=200)

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=fixed.effects.domain,main="Mean of fixed effects (log-scale)",nlevel=200,ylab="")

    dev.off()
  }

  ## Producing an 2x2 summary plot

  if (summaryplot){
    pdf(file=file.path(savingFolder,"results_with_data.pdf"),width=12,height=12)
    par(mfrow=c(2,2))

    ## Mean posterior field
    if (sum(leaveOutTransect)>0){
      fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$mean.field,main="Mean of latent field (log-scale)",nlevel=200,col=topo.colors(200))
    } else {
      fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$mean.field.samp,main="Mean of latent field (log-scale)",nlevel=200,col=topo.colors(200))
    }

    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)
    lines(countingDomain,col=6,lwd=3)
    legend("topleft",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")
    legend("bottomright",c("Observed","y=0","y>0"),col=c("white","black","black"),pch=c(0,0,15),bg="white")
    legend("bottomleft",c("Unmodelled","y=0","y>0"),col=c("white","red","red"),pch=c(0,0,15),bg="white")




    ## Sd of posterior field
    if (sum(leaveOutTransect)>0){
      frame()
    } else {
      fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$sd.field.samp,main="Sd of latent field (log-scale)",nlevel=200,col=topo.colors(200))
      rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
      rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
      rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
      rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)
      lines(countingDomain,col=6,lwd=3)

    }



    ## Posterior predictive dist
    if (sum(leaveOutTransect)>0){
      plotTo <- which.min((cumsum(finalResList$posteriorDistTransect)-0.99)^2)
      plot(finalResList$posteriorevalFullTransect[1:plotTo],finalResList$posteriorDistTransect[1:plotTo],main="Posterior Predictive dist FOR PHOTOS IN TRANSECT",type='h',xlab="seals",ylab="probability",col="grey") # Now uses the unsmoothed version for plotting instead
    } else {
      plotTo <- which.min((cumsum(finalResList$posteriorDist)-0.99)^2)
      plot(finalResList$posteriorEvalPoints[1:plotTo],finalResList$posteriorDist[1:plotTo],main="Posterior Predictive dist",type='h',xlab="seals",ylab="probability",col="grey") # Now uses the unsmoothed version for plotting instead
    }
    lines(rep(finalResList$posteriorMean,2),c(0,1000),col=2)
    lines(rep(finalResList$posteriorMedian,2),c(0,1000),col=3)
    lines(rep(finalResList$posteriorMode,2),c(0,1000),col=4)

    legend("topright",c(paste("mean =",round(finalResList$posteriorMean,2)),
                        paste("median =",round(finalResList$posteriorMedian,2)),
                        paste("mode =",round(finalResList$posteriorMode,2)),
                        paste("IQR =",round(finalResList$posteriorIQR,2)),
                        paste(round(results.CI.level*100),"% CI = (",round(finalResList$posteriorCI[1]),",",round(finalResList$posteriorCI[2]),")"),
                        paste("range = ", round(finalResList$mean.range.param,3))),
           lty=1,col=c(2:4,rep("white",3)))


    ## Just some parameters and variables

    # First splitting comment if it is too long:

    maxChar <- 50
    if (nchar(comment)>maxChar){
      splits <- gregexpr(pattern="=",comment)[[1]]
      splitHere <- max(splits[splits<maxChar])
      comment <- paste(substr(comment, 1, splitHere-1), "\n", substr(comment, splitHere, nchar(comment)), sep = "")
    }

    stStart <- sort(gregexpr('/',sealPhotoDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealPhotoDataFile)
    photoFrom <- substr(sealPhotoDataFile,stStart,stStop)

    stStart <- sort(gregexpr('/',sealTransectDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealTransectDataFile)
    transectsFrom <- substr(sealTransectDataFile,stStart,stStop)


    par(xpd=TRUE)
    frame()
    #if (dataType=="simulated"){
    #  text(0.5,1.15,paste("True # seals = ",sampPoisCounted$n,sep=""))
    #}
    fixed.effects.vec <- pp.res$summary.fixed$mean
    text(0.5,1.10,paste("Photos and transects from = ",photoFrom," and ",transectsFrom,sep=""))
    text(0.5,1.05,paste("SealType = ",sealType,sep=""))
    text(0.5,1.00,paste("use.covariates = ",use.covariates,sep=""))
    text(0.5,0.95,paste("additional.iid.term = ",additional.iid.term,sep=""))
    text(0.5,0.90,paste("covariates.type = ",covariates.type,sep=""))
    text(0.5,0.85,paste("spatial = ",spatial,sep=""))
    text(0.5,0.80,paste("covariate.fitting = ",covariate.fitting,sep=""))
    text(0.5,0.75,paste("convHullVar.convex = ",convHullVar.convex,sep=""))
    text(0.5,0.70,paste("convHullVar.concave = ",convHullVar.concave,sep=""))
    text(0.5,0.65,paste("convHullVar.resolution = ",paste(convHullVar.resolution,collapse=", "),sep=""))
    text(0.5,0.60,paste("meshVar.max.edge = ",paste(meshVar.max.edge,collapse=", "),sep=""))
    text(0.5,0.55,paste("meshVar.offset = ",meshVar.offset,sep=""))
    text(0.5,0.50,paste("meshVar.cutoff = ",paste(meshVar.cutoff,collapse=", "),sep=""))
    text(0.5,0.45,paste("Matern.alpha =",Matern.alpha,sep=""))
    text(0.5,0.40,paste("grid.pixelsize =",grid.pixelsize,sep=""))
    text(0.5,0.35,paste("INLA.theta.startval: ",paste(INLA.theta.startval,collapse=", "),sep=""))
    text(0.5,0.30,paste("parallelize.numCores",parallelize.numCores,sep=""))
    text(0.5,0.25,paste("poisson.maxEvals = ",poisson.maxEvals,sep=""))
    text(0.5,0.20,paste("Number of posterior samples = ",noSamp,sep=""))
    text(0.5,0.15,paste("Mean of fixed effects = ", paste(round(fixed.effects.vec,4),collapse=", "),sep=""))
    text(0.5,0.10,paste("DIC = ", round(finalResList$dic,4),sep=""))
    text(0.5,0.05,paste("WAIC = ", round(finalResList$waic,4),sep=""))
    text(0.5,0.00,paste("Marginal log-likelihood = ", round(finalResList$mlik,4),sep=""))
    text(0.5,-0.05,paste("Running time = ", round(time[3]/60,2), " minutes", sep=""))
    text(0.5,-0.10,paste("testing: ",testing,sep=""))
    text(0.5,-0.16,paste("comment: ",comment,sep=""))
    dev.off()



    pdf(file=file.path(savingFolder,"results.pdf"),width=12,height=12)
    par(mfrow=c(2,2))

    ## Mean posterior field
    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$mean.field.domain,main="Mean of latent field (log-scale)",nlevel=200)

    ## Sd of posterior field
    if (sum(leaveOutTransect)>0){
      frame()
    } else {
      fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$sd.field.domain.samp,main="Sd of latent field (log-scale)",nlevel=200)
    }

    ## Posterior predictive dist
    if (sum(leaveOutTransect)>0){
      plotTo <- which.min((cumsum(finalResList$posteriorDistTransect)-0.99)^2)
      plot(finalResList$posteriorevalFullTransect[1:plotTo],finalResList$posteriorDistTransect[1:plotTo],main="Posterior Predictive dist FOR PHOTOS IN TRANSECT",type='h',xlab="seals",ylab="probability",col="grey") # Now uses the unsmoothed version for plotting instead
    } else {
      plotTo <- which.min((cumsum(finalResList$posteriorDist)-0.99)^2)
      plot(finalResList$posteriorEvalPoints[1:plotTo],finalResList$posteriorDist[1:plotTo],main="Posterior Predictive dist",type='h',xlab="seals",ylab="probability",col="grey") # Now uses the unsmoothed version for plotting instead
    }
    lines(rep(finalResList$posteriorMean,2),c(0,1000),col=2)
    lines(rep(finalResList$posteriorMedian,2),c(0,1000),col=3)
    lines(rep(finalResList$posteriorMode,2),c(0,1000),col=4)

    legend("topright",c(paste("mean =",round(finalResList$posteriorMean,2)),
                        paste("median =",round(finalResList$posteriorMedian,2)),
                        paste("mode =",round(finalResList$posteriorMode,2)),
                        paste("IQR =",round(finalResList$posteriorIQR,2)),
                        paste(round(results.CI.level*100),"% CI = (",round(finalResList$posteriorCI[1]),",",round(finalResList$posteriorCI[2]),")"),
                        paste("range = ", round(finalResList$mean.range.param,3))),
           lty=1,col=c(2:4,rep("white",3)))

    ## Just some parameters and variables

    # First splitting comment if it is too long:

    maxChar <- 50
    if (nchar(comment)>maxChar){
      splits <- gregexpr(pattern="=",comment)[[1]]
      splitHere <- max(splits[splits<maxChar])
      comment <- paste(substr(comment, 1, splitHere-1), "\n", substr(comment, splitHere, nchar(comment)), sep = "")
    }

    stStart <- sort(gregexpr('/',sealPhotoDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealPhotoDataFile)
    photoFrom <- substr(sealPhotoDataFile,stStart,stStop)

    stStart <- sort(gregexpr('/',sealTransectDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealTransectDataFile)
    transectsFrom <- substr(sealTransectDataFile,stStart,stStop)


    par(xpd=TRUE)
    frame()
    #if (dataType=="simulated"){
    #  text(0.5,1.15,paste("True # seals = ",sampPoisCounted$n,sep=""))
    #}
    fixed.effects.vec <- pp.res$summary.fixed$mean
    text(0.5,1.10,paste("Photos and transects from = ",photoFrom," and ",transectsFrom,sep=""))
    text(0.5,1.05,paste("SealType = ",sealType,sep=""))
    text(0.5,1.00,paste("use.covariates = ",use.covariates,sep=""))
    text(0.5,0.95,paste("additional.iid.term = ",additional.iid.term,sep=""))
    text(0.5,0.90,paste("covariates.type = ",covariates.type,sep=""))
    text(0.5,0.85,paste("spatial = ",spatial,sep=""))
    text(0.5,0.80,paste("covariate.fitting = ",covariate.fitting,sep=""))
    text(0.5,0.75,paste("convHullVar.convex = ",convHullVar.convex,sep=""))
    text(0.5,0.70,paste("convHullVar.concave = ",convHullVar.concave,sep=""))
    text(0.5,0.65,paste("convHullVar.resolution = ",paste(convHullVar.resolution,collapse=", "),sep=""))
    text(0.5,0.60,paste("meshVar.max.edge = ",paste(meshVar.max.edge,collapse=", "),sep=""))
    text(0.5,0.55,paste("meshVar.offset = ",meshVar.offset,sep=""))
    text(0.5,0.50,paste("meshVar.cutoff = ",paste(meshVar.cutoff,collapse=", "),sep=""))
    text(0.5,0.45,paste("Matern.alpha =",Matern.alpha,sep=""))
    text(0.5,0.40,paste("grid.pixelsize =",grid.pixelsize,sep=""))
    text(0.5,0.35,paste("INLA.theta.startval: ",paste(INLA.theta.startval,collapse=", "),sep=""))
    text(0.5,0.30,paste("parallelize.numCores",parallelize.numCores,sep=""))
    text(0.5,0.25,paste("poisson.maxEvals = ",poisson.maxEvals,sep=""))
    text(0.5,0.20,paste("Number of posterior samples = ",noSamp,sep=""))
    text(0.5,0.15,paste("Mean of fixed effects = ", paste(round(fixed.effects.vec,4),collapse=", "),sep=""))
    text(0.5,0.10,paste("DIC = ", round(finalResList$dic,4),sep=""))
    text(0.5,0.05,paste("WAIC = ", round(finalResList$waic,4),sep=""))
    text(0.5,0.00,paste("Marginal log-likelihood = ", round(finalResList$mlik,4),sep=""))
    text(0.5,-0.05,paste("Running time = ", round(time[3]/60,2), " minutes", sep=""))
    text(0.5,-0.10,paste("testing: ",testing,sep=""))
    text(0.5,-0.16,paste("comment: ",comment,sep=""))
    dev.off()


  }

  print("All plotting to file completed")

}



#' Summary plot for continuous likelihood approximation approach
#'
#' Produces several different summary plots showing the results from INLA approach
#'
#' @param meshplot Logical, indicating whether a plot with the mesh should be produced
#' @param boundariesplot Logical, indicating whether plots with the used boundaries and domains should be produced
#' @param covariatesplot Logical, indicating whether various plots showing the covariates and their effects should be produced
#' @param summaryplot Logical, indicating whether a plot showing the final posterior distribution and mean and pointwise sd of the latent field should be produced
#' @param savingFolder String, indicating the complete path to where the plots are to be stored
#' @param sealPhotoDataFile String, indicating the file where seal photo data are stored
#' @param sealTransectDataFile String, indicating the file where seal transect data are stored
#' @param gridList List containing information about the grid, being the output of the function GridCreation
#' @param dataList List with various data variables
#' @param finalResList List gathering all final results
#' @param mesh Mesh object, being the output from inla.mesh.2d
#' @param countingDomain data.frame containing x- and y-coordinates of the counting domain polygon
#' @param logicalGridPointsInsideCountingDomain Logical vector, indicating which of the grid elements are within the counting domain
#' @param covNewGridval Matrix with all covariate values at grid points
#' @param fixed.effects.vec Numeric vector with the estimated mean fixed effects
#' @param sealType String, indicating which seal type to model "hooded" (default as it is quicker), or "harps"
#' @param use.covariates Logical, indicating whether covariates are used or not (see description!)
#' @param covariate.fitting String, indicating how to model covariates. "linear", quadratic (default) or "linAndLog", or FALSE for no covariates
#' @param covariates.type String equal to "band1" or "band2" indicating which of the bands from the satellite data should be used.
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param convHullVar.convex Numeric, corrresponding to the convex parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.concave Numeric, corrresponding to the concave parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.resolution Numeric vector of length 2, corrresponding to the resolution parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param meshVar.max.edge Numeric, corrresponding to the max.edge parameter of inla.mesh.2d. Smaller values gives smaller triangles outside triangles. See ?inla.mesh.2d for details
#' @param meshVar.offset Numeric, corrresponding to the offset parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param meshVar.cutoff Numeric, corrresponding to the cutoff parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param Matern.alpha Numeric, corresponding to the alpha parameter in the Matern covariance function (2 is the default)
#' @param grid.pixelsize Numeric, denoting the size (in km, in both x- and y-direction) of each pixel of the grid being used
#' @param INLA.theta.startval List containing the start values for the theta parameter in the INLA run. (NULL indicates that automatic start values should be used, and is the default)
#' @param parallelize.numCores Numeric, corresponding to the number of cores any parallelization should be run at
#' @param parallelize.noSplits Numeric, deciding how many sublists samp should be splitted into. Should be a multiple of parallelize.numCores for the highest efficiency. The larger number the less memory is used (and longer time).
#' @param poisson.maxEvals Numeric, corresponding to maximum number of points the Poisson distribution should be evaluated at (a much smaller number is typically used)
#' @param noSamp Numeric, indicating the number of samples from posterior model to use when computing the posterior predictive distribution. 5000 (default) typically sufficient.
#' @param time The number of seconds the function has been running
#' @param testing Logical, indicating whether the testing parts of this function should be used
#' @param comment String to add as a comment to the summary plot
#' @return Produces a set of plots and writes them as pdf to file. Nothing else is returned
#' @keywords plot
#' @export


SummaryPlotFuncContLikApprox <- function(meshplot = TRUE,
                                         boundariesplot = TRUE,
                                         covariatesplot = TRUE,
                                         summaryplot = TRUE,
                                         savingFolder,
                                         sealPhotoDataFile,
                                         sealTransectDataFile,
                                         dataList,
                                         results.CI.level = 0.95,
                                         observationDomain,
                                         gridList,
                                         finalResList,
                                         mesh,
                                         countingDomain,
                                         logicalGridPointsInsideCountingDomain,
                                         covNewGridval,
                                         fixed.effects.vec,
                                         sealType,
                                         use.covariates,
                                         covariates.type,
                                         covariate.fitting,
                                         additional.iid.term,
                                         convHullVar.convex,
                                         convHullVar.concave,
                                         convHullVar.resolution,
                                         meshVar.max.edge,
                                         meshVar.offset,
                                         meshVar.cutoff,
                                         Matern.alpha,
                                         grid.pixelsize,
                                         INLA.theta.startval,
                                         parallelize.noSplits,
                                         parallelize.numCores,
                                         poisson.maxEvals,
                                         noSamp,
                                         time,
                                         testing,
                                         comment
){

  ## Creating some helping variables
  dataList$obsPoints
  fixed.effects.mesh <- finalResList$intercept.mean + finalResList$covariate.mean*covNewGridval + finalResList$covariate2.mean*covNewGridval^2 + finalResList$covariateLog.mean*log(covNewGridval)
  fixed.effects.domain <- fixed.effects.mesh
  fixed.effects.domain[!logicalGridPointsInsideCountingDomain] = NA

  covNewGridvalDomain <- covNewGridval
  covNewGridvalDomain[!logicalGridPointsInsideCountingDomain] = NA


  ## Producing a mesh plot

  if (meshplot){
    ## Plotting the mesh and defined boundaries
    pdf(file=file.path(savingFolder,"mesh.pdf"),width=6,height=6)
    par(mar=c(0,0,0,0))
    plot(mesh, asp=1, main='')
    points(dataList$coordPhotoX[dataList$noObsPerPhoto==0],dataList$coordPhotoY[dataList$noObsPerPhoto==0],col="red",cex=0.5,pch=16)
    points(dataList$coordPhotoX[dataList$noObsPerPhoto>0],dataList$coordPhotoY[dataList$noObsPerPhoto>0],col="green",cex=0.5,pch=16)
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("y=0","y>0"),col=c("red","green"),pch=16)
    legend("topleft",c("Model boundary","counting domain"),col=c("blue",6),lty=c(1),lwd=3)
    dev.off()
  }

  ## Producing 2 x 2 boundaries plot

  if (boundariesplot){ ### New type of plot here!

    pdf(file=file.path(savingFolder,"boundaries.pdf"),width=12,height=12)
    par(mfrow=c(2,2))

    ## Plot of the tiles and give their corresponding obesrvation value  -- complete area
    plot(mesh$loc,type='n',main="Observations",xlab="x", ylab='y')
    rect(xleft = dataList$coordPhotoX-dataList$photoWidth*0.5,
         xright = dataList$coordPhotoX+dataList$photoWidth*0.5,
         ybottom = dataList$coordPhotoY-dataList$photoHeight*0.5,
         ytop = dataList$coordPhotoY+dataList$photoHeight*0.5,
         col = c("red","green")[(dataList$noObsPerPhoto>0)*1+1],
         border=c("red","green")[(dataList$noObsPerPhoto>0)*1+1])

    legend("bottomright",c("y=0","y>0"),col=c("red","green"),pch=16,bg="white")

    ## Plot of the tiles and give their corresponding weight value ( e ) -- complete area
    plot(mesh$loc,type='n',main="Weights/offset ( e )",xlab="x", ylab='y')
    points(dataList$intPoints$x, dataList$intPoints$y,col="purple",cex=0.75,pch=16)
    points(dataList$obsPoints$x, dataList$obsPoints$y,col="orange",cex=0.75,pch=16)

    legend("bottomright",c("e>0 & observed","e>0 & intPoint"),col=c("orange","purple"),pch=16,bg="white")

    ## Plot of the tiles and give their corresponding obesrvation value  -- zooms in on speciic area

    par(mar=c(2,2,1,1), mgp=2:0)
    plot(mesh,xlim=c(25,57), ylim=c(-20,5)-3,main="Observations",xlab="x", ylab='y')
    rect(xleft = dataList$coordPhotoX-dataList$photoWidth*0.5,
         xright = dataList$coordPhotoX+dataList$photoWidth*0.5,
         ybottom = dataList$coordPhotoY-dataList$photoHeight*0.5,
         ytop = dataList$coordPhotoY+dataList$photoHeight*0.5,
         col = c("red","green")[(dataList$noObsPerPhoto>0)*1+1],
         border=c("red","green")[(dataList$noObsPerPhoto>0)*1+1])
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("y=NA","y=0","y>0"),col=c("white","red","green"),pch=15,bg="white")
    legend("bottomleft",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")


    ## Plot of the tiles and give their corresponding weight value ( e ) -- zooms in on specific area
    plot(mesh,main="Weights/offset ( e )",xlim=c(25,57), ylim=c(-23,2),xlab="x", ylab='y')
    for (i in unique(observationDomain$SID)){
      theseObs <- observationDomain[observationDomain$SID==i,]
      n.theseObs <- nrow(theseObs)
      lines(theseObs$X[c(1:n.theseObs,1)],theseObs$Y[c(1:n.theseObs,1)],col='green')
    }
    points(dataList$intPoints$x, dataList$intPoints$y,col="purple",cex=1,pch=16)
    points(dataList$obsPoints$x, dataList$obsPoints$y,col="orange",cex=1,pch=16)
    lines(countingDomain,col=6,lwd=3)
    legend("bottomleft",c("Observation domain","Count domain"),col=c("green",6),lty=c(1,1),lwd=c(1,3),bg="white")
    legend("bottomright",c("e>0 & observed","e>0 & intPoint"),col=c("orange","purple"),pch=16,bg="white")

    dev.off()

  }

  ## Producing 4 1x2 covariate plots

  if (covariatesplot){
    pdf(file=file.path(savingFolder,"covariate_effects_mesh_with_obs.pdf"),width=13,height=6)
    par(mfrow=c(1,2))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=covNewGridval,main="Covariate values",nlevel=200,col=topo.colors(200))
    lines(countingDomain,col=6,lwd=3)
    points(dataList$coordPhotoX,dataList$coordPhotoY,cex=0.5,col=1)
    points(dataList$obsPoints[,1:2],cex=1,col="red")

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=fixed.effects.mesh,main="Mean of fixed effects (log-scale)",nlevel=200,col=topo.colors(200))
    lines(countingDomain,col=6,lwd=3)
    points(dataList$coordPhotoX,dataList$coordPhotoY,cex=0.5,col=1)
    points(dataList$obsPoints[,1:2],cex=1,col="red")
    legend("bottomright",c("y=0","y>0"),col=c("black","red"),pch=c(15,1),bg="white")
    legend("topleft",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_mesh_without_obs.pdf"),width=13,height=6)
    par(mfrow=c(1,2))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=covNewGridval,main="Covariate values",nlevel=200)
    lines(countingDomain,col=6,lwd=3)

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=fixed.effects.mesh,main="Mean of fixed effects (log-scale)",nlevel=200)
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_count_domain_with_obs.pdf"),width=13,height=6)
    par(mfrow=c(1,2))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=covNewGridvalDomain,main="Covariate values",nlevel=200,col=topo.colors(200))
    points(dataList$coordPhotoX,dataList$coordPhotoY,cex=0.5,col=1)
    points(dataList$obsPoints[,1:2],cex=1,col="red")

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=fixed.effects.domain,main="Mean of fixed effects (log-scale)",nlevel=200,col=topo.colors(200))
    points(dataList$coordPhotoX,dataList$coordPhotoY,cex=0.5,col=1)
    points(dataList$obsPoints[,1:2],cex=1,col="red")
    legend("bottomright",c("y=0","y>0"),col=c("black","red"),pch=c(15,1),bg="white")

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_count_domain_without_obs.pdf"),width=13,height=6)
    par(mfrow=c(1,2))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=covNewGridvalDomain,main="Covariate values",nlevel=200)

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=fixed.effects.domain,main="Mean of fixed effects (log-scale)",nlevel=200)

    dev.off()
  }

  ## Producing an 2x2 summary plot

  if (summaryplot){
    pdf(file=file.path(savingFolder,"results.pdf"),width=12,height=12)
    par(mfrow=c(2,2))

    ## Mean posterior field
    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$mean.field,main="Mean of latent field (log-scale)",nlevel=200)
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")

    ## Sd of posterior field
    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$sd.field.samp,main="Sd of latent field (log-scale)",nlevel=200)
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")


    ## Posterior predictive dist
    plot(finalResList$posteriorEvalPoints,finalResList$posteriorDist,main="Posterior Predictive dist",type='h',col="grey",xlab="seals",ylab="probability") # Now uses the unsmoothed version for plotting instead
    lines(rep(finalResList$posteriorMean,2),c(0,1000),col=2)
    lines(rep(finalResList$posteriorMedian,2),c(0,1000),col=3)
    lines(rep(finalResList$posteriorMode,2),c(0,1000),col=4)

    legend("topright",c(paste("mean =",round(finalResList$posteriorMean,2)),
                        paste("median =",round(finalResList$posteriorMedian,2)),
                        paste("mode =",round(finalResList$posteriorMode,2)),
                        paste("IQR =",round(finalResList$posteriorIQR,2)),
                        paste(round(results.CI.level*100),"% CI = (",round(finalResList$posteriorCI[1]),",",round(finalResList$posteriorCI[2]),")"),
                        paste("range = ", round(finalResList$mean.range.param,3))),
           lty=1,col=c(2:4,rep("white",3)))

    ## Just some parameters and variables

    # First splitting comment if it is too long:

    maxChar <- 50
    if (nchar(comment)>maxChar){
      splits <- gregexpr(pattern="=",comment)[[1]]
      splitHere <- max(splits[splits<maxChar])
      comment <- paste(substr(comment, 1, splitHere-1), "\n", substr(comment, splitHere, nchar(comment)), sep = "")
    }

    stStart <- sort(gregexpr('/',sealPhotoDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealPhotoDataFile)
    photoFrom <- substr(sealPhotoDataFile,stStart,stStop)

    stStart <- sort(gregexpr('/',sealTransectDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealTransectDataFile)
    transectsFrom <- substr(sealTransectDataFile,stStart,stStop)


    par(xpd=TRUE)
    frame()
    #if (dataType=="simulated"){
    #  text(0.5,1.15,paste("True # seals = ",sampPoisCounted$n,sep=""))
    #}

    text(0.5,1.10,paste("Photos and transects from = ",photoFrom," and ",transectsFrom,sep=""))
    text(0.5,1.05,paste("SealType = ",sealType,sep=""))
    text(0.5,1.00,paste("use.covariates = ",use.covariates,sep=""))
    text(0.5,0.95,paste("additional.iid.term = ",additional.iid.term,sep=""))
    text(0.5,0.90,paste("covariates.type = ",covariates.type,sep=""))
    text(0.5,0.85,paste("covariate.fitting = ",covariate.fitting,sep=""))
    text(0.5,0.80,paste("convHullVar.convex = ",convHullVar.convex,sep=""))
    text(0.5,0.75,paste("convHullVar.concave = ",convHullVar.concave,sep=""))
    text(0.5,0.70,paste("convHullVar.resolution = ",paste(convHullVar.resolution,collapse=", "),sep=""))
    text(0.5,0.65,paste("meshVar.max.edge = ",paste(meshVar.max.edge,collapse=", "),sep=""))
    text(0.5,0.60,paste("meshVar.offset = ",meshVar.offset,sep=""))
    text(0.5,0.55,paste("meshVar.cutoff =",meshVar.cutoff,sep=""))
    text(0.5,0.50,paste("Matern.alpha =",Matern.alpha,sep=""))
    text(0.5,0.45,paste("grid.pixelsize =",grid.pixelsize,sep=""))
    text(0.5,0.40,paste("INLA.theta.startval: ",paste(INLA.theta.startval,collapse=", "),sep=""))
    text(0.5,0.35,paste("parallelize.numCores",parallelize.numCores,sep=""))
    text(0.5,0.30,paste("poisson.maxEvals = ",poisson.maxEvals,sep=""))
    text(0.5,0.25,paste("Number of posterior samples = ",noSamp,sep=""))
    text(0.5,0.20,paste("Mean of fixed effects = ", paste(round(fixed.effects.vec,4),collapse=", "),sep=""))
    text(0.5,0.15,paste("DIC = ", round(finalResList$dic,4),sep=""))
    text(0.5,0.10,paste("WAIC = ", round(finalResList$waic,4),sep=""))
    text(0.5,0.05,paste("Marginal log-likelihood = ", round(finalResList$mlik,4),sep=""))
    text(0.5,0.00,paste("Running time = ", round(time[3]/60,2), " minutes", sep=""))
    text(0.5,-0.05,paste("testing: ",testing,sep=""))
    text(0.5,-0.10,paste("comment: ",comment,sep=""))

    dev.off()
  }

  print("All plotting to file completed")

}


#' Posterior summary statistics for GAM
#'
#' Computes basic summary statistics based on distribution
#'
#' @param posteriorevalFullPhoto Numeric vector of evaluation points used for all posteriorDistributions
#' @param posteriorDistPhoto Matrix where each column is the posterior distribution probability for one photo
#' @param posteriorevalFullTransect Numeric vector of evaluation points used for posteriorDistTransect
#' @param posteriorDistTransect Numeric vector of posterior distribution probabilitites for the whole transect in question
#' @param areaPerSubSampPhotoMatrix Numeric matrix with
#' @param photoCounts Numeric vector with the observed photo counts (not used in the modelling)
#' @param photoOrder Numeric vector specifying the plotting order of the photos (smallest to largest original x-coordinate)
#' @param leaveOutTransect The transect in question
#' @param inputVar List of all input variables used in the original run
#' @param savingFolder String with the path for where to save these comparison results
#' @return Nothing, just save the output to file
#' @import Hmisc
#' @export

ComparePhotoCountsAndPred <- function(posteriorevalFullPhoto,
                                      posteriorDistPhoto,
                                      posteriorevalFullTransect,
                                      posteriorDistTransect,
                                      areaPerSubSampPhotoMatrix,
                                      photoCounts,
                                      photoOrder,
                                      leaveOutTransect,
                                      inputVar,
                                      savingFolder,
                                      plot = TRUE){

  counts <- 0:(max(photoCounts,posteriorevalFullPhoto)+1)
  posteriorDistPhotoMat <- FhatPhotoMat <- matrix(0,ncol=length(counts),nrow=length(photoCounts))
  for (i in 1:length(photoCounts)){
    theseEvals <- which(counts %in% posteriorevalFullPhoto)
    posteriorDistPhotoMat[i,theseEvals] <- posteriorDistPhoto[,i]
    FhatPhotoMat[i,] <- cumsum(posteriorDistPhotoMat[i,])
  }


  sampWithinMat <- matrix(NA,ncol=2,nrow=length(photoCounts))
  sampVec <- rep(NA,length(photoCounts))
  for (i in 1:length(photoCounts)){
    evalPoint <- which(counts==photoCounts[i])
    FhatPhotoVec <- c(0,FhatPhotoMat[i,])
    sampWithinMat[i,] <- c(FhatPhotoVec[evalPoint],FhatPhotoVec[evalPoint+1])
    sampVec[i] <- runif(n=1,min = sampWithinMat[i,1],max= sampWithinMat[i,2])
  }

  meanVec <- as.vector(t(counts) %*% t(posteriorDistPhotoMat))

  quant.finder <- function(vec,q,evalvec){
    vec.q = vec-q
    evalvec[vec.q>0][1]
  }
  quantMat <- matrix(NA,ncol=5,nrow=length(photoCounts))
  colnames(quantMat) = c("quant005","quant025","median","quant075","quant095")
  for (i in 1:length(photoCounts)){
    quantMat[i,1] <- quant.finder(vec=FhatPhotoMat[i,],q=0.05,evalvec = counts)
    quantMat[i,2] <- quant.finder(vec=FhatPhotoMat[i,],q=0.25,evalvec = counts)
    quantMat[i,3] <- quant.finder(vec=FhatPhotoMat[i,],q=0.5,evalvec = counts)
    quantMat[i,4] <- quant.finder(vec=FhatPhotoMat[i,],q=0.75,evalvec = counts)
    quantMat[i,5] <- quant.finder(vec=FhatPhotoMat[i,],q=0.95,evalvec = counts)
  }
  predPhotoDf <- as.data.frame(quantMat)
  predPhotoDf$mean <- meanVec
  predPhotoDf$photoCounts <- photoCounts
  predPhotoDf$photoOrder <- photoOrder


  countsTrans <- 0:(max(posteriorevalFullTransect)+1)
  posteriorDistTransectNew <- rep(0,length(countsTrans))
  theseEvals <- which(countsTrans %in% posteriorevalFullTransect)
  posteriorDistTransectNew[theseEvals] <- posteriorDistTransect
  FhatTrans <- cumsum(posteriorDistTransectNew)

  meanTrans <- c(t(countsTrans)%*%posteriorDistTransectNew)
  medTrans <- quant.finder(vec=FhatTrans,q=0.5,evalvec = countsTrans)

  quant005Trans <- quant.finder(vec=FhatTrans,q=0.05,evalvec = countsTrans)
  quant025Trans <- quant.finder(vec=FhatTrans,q=0.25,evalvec = countsTrans)
  quant075Trans <- quant.finder(vec=FhatTrans,q=0.75,evalvec = countsTrans)
  quant095Trans <- quant.finder(vec=FhatTrans,q=0.95,evalvec = countsTrans)

  trans.list <- list(posteriorevalFullTransect=countsTrans,
                     posteriorDistTransect = posteriorDistTransectNew,
                     FhatTrans = FhatTrans,
                     meanTrans = meanTrans,
                     medTrans = medTrans,
                     quant095Trans = quant095Trans,
                     quant005Trans = quant005Trans,
                     quant075Trans = quant075Trans,
                     quant025Trans = quant025Trans)

  #################

  posteriorFieldAtPhoto = data.frame(median = apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=median))

  posteriorFieldAtPhoto$mean <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=mean)
  posteriorFieldAtPhoto$quant005 <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=quantile, probs=0.05)
  posteriorFieldAtPhoto$quant025 <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=quantile, probs=0.25)
  posteriorFieldAtPhoto$quant075 <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=quantile, probs=0.75)
  posteriorFieldAtPhoto$quant095 <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=quantile, probs=0.95)
  posteriorFieldAtPhoto$photoCounts <- photoCounts

  PropTrueCountIn050CIMEAN <- mean(posteriorFieldAtPhoto$quant025 <= posteriorFieldAtPhoto$photoCounts & posteriorFieldAtPhoto$quant075 >= posteriorFieldAtPhoto$photoCounts)
  PropTrueCountIn090CIMEAN <- mean(posteriorFieldAtPhoto$quant005 <= posteriorFieldAtPhoto$photoCounts & posteriorFieldAtPhoto$quant095 >= posteriorFieldAtPhoto$photoCounts)

  (rmse.meanPredMEAN <- sqrt(mean((posteriorFieldAtPhoto$photoCounts-posteriorFieldAtPhoto$mean)^2)))
  (mae.medianPredMEAN <- mean(abs(posteriorFieldAtPhoto$photoCounts-posteriorFieldAtPhoto$median)))




  if (plot){

  pdf(file=file.path(savingFolder,"photo_counts_and_predictions.pdf"),width=7,height=4)
  par(mar=c(3,3,2,2),mgp=2:0,cex=0.8)

  plot(1:length(predPhotoDf$photoCounts),predPhotoDf$median[photoOrder],type='n',ylim=c(0,max(c(predPhotoDf$quant095,predPhotoDf$photoCounts))),
       main = paste("Photo counts and prediction comparison for transect/CVFold ",paste(leaveOutTransect,collapse = "&"),sep=""),ylab="Predicted/observed counts",xlab="Photo in transect (left to right)")

  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$quant005[photoOrder],lty=1,lwd=1,col=2)
  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$quant095[photoOrder],lty=1,lwd=1,col=2)

  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$quant025[photoOrder],lty=1,lwd=1,col=4)
  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$quant075[photoOrder],lty=1,lwd=1,col=4)


  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$median[photoOrder],col=1,lwd=2)
  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$photoCounts[photoOrder],col=3,lty=2,lwd=1.5)

  legend("topleft",c("Post.pred. median","Post.pred 90% CI","Post.pred 50% CI","Observed counts"),col=c(1,2,4,3),lty=c(1,1,1,2),lwd=c(2,1,1,1.5))

  dev.off()


  PropTrueCountIn050CI <- mean(predPhotoDf$quant025 <= predPhotoDf$photoCounts & predPhotoDf$quant075 >= predPhotoDf$photoCounts)
  PropTrueCountIn090CI <- mean(predPhotoDf$quant005 <= predPhotoDf$photoCounts & predPhotoDf$quant095 >= predPhotoDf$photoCounts)

  (rmse.meanPred <- sqrt(mean((predPhotoDf$photoCounts-predPhotoDf$mean)^2)))
  (mae.medianPred <- mean(abs(predPhotoDf$photoCounts-predPhotoDf$median)))
  photoVec <- 1:length(predPhotoDf$photoCounts)


  pdf(file=file.path(savingFolder,"photo_counts_and_predictions2.pdf"),width=10,heigh=5)
  par(mfrow=c(1,1))
  plot(photoVec,predPhotoDf$photoCounts[photoOrder],type='n',ylim=c(0,max(predPhotoDf$quant095,predPhotoDf$photoCounts)),lwd=2,
       main="REG\n Comparison of true counts and out-of-sample posterior predictions per photo in transect/CVFold",
       xlab="Photo number (from left)",
       ylab="Seal count")
  k=0.25
  for(i in 1:length(photoVec))
  {
    a <- photoVec[photoOrder][i]
    b <- predPhotoDf$quant095[photoOrder][i]
    c <- predPhotoDf$quant005[photoOrder][i]
    polygon(c(a-k, a+k, a+k, a-k), c(c, c, b, b), col="gray90", border="gray90")
    b <- predPhotoDf$quant075[photoOrder][i]
    c <- predPhotoDf$quant025[photoOrder][i]
    polygon(c(a-k, a+k, a+k, a-k), c(c, c, b, b), col="gray70", border="gray70")

    d <- predPhotoDf$median[photoOrder][i]
    lines(c(a-k, a+k), c(d,d), lwd=2, col="red")

    e <- predPhotoDf$mean[photoOrder][i]
    lines(c(a-k, a+k), c(e,e), lwd=2, col="blue")

  }

  points(photoVec[photoOrder],predPhotoDf$photoCounts[photoOrder],lwd=1)
  legend("topright",c("True Count","Pred median","Pred mean","90 CI","50% CI"),col=c(1,"red","blue","grey90","grey70"),lty=c(1,1,1,1,1),lwd=c(NA,1,1,10,10),pch=c(1,NA,NA,NA,NA),bg="white")
  legend("topleft",c(paste("RMSE mean pred = ",round(rmse.meanPred,2),sep=""),
                     paste("MAE median pred = ",round(mae.medianPred,2),sep=""),
                     paste("Prop truth in 50% CI = ",round(PropTrueCountIn050CI,2),sep=""),
                     paste("Prop truth in 90% CI = ",round(PropTrueCountIn090CI,2),sep="")),
         pch=NA,col="white")
  dev.off()



  pdf(file=file.path(savingFolder,"photo_counts_and_posterior_mean.pdf"),width=10,heigh=5)
  par(mfrow=c(1,1))
  plot(photoVec,predPhotoDf$photoCounts[photoOrder],type='n',ylim=c(0,max(posteriorFieldAtPhoto$quant095,predPhotoDf$photoCounts)),lwd=2,
       main="REG\n Comparison of true counts and out-of-sample posterior mean per photo in transect/CVFold",
       xlab="Photo number (from left)",
       ylab="Mean seal count")
  k=0.25
  for(i in 1:length(photoVec))
  {
    a <- photoVec[photoOrder][i]
    b <- posteriorFieldAtPhoto$quant095[photoOrder][i]
    c <- posteriorFieldAtPhoto$quant005[photoOrder][i]
    polygon(c(a-k, a+k, a+k, a-k), c(c, c, b, b), col="gray90", border="gray90")
    b <- posteriorFieldAtPhoto$quant075[photoOrder][i]
    c <- posteriorFieldAtPhoto$quant025[photoOrder][i]
    polygon(c(a-k, a+k, a+k, a-k), c(c, c, b, b), col="gray70", border="gray70")

    d <- posteriorFieldAtPhoto$median[photoOrder][i]
    lines(c(a-k, a+k), c(d,d), lwd=2, col="red")

    e <- posteriorFieldAtPhoto$mean[photoOrder][i]
    lines(c(a-k, a+k), c(e,e), lwd=2, col="blue")

  }

  points(photoVec[photoOrder],posteriorFieldAtPhoto$photoCounts[photoOrder],lwd=1)
  legend("topright",c("True Count","Pred median","Pred mean","90 CI","50% CI"),col=c(1,"red","blue","grey90","grey70"),lty=c(1,1,1,1,1),lwd=c(NA,1,1,10,10),pch=c(1,NA,NA,NA,NA),bg="white")
  legend("topleft",c(paste("RMSE mean meanPred = ",round(rmse.meanPredMEAN,2),sep=""),
                     paste("MAE median meanPred = ",round(mae.medianPredMEAN,2),sep=""),
                     paste("Prop truth in 50% CI = ",round(PropTrueCountIn050CIMEAN,2),sep=""),
                     paste("Prop truth in 90% CI = ",round(PropTrueCountIn090CIMEAN,2),sep="")),
         pch=NA,col="white")
  dev.off()









  }




  save.list <- list(FhatMat=FhatPhotoMat,sampWithinMat=sampWithinMat,sampVec=sampVec,inputVar=inputVar,predPhotoDf=predPhotoDf,
                    posteriorevalFullPhoto = counts,posteriorDistPhoto = posteriorDistPhotoMat,
                    trans.list = trans.list) # inputVar included here for easy check of
  saveRDS(save.list,file = file.path(savingFolder,paste("photo_pred_comparisons_transect_",paste(leaveOutTransect,collapse = "&"),".rds",sep="")))
}


#' Comparison of photos counts and predictions for GAM
#'
#' Computes basic summary statistics based on distribution
#'
#' @param posteriorevalFullPhoto Numeric vector of evaluation points used for all posteriorDistributions
#' @param posteriorDistPhoto Matrix where each column is the posterior distribution probability for one photo
#' @param posteriorevalFullTransect Numeric vector of evaluation points used for posteriorDistTransect
#' @param posteriorDistTransect Numeric vector of posterior distribution probabilitites for the whole transect in question
#' @param photoCounts Numeric vector with the observed photo counts (not used in the modelling)
#' @param photoOrder Numeric vector specifying the plotting order of the photos (smallest to largest original x-coordinate)
#' @param leaveOutTransect The transect in question
#' @param inputVar List of all input variables used in the original run
#' @param savingFolder String with the path for where to save these comparison results
#' @return Nothing, just save the output to file
#' @import Hmisc
#' @export

ComparePhotoCountsAndPredGAM <- function(posteriorevalFullPhoto,
                                         posteriorDistPhoto,
                                         posteriorevalFullTransect,
                                         posteriorDistTransect,
                                         areaPerSubSampPhotoMatrix,
                                         photoCounts,
                                         photoOrder,
                                         leaveOutTransect,
                                         inputVar,
                                         savingFolder,
                                         plot= TRUE){

  ## For each photo
  counts <- 0:(max(c(photoCounts,unlist(posteriorevalFullPhoto)))+1)
  posteriorDistPhotoMat <- FhatPhotoMat <- matrix(0,ncol=length(counts),nrow=length(photoCounts))
  for (i in 1:length(photoCounts)){
    theseEvals <- which(counts %in% posteriorevalFullPhoto[[i]])
    posteriorDistPhotoMat[i,theseEvals] <- posteriorDistPhoto[[i]]
    FhatPhotoMat[i,] <- cumsum(posteriorDistPhotoMat[i,]) ### Adding an extra zero here
  }


  sampWithinMat <- matrix(NA,ncol=2,nrow=length(photoCounts))
  sampVec <- rep(NA,length(photoCounts))
  for (i in 1:length(photoCounts)){
    evalPoint <- which(counts==photoCounts[i])
    FhatPhotoVec <- c(0,FhatPhotoMat[i,])
    sampWithinMat[i,] <- c(FhatPhotoVec[evalPoint],FhatPhotoVec[evalPoint+1])
    sampVec[i] <- runif(n=1,min = sampWithinMat[i,1],max= sampWithinMat[i,2])
  }

  meanVec <- as.vector(t(counts) %*% t(posteriorDistPhotoMat))

  quant.finder <- function(vec,q,evalvec){
    vec.q = vec-q
    evalvec[vec.q>0][1]
  }

  quantMat <- matrix(NA,ncol=5,nrow=length(photoCounts))
  colnames(quantMat) = c("quant005","quant025","median","quant075","quant095")
  for (i in 1:length(photoCounts)){
    quantMat[i,1] <- quant.finder(vec=FhatPhotoMat[i,],q=0.05,evalvec = counts)
    quantMat[i,2] <- quant.finder(vec=FhatPhotoMat[i,],q=0.25,evalvec = counts)
    quantMat[i,3] <- quant.finder(vec=FhatPhotoMat[i,],q=0.5,evalvec = counts)
    quantMat[i,4] <- quant.finder(vec=FhatPhotoMat[i,],q=0.75,evalvec = counts)
    quantMat[i,5] <- quant.finder(vec=FhatPhotoMat[i,],q=0.95,evalvec = counts)
  }

  predPhotoDf <- as.data.frame(quantMat)
  predPhotoDf$mean <- meanVec
  predPhotoDf$photoCounts <- photoCounts
  predPhotoDf$photoOrder <- photoOrder

  countsTrans <- 0:(max(posteriorevalFullTransect)+1)
  posteriorDistTransectNew <- rep(0,length(countsTrans))
  theseEvals <- which(countsTrans %in% posteriorevalFullTransect)
  posteriorDistTransectNew[theseEvals] <- posteriorDistTransect
  FhatTrans <- cumsum(posteriorDistTransectNew)

  meanTrans <- c(t(countsTrans)%*%posteriorDistTransectNew)
  medTrans <- quant.finder(vec=FhatTrans,q=0.5,evalvec = countsTrans)

  quant005Trans <- quant.finder(vec=FhatTrans,q=0.05,evalvec = counts)
  quant025Trans <- quant.finder(vec=FhatTrans,q=0.25,evalvec = counts)
  quant075Trans <- quant.finder(vec=FhatTrans,q=0.75,evalvec = counts)
  quant095Trans <- quant.finder(vec=FhatTrans,q=0.95,evalvec = counts)


  trans.list <- list(posteriorevalFullTransect=countsTrans,
                     posteriorDistTransect = posteriorDistTransectNew,
                     FhatTrans = FhatTrans,
                     meanTrans = meanTrans,
                     medTrans = medTrans,
                     quant095Trans = quant095Trans,
                     quant005Trans = quant005Trans,
                     quant075Trans = quant075Trans,
                     quant025Trans = quant025Trans)

  posteriorFieldAtPhoto = data.frame(median = apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=median))

  posteriorFieldAtPhoto$mean <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=mean)
  posteriorFieldAtPhoto$quant005 <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=quantile, probs=0.05)
  posteriorFieldAtPhoto$quant025 <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=quantile, probs=0.25)
  posteriorFieldAtPhoto$quant075 <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=quantile, probs=0.75)
  posteriorFieldAtPhoto$quant095 <-  apply(areaPerSubSampPhotoMatrix,MARGIN=2,FUN=quantile, probs=0.95)
  posteriorFieldAtPhoto$photoCounts <- photoCounts

  PropTrueCountIn050CIMEAN <- mean(posteriorFieldAtPhoto$quant025 <= posteriorFieldAtPhoto$photoCounts & posteriorFieldAtPhoto$quant075 >= posteriorFieldAtPhoto$photoCounts)
  PropTrueCountIn090CIMEAN <- mean(posteriorFieldAtPhoto$quant005 <= posteriorFieldAtPhoto$photoCounts & posteriorFieldAtPhoto$quant095 >= posteriorFieldAtPhoto$photoCounts)

  (rmse.meanPredMEAN <- sqrt(mean((posteriorFieldAtPhoto$photoCounts-posteriorFieldAtPhoto$mean)^2)))
  (mae.medianPredMEAN <- mean(abs(posteriorFieldAtPhoto$photoCounts-posteriorFieldAtPhoto$median)))





  if (plot){

  pdf(file=file.path(savingFolder,"photo_counts_and_predictions.pdf"),width=7,height=4)
  par(mar=c(3,3,2,2),mgp=2:0,cex=0.8)

  plot(1:length(predPhotoDf$photoCounts),predPhotoDf$median[photoOrder],type='n',ylim=c(0,max(c(predPhotoDf$quant095,predPhotoDf$photoCounts))),
       main = paste("Photo counts and prediction comparison for transect/CVFold ",paste(leaveOutTransect,collapse = "&"),sep=""),ylab="Predicted/observed counts",xlab="Photo in transect (left to right)")

  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$quant005[photoOrder],lty=1,lwd=1,col=2)
  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$quant095[photoOrder],lty=1,lwd=1,col=2)

  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$quant025[photoOrder],lty=1,lwd=1,col=4)
  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$quant075[photoOrder],lty=1,lwd=1,col=4)


  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$median[photoOrder],col=1,lwd=2)
  lines(1:length(predPhotoDf$photoCounts),predPhotoDf$photoCounts[photoOrder],col=3,lty=2,lwd=1.5)

  legend("topleft",c("Post.pred. median","Post.pred 90% CI","Post.pred 50% CI","Observed counts"),col=c(1,2,4,3),lty=c(1,1,1,2),lwd=c(2,1,1,1.5))

  dev.off()


  ######

  PropTrueCountIn050CI <- mean(predPhotoDf$quant025 <= predPhotoDf$photoCounts & predPhotoDf$quant075 >= predPhotoDf$photoCounts)
  PropTrueCountIn090CI <- mean(predPhotoDf$quant005 <= predPhotoDf$photoCounts & predPhotoDf$quant095 >= predPhotoDf$photoCounts)

  (rmse.meanPred <- sqrt(mean((predPhotoDf$photoCounts-predPhotoDf$mean)^2)))
  (mae.medianPred <- mean(abs(predPhotoDf$photoCounts-predPhotoDf$median)))
  photoVec <- 1:length(predPhotoDf$photoCounts)


  pdf(file=file.path(savingFolder,"photo_counts_and_predictions2.pdf"),width=10,heigh=5)
  par(mfrow=c(1,1))
  plot(photoVec,predPhotoDf$photoCounts[photoOrder],type='n',ylim=c(0,max(predPhotoDf$quant095,predPhotoDf$photoCounts)),lwd=2,
       main="GAM\n Comparison of true counts and out-of-sample posterior predictions per photo in transect/CVFold",
       xlab="Photo number (from left)",
       ylab="Seal count")
  k=0.25
  for(i in 1:length(photoVec))
  {
    a <- photoVec[photoOrder][i]
    b <- predPhotoDf$quant095[photoOrder][i]
    c <- predPhotoDf$quant005[photoOrder][i]
    polygon(c(a-k, a+k, a+k, a-k), c(c, c, b, b), col="gray90", border="gray90")
    b <- predPhotoDf$quant075[photoOrder][i]
    c <- predPhotoDf$quant025[photoOrder][i]
    polygon(c(a-k, a+k, a+k, a-k), c(c, c, b, b), col="gray70", border="gray70")

    d <- predPhotoDf$median[photoOrder][i]
    lines(c(a-k, a+k), c(d,d), lwd=2, col="red")

    e <- predPhotoDf$mean[photoOrder][i]
    lines(c(a-k, a+k), c(e,e), lwd=2, col="blue")

  }

  points(photoVec[photoOrder],predPhotoDf$photoCounts[photoOrder],lwd=1)
  legend("topright",c("True Count","Pred median","Pred mean","90 CI","50% CI"),col=c(1,"red","blue","grey90","grey70"),lty=c(1,1,1,1,1),lwd=c(NA,1,1,10,10),pch=c(1,NA,NA,NA,NA),bg="white")
  legend("topleft",c(paste("RMSE mean pred = ",round(rmse.meanPred,2),sep=""),
                     paste("MAE median pred = ",round(mae.medianPred,2),sep=""),
                     paste("Prop truth in 50% CI = ",round(PropTrueCountIn050CI,2),sep=""),
                     paste("Prop truth in 90% CI = ",round(PropTrueCountIn090CI,2),sep="")),
         pch=NA,col="white")
  dev.off()



  pdf(file=file.path(savingFolder,"photo_counts_and_posterior_mean.pdf"),width=10,heigh=5)
  par(mfrow=c(1,1))
  plot(photoVec,predPhotoDf$photoCounts[photoOrder],type='n',ylim=c(0,max(posteriorFieldAtPhoto$quant095,predPhotoDf$photoCounts)),lwd=2,
       main="GAM\n Comparison of true counts and out-of-sample posterior mean per photo in transect/CVFold",
       xlab="Photo number (from left)",
       ylab="Mean seal count")
  k=0.25
  for(i in 1:length(photoVec))
  {
    a <- photoVec[photoOrder][i]
    b <- posteriorFieldAtPhoto$quant095[photoOrder][i]
    c <- posteriorFieldAtPhoto$quant005[photoOrder][i]
    polygon(c(a-k, a+k, a+k, a-k), c(c, c, b, b), col="gray90", border="gray90")
    b <- posteriorFieldAtPhoto$quant075[photoOrder][i]
    c <- posteriorFieldAtPhoto$quant025[photoOrder][i]
    polygon(c(a-k, a+k, a+k, a-k), c(c, c, b, b), col="gray70", border="gray70")

    d <- posteriorFieldAtPhoto$median[photoOrder][i]
    lines(c(a-k, a+k), c(d,d), lwd=2, col="red")

    e <- posteriorFieldAtPhoto$mean[photoOrder][i]
    lines(c(a-k, a+k), c(e,e), lwd=2, col="blue")

  }

  points(photoVec[photoOrder],posteriorFieldAtPhoto$photoCounts[photoOrder],lwd=1)
  legend("topright",c("True Count","Pred median","Pred mean","90 CI","50% CI"),col=c(1,"red","blue","grey90","grey70"),lty=c(1,1,1,1,1),lwd=c(NA,1,1,10,10),pch=c(1,NA,NA,NA,NA),bg="white")
  legend("topleft",c(paste("RMSE mean meanPred = ",round(rmse.meanPredMEAN,2),sep=""),
                     paste("MAE median meanPred = ",round(mae.medianPredMEAN,2),sep=""),
                     paste("Prop truth in 50% CI = ",round(PropTrueCountIn050CIMEAN,2),sep=""),
                     paste("Prop truth in 90% CI = ",round(PropTrueCountIn090CIMEAN,2),sep="")),
         pch=NA,col="white")
  dev.off()




  }


  save.list <- list(FhatMat=FhatPhotoMat,sampWithinMat=sampWithinMat,sampVec=sampVec,inputVar=inputVar,predPhotoDf,
                    posteriorevalFullPhoto = counts,posteriorDistPhoto = posteriorDistPhotoMat,
                    trans.list = trans.list) # inputVar included here for easy check of
  saveRDS(save.list,file = file.path(savingFolder,paste("photo_pred_comparisons_transect_",paste(leaveOutTransect,collapse = "&"),".rds",sep="")))
}




#' Posterior summary statistics
#'
#' Computes basic summary statistics based on distribution
#'
#' @param evalPoints Numeric vector of evaluation points for a distribution
#' @param dist Numeric vector with the probability at the corresponding evaluation points
#' @param results.CI.level Numeric, denoting the confidence/credibility degree to use in the final credibility interval for the total number of counts
#' @param posterior Logical, indicating whether this is a posterior distribution for which the name posterior should be used on the output
#' @return List with basic self-explanatory summary statistics
#' @import Hmisc
#' @export


SummaryStat <- function(evalPoints,
                        dist,
                        results.CI.level = 0.95,
                        posterior=TRUE){

  ## Computes the mean
  meanmean <- sum(evalPoints*dist)

  ## Computes the median, IQR and 95% CI
  CI.quantiles <- c((1-results.CI.level)/2, 1- (1-results.CI.level)/2)

  quantiles <- Hmisc::wtd.quantile(x=evalPoints,weights=dist,normwt=TRUE,na.rm=TRUE,probs=c(CI.quantiles[1],0.25,0.5,0.75,CI.quantiles[2]))
  med <- quantiles[3]
  IQR <- quantiles[4]- quantiles[2]
  CI <- c(quantiles[1], quantiles[5])

  ## Computes  the mode
  mode <- evalPoints[which.max(dist)]

  retret <- list()
  if (posterior){
    retret$posteriorMean <- meanmean
    retret$posteriorMedian <- med
    retret$posteriorMode <- mode
    retret$posteriorIQR <- IQR
    retret$posteriorCI <- CI
  } else {
    retret$mean <- meanmean
    retret$median <- med
    retret$mode <- mode
    retret$IQR <- IQR
    retret$CI <- CI
    }

  return(retret)
}


#' Computing the posterior predictive distribution of total area counts in the location of the PredPhotots
#'
#' Currently only to be used with a spatial model with a linear covariate and no additional iid term
#'
#' @param samp List with all posterior samples, one sample in each sublist (as obtained from inla.posterior.sample)
#' @param thisMeshPoint Numeric vector giving the mesh point corresponding to each of the pred photos
#' @param covAtPredLoc Numeric vector containing the covariate value at the location of the prediction photos
#' @param poisson.maxEvals Numeric giving the maximum number of poisson evaluations
#' @return List containing the resulting posterior predicitve distribution for all PredPoints
#' @keywords inla
#' @export



ComputePostPredDistforPredPhotos <- function(samp,
                                             thisMeshPoint,
                                             covAtPredLoc,
                                             areaPredPhotos,
                                             poisson.maxEvals,
                                             parallelize.numCores,
                                             parallelize.noSplits,
                                             covariate.fitting,
                                             mesh){
  noSamp <- length(samp)

  ## Finds the columns of the latent sample corresponding to the different sampled variables
  extractTypes <- c("rf","iid","intercept","covariate","nonlinear","spatialX","spatialY","spatialXY")

  ids <- lapply(extractTypes,function(x) grep(x,rownames(samp[[1]]$latent),fixed=TRUE))
  names(ids)=extractTypes

  if (covariate.fitting=="linear"){
    etaAtGrid <- function(s,thisMeshPoint,ids,covAtPredLoc,val.spatialX,val.spatialY){
      eta1AtMesh <- s$latent[ids$rf,1]
      eta1AtPredLoc <- eta1AtMesh[thisMeshPoint]
      eta2AtPredLoc <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covAtPredLoc
      etaCombinedAtPredLoc <- eta1AtPredLoc + eta2AtPredLoc
      return(etaCombinedAtPredLoc)
    }
  }
  if (covariate.fitting=="linearAndSpatial"){
    etaAtGrid <- function(s,thisMeshPoint,ids,covAtPredLoc,val.spatialX,val.spatialY){
      eta1AtMesh <- s$latent[ids$rf,1]
      eta1AtPredLoc <- eta1AtMesh[thisMeshPoint]
      eta2AtPredLoc <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covAtPredLoc
      eta3AtPredLoc <- s$latent[ids$spatialX[1],1]*val.spatialX + s$latent[ids$spatialY,1]*val.spatialY + s$latent[ids$spatialXY,1]*sqrt(val.spatialX^2+val.spatialY^2)
      etaCombinedAtPredLoc <- eta1AtPredLoc + eta2AtPredLoc + eta3AtPredLoc
      return(etaCombinedAtPredLoc)
    }
  }

  val.spatialX <- mesh$loc[thisMeshPoint,1]
  val.spatialY <- mesh$loc[thisMeshPoint,2]
  val.spatialXY <- sqrt(val.spatialX^2+val.spatialY^2)




    etaCombinedAtPredLocList <- sapply(samp,etaAtGrid,simplify=TRUE,thisMeshPoint = thisMeshPoint,ids = ids, covAtPredLoc = covAtPredLoc, val.spatialX = val.spatialX, val.spatialY = val.spatialY)

    areaPerPredPhoto <- exp(etaCombinedAtPredLocList)*areaPredPhotos
    areaAllPredPhoto <- colSums(areaPerPredPhoto)

    evalFullPredPhoto <- unique(round(seq(qpois(0.001,min(areaPerPredPhoto,na.rm=TRUE)),min(qpois(0.999,max(areaPerPredPhoto,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
    evalFullAllPredPhoto <- unique(round(seq(qpois(0.001,min(areaAllPredPhoto,na.rm=TRUE)),min(qpois(0.999,max(areaAllPredPhoto,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))

    splittedevalFullPredPhoto=split(evalFullPredPhoto,ceiling((1:length(evalFullPredPhoto))/length(evalFullPredPhoto)*parallelize.noSplits))
    splittedevalFullAllPredPhoto=split(evalFullAllPredPhoto,ceiling((1:length(evalFullAllPredPhoto))/length(evalFullAllPredPhoto)*parallelize.noSplits))




    # Then poisson evaluation for each photo
    doMC::registerDoMC(parallelize.numCores)
    export.var=c("areaPerPredPhoto","splittedevalFullPredPhoto")
    non.export.var=ls()[!(ls()%in%export.var)]
    uu=proc.time()
    posteriorDistPredPhoto <- foreach::foreach(i=1:length(splittedevalFullPredPhoto),.noexport=non.export.var,.verbose=FALSE,.inorder=TRUE) %dopar% {
      lapply(X=splittedevalFullPredPhoto[[i]],FUN=MeanPoissonDistMat,lambdaMat=t(areaPerPredPhoto))
    }

    print(paste("Finished computing the mean of the Poisson counts in each prediction photo in ",round((proc.time()-uu)[3])," seconds.",sep=""))

    posteriorDistPredPhoto <- matrix(unlist(posteriorDistPredPhoto),ncol=length(thisMeshPoint),byrow=T)
    posteriorDistPredPhoto <- posteriorDistPredPhoto%*%diag(1/colSums(posteriorDistPredPhoto))


    ## First poisson evaluation for the transect
    doMC::registerDoMC(parallelize.numCores)
    export.var=c("areaAllPredPhoto","splittedevalFullAllPredPhoto")
    non.export.var=ls()[!(ls()%in%export.var)]
    uu=proc.time()
    posteriorDistAllPredPhoto <- foreach::foreach(i=1:length(splittedevalFullAllPredPhoto),.noexport=non.export.var,.verbose=FALSE,.inorder=TRUE) %dopar% {
      sapply(X=splittedevalFullAllPredPhoto[[i]],FUN=MeanPoissonDist,lambdavec=areaAllPredPhoto)
    }
    print(paste("Finished computing the mean of the Poisson counts for the full transect in ",round((proc.time()-uu)[3])," seconds.",sep=""))
    posteriorDistAllPredPhoto <- unlist(posteriorDistAllPredPhoto)
    posteriorDistAllPredPhoto <- posteriorDistAllPredPhoto/sum(posteriorDistAllPredPhoto)


    retret <- list()
    #  retret$posteriorEvalPoints <- evalFull
    #  retret$posteriorDist <- posteriorDist
    retret$evalFullPredPhoto <- evalFullPredPhoto
    retret$posteriorDistPredPhoto <- posteriorDistPredPhoto
    retret$evalFullAllPredPhoto <- evalFullAllPredPhoto
    retret$posteriorDistAllPredPhoto <- posteriorDistAllPredPhoto
    retret$areaPerPredPhoto <- areaPerPredPhoto

    #  retret$mean.field.samp <- expField
    #  retret$sd.field.samp <- sdField
    #  retret$mean.field.domain.samp <- expFieldDomain
    #  retret$sd.field.domain.samp <- sdFieldDomain

    return(retret)
    print("Finished running the ComputePostPredDistforPredPhotos function")

}


#' Computing the posterior predictive distribution of total area counts
#'
#' Computes the posterior predictive distribution of total area counts based on a set of samples from the posterior field.
#' The function first splits the samples into several sublists and saves them to disk. Then computes the eta for each point in each grid along with the
#' sd and so on. Finally computes the dist of the mean of the Poisson counts corresponding to the posterior predictive dist of total area counts.
#'
#' @param samp List with all posterior samples, one sample in each sublist (as obtained from inla.posterior.sample)
#' @param spatial Logical, indicating whether the model fitted includes a spatial spde term
#' @param parallelize.noSplits Numeric, deciding how many sublists samp should be splitted into. Should be a multiple of parallelize.numCores for the highest efficiency. The larger number the less memory is used (and longer time)
#' @param parallelize.numCores Numeric, corresponding to the number of cores any parallelization should be run at
#' @param tempFolder Path to the folder where temporary files (samp subfiles and gridded eta subfiles)
#' @param use.covariates Logical, indicating whether covariates are used or not (see description!)
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param covariate.fitting String, indicating how to model covariates. "linear", quadratic (default) or "linAndLog", or FALSE for no covariates
#' @param gridList List containing information about the grid, being the output of the function GridCreation
#' @param covGridList List containing information about the covariates at the grid (also when use.covariates = FALSE), being the output of the function covAtNewGrid
#' @param nxy Numberic vector of size 2, giving the dimension in x- and y-direction for the grid
#' @param poisson.maxEvals Numeric, corresponding to maximum number of points the Poisson distribution should be evaluated at (a much smaller number is typically used)
#' @param noMeshPoints Numeric, corresponding to the number of points in the mesh being used
#' @param extraNonlinear.covGridList List of additional projection objects related to the nonlinear covariate effect when applicable
#' @return List containing the resulting posterior predicitve distribution for all evaluation points (also returned) in addition to the mean field and sd field computed from the samples)
#' @keywords inla
#' @export


ComputePostPredDist <- function(samp,
                                spatial,
                                parallelize.noSplits = parallelize.numCores,
                                parallelize.numCores,
                                tempFolder,
                                use.covariates,
                                additional.iid.term,
                                covariate.fitting,
                                gridList,
                                covGridList,
                                nxy,
                                poisson.maxEvals,
                                noMeshPoints,
                                extraNonlinear.covGridList,
                                delete.samp = T) {

  noSamp <- length(samp)

  ## Finds the columns of the latent sample corresponding to the different sampled variables
  extractTypes <- c("rf","iid","intercept","covariate","nonlinear","spatialX","spatialY","spatialXY")

  ids <- lapply(extractTypes,function(x) grep(x,rownames(samp[[1]]$latent),fixed=TRUE))
  names(ids)=extractTypes

  ## Splits the samp object into smaller lists, saves them to disk and them delete them from RAM

  splittedSamp=split(samp,ceiling((1:noSamp)/noSamp*parallelize.noSplits))

  for (i in 1:parallelize.noSplits){
    saveRDS(splittedSamp[[i]],file = paste(tempFolder,"/samp_",i,".rds",sep=""))
    print(paste("Wrote splitted samp file ",i," of ",parallelize.noSplits," to disk",sep=""))
  }
  if (delete.samp){
    rm("samp",envir =sys.frame(-1))
    rm("samp","splittedSamp")
    gc()
  }
  print("Finished saving the splitted samp files to the temp-folder")

  ## Creates a function for extracting the eta (logIntensity of each sample)

  if (spatial){
    if (!use.covariates){
      if (!additional.iid.term){
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
          eta1AtMesh <- s$latent[ids$rf,1]
          eta1AtGrid <- inla.mesh.project(project=projgrid,field=eta1AtMesh)

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      } else {
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
      }
    }
    if (use.covariates){
      if (!additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
            }
        }
        if (covariate.fitting=="linearAndSpatial"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)
            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval
            eta3AtGrid <- s$latent[ids$spatialX[1],1]*val.spatialX + s$latent[ids$spatialY,1]*val.spatialY + s$latent[ids$spatialXY,1]*sqrt(val.spatialX^2+val.spatialY^2)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      if (additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
    }
  }
    } else {
    if (!use.covariates){
      if (!additional.iid.term){
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
          eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
          }
        } else {
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
    }
    if (use.covariates){
      if (!additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linearAndSpatial"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))
            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval
            eta3AtGrid <- s$latent[ids$spatialX[1],1]*val.spatialX + s$latent[ids$spatialY,1]*val.spatialY + s$latent[ids$spatialXY,1]*sqrt(val.spatialX^2+val.spatialY^2)
            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
      if (additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
    }
  }



  doMC::registerDoMC(parallelize.numCores)

  export.var=c("etaAtGrid","ids","gridList","covGridList","tempFolder","noMeshPoints","extraNonlinear.covGridList") # Not functions here
  non.export.var=ls()[!(ls()%in%export.var)]

  uu=proc.time()

  parallelSampHandling <- foreach::foreach(i=1:parallelize.noSplits,.noexport=non.export.var,.packages="INLA",.verbose=FALSE,.inorder=FALSE) %dopar% {

    # Reading a part of the subsampled list
    thisSubSamp <- readRDS(file = paste(tempFolder,"/samp_",i,".rds",sep=""))


    val.spatialX <- matrix(rep(gridList$projgrid$x,nxy[2]),ncol=nxy[2],nrow=nxy[1])
    val.spatialY <- matrix(rep(gridList$projgrid$y,nxy[1]),ncol=nxy[2],nrow=nxy[1],byrow = T)
    val.spatialXY <- sqrt(val.spatialX^2+val.spatialY^2)


    # Computing the eta for all grid points in each of the list in the subsampled list and write it to file
    etaAtGridListSub=sapply(thisSubSamp,etaAtGrid,simplify = TRUE,projgrid = gridList$projgrid,covNewGridval = covGridList$covariateValues,ids=ids,noMeshPoints=noMeshPoints,extraNonlinear.covGridList=extraNonlinear.covGridList,
                            val.spatialX = val.spatialX,val.spatialY = val.spatialY)
    #saveRDS(etaAtGridListSub,file = paste(tempFolder,"/etaAtGridList_",i,".rds",sep="")) # Really no point in saving these really big files.

    # Computing the empirical mean field value at each gridpoint over the subsampled list -- strictly not need, but computed to check results.
    expContrib <- rowMeans(etaAtGridListSub)

    # Computing the squared empirical mean field value at each gridpoint over the subsampled list -- to be used for computing the mean/variance
    squaredExpContrib <- rowMeans(etaAtGridListSub^2)

    # Computing the integrated intensity, i.e. the expected number of Poisson distributed counts for each of the sampled fields in the subsamp
    areaPerSubSamp <- IntegrateVectorizedGrids(arrayGrid=etaAtGridListSub,
                                               logicalGridPointsInsideCountingDomain=gridList$logicalGridPointsInsideCountingDomain,
                                               truePixelSize=gridList$truePixelSize,
                                               scaleArea=gridList$scaleArea)
    print(paste("Finished parallelSamphandling for core",i,sep="")) # Tries to print to the terminal how far it has gotton

    ret0List <- list()
    ret0List$expContrib <- expContrib
    ret0List$squaredExpContrib <- squaredExpContrib
    ret0List$areaPerSubSamp <- areaPerSubSamp
    ret0List
  }
  print(paste("Finished all handling of the original samp files in ",round((proc.time()-uu)[3])," seconds.",sep=""))

  ## Re-arranges the output from the parallelization
  expContribList <- list()
  squaredExpContribList <- list()
  areaPerSubSampVec <- NULL
  for (i in 1:parallelize.noSplits){
    expContribList[[i]] <- parallelSampHandling[[i]]$expContrib
    squaredExpContribList[[i]] <- parallelSampHandling[[i]]$squaredExpContrib
    areaPerSubSampVec <- c(areaPerSubSampVec,parallelSampHandling[[i]]$areaPerSubSamp)
  }

  ## Computes E[X] and E[X^2] for X the posterior in each location of the grid, transform to a matrix and computes the corresponding pointwise sd of the field
  squaredExpField <- matrix(rowMeans(simplify2array(squaredExpContribList)),ncol=nxy[2],nrow=nxy[1])
  expField <- matrix(rowMeans(simplify2array(expContribList)),ncol=nxy[2],nrow=nxy[1]) # This is actually know, but "better" to compute it as it does not give NAs in sd computation below due to approx error.

  sdField <- matrix(sqrt(squaredExpField - expField^2),ncol=nxy[2],nrow=nxy[1])

  ## Computing the Poisson distribution having the sampled total intensities as the mean ##

  # Finding the evaluations points and splits them into several sub-vectors
  evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
# old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))

  splittedevalFull=split(evalFull,ceiling((1:length(evalFull))/length(evalFull)*parallelize.noSplits))

  doMC::registerDoMC(parallelize.numCores)

  export.var=c("areaPerSubSampVec","splittedevalFull")
  non.export.var=ls()[!(ls()%in%export.var)]

  uu=proc.time()
  posteriorDist <- foreach::foreach(i=1:parallelize.noSplits,.noexport=non.export.var,.verbose=FALSE,.inorder=TRUE) %dopar% {
    sapply(X=splittedevalFull[[i]],FUN=MeanPoissonDist,lambdavec=areaPerSubSampVec)
  }
  print(paste("Finished computing the mean of the Poisson counts in ",round((proc.time()-uu)[3])," seconds.",sep=""))

  posteriorDist <- unlist(posteriorDist)
  posteriorDist <- posteriorDist/sum(posteriorDist)

  expFieldDomain <- expField
  expFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA
  sdFieldDomain <- sdField
  sdFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA


  retret <- list()
  retret$posteriorEvalPoints <- evalFull
  retret$posteriorDist <- posteriorDist
  retret$mean.field.samp <- expField
  retret$sd.field.samp <- sdField
  retret$mean.field.domain.samp <- expFieldDomain
  retret$sd.field.domain.samp <- sdFieldDomain


  return(retret)
  print("Finished running the ComputePostPredDist function")

}


#' Computing the posterior predictive distribution for single photos and the sums of these photos in addition to the regular counting domain
#'
#' As for the ComputePostPredDist, but with the additional stuff for doing it also per photo
#'
#' @param samp List with all posterior samples, one sample in each sublist (as obtained from inla.posterior.sample)
#' @param spatial Logical, indicating whether the model fitted includes a spatial spde term
#' @param parallelize.noSplits Numeric, deciding how many sublists samp should be splitted into. Should be a multiple of parallelize.numCores for the highest efficiency. The larger number the less memory is used (and longer time)
#' @param parallelize.numCores Numeric, corresponding to the number of cores any parallelization should be run at
#' @param tempFolder Path to the folder where temporary files (samp subfiles and gridded eta subfiles)
#' @param use.covariates Logical, indicating whether covariates are used or not (see description!)
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param covariate.fitting String, indicating how to model covariates. "linear", quadratic (default) or "linAndLog", or FALSE for no covariates
#' @param gridList List containing information about the grid, being the output of the function GridCreation
#' @param covGridList List containing information about the covariates at the grid (also when use.covariates = FALSE), being the output of the function covAtNewGrid
#' @param nxy Numberic vector of size 2, giving the dimension in x- and y-direction for the grid
#' @param poisson.maxEvals Numeric, corresponding to maximum number of points the Poisson distribution should be evaluated at (a much smaller number is typically used)
#' @param noMeshPoints Numeric, corresponding to the number of points in the mesh being used
#' @param extraNonlinear.covGridList List of additional projection objects related to the nonlinear covariate effect when applicable
#' @param predPhotoGridPointList List of which grid point which define the photos when sampling
#' @return List containing the resulting posterior predicitve distribution for all evaluation points (also returned) in addition to the mean field and sd field computed from the samples)
#' @keywords inla
#' @export


PhotoPostPredDist <- function(samp,
                                spatial,
                                parallelize.noSplits = parallelize.numCores,
                                parallelize.numCores,
                                tempFolder,
                                use.covariates,
                                additional.iid.term,
                                covariate.fitting,
                                gridList,
                                covGridList,
                                nxy,
                                poisson.maxEvals,
                                noMeshPoints,
                                extraNonlinear.covGridList,
                                predPhotoGridPointList) {


  noSamp <- length(samp)

  ## Finds the columns of the latent sample corresponding to the different sampled variables
  extractTypes <- c("rf","iid","intercept","covariate","nonlinear","spatialX","spatialY","spatialXY")

  ids <- lapply(extractTypes,function(x) grep(x,rownames(samp[[1]]$latent),fixed=TRUE))
  names(ids)=extractTypes

  ## Splits the samp object into smaller lists, saves them to disk and them delete them from RAM

  splittedSamp=split(samp,ceiling((1:noSamp)/noSamp*parallelize.noSplits))

  for (i in 1:parallelize.noSplits){
    saveRDS(splittedSamp[[i]],file = paste(tempFolder,"/samp_",i,".rds",sep=""))
    print(paste("Wrote splitted samp file ",i," of ",parallelize.noSplits," to disk",sep=""))
  }
  rm("samp",envir =sys.frame(-1))
  rm("samp","splittedSamp")
  print("Finished saving the splitted samp files to the temp-folder")

  ## Creates a function for extracting the eta (logIntensity of each sample)

  if (spatial){
    if (!use.covariates){
      if (!additional.iid.term){
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
          eta1AtMesh <- s$latent[ids$rf,1]
          eta1AtGrid <- inla.mesh.project(project=projgrid,field=eta1AtMesh)

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      } else {
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
          eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
          eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      }
    }
    if (use.covariates){
      if (!additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linearAndSpatial"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)
            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval
            eta3AtGrid <- s$latent[ids$spatialX[1],1]*val.spatialX + s$latent[ids$spatialY,1]*val.spatialY + s$latent[ids$spatialXY,1]*sqrt(val.spatialX^2+val.spatialY^2)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (additional.iid.term){
          if (covariate.fitting=="linear"){
            etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
              eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
              eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

              eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

              etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
              return(etaCombinedAtGrid)
            }
          }
          if (covariate.fitting=="quadratic"){
            etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
              eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
              eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

              eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

              etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
              return(etaCombinedAtGrid)
            }
          }
          if (covariate.fitting=="linAndLog"){
            etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
              eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
              eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

              eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

              etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
              return(etaCombinedAtGrid)
            }
          }
          if (covariate.fitting=="nonlinear"){
            etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
              eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
              eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

              eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
              eta2AtMesh <- s$latent[ids$nonlinear,1]
              eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
              eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

              eta3AtGrid <- s$latent[ids$intercept,1]

              etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
              return(etaCombinedAtGrid)
            }
          }
        }
      }
    }
  } else {
    if (!use.covariates){
      if (!additional.iid.term){
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
          eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      } else {
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
          eta1AtMesh <- 0 + s$latent[ids$iid,1]
          eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      }
    }
    if (use.covariates){
      if (!additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linearAndSpatial"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))
            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval
            eta3AtGrid <- s$latent[ids$spatialX[1],1]*val.spatialX + s$latent[ids$spatialY,1]*val.spatialY + s$latent[ids$spatialXY,1]*sqrt(val.spatialX^2+val.spatialY^2)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
      if (additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList,val.spatialX,val.spatialY){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
    }
  }


  doMC::registerDoMC(parallelize.numCores)

  export.var=c("etaAtGrid","ids","gridList","covGridList","tempFolder","noMeshPoints","extraNonlinear.covGridList","predPhotoGridPointList") # Not functions here
  non.export.var=ls()[!(ls()%in%export.var)]

  uu=proc.time()

  parallelSampHandling <- foreach::foreach(i=1:parallelize.noSplits,.noexport=non.export.var,.packages="INLA",.verbose=FALSE,.inorder=FALSE) %dopar% {

    # Reading a part of the subsampled list
    thisSubSamp <- readRDS(file = paste(tempFolder,"/samp_",i,".rds",sep=""))



    val.spatialX <- matrix(rep(gridList$projgrid$x,nxy[2]),ncol=nxy[2],nrow=nxy[1])
    val.spatialY <- matrix(rep(gridList$projgrid$y,nxy[1]),ncol=nxy[2],nrow=nxy[1],byrow = T)
    val.spatialXY <- sqrt(val.spatialX^2+val.spatialY^2)


    # Computing the eta for all grid points in each of the list in the subsampled list and write it to file
    etaAtGridListSub=sapply(thisSubSamp,etaAtGrid,simplify = TRUE,projgrid = gridList$projgrid,
                            covNewGridval = covGridList$covariateValues,ids=ids,noMeshPoints=noMeshPoints,
                            extraNonlinear.covGridList=extraNonlinear.covGridList,
                            val.spatialX = val.spatialX,val.spatialY = val.spatialY)
      #saveRDS(etaAtGridListSub,file = paste(tempFolder,"/etaAtGridList_",i,".rds",sep="")) # Really no point in saving these really big files.

    # Computing the empirical mean field value at each gridpoint over the subsampled list -- strictly not need, but computed to check results.
#    expContrib <- rowMeans(etaAtGridListSub)

    # Computing the squared empirical mean field value at each gridpoint over the subsampled list -- to be used for computing the mean/variance
#    squaredExpContrib <- rowMeans(etaAtGridListSub^2)


    # Computing the integrated intensity, i.e. the expected number of Poisson distributed counts for each of the sampled fields in the subsamp
#    areaPerSubSamp <- IntegrateVectorizedGrids(arrayGrid=etaAtGridListSub,
#                                               logicalGridPointsInsideCountingDomain=gridList$logicalGridPointsInsideCountingDomain,
#                                               truePixelSize=gridList$truePixelSize,
#                                               scaleArea=gridList$scaleArea)

    areaPerSubSampListPhoto=lapply(predPhotoGridPointList,FUN=lapplyIntegrateVectorizedGrids,arrayGrid=etaAtGridListSub,truePixelSize=gridList$truePixelSize)

    #areaPerSubSampListPhoto=lapply(predPhotoGridPointList,FUN=lapplyIntegrateVectorizedGrids,arrayGrid=etaAtGridListSub,truePixelSize=gridList$truePixelSize,allPhotosinGrid=allPhotosinGrid)
    print(paste("Finished parallelSamphandling for core",i,sep="")) # Tries to print to the terminal how far it has gotton

    ret0List <- list()
#    ret0List$expContrib <- expContrib
#    ret0List$squaredExpContrib <- squaredExpContrib
#    ret0List$areaPerSubSamp <- areaPerSubSamp
    ret0List$areaPerSubSampListPhoto <- areaPerSubSampListPhoto
    ret0List
  }
  print(paste("Finished all handling of the original samp files in ",round((proc.time()-uu)[3])," seconds.",sep=""))

  ## Re-arranges the output from the parallelization
#  expContribList <- list()
#  squaredExpContribList <- list()
#  areaPerSubSampVec <- NULL
  areaPerSubSampPhotoMatrix <- NULL
  for (i in 1:parallelize.noSplits){
#    expContribList[[i]] <- parallelSampHandling[[i]]$expContrib
#    squaredExpContribList[[i]] <- parallelSampHandling[[i]]$squaredExpContrib
#    areaPerSubSampVec <- c(areaPerSubSampVec,parallelSampHandling[[i]]$areaPerSubSamp)
    areaPerSubSampPhotoMatrix <- rbind(areaPerSubSampPhotoMatrix,matrix(unlist(parallelSampHandling[[i]]$areaPerSubSampListPhoto),ncol=length(predPhotoGridPointList)))
  }

  ## Creating also the subsampled area intensity for the union of the photos
  areaPerSubSampTransectVec <- rowSums(areaPerSubSampPhotoMatrix)

  ## Computes E[X] and E[X^2] for X the posterior in each location of the grid, transform to a matrix and computes the corresponding pointwise sd of the field
#  squaredExpField <- matrix(rowMeans(simplify2array(squaredExpContribList)),ncol=nxy[2],nrow=nxy[1])
#  expField <- matrix(rowMeans(simplify2array(expContribList)),ncol=nxy[2],nrow=nxy[1]) # This is actually know, but "better" to compute it as it does not give NAs in sd computation below due to approx error.
#  sdField <- matrix(sqrt(squaredExpField - expField^2),ncol=nxy[2],nrow=nxy[1])

  ## Computing the Poisson distribution having the sampled total intensities as the mean ##

  # Finding the evaluations points and splits them into several sub-vectors # FOR FULL AREA
#  evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
#  splittedevalFull=split(evalFull,ceiling((1:length(evalFull))/length(evalFull)*parallelize.noSplits))

  # Finding the evaluations points and splits them into several sub-vectors # FOR PHOTOS
  evalFullPhoto <- unique(round(seq(qpois(0.001,min(areaPerSubSampPhotoMatrix,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampPhotoMatrix,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
  splittedevalFullPhoto=split(evalFullPhoto,ceiling((1:length(evalFullPhoto))/length(evalFullPhoto)*parallelize.noSplits))

  # Finding the evaluations points and splits them into several sub-vectors # FOR TRANSECT
  evalFullTransect <- unique(round(seq(qpois(0.001,min(areaPerSubSampTransectVec,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampTransectVec,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
  splittedevalFullTransect=split(evalFullTransect,ceiling((1:length(evalFullTransect))/length(evalFullTransect)*parallelize.noSplits))



#  ## First poisson evaluation for the counting domain
#  doMC::registerDoMC(parallelize.numCores)
#  export.var=c("areaPerSubSampVec","splittedevalFull")
#  non.export.var=ls()[!(ls()%in%export.var)]
#  uu=proc.time()
#  posteriorDist <- foreach::foreach(i=1:length(splittedevalFull),.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {
#    sapply(X=splittedevalFull[[i]],FUN=MeanPoissonDist,lambdavec=areaPerSubSampVec)
#  }
#  print(paste("Finished computing the mean of the Poisson counts in countingDomain in ",round((proc.time()-uu)[3])," seconds.",sep=""))
#  posteriorDist <- unlist(posteriorDist)
#  posteriorDist <- posteriorDist/sum(posteriorDist)


  # Then poisson evaluation for each photo
  doMC::registerDoMC(parallelize.numCores)
  export.var=c("areaPerSubSampPhotoMatrix","splittedevalFullPhoto")
  non.export.var=ls()[!(ls()%in%export.var)]
  uu=proc.time()
  posteriorDistPhoto <- foreach::foreach(i=1:length(splittedevalFullPhoto),.noexport=non.export.var,.verbose=FALSE,.inorder=TRUE) %dopar% {
    lapply(X=splittedevalFullPhoto[[i]],FUN=MeanPoissonDistMat,lambdaMat=areaPerSubSampPhotoMatrix)
  }
  print(paste("Finished computing the mean of the Poisson counts in each prediction photo in ",round((proc.time()-uu)[3])," seconds.",sep=""))

  posteriorDistPhoto <- matrix(unlist(posteriorDistPhoto),ncol=length(predPhotoGridPointList),byrow=T)
  posteriorDistPhoto <- posteriorDistPhoto%*%diag(1/colSums(posteriorDistPhoto))


  ## First poisson evaluation for the transect
  doMC::registerDoMC(parallelize.numCores)
  export.var=c("areaPerSubSampTransectVec","splittedevalFullTransect")
  non.export.var=ls()[!(ls()%in%export.var)]
  uu=proc.time()
  posteriorDistTransect <- foreach::foreach(i=1:length(splittedevalFullTransect),.noexport=non.export.var,.verbose=FALSE,.inorder=TRUE) %dopar% {
    sapply(X=splittedevalFullTransect[[i]],FUN=MeanPoissonDist,lambdavec=areaPerSubSampTransectVec)
  }
  print(paste("Finished computing the mean of the Poisson counts for the full transect in ",round((proc.time()-uu)[3])," seconds.",sep=""))
  posteriorDistTransect <- unlist(posteriorDistTransect)
  posteriorDistTransect <- posteriorDistTransect/sum(posteriorDistTransect)



#  expFieldDomain <- expField
#  expFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA
#  sdFieldDomain <- sdField
#  sdFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA


  retret <- list()
#  retret$posteriorEvalPoints <- evalFull
#  retret$posteriorDist <- posteriorDist
  retret$posteriorevalFullPhoto <- evalFullPhoto
  retret$posteriorDistPhoto <- posteriorDistPhoto
  retret$posteriorevalFullTransect <- evalFullTransect
  retret$posteriorDistTransect <- posteriorDistTransect
  retret$areaPerSubSampPhotoMatrix <- areaPerSubSampPhotoMatrix

#  retret$mean.field.samp <- expField
#  retret$sd.field.samp <- sdField
#  retret$mean.field.domain.samp <- expFieldDomain
#  retret$sd.field.domain.samp <- sdFieldDomain

  return(retret)
  print("Finished running the PhotoPostPredDist function")

}

#' Computing the posterior predictive distribution for single photos and the sums of these photos in addition to the regular counting domain
#'
#' As for the ComputePostPredDist, but with the additional stuff for doing it also per photo
#'
#' @param samp List with all posterior samples, one sample in each sublist (as obtained from inla.posterior.sample)
#' @param spatial Logical, indicating whether the model fitted includes a spatial spde term
#' @param parallelize.noSplits Numeric, deciding how many sublists samp should be splitted into. Should be a multiple of parallelize.numCores for the highest efficiency. The larger number the less memory is used (and longer time)
#' @param parallelize.numCores Numeric, corresponding to the number of cores any parallelization should be run at
#' @param tempFolder Path to the folder where temporary files (samp subfiles and gridded eta subfiles)
#' @param use.covariates Logical, indicating whether covariates are used or not (see description!)
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param covariate.fitting String, indicating how to model covariates. "linear", quadratic (default) or "linAndLog", or FALSE for no covariates
#' @param gridList List containing information about the grid, being the output of the function GridCreation
#' @param covGridList List containing information about the covariates at the grid (also when use.covariates = FALSE), being the output of the function covAtNewGrid
#' @param nxy Numberic vector of size 2, giving the dimension in x- and y-direction for the grid
#' @param poisson.maxEvals Numeric, corresponding to maximum number of points the Poisson distribution should be evaluated at (a much smaller number is typically used)
#' @param noMeshPoints Numeric, corresponding to the number of points in the mesh being used
#' @param extraNonlinear.covGridList List of additional projection objects related to the nonlinear covariate effect when applicable
#' @param predPhotoGridPointList List of which grid point which define the photos when sampling
#' @return List containing the resulting posterior predicitve distribution for all evaluation points (also returned) in addition to the mean field and sd field computed from the samples)
#' @keywords inla
#' @export


PhotoPostPredDistold <- function(samp,
                              spatial,
                              parallelize.noSplits = parallelize.numCores,
                              parallelize.numCores,
                              tempFolder,
                              use.covariates,
                              additional.iid.term,
                              covariate.fitting,
                              gridList,
                              covGridList,
                              nxy,
                              poisson.maxEvals,
                              noMeshPoints,
                              extraNonlinear.covGridList,
                              predPhotoGridPointList) {

  noSamp <- length(samp)

  ## Finds the columns of the latent sample corresponding to the different sampled variables
  extractTypes <- c("rf","iid","intercept","covariate","nonlinear")

  ids <- lapply(extractTypes,function(x) grep(x,rownames(samp[[1]]$latent),fixed=TRUE))
  names(ids)=extractTypes

  ## Splits the samp object into smaller lists, saves them to disk and them delete them from RAM

  splittedSamp=split(samp,ceiling((1:noSamp)/noSamp*parallelize.noSplits))

  for (i in 1:parallelize.noSplits){
    saveRDS(splittedSamp[[i]],file = paste(tempFolder,"/samp_",i,".rds",sep=""))
    print(paste("Wrote splitted samp file ",i," of ",parallelize.noSplits," to disk",sep=""))
  }
  rm("samp",envir =sys.frame(-1))
  rm("samp","splittedSamp")
  print("Finished saving the splitted samp files to the temp-folder")

  ## Creates a function for extracting the eta (logIntensity of each sample)

  if (spatial){
    if (!use.covariates){
      if (!additional.iid.term){
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
          eta1AtMesh <- s$latent[ids$rf,1]
          eta1AtGrid <- inla.mesh.project(project=projgrid,field=eta1AtMesh)

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      } else {
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
          eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
          eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      }
    }
    if (use.covariates){
      if (!additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- s$latent[ids$rf,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
      if (additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- s$latent[ids$rf,1] + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
    }
  } else {
    if (!use.covariates){
      if (!additional.iid.term){
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
          eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      } else {
        etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
          eta1AtMesh <- 0 + s$latent[ids$iid,1]
          eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

          eta2AtGrid <- s$latent[ids$intercept,1] # This is a single number, but that is fine here

          etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
          return(etaCombinedAtGrid)
        }
      }
    }
    if (use.covariates){
      if (!additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtGrid <- inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=rep(0,noMeshPoints))

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
      if (additional.iid.term){
        if (covariate.fitting=="linear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="quadratic"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*covNewGridval^2

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="linAndLog"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- s$latent[ids$intercept,1] + s$latent[ids$covariate[1],1]*covNewGridval + s$latent[ids$covariate[2],1]*log(covNewGridval)

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid
            return(etaCombinedAtGrid)
          }
        }
        if (covariate.fitting=="nonlinear"){
          etaAtGrid <- function(s,projgrid,covNewGridval,ids,noMeshPoints,extraNonlinear.covGridList){
            eta1AtMesh <- 0 + s$latent[ids$iid,1]
            eta1AtGrid <- INLA::inla.mesh.project(project=projgrid,field=eta1AtMesh)

            eta2AtGrid <- rep(NA,length(extraNonlinear.covGridList$thesePoints))
            eta2AtMesh <- s$latent[ids$nonlinear,1]
            eta2AtGrid0 <- INLA::inla.mesh.project(project=extraNonlinear.covGridList$projgrid.cov,field=eta2AtMesh)
            eta2AtGrid[extraNonlinear.covGridList$thesePoints] <- eta2AtGrid0

            eta3AtGrid <- s$latent[ids$intercept,1]

            etaCombinedAtGrid <- eta1AtGrid + eta2AtGrid + eta3AtGrid
            return(etaCombinedAtGrid)
          }
        }
      }
    }
  }



  doMC::registerDoMC(parallelize.numCores)

  export.var=c("etaAtGrid","ids","gridList","covGridList","tempFolder","noMeshPoints","extraNonlinear.covGridList","predPhotoGridPointList") # Not functions here
  non.export.var=ls()[!(ls()%in%export.var)]

  uu=proc.time()

  parallelSampHandling <- foreach::foreach(i=1:parallelize.noSplits,.noexport=non.export.var,.packages="INLA",.verbose=FALSE,.inorder=FALSE) %dopar% {

    # Reading a part of the subsampled list
    thisSubSamp <- readRDS(file = paste(tempFolder,"/samp_",i,".rds",sep=""))

    # Computing the eta for all grid points in each of the list in the subsampled list and write it to file
    etaAtGridListSub=sapply(thisSubSamp,etaAtGrid,simplify = TRUE,projgrid = gridList$projgrid,
                            covNewGridval = covGridList$covariateValues,ids=ids,noMeshPoints=noMeshPoints,
                            extraNonlinear.covGridList=extraNonlinear.covGridList)
    #saveRDS(etaAtGridListSub,file = paste(tempFolder,"/etaAtGridList_",i,".rds",sep="")) # Really no point in saving these really big files.

    # Computing the empirical mean field value at each gridpoint over the subsampled list -- strictly not need, but computed to check results.
    expContrib <- rowMeans(etaAtGridListSub)

    # Computing the squared empirical mean field value at each gridpoint over the subsampled list -- to be used for computing the mean/variance
    squaredExpContrib <- rowMeans(etaAtGridListSub^2)


    # Computing the integrated intensity, i.e. the expected number of Poisson distributed counts for each of the sampled fields in the subsamp
    areaPerSubSamp <- IntegrateVectorizedGrids(arrayGrid=etaAtGridListSub,
                                               logicalGridPointsInsideCountingDomain=gridList$logicalGridPointsInsideCountingDomain,
                                               truePixelSize=gridList$truePixelSize,
                                               scaleArea=gridList$scaleArea)

    areaPerSubSampListPhoto=lapply(predPhotoGridPointList,FUN=lapplyIntegrateVectorizedGrids,arrayGrid=etaAtGridListSub,truePixelSize=gridList$truePixelSize)
    print(paste("Finished parallelSamphandling for core",i,sep="")) # Tries to print to the terminal how far it has gotton

    ret0List <- list()
    ret0List$expContrib <- expContrib
    ret0List$squaredExpContrib <- squaredExpContrib
    ret0List$areaPerSubSamp <- areaPerSubSamp
    ret0List$areaPerSubSampListPhoto <- areaPerSubSampListPhoto
    ret0List
  }
  print(paste("Finished all handling of the original samp files in ",round((proc.time()-uu)[3])," seconds.",sep=""))

  ## Re-arranges the output from the parallelization
  expContribList <- list()
  squaredExpContribList <- list()
  areaPerSubSampVec <- NULL
  areaPerSubSampPhotoMatrix <- NULL
  for (i in 1:parallelize.noSplits){
    expContribList[[i]] <- parallelSampHandling[[i]]$expContrib
    squaredExpContribList[[i]] <- parallelSampHandling[[i]]$squaredExpContrib
    areaPerSubSampVec <- c(areaPerSubSampVec,parallelSampHandling[[i]]$areaPerSubSamp)
    areaPerSubSampPhotoMatrix <- rbind(areaPerSubSampPhotoMatrix,matrix(unlist(parallelSampHandling[[i]]$areaPerSubSampListPhoto),ncol=length(predPhotoGridPointList)))
  }

  ## Creating also the subsampled area intensity for the union of the photos
  areaPerSubSampTransectVec <- rowSums(areaPerSubSampPhotoMatrix)
  ############################# CONTINUE HERE!!!!!

  ## Computes E[X] and E[X^2] for X the posterior in each location of the grid, transform to a matrix and computes the corresponding pointwise sd of the field
  squaredExpField <- matrix(rowMeans(simplify2array(squaredExpContribList)),ncol=nxy[2],nrow=nxy[1])
  expField <- matrix(rowMeans(simplify2array(expContribList)),ncol=nxy[2],nrow=nxy[1]) # This is actually know, but "better" to compute it as it does not give NAs in sd computation below due to approx error.
  sdField <- matrix(sqrt(squaredExpField - expField^2),ncol=nxy[2],nrow=nxy[1])

  ## Computing the Poisson distribution having the sampled total intensities as the mean ##

  # Finding the evaluations points and splits them into several sub-vectors # FOR FULL AREA
  evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
  splittedevalFull=split(evalFull,ceiling((1:length(evalFull))/length(evalFull)*parallelize.noSplits))

  # Finding the evaluations points and splits them into several sub-vectors # FOR PHOTOS
  evalFullPhoto <- unique(round(seq(qpois(0.001,min(areaPerSubSampPhotoMatrix,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampPhotoMatrix,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
  splittedevalFullPhoto=split(evalFullPhoto,ceiling((1:length(evalFullPhoto))/length(evalFullPhoto)*parallelize.noSplits))

  # Finding the evaluations points and splits them into several sub-vectors # FOR TRANSECT
  evalFullTransect <- unique(round(seq(qpois(0.001,min(areaPerSubSampTransectVec,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampTransectVec,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
  splittedevalFullTransect=split(evalFullTransect,ceiling((1:length(evalFullTransect))/length(evalFullTransect)*parallelize.noSplits))



  ## First poisson evaluation for the counting domain
  doMC::registerDoMC(parallelize.numCores)
  export.var=c("areaPerSubSampVec","splittedevalFull")
  non.export.var=ls()[!(ls()%in%export.var)]
  uu=proc.time()
  posteriorDist <- foreach::foreach(i=1:length(splittedevalFull),.noexport=non.export.var,.verbose=FALSE,.inorder=TRUE) %dopar% {
    sapply(X=splittedevalFull[[i]],FUN=MeanPoissonDist,lambdavec=areaPerSubSampVec)
  }
  print(paste("Finished computing the mean of the Poisson counts in countingDomain in ",round((proc.time()-uu)[3])," seconds.",sep=""))
  posteriorDist <- unlist(posteriorDist)
  posteriorDist <- posteriorDist/sum(posteriorDist)


  # Then poisson evaluation for each photo
  doMC::registerDoMC(parallelize.numCores)
  export.var=c("areaPerSubSampPhotoMatrix","splittedevalFullPhoto")
  non.export.var=ls()[!(ls()%in%export.var)]
  uu=proc.time()
  posteriorDistPhoto <- foreach::foreach(i=1:length(splittedevalFullPhoto),.noexport=non.export.var,.verbose=FALSE,.inorder=TRUE) %dopar% {
    lapply(X=splittedevalFullPhoto[[i]],FUN=MeanPoissonDistMat,lambdaMat=areaPerSubSampPhotoMatrix)
  }
  print(paste("Finished computing the mean of the Poisson counts in each prediction photo in ",round((proc.time()-uu)[3])," seconds.",sep=""))

  posteriorDistPhoto <- matrix(unlist(posteriorDistPhoto),ncol=length(predPhotoGridPointList),byrow=T)
  posteriorDistPhoto <- posteriorDistPhoto%*%diag(1/colSums(posteriorDistPhoto))


  ## First poisson evaluation for the transect
  doMC::registerDoMC(parallelize.numCores)
  export.var=c("areaPerSubSampTransectVec","splittedevalFullTransect")
  non.export.var=ls()[!(ls()%in%export.var)]
  uu=proc.time()
  posteriorDistTransect <- foreach::foreach(i=1:length(splittedevalFullTransect),.noexport=non.export.var,.verbose=FALSE,.inorder=TRUE) %dopar% {
    sapply(X=splittedevalFullTransect[[i]],FUN=MeanPoissonDist,lambdavec=areaPerSubSampTransectVec)
  }
  print(paste("Finished computing the mean of the Poisson counts for the full transect in ",round((proc.time()-uu)[3])," seconds.",sep=""))
  posteriorDistTransect <- unlist(posteriorDistTransect)
  posteriorDistTransect <- posteriorDistTransect/sum(posteriorDistTransect)



  expFieldDomain <- expField
  expFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA
  sdFieldDomain <- sdField
  sdFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA


  retret <- list()
  retret$posteriorEvalPoints <- evalFull
  retret$posteriorDist <- posteriorDist
  retret$posteriorevalFullPhoto <- evalFullPhoto
  retret$posteriorDistPhoto <- posteriorDistPhoto
  retret$posteriorevalFullTransect <- evalFullTransect
  retret$posteriorDistTransect <- posteriorDistTransect

  retret$mean.field.samp <- expField
  retret$sd.field.samp <- sdField
  retret$mean.field.domain.samp <- expFieldDomain
  retret$sd.field.domain.samp <- sdFieldDomain

  return(retret)
  print("Finished running the PhotoPostPredDist function")

}

#' Reversed negative binomial sampler
#'
#' @param muVec vector or means
#' @param est.theta numeric theta value
#' @return one sample for every value of mu
#' @export

rnbinomSum <- function(muVec,
                       est.theta,
                       subSampPerSamp){
  samps <- matrix(rnbinom(n=length(muVec)*subSampPerSamp, size=est.theta, mu=muVec),ncol=subSampPerSamp)
  return(colSums(samps))
}

#' Reversed negative binomial sampler
#'
#' @param muVec vector or means
#' @param est.theta numeric theta value
#' @return one sample for every value of mu
#' @export

rnbinomPredPhoto <- function(muVec,
                       est.theta,
                       subSampPerSamp){
  samps <- rnbinom(n=length(muVec)*subSampPerSamp, size=est.theta, mu=muVec)
  return(samps)
}


#' Reversed poisson sampler
#'
#' @param muVec vector or means
#' @param subSampPerSamp ...
#' @return one sample for every value of mu
#' @export

rpoisSum <- function(muVec,
                     subSampPerSamp){
  fullsamps <- rpois(n=subSampPerSamp,lambda=sum(muVec))
  return(fullsamps)
}

#' Reversed poisson sampler for pred photo
#'
#' @param muVec vector or means
#' @param subSampPerSamp ...
#' @return one sample for every value of mu
#' @export

rpoisPredPhoto <- function(muVec,
                     subSampPerSamp){
  fullsamps <- rpois(n=length(muVec)*subSampPerSamp,lambda=muVec)
  return(fullsamps)
}


#' Sample samplePostPredDistGAMPoisson
#'
#' DO not bother
#'
#' @export


samplePostPredDistGAMPoissonPredPhoto <- function(arrayGrid,
                                                  subSampPerSamp,
                                                  areaPredPhotos){
  #  =NA # Inserts NA for the grid points which should not be counted
  return(apply(X = exp(arrayGrid)*areaPredPhotos,MARGIN=1,
               FUN = rpoisPredPhoto,subSampPerSamp = subSampPerSamp))
}



#' Sample postpredDistGAM
#'
#' DO not bother
#'
#' @export

samplePostPredDistGAMPredPhoto <- function(arrayGrid,
                                           est.theta,
                                           subSampPerSamp,
                                           areaPredPhotos){
  #  =NA # Inserts NA for the grid points which should not be counted
  return(apply(X = exp(arrayGrid)*areaPredPhotos,MARGIN=1,
               FUN = rnbinomPredPhoto,est.theta=est.theta,subSampPerSamp = subSampPerSamp))
}



#' Sample postpredDistGAM
#'
#' Sample from the posterior predictive distribution for the complete domain
#'
#' @param arrayGrid Matrix where the first dimension gives the complete vectorized (1-dimensional) grid ( c() applied to the matrix of grid values) and teh second dimension corresponds to differnet samples
#' @param logicalGridPointsInsideCountingDomain Logical vector, indicating which of the grid elements are within the counting domain
#' @param truePixelSize Numeric vector of dim 2 specifying the size of each pixel (in x- and y-direction)
#' @param scaleArea Numeric specifying how much the resulting integreated field should be adjust to correct for the inaccurate size of the counting domain when gridding
#' @param est.theta numeric theta value
#' @return Samples from predictive disttribution of fitted GAM model
#' @keywords inla
#' @export

samplePostPredDistGAM <- function(arrayGrid,
                                  logicalGridPointsInsideCountingDomain,
                                  truePixelSize,
                                  scaleArea,
                                  est.theta,
                                  subSampPerSamp){
  #  =NA # Inserts NA for the grid points which should not be counted
  return(apply(X = exp(arrayGrid[logicalGridPointsInsideCountingDomain,])*(truePixelSize[1]*truePixelSize[2])*scaleArea,MARGIN=2,
               FUN = rnbinomSum,est.theta=est.theta,subSampPerSamp = subSampPerSamp))
}

#' Sample postpredDistGAM
#'
#' Sample from the posterior predictive distribution for the complete domain
#'
#' @param arrayGrid Matrix where the first dimension gives the complete vectorized (1-dimensional) grid ( c() applied to the matrix of grid values) and teh second dimension corresponds to differnet samples
#' @param logicalGridPointsInsideCountingDomain Logical vector, indicating which of the grid elements are within the counting domain
#' @param truePixelSize Numeric vector of dim 2 specifying the size of each pixel (in x- and y-direction)
#' @param scaleArea Numeric specifying how much the resulting integreated field should be adjust to correct for the inaccurate size of the counting domain when gridding
#' @param est.theta numeric theta value
#' @return Samples from predictive disttribution of fitted GAM model
#' @keywords inla
#' @export

samplePostPredDistGAMPoisson <- function(arrayGrid,
                                  logicalGridPointsInsideCountingDomain,
                                  truePixelSize,
                                  scaleArea,
                                  subSampPerSamp){
  #  =NA # Inserts NA for the grid points which should not be counted
  return(apply(X = exp(arrayGrid[logicalGridPointsInsideCountingDomain,])*(truePixelSize[1]*truePixelSize[2])*scaleArea,MARGIN=2,
               FUN = rpoisSum,subSampPerSamp = subSampPerSamp))
}


#' Sample postpredDistGAM
#'
#' Sample from the posterior predictive distribution for the complete domain
#'
#' @param arrayGrid Matrix where the first dimension gives the complete vectorized (1-dimensional) grid ( c() applied to the matrix of grid values) and teh second dimension corresponds to differnet samples
#' @param logicalGridPointsInsideCountingDomain Logical vector, indicating which of the grid elements are within the counting domain
#' @param truePixelSize Numeric vector of dim 2 specifying the size of each pixel (in x- and y-direction)
#' @param scaleArea Numeric specifying how much the resulting integreated field should be adjust to correct for the inaccurate size of the counting domain when gridding
#' @return Samples from predictive disttribution of fitted GAM model
#' @keywords inla
#' @export

meanPostPredDistGAMBoth <- function(arrayGrid,
                                       logicalGridPointsInsideCountingDomain,
                                       truePixelSize,
                                       scaleArea){

  bb=arrayGrid[logicalGridPointsInsideCountingDomain,]
  cc=exp(bb)
  dd=colSums(cc)
  ee=dd*(truePixelSize[1]*truePixelSize[2])*scaleArea
  return(ee)
}





#' Wrapper for IntegrateVectorizedGrids
#'
#' Integrate over a sampled field given as a matrix within a specified counting region
#'
#' @param listi list with logicalGridPointsInsideCountingDomain and scaleArea stuff
#' @param arrayGrid Matrix where the first dimension gives the complete vectorized (1-dimensional) grid ( c() applied to the matrix of grid values) and teh second dimension corresponds to differnet samples
#' @param truePixelSize Numeric vector of dim 2 specifying the size of each pixel (in x- and y-direction)
#' @return Numeric vector with the expected number of Poisson distributed counts given each of the sampled fields
#' @keywords inla
#' @export


lapplyIntegrateVectorizedGrids <- function(listi,arrayGrid,truePixelSize){
  IntegrateVectorizedGrids(arrayGrid=arrayGrid,
                           logicalGridPointsInsideCountingDomain=listi$logical,
                           truePixelSize=truePixelSize,
                           scaleArea=listi$scaleArea)
}

#' Wrapper for samplePostPredDistGAM
#'
#' Do not bother...
#'
#' @export


lapplysamplePostPredDistGAM <- function(listi,
                                        arrayGrid,
                                        truePixelSize,
                                        allPhotosinGrid = rep(TRUE,nrow(arrayGrid)),
                                        est.theta,
                                        subSampPerSamp){
  ret=samplePostPredDistGAM(arrayGrid=arrayGrid,
                            logicalGridPointsInsideCountingDomain=listi$logical[allPhotosinGrid],
                            truePixelSize=truePixelSize,
                            scaleArea=listi$scaleArea,
                            est.theta = est.theta,
                            subSampPerSamp = subSampPerSamp)
  return(c(ret))
}


#' Wrapper for samplePostPredDistGAM
#'
#' Do not bother...
#'
#' @export


lapplysamplePostPredDistGAMPoisson <- function(listi,
                                        arrayGrid,
                                        truePixelSize,
                                        allPhotosinGrid = rep(TRUE,nrow(arrayGrid)),
                                        subSampPerSamp){
  ret=samplePostPredDistGAMPoisson(arrayGrid=arrayGrid,
                                   logicalGridPointsInsideCountingDomain=listi$logical[allPhotosinGrid],
                                   truePixelSize=truePixelSize,
                                   scaleArea=listi$scaleArea,
                                   subSampPerSamp = subSampPerSamp)
  return(c(ret))
}


#' Wrapper for samplePostPredDistGAM
#'
#' Do not bother...
#'
#' @export


lapplymeanPostPredDistGAMBoth <- function(listi,
                                          arrayGrid,
                                          truePixelSize,
                                          allPhotosinGrid = rep(TRUE,nrow(arrayGrid))){
  ret=meanPostPredDistGAMBoth(arrayGrid=arrayGrid,
                              logicalGridPointsInsideCountingDomain=listi$logical[allPhotosinGrid],
                              truePixelSize=truePixelSize,
                              scaleArea=listi$scaleArea)
  return(ret)
}



#' Integrate over sampled filed on a grid
#'
#' Integrate over a sampled field given as a matrix within a specified counting region
#'
#' @param arrayGrid Matrix where the first dimension gives the complete vectorized (1-dimensional) grid ( c() applied to the matrix of grid values) and teh second dimension corresponds to differnet samples
#' @param logicalGridPointsInsideCountingDomain Logical vector, indicating which of the grid elements are within the counting domain
#' @param truePixelSize Numeric vector of dim 2 specifying the size of each pixel (in x- and y-direction)
#' @param scaleArea Numeric specifying how much the resulting integreated field should be adjust to correct for the inaccurate size of the counting domain when gridding
#' @return Numeric vector with the expected number of Poisson distributed counts given each of the sampled fields
#' @keywords inla
#' @export

IntegrateVectorizedGrids <- function(arrayGrid,
                                 logicalGridPointsInsideCountingDomain,
                                 truePixelSize,
                                 scaleArea){
#  =NA # Inserts NA for the grid points which should not be counted
  return(colSums(exp(arrayGrid[logicalGridPointsInsideCountingDomain,]),na.rm=T)*(truePixelSize[1]*truePixelSize[2])*scaleArea) # Numerical integration over the filed assuming constant field value in each pixed. And multiplies with the scaling factor
}


#' Probability distribution computation for a mean of Poisson counts
#'
#' Computing the probability distribution for the mean of a set of Poisson distributed counts
#'
#' Note: The distribution is not normalized
#'
#' @param eval Numeric vector with points where the distribution is evaluated
#' @param lambdavec Numeric vector with lambda values (characterizing the Poisson counts)
#' @return Numeric vector with the expected number of Poisson distributed counts given each of the sampled fields
#' @keywords inla
#' @export

MeanPoissonDist <- function(eval,lambdavec){
  return(mean(dpois(eval,lambda=lambdavec)))
}

#' Reversed Probability distribution computation for a mean of Poisson counts
#'
#' Computing the probability distribution for the mean of a set of Poisson distributed counts
#'
#' Note: The distribution is not normalized
#'
#' @param eval Numeric vector with points where the distribution is evaluated
#' @param lambdavec Numeric vector with lambda values (characterizing the Poisson counts)
#' @return Numeric vector with the expected number of Poisson distributed counts given each of the sampled fields
#' @keywords inla
#' @export

MeanPoissonDistRev <- function(lambdavec,eval){
  return(mean(dpois(eval,lambda=lambdavec)))
}


#' Probability distribution computation for a mean of Poisson counts with a matrix of lambdas
#'
#' Computing the probability distribution for the mean of a set of Poisson distributed counts
#'
#' Note: The distribution is not normalized
#'
#' @param eval Numeric vector with points where the distribution is evaluated
#' @param lambdavec Numeric vector with lambda values (characterizing the Poisson counts)
#' @return Numeric vector with the expected number of Poisson distributed counts given each of the sampled fields
#' @keywords inla
#' @export

MeanPoissonDistMat <- function(eval,lambdaMat){
  return(apply(lambdaMat,MARGIN=2,FUN=MeanPoissonDistRev,eval=eval))
}


#' Probability distribution computation for a mean of Negative binomial counts
#'
#' Computing the probability distribution for the mean of a set of Negative binomial distributed counts
#'
#' Note: The distribution is not normalized
#'
#' @param eval Numeric vector with points where the distribution is evaluated
#' @param muvec Numeric vector with lambda values (characterizing the Negative binomial counts)
#' @param est.theta The fixed theta parameter
#' @return Numeric vector with the expected number of Negative binomial distributed counts given each of the sampled fields
#' @export

MeanNegbinDist <- function(eval,muvec,est.theta){
  return(mean(dnbinom(eval,mu=muvec,size = est.theta)))
}

#' Reversed probability distribution computation for a mean of Negative binomial counts
#'
#' Computing the probability distribution for the mean of a set of Negative binomial distributed counts
#'
#' Note: The distribution is not normalized
#'
#' @param muvec Numeric vector with lambda values (characterizing the Negative binomial counts)
#' @param eval Numeric vector with points where the distribution is evaluated
#' @param est.theta The fixed theta parameter
#' @return Numeric vector with the expected number of Negative binomial distributed counts given each of the sampled fields
#' @export

MeanNegbinDistRev <- function(muvec,eval,est.theta){
  return(mean(dnbinom(eval,mu=muvec,size = est.theta)))
}


#' Probability distribution computation for a mean of Negative binomial counts with a matrix of lambdas
#'
#' Computing the probability distribution for the mean of a set of Negative binomial distributed counts
#'
#' Note: The distribution is not normalized
#'
#' @param eval Numeric vector with points where the distribution is evaluated
#' @param lambdavec Numeric vector with lambda values (characterizing the Negative binomial counts)
#' @param est.theta The fixed theta parameter
#' @return Numeric vector with the expected number of Negative binomial distributed counts given each of the sampled fields
#' @keywords
#' @export

MeanNegbinDistMat <- function(eval,muMat,est.theta){
  return(apply(muMat,MARGIN=2,FUN=MeanNegbinDistRev,eval=eval,est.theta=est.theta))
}



#' Transform covariates to new grid
#'
#' Transform a covariate grid of class 'im' to a new grid of covariate values, along with only those values being within
#' a specified counting domain. Note: This function should also be used if covariates are not present to simplify
#' later computations.
#'
#' @param pp.res inla object (being the output of running the inla function)
#' @param inlaPrepList list of data provided to the inla function  (being the output from either PrepareINLAFuncContLikApprox or PrepareINLAFunc)
#' @param projgrid inla.mesh.projector object (being the output from running inla.mesh.projector)
#' @param use.covariates Logical, indicating whether covariates are used or not (see description!)
#' @param covariate.fitting String, indicating how to model covariates. "linear", quadratic (default) or "linAndLog", or FALSE for no covariates
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param covariateValues, Matrix giving the covariate values on the grid
#' @param logicalGridPointsInsideCountingDomain Logical vector, indicating which of the grid elements are within the counting domain
#' @param nxy Numberic vector of size 2, giving the dimension in x- and y-direction for the grid
#' @param extraNonlinear.covGridList List with additional projection object and such for nonlinear covariate effect, when applicable
#' @return List with several results requiring minimal interaction with inla object
#' @keywords inla
#' @export


BasicResultExtraction <- function(pp.res,
                                  inlaPrepList,
                                  projgrid,
                                  use.covariates = TRUE,
                                  covariate.fitting = "quadratic",
                                  additional.iid.term = TRUE,
                                  covariateValues,
                                  logicalGridPointsInsideCountingDomain,
                                  nxy,
                                  extraNonlinear.covGridList) {

  ## Extracting model fitting measures
  resultList <- list()
  resultList$dic <- pp.res$dic$dic
  resultList$waic <- pp.res$waic$waic
  resultList$mlik <- pp.res$mlik[1]

  spde <- inlaPrepList$spde # Gives null if spatial was FALSE when creating inlaPrepList

  spatial <- !is.null(spde)

  if(spatial){
    ## Extracting results for the underlying random field
    pp.rf <- inla.spde2.result(pp.res, 'rf', spde) # 'rf' here corresponds to the random field. Only used to estimate the range parameter below

    ## Computes the Estimates the range parameter
    resultList$mean.range.param = inla.emarginal(function(x) x,pp.rf$marginals.range.nominal[[1]]) #expectation of range = distance at which correlation is about 0.1

    resultList$rf.mean.grid <-  inla.mesh.project(projgrid,pp.res$summary.random$rf$mean) # extracts the pointwise mean of the posterior *spatial* random field
  } else {
    resultList$mean.range.param <- 0

    resultList$rf.mean.grid <-  inla.mesh.project(projgrid,rep(0,ncol(projgrid$proj$A))) # Puts 0 everywhere ( ncol(projgrid$proj$A) are the number of mesh nodes)

  }

  val.spatialX <- matrix(rep(projgrid$x,nxy[2]),ncol=nxy[2],nrow=nxy[1])
  val.spatialY <- matrix(rep(projgrid$y,nxy[1]),ncol=nxy[2],nrow=nxy[1],byrow = T)
  val.spatialXY <- sqrt(val.spatialX^2+val.spatialY^2)


  if (additional.iid.term){
    resultList$iid.mean.grid <- inla.mesh.project(projgrid,pp.res$summary.random$iid$mean)
  } else {
    resultList$iid.mean.grid <- matrix(0,ncol=nxy[2],nrow=nxy[1])
  }

  resultList$intercept.mean <- pp.res$summary.fixed$mean[1]
  resultList$covariate.mean <- 0
  resultList$covariate2.mean <- 0
  resultList$covariateLog.mean <- 0
  resultList$nonlinearCovgrid.mean <- 0
  resultList$spatialX.mean <- 0
  resultList$spatialY.mean <- 0
  resultList$spatialXY.mean <- 0


  if (use.covariates){
    if (covariate.fitting=="linear"){
      resultList$covariate.mean <- pp.res$summary.fixed$mean[2]
    }
    if (covariate.fitting=="quadratic"){
      resultList$covariate.mean <- pp.res$summary.fixed$mean[2]
      resultList$covariate2.mean <- pp.res$summary.fixed$mean[3]
    }
    if (covariate.fitting=="linAndLog"){
      resultList$covariate.mean <- pp.res$summary.fixed$mean[2]
      resultList$covariateLog.mean <- pp.res$summary.fixed$mean[3]
    }
    if (covariate.fitting=="nonlinear"){
    meanGrid <- inla.mesh.project(extraNonlinear.covGridList$projgrid.cov,pp.res$summary.random$nonlinear$mean)
    resultList$nonlinearCovgrid.mean <- extraNonlinear.covGridList$thesePoints*NA
    resultList$nonlinearCovgrid.mean[extraNonlinear.covGridList$thesePoints] <- meanGrid
    }
    if (covariate.fitting=="linearAndSpatial"){
      resultList$covariate.mean <- pp.res$summary.fixed$mean[2]
      resultList$spatialX.mean <- pp.res$summary.fixed$mean[3]
      resultList$spatialY.mean <- pp.res$summary.fixed$mean[4]
      resultList$spatialXY.mean <- pp.res$summary.fixed$mean[5]
    }
  }

  if (sum(covariateValues,na.rm = T)==0){
    logcovar  <- 0
  } else {
    logcovar <- log(covariateValues)
  }

  resultList$mean.field <- resultList$rf.mean.grid + resultList$iid.mean.grid + resultList$intercept.mean +
                           resultList$covariate.mean*covariateValues + resultList$covariate2.mean*covariateValues^2 +
                           resultList$covariateLog.mean*logcovar + resultList$nonlinearCovgrid.mean +
    resultList$spatialX.mean*val.spatialX + resultList$spatialY.mean*val.spatialY + resultList$spatialXY.mean*val.spatialXY

  # Mean field with values only within the counting domain
  resultList$mean.field.domain <- resultList$mean.field
  resultList$mean.field.domain[!logicalGridPointsInsideCountingDomain] = NA

  return(resultList)
}



#' Transform covariates to new grid
#'
#' Transform a covariate grid of class 'im' to a new grid of covariate values, along with only those values being within
#' a specified counting domain. Note: This function should also be used if covariates are not present to simplify
#' later computations.
#'
#' @param use.covariates Logical, indicating whether covariates are used or not (see description!)
#' @param covGrid, im object representing the covariate values on a dense grid where the counts live.
#' @param gridvalX Numeric vector with the x-values for the grid
#' @param gridvalY Numeric vecotr with the y-values for the grid
#' @param modelledGridVal Matrix with same dimension as the new grid, where 1 indicates that the grid point is modelled, NA means it is outside the modelling domain.
#' @param logicalGridPointsInsideCountingDomain Logical vector, indicating which of the grid elements are within the counting domain
#' @return List with the covariate values at the new grid (with respectively NA outside the modelling domain and outside the coutning domain)
#' @keywords inla
#' @import fields
#' @export

covAtNewGrid <- function(use.covariates,
                         covGrid,
                         gridvalX,
                         gridvalY,
                         modelledGridVal,
                         logicalGridPointsInsideCountingDomain){

  if (use.covariates){
    indGridX <- covGrid$xcol
    indGridY <- covGrid$yrow

    alloldGridX <- rep(indGridX,times=covGrid$dim[2])
    alloldGridY <- rep(indGridY,each=covGrid$dim[1])

    covNewGrid <- fields::as.image(Z=c(t(covGrid$v)), x = cbind(x=alloldGridX,y=alloldGridY), grid = list(x=gridvalX,y=gridvalY),na.rm=T)
    covNewGridval=modelledGridVal*covNewGrid$z
  } else {
    covNewGridval <- modelledGridVal*0
  }


covariateValuesCountingDomain <- covNewGridval
covariateValuesCountingDomain[!logicalGridPointsInsideCountingDomain] <- NA

retList <- list()
retList$covariateValues <- covNewGridval
retList$covariateValuesCountingDomain <- covariateValuesCountingDomain
return(retList)
}




#' Grid creation
#'
#' Creating a grid spanning the complete mesh with a specified pixel size and a projection from the mesh to the grid
#'
#' @param mesh Mesh object, being the output from inla.mesh.2d
#' @param countingDomain data.frame containing x- and y-coordinates of the counting domain polygon
#' @param areaCountDomain Numeric giving the area of the counting domain (in km)
#' @param grid.pixelsize Numeric, denoting the size (in km, in both x- and y-direction) of each pixel of the grid being used
#' @param gridSpan String equal to either "mesh" or "countingDomain" specifying what area the grid should span.
#' @return List with the X and Y-values for the grid, the mesh-to-grid projection object (output form inla.mesh.projector), and other related variables
#' @keywords inla
#' @import sp
#' @import INLA
#' @export

GridCreation <- function(mesh,
                         countingDomain,
                         areaCountDomain,
                         grid.pixelsize = 0.2,
                         gridSpan = "mesh"){

  ## Defines a projection from the mesh to a regular grid for sampling and plott of the posterior field
  if (gridSpan == "countingDomain"){
    rangeX <- range(countingDomain$x)
    rangeY <- range(countingDomain$y)
  } else {
    rangeX <- range(mesh$loc[,1]) # range of x-coord of grid
    rangeY <- range(mesh$loc[,2]) # range of y-coord of grid
  }


  nxy <- round(c(diff(rangeX),diff(rangeY))/grid.pixelsize) # The number of points of the grid in x- and y-direction
  projgrid <- inla.mesh.projector(mesh,xlim=rangeX,ylim=rangeY,dims=nxy) # Defines the projection
  truePixelSize <- c(projgrid$x[2]-projgrid$x[1],projgrid$y[2]-projgrid$y[1]) # The actual pixelSize

  ## Listing the x- and y-values of each grid point
  allX <- rep(projgrid$x,times=nxy[2])
  allY <- rep(projgrid$y,each=nxy[1])

  ## Checks which of the grid points which have centers within the count domain
  gridPointsInsideCountingDomain <- point.in.polygon(point.x=allX, point.y=allY,pol.x=countingDomain$x,pol.y=countingDomain$y)
  logicalGridPointsInsideCountingDomain <- as.logical(gridPointsInsideCountingDomain)
  noPixels <- sum(logicalGridPointsInsideCountingDomain)

  ## Finds the corresponding areas covered by the pixel which are counted
  pixelArea <- noPixels*(truePixelSize[1]*truePixelSize[2])

  ## How much should the estimates of the number of seals based on counting from the pixels be adjusted?
  scaleArea <- areaCountDomain/pixelArea

  modelledGridVal <- inla.mesh.project(projgrid,rep(1,mesh$n))

  retList <- list()
  retList$nxy <- nxy
  retList$gridvalX <- projgrid$x
  retList$gridvalY <- projgrid$y
  retList$truePixelSize <- truePixelSize
  retList$projgrid <- projgrid
  retList$logicalGridPointsInsideCountingDomain <- logicalGridPointsInsideCountingDomain
  retList$scaleArea <- scaleArea
  retList$modelledGridVal <- modelledGridVal
  return(retList)
}


#' Grid Points in prediction photos
#'
#' Finds which grid points are to be used to represent the value in each of the photos.
#'
#' @param gridList List with grid objects, being the output from GridCreation
#' @param predPhotos data.frame containing x- and y-coordinates of the four corners of each of the pictures
#' @return List with which gridpoints are within each of the pictures and how much their integrated intensity should be scaled to match the size of the picture
#' @keywords inla
#' @import sp
#' @export

GridPointsInPredPhotos <- function(gridList,
                                   predPhotos){

  truePixelSize <- gridList$truePixelSize
  nxy <- gridList$nxy

  ## Listing the x- and y-values of each grid point
  allX <- rep(gridList$gridvalX,times=nxy[2])
  allY <- rep(gridList$gridvalY,each=nxy[1])

  #plot(c(predPhotos[1,c(1,2,2,1,1)]),c(predPhotos[1,c(3,3,4,4,3)]),type='l'),xlim=c(-7,20),ylim=c(33.8,35))
  #points(allX,allY)
  ## Checks which of the grid points which have centers within the count domain
  gridPointsInPredPhotosList <- list()
  for (i in 1:nrow(predPhotos)){
    picArea <- abs(predPhotos[i,2]-predPhotos[i,1])*abs(predPhotos[i,4]-predPhotos[i,3])

    gridPointsInPredPhotosList[[i]] <- list()
    gridPointsInPredPhotosList[[i]]$logical <- sp::point.in.polygon(point.x=allX, point.y=allY,pol.x=predPhotos[i,c(1,2,2,1,1)],pol.y=predPhotos[i,c(3,3,4,4,3)])
    gridPointsInPredPhotosList[[i]]$logical <- as.logical(gridPointsInPredPhotosList[[i]]$logical)
    gridPointsInPredPhotosList[[i]]$noPixels <- sum(gridPointsInPredPhotosList[[i]]$logical)
    gridPointsInPredPhotosList[[i]]$pixelArea <- gridPointsInPredPhotosList[[i]]$noPixels*(truePixelSize[1]*truePixelSize[2])
    gridPointsInPredPhotosList[[i]]$scaleArea <- picArea/gridPointsInPredPhotosList[[i]]$pixelArea
 #   points(allX[gridPointsInPredPhotosList[[i]]$logical],allY[gridPointsInPredPhotosList[[i]]$logical],col=2)
  }
  return(gridPointsInPredPhotosList)
}

#' Mesh points in prediction photos
#'
#' Finds which mesh points are to be used to represent the value in each of the photos.
#'
#' @param gridList The mesh
#' @param predPhotos data.frame containing x- and y-coordinates of the four corners of each of the pictures
#' @return List with which mesh points are within each of the pictures and how much their integrated intensity should be scaled to match the size of the picture
#' @keywords inla
#' @import sp
#' @export

MeshPointsInPredPhotos <- function(mesh,
                                   predPhotos){

  #truePixelSize <- gridList$truePixelSize
  #nxy <- gridList$nxy

  ## Listing the x- and y-values of each grid point
  allX <- mesh$loc[,1]
  allY <- mesh$loc[,2]

  #plot(c(predPhotos[1,c(1,2,2,1,1)]),c(predPhotos[1,c(3,3,4,4,3)]),type='l'),xlim=c(-7,20),ylim=c(33.8,35))
  #points(allX,allY)
  ## Checks which of the grid points which have centers within the count domain
  meshPointsInPredPhotosList <- list()
  for (i in 1:nrow(predPhotos)){
    picArea <- abs(predPhotos[i,2]-predPhotos[i,1])*abs(predPhotos[i,4]-predPhotos[i,3])

    meshPointsInPredPhotosList[[i]] <- list()
    meshPointsInPredPhotosList[[i]]$logical <- sp::point.in.polygon(point.x=allX, point.y=allY,pol.x=predPhotos[i,c(1,2,2,1,1)],pol.y=predPhotos[i,c(3,3,4,4,3)])
    meshPointsInPredPhotosList[[i]]$logical <- as.logical(meshPointsInPredPhotosList[[i]]$logical)
    meshPointsInPredPhotosList[[i]]$noMeshPoints <- sum(meshPointsInPredPhotosList[[i]]$logical)
    meshPointsInPredPhotosList[[i]]$picArea <- picArea
    #   points(allX[gridPointsInPredPhotosList[[i]]$logical],allY[gridPointsInPredPhotosList[[i]]$logical],col=2)
  }
  return(meshPointsInPredPhotosList)
}




  #' Creating the grid for the nonlinear covariates
  #'
  #' Extracts the nonlinear effect of a given set of covariate values
  #'
  #' @param covMesh Mesh object representing the mesh for the covariate when applicable, being the output from inla.mesh.1d
  #' @param covariateValues, Matrix giving the covariate values on the grid
  #' @return List with the projection object and which values they represent
  #' @keywords inla
  #' @import INLA
  #' @export


covGridCreationNonlinear <- function(covMesh,
                                       covariateValues){

  nonNAcovariateValues <- covariateValues[!is.na(covariateValues)]

  projgrid.cov <- inla.mesh.projector(covMesh,loc = nonNAcovariateValues)
  thesePoints <- !is.na(covariateValues)

  retList <- list()
  retList$projgrid.cov <- projgrid.cov
  retList$thesePoints <- thesePoints
  return(retList)
}



#' Nonlinear covariate effects at grid
#'
#' Extracts the nonlinear effect of a given set of covariate values
#'
#' @param covMesh Mesh object representing the mesh for the covariate when applicable, being the output from inla.mesh.1d
#' @param covariateValues, Matrix giving the covariate values on the grid
#' @param pp.res inla object (being the output of running the inla function)
#' @return The nonlinear effect of the provided covariateValues at the same format as the covariateValues (vector or matrix)
#' @keywords inla
#' @import INLA
#' @export


covGridVal.nonlinear <- function(covMesh,
                                 covariateValues,
                                 pp.res){

  nonNAcovariateValues <- covariateValues[!is.na(covariateValues)]

  covEvalProjgrid <- inla.mesh.projector(covMesh,loc = nonNAcovariateValues)
  covEvalMeanVec <- inla.mesh.project(covEvalProjgrid,pp.res$summary.random$nonlinear$mean)

  covariateNonlinearEffects <- covariateValues*0

  covariateNonlinearEffects[!is.na(covariateValues)] <- covEvalMeanVec
  return(covariateNonlinearEffects)
}



#' Prepare for running INLA functio
#'
#' Gathers and stacks the input to the INLA function
#'
#' @param mesh Mesh object, being the output from inla.mesh.2d
#' @param covMesh Mesh object representing the mesh for the covariate when applicable, being the output from inla.mesh.1d
#' @param obsLoc Matrix with at least two columns specifying, respectively, the x- and y-coordinates of the observations
#' @param y.pp Numeric vector, indicating number of counts at each of the locations
#' @param e.pp Numeric vector, indicating the weight to be assigned to each location
#' @param covGrid, im object representing the covariate values on a dense grid where the counts live.
#' @param spatial Logical indicating whether a spatial spde model should be used
#' @param covariate.fitting String, indicating how to model covariates. "no", "linear", "quadratic" (default), "linAndLog" or "nonlinear"
#' @param covariates.type String equal to "band1" or "band2" indicating which of the bands from the satellite data should be used.
#' @param additional.iid.term Logical, indicating whether to include an additional iid (Gaussian) term in the latent field specification. FALSE is default
#' @param Matern.alpha Numeric, corresponding to the alpha parameter in the Matern covariance function (2 is the default)
#' @param INLA.theta.startval List containing the start values for the theta parameter in the INLA run. (NULL indicates that automatic start values should be used, and is the default)
#' @return List with all variables necessary to run the INLA function, in addition to the spde object
#' @keywords inla
#' @import spatstat
#' @import INLA
#' @export

PrepareINLAFunc <- function(mesh,
                            covMesh,
                            obsLoc,
                            y.pp,
                            e.pp,
                            covGrid,
                            spatial = TRUE,
                            covariate.fitting = "quadratic",
                            additional.iid.term = FALSE,
                            Matern.alpha = 2,
                            covariates.type = "band1",
                            INLA.theta.startval = NULL,
                            INLA.constr = TRUE) {

  if (spatial){
    spde <- inla.spde2.matern(mesh=mesh, alpha=Matern.alpha, constr = INLA.constr) # Constructing the Matern SPDE object
    A.pp <- inla.spde.make.A(mesh,loc=obsLoc) # Projection matrix to linear predictor for all mesh points
  } else {
    spde <- NULL
    A.pp <- NULL # Not sure if this is needed...
  }

  if(!is.null(covGrid)){
    nearestPixelObsPoints <- spatstat::nearest.pixel(obsLoc[,1], obsLoc[,2],covGrid)
    covAtObsPoints <- covGrid[Reduce('cbind', nearestPixelObsPoints)]
    covAtObsPoints2 <- covAtObsPoints^2
    covAtObsPointsLog <- log(covAtObsPoints)
  }

  if (covariate.fitting == "nonlinear"){
    covSpde <- inla.spde2.matern(mesh = covMesh, alpha = Matern.alpha)
    covA.pp <- inla.spde.make.A(mesh = covMesh, loc = covAtObsPoints)
  } else {
    covSpde <- NULL
    covA.pp <- NULL
  }


  ## Setting up the inla stack:

  direct.A.list <- 1

  if (covariate.fitting=="no"){
    direct.effects.list <- list(intercept=rep(1,length(covAtObsPoints)))
    direct.formula <- "0 + intercept"
  }
  if (covariate.fitting=="linear"){
    direct.effects.list <- list(intercept=1, covariate = covAtObsPoints)
    direct.formula <- "0 + intercept + covariate"
  }
  if (covariate.fitting=="quadratic"){
    direct.effects.list <- list(intercept=1, covariate = covAtObsPoints,
                                covariate2 = covAtObsPoints2)
    direct.formula <- "0 + intercept + covariate + covariate2"
  }
  if (covariate.fitting=="linAndLog"){
    direct.effects.list <- list(intercept=1, covariate = covAtObsPoints,
                                covariateLog = covAtObsPointsLog)
    direct.formula <- "0 + intercept + covariate + covariateLog"
  }
  if (covariate.fitting=="nonlinear"){
    direct.effects.list <- list(intercept=rep(1,length(covAtObsPoints)))
    direct.formula <- "0 + intercept"
  }
  if (covariate.fitting=="linearAndSpatial"){
    direct.effects.list <- list(intercept=1, covariate = covAtObsPoints, spatialX = obsLoc[,1],spatialY = obsLoc[,2], spatialXY = sqrt(obsLoc[,1]^2+obsLoc[,2]^2))
    direct.formula <- "0 + intercept + covariate + spatialX + spatialY + spatialXY"
  }
  if (covariate.fitting=="onlySpatial"){
    direct.effects.list <- list(intercept=1, spatialX = obsLoc[,1],spatialY = obsLoc[,2], spatialXY = sqrt(obsLoc[,1]^2+obsLoc[,2]^2))
    direct.formula <- "0 + intercept + spatialX + spatialY + spatialXY"
  }


  # Spatial variables

  if (spatial){
    spatial.formula <- "f(rf,model=spde)"
    spatial.A.list <- A.pp
    spatial.effects.list <- list(rf=1:mesh$n)
  } else {
    spatial.formula = spatial.A.list = spatial.effects.list = NULL
  }

  # Additional iid term

  if (additional.iid.term){
    iid.formula <- "f(iid,model='iid')"
    iid.A.list <- A.pp
    iid.effects.list <- list(iid=1:nrow(obsLoc))
  } else {
    iid.formula = iid.A.list = iid.effects.list = NULL
  }


  if (covariate.fitting == "nonlinear"){
    nonlinear.formula <- "f(nonlinear,model=covSpde)"
    nonlinear.A.list <- covA.pp
    nonlinear.effects.list <- list(nonlinear = 1:covSpde$f$n)
  } else {
    nonlinear.formula = nonlinear.A.list = nonlinear.effects.list = NULL
  }



  A.list <- list(direct.A.list,spatial.A.list,iid.A.list,nonlinear.A.list)
  A.list <- A.list[!sapply(A.list,is.null)]

  effects.list <- list(direct.effects.list,spatial.effects.list,iid.effects.list,nonlinear.effects.list)
  effects.list <- effects.list[!sapply(effects.list,is.null)]

  formula.list <- list(direct.formula,spatial.formula,iid.formula,nonlinear.formula)
  formula.list <- formula.list[!sapply(formula.list,is.null)]

  stk.pp <- inla.stack(data=list(y=y.pp, e=e.pp),
                       A=A.list,
                       tag='pp',
                       effects=effects.list)

  formula = as.formula(paste('y ~ ',paste(unlist(formula.list),collapse=' + '),sep=''))
  control.compute.list <- list(config = TRUE, dic = TRUE, waic = TRUE, mlik = TRUE) # control.compute=list(config = TRUE) is needed for posterior sampling of the posterior random field


  if (all(is.null(INLA.theta.startval))){
    control.mode.list <- inla.set.control.mode.default()
  } else {
    control.mode.list <- list(theta = INLA.theta.startval,restart = TRUE) # Use the start values for the (internally specified) theta parameters here
  }

  retList <- list()
  retList$stk.pp <- stk.pp
  retList$formula <- formula
  retList$control.compute.list <- control.compute.list
  retList$control.mode.list <- control.mode.list
  retList$spde <- spde
  retList$covSpde <- covSpde
  return(retList)



  #
  #
  #
  # basic.formula.string <- 'y ~ 0 + intercept +f(rf,model=spde)'
  # A.list <- list(A.pp,A.pp)
  # extractTypes <- c("intercept","rf")
  # extractF <- c("rf")
  #
  # if (!use.covariates){
  #   effects.list <- list(list(intercept=rep(1,mesh$n)), list(rf=1:mesh$n))
  # }
  # if (use.covariates){
  #   extractTypes <- c(extractTypes,"covariate")
  #   if (covariate.fitting=="linear"){
  #     effects.list <- list(list(intercept=rep(1,mesh$n), covariate = covAtMeshLoc), list(rf=1:mesh$n))
  #     additional.formula <- "covariate"
  #   }
  #   if (covariate.fitting=="quadratic"){
  #     effects.list <- list(list(intercept=rep(1,mesh$n), covariate = covAtMeshLoc, covariate2 = covAtMeshLoc2), list(rf=1:mesh$n))
  #     additional.formula <- c("covariate","covariate2")
  #   }
  #   if (covariate.fitting=="linAndLog"){
  #     effects.list <- list(list(intercept=rep(1,mesh$n), covariate = covAtMeshLoc, covariateLog = covAtMeshLocLog), list(rf=1:mesh$n))
  #     additional.formula <- c("covariate","covariateLog")
  #   }
  # }
  # if (additional.iid.term){
  #   effects.list[[length(effects.list)+1]] <- list(iid=1:mesh$n)
  #   additional.formula <- c(additional.formula,"f(iid,model='iid')")
  #   A.list[[length(A.list)+1]] <- A.pp
  #   extractTypes <- c(extractTypes,"iid")
  #   extractF <- c(extractF,"iid")
  # }
  #

}


#' Merge polygons
#'
#' Merge several polygons into their union
#'
#' @param photos Data frame with 4 columns, the x-left, x-right, y-bottom, y-top coordinates of the photos
#' @return A data frame with the merged polygons of the same format as the PolySet of the PBSmapping package
#' @keywords polygon
#' @import PBSmapping
#' @export


MergePolygons <- function(photos,plot=F){

  photoPoly <- NULL
  for (i in 1:nrow(photos)){
    photoPoly <- rbind(photoPoly,cbind(PID=i, POS=1:4, X=c(photos[i,1],photos[i,2])[c(1,2,2,1)], Y=rep(c(photos[i,3],photos[i,4]),each=2)))
  }

  mergedPoly <- PBSmapping::joinPolys(PBSmapping::as.PolySet(photoPoly),operation="UNION")
  if(plot){
    PBSmapping::plotPolys(mergedPoly)
  }
  mergedPoly <- as.data.frame(mergedPoly)
  return(mergedPoly)
}


#' Create counting domain polygon
#'
#' Creates a polygon corresponding to the area where seals might live and we are going to count
#'
#' @param transectStartCoordX Numeric vector with X-coordinate for which each transect is starting, in km
#' @param transectStartCoordY Numeric vector with Y-coordinate for which each transect is starting, in km
#' @param transectEndCoordX Numeric vector with X-coordinate for which each transect is ending, in km
#' @param transectEndCoordY Numeric vector with X-coordinate for which each transect is ending, in km
#' @param coordPhotoX Numeric vector with X-coordinate for the center point of each photo, in kilometers
#' @param coordPhotoY Numeric vector with Y-coordinate for the center point of each photo, in kilometers
#' @param photoWidth Numeric vector with the width of each photo, in km
#' @param photoHeight Numeric vector with the height of each photo, in km
#' @param transectYSpan Numeric, specifying how far transects are spanning in Y-direction, given in kilometers (for this data set this is 1.5Nm = 1.5*1.852 km)
#' @param transectXSpan Numeric, specifying how far transects are spanning in X-direction, given in kilometers
#' @param theseTransInCountDomain Vector of the transects to include in the counting domain (should be adjecent transects to fully sense), default is all (=length(transectStartCoordX))
#' @return A polygon of the counting domain given as a data.frame with the x- and y-coordinates of the complete polygon
#' @keywords polygon
#' @import PBSmapping
#' @export


CreateCountDomainPolygon <- function(transectStartCoordX,
                                     transectStartCoordY,
                                     transectEndCoordX,
                                     transectEndCoordY,
                                     coordPhotoX,
                                     coordPhotoY,
                                     photoWidth,
                                     photoHeight,
                                     transectYSpan = 1.5*1.852,
                                     transectXSpan = 0,
                                     theseTransInCountDomain = length(transectStartCoordX),
                                     cornerpointPolyInitializer = T)
                                     {

  transectCenters <- cbind(x=rowMeans(cbind(transectStartCoordX,transectEndCoordX)),
                          y=rowMeans(cbind(transectStartCoordY,transectEndCoordY))) # The center of each transect (in x and y-coordinate)
  noTransects <- length(transectStartCoordY)

  coordPhoto <- cbind(coordPhotoX,coordPhotoY)
  noPhoto <- length(coordPhotoX)

  photoTransect <- rep(NA,noPhoto) # Which transect does each photo belong to
  for (i in 1:noPhoto){
    distMat <- (t(coordPhoto[i,]-t(transectCenters)))^2
    photoTransect[i] <- which.min(distMat[,1]/(10^10)+distMat[,2])
  }

  ## Finds the coordinates for the end points of each transect
  transectCoord <- data.frame(leftPhoto = rep(NA,noTransects),rightPhoto = rep(NA,noTransects))

  # Finds ends pictures for each transect
  for (i in theseTransInCountDomain){
    theseTransects <- which(photoTransect==i)
    thisMinTransect0 <- which.min(coordPhotoX[theseTransects])
    thisMaxTransect0 <- which.max(coordPhotoX[theseTransects])
    thisMinTransect <- theseTransects[thisMinTransect0]
    thisMaxTransect <- theseTransects[thisMaxTransect0]

    transectCoord$leftPhoto[i] <- thisMinTransect
    transectCoord$rightPhoto[i] <- thisMaxTransect
  }

  # First creating a polygon which connects all cornerpoints
  poly.x <- c(rep(coordPhotoX[transectCoord$leftPhoto]-0.5*photoWidth[transectCoord$leftPhoto]-transectXSpan,each=2),
              rev(rep(coordPhotoX[transectCoord$rightPhoto]+0.5*photoWidth[transectCoord$rightPhoto] + transectXSpan,each=2)))
  poly.y <- NA*poly.x
  poly.y[seq(1,by=2,length.out=noTransects)] <- coordPhotoY[transectCoord$leftPhoto]+0.5*photoHeight[transectCoord$leftPhoto[i]] + transectYSpan
  poly.y[seq(2,by=2,length.out=noTransects)] <- coordPhotoY[transectCoord$leftPhoto]-0.5*photoHeight[transectCoord$leftPhoto[i]] - transectYSpan
  poly.y[2*noTransects+seq(1,by=2,length.out=noTransects)] <- rev(coordPhotoY[transectCoord$rightPhoto]-0.5*photoHeight[transectCoord$rightPhoto[i]] - transectYSpan)
  poly.y[2*noTransects+seq(2,by=2,length.out=noTransects)] <- rev(coordPhotoY[transectCoord$rightPhoto]+0.5*photoHeight[transectCoord$rightPhoto[i]] + transectYSpan)
  poly.x <- c(poly.x,poly.x[1])
  poly.y <- c(poly.y,poly.y[1])

  poly.x <- poly.x[!is.na(poly.x)]
  poly.y <- poly.y[!is.na(poly.y)]

  basisPoly <- data.frame(PID=rep(1, length(poly.x)), POS=1:length(poly.x),X=poly.x,Y=poly.y)


  transectPolyList <- list()
  for (i in theseTransInCountDomain){

    xleft <- coordPhotoX[transectCoord$leftPhoto[i]]-0.5*photoWidth[transectCoord$leftPhoto[i]]-transectXSpan
    xright <- coordPhotoX[transectCoord$rightPhoto[i]]+0.5*photoWidth[transectCoord$rightPhoto[i]]+transectXSpan
    ybottomL <- coordPhotoY[transectCoord$leftPhoto[i]]-0.5*photoHeight[transectCoord$leftPhoto[i]]-transectYSpan
    ybottomR <- coordPhotoY[transectCoord$rightPhoto[i]]-0.5*photoHeight[transectCoord$rightPhoto[i]]-transectYSpan
    ytopL <- coordPhotoY[transectCoord$leftPhoto[i]]+0.5*photoHeight[transectCoord$leftPhoto[i]]+transectYSpan
    ytopR <- coordPhotoY[transectCoord$rightPhoto[i]]+0.5*photoHeight[transectCoord$rightPhoto[i]]+transectYSpan

    transectPolyList[[i]] <- data.frame(PID=rep(1, 4), POS=1:4, X=c(xleft,xright)[c(1,2,2,1)], Y=c(ybottomL,ybottomR,ytopR,ytopL))
  }


  if(cornerpointPolyInitializer){
    fullPoly <- basisPoly
  } else {
    fullPoly <- transectPolyList[[1]]
  }

  for (i in theseTransInCountDomain){
    fullPoly <- PBSmapping::joinPolys(polysA=transectPolyList[[i]],polysB=fullPoly,operation="UNION")
  }

  nPointsFullPoly <- dim(fullPoly)[1]

  ret <- data.frame(x=fullPoly$X[c(1:nPointsFullPoly,1)],y=fullPoly$Y[c(1:nPointsFullPoly,1)])

  return(ret)
}



#' Mesh creation for rectangle matching
#'
#' Creates a 2D mesh with the inla.mesh.2d function which have a Voronoi triangulation which matches rectanglges given as input.
#'
#' @param rectangleCentersX Vector with X-coordinate of the centerpoint of each of the obligatory rectangles
#' @param rectangleCentersY Vector with X-coordinate of the centerpoint of each of the obligatory rectangles
#' @param rectangleWidth Vector with the width of each of the obligatory rectangles
#' @param rectangleHeight Vector with the height of each of the obligatory rectanglges
#' @param convHullVar.convex convex parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.concave concave parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param convHullVar.resolution resolution parameter of inla.nonvonvex.hull. See ?inla.nonconvex.hull for details
#' @param meshVar.max.edge max.edge parameter of inla.mesh.2d. Smaller values gives smaller triangles outside triangles. See ?inla.mesh.2d for details
#' @param meshVar.offset offset parameter of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param meshVar.cutoff vector with cutoff parameters of inla.mesh.2d. See ?inla.mesh.2d for details
#' @param y.cutoff.boundary numeric vector deciding the y-values which should be the coundary for the  different cutoff-values. NULL means the first is applied everywhere
#' @param countingDomainExpanded Expanded coutning domain used to define the inner boundary for the mesh
#' @return An INLA mesh object
#' @keywords mesh
#' @import INLA
#' @export


MeshCreationMatchingRectangles <- function(rectangleCentersX,
                                           rectangleCentersY,
                                           rectangleWidth,
                                           rectangleHeight,
                                           convHullVar.convex = -0.15,
                                           convHullVar.concave = convHullVar.convex,
                                           convHullVar.resolution = c(120,120),
                                           meshVar.max.edge = c(2,10),
                                           meshVar.offset = 6,
                                           meshVar.cutoff = 0.2,
                                           y.cutoff.boundary = NULL,
                                           countingDomainExpanded){

  obligMeshLoc <- list()
  obligMeshLoc$x=matrix(NA,ncol=9,nrow=length(rectangleCentersX))
  obligMeshLoc$y=matrix(NA,ncol=9,nrow=length(rectangleCentersX))

  transmat <- cbind(c(-1,0,1,-1,0,1,-1,0,1),c(-1,-1,-1,0,0,0,1,1,1)) # Transformation matrix defining the 9 directions

  # creates a matrix with all the different locations
  for (i in 1:9){
    obligMeshLoc$x[,i]=rectangleCentersX+transmat[i,1]*rectangleWidth
    obligMeshLoc$y[,i]=rectangleCentersY+transmat[i,2]*rectangleHeight
  }

  # Gather all the coordinates in vector form

  obligMeshLoc$x=c(obligMeshLoc$x)
  obligMeshLoc$y=c(obligMeshLoc$y)

#  set.seed(123)
#  add.noise <- cbind(runif(length(obligMeshLoc$x),-0.01,0.01),runif(length(obligMeshLoc$x),-0.01,0.01),0)

  obligMeshLoc=as.data.frame(obligMeshLoc)

#  obligMeshLoc <- obligMeshLoc + add.noise

  #### Creating the mesh based on the positions in obligMeshLoc, with cutoff to remove close duplicates ####

  # The function makes larger trinalges outside these observations automatically.
  domain.outer  <- inla.nonconvex.hull(points=cbind(rectangleCentersX,rectangleCentersY),
                                       convex=convHullVar.convex+0.03,
                                       concave=convHullVar.concave+0.03,
                                       resolution=convHullVar.resolution)

  domain.outer.final  <- inla.nonconvex.hull(points=cbind(rectangleCentersX,rectangleCentersY),
                                       convex=convHullVar.convex,
                                       concave=convHullVar.concave,
                                       resolution=convHullVar.resolution)



  domain.inner <- INLA::inla.mesh.segment(loc = as.matrix(countingDomainExpanded))



  initial.meshList <- list()
  for (i in 1:length(meshVar.cutoff)){
    initial.meshList[[i]] <- inla.mesh.2d(loc=obligMeshLoc,
                                          boundary = list(domain.inner,domain.outer),
                                          max.edge=meshVar.max.edge,
                                          offset=meshVar.offset,
                                          cutoff=meshVar.cutoff[i])

  }
  if (is.null(y.cutoff.boundary)){y.cutoff.boundary=10^30}

  y.cutoff.val <- c(-10^50,y.cutoff.boundary,10^50)

  newMeshLoc <- NULL
  for (i in 1:length(meshVar.cutoff)){
    keep.these.meshLoc <- (initial.meshList[[i]]$loc[,2]>y.cutoff.val[i] & initial.meshList[[i]]$loc[,2]<y.cutoff.val[i+1])

    newMeshLoc <- rbind(newMeshLoc,initial.meshList[[i]]$loc[keep.these.meshLoc,])
  }


  meshLocInInner <- (rowSums(newMeshLoc)%in%rowSums(domain.inner$loc))
  meshLocInOuter <- (rowSums(newMeshLoc)%in%rowSums(domain.outer$loc))

  #adjustTheseMeshLoc <- which(!as.logical(meshLocInInner+meshLocInOuter))
  outerBoundaryMeshLoc <- which(as.logical(meshLocInOuter))
  innerBoundaryMeshLoc <- which(as.logical(meshLocInInner))


  D <- as.matrix(dist(newMeshLoc[,1:2]))
  D.OuterBoundary <- D[,outerBoundaryMeshLoc]
  minDistOuterBoundary = apply(X=D.OuterBoundary,MARGIN=1,FUN=min)
  D.InnerBoundary <- D[,innerBoundaryMeshLoc]
  minDistInnerBoundary = apply(X=D.InnerBoundary,MARGIN=1,FUN=min)


  adjustTheseMeshLoc <- which(!as.logical((minDistOuterBoundary<5)+(minDistInnerBoundary<1.5)))


  set.seed(123)

  add.noise <- matrix(0,ncol=3,nrow=nrow(newMeshLoc))

  add.noise[adjustTheseMeshLoc,] <- cbind(runif(length(adjustTheseMeshLoc),-0.005,0.005),runif(length(adjustTheseMeshLoc),-0.005,0.005),0)

  newMeshLoc.adjusted <- newMeshLoc + add.noise

  newMeshLoc.adjusted <- newMeshLoc.adjusted[!as.logical(meshLocInOuter),]

  mesh <- inla.mesh.2d(loc=newMeshLoc.adjusted,
                       max.edge=meshVar.max.edge[2],
                       boundary=domain.outer.final,
                       cutoff=min(meshVar.cutoff))

  if(mesh$n!=nrow(unique(mesh$loc))){
    print("Mesh not fitted well at attempt!")
#    set.seed(123)
#    add.noise <- cbind(runif(nrow(mesh$loc),-0.01,0.01),runif(nrow(mesh$loc),-0.01,0.01),0)

    mesh <- inla.mesh.2d(loc=mesh$loc,
                         boundary=inla.mesh.boundary(mesh)[[1]],
                         max.edge=max.edge,
                         offset=offset,
                         cutoff=0)
    }

  return(mesh)
}



#### Specifies the corners of the pictures, which transect they belong to and the end points of each transect ####

#' Transform regular coordinates to SpatialPolygons
#'
#' Specify a function for transforming from regular coordinates to SpatialPolygons as used in the sp-package.
#' Copied from INLAs SPDE-tutorial
#'
#' @param coo Numeric matrix of dim 2, specifying the x- and y-coordinates to transform
#' @keywords Coordinates tranformation
#' @import sp
#' @export

coo2sp <- function(coo) {
  n <- nrow(coo)
  if (any(coo[1,]!=coo[n,]))
    coo <- coo[c(1:n,1),]
  sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coo)), '0')))
}

#### Specifies the corners of the pictures, which transect they belong to and the end points of each transect ####

#' Transform regular coordinates to SpatialPolygons
#'
#' Specify a function for transforming from regular coordinates to SpatialPolygons as used in the sp-package.
#' Copied from INLAs SPDE-tutorial
#'
#' @param coo Numeric matrix of dim 2, specifying the x- and y-coordinates to transform
#' @param ID The polygon ID each of the coordinates belongs to.
#' @keywords Coordinates tranformation
#' @import sp
#' @export

multiCoo2sp <- function(coo,ID) {
  ids <- unique(ID)
  polygonList <- list()
  for (i in 1:length(ids)){
    thiscoo <- coo[ID==ids[i],]
    n <- nrow(thiscoo)
    if (any(thiscoo[1,]!=thiscoo[n,])){
      thiscoo <- thiscoo[c(1:n,1),]
    }
    polygonList[[i]] <- sp::Polygon(thiscoo,hole=F)
  }
  ret <- sp::SpatialPolygons(list(sp::Polygons(polygonList, '0')))
  return(ret)
}


#' Computing weights based on Voronoi tessellation
#'
#' Creating a Voronoi tessellation from a set of points and calculates the area of each tile which are within a given domain
#'
#' @param locationsCoordX Numeric vector, specifying the x-coordinates of the locations where where tesslations are to be computes
#' @param locationsCoordY Numeric vector, specifying the y-coordinates of the locations where where tesslations are to be computes
#' @param observationDomain Data frame with 5 columns, where the 2 last are the coordinates of each of the polygons defining the observation domain, and the 2nd is the
#' polygon ID each coordinate belong to
#' @return A list whose first element (tiles) gives the actual tiles and second element giving the area of these intersecting with the specified domain
#' @keywords Voronoi tessallation
#' @import INLA
#' @import sp
#' @import deldir
#' @import rgeos
#' @import parallel
#' @export

CreateVoronoiTessellation <- function(locationsCoordX,
                                      locationsCoordY,
                                      observationDomain,
                                      parallelize.numCores){
  ## Builing Voronoi triangulated polygons for the complete mesh
  uu <- proc.time()
  dd <- deldir::deldir(locationsCoordX, locationsCoordY)
  tiles <- deldir::tile.list(dd)  # These tiles defines all the polygons which we are going to use in the modeling
  time <- proc.time()-uu
  print(paste("deldir call took",round(time[3])," seconds.",sep=""))

  # Creating polygon with the complete study region
  uu <- proc.time()
  pl.study <- multiCoo2sp(coo = observationDomain[,4:5], ID = observationDomain[,2]) # Defines the region which are going to contribute to the likelihood (have nonzero weight e.pp below)
  time <- proc.time()-uu
  print(paste("multiCoo2sp call took",round(time[3])," seconds.",sep=""))


  # Function for checking the area covered by the new polygons within the region with likelihood contribution
  PFunc <- function(p,pl.study) {
    pl <- coo2sp(cbind(p$x, p$y))
    if (rgeos::gIntersects(pl, pl.study)){
      return(rgeos::gArea(rgeos::gIntersection(pl, pl.study)))
    }
    else {
      return(0)
    }
  }
  uu <- proc.time()
  tileSize <- parallel::mclapply(X = tiles, FUN = PFunc, pl.study = pl.study, mc.cores = parallelize.numCores)
  time <- proc.time()-uu
  print(paste("mclapply-Pfunc call took",round(time[3])," seconds.",sep=""))

  tileSize <- unlist(tileSize)

  retList <- list()
  retList$tiles <- tiles
  retList$tileSize <- tileSize
  return(retList)
}





#### GAM functions

#' Do not both to write help function
#'
#' @export


BasicGridCreation <- function(countingDomain,
                              grid.pixelsize,
                              areaCountDomain){


  rangeX <- range(countingDomain$x)
  rangeY <- range(countingDomain$y)
  nxy <- round(c(diff(rangeX),diff(rangeY))/grid.pixelsize) # The number of points of the grid in x- and y-direction
  ## Listing the x- and y-values of each grid point

  eachX <- seq(rangeX[1],rangeX[2],length.out = nxy[1])
  eachY <- seq(rangeY[1],rangeY[2],length.out = nxy[2])

  allX <- rep(eachX,times=nxy[2])
  allY <- rep(eachY,each=nxy[1])

  ## Checks which of the grid points which have centers within the count domain
  gridPointsInsideCountingDomain <- point.in.polygon(point.x=allX, point.y=allY,pol.x=countingDomain$x,pol.y=countingDomain$y)
  logicalGridPointsInsideCountingDomain <- as.logical(gridPointsInsideCountingDomain)
  noPixels <- sum(logicalGridPointsInsideCountingDomain)

  ## Finds the corresponding areas covered by the pixel which are counted
  pixelArea <- noPixels*(grid.pixelsize^2)

  ## How much should the estimates of the number of seals based on counting from the pixels be adjusted?
  scaleArea <- areaCountDomain/pixelArea


  modelledGridVal <- matrix(NA,ncol=nxy[2],nrow=nxy[1])
  modelledGridVal[logicalGridPointsInsideCountingDomain] <- 1

  retList <- list()
  retList$nxy <- nxy
  retList$gridvalX <- eachX
  retList$gridvalY <- eachY
  retList$modelledGridVal <- modelledGridVal
  retList$scaleArea <- scaleArea
  retList$logicalGridPointsInsideCountingDomain <- logicalGridPointsInsideCountingDomain
  retList$truePixelSize=c(grid.pixelsize,grid.pixelsize)
  return(retList)
}


#' Do not both to write help function
#'
#' @export


GAMImputePred <- function(allX,
                          allY,
                          use.covariates,
                          covGrid,
                          GAMfit){

  if (use.covariates){
    nearestPixelImputePoints <- spatstat::nearest.pixel(allX, allY,covGrid)
    imputeCov <- covGrid[Reduce('cbind', nearestPixelImputePoints)]
    imputeCov2 <- imputeCov^2
    imputeCovLog <- log(imputeCov)

    imputeData <- data.frame(x=allX, y = allY, covariate = imputeCov, covariate2 = imputeCov2, covariateLog = imputeCovLog, logAreakm = 0,
                             spatialX = allX, spatialY = allY, spatialXY = sqrt(allX^2+allY^2))
  } else {
    imputeData <- data.frame(x=allX, y = allY, logAreakm = 0)
  }

  pred <- predict(GAMfit,imputeData,se.fit=T)
  imputedPred <- pred$fit
  imputedSe <- pred$se.fit


  return(list(imputedPred=imputedPred,imputedSe=imputedSe,imputeData=imputeData))
}

#' Do not both to write help function
#'
#' @export


GAMImpute <- function(gridvalX,
                      gridvalY,
                      use.covariates,
                      covGrid,
                      GAMfit,
                      logicalGridPointsInsideCountingDomain){

  nxy <- c(length(gridvalX),length(gridvalY))

  allX <- rep(gridvalX,times=nxy[2])
  allY <- rep(gridvalY,each=nxy[1])

  if (use.covariates){
    nearestPixelImputePoints <- spatstat::nearest.pixel(allX, allY,covGrid)
    imputeCov <- covGrid[Reduce('cbind', nearestPixelImputePoints)]
    imputeCov2 <- imputeCov^2
    imputeCovLog <- log(imputeCov)

    imputeData <- data.frame(x=allX, y = allY, covariate = imputeCov, covariate2 = imputeCov2, covariateLog = imputeCovLog, logAreakm = 0,
                             spatialX = allX, spatialY = allY, spatialXY = sqrt(allX^2+allY^2))

    imputeData.rf <- data.frame(x=allX, y = allY, covariate = 0, covariate2 = 0, covariateLog = 0, logAreakm = 0,
                                spatialX = 0, spatialY = 0, spatialXY = 0)

  } else {
    imputeData <- data.frame(x=allX, y = allY, logAreakm = 0)
    imputeData.rf <- imputeData
  }
  imputeData0 <- imputeData[logicalGridPointsInsideCountingDomain,]
  imputeData.rf0 <- imputeData.rf[logicalGridPointsInsideCountingDomain,]


  imputedPred <- imputedSe <- imputedPred.rf <- matrix(NA,ncol=nxy[2],nrow=nxy[1])
  pred <- predict(GAMfit,imputeData0,se.fit=T)
  pred.rf <- predict(GAMfit,imputeData.rf0,se.fit=F)


  imputedPred[logicalGridPointsInsideCountingDomain] <- pred$fit
  imputedSe[logicalGridPointsInsideCountingDomain] <- pred$se.fit

  imputedPred.rf[logicalGridPointsInsideCountingDomain] <- pred.rf


  return(list(imputedPred=imputedPred,imputedSe=imputedSe,imputeData=imputeData, imputedPred.rf = imputedPred.rf))
}


#' Do not both to write help function
#'
#' @export


GridPointsInPredPhotosGAM <- function(gridList,
                                      predPhotos,
                                      grid.pixelsize){

  nxy <- gridList$nxy

  ## Listing the x- and y-values of each grid point
  allX <- rep(gridList$gridvalX,times=nxy[2])
  allY <- rep(gridList$gridvalY,each=nxy[1])

  #  plot(c(predPhotos[1,c(1,2,2,1,1)]),c(predPhotos[1,c(3,3,4,4,3)]),type='l',xlim=c(-7,20),ylim=c(33.8,35))
  #  points(allX,allY)
  ## Checks which of the grid points which have centers within the count domain
  gridPointsInPredPhotosList <- list()
  for (i in 1:nrow(predPhotos)){
    picArea <- abs(predPhotos[i,2]-predPhotos[i,1])*abs(predPhotos[i,4]-predPhotos[i,3])

    gridPointsInPredPhotosList[[i]] <- list()
    gridPointsInPredPhotosList[[i]]$logical <- sp::point.in.polygon(point.x=allX, point.y=allY,pol.x=predPhotos[i,c(1,2,2,1,1)],pol.y=predPhotos[i,c(3,3,4,4,3)])
    gridPointsInPredPhotosList[[i]]$logical <- as.logical(gridPointsInPredPhotosList[[i]]$logical)
    gridPointsInPredPhotosList[[i]]$noPixels <- sum(gridPointsInPredPhotosList[[i]]$logical)
    gridPointsInPredPhotosList[[i]]$pixelArea <- gridPointsInPredPhotosList[[i]]$noPixels*(grid.pixelsize^2)
    gridPointsInPredPhotosList[[i]]$scaleArea <- picArea/gridPointsInPredPhotosList[[i]]$pixelArea
    #   points(allX[gridPointsInPredPhotosList[[i]]$logical],allY[gridPointsInPredPhotosList[[i]]$logical],col=2)
  }
  return(gridPointsInPredPhotosList)
}


#' Do not both to write help function
#'
#' @export


postPredDistGAMFunc <- function(GAMfit,
                                noSamp,
                                imputedList,
                                allPhotosinGrid){

  Rbeta <- MASS::mvrnorm(n = noSamp, coef(GAMfit), vcov(GAMfit))
  Xp <- predict(GAMfit, newdata = imputedList$imputeData[allPhotosinGrid,], type = "lpmatrix")
  sampInPhotoGridPoints <- Xp %*% t(Rbeta)

  return(sampInPhotoGridPoints)
}


#' Do not both to write help function
#'
#' @export

PhotoPostPredDistGAM <- function(samp,
                                 parallelize.noSplits = parallelize.numCores,
                                 parallelize.numCores,
                                 tempFolder,
                                 gridList,
                                 predPhotoGridPointList,
                                 est.theta,
                                 allPhotosinGrid,
                                 subSampPerSamp,
                                 fam) {

  ## Splits the samp object into smaller lists, saves them to disk and them delete them from RAM

  noSamp <- ncol(samp)
  splittedSamp=split(1:noSamp,ceiling((1:noSamp)/noSamp*parallelize.noSplits))

  for (i in 1:parallelize.noSplits){
    saveRDS(samp[,splittedSamp[[i]]],file = paste(tempFolder,"/samp_",i,".rds",sep=""))
    print(paste("Wrote splitted samp file ",i," of ",parallelize.noSplits," to disk",sep=""))
  }
  rm("samp",envir =sys.frame(-1))
  rm("samp","splittedSamp")
  print("Finished saving the splitted samp files to the temp-folder")


  doMC::registerDoMC(parallelize.numCores)

  export.var=c("gridList","tempFolder","predPhotoGridPointList","allPhotosinGrid","est.theta","subSampPerSamp","fam") # Not functions here
  non.export.var=ls()[!(ls()%in%export.var)]

  uu=proc.time()

  parallelSampHandling <- foreach::foreach(i=1:parallelize.noSplits,.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {

    # Reading a part of the subsampled list
    thisSubSamp <- readRDS(file = paste(tempFolder,"/samp_",i,".rds",sep=""))

    # Computing the eta for all grid points in each of the list in the subsampled list and write it to file
    etaAtGridListSub <- thisSubSamp

    # Computing the empirical mean field value at each gridpoint over the subsampled list -- strictly not need, but computed to check results.
    ##OFFV expContrib <- rowMeans(etaAtGridListSub)

    # Computing the squared empirical mean field value at each gridpoint over the subsampled list -- to be used for computing the mean/variance
    ##OFFV  squaredExpContrib <- rowMeans(etaAtGridListSub^2)


    ### ONLY FOR TH FULL VERSION
    # if(sum(leaveOutTransect)>0){
    #  areaPerSubSamp <- NULL
    #} else {
    #  areaPerSubSamp <- IntegrateVectorizedGrids(arrayGrid=etaAtGridListSub,
    #                                             logicalGridPointsInsideCountingDomain=gridList$logicalGridPointsInsideCountingDomain,
    #                                             truePixelSize=gridList$truePixelSize,
    #                                             scaleArea=gridList$scaleArea)
    #}


    # Computing the integrated intensity, i.e. the expected number of Poisson distributed counts for each of the sampled fields in the subsamp

    if (fam=="negbin"){
      areaPerSubSampListPhoto=lapply(predPhotoGridPointList,### MJ: FIXHERE
                                     FUN=lapplysamplePostPredDistGAM,
                                     arrayGrid=etaAtGridListSub,
                                     truePixelSize=gridList$truePixelSize,
                                     allPhotosinGrid=allPhotosinGrid,
                                     est.theta=est.theta,
                                     subSampPerSamp = subSampPerSamp)
    } else {
      areaPerSubSampListPhoto=lapply(predPhotoGridPointList, ### MJ: FIXHERE
                                     FUN=lapplysamplePostPredDistGAMPoisson,
                                     arrayGrid=etaAtGridListSub,
                                     truePixelSize=gridList$truePixelSize,
                                     allPhotosinGrid=allPhotosinGrid,
                                     subSampPerSamp = subSampPerSamp)
    }

    # This is OK
    areaPerSubSampListPhotoMEAN=lapply(predPhotoGridPointList,
                                       FUN=lapplymeanPostPredDistGAMBoth,
                                       arrayGrid=etaAtGridListSub,
                                       truePixelSize=gridList$truePixelSize,
                                       allPhotosinGrid=allPhotosinGrid)



    print(paste("Finished parallelSamphandling for core",i,sep="")) # Tries to print to the terminal how far it has gotton

    ret0List <- list()
    #OFFV ret0List$expContrib <- expContrib
    #OFFVret0List$squaredExpContrib <- squaredExpContrib
    #OFFVret0List$areaPerSubSamp <- areaPerSubSamp
    ret0List$areaPerSubSampListPhoto <- areaPerSubSampListPhoto
    ret0List$areaPerSubSampListPhotoMEAN <- areaPerSubSampListPhotoMEAN
    ret0List
  }
  print(paste("Finished all handling of the original samp files in ",round((proc.time()-uu)[3])," seconds.",sep=""))

  ## Re-arranges the output from the parallelization
  #OFFV  expContribList <- list()
  #OFFVsquaredExpContribList <- list()
  #OFFVareaPerSubSampVec <- NULL
  areaPerSubSampPhotoMatrix <- NULL
  areaPerSubSampPhotoMatrixMEAN <- NULL
  for (i in 1:parallelize.noSplits){
    #OFFV    expContribList[[i]] <- parallelSampHandling[[i]]$expContrib
    #OFFV    squaredExpContribList[[i]] <- parallelSampHandling[[i]]$squaredExpContrib
    #OFFV    areaPerSubSampVec <- c(areaPerSubSampVec,parallelSampHandling[[i]]$areaPerSubSamp)
    areaPerSubSampPhotoMatrix <- rbind(areaPerSubSampPhotoMatrix,matrix(unlist(parallelSampHandling[[i]]$areaPerSubSampListPhoto),ncol=length(predPhotoGridPointList)))
    areaPerSubSampPhotoMatrixMEAN <- rbind(areaPerSubSampPhotoMatrixMEAN,matrix(unlist(parallelSampHandling[[i]]$areaPerSubSampListPhotoMEAN),ncol=length(predPhotoGridPointList)))

  }


  ## Creating also the subsampled area intensity for the union of the photos
  areaPerSubSampTransectVec <- rowSums(areaPerSubSampPhotoMatrix)


  posthistTransect <- hist(areaPerSubSampTransectVec,breaks=50,plot=F)
  tab.areaPerSubSampTransectVec <- table(areaPerSubSampTransectVec)
  evalFullTransect <- as.numeric(names(tab.areaPerSubSampTransectVec))
  posteriorDistTransect <- as.numeric(tab.areaPerSubSampTransectVec)/sum(as.numeric(tab.areaPerSubSampTransectVec))

  posthistPhotoList <- list()
  evalFullPhotoList <- list()
  posteriorDistPhotoList <- list()
  for (i in 1:ncol(areaPerSubSampPhotoMatrix)){
    posthistPhotoList[[i]] <- hist(areaPerSubSampPhotoMatrix[,i],breaks=50,plot=F)
    tab.areaPerSubSampPhotoVec <- table(areaPerSubSampPhotoMatrix[,i])
    evalFullPhotoList[[i]] <- as.numeric(names(tab.areaPerSubSampPhotoVec))
    posteriorDistPhotoList[[i]] <- as.numeric(tab.areaPerSubSampPhotoVec)/sum(as.numeric(tab.areaPerSubSampPhotoVec))
  }

  # densityeval <- 512*2
  #
  # densTrans <- density(areaPerSubSampTransectVec,n=densityeval)
  # evalFullTransect <- densTrans$x
  # posteriorDistTransect <- densTrans$y/sum(densTrans$y)
  #
  # evalFullPhoto <- matrix(NA,nrow=densityeval,ncol=ncol(areaPerSubSampPhotoMatrix))
  # posteriorDistPhoto <- matrix(NA,nrow=densityeval,ncol=ncol(areaPerSubSampPhotoMatrix))
  # for (i in 1:ncol(areaPerSubSampPhotoMatrix)){
  #   densPhoto <- density(areaPerSubSampPhotoMatrix[,i],n=densityeval)
  #   evalFullPhoto[,i] <- densPhoto$x
  #   posteriorDistPhoto[,i] <- densPhoto$y/sum(densPhoto$y)
  # }

  ### ONLY FOR THE FULL VERISON
  ### Computes E[X] and E[X^2] for X the posterior in each location of the grid, transform to a matrix and computes the corresponding pointwise sd of the field
  #squaredExpField <- matrix(rowMeans(simplify2array(squaredExpContribList)),ncol=nxy[2],nrow=nxy[1])
  #expField <- matrix(rowMeans(simplify2array(expContribList)),ncol=nxy[2],nrow=nxy[1]) # This is actually know, but "better" to compute it as it does not give NAs in sd computation below due to approx error.
  #sdField <- matrix(sqrt(squaredExpField - expField^2),ncol=nxy[2],nrow=nxy[1])

  ## Computing the Poisson distribution having the sampled total intensities as the mean ##

  # Finding the evaluations points and splits them into several sub-vectors # FOR FULL AREA
  #OFFVevalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
  #OFFV splittedevalFull=split(evalFull,ceiling((1:length(evalFull))/length(evalFull)*parallelize.noSplits))

  # Finding the evaluations points and splits them into several sub-vectors # FOR PHOTOS
#  evalFullPhoto <- unique(round(seq(qnbinom(0.001,size = est.theta, mu = min(areaPerSubSampPhotoMatrix,na.rm=TRUE)),
#                                    min(qnbinom(0.999,size = est.theta, mu = max(areaPerSubSampPhotoMatrix,na.rm=TRUE)),negbin.maxEvals),
#                                    length.out =negbin.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
#  splittedevalFullPhoto=split(evalFullPhoto,ceiling((1:length(evalFullPhoto))/length(evalFullPhoto)*parallelize.noSplits))

  # Finding the evaluations points and splits them into several sub-vectors # FOR TRANSECT
#  evalFullTransect <- unique(round(seq(qnbinom(0.001,size = est.theta, mu = min(areaPerSubSampTransectVec,na.rm=TRUE)),
#                                       min(qnbinom(0.999,size = est.theta, mu = max(areaPerSubSampTransectVec,na.rm=TRUE)),negbin.maxEvals),
#                                       length.out =negbin.maxEvals)))

  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
#  splittedevalFullTransect=split(evalFullTransect,ceiling((1:length(evalFullTransect))/length(evalFullTransect)*parallelize.noSplits))



  # ## First poisson evaluation for the counting domain
  # doMC::registerDoMC(parallelize.numCores)
  # export.var=c("areaPerSubSampVec","splittedevalFull")
  # non.export.var=ls()[!(ls()%in%export.var)]
  # uu=proc.time()
  # posteriorDist <- foreach::foreach(i=1:length(splittedevalFull),.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {
  #   sapply(X=splittedevalFull[[i]],FUN=MeanPoissonDist,lambdavec=areaPerSubSampVec)
  # }
  # print(paste("Finished computing the mean of the Poisson counts in countingDomain in ",round((proc.time()-uu)[3])," seconds.",sep=""))
  # posteriorDist <- unlist(posteriorDist)
  # posteriorDist <- posteriorDist/sum(posteriorDist)

#
#   ## Then poisson evaluation for each photo
#   doMC::registerDoMC(parallelize.numCores)
#   export.var=c("areaPerSubSampPhotoMatrix","splittedevalFullPhoto","est.theta")
#   non.export.var=ls()[!(ls()%in%export.var)]
#   uu=proc.time()
#   posteriorDistPhoto <- foreach::foreach(i=1:length(splittedevalFullPhoto),.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {
#     lapply(X=splittedevalFullPhoto[[i]],FUN=MeanNegbinDistMat,muMat=areaPerSubSampPhotoMatrix,est.theta=est.theta)
#   }
#   print(paste("Finished computing the mean of the Negbin counts in each prediction photo in ",round((proc.time()-uu)[3])," seconds.",sep=""))
#
#   posteriorDistPhoto <- matrix(unlist(posteriorDistPhoto),ncol=length(predPhotoGridPointList),byrow=T)
#   posteriorDistPhoto <- posteriorDistPhoto%*%diag(1/colSums(posteriorDistPhoto))
#
#
#   ## First poisson evaluation for the transect
#   doMC::registerDoMC(parallelize.numCores)
#   export.var=c("areaPerSubSampTransectVec","splittedevalFullTransect","est.theta")
#   non.export.var=ls()[!(ls()%in%export.var)]
#   uu=proc.time()
#   posteriorDistTransect <- foreach::foreach(i=1:length(splittedevalFullTransect),.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {
#     sapply(X=splittedevalFullTransect[[i]],FUN=MeanNegbinDist,muvec=areaPerSubSampTransectVec,est.theta=est.theta)
#   }
#   print(paste("Finished computing the mean of the Negbin counts for the full transect in ",round((proc.time()-uu)[3])," seconds.",sep=""))
#   posteriorDistTransect <- unlist(posteriorDistTransect)
#   posteriorDistTransect <- posteriorDistTransect/sum(posteriorDistTransect)
#

  #OFFVexpFieldDomain <- expField
  #OFFVexpFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA
  #OFFVsdFieldDomain <- sdField
  #OFFVsdFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA


  retret <- list()
  #  retret$densTrans <- densTrans
  # retret$densPhotoList <- densPhotoList
  # retret$posteriorEvalPoints <- evalFull
  #  retret$posteriorDist <- posteriorDist
  retret$posteriorevalFullPhotoList <- evalFullPhotoList
  retret$posteriorDistPhotoList <- posteriorDistPhotoList
  retret$posteriorevalFullTransect <- evalFullTransect
  retret$posteriorDistTransect <- posteriorDistTransect
  retret$posthistTransect <- posthistTransect
  retret$posthistPhotoList <- posthistPhotoList
  retret$areaPerSubSampPhotoMatrix <- areaPerSubSampPhotoMatrixMEAN # Note that this is not the same as areaPerSubSampPhotoMatrix used internally in this function... (but the same as used for the INLA approach)

  #OFFVretret$mean.field.samp <- expField
  #OFFVretret$sd.field.samp <- sdField
  #OFFVretret$mean.field.domain.samp <- expFieldDomain
  #OFFVretret$sd.field.domain.samp <- sdFieldDomain

  return(retret)
  print("Finished running the PhotoPostPredDist function")

}

#' Do not both to write help function
#'
#' @export


FullPostPredDistGAM <- function(samp,
                                parallelize.noSplits = parallelize.numCores,
                                parallelize.numCores,
                                tempFolder,
                                gridList,
                                predPhotoGridPointList,
                                est.theta,
                                subSampPerSamp,
                                fam,
                                delete.samp = T) {

  ## Splits the samp object into smaller lists, saves them to disk and them delete them from RAM

  noSamp <- ncol(samp)
  splittedSamp=split(1:noSamp,ceiling((1:noSamp)/noSamp*parallelize.noSplits))

  for (i in 1:parallelize.noSplits){
    saveRDS(samp[,splittedSamp[[i]]],file = paste(tempFolder,"/samp_",i,".rds",sep=""))
    print(paste("Wrote splitted samp file ",i," of ",parallelize.noSplits," to disk",sep=""))
  }
  if(delete.samp){
    rm("samp",envir =sys.frame(-1))
    rm("samp","splittedSamp")
  }
  print("Finished saving the splitted samp files to the temp-folder")


  doMC::registerDoMC(parallelize.numCores)

  export.var=c("gridList","tempFolder","predPhotoGridPointList","est.theta","subSampPerSamp","fam") # Not functions here
  non.export.var=ls()[!(ls()%in%export.var)]

  uu=proc.time()

  parallelSampHandling <- foreach::foreach(i=1:parallelize.noSplits,.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {

    # Reading a part of the subsampled list
    thisSubSamp <- readRDS(file = paste(tempFolder,"/samp_",i,".rds",sep=""))

    # Computing the eta for all grid points in each of the list in the subsampled list and write it to file
    etaAtGridListSub <- thisSubSamp

    # Computing the empirical mean field value at each gridpoint over the subsampled list -- strictly not need, but computed to check results.
    expContrib <- rowMeans(etaAtGridListSub)

    # Computing the squared empirical mean field value at each gridpoint over the subsampled list -- to be used for computing the mean/variance
    squaredExpContrib <- rowMeans(etaAtGridListSub^2)

    if (fam=="negbin"){
      areaPerSubSamp <- samplePostPredDistGAM(arrayGrid=etaAtGridListSub,
                                              logicalGridPointsInsideCountingDomain=gridList$logicalGridPointsInsideCountingDomain,
                                              truePixelSize=gridList$truePixelSize,
                                              scaleArea=gridList$scaleArea,
                                              est.theta = est.theta,
                                              subSampPerSamp = subSampPerSamp)
    } else {
      areaPerSubSamp <- samplePostPredDistGAMPoisson(arrayGrid=etaAtGridListSub,
                                                     logicalGridPointsInsideCountingDomain=gridList$logicalGridPointsInsideCountingDomain,
                                                     truePixelSize=gridList$truePixelSize,
                                                     scaleArea=gridList$scaleArea,
                                                     subSampPerSamp = subSampPerSamp)
    }



    # Computing the integrated intensity, i.e. the expected number of Poisson distributed counts for each of the sampled fields in the subsamp

    print(paste("Finished parallelSamphandling for core",i,sep="")) # Tries to print to the terminal how far it has gotton

    ret0List <- list()
    ret0List$expContrib <- expContrib
    ret0List$squaredExpContrib <- squaredExpContrib
    ret0List$areaPerSubSamp <- c(areaPerSubSamp)
    ret0List
  }
  print(paste("Finished all handling of the original samp files in ",round((proc.time()-uu)[3])," seconds.",sep=""))

  ## Re-arranges the output from the parallelization
  expContribList <- list()
  squaredExpContribList <- list()
  areaPerSubSampVec <- NULL
  for (i in 1:parallelize.noSplits){
    expContribList[[i]] <- parallelSampHandling[[i]]$expContrib
    squaredExpContribList[[i]] <- parallelSampHandling[[i]]$squaredExpContrib
    areaPerSubSampVec <- c(areaPerSubSampVec,parallelSampHandling[[i]]$areaPerSubSamp)
  }

  # densityeval <- 512*2
  #
  # dens <- density(areaPerSubSampVec,n=densityeval)
  # evalFull <- dens$x
  # posteriorDist <- dens$y/sum(dens$y)

  ### Computes E[X] and E[X^2] for X the posterior in each location of the grid, transform to a matrix and computes the corresponding pointwise sd of the field
  squaredExpField <- matrix(rowMeans(simplify2array(squaredExpContribList)),ncol=gridList$nxy[2],nrow=gridList$nxy[1])
  expField <- matrix(rowMeans(simplify2array(expContribList)),ncol=gridList$nxy[2],nrow=gridList$nxy[1]) # This is actually know, but "better" to compute it as it does not give NAs in sd computation below due to approx error.
  sdField <- matrix(sqrt(squaredExpField - expField^2),ncol=gridList$nxy[2],nrow=gridList$nxy[1])

  ## Computing the Poisson distribution having the sampled total intensities as the mean ##

  posthist <- hist(areaPerSubSampVec,breaks=50,plot=F)
  tab.areaPerSubSampVec <- table(areaPerSubSampVec)
  evalFull <- as.numeric(names(tab.areaPerSubSampVec))
  posteriorDist <- as.numeric(tab.areaPerSubSampVec)/sum(as.numeric(tab.areaPerSubSampVec))

  # Finding the evaluations points and splits them into several sub-vectors # FOR FULL AREA
#  evalFull <- unique(round(seq(qnbinom(0.001,size = est.theta, mu = min(areaPerSubSampVec,na.rm=TRUE)),
#                               min(qnbinom(0.999,size= est.theta, mu = max(areaPerSubSampVec,na.rm=TRUE)),negbin.maxEvals),
#                               length.out =negbin.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
#  splittedevalFull=split(evalFull,ceiling((1:length(evalFull))/length(evalFull)*parallelize.noSplits))


  # Finding the evaluations points and splits them into several sub-vectors # FOR PHOTOS
  #  evalFullPhoto <- unique(round(seq(qpois(0.001,min(areaPerSubSampPhotoMatrix,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampPhotoMatrix,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
  #  splittedevalFullPhoto=split(evalFullPhoto,ceiling((1:length(evalFullPhoto))/length(evalFullPhoto)*parallelize.noSplits))

  # Finding the evaluations points and splits them into several sub-vectors # FOR TRANSECT
  #  evalFullTransect <- unique(round(seq(qpois(0.001,min(areaPerSubSampTransectVec,na.rm=TRUE)),min(qpois(0.999,max(areaPerSubSampTransectVec,na.rm=TRUE)),poisson.maxEvals),length.out =poisson.maxEvals)))
  # old # evalFull <- unique(round(seq(qpois(0.001,min(areaPerSubSampVec,na.rm=TRUE)),qpois(0.999,max(areaPerSubSampVec,na.rm=TRUE)),length.out =poisson.maxEvals)))
  #  splittedevalFullTransect=split(evalFullTransect,ceiling((1:length(evalFullTransect))/length(evalFullTransect)*parallelize.noSplits))



#  ## First poisson evaluation for the counting domain
#  doMC::registerDoMC(parallelize.numCores)
#  export.var=c("areaPerSubSampVec","splittedevalFull","est.theta")
#  non.export.var=ls()[!(ls()%in%export.var)]
#  uu=proc.time()
#  posteriorDist <- foreach::foreach(i=1:length(splittedevalFull),.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {
#    sapply(X=splittedevalFull[[i]],FUN=MeanNegbinDist,muvec=areaPerSubSampVec,est.theta=est.theta)
#  }
#  print(paste("Finished computing the mean of the Negbin counts in countingDomain in ",round((proc.time()-uu)[3])," seconds.",sep=""))
#  posteriorDist <- unlist(posteriorDist)
#  posteriorDist <- posteriorDist/sum(posteriorDist)



  # # Then poisson evaluation for each photo
  # doMC::registerDoMC(parallelize.numCores)
  # export.var=c("areaPerSubSampPhotoMatrix","splittedevalFullPhoto")
  # non.export.var=ls()[!(ls()%in%export.var)]
  # uu=proc.time()
  # posteriorDistPhoto <- foreach::foreach(i=1:length(splittedevalFullPhoto),.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {
  #   lapply(X=splittedevalFullPhoto[[i]],FUN=MeanPoissonDistMat,lambdaMat=areaPerSubSampPhotoMatrix)
  # }
  # print(paste("Finished computing the mean of the Poisson counts in each prediction photo in ",round((proc.time()-uu)[3])," seconds.",sep=""))
  #
  # posteriorDistPhoto <- matrix(unlist(posteriorDistPhoto),ncol=length(predPhotoGridPointList),byrow=T)
  # posteriorDistPhoto <- posteriorDistPhoto/colSums(posteriorDistPhoto)
  #
  #
  # ## First poisson evaluation for the transect
  # doMC::registerDoMC(parallelize.numCores)
  # export.var=c("areaPerSubSampTransectVec","splittedevalFullTransect")
  # non.export.var=ls()[!(ls()%in%export.var)]
  # uu=proc.time()
  # posteriorDistTransect <- foreach::foreach(i=1:length(splittedevalFullTransect),.noexport=non.export.var,.verbose=FALSE,.inorder=FALSE) %dopar% {
  #   sapply(X=splittedevalFullTransect[[i]],FUN=MeanPoissonDist,lambdavec=areaPerSubSampTransectVec)
  # }
  # print(paste("Finished computing the mean of the Poisson counts for the full transect in ",round((proc.time()-uu)[3])," seconds.",sep=""))
  # posteriorDistTransect <- unlist(posteriorDistTransect)
  # posteriorDistTransect <- posteriorDistTransect/sum(posteriorDistTransect)
  #

  expFieldDomain <- expField
  expFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA
  sdFieldDomain <- sdField
  sdFieldDomain[!gridList$logicalGridPointsInsideCountingDomain] <- NA

  retret <- list()

  retret$posteriorEvalPoints <- evalFull
  retret$posteriorDist <- posteriorDist
  retret$posteriorhist <- posthist
  retret$mean.field.samp <- expField
  retret$sd.field.samp <- sdField
  retret$mean.field.domain.samp <- expFieldDomain
  retret$sd.field.domain.samp <- sdFieldDomain

  return(retret)
  print("Finished running the PhotoPostPredDist function")

}

#' Do not both to write help function
#'
#' @export


BasicResultExtractionGAM <- function(GAMfit,
                                     imputedList,
                                     use.covariates = TRUE,
                                     covariate.fitting = "quadratic",
                                     logicalGridPointsInsideCountingDomain,
                                     gridList) {

  GAM.summary <- summary(GAMfit)

  ## Extracting model fitting measures
  resultList <- list()
  resultList$deviance <- GAMfit$deviance
  resultList$intercept.mean <- GAM.summary$p.coeff[1]
  resultList$covariate.mean <- 0
  resultList$covariate2.mean <- 0
  resultList$covariateLog.mean <- 0
  resultList$spatialX.mean <- 0
  resultList$spatialY.mean <- 0
  resultList$spatialXY.mean <- 0


  resultList$rf.mean.grid <- imputedList$imputedPred.rf

  if (use.covariates){
    if (covariate.fitting=="linear"){
      resultList$covariate.mean <- GAM.summary$p.coeff[2]
    }
    if (covariate.fitting=="quadratic"){
      resultList$covariate.mean <- GAM.summary$p.coeff[2]
      resultList$covariate2.mean <- GAM.summary$p.coeff[3]
    }
    if (covariate.fitting=="linAndLog"){
      resultList$covariate.mean <- GAM.summary$p.coeff[2]
      resultList$covariateLog.mean <- GAM.summary$p.coeff[3]
    }
    if (covariate.fitting=="linearAndSpatial"){
      resultList$covariate.mean <- GAM.summary$p.coeff[2]
      resultList$spatialX.mean <- GAM.summary$p.coeff[3]
      resultList$spatialY.mean <- GAM.summary$p.coeff[4]
      resultList$spatialXY.mean <- GAM.summary$p.coeff[5]
    }

  }

  resultList$mean.field <- imputedList$imputedPred

# I think covariates etc is already accounted for...
#  if (use.covariates){
#    resultList$mean.field <- imputedList$imputedPred + resultList$intercept.mean +
#      resultList$covariate.mean*imputedList$imputeData$covariate +
#      resultList$covariate2.mean*imputedList$imputeData$covariate2 +
#      resultList$covariateLog.mean*imputedList$imputeData$covariateLog
#  } else {
#    resultList$mean.field <- imputedList$imputedPred + resultList$intercept.mean
#  }

  # Mean field with values only within the counting domain
  resultList$mean.field.domain <- resultList$mean.field
  resultList$mean.field.domain[!logicalGridPointsInsideCountingDomain] = NA

  resultList$fullAreaBasicPred= sum(exp(imputedList$imputedPred),na.rm=TRUE)*gridList$truePixelSize[1]*gridList$truePixelSize[2]*gridList$scaleArea


  return(resultList)
}


#' Do not both to write help function
#'
#' @export

SummaryPlotFuncGAM <- function(covariatesplot = TRUE,
                               summaryplot = TRUE,
                               savingFolder,
                               sealPhotoDataFile,
                               sealTransectDataFile,
                               dataList,
                               orgPhotos,
                               modPhotos,
                               results.CI.level = 0.95,
                               gridList,
                               finalResList,
                               countingDomain,
                               logicalGridPointsInsideCountingDomain,
                               covNewGridval,
                               GAMfit,
                               sealType,
                               use.covariates,
                               covariates.type,
                               covariate.fitting,
                               grid.pixelsize,
                               parallelize.noSplits,
                               parallelize.numCores,
                               noSamp,
                               subSampPerSamp,
                               time,
                               testing,
                               comment,
                               leaveOutTransect,
                               fam){


  if (use.covariates){
    covNewGridvalDomain <- covNewGridval
    covNewGridvalDomain[!logicalGridPointsInsideCountingDomain] = NA
  } else {
    covNewGridval <- 0
    covNewGridvalDomain <- 0
  }

  fixed.effects.grid <- finalResList$intercept.mean + finalResList$covariate.mean*covNewGridval + finalResList$covariate2.mean*covNewGridval^2 + finalResList$covariateLog.mean*log(covNewGridval)
  fixed.effects.domain <- fixed.effects.grid
  fixed.effects.domain[!logicalGridPointsInsideCountingDomain] = NA



  ## Producing 4 1x2 covariate plots

  if (covariatesplot){
    par(mar=c(2,2,2,2), mgp=2:0)

    ## For covariate_effect_plot.pdf
    rangeCov <- range(covNewGridval,na.rm=T)
    covVal <- seq(rangeCov[1],rangeCov[2],length.out = 500)
    linearEffect <- finalResList$intercept.mean + finalResList$covariate.mean*covVal + finalResList$covariate2.mean*covVal^2 + finalResList$covariateLog.mean*log(covVal^2)


    pdf(file=file.path(savingFolder,"covariate_effect_plot.pdf"),width=7,height=4)
    plot(covVal,linearEffect, type='l',main="Covariate effect",xlab="Satellite image value",ylab="Log-intensity effect")
    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_grid_with_obs.pdf"),width=14,height=6)
    par(mfrow=c(1,2),oma=c(0,0,0,1.5))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=matrix(covNewGridval,ncol=gridList$nxy[2]),main="Covariate values",nlevel=200,col=topo.colors(200))
    lines(countingDomain,col=6,lwd=3)

    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)
    legend("bottomright",c("Observed","y=0","y>0"),col=c("white","black","black"),pch=c(0,0,15),bg="white")
    legend("bottomleft",c("Unmodelled","y=0","y>0"),col=c("white","red","red"),pch=c(0,0,15),bg="white")
    legend("topleft",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")


    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=matrix(fixed.effects.grid,ncol=gridList$nxy[2]),main="Mean of fixed effects (log-scale)",nlevel=200,col=topo.colors(200),ylab="")
    lines(countingDomain,col=6,lwd=3)
    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_grid_without_obs.pdf"),width=14,height=6)
    par(mfrow=c(1,2),oma=c(0,0,0,1.5))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=matrix(covNewGridval,ncol=gridList$nxy[2]),main="Covariate values",nlevel=200)
    lines(countingDomain,col=6,lwd=3)

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=matrix(fixed.effects.grid,ncol=gridList$nxy[2]),main="Mean of fixed effects (log-scale)",nlevel=200,ylab="")
    lines(countingDomain,col=6,lwd=3)
    legend("bottomright",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_count_domain_with_obs.pdf"),width=14,height=6)
    par(mfrow=c(1,2),oma=c(0,0,0,1.5))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=matrix(covNewGridvalDomain,ncol=gridList$nxy[2]),main="Covariate values",nlevel=200,col=topo.colors(200))
    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)
    legend("bottomright",c("Observed","y=0","y>0"),col=c("white","black","black"),pch=c(0,0,15),bg="white")
    legend("bottomleft",c("Unmodelled","y=0","y>0"),col=c("white","red","red"),pch=c(0,0,15),bg="white")


    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=matrix(fixed.effects.domain,ncol=gridList$nxy[2]),main="Mean of fixed effects (log-scale)",nlevel=200,col=topo.colors(200),ylab="")
    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)

    dev.off()


    pdf(file=file.path(savingFolder,"covariate_effects_count_domain_without_obs.pdf"),width=14,height=6)
    par(mfrow=c(1,2),oma=c(0,0,0,1.5))

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=matrix(covNewGridvalDomain,ncol=gridList$nxy[2]),main="Covariate values",nlevel=200)

    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=matrix(fixed.effects.domain,ncol=gridList$nxy[2]),main="Mean of fixed effects (log-scale)",nlevel=200,ylab="")

    dev.off()
  }

  ## Producing an 2x2 summary plot

  if (summaryplot){
    pdf(file=file.path(savingFolder,"results_with_data.pdf"),width=12,height=12)
    par(mfrow=c(2,2))

    ## Mean posterior field
    if (sum(leaveOutTransect)>0){
      fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$mean.field,main="Mean of latent field (log-scale)",nlevel=200,col=topo.colors(200))
    } else {
      fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$mean.field.samp,main="Mean of latent field (log-scale)",nlevel=200,col=topo.colors(200))

      }
    rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
    rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
    rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
    rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)
    lines(countingDomain,col=6,lwd=3)
    legend("toplef",c("Count domain"),col=c(6),lty=c(1),lwd=3,bg="white")
    legend("bottomright",c("Observed","y=0","y>0"),col=c("white","black","black"),pch=c(0,0,15),bg="white")
    legend("bottomleft",c("Unmodelled","y=0","y>0"),col=c("white","red","red"),pch=c(0,0,15),bg="white")





    ## Sd of posterior field
    if (sum(leaveOutTransect)>0){
      frame()
    } else {
      fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$sd.field.samp,main="Sd of latent field (log-scale)",nlevel=200,col=topo.colors(200))
      rect(orgPhotos[,1],orgPhotos[,3],orgPhotos[,2],orgPhotos[,4],border="red",lwd=0.5)
      rect(orgPhotos[dataList$org$noObsPerPhoto>0,1],orgPhotos[dataList$org$noObsPerPhoto>0,3],orgPhotos[dataList$org$noObsPerPhoto>0,2],orgPhotos[dataList$org$noObsPerPhoto>0,4],border="red",col="red",lwd=0.5)
      rect(modPhotos[,1],modPhotos[,3],modPhotos[,2],modPhotos[,4],border=1,lwd=0.5)
      rect(modPhotos[dataList$mod$noObsPerPhoto>0,1],modPhotos[dataList$mod$noObsPerPhoto>0,3],modPhotos[dataList$mod$noObsPerPhoto>0,2],modPhotos[dataList$mod$noObsPerPhoto>0,4],col=1,lwd=0.5)
      lines(countingDomain,col=6,lwd=3)

    }

    ## Posterior predictive dist
    if (sum(leaveOutTransect)>0){
      plotTo <- which.min((cumsum(finalResList$posteriorDistTransect)-0.99)^2)
      xlim <- range(finalResList$posteriorevalFullTransect[1:plotTo])
      plot(finalResList$posthistTransect,main="Posterior Predictive dist FOR PHOTOS IN TRANSECT",xlab="seals",ylab="probability",col="grey",xlim=xlim,freq=F)
    } else {
      plotTo <- which.min((cumsum(finalResList$posteriorDist)-0.99)^2)
      xlim <- range(finalResList$posteriorEvalPoints[1:plotTo])
      plot(finalResList$posteriorhist,main="Posterior Predictive dist",xlab="seals",ylab="probability",col="grey",xlim=xlim,freq=F)
    }

    lines(rep(finalResList$posteriorMean,2),c(0,1000),col=2)
    lines(rep(finalResList$posteriorMedian,2),c(0,1000),col=3)
    lines(rep(finalResList$posteriorMode,2),c(0,1000),col=4)

    legend("topright",c(paste("mean =",round(finalResList$posteriorMean,2)),
                        paste("median =",round(finalResList$posteriorMedian,2)),
                        paste("mode =",round(finalResList$posteriorMode,2)),
                        paste("IQR =",round(finalResList$posteriorIQR,2)),
                        paste(round(results.CI.level*100),"% CI = (",round(finalResList$posteriorCI[1]),",",round(finalResList$posteriorCI[2]),")")),
           lty=1,col=c(2:4,rep("white",2)))

    ## Just some parameters and variables

    # First splitting comment if it is too long:

    maxChar <- 50
    if (nchar(comment)>maxChar){
      splits <- gregexpr(pattern="=",comment)[[1]]
      splitHere <- max(splits[splits<maxChar])
      comment <- paste(substr(comment, 1, splitHere-1), "\n", substr(comment, splitHere, nchar(comment)), sep = "")
    }

    stStart <- sort(gregexpr('/',sealPhotoDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealPhotoDataFile)
    photoFrom <- substr(sealPhotoDataFile,stStart,stStop)

    stStart <- sort(gregexpr('/',sealTransectDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealTransectDataFile)
    transectsFrom <- substr(sealTransectDataFile,stStart,stStop)


    par(xpd=TRUE)
    frame()
    #if (dataType=="simulated"){
    #  text(0.5,1.15,paste("True # seals = ",sampPoisCounted$n,sep=""))
    #}
    fixed.effects.vec <- summary(GAMfit)$p.coef
    text(0.5,1.10,paste("Photos and transects from = ",photoFrom," and ",transectsFrom,sep=""))
    text(0.5,1.05,paste("SealType = ",sealType,sep=""))
    text(0.5,1.00,paste("use.covariates = ",use.covariates,sep=""))
    text(0.5,0.90,paste("covariates.type = ",covariates.type,sep=""))
    text(0.5,0.80,paste("covariate.fitting = ",covariate.fitting,sep=""))
    text(0.5,0.70,paste("fullAreaBasicPred = ", round(finalResList$fullAreaBasicPred,2),sep=""))
    text(0.5,0.60,paste("noSamp = ", noSamp,sep=""))
    text(0.5,0.55,paste("family = ", fam,sep=""))
    text(0.5,0.50,paste("subSampPerSamp = ", subSampPerSamp,sep=""))
    text(0.5,0.40,paste("grid.pixelsize =",grid.pixelsize,sep=""))
    text(0.5,0.30,paste("parallelize.numCores",parallelize.numCores,sep=""))
    text(0.5,0.20,paste("Number of posterior samples = ",noSamp,sep=""))
    text(0.5,0.15,paste("Mean of fixed effects = ", paste(round(fixed.effects.vec,4),collapse=", "),sep=""))
    text(0.5,0.10,paste("deviance = ", round(finalResList$deviance,4),sep=""))
    text(0.5,-0.05,paste("Running time = ", round(time[3]/60,2), " minutes", sep=""))
    text(0.5,-0.15,paste("comment: ",comment,sep=""))
    dev.off()



    pdf(file=file.path(savingFolder,"results.pdf"),width=12,height=12)
    par(mfrow=c(2,2))

    ## Mean posterior field
    fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$mean.field.domain,main="Mean of latent field (log-scale)",nlevel=200)

    ## Sd of posterior field
    if (sum(leaveOutTransect)>0){
      frame()
    } else {
      fields::image.plot(x=gridList$gridvalX, y=gridList$gridvalY, z=finalResList$sd.field.domain.samp,main="Sd of latent field (log-scale)",nlevel=200)

    }







    ## Posterior predictive dist
    if (sum(leaveOutTransect)>0){
      plotTo <- which.min((cumsum(finalResList$posteriorDistTransect)-0.99)^2)
      xlim <- range(finalResList$posteriorevalFullTransect[1:plotTo])
      plot(finalResList$posthistTransect,main="Posterior Predictive dist FOR PHOTOS IN TRANSECT",xlab="seals",ylab="probability",col="grey",xlim=xlim,freq=F)
    } else {
      plotTo <- which.min((cumsum(finalResList$posteriorDist)-0.99)^2)
      xlim <- range(finalResList$posteriorEvalPoints[1:plotTo])
      plot(finalResList$posteriorhist,main="Posterior Predictive dist",xlab="seals",ylab="probability",col="grey",xlim=xlim,freq=F)
    }

    lines(rep(finalResList$posteriorMean,2),c(0,1000),col=2)
    lines(rep(finalResList$posteriorMedian,2),c(0,1000),col=3)
    lines(rep(finalResList$posteriorMode,2),c(0,1000),col=4)

    legend("topright",c(paste("mean =",round(finalResList$posteriorMean,2)),
                        paste("median =",round(finalResList$posteriorMedian,2)),
                        paste("mode =",round(finalResList$posteriorMode,2)),
                        paste("IQR =",round(finalResList$posteriorIQR,2)),
                        paste(round(results.CI.level*100),"% CI = (",round(finalResList$posteriorCI[1]),",",round(finalResList$posteriorCI[2]),")")),
           lty=1,col=c(2:4,rep("white",2)))

    ## Just some parameters and variables

    # First splitting comment if it is too long:

    maxChar <- 50
    if (nchar(comment)>maxChar){
      splits <- gregexpr(pattern="=",comment)[[1]]
      splitHere <- max(splits[splits<maxChar])
      comment <- paste(substr(comment, 1, splitHere-1), "\n", substr(comment, splitHere, nchar(comment)), sep = "")
    }

    stStart <- sort(gregexpr('/',sealPhotoDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealPhotoDataFile)
    photoFrom <- substr(sealPhotoDataFile,stStart,stStop)

    stStart <- sort(gregexpr('/',sealTransectDataFile)[[1]],decreasing = T)[2]+1
    stStop <- nchar(sealTransectDataFile)
    transectsFrom <- substr(sealTransectDataFile,stStart,stStop)


    par(xpd=TRUE)
    frame()
    #if (dataType=="simulated"){
    #  text(0.5,1.15,paste("True # seals = ",sampPoisCounted$n,sep=""))
    #}
    fixed.effects.vec <- summary(GAMfit)$p.coef
    text(0.5,1.10,paste("Photos and transects from = ",photoFrom," and ",transectsFrom,sep=""))
    text(0.5,1.05,paste("SealType = ",sealType,sep=""))
    text(0.5,1.00,paste("use.covariates = ",use.covariates,sep=""))
    text(0.5,0.90,paste("covariates.type = ",covariates.type,sep=""))
    text(0.5,0.80,paste("covariate.fitting = ",covariate.fitting,sep=""))
    text(0.5,0.70,paste("fullAreaBasicPred = ", round(finalResList$fullAreaBasicPred,2),sep=""))
    text(0.5,0.60,paste("noSamp = ", noSamp,sep=""))
    text(0.5,0.55,paste("family = ", fam,sep=""))
    text(0.5,0.50,paste("subSampPerSamp = ", subSampPerSamp,sep=""))
    text(0.5,0.40,paste("grid.pixelsize =",grid.pixelsize,sep=""))
    text(0.5,0.30,paste("parallelize.numCores",parallelize.numCores,sep=""))
    text(0.5,0.20,paste("Number of posterior samples = ",noSamp,sep=""))
    text(0.5,0.15,paste("Mean of fixed effects = ", paste(round(fixed.effects.vec,4),collapse=", "),sep=""))
    text(0.5,0.10,paste("deviance = ", round(finalResList$deviance,4),sep=""))
    text(0.5,-0.05,paste("Running time = ", round(time[3]/60,2), " minutes", sep=""))
    text(0.5,-0.15,paste("comment: ",comment,sep=""))
    dev.off()

  }

  print("All plotting to file completed")

}









