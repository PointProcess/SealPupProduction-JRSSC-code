rm(list=ls())


library(SealCountINLA)
library(foreach)
library(INLA)
library(mgcv)
library(fields)
library(sp)


#setwd("~/Prosjekter/FRINATEK-Point/Git_FRINATEK-Point/CoxProcessesINLA/seals/SealCountINLA")

#source("R/BasicFunctions.R")
#source("R/PoissonRegressionFunction.R")

#
# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/REGSpatialCVFold10/"
# CVFold <- list(NA,10,0)
#
# for (i in CVFold[[2]]:1){
#   CVFold[[1]] <- i
#   add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
#   INLAPPSealsPoissonReg(resultsBaseFolder = resultsBaseFolder,
#                         sealType = "harps",
#                         noSamp=10000,
#                         additional_comment = add_comment,
#                         parallelize.numCores = 10,
#                         parallelize.noSplits = 20,
#                         covariate.fitting="linearAndSpatial",
#                         spatial=TRUE,
#                         save.data=FALSE,
#                         Matern.alpha = 2,
#                         transectAsCountDomain = FALSE,
#                         leaveOutTransect = FALSE,
#                         testing = FALSE,
#                         INLA.constr = TRUE,
#                         CVFold = CVFold,
#                         grid.pixelsize=2)
# }
#
#
#
#
# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/REGSpatialLeaveOutAsCV/"
# CVFold <- list(NA,27,1)
#
# for (i in CVFold[[2]]:1){
#   CVFold[[1]] <- i
#   add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
#   INLAPPSealsPoissonReg(resultsBaseFolder = resultsBaseFolder,
#                         sealType = "harps",
#                         noSamp=10000,
#                         additional_comment = add_comment,
#                         parallelize.numCores = 10,
#                         parallelize.noSplits = 20,
#                         covariate.fitting="linearAndSpatial",
#                         spatial=TRUE,
#                         save.data=FALSE,
#                         Matern.alpha = 2,
#                         transectAsCountDomain = FALSE,
#                         leaveOutTransect = FALSE,
#                         testing = FALSE,
#                         INLA.constr = TRUE,
#                         CVFold = CVFold,
#                         grid.pixelsize=2)
# }
#
#
#
# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/REGSpatialLeaveOutAsCV/"
# CVFold <- list(NA,27,1)
#
# for (i in 24:1){
#   CVFold[[1]] <- i
#   add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
#   INLAPPSealsPoissonReg(resultsBaseFolder = resultsBaseFolder,
#                         sealType = "hooded",
#                         noSamp=10000,
#                         additional_comment = add_comment,
#                         parallelize.numCores = 10,
#                         parallelize.noSplits = 20,
#                         covariate.fitting="linearAndSpatial",
#                         spatial=TRUE,
#                         save.data=FALSE,
#                         Matern.alpha = 2,
#                         transectAsCountDomain = FALSE,
#                         leaveOutTransect = FALSE,
#                         testing = FALSE,
#                         INLA.constr = TRUE,
#                         CVFold = CVFold,
#                         grid.pixelsize=2)
# }
#
# ##### GAM
#
# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/GAMSpatialCVFold10/"
# CVFold <- list(NA,10,0)
#
# for (i in 1:CVFold[[2]]){
#   CVFold[[1]] <- i
#   add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
#   SealCountINLA::GAMFittingFunc(resultsBaseFolder = resultsBaseFolder,
#                                 sealType = "harps",
#                                 noSamp=10000,
#                                 subSampPerSamp = 5,
#                                 additional_comment = add_comment,
#                                 parallelize.numCores = 10,
#                                 parallelize.noSplits = 20,
#                                 covariate.fitting="linearAndSpatial",
#                                 grid.pixelsize = 0,
#                                 save.data=FALSE,
#                                 transectAsCountDomain = FALSE,
#                                 leaveOutTransect = FALSE,
#                                 CVFold = CVFold)
# }

resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/GAMSpatialCVFold10/"
CVFold <- list(NA,10,0)

for (i in 1:CVFold[[2]]){
  CVFold[[1]] <- i
  add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
  SealCountINLA::GAMFittingFunc(resultsBaseFolder = resultsBaseFolder,
                                sealType = "hooded",
                                noSamp=10000,
                                subSampPerSamp = 5,
                                additional_comment = add_comment,
                                parallelize.numCores = 10,
                                parallelize.noSplits = 20,
                                covariate.fitting="linearAndSpatial",
                                grid.pixelsize = 0,
                                save.data=FALSE,
                                transectAsCountDomain = FALSE,
                                leaveOutTransect = FALSE,
                                CVFold = CVFold)
}

resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/GAMPoissonSpatialCVFold10/"
CVFold <- list(NA,10,0)

for (i in 1:CVFold[[2]]){
  CVFold[[1]] <- i
  add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
  SealCountINLA::GAMFittingFunc(resultsBaseFolder = resultsBaseFolder,
                                sealType = "harps",
                                noSamp=10000,
                                subSampPerSamp = 5,
                                additional_comment = add_comment,
                                parallelize.numCores = 10,
                                parallelize.noSplits = 20,
                                covariate.fitting="linearAndSpatial",
                                grid.pixelsize = 0,
                                save.data=FALSE,
                                transectAsCountDomain = FALSE,
                                leaveOutTransect = FALSE,
                                CVFold = CVFold,
                                fam="poisson")
}

resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/GAMPoissonSpatialCVFold10/"
CVFold <- list(NA,10,0)

for (i in 1:CVFold[[2]]){
  CVFold[[1]] <- i
  add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
  SealCountINLA::GAMFittingFunc(resultsBaseFolder = resultsBaseFolder,
                                sealType = "hooded",
                                noSamp=10000,
                                subSampPerSamp = 5,
                                additional_comment = add_comment,
                                parallelize.numCores = 10,
                                parallelize.noSplits = 20,
                                covariate.fitting="linearAndSpatial",
                                grid.pixelsize = 0,
                                save.data=FALSE,
                                transectAsCountDomain = FALSE,
                                leaveOutTransect = FALSE,
                                CVFold = CVFold,
                                fam="poisson")
}

###########################


resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/GAMSpatialLeaveOutAsCV/"
CVFold <- list(NA,27,1)

for (i in 1:CVFold[[2]]){
  CVFold[[1]] <- i
  add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
  SealCountINLA::GAMFittingFunc(resultsBaseFolder = resultsBaseFolder,
                                sealType = "harps",
                                noSamp=10000,
                                subSampPerSamp = 5,
                                additional_comment = add_comment,
                                parallelize.numCores = 10,
                                parallelize.noSplits = 20,
                                covariate.fitting="linearAndSpatial",
                                grid.pixelsize = 0,
                                save.data=FALSE,
                                transectAsCountDomain = FALSE,
                                leaveOutTransect = FALSE,
                                CVFold = CVFold)
}

resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/GAMSpatialLeaveOutAsCV/"
CVFold <- list(NA,27,1)

for (i in 1:CVFold[[2]]){
  CVFold[[1]] <- i
  add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
  SealCountINLA::GAMFittingFunc(resultsBaseFolder = resultsBaseFolder,
                                sealType = "hooded",
                                noSamp=10000,
                                subSampPerSamp = 5,
                                additional_comment = add_comment,
                                parallelize.numCores = 10,
                                parallelize.noSplits = 20,
                                covariate.fitting="linearAndSpatial",
                                grid.pixelsize = 0,
                                save.data=FALSE,
                                transectAsCountDomain = FALSE,
                                leaveOutTransect = FALSE,
                                CVFold = CVFold)
}

resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/GAMPoissonSpatialLeaveOutAsCV/"
CVFold <- list(NA,27,1)

for (i in 1:CVFold[[2]]){
  CVFold[[1]] <- i
  add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
  SealCountINLA::GAMFittingFunc(resultsBaseFolder = resultsBaseFolder,
                                sealType = "harps",
                                noSamp=10000,
                                subSampPerSamp = 5,
                                additional_comment = add_comment,
                                parallelize.numCores = 10,
                                parallelize.noSplits = 20,
                                covariate.fitting="linearAndSpatial",
                                grid.pixelsize = 0,
                                save.data=FALSE,
                                transectAsCountDomain = FALSE,
                                leaveOutTransect = FALSE,
                                CVFold = CVFold,
                                fam="poisson")
}

resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/GAMPoissonSpatialLeaveOutAsCV/"
CVFold <- list(NA,27,1)

for (i in 1:CVFold[[2]]){
  CVFold[[1]] <- i
  add_comment <- paste("CVFold",CVFold[[2]],"_",i,sep="")
  SealCountINLA::GAMFittingFunc(resultsBaseFolder = resultsBaseFolder,
                                sealType = "hooded",
                                noSamp=10000,
                                subSampPerSamp = 5,
                                additional_comment = add_comment,
                                parallelize.numCores = 10,
                                parallelize.noSplits = 20,
                                covariate.fitting="linearAndSpatial",
                                grid.pixelsize = 0,
                                save.data=FALSE,
                                transectAsCountDomain = FALSE,
                                leaveOutTransect = FALSE,
                                CVFold = CVFold,
                                fam="poisson")
}















#
#
# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/REGLeaveOutLinear_sumtozero/"
#
# for (i in 1:27){
#   add_comment <- paste("leaveOutTrans_",i,sep="")
#   INLAPPSealsPoissonReg(resultsBaseFolder = resultsBaseFolder,
#                         sealType = "hooded",
#                         noSamp=10000,
#                         additional_comment = paste("leaveOutTrans_",i,sep=""),
#                         parallelize.numCores = 10,
#                         parallelize.noSplits = 20,
#                         covariate.fitting="linear",
#                         grid.pixelsize = 0.04,
#                         spatial=TRUE,
#                         save.data=FALSE,
#                         Matern.alpha = 2,
#                         transectAsCountDomain = TRUE,
#                         leaveOutTransect = i,
#                         testing = FALSE,
#                         INLA.constr = TRUE)
# }
#
#
#
# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/REGLeaveOutLinearNewMeshAtLeaveOutTrans_sumtozero/"
#
# for (i in 1:27){
#   add_comment <- paste("leaveOutTrans_",i,sep="")
#   INLAPPSealsPoissonReg(resultsBaseFolder = resultsBaseFolder,
#                         sealType = "hooded",
#                         noSamp=10000,
#                         additional_comment = paste("leaveOutTrans_",i,sep=""),
#                         parallelize.numCores = 10,
#                         parallelize.noSplits = 20,
#                         covariate.fitting="linear",
#                         grid.pixelsize = 0.04,
#                         spatial=TRUE,
#                         save.data=FALSE,
#                         Matern.alpha = 2,
#                         transectAsCountDomain = TRUE,
#                         leaveOutTransect = i,
#                         testing = TRUE,
#                         INLA.constr = TRUE)
# }
#
#
#
#
#
#
# #
#
#
#
#
#
#
# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/REGLeaveOutLinearNewMeshAtLeaveOutTrans_meshVar.max.edge_1_10/"
#
# for (i in 1:27){
#   add_comment <- paste("leaveOutTrans_",i,sep="")
#   INLAPPSealsPoissonReg(resultsBaseFolder = resultsBaseFolder,
#                         sealType = "hooded",
#                         noSamp=10000,
#                         additional_comment = paste("leaveOutTrans_",i,sep=""),
#                         parallelize.numCores = 10,
#                         parallelize.noSplits = 20,
#                         covariate.fitting="linear",
#                         grid.pixelsize = 0.04,
#                         spatial=TRUE,
#                         save.data=FALSE,
#                         Matern.alpha = 2,
#                         transectAsCountDomain = TRUE,
#                         leaveOutTransect = i,
#                         testing = TRUE,
#                         meshVar.max.edge = c(1,10))
# }
#
#
# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/REGLeaveOutLinear_meshVar.max.edge_1_10/"
#
# for (i in 1:27){
#   add_comment <- paste("leaveOutTrans_",i,sep="")
#   INLAPPSealsPoissonReg(resultsBaseFolder = resultsBaseFolder,
#                         sealType = "hooded",
#                         noSamp=10000,
#                         additional_comment = paste("leaveOutTrans_",i,sep=""),
#                         parallelize.numCores = 10,
#                         parallelize.noSplits = 20,
#                         covariate.fitting="linear",
#                         grid.pixelsize = 0.04,
#                         spatial=TRUE,
#                         save.data=FALSE,
#                         Matern.alpha = 2,
#                         transectAsCountDomain = TRUE,
#                         leaveOutTransect = i,
#                         testing = FALSE,
#                         meshVar.max.edge = c(1,10))
# }
#
#
#



# resultsBaseFolder <- "/nr/samba/user/jullum/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/test/"
# i=12
#
# add_comment <- "testing_new_mesh_at_prediction_transect"
# #add_comment <- "just_testing"
#
# INLAPPSealsPoissonReg(resultsBaseFolder = resultsBaseFolder,
#                       sealType = "hooded",
#                       noSamp=10000,
#                       additional_comment = add_comment,
#                       parallelize.numCores = 10,
#                       parallelize.noSplits = 20,
#                       covariate.fitting="linear",
#                       grid.pixelsize = 0.04,
#                       spatial=TRUE,
#                       save.data=FALSE,
#                       Matern.alpha = 2,
#                       leaveOutTransect = i,
#                       transectAsCountDomain = TRUE,
#                       testing=TRUE)
#
