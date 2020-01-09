# Counts per transect as contributed by Tor Arne Øigård

processed_folder <- here::here("inst","extdata", "processed")
original_folder <- here::here("inst","extdata", "original")

outfile <- "OigardTransectCounts.rds"

noTransects <- 27
tab1 <- matrix(c(1, 72.01, 16.52, 17.29, 75, 14, 1,
                 2, 71.58, 17.04, 16.51, 25, 0, 0,
                 3,71.55, 16.27, 17.43, 150, 1, 0,
                 4, 71.52, 17.48, 16.15, 188, 132,6,
                 5, 71.49, 16.35, 17.48, 146, 62,10,
                 6, 71.46, 17.44, 16.47, 118, 69,7,
                 7, 71.43, 16.53, 17.47, 111, 238,7,
                 8, 71.40, 17.51, 16.54, 119, 271,10,
                 9, 71.37, 16.38, 17.55, 155, 453,56,
                 10, 71.34, 17.45, 16.37, 140, 955,75,
                 11, 71.31, 16.36, 17.42, 144, 343,164,
                 12, 71.28, 17.50, 16.36, 157, 329,104,
                 13, 71.25, 16.33, 17.54, 169, 88,67,
                 14, 71.22, 18.22, 16.49, 195, 242,88,
                 15, 71.19, 18.29, 17.30, 121, 740,19,
                 16, 71.16, 18.29, 17.53, 76 ,394,22,
                 17, 71.13, 17.57, 18.24, 54, 136,7,
                 18, 71.10, 18.24, 17.51 ,67, 213,12,
                 19, 71.07, 17.56, 18.28 ,68, 74,8,
                 20, 71.04, 18.25, 17.57, 61, 158,19,
                 21, 71.01, 18.01, 18.36, 75, 117,20,
                 22, 70.58, 18.33, 18.09, 50, 86,11,
                 23, 70.55, 18.13, 18.46, 76, 116,37,
                 24, 70.52, 18.04, 18.38, 73, 407,15,
                 25, 70.49, 18.34, 18.00, 75, 309,11,
                 26, 70.46, 18.17, 18.38, 58, 87,1,
                 27, 70.43, 18.31, 18.15, 46, 0,0),ncol=7,byrow=TRUE)

# Converting minutes to decimals for lon and lat
tab1[,2]=floor(tab1[,2])+(tab1[,2]-floor(tab1[,2]))/60*100
tab1[,3]=floor(tab1[,3])+(tab1[,3]-floor(tab1[,3]))/60*100
tab1[,4]=floor(tab1[,4])+(tab1[,4]-floor(tab1[,4]))/60*100

tab1[,3]=-tab1[,3]
tab1[,4]=-tab1[,4]


colnames(tab1) <- c("","lat","lonStart","lonEnd","noPhotos","noHarps","noHooded")
tab1Df <- as.data.frame(tab1[,-1])

dir.create(processed_folder,recursive = T)

saveRDS(tab1Df,file = file.path(processed_folder,outfile))
