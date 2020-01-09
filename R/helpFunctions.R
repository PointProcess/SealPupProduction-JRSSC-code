#' Default function values to global variables
#'
#' Function for setting default function arguments as global variables
#'
#' @param fun The function, whose defualt arguments are to be set as global arguments
#' @param noDefaultVal Which "value" to use for the arguments without a default value
#' @return Creates variables with the default values as global variables
#' @export

DefToGlob <- function(fun,noDefaultVal = NULL){
  formalsFun <- formals(fun)
  formalNames <- names(formalsFun)
  nargs <- length(formalNames)

  for (i in 1:nargs){
    this <- eval(parse(text=paste("formalsFun$",formalNames[i],sep="")))
    this <- tryCatch(eval(this), error = function(x) noDefaultVal)
    assign(x=formalNames[i],value = this, envir = .GlobalEnv)
  }
  }

## Defines function for transfering longitude-latitude to x and y in km.
#'
#' Contributed by Tor Arne Øigård
#'
#' @param lon longitude
#' @param lat latitude
#' @param lon0 longitude base value
#' @param lat0 latitude base value
#' @export
lb2xykm <- function(lon,lat,lon0,lat0){
  n = length(lon);
  l0 = lon0/180*pi;
  b0 = lat0/180*pi;
  R=6360/1.852;
  l = lon/180*pi;
  b = lat/180*pi;

  x=array(0,length(l));
  y=array(0,length(b));
  cb0=cos(b0);
  sb0=sin(b0);
  A = rbind(c(0,1,0),c(-sb0,0,cb0),c(cb0,0,sb0))
  cl=cos(l-l0);
  sl=sin(l-l0);
  cb=cos(b);
  sb=sin(b);
  xg=rbind(cb*cl,cb*sl,sb)
  xp=A%*%xg;
  x=xp[1,]*R;
  y=xp[2,]*R;
  xkm=x*1.852
  ykm=y*1.852
  NM = data.frame(x=xkm,y=ykm)
  return(NM)
}


