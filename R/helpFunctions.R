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

