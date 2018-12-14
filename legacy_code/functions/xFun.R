### ------------------------------------------------------------------------ ###
### xFun: harvest control rule ####
### ------------------------------------------------------------------------ ###

xFun <- function(...){
  
  args <- list(...)
  
  ### if no method specified, return nothing
  if (is.null(args$method)) {
    
    return(list(ctrl = NULL, refpts = NULL))
    
  } else {
    method <- args$method
  }
  
  args$method <- NULL
  
  # Check inputs
  #if (!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
  
  # dispatch
  ctrl <- do.call(method, args)
  
  # check outputs
  if (!is(ctrl$refpts, "FLPar")) {
    stop("The HCR parametrization must return and object of class FLPar")
  }
  
  # return
  ctrl  
  
}


### turn reference points into FLPar
x_ICES <- function(...) {
  
  args <- list(...)
  
  ### coerce into FLPar
  refpts <- FLPar(args)
  
  return(list(refpts = refpts))
  
}