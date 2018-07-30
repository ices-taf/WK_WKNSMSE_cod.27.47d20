### ------------------------------------------------------------------------ ###
### hFun: harvest control rule ####
### ------------------------------------------------------------------------ ###
### set target according to HCR


hFun <- function(...){
  
  args <- list(...)
  
  ### if no method specified, return nothing
  if (is.null(args$method)) {
    
    return(NULL)
    
  } else {
    method <- args$method
  }
  
  args$method <- NULL
  
  # Check inputs
  if (!is(args$stk0,"FLS")) stop("stk argument must be an FLStock")
  
  # dispatch
  ctrl <- do.call(method, args)
  
  # check outputs
  if (!is(ctrl$ctrl, "fwdControl")) {
    stop("The HCR must return and object of class fwdControl")
  }
  
  # return
  ctrl  
  
}


### ICES default HCR
### target Ftrgt
### if SSB < Btrigger, target reduced: Ftrgt * SSB/Btrigger
h_ICES <- function(stk0, refpts, ay, it, tracking, ...) {
  
  ### target F
  ### reduce F if SSB below Btrigger
  
  ### F target
  Ftrgt <- refpts["Ftrgt"]
  
  ### SSB status (last year only) relative to Btrigger
  status <- tail(ssb(stk0)) / refpts["Btrigger"]
  
  ### keep only values <1
  status <- ifelse(status < 1, status, 1)
  
  ### new F target: Ftrgt reduced if SSB below Btrigger
  Ftrgt <- Ftrgt * status
  
  ### create ctrl object
  ctrl <- getCtrl(values = Ftrgt, quantity = "f", years = ay + 1, it = it)
  
  ### save in tracking
  tracking["advice", ac(ay)] <- ctrl@trgtArray[ac(ay + 1), "val", ]
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}