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
h_ICES <- function(stk0, refpts, ay, it, tracking, 
                   option = "A", ### WKNSMSE options
                   ...) {
  
  ### target F
  ### reduce F if SSB below Btrigger
  
  ### get reference values
  Ftrgt <- propagate(FLPar(refpts["Ftrgt"]), it)
  Btrigger <- propagate(FLPar(refpts["Btrigger"]), it)
  if ("Blim" %in% dimnames(refpts)$params) {
    Blim <- propagate(FLPar(refpts["Blim"]), it)
  } else {
    Blim <- propagate(FLPar(0), it)
  } 
  
  ### SSB status (last year only) relative to Btrigger
  status_Btrigger <- tail(ssb(stk0), 1) / Btrigger 
  ### SSB status (last year only) relative to Blim
  status_Blim <- tail(ssb(stk0), 1) / Blim
  
  ### positions (iterations) where SSB is below Btrigger
  pos_Btrigger <- which(status_Btrigger < 1)
  ### below Blim
  pos_Blim <- which(status_Blim < 1)
  
  ### default ICES HCR (option A):
  ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
  if (option == "A") {
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    
  } else if (option == "B") {
  ### option B:
  ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
  ### if SSB<Blim => F = Ftrgt * 0.25
    ### SSB < Btrigger
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    ### set to 0.25 if SSB < Blim
    mult[,,,,, pos_Blim] <- 0.25
    
  } else if (option == "C") {
  ### option C:
  ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
  ### if SSB<Blim => F = max(Ftrgt * 0.25, Ftrgt * SSB/Btrigger)
    ### SSB < Btrigger
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    ### if SSB < Blim, use maximum of SSB/Btrigger or 0.25
    ### i.e. limit ratio to 0.25
    mult[,,,,, pos_Blim] <- c(ifelse(mult[,,,,, pos_Blim] < 0.25, 0.25, 
                                     mult[,,,,, pos_Blim]))
    
  } else {
  ### unknown options
    mult <- 1
  }
  
  ### new target
  Ftrgt <- Ftrgt * mult
  
  ### create ctrl object
  ctrl <- getCtrl(values = Ftrgt, quantity = "f", years = ay + 1, it = it)
  
  ### save in tracking
  tracking["advice", ac(ay)] <- ctrl@trgtArray[ac(ay + 1), "val", ]
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}
