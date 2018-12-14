### ------------------------------------------------------------------------ ###
### jFun: fleet dynamics/behaviour ####
### ------------------------------------------------------------------------ ###
### 


jFun <- function(...){
  
  args <- list(...)
  
  ### if no method specified, return nothing
  if (is.null(args$method)) {
    
    return(NULL)
    
  } else {
    
    method <- args$method
    
  }
  
  args$method <- NULL
  
  ### Check inputs
  if (!is(args$ctrl, "fwdControl")) stop("ctrl argument must be fwdControl")
  
  ### dispatch
  out <- do.call(method, args)
  
  ### check outputs
  if (!is(out$ctrl, "fwdControl")) {
    stop("jFun must return an object of class fwdControl")
  }
  
  # return
  return(out)
  
}


### no method defined yet


