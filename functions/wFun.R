### ------------------------------------------------------------------------ ###
### wFun: technical measures ####
### ------------------------------------------------------------------------ ###
### implement technical measures, e.g. TAC constraints


wFun <- function(...){
  
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
    stop("wFun must return an object of class fwdControl")
  }
  
  # return
  return(out)
  
}


TAC_constraint <- function(ctrl, tracking, ay, upper = Inf, lower = -Inf,
                           ...) {
  
  ### check if catch targeted
  if (ctrl@target$quantity != "catch") stop("target has to be catch") 
  
  ### target year
  yr_target <- ctrl@target$year
  
  ### get catch target
  catch_target <- ctrl@trgtArray[, "val", ]
  
  ### get previous target
  catch_prev <- tracking["Implementation", ac(yr_target - 1), drop = TRUE]
  
  ### change in advice in %
  change <- (catch_target / catch_prev) * 100
  
  ### limit changes
  changes_new <- change
  changes_new <- ifelse(changes_new > upper, upper, changes_new) ### upper limit
  changes_new <- ifelse(changes_new < lower, lower, changes_new) ### lower limit
  
  ### find positions which exceed limits
  pos <- which(changes_new != change)
  
  ### modify them
  catch_target[pos] <- catch_prev[pos] * changes_new[pos]/100
  
  ### save new targets
  tracking["Implementation", ac(yr_target), ] <- catch_target
  ctrl@trgtArray[, "val", ] <- catch_target
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}
