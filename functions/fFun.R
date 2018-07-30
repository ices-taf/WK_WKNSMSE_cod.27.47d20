### ------------------------------------------------------------------------ ###
### fFun: assessment ####
### ------------------------------------------------------------------------ ###

fFun <- function(...){#browser()
  # in: stk = FLStock, idx = FLIndices
  #	method = character for the wrapper function on the assessment.
  # out: list with elements 'stk' updated with the assessment and 'convergence' diagnostics
  args <- list(...)
  
  ### if no method specified, return nothing
  if (is.null(args$method)) {
    
    return(NULL)
    
  } else {
    method <- args$method
  }
  
  args$method <- NULL

  # checks 
  if (!is(args$stk0, "FLS")) stop("stk must be of class FLStock")
  if (!is(args$idx0, "FLIndices")) stop("idx must be of class FLIndices")
  
  # dispatch
  out <- do.call(method, args)
  
  if (!is(out$stk0, "FLS")) stop("stk must be of class FLStock")
  
  out
}

### wrapper for calling SAM
### this makes used of the SAM wrapper in FLfse,
### which in turn calls SAM from the package stockassessment
SAM_wrapper <- function(stk0, idx0, ay, forecast = FALSE, 
                        fwd_trgt = "fsq", ### what to target in forecast
                        fwd_yrs = 1, ### number of years to add
                        fwd_yrs_average = -3:0, ### years used for averages
                        fwd_yrs_rec_start = NULL, ### recruitment 
                        fwd_yrs_sel = -3:-1, ### selectivity
                        fwd_yrs_lf_remove = -2:-1,
                        fwd_splitLD = TRUE,
                        parallel = FALSE,
                        ...){
  
  args <- list(...)
  
  ### fit SAM to provided data
  fit <- FLR_SAM(stk = stk0, idx = idx0, conf = args$conf, 
                 DoParallel = parallel)
  
  ### convert into FLStock
  stk0 <- SAM2FLStock(object = fit, stk = stk0) 
  
  ### perform forecast to get SSB ay+1
  if (isTRUE(forecast)) {
    
    fit1 <- fit[[2]]
    
    ### overwrite landing fraction with last year, if requested
    if (!is.null(fwd_yrs_lf_remove)) {
      ### index for years to remove/overwrite
      idx_remove <- nrow(fit1$data$landFrac) + args$fwd_yrs_lf_remove
      ### overwrite
      fit1$data$landFrac[idx_remove, ] <- fit1$data$landFrac[rep(nrow(fit1$data$landFrac), length(idx_remove)), ]
    }
    
    ### check how to do forecast
    ### currently, can only do F status quo
    if (fwd_trgt != "fsq") stop("only fsq supported in forecast")
    
    ### years for average values
    ave.years <- max(fit1$data$years) + fwd_yrs_average
    ### years for sampling of recruitment years
    rec.years <- seq(from = fwd_yrs_rec_start, max(fit1$data$years))
    ### years where selectivity is not used for mean in forecast
    overwriteSelYears <- max(fit1$data$years) + fwd_yrs_sel
    
    ### TODO:
    ### stockassessment::forecast needs to be changes to return stock 
    ### details, so far only summary stats SSB, catch, fbar, ... are returned
    fc <- stockassessment::forecast(fit = fit1, fscale = rep(1, fwd_yrs),
                                    ave.years = ave.years,
                                    rec.years = rec.years,
                                    overwriteSelYears = overwriteSelYears,
                                    splitLD = fwd_splitLD)
    
    ### ISSUES: forecast fails if SAM did not converge, creates error and stops
    
  }
  
  ### save convergence for all iterations
  tracking["convergence", ac(ay)] <- sapply(fit, function(x) x$opt$convergence)
  ### add perceived F
  tracking["Fperc", ac(ay)] <- fbar(stk0)[, ac(ay)]
  
  ### save model fit (list) as attribute in stk0
  attr(stk0, "fit") <- fit
  
  ### return assessed stock, tracking & model fit 
  ### (model fit required later for TAC calculation)
  return(list(stk0 = stk0, tracking = tracking))
  
}