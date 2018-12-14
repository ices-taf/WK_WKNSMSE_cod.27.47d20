### ------------------------------------------------------------------------ ###
### kFun: Management implementation function ####
### ------------------------------------------------------------------------ ###
### Evaluate the chosen Management implementation function
### Evaluate the chosen HCR function using the current stock perception and 
### a control object


kFun <- function(...){
  
  args <- list(...)
  
  ### if no method specified, return nothing
  if (is.null(args$method)) {
    
    return(NULL)
    
  } else {
    
    method <- args$method
    
  }
  
  args$method <- NULL
  
  # Check inputs
  if (!is(args$stk0, "FLS")) stop("stk0 argument must be an FLStock")
  
  # dispatch
  out <- do.call(method, args)
  
  # check outputs
  if (!is(out$ctrl, "fwdControl")) {
    stop("The HCR must return an object of class fwdControl")
  }
  
  # return
  return(out)
  
}

k_SAM_forecast <- function(tracking, stk0, ay, ctrl, it,
                           fwd_trgt = c("fsq", "hcr"), ### target in forecast
                           fwd_yrs = 2, ### number of years to add
                           fwd_yrs_average = -3:0, ### years used for averages
                           fwd_yrs_rec_start = NULL, ### recruitment 
                           fwd_yrs_sel = -3:-1, ### selectivity
                           fwd_yrs_lf_remove = -2:-1,
                           fwd_splitLD = TRUE, 
                           ...) {
  
  ### retrieve SAM model fit (list)
  fit <- attr(stk0, "fit")
  
  ### check class of model fit(s)
  if (!class(fit) %in% c("sam", "sam_list")) 
    stop("attr(stk0, \"fit\") has to be class sam or sam_list")
  
  ### if single fit, turn into list
  if (is(fit, "sam")) fit <- list(fit)
  
  ### go through all model fits
  fc <- foreach(fit_i = fit, iter_i = seq_along(fit), 
                .errorhandling = "pass") %do% {
    
    ### overwrite landing fraction with last year, if requested
    if (!is.null(fwd_yrs_lf_remove)) {
      ### index for years to remove/overwrite
      idx_remove <- nrow(fit_i$data$landFrac) + fwd_yrs_lf_remove
      ### overwrite
      fit_i$data$landFrac[idx_remove, ] <- fit_i$data$landFrac[rep(nrow(fit_i$data$landFrac), length(idx_remove)), ]
    }
    
    ### check how to do forecast
    ### currently, can only do F status quo and F target from ctrl object
    fscale <- ifelse(fwd_trgt == "fsq", 1, NA) ### scaled Fsq
    ### target F values
    fval <- ifelse(fwd_trgt == "hcr", ctrl@trgtArray[, "val", iter_i], NA) 
    
    ### years for average values
    ave.years <- max(fit_i$data$years) + fwd_yrs_average
    ### years for sampling of recruitment years
    if (is.null(fwd_yrs_rec_start)) {
      rec.years <- fit_i$data$years ### use all years, if not defined
    } else {
      rec.years <- seq(from = fwd_yrs_rec_start, max(fit_i$data$years))
    }
    
    ### years where selectivity is not used for mean in forecast
    overwriteSelYears <- max(fit_i$data$years) + fwd_yrs_sel
    
    ### forecast 
    fc <- stockassessment::forecast(fit = fit_i, fscale = fscale, fval = fval,
                                    ave.years = ave.years,
                                    rec.years = rec.years,
                                    overwriteSelYears = overwriteSelYears,
                                    splitLD = fwd_splitLD)
    
    ### return catch in ay + 1
    return(attr(fc, "tab")[ac(ay + 1), "catch:median"])
    
  }
  ### if forecast fails, error message returned
  ### replace error message with NA
  catch_target <- unlist(ifelse(sapply(fc, is, "error"), NA, fc))

  ### create ctrl object
  ctrl <- getCtrl(values = catch_target, quantity = "catch", 
                  years = ay + 1, it = it)
  
  ### save result in tracking
  tracking["Implementation", ac(ay + 1)] <- catch_target
  
  return(list(ctrl = ctrl, tracking = tracking))

}


