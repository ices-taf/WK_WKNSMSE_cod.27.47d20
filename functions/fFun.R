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
                        #newtonsteps = 3, 
                        conf = NULL,
                        tracking = NULL,
                        ...){
  
  args <- list(...)
  
  ### fit SAM to provided data
  fit <- FLR_SAM(stk = stk0, idx = idx0, conf = conf, 
                 DoParallel = parallel, ...)
  
  ### convert into FLStock
  stk0 <- SAM2FLStock(object = fit, stk = stk0) 
  
  ### perform forecast to get SSB ay+1
  if (isTRUE(forecast)) {
    
    ### check how to do forecast
    ### currently, can only do F status quo
    if (fwd_trgt != "fsq") stop("only fsq supported in forecast")
    
    ### years for average values
    ave.years <- range(stk0)[["maxyear"]] + fwd_yrs_average
    ### years for sampling of recruitment years
    rec.years <- seq(from = fwd_yrs_rec_start, to = range(stk0)[["maxyear"]] - 1)
    ### years where selectivity is not used for mean in forecast
    overwriteSelYears <- range(stk0)[["maxyear"]] + fwd_yrs_sel
    
    lst_yr <- range(stk0)[["maxyear"]]
    
    ### extend stk0
    stk0 <- window(stk0, end = range(stk0)[["maxyear"]] + fwd_yrs)
    
    ### modify fwd_yrs in case last data year does not include catch
    if (all(is.na(catch(stk0)[, ac(lst_yr)]))) fwd_yrs <- fwd_yrs + 1
    
    ### forecast years
    yrs <- seq(to = dims(stk0)$maxyear, length.out = fwd_yrs)
    
    ### coerce fit into list if only 1 iteration
    if (is(fit, "sam")) {
      fit <- list(fit)
      class(fit) <- "sam_list"
    }
    
    ### do forecast for all iterations
    fc <- lapply(fit, function(fit_i) {
      
      ### overwrite landing fraction with last year, if requested
      if (!is.null(fwd_yrs_lf_remove)) {
        ### index for years to remove/overwrite
        idx_remove <- nrow(fit_i$data$landFrac) + args$fwd_yrs_lf_remove
        ### overwrite
        fit_i$data$landFrac[idx_remove, ] <- fit_i$data$landFrac[rep(nrow(fit_i$data$landFrac), length(idx_remove)), ]
      }
      
      ### run forecast
      fc_i <- stockassessment::forecast(fit = fit_i, fscale = rep(1, fwd_yrs),
                                        ave.years = ave.years,
                                        rec.years = rec.years,
                                        overwriteSelYears = overwriteSelYears,
                                        splitLD = fwd_splitLD)
      
      ### get numbers at age for all forecast years
      numbers <- lapply(seq(fwd_yrs), function(x) {
        ### index for numbers at age
        idx <- seq(length(fit_i$conf$keyLogFsta[1, ]))
        ### get simulated numbers
        n <- exp(fc_i[[x]]$sim[, idx])
        ### median
        apply(n, 2, median)
      })
      numbers <- do.call(cbind, numbers)
      
      return(list(stock.n = numbers))
      
    })
    
    ### get numbers
    fc_stock.n <- lapply(fc, "[[", "stock.n")
    
    ### insert stock numbers
    stock.n(stk0)[, ac(yrs)] <- unlist(fc_stock.n)
    
    ### extend stock characteristics required for calculation of SSB, 
    ### weights, etc.
    
    ### find years to fill (do not fill years, if there is already data inside)
    yrs_fill <- setdiff(yrs, lst_yr)
    
    stock.wt(stk0)[, ac(yrs_fill)] <- yearMeans(stock.wt(stk0)[, ac(ave.years)])
    m(stk0)[, ac(yrs_fill)] <- yearMeans(m(stk0)[, ac(ave.years)])
    mat(stk0)[, ac(yrs_fill)] <- yearMeans(mat(stk0)[, ac(ave.years)])
    
    harvest.spwn(stk0)[, ac(yrs_fill)] <- yearMeans(harvest.spwn(stk0)[, 
                                                                ac(ave.years)])
    m.spwn(stk0)[, ac(yrs_fill)] <- yearMeans(m.spwn(stk0)[, ac(ave.years)])
    harvest(stk0)[, ac(yrs_fill)] <- yearMeans(harvest(stk0)[, ac(ave.years)])
    
    #ssb(stk0)
    
    ### PLEASE NOTE:
    ### SSB value slightly different from SSB value generated from SAM:
    ### SAM calculates SSB per simulation and gives median
    ### here: calculate median of numbers at age and calculate SSB from
    ###       median numbers
    
    ### ISSUES: forecast fails if SAM did not converge, creates error and stops
    ### 
    
  }
  
  ### save convergence for all iterations
  tracking["convergence", ac(ay)] <- sapply(fit, function(x) x$opt$convergence)
  ### add perceived F
  tracking["MP.f", ac(ay)] <- fbar(stk0)[, ac(ay - 1)]
  ### and SSB
  tracking["MP.ssb", ac(ay)] <- tail(ssb(stk0))
  
  ### save model fit (list) as attribute in stk0
  attr(stk0, "fit") <- fit
  
  ### return assessed stock, tracking & model fit 
  ### (model fit required later for TAC calculation)
  return(list(stk0 = stk0, tracking = tracking))
  
}