### ------------------------------------------------------------------------ ###
### oFun: observations ####
### ------------------------------------------------------------------------ ###

oFun <- function(...){
	# in: stk = FLStock, idx = FLIndices
	#	in: method = character for the wrapper function on the assessment.
	# out: list with elements stk and idx sampled from OM
	args <- list(...)
	
	### select method or use default
	if (is.null(args$method)) {
	  method <- "perfectInfo.wrapper"
	} else {
	  method <- args$method
	}

	args$method <- NULL
	
	### checks 
	if (!is(args$stk, "FLS")) stop("stk must be of class FLStock")
	
	### dispatch
	out <- do.call(method, args)
	
	if (!is(out$stk, "FLS")) stop("stk must be of class FLStock")
	if (!is(out$idx, "FLIndices")) stop("idx must be of class FLIndices")
	
	out
}

### ------------------------------------------------------------------------ ###
### function used for WKNSMSE ####
### ------------------------------------------------------------------------ ###

o_WKNSMSE <- function(stk, idx, ay, tracking, idx_timing = -1,
                      catch_timing = -1, ...) {
  
  ### assume catch is known without error for now
  stk0 <- stk
  
  ### cut of years
  ### workaround for NS cod: survey until intermediate year, but catch stops
  ### 1 year earlier
  ### slots such as natural mortality, maturity need to be kept, otherwise
  ### SAM will fall over
  if (any(idx_timing > catch_timing)) {
    ### keep stock until last survey data year
    stk0 <- window(stk0, end = ay + max(idx_timing))
    ### find years to remove
    yrs_remove <- (ay + catch_timing + 1):ay
    ### remove catch data
    catch(stk0)[, ac(yrs_remove)] <- NA
    catch.n(stk0)[, ac(yrs_remove)] <- NA
    catch.wt(stk0)[, ac(yrs_remove)] <- NA
    landings(stk0)[, ac(yrs_remove)] <- NA
    landings.n(stk0)[, ac(yrs_remove)] <- NA
    landings.wt(stk0)[, ac(yrs_remove)] <- NA
    discards(stk0)[, ac(yrs_remove)] <- NA
    discards.n(stk0)[, ac(yrs_remove)] <- NA
    discards.wt(stk0)[, ac(yrs_remove)] <- NA
  } else {
    stk0 <- window(stk0, end = ay + catch_timing)
  }
  
  ### calculate index values
  idx <- calc_survey(stk = stk, idx = idx)
  idx0 <- idx ### observed survey
  
  ### timing of survey
  ### 0: intermediate/assessment year
  ### <0: fewer years available & vice versa
  if (length(idx_timing) < length(idx)) { 
    idx_timing <- rep(idx_timing, length(idx))
  }
  ### restrict years for indices based on timing
  idx0 <- lapply(seq_along(idx), function(x) {
    window(idx[[x]], end = ay + idx_timing[x])
  })
  idx0 <- FLIndices(idx0) ### restore class
  names(idx0) <- names(idx)
  
  list(stk0 = stk0, idx0 = idx0, idx = idx, 
       tracking = tracking)
  
}





perfectInfo.wrapper <- function(stk, observations, vy0, ay, tracking){
	dataYears <- vy0
	assessmentYear <- ac(ay)
	stk0 <- stk[,dataYears] # Only data years
	catch.n(stk0) <- (catch.n(stk0) + 1)
	idx <- observations$idx
	idx0 <- lapply(idx, function(x) x[,dataYears])
	# Generate full index up to projection year 
	# using catchability from first year, which means no error
	for (idx_count in 1:length(idx)) {
		index(idx[[idx_count]])[,assessmentYear] <- stock.n(stk)[,assessmentYear]*index.q(idx[[idx_count]])[,1]
	}
	list(stk = stk0, idx = idx0, observations = list(idx = idx), 
	     tracking = tracking)
}

