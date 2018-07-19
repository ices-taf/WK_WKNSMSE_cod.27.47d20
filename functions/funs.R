### ------------------------------------------------------------------------ ###
### extract uncertainty from SAM object ####
### ------------------------------------------------------------------------ ###
### function for creating iterations based on estimation of uncertainty in SAM



SAM_uncertainty <- function(fit, n = 1000) {
  
  if (!is(fit, "sam")) stop("fit has to be class \"sam\"")
  
  ### index for fishing mortality ages
  idxF <- fit$conf$keyLogFsta[1, ] + dim(stk)[1] + 1
  idxF <- idxF[idxF != 0] ### remove 0s
  
  ### index for F variances (usually some ages are bound)
  idxNVar <- fit$conf$keyVarLogN 
  
  ### get ages used for calculating fbar
  bAges <- fit$conf$fbarRange
  bAges <- do.call(':', as.list(bAges))
  
  ### index for stock numbers ages
  #idxN <- 1:ncol(natural.mortality)
  idxN <- seq(min(fit$data$minAgePerFleet), max(fit$data$maxAgePerFleet))
  ### index for observation variances
  idxObs <- fit$conf$keyVarObs # starts at 0
  
  ##Resample estimated values to get N, F and q 
  
  ### calculate standard deviations of model parameters
  sds <- TMB::sdreport(obj = fit$obj, par.fixed = fit$opt$par, 
                       getJointPrecision = TRUE)
  ### extract values for parameters
  est <- c(sds$par.fixed, sds$par.random)
  ### estimate covariance?
  cov <- solve(sds$jointPrecision)
  
  ### create random values based on estimation and covariance
  sim.states <- mvrnorm(n, est, cov)
  ### matrix, columns are values, rows are requested samples
  table(colnames(sim.states))
  ### contains, among others, logF, logN...
  
  ### combine SAM estimate and random samples
  dat <- rbind(est, sim.states)
  
  ### ---------------------------------------------------------------------- ###
  ### stock ####
  ### ---------------------------------------------------------------------- ###
  
  ### stock characteristics
  min_age <- min(fit$data$minAgePerFleet[fit$data$fleetTypes == 0])
  max_age <- max(fit$data$maxAgePerFleet[fit$data$fleetTypes == 0])
  years <- fit$data$years
  ### FLQuant template for stock
  stk_template <- FLQuant(dimnames = list(age = min_age:max_age, year = years,
                                          iter = 1:(n + 1)))
  
  ### numbers at age
  stock.n <- stk_template
  stock.n[] <- exp(t(dat[, colnames(dat) == "logN"]))
  
  ### F at age
  harvest <- stk_template
  harvest[] <- exp(t(dat[, colnames(dat) == "logF"]))
  
  ### ---------------------------------------------------------------------- ###
  ### surveys ####
  ### ---------------------------------------------------------------------- ###
  
  ### survey specs
  idx_surveys <- which(fit$data$fleetTypes > 0) ### which observation are surveys
  ### age range of surveys
  survey_ages <- lapply(seq_along(idx_surveys), function(x) {
    seq(fit$data$minAgePerFleet[idx_surveys][x],
        fit$data$maxAgePerFleet[idx_surveys][x])
  })
  ### index for estimated parameters
  sum(colnames(dat) == "logFpar") ### there are 9 parameters
  ### I assume: 5 for Q1, 4 for Q3, as they have this many ages...
  survey_ages_idx <- split(seq(length(unlist(survey_ages))), 
                           rep(seq(survey_ages), sapply(survey_ages, length))) 
  
  ### get catchability at age (time-invariant) samples
  catchability <- lapply(seq_along(idx_surveys), function(x) {
    
    ### create FLQuant template
    tmp <- FLQuant(dimnames = list(age = survey_ages[[x]],
                                   year = "all", iter = 1:(n + 1)))
    ### fill with catchability values
    tmp[] <- exp(t(dat[, colnames(dat) == "logFpar"][,survey_ages_idx[[x]]]))
    
    return(tmp)
    
  })
  
  ### ---------------------------------------------------------------------- ###
  ### standard deviation - catch ####
  ### ---------------------------------------------------------------------- ###
  ### time-invariant
  
  ### template
  catch_sd <- FLQuant(dimnames = list(age = dimnames(stock.n)$age, year = "all",
                                      iter = 1:(n + 1)))
  
  ### index for catch sd (some ages are linked)
  catch_sd_idx <- idxObs[1, ][idxObs[1, ] > -1] + 1
  
  ### extract values
  catch_sd[] <- exp(t(dat[, colnames(dat) == "logSdLogObs"][, catch_sd_idx]))
  
  ### ---------------------------------------------------------------------- ###
  ### standard deviation - surveys ####
  ### ---------------------------------------------------------------------- ###
  ### time-invariant
  
  ### get catchability at age (time-invariant) samples
  survey_sd <- lapply(seq_along(idx_surveys), function(x) {
    
    ### create FLQuant template
    tmp <- FLQuant(dimnames = list(age = survey_ages[[x]],
                                   year = "all", iter = 1:(n + 1)))
    ### index for sd (some ages are linked)
    idx <- idxObs[idx_surveys[x], ]
    idx <- idx[idx > -1] + 1
    
    ### fill with catchability values
    tmp[] <- exp(t(dat[, colnames(dat) == "logSdLogObs"][, idx]))
    
    return(tmp)
    
  })
  
  return(list(stock.n = stock.n, harvest = harvest,
              survey_catchability = catchability, catch_sd = catch_sd,
              survey_sd = survey_sd))
  
}

### ------------------------------------------------------------------------ ###
### calculate survey index/indices from stock ####
### ------------------------------------------------------------------------ ###
#' calculate survey index/indices from FLStock
#'
#' This function calculates survey indices from the numbers at age of an 
#' FLStock object
#'
#' @param stk Object of class \linkS4class{FLStock} with stock and fishery data.
#' @param idx Object of class \linkS4class{FLIndices} or \linkS4class{FLIndex}.
#' @param use_q If \code{TRUE} index numbers are calculated by multiplying
#'   numbers at age in the stock with the catchability.
#' @param use_time If \code{TRUE} observed numbers in the survey are corrected
#'   for fishing and natural mortality.
#'
#' @return An object of class \code{FLIndex} or \code{FLIndex} with the survey
#'   index stored in the \code{index} slot.
#'
#' @export

setGeneric("calc_survey", function(stk, idx, use_q = TRUE, use_time = TRUE) {
  standardGeneric("calc_survey")
})

### stk = FLStock, idx = FLIndices
#' @rdname calc_survey
setMethod(f = "calc_survey",
          signature = signature(stk = "FLStock", idx = "FLIndices"),
          definition = function(stk, idx, use_q = TRUE, use_time = TRUE) {
  
  ### apply function to every element of FLIndices
  lapply(X = idx, FUN = calc_survey_ind, stk = stk, use_q = TRUE, use_time = TRUE)
            
})
### stk = FLStock, idx = FLIndex
#' @rdname calc_survey
setMethod(f = "calc_survey",
          signature = signature(stk = "FLStock", idx = "FLIndex"),
          definition = function(stk, idx, use_q = TRUE, use_time = TRUE) {
            
  calc_survey(stk = stk, idx = idx, use_q = TRUE, use_time = TRUE)
            
})



calc_survey_ind <- function(stk, idx, 
                            use_q = TRUE, ### catchability
                            use_time = TRUE ### timing of survey
) {
  
  ### find ranges for years, ages & iters
  ages <- intersect(dimnames(index(idx))$age, dimnames(stock.n(stk))$age)
  years <- intersect(dimnames(index(idx))$year, dimnames(stock.n(stk))$year)
  iter <- intersect(dimnames(index(idx))$iter, dimnames(stock.n(stk))$iter)
  
  ### timing of survey
  if (isTRUE(use_time)) {
    ### use mean of fishing period
    time <- mean(range(idx)[c("startf", "endf")])
  } else {
    ### otherwise assume beginning of year
    time <- 0
  }
  
  ### extract stock numbers for requested/available dimensions
  index.n <- stock.n(stk)[ac(ages), ac(years),,,, ac(iter)]
  ### get Z = M & F
  Z <- m(stk)[ac(ages), ac(years),,,, ac(iter)] +
    harvest(stk)[ac(ages), ac(years),,,, ac(iter)]
  
  ### estimate stock numbers at time of survey
  index.n <- index.n * exp(-time * Z)
  
  ### add catchability, if requested
  if (isTRUE(use_q)) {
    index.n <- index.n * index.q(idx)[ac(ages), ac(years),,,, ac(iter)]
  }
  
  ### insert values into index
  index(idx)[ac(ages), ac(years),,,, ac(iter)] <- index.n
  
  return(idx)
  
}
