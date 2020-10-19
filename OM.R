### ------------------------------------------------------------------------ ###
### create FLStock for cod ####
### ------------------------------------------------------------------------ ###
### base on SAM assessment

### check versions of required R packages
if (packageVersion("FLCore") < "2.6.11.9001") 
  stop("please update FLCore")
if (packageVersion("FLfse") < "0.0.0.9003") 
  stop("please update FLfse")
if (packageVersion("stockassessment") < "0.8.1") 
  stop("please update stockassessment")
if (packageVersion("mse") < "0.9.1") 
  stop("please update stockassessment")

### load packages
library(FLfse)
library(stockassessment)
library(ggplotFL)
library(FLAssess)
library(mse)
### load files from package mse for easier debugging
# devtools::load_all("../mse/")
library(FLash)
library(tidyr)
library(dplyr)
library(doParallel)

source("a4a_mse_WKNSMSE_funs.R")

dir.create(path = "input/cod4", recursive = TRUE)
dir.create(path = "output/runs/cod4", recursive = TRUE)

### create plots and print to screen?
verbose <- TRUE


### the original WKNSMSE was run with R 3.5
### for exact reproducibility in R 3.6, the random number generation must be
### be changed
if (getRversion() >= 3.6) RNGkind(sample.kind = "Rounding")

### ------------------------------------------------------------------------ ###
### simulation specifications ####
### ------------------------------------------------------------------------ ###

### number of iterations/replicates
n <- 1000
### number of years
n_years <- 20
### last data year
yr_data <- 2018

### ------------------------------------------------------------------------ ###
### fit SAM ####
### ------------------------------------------------------------------------ ###
### use input data provided in FLfse
### recreates the WGNSSK2018 cod assessment
fit <- FLR_SAM(stk = cod4_stk, idx = cod4_idx, conf = cod4_conf_sam)

if (isTRUE(verbose)) {
  is(fit)
  fit
  plot(fit)
}


### ------------------------------------------------------------------------ ###
### remove catch multiplier for cod ####
### ------------------------------------------------------------------------ ###
### For the cod SAM assessment a catch multiplier is estimated for the years
### 1993-2005.
### For the simulation, we correct the catch with this multiplier and remove
### estimation of the catch multiplier 
### -> saves time in simulation (less parameters to estimate)

### get catch multiplier
ages <- fit$conf$minAge:fit$conf$maxAge
yrs <- fit$conf$keyScaledYears
catch_mult <- FLQuant(
  matrix(data = fit$pl$logScale[(fit$conf$keyParScaledYA + 1)], 
         ncol = fit$conf$noScaledYears,
         nrow = length(fit$conf$minAge:fit$conf$maxAge),
         byrow = TRUE), 
  dimnames = list(year = fit$conf$keyScaledYears, 
                  age = fit$conf$minAge:fit$conf$maxAge))
catch_mult <- exp(catch_mult)

cod4_stk2 <- cod4_stk
### correct catch.n
catch.n(cod4_stk2)[ac(ages), ac(yrs)] <- catch.n(cod4_stk2)[ac(ages), ac(yrs)] *
  catch_mult
### split into landings and discards, based on landing fraction
landings.n(cod4_stk2)[ac(ages), ac(yrs)] <- catch.n(cod4_stk2)[ac(ages), ac(yrs)] *
  (landings.n(cod4_stk)[ac(ages), ac(yrs)] / catch.n(cod4_stk)[ac(ages), ac(yrs)])
discards.n(cod4_stk2)[ac(ages), ac(yrs)] <- catch.n(cod4_stk2)[ac(ages), ac(yrs)] *
  (1 - landings.n(cod4_stk)[ac(ages), ac(yrs)] / 
     catch.n(cod4_stk)[ac(ages), ac(yrs)])
### update stock
catch(cod4_stk2)[, ac(yrs)] <- computeCatch(cod4_stk2)[, ac(yrs)]
landings(cod4_stk2)[, ac(yrs)] <- computeLandings(cod4_stk2)[, ac(yrs)]
discards(cod4_stk2)[, ac(yrs)] <- computeDiscards(cod4_stk2)[, ac(yrs)]

### fit SAM to "corrected" catches
cod4_conf_sam_no_mult <- cod4_conf_sam[!names(cod4_conf_sam) %in% 
                                         c("noScaledYears", "keyScaledYears",
                                           "keyParScaledYA")]
fit2 <- FLR_SAM(stk = cod4_stk2, idx = cod4_idx, 
                conf = cod4_conf_sam_no_mult)
### compare results
if (isTRUE(verbose)) summary(fit2) / summary(fit)
### estimates and log likelihood identical, only bounds smaller

### fit SAM as it is done during the MSE simulation
fit_est <- FLR_SAM(stk = cod4_stk2, idx = cod4_idx, 
                   conf = cod4_conf_sam_no_mult, 
                   newtonsteps = 0, rel.tol = 0.001)
### extract model parameters and use them in the simulation as starting values
sam_initial <- sam_getpar(fit_est)
### ------------------------------------------------------------------------ ###
### create FLStock ####
### ------------------------------------------------------------------------ ###
### create template with 1 iteration
### cod4_stk2 is used as template, i.e. the input values (catch) include
### the catch multiplier, 
### the results (stock numbers & harvest) are used from the real WGNSSK fit
stk <- SAM2FLStock(object = fit, stk = cod4_stk2)
if (isTRUE(verbose)) summary(stk)

### set units
units(stk)[1:17] <- as.list(c(rep(c("t", "1000", "kg"), 4),
                              "", "", "f", "", ""))
if (isTRUE(verbose)) plot(stk)

### save for later comparison
stk_orig <- stk

### ------------------------------------------------------------------------ ###
### add uncertainty ####
### ------------------------------------------------------------------------ ###
### first approach: use variance-covariance


### add iteration dimension
stk <- FLCore::propagate(stk, n)
dim(stk)

### add uncertainty estimated by SAM as iterations
set.seed(1)
uncertainty <- SAM_uncertainty(fit = fit, n = n, print_screen = FALSE, 
                               idx_cov = TRUE, catch_est = TRUE)
### add noise to stock
stock.n(stk)[] <- uncertainty$stock.n
stock(stk)[] <- computeStock(stk)
### add noise to F
harvest(stk)[] <- uncertainty$harvest

### add noise to catch numbers
catch.n(stk)[, dimnames(stk)$year[-dims(stk)$year]] <- uncertainty$catch_n
catch(stk) <- computeCatch(stk)

if (isTRUE(verbose)) plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

### maximum observed F
max(fbar(stk))
# 1.359563 in year 1999 with 10,000 iterations
max(harvest(stk))
# 2.086832 in year 2001 for age 6 (plusgroup)

### get estimated catch numbers
catch_n <- uncertainty$catch_n
catch_n[dimnames(catch_mult)$age, dimnames(catch_mult)$year] <- 
  catch_n[dimnames(catch_mult)$age, dimnames(catch_mult)$year] * catch_mult


### ------------------------------------------------------------------------ ###
### check MCMC approach ####
### ------------------------------------------------------------------------ ###
# library(tmbstan) # cran package
# 
# ### create template stock for storing results
# MCMC_iter <- n
# MCMC_warmup <- 1000
# stk_MCMC <- propagate(stk, MCMC_iter)
# 
# ### run MCMC
# system.time(mcmc <- tmbstan(fit$obj, chains = 1, iter = MCMC_warmup + MCMC_iter, 
#                             warmup = MCMC_warmup, 
#                             seed = 1, control = list(max_treedepth = 15)))
# ### extract 
# mc <- extract(mcmc, inc_warmup = FALSE, permuted = FALSE)
# 
# table(gsub(x = dimnames(mc)$parameters, pattern = "\\[[0-9]{1,}\\]", 
#            replacement = ""))
# # itrans_rho         logF      logFpar         logN     logScale logSdLogFsta 
# #          1          336            9          336           13            2 
# # logSdLogN  logSdLogObs         lp__ 
# #         2            7            1 
# ### gives N and F @age
# ### scale for 
# 
# 
# ### find positions for results in MCMC
# nms <- dimnames(mc)$parameters
# F_pos <- grep(x = nms, pattern = "logF\\[[0-9]{1,3}\\]$")
# N_pos <- grep(x = nms, pattern = "logN\\[[0-9]{1,3}\\]$")
# factor_pos <- grep(x = nms, pattern = "logScale\\[[0-9]{1,2}\\]$")
# 
# ### in the MCMC results age and years are mixed within the same row,
# ### each row represents one iteration
# ### this needs to be reformatted to be useful...
# 
# ### fishing mortality:
# harvest(stk_MCMC)[] <- aperm(array(data = c(exp(mc[,, F_pos])),
#                                    dim = dim(harvest(stk_MCMC))[c(6, 1:5)]),
#                              perm = c(2:6, 1))
# 
# ### stock numbers at age
# stock.n(stk_MCMC)[] <- aperm(array(data = c(exp(mc[,, N_pos])),
#                                    dim = dim(stock.n(stk_MCMC))[c(6, 1:5)]),
#                              perm = c(2:6, 1))
# stock(stk_MCMC) <- computeStock(stk_MCMC)
# 
# ### catch factor
# ### get ages and years
# ages <- fit$conf$minAge:fit$conf$maxAge
# yrs <- fit$conf$keyScaledYears
# ### get catch factors
# catch_factor <- lapply(split(exp(mc[,, factor_pos]), seq(MCMC_iter)), 
#                        function(x) {
#   x[t(fit$conf$keyParScaledYA + 1)]
# })
# catch_factor <- unlist(catch_factor)
# ### coerce into FLQuant
# catch_factor <- FLQuant(catch_factor,
#                       dimnames = dimnames(catch.n(stk_MCMC[ac(ages), ac(yrs)])))
# ### multiply catch numbers
# catch.n(stk_MCMC)[ac(ages), ac(yrs)] <- catch_factor * 
#   catch.n(stk_MCMC)[ac(ages), ac(yrs)]
# ### split into landings and discards, based on landing fraction
# lfrac <- propagate((landings.n(stk)[ac(ages), ac(yrs)] / 
#                       catch.n(stk)[ac(ages), ac(yrs)]), MCMC_iter)
# landings.n(stk_MCMC)[ac(ages), ac(yrs)] <- catch.n(stk_MCMC)[ac(ages), ac(yrs)] *
#   lfrac
# discards.n(stk_MCMC)[ac(ages), ac(yrs)] <- catch.n(stk_MCMC)[ac(ages), ac(yrs)] *
#   (1 - lfrac)
# ### update stock
# catch(stk_MCMC)[, ac(yrs)] <- computeCatch(stk_MCMC)[, ac(yrs)]
# landings(stk_MCMC)[, ac(yrs)] <- computeLandings(stk_MCMC)[, ac(yrs)]
# discards(stk_MCMC)[, ac(yrs)] <- computeDiscards(stk_MCMC)[, ac(yrs)]
# 
# ### plot MCMC stock
# plot(stk_MCMC)
# ### compare with original SAM fit
# plot(FLStocks(MCMC = stk_MCMC, original = stk_orig))
# ### compare with first uncertainty approach
# plot(FLStocks(MCMC = stk_MCMC, original = stk_orig, Cov = stk))

### ------------------------------------------------------------------------ ###
### try SAM internal "simstudy" ####
### ------------------------------------------------------------------------ ###
### "Simulate data from fitted model and re-estimate from each run"

# set.seed(0)
# system.time(fits <- simstudy(fit = fit, nsim = n))
# class(fits) <- "sam_list"
# 
# stk_sim <- SAM2FLStock(object = fits, stk = cod4_stk2)
# plot(FLStocks(Cov = stk, simstudy = stk_sim,
#               original = stk_orig))

### ------------------------------------------------------------------------ ###
### extend stock for MSE simulation ####
### ------------------------------------------------------------------------ ###

### special case for NS cod
### maturity data available for 2018 (based on the IBTS Q1)
### stock weights, M etc in 2018 based on a three year average to enable 
### calculation of SSB
### although SAM estimates F in 2018, this is not reported or taken forward into
### forcasts by the WG
# stk_stf2017 <- stf(window(stk, end = 2017), n_years + 1)
stk_stf2018 <- stf(window(stk, end = 2018), n_years)
### use all available data
stk_stf <- stk_stf2018

### ------------------------------------------------------------------------ ###
### biological data for OM ####
### ------------------------------------------------------------------------ ###
### Resample weights, maturity and natural mortality from the last 5 years 
### (2013-2017)
### set up an array with one resampled year for each projection year 
### (including intermediate year) and replicate
### use the same resampled year for all biological parameters
### this is the approach used in eqsim for North Sea cod
set.seed(2)

### use last five data years to sample biological parameters
sample_yrs <- 2013:2017
### get year position of sample years
sample_yrs_pos <- which(dimnames(stk_stf)$year %in% sample_yrs)

### create samples for biological data (weights, etc.)
### the historical biological parameters are identical for all iterations
### and consequently do not need to be treated individually
### (but keep age structure)
### create vector with resampled years
bio_samples <- sample(x = sample_yrs_pos, 
                      size = (n_years + 1) * n, replace = TRUE)
### do the same for selectivity
sel_samples <- sample(x = sample_yrs_pos, 
                      size = (n_years + 1) * n, replace = TRUE)
### years to be populated
bio_yrs <- which(dimnames(stk_stf)$year %in% 2018:dims(stk_stf)$maxyear)


### insert values
catch.wt(stk_stf)[, bio_yrs] <- c(catch.wt(stk)[, bio_samples,,,, 1])
stock.wt(stk_stf)[, bio_yrs] <- c(stock.wt(stk)[, bio_samples,,,, 1])
landings.wt(stk_stf)[, bio_yrs] <- c(landings.wt(stk)[, bio_samples,,,, 1])
discards.wt(stk_stf)[, bio_yrs] <- c(discards.wt(stk)[, bio_samples,,,, 1])
m(stk_stf)[, bio_yrs] <- c(m(stk)[, bio_samples,,,, 1])
mat(stk_stf)[, bio_yrs] <- c(mat(stk)[, bio_samples,,,, 1])
### maturity data for 2018 exists, re-insert real data
mat(stk_stf)[, ac(2018)] <- mat(stk_orig)[, ac(2018)]
### use different samples for selectivity
harvest(stk_stf)[, bio_yrs] <- c(harvest(stk)[, sel_samples,,,, 1])

if (isTRUE(verbose)) plot(stk_stf)

### ------------------------------------------------------------------------ ###
### stock recruitment ####
### ------------------------------------------------------------------------ ###
### fit hockey-stick model
### get residuals from smoothed residuals

### use only data from 1997 and later
sr <- as.FLSR(window(stk_stf, start = 1997), model = "segreg")
### fit model individually to each iteration and suppress output to screen
suppressWarnings(. <- capture.output(sr <- fmle(sr)))
### run in parallel
# library(doParallel)
# cl <- makeCluster(10)
# registerDoParallel(cl)
# sr <- fmle_parallel(sr, cl)
# ### run again for failed iterations
# pos_error <- which(is.na(params(sr)["a"]))
# sr_corrected <- fmle(FLCore::iter(sr, pos_error))
# sr[,,,,, pos_error] <- sr_corrected[]
# params(sr)[, pos_error] <- params(sr_corrected)

if (isTRUE(verbose)) {
  
  plot(sr)
  ### check breakpoints
  summary(params(sr)["b"])
  
  ### plot model and data
  as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb)) %>%
    mutate(age = NULL,
           year = ifelse(qname == "SSB", year + 1, year)) %>%
    tidyr::spread(key = qname, value = data) %>%
    ggplot() +
    geom_point(aes(x = SSB, y = rec, group = iter), 
               alpha = 0.5, colour = "grey", shape = 1) +
    geom_line(aes(x = SSB, y = fitted, group = iter)) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  
  ### Check extent of autocorrelation
  # Not significant, so no need to account for it in this OM
  acf(window(stock.n(stk_orig)[1], start = 1998))
  
  ### Check method proposed for generating recruitment compares with past recruitment estimates
  test <- as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb))
  test <- mutate(test, age = NULL, year = ifelse(qname == "SSB", year + 1, year))
  test <- tidyr::spread(test, key = qname, value = data)
  test <- test[complete.cases(test),]
  test$res <- rep(NA, nrow(test))
  
  # Generate residuals for future recruitments
  foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
          .errorhandling = "pass") %do% {
            
            set.seed(iter_i^2)
            
            ### get residuals for current iteration
            res_i <- c(FLCore::iter(residuals(sr), iter_i))
            res_i <- res_i[!is.na(res_i)]
            
            ### calculate kernel density of residuals
            density <- density(x = res_i)
            ### sample residuals
            mu <- sample(x = res_i, size = length(res_i), replace = TRUE)
            ### "smooth", i.e. sample from density distribution
            test$res[test$iter==iter_i] <- rnorm(n = length(res_i), mean = mu, sd = density$bw)
            
          }
  
  # Generate future recruitments from past SSBs and generated residuals
  test$future <- test$fitted * exp(test$res)
  
  # 10 randomly selected iters for plotting
  # should probably increase the number later
  i_samp <- sample(seq(dim(sr)[6]), 10, replace=FALSE)
  
  # Plot past and future stock recruit pairs for selected iters
  ggplot(test[is.element(test$iter, i_samp),]) +
    geom_point(aes(x = SSB, y = rec), 
               alpha = 0.5, colour = "red", shape = 19) +
    geom_point(aes(x = SSB, y = future), 
               alpha = 0.5, colour = "black", shape = 19) +
    geom_line(aes(x = SSB, y = fitted)) +
    facet_wrap(~iter) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  
  # Empirical cumulative distributions for the same iters
  ggplot(test[is.element(test$iter, i_samp),]) +
    stat_ecdf(aes(rec), geom = "step", colour = "red") +
    stat_ecdf(aes(future), geom = "step", colour = "black") +
    facet_wrap(~iter) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  
  # Combine previous two plots over iters
  # Stock recruit pairs
  ggplot(test[is.element(test$iter, i_samp),]) +
    geom_point(aes(x = SSB, y = rec), 
               alpha = 0.5, colour = "red", shape = 19) +
    geom_point(aes(x = SSB, y = future), 
               alpha = 0.5, colour = "black", shape = 19) +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  
  # Empirical cumulative distribution
  ggplot(test[is.element(test$iter, i_samp),]) +
    stat_ecdf(aes(rec), geom = "step", colour = "red") +
    stat_ecdf(aes(future), geom = "step", colour = "black") +
    theme_bw() + xlim(0, NA) + ylim(0, NA)
  
  rm(test, i_samp)
}

### generate residuals for MSE
### years with missing residuals
# NW: dimnames produces NULL for me
# yrs_res <- dimnames(sr)$year[which(is.na(iterMeans(rec(sr))))]
yrs_res <- colnames(rec(sr))[which(is.na(iterMeans(rec(sr))))]

### go through iterations and create residuals
### use kernel density to create smooth distribution of residuals
### and sample from this distribution
res_new <- foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
                   .errorhandling = "pass") %dopar% {
                     
  set.seed(iter_i)
  
  ### get residuals for current iteration
  res_i <- c(FLCore::iter(residuals(sr), iter_i))
  res_i <- res_i[!is.na(res_i)]
  
  ### calculate kernel density of residuals
  density <- density(x = res_i)
  ### sample residuals
  mu <- sample(x = res_i, size = length(yrs_res), replace = TRUE)
  ### "smooth", i.e. sample from density distribution
  res_new <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)

  return(res_new)
  
}
summary(exp(unlist(res_new)))
### insert into model
residuals(sr)[, yrs_res] <- unlist(res_new)
### exponeniate residuals to get factor
residuals(sr) <- exp(residuals(sr))
sr_res <- residuals(sr)

if (isTRUE(verbose)) plot(sr_res)

### ------------------------------------------------------------------------ ###
### process noise ####
### ------------------------------------------------------------------------ ###
### create FLQuant with process noise
### this will be added to the values obtained from fwd() in the MSE

### create noise for process error
set.seed(3)
proc_res <- stock.n(stk_stf) %=% 0 ### template FLQuant
proc_res[] <- stats::rnorm(n = length(proc_res), mean = 0, 
                           sd = uncertainty$proc_error)
### the proc_res values are on a normale scale,
### exponentiate to get log-normal 
proc_res <- exp(proc_res)
### proc_res is a factor by which the numbers at age are multiplied

### for historical period, numbers already include process error from SAM
### -> remove deviation
proc_res[, dimnames(proc_res)$year <= 2017] <- 1

### remove deviation for first age class (recruits)
proc_res[1, ] <- 1

### try saving in stock recruitment model ... 
### this gets passed on to the projection module
fitted(sr) <- proc_res

if (isTRUE(verbose)) plot(proc_res)


### ------------------------------------------------------------------------ ###
### stf for 2018: assume catch advice is taken ####
### ------------------------------------------------------------------------ ###
c2018 <- 53058
ctrl <- fwdControl(data.frame(year = 2018, quantity = "catch", 
                              val = c2018))

### project forward for intermediate year (2018)
stk_int <- stk_stf
stk_int[] <- fwd(stk_stf, ctrl = ctrl, sr = sr, sr.residuals = sr_res,
                 sr.residuals.mult = TRUE, maxF = 5)[]
### add process noise
stock.n(stk_int) <- stock.n(stk_int) * proc_res
stock(stk_int)[] <- computeStock(stk_int)

### create stock for MSE simulation
stk_fwd <- stk_stf
### insert values for 2018
stk_fwd[, ac(2018)] <- stk_int[, ac(2018)]
### insert stock number for 2019 in order to calculate SSB at beginning of 
### 2019
stock.n(stk_fwd)[, ac(2019)] <- stock.n(stk_int)[, ac(2019)]
stock(stk_fwd)[, ac(2019)] <- computeStock(stk_fwd[, ac(2019)])

#all.equal(window(stk_fwd, end = 2018), window(stk_stf, end = 2018))

### ------------------------------------------------------------------------ ###
### biological data for OEM ####
### ------------------------------------------------------------------------ ###

### base on OM stock
stk_oem <- stk_fwd

### projection years
proj_yrs <- 2018:range(stk_oem)[["maxyear"]]

### use means of sampled values for projection period
catch.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(catch.wt(stk_oem)[, ac(sample_yrs)])
landings.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(landings.wt(stk_oem)[, ac(sample_yrs)])
discards.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(discards.wt(stk_oem)[, ac(sample_yrs)])
stock.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(stock.wt(stk_oem)[, ac(sample_yrs)])
m(stk_oem)[, ac(proj_yrs)] <- yearMeans(m(stk_oem)[, ac(sample_yrs)])
### maturity starts one year later because there is data for 2018
mat(stk_oem)[, ac(proj_yrs[-1])] <- yearMeans(mat(stk_oem)[, ac(sample_yrs)])

### remove stock assessment results
stock.n(stk_oem)[] <- stock(stk_oem)[] <- harvest(stk_oem)[] <- NA


### ------------------------------------------------------------------------ ###
### indices ####
### ------------------------------------------------------------------------ ###
### use real FLIndices object as template (included in FLfse)
idx <- cod4_idx
### extend for simulation period
idx <- window(idx, end = yr_data + n_years)
### add iterations
idx <- lapply(idx, propagate, n)

### insert catchability
for (idx_i in seq_along(idx)) {
  
  ### set catchability for projection
  index.q(idx[[idx_i]])[] <- uncertainty$survey_catchability[[idx_i]]
  
}
### create copy of index with original values
idx_raw <- lapply(idx ,index)
### calculate index values
idx <- calc_survey(stk = stk_fwd, idx = idx)

### create deviances for indices
### first, get template
idx_dev <- lapply(idx, index)
### create random noise based on sd
set.seed(4)
for (idx_i in seq_along(idx_dev)) {
  ### insert sd
  idx_dev[[idx_i]][] <- uncertainty$survey_sd[[idx_i]]
  ### noise
  idx_dev[[idx_i]][] <- stats::rnorm(n = length(idx_dev[[idx_i]]),
                                     mean = 0, sd = idx_dev[[idx_i]])
  ### exponentiate to get from normal to log-normal scale
  idx_dev[[idx_i]] <- exp(idx_dev[[idx_i]])
}

### modify residuals for historical period so that index values passed to 
### stock assessment are the ones observed in reality
### IBTS Q1, values up to 2018
idx_dev$IBTS_Q1_gam[, dimnames(idx_dev$IBTS_Q1_gam)$year <= 2018] <- 
  idx_raw$IBTS_Q1_gam[, dimnames(idx_raw$IBTS_Q1_gam)$year <= 2018] /
  index(idx$IBTS_Q1_gam)[, dimnames(idx$IBTS_Q1_gam@index)$year <= 2018]
### IBTS Q3, values up to 2017
idx_dev$IBTS_Q3_gam[, dimnames(idx_dev$IBTS_Q3_gam)$year <= 2017] <- 
  idx_raw$IBTS_Q3_gam[, dimnames(idx_raw$IBTS_Q3_gam)$year <= 2017] /
  index(idx$IBTS_Q3_gam)[, dimnames(idx$IBTS_Q3_gam@index)$year <= 2017]

if (isTRUE(verbose)) {

  ### compare simulated to original survey(s)
  as.data.frame(FLQuants(cod4_q1 = index(cod4_idx$IBTS_Q1_gam), 
                         cod4_q3 = index(cod4_idx$IBTS_Q3_gam),
                         sim_q1 = (index(idx$IBTS_Q1_gam)),
                         sim_q3 = (index(idx$IBTS_Q3_gam))
  )) %>%
    mutate(survey = ifelse(grepl(x = qname, pattern = "*_q1$"), "Q1", "Q3"),
           source = ifelse(grepl(x = qname, pattern = "^sim*"), "sim", "data")) %>%
    filter(year <= 2019) %>%
    ggplot(aes(x = year, y = data, colour = source)) +
    facet_grid(paste("age", age) ~ paste("IBTS", survey), scales = "free_y") +
    stat_summary(fun.y = quantile, fun.args = 0.25, geom = "line",
                 alpha = 0.5) +
    stat_summary(fun.y = quantile, fun.args = 0.75, geom = "line",
                 alpha = 0.5) +
    stat_summary(fun.y = median, geom = "line") +
    theme_bw()

}

### check survey
# idx0 <- calc_survey(stk = stk_fwd, idx = idx)
# idx0 <- lapply(seq_along(idx0), function(idx_i) {
#   idx_tmp <- idx0[[idx_i]]
#   index(idx_tmp) <- index(idx_tmp) * idx_dev[[idx_i]]
#   return(idx_tmp)
# })
# plot(index(idx0[[2]]))

### ------------------------------------------------------------------------ ###
### catch noise ####
### ------------------------------------------------------------------------ ###
### take estimates from sam: uncertainty$catch_sd is "logSdLogObs"
### assume catch observed by SAM in projection is log-normally distributed
### around operating model catch

### create noise for catch
set.seed(5)
catch_res <- catch.n(stk_fwd) %=% 0 ### template FLQuant
catch_res[] <- stats::rnorm(n = length(catch_res), mean = 0, 
                            sd = uncertainty$catch_sd)
### the catch_res values are on a normale scale,
### exponentiate to get log-normal 
catch_res <- exp(catch_res)
### catch_res is a factor by which the numbers at age are multiplied

### for historical period, pass on real observed catch
### -> remove deviation
catch_res[, dimnames(catch_res)$year <= 2017] <- 1

if (isTRUE(verbose)) plot(catch_res, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

### ------------------------------------------------------------------------ ###
### check SAM ####
### ------------------------------------------------------------------------ ###
# 
# stk_tmp <- window(stk_stf, end = 2018)
# catch.wt(stk_tmp)[,ac(2018)] <- landings.wt(stk_tmp)[,ac(2018)] <-
#   discards.wt(stk_tmp)[,ac(2018)] <- NA
# stk_tmp <- stk_tmp[,,,,, 1:10]
# idx_tmp <- window(idx, end = 2018)
# idx_tmp[[2]] <- window(idx_tmp[[2]], end = 2017)
# idx_tmp <- lapply(idx_tmp, FLCore::iter, 1:10)
# 
# fit3 <- FLR_SAM(stk = stk_tmp, 
#                idx = idx_tmp, conf = cod4_conf_sam)
# stk3 <- SAM2FLStock(fit3)
# plot(iter(FLStocks(cod4 = stk, sim = stk3), 1))
# 
# 

### ------------------------------------------------------------------------ ###
### save OM ####
### ------------------------------------------------------------------------ ###

### path
input_path <- paste0("input/cod4/", n, "_", n_years, "/")
dir.create(input_path)
### stock
saveRDS(stk_fwd, file = paste0(input_path, "stk.rds"))
### stock recruitment
saveRDS(sr, file = paste0(input_path, "sr.rds"))
### recruitment residuals
saveRDS(sr_res, file = paste0(input_path, "sr_res.rds"))
### surveys
saveRDS(idx, file = paste0(input_path, "idx.rds"))
saveRDS(idx_dev, file = paste0(input_path, "idx_dev.rds"))
### catch noise
saveRDS(catch_res, file = paste0(input_path, "catch_res.rds"))
### process error
saveRDS(proc_res, file = paste0(input_path, "proc_res.rds"))
### observed stock
saveRDS(stk_oem, file = paste0(input_path, "stk_oem.rds"))
### sam initial parameters
saveRDS(sam_initial, file = paste0(input_path, "sam_initial.rds"))
### sam configuration
saveRDS(cod4_conf_sam_no_mult, file = paste0(input_path, "cod4_conf_sam_no_mult"))
### catch numbers
saveRDS(catch_n, file = paste0(input_path, "catch_n.rds"))
save.image(file = paste0(input_path, "image.RData"))

# stk_fwd <- readRDS(file = paste0(input_path, "stk.rds"))
# sr <- readRDS(file = paste0(input_path, "sr.rds"))
# sr_res <- readRDS(file = paste0(input_path, "sr_res.rds"))
# idx <- readRDS(file = paste0(input_path, "idx.rds"))
# idx_dev <- readRDS(file = paste0(input_path, "idx_dev.rds"))
# catch_res <- readRDS(file = paste0(input_path, "catch_res.rds"))
# proc_res <- readRDS(file = paste0(input_path, "proc_res.rds"))
# stk_oem <- readRDS(file = paste0(input_path, "stk_oem.rds"))
# sam_initial <- readRDS(file = paste0(input_path, "sam_initial.rds"))
# cod4_conf_sam_no_mult <- readRDS(file = paste0(input_path, 
#                                                "cod4_conf_sam_no_mult"))

### ------------------------------------------------------------------------ ###
### prepare objects for new a4a standard mse package ####
### ------------------------------------------------------------------------ ###
### https://github.com/flr/mse

### save workspace to start from here
# save.image(file = "input/cod4/image_10.RData")
# load(file = "input/cod4/image_10.RData")

### reference points
refpts_mse <- list(Btrigger = 150000,
                   Ftrgt = 0.31,
                   Fpa = 0.39,
                   Bpa = 150000,
                   Blim = 107000)
### some specifications for short term forecast with SAM
cod4_stf_def <- list(fwd_yrs_average = -3:0,
                     fwd_yrs_rec_start = 1998,
                     fwd_yrs_sel = -3:-1,
                     fwd_yrs_lf_remove = -2:-1,
                     fwd_splitLD = TRUE)

### some arguments (passed to mp())
args <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
                y0 = dims(stk_fwd)$minyear, ### first data year
                iy = yr_data, ### first simulation (intermediate) year
                nsqy = 3, ### not used, but has to provided
                nblocks = 1, ### block for parallel processing
                seed = 1 ### random number seed before starting MSE
)

### operating model
om <- FLom(stock = stk_fwd, ### stock 
           sr = sr, ### stock recruitment and precompiled residuals
           projection = mseCtrl(method = fwd_WKNSMSE, 
                                args = list(maxF = 2,
                                            ### process noise on stock.n
                                            proc_res = "fitted"
                                ))
)

### observation (error) model
oem <- FLoem(method = oem_WKNSMSE,
             observations = list(stk = stk_oem, idx = idx), 
             deviances = list(stk = FLQuants(catch.dev = catch_res), 
                              idx = idx_dev),
             args = list(idx_timing = c(0, -1),
                         catch_timing = -1,
                         use_catch_residuals = TRUE, 
                         use_idx_residuals = TRUE,
                         use_stk_oem = TRUE))
### implementation error model (banking and borrowing)
# iem <- FLiem(method = iem_WKNSMSE, 
#              args = list(BB = TRUE))

### default management
ctrl_obj <- mpCtrl(list(
  est = mseCtrl(method = SAM_wrapper,
                     args = c(### short term forecast specifications
                       forecast = TRUE, 
                       fwd_trgt = "fsq", fwd_yrs = 1, 
                       cod4_stf_def,
                       ### speeding SAM up
                       newtonsteps = 0, rel.tol = 0.001,
                       par_ini = list(sam_initial),
                       track_ini = TRUE, ### store ini for next year
                       ### SAM model specifications
                       conf = list(cod4_conf_sam_no_mult),
                       parallel = FALSE ### TESTING ONLY
                     )),
  phcr = mseCtrl(method = phcr_WKNSMSE,
                      args = refpts_mse),
  hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
  isys = mseCtrl(method = is_WKNSMSE, 
                    args = c(hcrpars = list(refpts_mse),
                             ### for short term forecast
                             fwd_trgt = list(c("fsq", "hcr")), fwd_yrs = 2,
                             cod4_stf_def#,
                             ### TAC constraint
                             #TAC_constraint = TRUE,
                             #lower = -Inf, upper = Inf,
                             #Btrigger_cond = FALSE,
                             ### banking and borrowing 
                             #BB = TRUE,
                             #BB_check_hcr = FALSE,
                             #BB_check_fc = TRUE,
                             #BB_rho = list(c(-0.1, 0.1))
                    ))#,
  #tm = NULL
))
### additional tracking metrics
tracking_add <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")

### save mse objects
input <- list(om = om, oem = oem, ctrl = ctrl_obj,
              args = args, tracking = tracking_add)
saveRDS(object = input, 
        file = paste0(input_path, "base_run.rds"))
# input <- readRDS(paste0(input_path, "/base_run.rds"))

### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###


### run MSE
### WARNING: takes a while...
### check normal execution
# res1 <- mp(om = input$om,
#            oem = input$oem,
#            #iem = iem,
#            ctrl = input$ctrl,
#            args = input$args,
#            tracking = input$tracking)

debugonce(mp)
debugonce(goFish)
res1 <- mp(om = input$om,
           oem = input$oem,
           #iem = iem,
           ctrl = input$ctrl,
           args = input$args,
           tracking = input$tracking)



### create Rmarkdown file
# knitr::spin(hair = "OM.R", format = "Rmd", precious = TRUE, comment = c('^### ------------------------------------------------------------------------ ###$', '^### ------------------------------------------------------------------------ ###$'))