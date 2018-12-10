### ------------------------------------------------------------------------ ###
### create FLStock for cod ####
### ------------------------------------------------------------------------ ###
### base on SAM assessment

### load packages
library(FLfse)
library(ggplotFL)
library(FLAssess)
# library(mse)
### load files from package mse for easier debugging
devtools::load_all("../mse/")
library(FLash)
library(tidyr)
library(dplyr)

### source the scripts from functions folder
# invisible(lapply(list.files(path = "functions/", pattern = "*.R$", 
#                             full.names = TRUE), source))
source("a4a_mse_WKNSMSE_funs.R")

dir.create(path = "input/cod4", recursive = TRUE)
dir.create(path = "output/runs", recursive = TRUE)

### ------------------------------------------------------------------------ ###
### simulation specifications ####
### ------------------------------------------------------------------------ ###

### number of iterations/replicates
n <- 10
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

is(fit)
fit
plot(fit)

### extract model parameters and use them in the simulation as starting values
sam_initial <- sam_getpar(fit)
sam_initial$logScale <- numeric(0)

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
summary(fit2) / summary(fit)
### estimates and log likelihood identical, only bounds smaller

### ------------------------------------------------------------------------ ###
### create FLStock ####
### ------------------------------------------------------------------------ ###
### create template with 1 iteration
### cod4_stk2 is used as template, i.e. the input values (catch) include
### the catch multiplier, 
### the results (stock numbers & harvest) are used from the real WGNSSK fit
stk <- SAM2FLStock(object = fit, stk = cod4_stk2)
summary(stk)

### set units
units(stk)[1:17] <- as.list(c(rep(c("t", "1000", "kg"), 4),
                              "", "", "f", "", ""))
plot(stk)

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
set.seed(2)
uncertainty <- SAM_uncertainty(fit = fit, n = n, print_screen = FALSE)
### add noise to stock
stock.n(stk)[] <- uncertainty$stock.n
stock(stk)[] <- computeStock(stk)
### add noise to F
harvest(stk)[] <- uncertainty$harvest

### catch noise added later


plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

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
### stock weights, M etc in 2018 based on a three year average to enable calculation of SSB
### although SAM estimates F in 2018, this is not reported or taken forward into forcasts by the WG
# stk_stf2017 <- stf(window(stk, end = 2017), n_years + 1)
stk_stf2018 <- stf(window(stk, end = 2018), n_years)
### use all available data
stk_stf <- stk_stf2018

### Resample weights, maturity and natural mortality from the last 5 years (2013-2017)
# set up an array with one resampled year for each projection year (including intermediate year) and replicate
# use the same resampled year for all biological parameters
# this is the approach used in eqsim for North Sea cod
r.bio <- array(sample(2013:2017, (n_years + 1) * n, TRUE), c(n_years + 1, n))

# loop to fill in replicates
# can probably be coded more elegantly than this!
for (iter in 1:n){
  
  # 2018 weights and natural mortality are averages in the assessment
  # so use a resampled value for 2018
  catch.wt(stk_stf)[, ac(2018:(2018 + n_years)),,,,iter] <- catch.wt(stk)[, ac(r.bio[,iter]),,,,1]
  discards.wt(stk_stf)[, ac(2018:(2018 + n_years)),,,,iter] <- discards.wt(stk)[, ac(r.bio[,iter]),,,,1]
  landings.wt(stk_stf)[, ac(2018:(2018 + n_years)),,,,iter] <- landings.wt(stk)[, ac(r.bio[,iter]),,,,1]
  stock.wt(stk_stf)[, ac(2018:(2018 + n_years)),,,,iter] <- stock.wt(stk)[, ac(r.bio[,iter]),,,,1]
  m(stk_stf)[, ac(2018:(2018 + n_years)),,,,iter] <- m(stk)[, ac(r.bio[,iter]),,,,1]
  
  # 2018 maturiy is based on actual data, so use this value for 2018
  # and resample from 2019 onwards
  mat(stk_stf)[, ac(2019:(2018 + n_years)),,,,iter] <- mat(stk)[, ac(r.bio[2:nrow(r.bio),iter]),,,,1]
}

# This may get overwritten later, but usual approach is to assume last data year F in intermediate year
harvest(stk_stf)[, ac(2018)] <- harvest(stk)[, ac(2017)]

plot(stk_stf)

### ------------------------------------------------------------------------ ###
### stock recruitment ####
### ------------------------------------------------------------------------ ###
### fit hockey-stick model
### get residuals from smoothed residuals

### use only data from 1997 and later
sr <- as.FLSR(window(stk_stf, start = 1997), model = "segreg")
### fit model individually to each iteration and suppress output to screen
suppressWarnings(. <- capture.output(sr <- fmle(sr)))

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

# Check extent of autocorrelation (need to check what goes into this)
#acf(window(stock.n(stk_orig)[1], start = 1998))
#acf(rec(sr)[!is.na(rec(sr))])

### years with missing residuals
# NW: dimnames produces NULL for me
# yrs_res <- dimnames(sr)$year[which(is.na(iterMeans(rec(sr))))]
yrs_res <- colnames(rec(sr))[which(is.na(iterMeans(rec(sr))))]

### go through iterations and create residuals
### use kernel density to create smooth distribution of residuals
### and sample from this distribution
res_new <- foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
                   .errorhandling = "pass") %do% {
                     
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
residuals(sr) <- exp(residuals(sr))
sr_res <- residuals(sr)

plot(sr_res)


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

### create stock for MSE simulation
stk_fwd <- stk_stf
### insert values for 2018
stk_fwd[, ac(2018)] <- stk_int[, ac(2018)]
### insert stock number for 2019 in order to calculate SSB at beginning of 
### 2019
stock.n(stk_fwd)[, ac(2019)] <- stock.n(stk_int)[, ac(2019)]

#all.equal(window(stk_fwd, end = 2018), window(stk_stf, end = 2018))

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
### calculate index values
idx <- calc_survey(stk = stk_fwd, idx = idx)

### create deviances for indices
### first, get template
idx_dev <- lapply(idx, index)
### create random noise based on sd
set.seed(2)
for (idx_i in seq_along(idx_dev)) {
  ### insert sd
  idx_dev[[idx_i]][] <- uncertainty$survey_sd[[idx_i]]
  ### noise
  idx_dev[[idx_i]][] <- stats::rnorm(n = length(idx_dev[[idx_i]]),
                                   mean = 0, sd = idx_dev[[idx_i]])
  #idx_dev[[idx_i]] <- exp(idx_dev[[idx_i]]) # Already exp(logSdLogObs)?
}


### compare to original survey(s)
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

### ------------------------------------------------------------------------ ###
### catch noise ####
### ------------------------------------------------------------------------ ###
### take estimates from sam: uncertainty$catch_sd is "logSdLogObs"
### assume catch observed by SAM in projection is log-normally distributed
### around operating model catch

### create noise for catch
set.seed(2)
catch_res <- catch.n(stk_fwd) %=% 0 ### template FLQuant
catch_res[] <- stats::rnorm(n = length(catch_res), mean = 0, 
                            sd = uncertainty$catch_sd)
#catch_res <- exp(catch_res) # Already exp(logSdLogObs)?
plot(catch_res)


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

### stock
saveRDS(stk_fwd, file = "input/cod4/stk.rds")
### stock recruitment
saveRDS(sr, file = "input/cod4/sr.rds")
### recruitment residuals
saveRDS(sr_res, file = "input/cod4/sr_res.rds")
### surveys
saveRDS(idx, file = "input/cod4/idx.rds")
### catch noise
saveRDS(catch_res, file = "input/cod4/catch_res.rds")

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
                   Bpa = 150000)
### some specifications for short term forecast with SAM
cod4_stf_def <- list(fwd_yrs_average = -3:0,
                     fwd_yrs_rec_start = 1998,
                     fwd_yrs_sel = -3:-1,
                     fwd_yrs_lf_remove = -2:-1,
                     fwd_splitLD = TRUE)

### some arguments (passed to mp())
genArgs <- list(fy = yr_data + n_years - 1, ### final simulation year
                y0 = yr_data, ### first simulation year
                iy = yr_data,
                nsqy = 3, ### not used, but has to provided
                nblocks = 1, ### block for parallel processing
                seed = 1 ### random number seed before starting MSE
)

### operating model
om <- FLom(stock = stk_fwd, ### stock 
           sr = sr ### stock recruitment and precompiled residuals
           )

### observation (error) model
oem <- FLoem(method = oem_WKNSMSE,
             observations = list(stk = stk_fwd, idx = idx), 
             deviances = list(stk = FLQuants(catch.dev = catch_res), 
                              idx = idx_dev),
             args = list(idx_timing = c(0, -1),
                         catch_timing = -1,
                         use_catch_residuals = TRUE, 
                         use_idx_residuals = TRUE))

### default management
ctrl_obj <- mpCtrl(list(
  ctrl.est = mseCtrl(method = SAM_wrapper,
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
  ctrl.phcr = mseCtrl(method = phcr_WKNSMSE,
                      args = refpts_mse),
  ctrl.hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
  ctrl.is = mseCtrl(method = is_WKNSMSE, 
                    args = c(hcrpars = list(refpts_mse),
                             ### for short term forecast
                             fwd_trgt = c("fsq", "hcr"), fwd_yrs = 2,
                             cod4_stf_def#,
                             ### TAC constraint
                             #TAC_constraint = TRUE,
                             #lower = -Inf, upper = Inf,
                             #Btrigger_cond = FALSE,
                             ### banking and borrowing 
                             #BB = TRUE,
                             #BB_conditional = TRUE,
                             #BB_rho = list(c(-0.1, 0.1))
                             )),
  ctrl.tm = NULL
))
### additional tracking metrics
tracking_add <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")

### try running mse
# library(doParallel)
# cl <- makeCluster(10)
# registerDoParallel(cl)
# debugonce(mp)
# debugonce(is_WKNSMSE)
# res1 <- mp(om = om,
#            oem = oem,
#            ctrl.mp = ctrl_obj,
#            genArgs = genArgs,
#            tracking = tracking_add)

# ### without BB
# ctrl_obj2 <- ctrl_obj
# ctrl_obj2$ctrl.is@args$BB <- FALSE
# res2 <- mp(om = om,
#            oem = oem,
#            ctrl.mp = ctrl_obj2,
#            genArgs = genArgs,
#            tracking = tracking_add)

### run MSE
### WARNING: takes a while...
### check normal execution
res1 <- mp(om = om,
           oem = oem,
           ctrl.mp = ctrl_obj,
           genArgs = genArgs,
           tracking = tracking_add)
### check mpParallel function
resp1 <- mpParallel(om = om,
            oem = oem,
            ctrl.mp = ctrl_obj,
            genArgs = genArgs,
            tracking = tracking_add)

### split into 2 parts
genArgs$nblocks <- 2
resp2 <- mpParallel(om = om,
                    oem = oem,
                    ctrl.mp = ctrl_obj,
                    genArgs = genArgs,
                    tracking = tracking_add)
### execute in parallel
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
### load packages and additional functions into workers
clusterEvalQ(cl = cl, expr = {
  library(mse)
  library(FLash)
  library(FLfse)
  library(stockassessment)
  library(foreach)
  library(doRNG)
  source("a4a_mse_WKNSMSE_funs.R")
})
### run MSE
resp3 <- mpParallel(om = om,
                    oem = oem,
                    ctrl.mp = ctrl_obj,
                    genArgs = genArgs,
                    tracking = tracking_add)
### try reproducible parallel execution
library(doRNG)
registerDoRNG(123) 
resp4 <- mpParallel(om = om,
                    oem = oem,
                    ctrl.mp = ctrl_obj,
                    genArgs = genArgs,
                    tracking = tracking_add)
registerDoRNG(123) 
resp5 <- mpParallel(om = om,
                    oem = oem,
                    ctrl.mp = ctrl_obj,
                    genArgs = genArgs,
                    tracking = tracking_add)

### create Rmarkdown file
# knitr::spin(hair = "OM.R", format = "Rmd", precious = TRUE, comment = c('^### ------------------------------------------------------------------------ ###$', '^### ------------------------------------------------------------------------ ###$'))