### ------------------------------------------------------------------------ ###
### create FLStock for cod ####
### ------------------------------------------------------------------------ ###
### base on SAM assessment

### load packages
library(FLfse)
library(ggplotFL)
library(FLAssess)
library(FLash)
library(tidyr)
library(dplyr)

### source the scripts from functions folder
invisible(lapply(list.files(path = "functions/", pattern = "*.R$", 
                            full.names = TRUE), source))

dir.create(path = "input/cod4", recursive = TRUE)
dir.create(path = "output/runs", recursive = TRUE)

### ------------------------------------------------------------------------ ###
### simulation specifications ####
### ------------------------------------------------------------------------ ###

### number of iterations
n <- 100
### number of years
n_years <- 30

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
fit2 <- FLR_SAM(stk = cod4_stk2, idx = cod4_idx, 
                conf = cod4_conf_sam[!names(cod4_conf_sam) %in% 
                                       c("noScaledYears", "keyScaledYears",
                                         "keyParScaledYA")])
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
units(stk)[1:17] <- as.list(c(rep(c("tonnes", "thousands", "kg"), 4),
                              "NA", "NA", "f", "NA", "NA"))
plot(stk)

### save for later comparison
stk_orig <- stk

### ------------------------------------------------------------------------ ###
### add uncertainty ####
### ------------------------------------------------------------------------ ###
### first approach: use variance-covariance


### add iteration dimension
stk <- propagate(stk, n)
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
library(tmbstan) # cran package

### create template stock for storing results
MCMC_iter <- n
MCMC_warmup <- 1000
stk_MCMC <- propagate(stk, MCMC_iter)

### run MCMC
system.time(mcmc <- tmbstan(fit$obj, chains = 1, iter = MCMC_warmup + MCMC_iter, 
                            warmup = MCMC_warmup, 
                            seed = 1, control = list(max_treedepth = 15)))
### extract 
mc <- extract(mcmc, inc_warmup = FALSE, permuted = FALSE)

table(gsub(x = dimnames(mc)$parameters, pattern = "\\[[0-9]{1,}\\]", 
           replacement = ""))
# itrans_rho         logF      logFpar         logN     logScale logSdLogFsta 
#          1          336            9          336           13            2 
# logSdLogN  logSdLogObs         lp__ 
#         2            7            1 
### gives N and F @age
### scale for 


### find positions for results in MCMC
nms <- dimnames(mc)$parameters
F_pos <- grep(x = nms, pattern = "logF\\[[0-9]{1,3}\\]$")
N_pos <- grep(x = nms, pattern = "logN\\[[0-9]{1,3}\\]$")
factor_pos <- grep(x = nms, pattern = "logScale\\[[0-9]{1,2}\\]$")

### in the MCMC results age and years are mixed within the same row,
### each row represents one iteration
### this needs to be reformatted to be useful...

### fishing mortality:
harvest(stk_MCMC)[] <- aperm(array(data = c(exp(mc[,, F_pos])),
                                   dim = dim(harvest(stk_MCMC))[c(6, 1:5)]),
                             perm = c(2:6, 1))

### stock numbers at age
stock.n(stk_MCMC)[] <- aperm(array(data = c(exp(mc[,, N_pos])),
                                   dim = dim(stock.n(stk_MCMC))[c(6, 1:5)]),
                             perm = c(2:6, 1))
stock(stk_MCMC) <- computeStock(stk_MCMC)

### catch factor
### get ages and years
ages <- fit$conf$minAge:fit$conf$maxAge
yrs <- fit$conf$keyScaledYears
### get catch factors
catch_factor <- lapply(split(exp(mc[,, factor_pos]), seq(MCMC_iter)), 
                       function(x) {
  x[t(fit$conf$keyParScaledYA + 1)]
})
catch_factor <- unlist(catch_factor)
### coerce into FLQuant
catch_factor <- FLQuant(catch_factor,
                      dimnames = dimnames(catch.n(stk_MCMC[ac(ages), ac(yrs)])))
### multiply catch numbers
catch.n(stk_MCMC)[ac(ages), ac(yrs)] <- catch_factor * 
  catch.n(stk_MCMC)[ac(ages), ac(yrs)]
### split into landings and discards, based on landing fraction
lfrac <- propagate((landings.n(stk)[ac(ages), ac(yrs)] / 
                      catch.n(stk)[ac(ages), ac(yrs)]), MCMC_iter)
landings.n(stk_MCMC)[ac(ages), ac(yrs)] <- catch.n(stk_MCMC)[ac(ages), ac(yrs)] *
  lfrac
discards.n(stk_MCMC)[ac(ages), ac(yrs)] <- catch.n(stk_MCMC)[ac(ages), ac(yrs)] *
  (1 - lfrac)
### update stock
catch(stk_MCMC)[, ac(yrs)] <- computeCatch(stk_MCMC)[, ac(yrs)]
landings(stk_MCMC)[, ac(yrs)] <- computeLandings(stk_MCMC)[, ac(yrs)]
discards(stk_MCMC)[, ac(yrs)] <- computeDiscards(stk_MCMC)[, ac(yrs)]

### plot MCMC stock
plot(stk_MCMC)
### compare with original SAM fit
plot(FLStocks(MCMC = stk_MCMC, original = stk_orig))
### compare with first uncertainty approach
plot(FLStocks(MCMC = stk_MCMC, original = stk_orig, Cov = stk))

### ------------------------------------------------------------------------ ###
### try SAM internal "simstudy" ####
### ------------------------------------------------------------------------ ###
### "Simulate data from fitted model and re-estimate from each run"

set.seed(0)
system.time(fits <- simstudy(fit = fit, nsim = n))
class(fits) <- "sam_list"

stk_sim <- SAM2FLStock(object = fits, stk = cod4_stk2)
plot(FLStocks(Cov = stk, simstudy = stk_sim,
              original = stk_orig))

### ------------------------------------------------------------------------ ###
### extend stock for MSE simulation ####
### ------------------------------------------------------------------------ ###

### special case for NS cod
### catch & catch weights until 2017,
### maturity, stock weights, stock & F etc until 2018
stk_stf2017 <- stf(window(stk, end = 2017), n_years + 1)
stk_stf2018 <- stf(window(stk, end = 2018), n_years)
### 
stk_stf <- stk_stf2018
### maturity data 2018 is already average of previous years
### use this in projection
mat(stk_stf)[, ac((2018 + 1):(2018 + n_years))] <- mat(stk_stf[, ac(2018)])

### last 3 data years for catch weights
catch.wt(stk_stf)[, ac(2018:(2018 + n_years))] <- 
  apply(catch.wt(stk)[, ac(2015:2017)], c(1, 6), mean)
landings.wt(stk_stf)[, ac(2018:(2018 + n_years))] <- 
  apply(landings.wt(stk)[, ac(2015:2017)], c(1, 6), mean)
discards.wt(stk_stf)[, ac(2018:(2018 + n_years))] <- 
  apply(discards.wt(stk)[, ac(2015:2017)], c(1, 6), mean)

plot(stk_stf)

### ------------------------------------------------------------------------ ###
### stock recruitment ####
### ------------------------------------------------------------------------ ###

# Truncated stock objects for defining different recruitment regimes
# Note EQSIM would end with the data rather than intermediate year 
#stk_trun<-window(stk, start=1987, end=2017)
stk_trun2<-window(stk, start=1997, end=2018)

### create stock recruitment model: segmented regression (hockey-stick)
sr <- as.FLSR(stk_trun2, model = "segreg")
### fit model individually to each iteration and suppress output to screen
suppressWarnings(. <- capture.output(sr <- fmle(sr)))

### plot models
plot(sr)

as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb)) %>%
  mutate(age = NULL,
         year = ifelse(qname == "SSB", year + 1, year)) %>%
  spread(key = qname, value = data) %>%
  ggplot() +
  geom_point(aes(x = SSB, y = rec, group = iter), 
             alpha = 0.5, colour = "grey", shape = 1) +
  geom_line(aes(x = SSB, y = fitted, group = iter)) +
  theme_bw()


### create residuals for projection
set.seed(1)
sr_res <- window(residuals(sr), end = dims(residuals(sr))$maxyear + n_years)
### calculate based on sd from residuals of model fit
sr_res_new <- rnorm(n, 
  FLQuant(0, dimnames = list(year = (yr_data + 1):(yr_data + n_years))),
  mean(c(apply(residuals(sr), 6, sd))))
sr_res[, ac((yr_data + 1):(yr_data + n_years))] <- sr_res_new
### exponentiate
sr_res <- exp(sr_res)

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

### get catchability for surveys and add random noise
idx <- lapply(seq_along(idx), function(x) {
  
  ### select index
  idx_i <- idx[[x]]
  
  ### age range
  ages_i <- range(idx_i)[["min"]]:range(idx_i)[["max"]]
  ### year range
  years_i <- range(idx_i)[["minyear"]]:range(idx_i)[["maxyear"]]
  ### survey timing
  time_i <- mean(range(idx_i)[c("startf", "endf")])
  
  ### get catchability from SAM
  q <- uncertainty$survey_catchability[[x]]
  
  ### get survey sd & adapt dimensions
  idx_sd <- window(uncertainty$survey_sd[[x]], start = min(years_i), 
                   end = max(years_i))
  idx_sd[] <- uncertainty$survey_sd[[x]]
  
  ### create random noise, based on survey sd at age
  set.seed(2)
  noise <- idx_sd %=% 0
  noise[] <- stats::rnorm(n = length(idx_sd), mean = 0, sd = idx_sd)
  noise <- exp(noise)
  
  ### remove noise from first iteration
  noise[,,,,, 1] <- 1
  
  ### add noise to catchability
  index.q(idx_i)[] <- q
  index.q(idx_i) <- index.q(idx_i) * noise
  
  return(idx_i)
  
})
names(idx) <- names(cod4_idx)
idx <- FLIndices(idx)

### initialize index values
idx <- calc_survey(stk = stk_fwd, idx = idx)


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

### create noise for catch
set.seed(2)
catch_res <- catch.n(stk_fwd) %=% 0 ### template FLQuant
catch_res[] <- stats::rnorm(n = length(catch_res), mean = 0, 
                            sd = uncertainty$catch_sd)
catch_res <- exp(catch_res)
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
### remove fishy first iteration ####
### ------------------------------------------------------------------------ ###
### lengthy code in order to get iteration dimnames 1:n and not 2:n+1

# stk_fwd2 <- stk_fwd[,,,,, 1:n]
# stk_fwd2[] <- stk_fwd[,,,,, 2:(n + 1)]
# 
# sr2 <- FLCore::iter(sr, 1:n)
# sr2[] <- FLCore::iter(sr, 2:(n + 1))
# 
# sr_res2 <- FLCore::iter(sr_res, 1:n)
# sr_res2[] <- FLCore::iter(sr_res, 2:(n + 1))
# 
# idx2 <- idx
# idx2[[1]] <- FLCore::iter(idx[[1]], 1:n)
# idx2[[1]][] <- FLCore::iter(idx[[1]], 2:(n + 1))
# idx2[[2]] <- FLCore::iter(idx[[2]], 1:n)
# idx2[[2]][] <- FLCore::iter(idx[[2]], 2:(n + 1))

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

# knitr::spin(hair = "OM.R", format = "Rmd", precious = TRUE, comment = c('^### ------------------------------------------------------------------------ ###$', '^### ------------------------------------------------------------------------ ###$'))