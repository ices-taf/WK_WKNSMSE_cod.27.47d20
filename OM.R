### ------------------------------------------------------------------------ ###
### create FLStock for cod ####
### ------------------------------------------------------------------------ ###
### base on SAM assessment

### load packages
library(FLfse)
library(ggplotFL)
library(FLAssess)

### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))

### ------------------------------------------------------------------------ ###
### simulation specifications ####
### ------------------------------------------------------------------------ ###

### number of iterations
n <- 10
### number of years
n_years <- 30

### last data year
yr_data <- 2018


### ------------------------------------------------------------------------ ###
### fit SAM ####
### ------------------------------------------------------------------------ ###
### use input data provided in FLfse
fit <- FLR_SAM(stk = cod4_stk, idx = cod4_idx, conf = cod4_conf_sam)
is(fit)
fit

### ------------------------------------------------------------------------ ###
### create FLStock ####
### ------------------------------------------------------------------------ ###
### the stock contains n+1 iterations
### the first one is the estimate from SAM, the remaining ones are resampled

### creat template with 1 iteration
stk <- SAM2FLStock(object = fit)
summary(stk)

### add iteration dimension
stk <- propagate(stk, n + 1)
dim(stk)

### add uncertainty estimated by SAM as iterations
set.seed(1)
uncertainty <- SAM_uncertainty(fit = fit, n = n)
### add noise to stock
stock.n(stk) <- uncertainty$stock.n
stock(stk) <- computeStock(stk)
### add noise to F
harvest(stk) <- uncertainty$harvest

###
### assume perfect knowledge of catch for now
### ...

### set units
units(stk)[1:17] <- as.list(c(rep(c("tonnes", "thousands", "kg"), 4),
                              "NA", "NA", "f", "NA", "NA"))

### ------------------------------------------------------------------------ ###
### extend stock for MSE simulation ####
### ------------------------------------------------------------------------ ###

### special case for NS cod
### catch & catch weights until 2017,
### maturity, stock weights, stock & F etc until 2018
stk_stf2017 <- stf(window(stk, end = 2017), n_years - 1)
stk_stf2018 <- stf(window(stk, end = 2018), n_years)
### 
stk_stf <- stk_stf2018
### last 2 data years for catch weights
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

### create stock recruitment model: segmented regression (hockey-stick)
sr <- as.FLSR(stk, model = "segreg")
### fit model individually to each iteration
sr <- fmle(sr)

plot(sr)
plot(sr@fitted ~ sr@ssb)

### create residuals for projection
set.seed(1)
sr_res <- window(residuals(sr), end = dims(residuals(sr))$maxyear + n_years)
### calculate based on sd from residuals of model fit
sr_res_new <- rnorm(n + 1, 
  FLQuant(0, dimnames = list(year = (yr_data + 1):(yr_data + n_years))),
  mean(c(apply(residuals(sr), 6, sd))))
sr_res[, ac((yr_data + 1):(yr_data + n_years))] <- sr_res_new
### remove residuals from first (deterministic) iteration
sr_res[, ac((yr_data + 1):(yr_data + n_years)),,,, 1] <- 0

plot(sr_res)

### ------------------------------------------------------------------------ ###
### indices ####
### ------------------------------------------------------------------------ ###
### use real FLIndices object as template (included in FLfse)
idx <- cod4_idx
### extend for simulation period
idx <- window(idx, end = yr_data + n_years)
### add iterations
idx <- lapply(idx, propagate, n + 1)

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
idx <- calc_survey(stk = stk_stf, idx = idx)


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
# fit2 <- FLR_SAM(stk = stk_tmp, 
#                idx = idx_tmp, conf = cod4_conf_sam)
# stk2 <- SAM2FLStock(fit2)
# plot(stk2)
# plot(stk)
# 
# 

### ------------------------------------------------------------------------ ###
### save OM ####
### ------------------------------------------------------------------------ ###


### stock
saveRDS(stk_stf, file = "input/cod4/stk.rds")
### stock recruitment
saveRDS(sr, file = "input/cod4/sr.rds")
### recruitment residuals
saveRDS(sr_res, file = "input/cod4/sr_res.rds")
### surveys
saveRDS(idx, file = "input/cod4/idx.rds")
