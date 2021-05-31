### example run with the mse/FLasher branch and 10 iterations

library(mse)
library(FLasher)
library(FLfse)
library(stockassessment)
library(doParallel)

source("a4a_mse_WKNSMSE_funs.R")

### load input data
### (result of running OM.R)
input <- readRDS("input/cod4/10_20/base_run.rds")
### run mp with 10 iterations
### takes a few minutes
set.seed(1)
res <- do.call(mp, input)

