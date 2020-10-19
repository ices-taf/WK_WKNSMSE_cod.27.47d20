### ------------------------------------------------------------------------ ###
### R script to run WKNSMSE cod MSE on HPC ####
### ------------------------------------------------------------------------ ###
### This is designed to be called by a job submission script
### run_mse.qsub for systems using PBS and the qsub commands
### run_mse.bsub for system using LSF and the bsub commands


### ------------------------------------------------------------------------ ###
### load arguments from job script ####
### ------------------------------------------------------------------------ ###

### load arguments
args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  
  ### parallelisation environment
  if (!exists("par_env")) par_env <- 2
  if (!exists("n_workers")) n_workers <- 1
  ### some basic arguments
  if (!exists("saveMP")) saveMP <- TRUE
  if (!exists("calc_stats")) calc_stats <- TRUE
  
} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
library(FLfse)
library(stockassessment)
library(ggplotFL)
#library(FLAssess)
library(mse)
### load files from package mse for easier debugging
#devtools::load_all("../mse/")
library(FLash)
library(tidyr)
library(dplyr)
library(foreach)
library(doRNG)

### load additional functions
source("a4a_mse_WKNSMSE_funs.R")

### the original WKNSMSE was run with R 3.5
### for exact reproducibility in R 3.6, the random number generation must be
### be changed
if (getRversion() >= 3.6) RNGkind(sample.kind = "Rounding")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###
### par_env=1 -> MPI (Rmpi, DoMPI)
### par_env=2 -> DoParallel

if (par_env == 1) {
  
  library(doMPI)
  cl <- startMPIcluster()
  registerDoMPI(cl)
  cl_length <- cl$workerCount
  
} else if (par_env == 2) {
  
  library(doParallel)
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  cl_length <- length(cl)
  
}

### load packages and functions into workers
. <- foreach(i = seq(cl_length)) %dopar% {
  #devtools::load_all("../mse/")
  library(mse)
  library(FLash)
  library(FLfse)
  library(stockassessment)
  library(foreach)
  library(doRNG)
  source("a4a_mse_WKNSMSE_funs.R")
  if (getRversion() >= 3.6) RNGkind(sample.kind = "Rounding")
}

### set random seed for reproducibility
library(doRNG)
registerDoRNG(123)

### ------------------------------------------------------------------------ ###
### load data for MSE ####
### ------------------------------------------------------------------------ ###

### data path
OM_alt <- "cod4"
if (exists("OM")) {
  if (OM > 0) OM_alt <- paste0("cod4_alt", OM)
}
path_data <- paste0("input/", OM_alt, "/", iters, "_", years, "/")

### load input objects
input <- readRDS(paste0(path_data, "base_run.rds"))

### modify input for running in parallel
input$args$nblocks <- nblocks

### ------------------------------------------------------------------------ ###
### set up HCR & options ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### set HCR option: A, B, C
if (exists("HCRoption")) {
  
  input$ctrl$hcr@args$option <- switch(HCRoption, 
                                               "1" = "A", 
                                               "2" = "B", 
                                               "3" = "C",
                                               "4" = "A",
                                               "5" = "B",
                                               "6" = "C")
  
  cat(paste0("\nSetting custom HCR option: HCRoption = ", HCRoption, 
             " => HCR ", input$ctrl$hcr@args$option, "\n\n"))
  
} else {
  
  cat(paste0("\nUsing default HCR option: HCR ", 
             input$ctrl$hcr@args$option, "\n\n"))
  HCRoption <- 0
  
}


### ------------------------------------------------------------------------ ###
### set HCR parameters 
if (HCRoption %in% 1:6) {

  ### 1-451
  hcr_vals <- expand.grid(
    Ftrgt = seq(from = 0.1, to = 0.5, by = 0.01),
    Btrigger = seq(from = 110000, to = 210000, by = 10000))
  ### additional combinations after finding yield maximum
  ### 452-456
  comb_max <- switch(HCRoption, 
                     "1" = c(170000, 0.38), 
                     "2" = c(160000, 0.38), 
                     "3" = c(170000, 0.38),
                     "4" = c(190000, 0.40),
                     "5" = c(130000, 0.36),
                     "6" = c(140000, 0.36))
  hcr_vals <- rbind(hcr_vals,
                    expand.grid(Ftrgt = c(comb_max[2], comb_max[2]*0.9, 
                                          comb_max[2]*1.1, 0.198, 0.46),
                                Btrigger = comb_max[1]))
  
}

### ------------------------------------------------------------------------ ###
### loop through several options? ####
### ------------------------------------------------------------------------ ###


### default HCR
if (!exists("HCR_comb")) {
  HCR_comb <- 186
}
  
. <- foreach(HCR_comb_i = HCR_comb) %dopar% {
  
  ### set Btrigger
  Btrigger <- hcr_vals[HCR_comb_i, "Btrigger"]
  input$ctrl$phcr@args$Btrigger <- Btrigger
  input$ctrl$isys@args$hcrpars$Btrigger <- Btrigger
  
  ### set Ftrgt
  Ftrgt <- hcr_vals[HCR_comb_i, "Ftrgt"]
  input$ctrl$phcr@args$Ftrgt <- Ftrgt
  input$ctrl$isys@args$hcrpars$Ftrgt <- Ftrgt
  
  cat(paste0("\nSetting custom Btrigger/Ftrgt values.\n",
             "Using HCR_comb_i = ", HCR_comb_i, "\n",
             "Ftrgt = ", Ftrgt, "\n",
             "Btrigger = ", Btrigger, "\n\n"))
    
  ### try uniform grid search with iteration-specific Ftrgt/Btrgigger combination
  if (exists("grid_search")) {
    if (isTRUE(as.logical(grid))) {
      hcr_vals <- expand.grid(Btrigger = seq(from = 100000, to = 200000, 
                                             length.out = 25),
                              Ftrgt = seq(from = 0.1, to = 0.49, length.out = 40))
      ### set Btrigger
      Btrigger <- hcr_vals$Btrigger
      input$ctrl$phcr@args$Btrigger <- Btrigger
      input$ctrl$isys@args$hcrpars$Btrigger <- Btrigger
      ### set Ftrgt
      Ftrgt <- hcr_vals$Ftrgt
      input$ctrl$phcr@args$Ftrgt <- Ftrgt
      input$ctrl$isys@args$hcrpars$Ftrgt <- Ftrgt
      cat(paste0("\nTrying uniform Btrigger/Ftrgt combinations.\n"))
      
    }
  }
  
  ### ------------------------------------------------------------------------ ###
  ### TAC constraint
  input$ctrl$isys@args$TAC_constraint <- FALSE
  ### check conditions
  ### either manually requested or as part of HCR options 4-6 
  if (exists("TAC_constraint")) {
    if (isTRUE(as.logical(TAC_constraint))) {
      input$ctrl$isys@args$TAC_constraint <- TRUE
    }
  }
  if (HCRoption %in% 4:6) {
      input$ctrl$isys@args$TAC_constraint <- TRUE
  }
  ### implement
  if (isTRUE(input$ctrl$isys@args$TAC_constraint)) {
      
      input$ctrl$isys@args$lower <- 80
      input$ctrl$isys@args$upper <- 125
      input$ctrl$isys@args$Btrigger_cond <- TRUE
      
      cat(paste0("\nImplementing TAC constraint.\n\n"))
      
  } else {
      
      cat(paste0("\nTAC constraint NOT implemented.\n\n"))
      
  }
  ### manual overwrite/removal of BB if requested
  if (exists("TAC_constraint")) {
    if (isTRUE(TAC_constraint == -1)) {
      input$ctrl$isys@args$TAC_constraint <- FALSE
    }
  }
  
  ### ------------------------------------------------------------------------ ###
  ### banking & borrowing
  input$ctrl$isys@args$BB <- FALSE
  input$iem <- NULL
  
  ### check conditions
  ### either manually requested or as part of HCR options 4-6 
  if (exists("BB")) {
    if (isTRUE(BB == 1)) {
      
      input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE))
      input$ctrl$isys@args$BB <- TRUE
      input$ctrl$isys@args$BB_check_hcr <- TRUE
      input$ctrl$isys@args$BB_check_fc <- TRUE
      input$ctrl$isys@args$BB_rho <- c(-0.1, 0.1)
      
    }
  
  }
  
  if (HCRoption %in% 4:6) {
      
    input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE))
    input$ctrl$isys@args$BB <- TRUE
    input$ctrl$isys@args$BB_rho <- c(-0.1, 0.1)
    input$ctrl$isys@args$BB_check_hcr <- FALSE
    input$ctrl$isys@args$BB_check_fc <- FALSE
    
    if (HCRoption %in% 4) {
      
      input$ctrl$isys@args$BB_check_hcr <- TRUE
      
    } else if (HCRoption %in% 5:6) {
      
      input$ctrl$isys@args$BB_check_fc <- TRUE
      
    }
    
  }
  ### manual overwrite/removal of BB
  if (exists("BB")) {
    if (isTRUE(BB == -1)) {
      
      input$iem <- NULL
      input$ctrl$isys@args$BB <- FALSE
      
    }
    
  }
  
  
  if (!is.null(input$iem)) {
      
    cat(paste0("\nImplementing banking and borrowing.\n\n"))
    
  } else {
    
    cat(paste0("\nBanking and borrowing NOT implemented.\n\n"))
    
  }
  
  
  ### ------------------------------------------------------------------------ ###
  ### run MSE ####
  ### ------------------------------------------------------------------------ ###
  
  ### run MSE
  res1 <- mp(om = input$om,
             oem = input$oem,
             iem = input$iem,
             ctrl = input$ctrl,
             args = input$args,
             tracking = input$tracking)
  
  ### save results
  path_out <- paste0("output/runs/cod4/", iters, "_", years)
  dir.create(path = path_out, recursive = TRUE)
  file_out <- paste0(OM_alt, "_",
                     "HCR-", input$ctrl$hcr@args$option[1],
                     "_Ftrgt-", input$ctrl$phcr@args$Ftrgt[1],
                     "_Btrigger-", input$ctrl$phcr@args$Btrigger[1],
                     "_TACconstr-", input$ctrl$isys@args$TAC_constraint[1],
                     "_BB-", input$ctrl$isys@args$BB[1]
  )
  
  if (isTRUE(saveMP))
    saveRDS(object = res1, paste0(path_out, "/MP_", file_out, ".rds"))
  
  ### ---------------------------------------------------------------------- ###
  ### stats ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(calc_stats)) {
    
    res_stats <- mp_stats(input = input, res = res1, OM = OM_alt)
    saveRDS(object = res_stats, paste0(path_out, "/stats_", file_out, ".rds"))
    
    
  }
  
}

### ------------------------------------------------------------------------ ###
### combine and plot ####
### ------------------------------------------------------------------------ ###

# ### get stock before simulation
# stk <- input$om@stock
# ### add simulated data
# stk[, dimnames(res0@stock)$year] <- res1@stock
# ### save
# saveRDS(object = stk, file = paste0("output/runs/cod4/", iters, "_", years,
#                                     "_base_full_stk.rds"))
# ### plot
# plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) + 
#   xlab("year") + geom_vline(xintercept = 2018.5) +
#   geom_hline(data = data.frame(qname = "SSB", data = 107000),
#              aes(yintercept = data), linetype = "dashed") +
#   geom_hline(data = data.frame(qname = "SSB", data = 150000),
#              aes(yintercept = data), linetype = "solid") +
#   geom_hline(data = data.frame(qname = "F", data = 0.54),
#              aes(yintercept = data), linetype = "dashed") +
#   geom_hline(data = data.frame(qname = "F", data = 0.31),
#              aes(yintercept = data), linetype = "solid") +
#   theme_bw()
# ggsave(filename = paste0("output/runs/cod4/", iters, "_", years,
#                          "_base_full_stk.png"), 
#        width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### terminate ####
### ------------------------------------------------------------------------ ###

### close R
# mpi.finalize()
### mpi.finalize() or mpi.quit() hang...
### -> kill R, the MPI processes stop afterwards

### try killing current job...
# if (par_env == 1 & exists("kill")) {
#   system("bkill $LSB_JOBID")
# }

quit(save = "no")


