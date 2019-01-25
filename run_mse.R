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
  if (!exists("par_env")) par_env <- 1
  if (!exists("n_workers")) n_workers <- 1
  
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
library(FLAssess)
library(mse)
### load files from package mse for easier debugging
#devtools::load_all("../mse/")
library(FLash)
library(tidyr)
library(dplyr)

### load additional functions
source("a4a_mse_WKNSMSE_funs.R")

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
}

### set random seed for reproducibility
library(doRNG)
registerDoRNG(123)

### ------------------------------------------------------------------------ ###
### load data for MSE ####
### ------------------------------------------------------------------------ ###

### data path
path_data <- paste0("input/cod4/", iters, "_", years, "/")

### load input objects
input <- readRDS(paste0(path_data, "base_run.rds"))

### modify input for running in parallel
input$genArgs$nblocks <- nblocks

### ------------------------------------------------------------------------ ###
### set up HCR & options ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### set HCR option: A, B, C
if (exists("HCRoption")) {
  
  input$ctrl.mp$ctrl.hcr@args$option <- switch(HCRoption, 
                                               "1" = "A", 
                                               "2" = "B", 
                                               "3" = "C",
                                               "4" = "A",
                                               "5" = "B",
                                               "6" = "C")
  
  cat(paste0("\nSetting custom HCR option: HCRoption = ", HCRoption, 
             " => HCR ", input$ctrl.mp$ctrl.hcr@args$option, "\n\n"))
  
} else {
  
  cat(paste0("\nUsing default HCR option: HCR ", 
             input$ctrl.mp$ctrl.hcr@args$option, "\n\n"))
  HCRoption <- 0
  
}

### ------------------------------------------------------------------------ ###
### set HCR parameters 

if (HCRoption == 1) {
  ### create Btrigger & Ftrgt combinations
  hcr_vals <- expand.grid(
    Btrigger = seq(from = 110000, to = 190000, length.out = 5),
    Ftrgt = c(0.1, 0.2, 0.3, 0.4, 0.5))
  ### 1-25
  hcr_vals <- rbind(hcr_vals,
    expand.grid(
      Btrigger = seq(from = 110000, to = 190000, length.out = 5),
      Ftrgt = c(0.32, 0.34, 0.36, 0.38)))
  ### 26-45
  hcr_vals <- rbind(hcr_vals,
    expand.grid(
      Btrigger = c(130000, 150000, 170000),
      Ftrgt = c(0.39, 0.37)))
  ### 46-51
  hcr_vals <- rbind(hcr_vals,
    data.frame(Btrigger = c(110000), Ftrgt = c(0.35)))
  ### 52
  hcr_vals <- rbind(hcr_vals,
    expand.grid(Btrigger = 120000, Ftrgt = c(0.35, 0.36, 0.37)),
    expand.grid(Btrigger = 140000, Ftrgt = c(0.36, 0.37, 0.38)),
    expand.grid(Btrigger = 160000, Ftrgt = c(0.37, 0.38, 0.39, 0.4))
  )
  ### 53-62
  ### assume best option is Btrigger=originial Btrigger = 150,000
  ### additional runs for Ftrgt=0.37: 0.9*Ftrgt & 1.1*Ftrgt
  ### FMSYlower/upper
  ### and original Ftrgt=0.31
  hcr_vals <- rbind(hcr_vals,
    expand.grid(Btrigger = 150000, 
                Ftrgt = c(0.37*0.9, 0.37*1.1, 0.198, 0.46, 0.31))
  )
  ### 63-67
  ### find where risk surpasses 5% at high Btrigger values
  hcr_vals <- rbind(hcr_vals,
    expand.grid(Btrigger = 180000, Ftrgt = c(0.39, 0.4, 0.41, 0.42)),
    expand.grid(Btrigger = 190000, Ftrgt = c(0.41, 0.42, 0.43, 0.44))
  )
  ### 68-75

} else if (HCRoption %in% 2:6) {
  hcr_vals <- expand.grid(
    Btrigger = seq(from = 110000, to = 190000, length.out = 5),
    Ftrgt = c(0.1, 0.2, 0.3, 0.35, 0.37, 0.4, 0.5))
  ### 1-35
}

### implement
if (exists("HCR_comb")) {
  
  ### set Btrigger
  Btrigger <- hcr_vals[HCR_comb, "Btrigger"]
  input$ctrl.mp$ctrl.phcr@args$Btrigger <- Btrigger
  input$ctrl.mp$ctrl.is@args$hcrpars$Btrigger <- Btrigger
  
  ### set Ftrgt
  Ftrgt <- hcr_vals[HCR_comb, "Ftrgt"]
  input$ctrl.mp$ctrl.phcr@args$Ftrgt <- Ftrgt
  input$ctrl.mp$ctrl.is@args$hcrpars$Ftrgt <- Ftrgt
  
  cat(paste0("\nSetting custom Btrigger/Ftrgt values.\n",
             "Using HCR_comb = ", HCR_comb, "\n",
             "Ftrgt = ", Ftrgt, "\n",
             "Btrigger = ", Btrigger, "\n\n"))
  
} else {
  
  cat(paste0("\nUsing default Btrigger/Ftrgt values.\n",
             "Ftrgt = ", input$ctrl.mp$ctrl.phcr@args$Ftrgt, "\n",
             "Btrigger = ", input$ctrl.mp$ctrl.phcr@args$Btrigger, "\n\n"))
  
}

### ------------------------------------------------------------------------ ###
### TAC constraint
input$ctrl.mp$ctrl.is@args$TAC_constraint <- FALSE
### check conditions
### either manually requested or as part of HCR options 4-6 
if (exists("TAC_constraint")) {
  if (isTRUE(as.logical(TAC_constraint))) {
    input$ctrl.mp$ctrl.is@args$TAC_constraint <- TRUE
  }
}
if (HCRoption %in% 4:6) {
    input$ctrl.mp$ctrl.is@args$TAC_constraint <- TRUE
}
### implement
if (isTRUE(input$ctrl.mp$ctrl.is@args$TAC_constraint)) {
    
    input$ctrl.mp$ctrl.is@args$lower <- 80
    input$ctrl.mp$ctrl.is@args$upper <- 125
    input$ctrl.mp$ctrl.is@args$Btrigger_cond <- TRUE
    
    cat(paste0("\nImplementing TAC constraint.\n\n"))
    
} else {
    
    cat(paste0("\nTAC constraint NOT implemented.\n\n"))
    
}

### ------------------------------------------------------------------------ ###
### banking & borrowing
input$ctrl.mp$ctrl.is@args$BB <- FALSE
input$iem <- NULL

### check conditions
### either manually requested or as part of HCR options 4-6 
if (exists("BB")) {
  if (isTRUE(as.logical(TAC_constraint))) 
    
    input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE))
    input$ctrl.mp$ctrl.is@args$BB <- TRUE
    input$ctrl.mp$ctrl.is@args$BB_conditional <- TRUE
    input$ctrl.mp$ctrl.is@args$BB_rho <- c(-0.1, 0.1)

}
if (HCRoption %in% 4:6) {
    
  input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE))
  input$ctrl.mp$ctrl.is@args$BB <- TRUE
  input$ctrl.mp$ctrl.is@args$BB_rho <- c(-0.1, 0.1)
  
  if (HCRoption %in% 4) {
    
    input$ctrl.mp$ctrl.is@args$BB_conditional <- FALSE
    
  } else if (HCRoption %in% 5:6) {
    
    input$ctrl.mp$ctrl.is@args$BB_conditional <- TRUE
    
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
           ctrl.mp = input$ctrl.mp,
           genArgs = input$genArgs,
           tracking = input$tracking)

### save results
path_out <- paste0("output/runs/cod4/", iters, "_", years)
dir.create(path = path_out, recursive = TRUE)
file_out <- paste0("HCR-", input$ctrl.mp$ctrl.hcr@args$option,
                   "_Ftrgt-", input$ctrl.mp$ctrl.phcr@args$Ftrgt,
                   "_Btrigger-", input$ctrl.mp$ctrl.phcr@args$Btrigger,
                   "_TACconstr-", input$ctrl.mp$ctrl.is@args$TAC_constraint,
                   "_BB-", input$ctrl.mp$ctrl.is@args$BB
)

saveRDS(object = res1, paste0(path_out, "/", file_out, ".rds"))

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
if (par_env == 1 & exists("kill")) {
  system("bkill $LSB_JOBID")
}

quit(save = "no")


