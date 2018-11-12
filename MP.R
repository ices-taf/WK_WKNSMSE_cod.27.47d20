### ------------------------------------------------------------------------ ###
### MSE for North Sea cod for WKNSMSE
### based on MSE for ICES WKLIFE VII & VIII
### based on the a4a standard MSE developed at JRC
### by Ernesto Jardim, Iago Mosqueira, Finlay Scott, et al.
### ------------------------------------------------------------------------ ###
### author: Simon Fischer (Cefas), simon.fischer@cefas.co.uk
### ------------------------------------------------------------------------ ###
### created 07/2018
### last modifications:
### 2018-07 Simon Fischer
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### load arguments passed to R ####
### ------------------------------------------------------------------------ ###
### load arguments
args <- commandArgs(TRUE)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  
  ### set default values
  ### number of cores, i.e. processes to spawn
  if (!isTRUE(exists("n_cores"))) {
    stop("n_cores need to be passed to R!")
  } else {
    n_cores <- n_cores - 1 ### slaves, exluding master
  }
  ### parallelization architecture
  if (!isTRUE(exists("cluster_type"))) cluster_type <- 2
  ### split each scenario into n parts?
  if (!isTRUE(exists("n_parts"))) n_parts <- 1
  ### scenarios to be simulated
  if (!isTRUE(exists("scn_start")) | !isTRUE(exists("scn_end"))) {
    scns <- TRUE
  } else {
    scns <- scn_start:scn_end
  }

} else {
  n_parts <- 1 ### no split
  scns <- TRUE ### run all scenarios
  cluster_type <- NULL 
}

### ------------------------------------------------------------------------ ###
### set up parallel computing environment ####
### ------------------------------------------------------------------------ ###
### needs to be load first, otherwise the MPI breaks down ...

### for local in-node/PC parallelization
if (isTRUE(cluster_type == 1)) {

  library(doParallel)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  getDoParWorkers()
  getDoParName()

} else if (isTRUE(cluster_type == 2)) {

  library(doMPI)
  cl <- startMPIcluster(count = n_cores) # one less than requested in the .job file
  registerDoMPI(cl)

}

### ------------------------------------------------------------------------ ###
### packages ####
### ------------------------------------------------------------------------ ###
required_pckgs <- c("FLash", "stockassessment", "FLfse", "foreach")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})

### ------------------------------------------------------------------------ ###
### functions ####
### ------------------------------------------------------------------------ ###

### source the scripts from functions folder
invisible(lapply(list.files(path = "functions/", pattern = "*.R$", 
                            full.names = TRUE), source))

### ------------------------------------------------------------------------ ###
### load scenario definitions ####
### ------------------------------------------------------------------------ ###
source("MP_scenarios.R")

### ------------------------------------------------------------------------ ###
### start simulations ####
### ------------------------------------------------------------------------ ###
### "loops" through scenarios
### nested foreach loop, splits scenarios into smaller junks, if requested
Sys.time()

### "loop" through scenarios and parts
res <- foreach(scn = seq_along(ctrl.mps)[scns], .packages = required_pckgs,
               .export = ls(), .errorhandling = "pass") %:%
         foreach(part = 1:n_parts, .errorhandling = "pass") %dopar% {
  ### load external functions
  invisible(lapply(paste0("functions/", load_files), source))

# part <- 1
# scn <- 2

  
  ### set seed depending on part of scenario
  set.seed(part)
  
  ### -------------------------------------------------------------------- ###
  ### load scenario specifications and objects ####
  ### -------------------------------------------------------------------- ###
  ### load control object
  ctrl.mp <- ctrl.mps[[scn]]
  ### load objects (stock, recruitment, ...)
  load_OMs(ctrl.mp)
  
  ### ---------------------------------------------------------------------- ###
  ### split/subset objects, if required ####
  ### ---------------------------------------------------------------------- ###
  if (n_parts > 1) {

    ### check if splitting feasible
    if (it %% n_parts != 0) stop("'it' cannot be split into 'n_parts'!")
    ### get desired iterations
    it_part <- split(1:it, cut(1:it, n_parts))[[part]]
    it <- it/n_parts
    
    ### subset OM to iterations for current part
    stk <- FLCore::iter(stk, it_part)
    idx <- lapply(idx, FLCore::iter, iter = it_part)
    sr <- FLCore::iter(sr, it_part)
    sr_res <- FLCore::iter(sr_res, it_part)

  }
  
  ### ---------------------------------------------------------------------- ###
  ### create object for tracking ####
  ### ---------------------------------------------------------------------- ###
  tracking <- FLQuant(NA,
    dimnames = list(metric = c("MP.f", "MP.ssb", "convergence", "advice", 
                               "Implementation", "IEM", "FleetDyn",
                               "OM.f", "OM.ssb", "OM.catch"),
                    year = dimnames(catch(stk))$year,
                    iter = 1:it))
  ### start tracking
  tracking["Implementation", ac(iy)] <- catch(stk)[, ac(iy)]
  
  ### ---------------------------------------------------------------------- ###
  ### loop through simulation years ####
  ### ---------------------------------------------------------------------- ###
  for (ay in an(vy[-length(vy)])) {
  #ay = an(vy[1])
    
    gc()
    cat(ay, "> ")
    
    ### -------------------------------------------------------------------- ###
    ### OEM ####
    ### -------------------------------------------------------------------- ###
    ### observations & observation error
    
    ### -------------------------------------------------------------------- ###
    ### oFun - observations
    
    ### create list with parameters for o()
    ctrl.oem <- create_ctrl(ctrl.mp, name = "ctrl.oem", stk = stk, idx = idx,
                            catch_res = if (exists("catch_res")) catch_res,
                            ay = ay, tracking = tracking)
    ### call observation function
    o.out <- do.call("oFun", ctrl.oem)
    ### retrieve observations
    load_elements(o.out) ### loads observed stk0 & idx0, and idx & tracking
    
    ### -------------------------------------------------------------------- ###
    ### MP ####
    ### -------------------------------------------------------------------- ###
    
    ### -------------------------------------------------------------------- ###
    ### fFun: Assessment/Estimator of stock statistics
    
    ### create list with parameters
    ctrl.f <- create_ctrl(ctrl.mp, name = "ctrl.f", stk0 = stk0, idx0 = idx0, 
                          tracking = tracking, ay = ay, parallel = TRUE)
    ### WARNING: parallel=TRUE only for local testing
    ### this is inefficient as only the assessment is parallelized
    ### but everything else runs sequentially
    
    ### run assessment
    f.out <- do.call("fFun", ctrl.f)
    ### retrieve assessment results
    load_elements(f.out) ### loads assessed stk0, tracking (and model fit)
    
    ### -------------------------------------------------------------------- ###
    ### xFun: HCR parametrization
    
    ### create list with parameters
    ctrl.x <- create_ctrl(ctrl.mp, name = "ctrl.x")
    ### get HCR parameters
    x.out <- do.call("xFun", ctrl.x)
    ### load parameters
    load_elements(x.out) ### refpts
    
    ### -------------------------------------------------------------------- ###
    ### hFun(): apply HCR
    
    ### create list with parameters
    ctrl.h <- create_ctrl(ctrl.mp, name = "ctrl.h", stk0 = stk0, ay = ay,
                          tracking = tracking, refpts = refpts, it = it)
    ### apply HCR
    h.out <- do.call("hFun", ctrl.h)
    ### extract (ctrl & tracking)
    load_elements(h.out)
    
    ### -------------------------------------------------------------------- ###
    ### kFun(): Management Implementation (short term forecast)
    
    ### create list with parameters
    ctrl.k <- create_ctrl(ctrl.mp = ctrl.mp, name = "ctrl.k", stk0 = stk0, 
                          ay = ay, tracking = tracking, ctrl = ctrl,
                          refpts = refpts, it = it)
    
    ### implement HCR
    k.out <- do.call("kFun", ctrl.k)
    ### extract (ctrl & tracking)
    load_elements(k.out)
    
    ### -------------------------------------------------------------------- ###
    ### wFun: Technical measures
    
    ### create list with parameters
    ctrl.w <- create_ctrl(ctrl.mp = ctrl.mp, name = "ctrl.w", 
                          ay = ay, tracking = tracking, ctrl = ctrl)
    ### apply measures
    w.out <- do.call("wFun", ctrl.w)
    ### extract (ctrl & tracking)
    load_elements(w.out)
    
    ### -------------------------------------------------------------------- ###
    ### IEM ####
    ### -------------------------------------------------------------------- ###
    
    ### -------------------------------------------------------------------- ###
    ### lFun: implementation error
    
    ### create list with parameters
    ctrl.l <- create_ctrl(ctrl.mp = ctrl.mp, name = "ctrl.l", 
                          ay = ay, tracking = tracking, ctrl = ctrl)
    ### apply measures
    l.out <- do.call("lFun", ctrl.l)
    ### extract 
    load_elements(l.out)
    
    ### -------------------------------------------------------------------- ###
    ### OM ####
    ### -------------------------------------------------------------------- ###
    ### jFun: fleet dynamics/behaviour
    
    ### create list with parameters
    ctrl.j <- create_ctrl(ctrl.mp = ctrl.mp, name = "ctrl.j", 
                          ay = ay, tracking = tracking, ctrl = ctrl)
    ### apply measures
    j.out <- do.call("lFun", ctrl.j)
    ### extract 
    load_elements(j.out)
    
    ### -------------------------------------------------------------------- ###
    ### stock dynamics and OM projections

    ### project
    stk[] <- fwd(stk, ctrl = ctrl, sr = sr, sr.residuals = sr_res,
                 sr.residuals.mult = TRUE, maxF = 5)[]
    
    print(plot(stk))
    
    ### insert results into tracking object
    tracking["OM.f", ac(ay + 1)] <- fbar(stk)[, ac(ay + 1)]
    tracking["OM.ssb", ac(ay + 1)] <- ssb(stk)[, ac(ay + 1)]
    tracking["OM.catch", ac(ay + 1)] <- catch(stk)[, ac(ay + 1)]
    
  } ### end of year loop
  cat("\n")
  
  ### save tracking & scenario definition as attributes in stk
  attr(stk, "tracking") <- tracking
  attr(stk, "ctrl.mp") <- ctrl.mp
  ### save stk to disk
  saveRDS(stk, file = paste0("output/runs/", scn, "_", part, ".rds"))
  
} ### end of scenario loop

Sys.time()

### ------------------------------------------------------------------------ ###
### close down parallel workers ####
### ------------------------------------------------------------------------ ###

if (exists("cluster_type")) {
  if (cluster_type == 1) {
    stopCluster(cl)
  } else if (cluster_type == 2) {
    closeCluster(cl)
    mpi.quit()
  }
}

### quit R, if not already closed by previous shutdown signals
quit("no")

