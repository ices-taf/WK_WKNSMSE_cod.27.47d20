WK\_WKNSMSE\_cod.27.47d20
================

*note!*

This repo contains submodules. To update the submodules to the *commit referenced in the github repository* run:

    git submodule update

And to *update the reference* to the most recent commit in the tools repo, go into that directory and perform a git pull

    cd wk_WKNSMSE_tools
    git pull origin master

and then navigate back to the repository root directory and commit and push the update

    cd ..
    git commit -am"updated submodule refs
    git push

Some useful links: \* <https://git-scm.com/book/en/v2/Git-Tools-Submodules> \* <https://blog.github.com/2016-02-01-working-with-submodules/>

a4a cod4 MSE
============

Introduction
------------

This repository contains a simple Management Strategy Evaulation (MSE) template for North Sea cod. It is based on the Fisheries Library in R ([FLR](http://www.flr-project.org/)) Assessment for All (a4a) standard MSE framework developed during the Workshop on development of MSE algorithms with R/FLR/a4a ([Jardim et al., 2017](https://ec.europa.eu/jrc/en/publication/assessment-all-initiativea4a-workshop-development-mse-algorithms-rflra4a)).

Repository structure
--------------------

The root folder contains the following R scripts:

-   `OM.R`: create operating model(s)
-   `MP_scenarios.R`: define scenarios for MSE
-   `MP.R` script which runs the MSE simulation

and the following subdirectories:

-   `functions/`: contains additional functions sourced in `OM.R` and `MP.R`
-   `input/`: created by running `OM.R`, stores the operating models (stock, indices, stock recruitment model and residuals)
-   `output/`: results from running the MSE , created by `MP.R`

Prerequisites
-------------

It is recommended to use R version 3.5.x or higher and the latest versions of R packages (as of 11/2018). The framework uses FLR and requires the following FLR packages:

-   `FLCore`
-   `FLash`
-   `FLAssess`
-   `ggplotFL`

They can be installed from www.flr-project.org with

``` r
install.packages(pkgs = c("FLCore", "FLash", "FLAssess", "ggplotFL"), 
                 repos = "http://flr-project.org/R")
```

or from GitHub with the `devtools` package:

``` r
devtools::install_github(repo = "flr/FLCore")
devtools::install_github(repo = "flr/FLash")
devtools::install_github(repo = "flr/FLAssess")
devtools::install_github(repo = "flr/ggplotFL")
```

The cod MSE uses SAM as stockassessment method and requires the `stockassessment` R package from GitHub:

``` r
install.packages("TMB") ### required to use stockassessment
devtools::install_github("fishfollower/SAM/stockassessment")
```

In order to use SAM within FLR, the functionality of the inofficial package `FLfse` is used and can be installed with:

``` r
devtools::install_github("shfischer/FLfse/FLfse")
```

This package also contains the latest data for cod and can recreate the assessment as conducted during wGNSSK2018.

Furthermore, the packages `foreach` and `DoParallel` are required:

``` r
install.packages(c("foreach", "DoParallel")) 
```

Running the code
----------------

The code is designed to work on personal computers as well as on high performance computing clusters. The simulation is started by executing `MP.R` and this script can take arguments passed on to it, e.g. from submitting a job to an LSF system on a cluster.

``` r
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
```

Parallelization
---------------

The simulation uses `foreach` which allows parallelization independent of the underlying architecture. On individual computers (or computing nodes) `DoParallel` is used and `DoMPI` on computing clusters which can spread the simulation or multiple computing nodes.

The parallel environment is set up in `MP.R` with:

``` r
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
```

In the `MP.R` script, the following nested `foreach` loop is used:

``` r
### "loop" through scenarios and parts
res <- foreach(scn = seq_along(ctrl.mps)[scns], .packages = required_pckgs,
               .export = ls(), .errorhandling = "pass") %:%
         foreach(part = 1:n_parts, .errorhandling = "pass") %dopar% {
  
  ### MSE code
}
```

With this construct the simulation is split into individual scenarios and the individual scenarios can be split further into several parts. This reduces the overall runtime if many computer cores are available.

Scenarios
---------

The scenarios for the MSE simulation are defined in the script `MP_scenarios.R`. The scenario definition are stored in a list `ctrl.mps`. Each element in this list defines one scenario. For each scenario, there is a element `ctrl.om` which defines which operating model is used and the second element is `ctrl.def` which defines the year range of the simulation. The following elements each correspond to one of the modules of the MSE framwork and select a method as well as some additional arguments passed to this method. Skipping a module will lead to the MSE simulation skipping this module.

Structure of `MP.R`
-------------------

The structure of the `MP.R` script follows the a4a standard MSE.

In the beginning, the environment is set up by loading packages and additional functions. The scenario definitions are loaded by sourcing `MP_scenarios.R`.

The scenarios are executed in a nested `foreach` loop. Within each of these iterations, the operating model is loaded first with the following elements:

-   `stk` is an object of class `FLStock` with the biological stock
-   `idx` is an object of class `FLIndices` with the survey indices
-   `sr` is an object of class `FLSR` storing the stock recruitment model
-   `sr_res` stores the precompiled recruitment residuals in an `FLQuant`

If the scenario is split into several parts, the operating model is subset to the currently executed part by selecting the appropriate iterations.

The object `tracking` is used to track some of the management quantities.

The simulation of the current scenario is then started by going through all the years in a loop. In each year of this loop, the following modules are executed:

-   `oFun`: The observation model. This creates the observations as seen later by the management procedure. This includes the calculation of the survey index/indices, catches, etc. and adds uncertainty to them, if requested.
-   `fFun`: The stock assessment. In here, the stock assessment (SAM) is conducted, based on the data from the observation model.
-   `xFun`: HCR parametrization. This loads the parameters used in the harvest control rule (Btrigger, Ftrgt).
-   `hFun`: HCR. Application of HCR and set the target fishing mortality.
-   `kFun`: Management Implementation. Conducts a short-term forecast within SAM and converts the target fishing mortality into a catch value.
-   `wFun`: Technical measures. Implementation of technical measures, here: TAC constraint.
-   `lFun`: Implementation error. Deviation between advised and implemented catch. Should be precompiled. NOT implemented yet...
-   `jFun`: Fleet dynamics/behaviour. NOT implemented.

For each of these modules, first, a control object is created which contains all the data required to run the module. This is then passed to the method the results are loaded into the current environment.

After going through all modules, the stock is projected forward by one year and the next year is started.
