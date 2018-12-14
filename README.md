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
-   `a4a_mse_WKNSMSE_funs.R`: functions used in OM.R and for the MSE

and the following subdirectories:

-   `legacy_code/`: Legacy code with original MSE approach and original a4a mse framework. This is no longer supported.

Prerequisites
-------------

It is recommended to use R version 3.5.x or higher and the latest versions of R packages (as of 12/2018). The framework uses FLR and requires the following FLR packages:

-   `FLCore`
-   `FLash`
-   `FLAssess`
-   `ggplotFL`
-   `mse`

They can be installed from www.flr-project.org with

``` r
install.packages(pkgs = c("FLCore", "FLash", "FLAssess", "ggplotFL"), 
                 repos = "http://flr-project.org/R")
### package mse not available here
```

or from GitHub (recommended) with the `devtools` package:

``` r
devtools::install_github(repo = "flr/FLCore")
devtools::install_github(repo = "flr/FLash")
devtools::install_github(repo = "flr/FLAssess")
devtools::install_github(repo = "flr/ggplotFL")
devtools::install_github(repo = "flr/mse")
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

The cod MSE uses the <https://github.com/flr/mse> framework. An example of is at the end of the `OM.R` script.

... to be updated
