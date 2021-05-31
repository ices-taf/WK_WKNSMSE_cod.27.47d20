WK\_WKNSMSE\_cod.27.47d20
================

# WKNSMSE cod4 MSE - updated for the FLasher branch of flr/mse

## Introduction

This repository contains the Management Strategy Evaluation (MSE) for
North Sea cod (Cod (*Gadus morhua*) in Subarea 4, Division 7.d, and
Subdivision 20 (North Sea, eastern English Channel, Skagerrak),
cod.27.47d20) conducted in 2018/2019 for ICES WKNSMSE
(<https://doi.org/10.17895/ices.pub.5090>).

The simulation is based on the Fisheries Library in R
([FLR](http://www.flr-project.org/)) and makes use of the Assessment for
All (a4a) standard MSE framework ([`FLR/mse`](github.com/FLR/mse))
developed during the Workshop on development of MSE algorithms with
R/FLR/a4a ([Jardim et
al., 2017](https://ec.europa.eu/jrc/en/publication/assessment-all-initiativea4a-workshop-development-mse-algorithms-rflra4a)).

## Disclaimer

The content of this repository is licensed under the GNU General Public
License v3.0, see [`LICENSE`](LICENSE) for more information. The content
in this repository is provided as is, without any support by the
contributors listed in this GitHub repository.

## Repository structure

The root folder contains the following R scripts:

  - `OM.R`: This script creates the baseline operating model (`OM.md` is
    the corresponding R Markdown file),
  - `OM_alt1.R`, `OM_alt2.R`and `OM_alt3.R` are scripts for creating
    alternative operating models,
  - `a4a_mse_WKNSMSE_funs.R` contains functions and methods used for the
    creation of the operating models and for running the MSE,
  - `run_mse.R` is an R script for running MSE scenarios and is called
    from a job submission script
  - `run_mse.bsub` and `run_mse.qsub` are job submission scripts which
    are used on a high performance computing cluster and call
    `run_mse.R`
  - `run_mse_analyse.R` is for analysing the MSE results

and the following subdirectories:

  - `legacy_code/`: Legacy code with original MSE approach and original
    a4a mse framework. This is no longer supported.

## R, R packages and version info

The MSE simulation was run on a high performance computing cluster:

``` r
R version 4.1.0 (2021-05-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)
```

The framework uses FLR and requires the following FLR packages:

  - `FLCore` 2.6.15
  - `FLash` 2.5.11
  - `FLBRP` 2.5.3
  - `FLAssess` 2.6.3
  - `ggplotFL` 2.6.5
  - `mse` 2.0.3

The FLR package versions as used for the simulation can be installed
with `devtools`:

``` r
devtools::install_github(repo = "flr/FLCore")
devtools::install_github(repo = "shfischer/FLasher", ref = "64a539786eeca08fb273302be3a920dd176dc158")
devtools::install_github(repo = "flr/FLBRP")
devtools::install_github(repo = "flr/FLAssess")
devtools::install_github(repo = "flr/ggplotFL")
devtools::install_github(repo = "shfischer/mse", ref = "a60f7dfc0637db8da17e829ec880f33a789914bc")
```

The cod MSE uses SAM as stockassessment method and requires the
[`stockassessment`](https://github.com/fishfollower/SAM/) R package from
GitHub:

  - `stockassessment` 0.11.0

Install this version with

``` r
install.packages("TMB") ### required to use stockassessment
devtools::install_github("fishfollower/SAM/stockassessment")
```

In order to use SAM within FLR, the functionality of the inofficial
package `FLfse` is used:

  - `FLfse` 0.0.0.9008

<!-- end list -->

``` r
devtools::install_github("shfischer/FLfse/FLfse")
```

This package also contains the latest data for cod and can recreate the
assessment as conducted during wGNSSK 2018.

Furthermore, some more R packages available from CRAN are required:

``` r
install.packages(c("foreach", "DoParallel", "dplyr", "tidyr", "data.table")) 
```
