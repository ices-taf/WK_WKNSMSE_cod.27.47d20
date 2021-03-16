WK\_WKNSMSE\_cod.27.47d20
================

# a4a cod4 MSE

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
  - `OM_alt1.R`, `OM_alt2.R`and `OM_alt1.R` are scripts for creating
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
sessionInfo()
R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.7 (Final)
```

The framework uses FLR and requires the following FLR packages:

  - `FLCore` 2.6.11.9001
  - `FLash` 2.5.11
  - `FLBRP` 2.5.3
  - `FLAssess` 2.6.3
  - `ggplotFL` 2.6.5
  - `mse` 0.9.1

The FLR package versions as used for the simulation can be installed
with `devtools`:

``` r
devtools::install_github(repo = "flr/FLCore", ref = "d55bc6570c0134c6bea6c3fc44be20378691e042")
devtools::install_github(repo = "flr/FLash", ref = "7c47560cf57627068259404bb553f2b644682726")
devtools::install_github(repo = "flr/FLBRP", ref = "5644cfccefb0ec3965b1d028090bbf75b1e59da2")
devtools::install_github(repo = "flr/FLAssess", ref = "f1e5acb98c106bcdfdc81034f1583f76bb485514")
devtools::install_github(repo = "flr/ggplotFL", ref = "e9e0d74e872815c1df3f172522da35ade5c70638")
devtools::install_github(repo = "flr/mse", ref = "e39ddd75cdb2bb693601e31428404d48ea810308")
```

The cod MSE uses SAM as stockassessment method and requires the
[`stockassessment`](https://github.com/fishfollower/SAM/) R package from
GitHub:

  - `stockassessment` 0.8.1

Install this version with

``` r
install.packages("TMB") ### required to use stockassessment
devtools::install_github("fishfollower/SAM/stockassessment", ref = "362d4b661742ed418ce1af890ced05b374f6d6ed")
```

In order to use SAM within FLR, the functionality of the inofficial
package `FLfse` is used:

  - `FLfse` 0.0.0.9003

<!-- end list -->

``` r
devtools::install_github("shfischer/FLfse/FLfse", ref = "c561f5bf28cbad0f711ef53a49bde7e9868dc257")
```

This package also contains the latest data for cod and can recreate the
assessment as conducted during wGNSSK 2018.

Furthermore, some more R packages available from CRAN are required:

``` r
install.packages(c("foreach", "DoParallel", "dplyr", "tidyr", "data.table")) 
```
