### ------------------------------------------------------------------------ ###
### scenario defintion ####
### ------------------------------------------------------------------------ ###
### this script creates a list called ctrl.mps
### each element in the list defines a particular scenario

### a dummy template for NS cod
ctrl.mp1 <- list(
  ### definition of operating model
  ctrl.om = list(name = "cod4_non_mult", stk = "stk", sr = "sr", sr_res = "sr_res",
                 idx = "idx", catch_res = "catch_res"),
  ### simulation specifications
  ctrl.def = list(
    fy = 2048, # final year
    y0 = 1963, # initial data year
    dy = 2018, # final data year
    iy = 2018, # initial year of projection (also intermediate)
    ny = 32, # number of years to project from intial year
    vy = ac(2018:2048) # vector of years to be projected
  ),
  ctrl.oem = list(method = "o_WKNSMSE", idx_timing = c(0, -1),
                  catch_timing = -1),
  ctrl.f = list(method = "SAM_wrapper", conf = cod4_conf_sam,
                forecast = TRUE, ### do a forecast, .e.g. to get SSB 
                fwd_trgt = "fsq", ### what to target in forecast
                fwd_yrs = 1, ### number of years to add
                fwd_yrs_average = -3:0, ### years used for averages
                fwd_yrs_rec_start = 1998, ### recruitment 
                fwd_yrs_sel = -3:-1, ### selectivity
                fwd_yrs_lf_remove = -2:-1, ### remove years from landing fraction
                fwd_splitLD = TRUE),
  ctrl.x = list(method = "x_ICES", Btrigger = 150000, Ftrgt = 0.31),
  ctrl.h = list(method = "h_ICES", option = "A"),
  ctrl.k = list(method = "k_SAM_forecast",
                fwd_trgt = c("fsq", "hcr"), ### what to target in forecast
                fwd_yrs = 2, ### number of years to add
                fwd_yrs_average = -3:0, ### years used for averages
                fwd_yrs_rec_start = 1998, ### first recruitment sample year
                fwd_yrs_sel = -3:-1, ### selectivity years
                fwd_yrs_lf_remove = -2:-1, ### remove years from landing fraction
                fwd_splitLD = TRUE),
  ctrl.w = list(method = "TAC_constraint", upper = Inf, lower = -Inf)
)

### without catch multiplier
ctrl.mp2 <- ctrl.mp1
ctrl.mp2$ctrl.om$name <- "cod4"
ctrl.mp2$ctrl.f$conf <- cod4_conf_sam[!names(cod4_conf_sam) %in% 
                                        c("noScaledYears", "keyScaledYears",
                                          "keyParScaledYA")]
ctrl.mp2$ctrl.f$newtonsteps <- 0

ctrl.mps <- list(
  no_mult = ctrl.mp1,
  cod4 = ctrl.mp2
)
