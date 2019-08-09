library(mse)
source("a4a_mse_WKNSMSE_funs.R")
### get normal input data for base OM
input <- readRDS("input/cod4/1000_20/base_run.rds")
### remove observation model
input$oem <- FLoem(observations = list(stk = FLQuant(0)), 
                   deviances = list(stk = FLQuant(0)))
### fixed F=0 target
input$ctrl.mp <- mpCtrl(list(ctrl.hcr = mseCtrl(method = fixedF.hcr,
                                                args = list(ftrg = 0))))
### use one block, no parallelisation
### (this is ok here, because the only stochastic part -SAM forecast- is not used)
input$genArgs$nblocks <- 1

### run as usual
### (should not take more than a few minutes on one core)
res1 <- mp(om = input$om,
           oem = input$oem,
           iem = input$iem,
           ctrl.mp = input$ctrl.mp,
           genArgs = input$genArgs,
           tracking = input$tracking)

### save
saveRDS(object = res1, file = "output/runs/cod4/1000_20/cod4_F0.rds")

