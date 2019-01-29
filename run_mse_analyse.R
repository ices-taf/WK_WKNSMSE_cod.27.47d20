### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)

library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

### load additional functions
source("a4a_mse_WKNSMSE_funs.R")

### ------------------------------------------------------------------------ ###
### get files ####
### ------------------------------------------------------------------------ ###

path_res <- "output/runs/cod4/1000_20/"
files_res <- data.frame(file = list.files(path_res, pattern = "*.rds"), 
                        stringsAsFactors = FALSE)
files_res <- files_res[files_res$file != "stats.rds",, drop = FALSE]

files_res$Ftrgt <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                     split = "_", fixed = TRUE),
                                            "[[", 2), gsub,
                                     pattern = "Ftrgt-", replacement = ""))
files_res$Btrigger <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                     split = "_", fixed = TRUE),
                                            "[[", 3), gsub,
                                     pattern = "Btrigger-", replacement = ""))
files_res$HCR <- sapply(lapply(strsplit(files_res$file, 
                                        split = "_", fixed = TRUE),
                               "[[", 1), gsub,
                        pattern = "HCR-", replacement = "")
files_res$TACconstr <- as.logical(sapply(lapply(strsplit(files_res$file, 
                                                      split = "_", fixed = TRUE),
                                               "[[", 4), gsub,
                                      pattern = "TACconstr-", replacement = ""))
files_res$BB <- as.logical(
  sapply(lapply(lapply(strsplit(files_res$file, split = "_", fixed = TRUE), 
                       "[[", 5), gsub, pattern = "BB-", replacement = ""), 
         gsub, pattern = ".rds", replacement = ""))

stats <- readRDS(paste0(path_res, "stats.rds"))
stats_new <- merge(stats, files_res, all = TRUE)
stats_new <- stats_new[!stats_new$file %in% stats$file, ]

res_list <- vector(mode = "list", length = nrow(stats_new))
res_list <- foreach(i = seq(nrow(stats_new))) %dopar% {
  readRDS(paste0(path_res, stats_new$file[i]))
}

### ------------------------------------------------------------------------ ###
### calculate summary statistics ####
### ------------------------------------------------------------------------ ###
### calculate for short- (year 1-5), medium- (year 6-10) and 
### long-term (year 11-20)
### risk 1: proportion of stock below Blim, average over iterations and period
### risk 3: maximum of annual proportions
### catch: mean catch in period (mean over years and iterations)
### iav: inter-annual variation of catch, average over years and iterations


### catch
stats_new$catch_long <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(catch(x@stock), start = 2029))
}
stats_new$catch_short <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(catch(x@stock), start = 2019, end = 2023))
}
stats_new$catch_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(catch(x@stock), start = 2024, end = 2028))
}
### risks
stats_new$risk1_full <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2019) < 107000)
}
stats_new$risk1_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2029) < 107000)
}
stats_new$risk1_short <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2019, end = 2023) < 107000)
}
stats_new$risk1_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2024, end = 2028) < 107000)
}
stats_new$risk3_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2029) < 107000))
}
stats_new$risk3_short <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2019, end = 2023) < 107000))
}
stats_new$risk3_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2024, end = 2028) < 107000))
}
### inter-annual variation of catch
stats_new$iav_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c", .export = "iav") %dopar% {
  iav(object = catch(window(stock(x), start = 2028)), summary_per_iter = mean,
      summary = mean)
}
stats_new$iav_short <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c", .export = "iav") %dopar% {
  iav(object = catch(window(stock(x), start = 2018, end = 2023)), 
      summary_per_iter = mean, summary = mean)
}
stats_new$iav_medium <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c", .export = "iav") %dopar% {
  iav(object = catch(window(stock(x), start = 2023, end = 2028)), 
      summary_per_iter = mean, summary = mean)
}

stats <- rbind(stats, stats_new)
saveRDS(object = stats, file = paste0(path_res, "stats.rds"))

### ------------------------------------------------------------------------ ###
### plot ####
### ------------------------------------------------------------------------ ###


### plot function
grid <- function(dat, HCR = "A",
                 time = c("long", "short", "medium"),
                 add_risk1 = FALSE) {
  
  ### catch
  dat$catch <- dat[, paste0("catch_", time)]
  dat$risk <- dat[, paste0("risk3_", time)]
  dat$iav <- dat[, paste0("iav_", time)]
  dat$risk1 <- dat[, paste0("risk1_", time)]
  p1 <- ggplot(data = dat, 
                aes(x = Btrigger, y = Ftrgt, fill = catch)) +
    geom_raster() +
    scale_fill_gradient(paste0(time, "-term\ncatch"), low = "red", 
                        high = "green") +
    geom_text(aes(label = round(catch), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw()
  ### risk
  p2 <- ggplot(data = dat, 
                aes(x = Btrigger, y = Ftrgt, fill = risk)) +
    geom_raster(alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 3"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw()
  ### iav
  p3 <- ggplot(data = dat, 
                aes(x = Btrigger, y = Ftrgt, fill = iav)) +
    geom_raster() +
    geom_text(aes(label = round(iav, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05",
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    scale_fill_gradient(paste0(time,"-term\ninter-annual\nvariability of catch"),
                        low = "green", high = "red") +
    theme_bw()
  ### risk1
  p4 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt, fill = risk1)) +
    geom_raster(alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 1"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk1, 3), colour = risk1 <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw()
  
  if (isTRUE(add_risk1)) {
    plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv")
  } else {
    plot_grid(p1, p2, p3, nrow = 2, ncol = 2, align = "hv")
  }
  
  
}

### A
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "long", add_risk1 = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/grid_A_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/grid_A_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "A", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/grid_A_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### B
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "long", add_risk1 = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/grid_B_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/grid_B_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "B", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/grid_B_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### C
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "long", add_risk1 = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/grid_C_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/grid_C_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE) %>%
       filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)), 
     HCR = "C", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/grid_C_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")




