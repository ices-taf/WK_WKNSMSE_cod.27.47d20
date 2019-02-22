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
files_res <- files_res[!grepl(x = files_res$file, pattern = "*stats*"),, 
                       drop = FALSE]
files_res <- files_res[!grepl(x = files_res$file, pattern = "F0.rds"),, 
                       drop = FALSE]

# files_rename <- files_res$file
# files_rename <- files_rename[grep(x = files_rename, pattern = "^HCR*")]
# file.rename(from = paste0("output/runs/cod4/1000_20/", files_rename),
#             to = paste0("output/runs/cod4/1000_20/cod4_", files_rename))
files_res$OM <- sapply(strsplit(x = files_res$file, split = "\\_HCR"), "[[", 1)
files_res$Ftrgt <- as.numeric(
  gsub(x = regmatches(x = files_res$file, 
                      m = regexpr(text = files_res$file, 
                                  pattern = "Ftrgt-0.[0-9]{1,}")),
       pattern = "Ftrgt-", replacement = ""))
files_res$Btrigger <- as.numeric(
  gsub(x = regmatches(x = files_res$file, 
                      m = regexpr(text = files_res$file,
            pattern = "Btrigger-[0-9]{1,}[e+]{0,}[0-9]{0,}")),
       pattern = "Btrigger-", replacement = ""))
files_res$HCR <- gsub(x = regmatches(x = files_res$file, 
                      m = regexpr(text = files_res$file, 
                                  pattern = "HCR-[A-Z]{1,}")),
                      pattern = "HCR-", replacement = "")
files_res$TACconstr <- as.logical(
  gsub(x = regmatches(x = files_res$file, 
                      m = regexpr(text = files_res$file, 
                                  pattern = "TACconstr-TRUE|TACconstr-FALSE")),
       pattern = "TACconstr-", replacement = ""))
files_res$BB <-  as.logical(
  gsub(x = regmatches(x = files_res$file,
                      m = regexpr(text = files_res$file, 
                                  pattern = "BB-TRUE|BB-FALSE")),
                  pattern = "BB-", replacement = ""))


stats <- readRDS(paste0(path_res, "stats.rds"))
stats_new <- merge(stats, files_res, all = TRUE)
### set Blim depending on OM
stats$Blim <- sapply(stats$OM, function(x) {
  switch(x, "cod4" = 107000,"cod4_alt1" = 107000, "cod4_alt2" = 108000,
         "cod4_alt3" = 107000)})
### keep only new files
stats_new <- stats_new[!stats_new$file %in% stats$file, ]

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
### catch median
stats_new$catch_median_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2029))
}
stats_new$catch_median_short <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2019, end = 2023))
}
stats_new$catch_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                  .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2024, end = 2028))
}
### risks
stats_new$risk1_full <- foreach(x = res_list, 
                                Blim = stats_new$Blim[seq(nrow(stats_new))],
                                .packages = "FLCore",
                                .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2019) < Blim)
}
stats_new$risk1_long <- foreach(x = res_list, 
                                Blim = stats_new$Blim[seq(nrow(stats_new))],
                                .packages = "FLCore",
                                .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2029) < Blim)
}
stats_new$risk1_short <- foreach(x = res_list, 
                                 Blim = stats_new$Blim[seq(nrow(stats_new))],
                                 .packages = "FLCore",
                                 .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2019, end = 2023) < Blim)
}
stats_new$risk1_medium <- foreach(x = res_list, 
                                  Blim = stats_new$Blim[seq(nrow(stats_new))],
                                  .packages = "FLCore",
                                  .combine = "c") %dopar% {
  mean(window(ssb(x@stock), start = 2024, end = 2028) < Blim)
}
stats_new$risk3_long <- foreach(x = res_list, 
                                Blim = stats_new$Blim[seq(nrow(stats_new))],
                                .packages = "FLCore",
                                .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2029) < Blim))
}
stats_new$risk3_short <- foreach(x = res_list, 
                                 Blim = stats_new$Blim[seq(nrow(stats_new))],
                                 .packages = "FLCore",
                                 .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2019, end = 2023) < Blim))
}
stats_new$risk3_medium <- foreach(x = res_list, 
                                  Blim = stats_new$Blim[seq(nrow(stats_new))],
                                  .packages = "FLCore",
                                  .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2024, end = 2028) < Blim))
}
### inter-annual variation of catch
stats_new$iav_long <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  iav(object = catch(window(stock(x), start = 2028)), summary_all = median)
}
stats_new$iav_short <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  iav(object = catch(window(stock(x), start = 2018, end = 2023)), 
      summary_all = median)
}
stats_new$iav_medium <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  iav(object = catch(window(stock(x), start = 2023, end = 2028)), 
      summary_all = mean)
}
### SSB
stats_new$ssb_median_long <- foreach(x = res_list, .packages = "FLCore",
                                       .combine = "c") %dopar% {
  median(window(ssb(x@stock), start = 2029))
}
stats_new$ssb_median_short <- foreach(x = res_list, .packages = "FLCore",
                                        .combine = "c") %dopar% {
  median(window(ssb(x@stock), start = 2019, end = 2023))
}
stats_new$ssb_median_medium <- foreach(x = res_list, .packages = "FLCore",
                                         .combine = "c") %dopar% {
  median(window(ssb(x@stock), start = 2024, end = 2028))
}
### time to recovery
MSYBtrigger <- 150000
stats_new$recovery_proportion <- foreach(x = res_list, .packages = "FLCore",
                                   .combine = "c") %dopar% {
  mean(apply(window(ssb(x@stock), start = 2019) >= MSYBtrigger, 6, max))
}
stats_new$recovery_time <- foreach(x = res_list, .packages = "FLCore",
                                   .combine = "c") %dopar% {
  median(apply(window(ssb(x@stock), start = 2019)@.Data >= MSYBtrigger, 6, 
               function(x) {
    if (any(x)) {
      which(x)[1]
    } else {
      Inf
    }
  }))
}
### SAM convergence
stats_new$conv_failed <- foreach(x = res_list, .packages = "FLCore",
                                .combine = "c") %dopar% {
  sum(x@tracking["conv.est", ac(2018:2037)] != 0)
}
all(stats_new$conv_failed == 0)
### check F maxed (2)
stats_new$F_maxed <- foreach(x = res_list, .packages = "FLCore",
                                 .combine = "c") %dopar% {
  sum(window(fbar(stock(x)), start = 2019) >= 2)
}
# pos <- which(stats_new$F_maxed != 0)
# stats_new[pos, ]
# plot(stock(res_list[[pos]]))
# pos_iter <- which(apply(fbar(stock(res_list[[pos]])), c(1, 6), max) >= 2)
# plot(stock(res_list[[pos]])[,,,,, pos_iter])
### F maxed in THREE scenarios, 
### in each of them in ONE iteration and ONE time only
### HCR B, Btrigger = 110000, Ftrgt = 0.5, iteration 29
### HCR AD, Btrigger = 170000, Ftrgt = 0.5, iteration 442
### HCR AD, Btrigger = 190000, Ftrgt = 0.5, iteration 442
### HCR BE, Btrigger = 190000, Ftrgt = 0.5, iteration 442
### HCR CE, Btrigger = 190000, Ftrgt = 0.5, iteration 442
### OM_alt2 HCR A, Btrigger = 190000, Ftrgt = 0.4, iteration 756

stats <- rbind(stats, stats_new)
stats <- stats[order(stats$file), ]
saveRDS(object = stats, file = paste0(path_res, "stats.rds"))
write.csv(x = stats, file = paste0("output/stats.csv"), row.names = FALSE)

### ------------------------------------------------------------------------ ###
### plot ####
### ------------------------------------------------------------------------ ###


### plot function
grid <- function(dat, HCR = "A",
                 time = c("long", "short", "medium"),
                 add_risk1 = FALSE, highlight_max = FALSE) {
  
  ### catch
  dat$catch <- dat[, paste0("catch_median_", time)]
  dat$risk <- dat[, paste0("risk3_", time)]
  dat$iav <- dat[, paste0("iav_", time)]
  dat$ssb <- dat[, paste0("ssb_median_", time)]
  dat$risk1 <- dat[, paste0("risk1_", time)]
  dat$Btrigger <- dat$Btrigger / 1000
  ### find yield maximum
  dat_max <- dat %>% 
    filter(risk <= 0.05) %>%
    filter(catch == max(catch)) %>%
    select(Ftrgt, Btrigger)
  # p1 <- ggplot(data = dat, 
  #               aes(x = Btrigger, y = Ftrgt, fill = catch)) +
  #   geom_raster() +
  #   scale_fill_gradient(paste0(time, "-term\ncatch (median)"), low = "red", 
  #                       high = "green") +
  #   geom_text(aes(label = round(catch), colour = risk <= 0.05),
  #             size = 2) +
  #   scale_colour_manual("risk <= 0.05", 
  #                       values = c("FALSE" = "red", "TRUE" = "black")) +
  #   theme_bw()
  p1 <- ggplot() +
    geom_raster(data = dat %>% 
                  filter(risk <= 0.05) %>%
                  filter(catch >= 0.95 * max(catch)),
                aes(x = Btrigger, y = Ftrgt, fill = catch)) +
    scale_fill_gradient(paste0("yield maximum\narea [t]"), low = "red",
                        high = "green") +
    geom_text(data = dat, 
              aes(x = Btrigger, y = Ftrgt, 
                  label = round(catch), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, "-term catch [t]")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### risk
  p2 <- ggplot(data = dat, 
                aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = risk), alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 3"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0(time, "-term risk 3")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### iav
  p3 <- ggplot(data = dat, 
                aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = iav)) +
    geom_text(aes(label = round(iav, 3), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05",
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    scale_fill_gradient(paste0("median\n", time,
                               "-term\ninter-annual\ncatch variability"),
                        low = "green", high = "red") +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, 
                        "-term inter-annual\ catch variability")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### SSB
  p4 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = ssb), alpha = 0.5) +
    geom_text(aes(label = round(ssb), colour = risk <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05",
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    scale_fill_gradient(paste0("median\n", time,
                               "-term\nSSB [t]"),
                        low = "red", high = "green") +
    theme_bw() +
    facet_wrap(~ paste0("median ", time, 
                        "-term SSB [t]")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### risk1
  p5 <- ggplot(data = dat, 
               aes(x = Btrigger, y = Ftrgt)) +
    geom_raster(aes(fill = risk1), alpha = 0.75) +
    scale_fill_gradient(paste0(time, "-term\nrisk 1"), 
                        low = "green", high = "red") +
    geom_text(aes(label = round(risk1, 3), colour = risk1 <= 0.05),
              size = 2) +
    scale_colour_manual("risk <= 0.05", 
                        values = c("FALSE" = "red", "TRUE" = "black")) +
    theme_bw() +
    facet_wrap(~ paste0(time, "-term risk 1")) +
    scale_x_continuous(breaks = c(seq(from = 110, to = 190, by = 10))) +
    labs(x = expression(B[trigger]~"[1000t]"),
         y = expression(F[trgt]))
  ### highlight maximum
  if (isTRUE(highlight_max)) {
    p_add <- geom_tile(data = dat_max, aes(x = Btrigger, y = Ftrgt),
                width = 10, height = 0.01, linetype = "solid",
                alpha = 0, colour = "black", size = 0.3)
    p1 <- p1 + p_add
    p2 <- p2 + p_add
    p3 <- p3 + p_add
    p4 <- p4 + p_add
    p5 <- p5 + p_add
  }
  
  if (isTRUE(add_risk1)) {
    plot_grid(p1, p2, p3, p5, nrow = 2, ncol = 2, align = "hv")
  } else {
    plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv")
  }
  
  
}

### A
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "A", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_A_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "A", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_A_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "A", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_A_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### B
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "B", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_B_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "B", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_B_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "B", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_B_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### C
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "C", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_C_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "C", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_C_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == FALSE & TACconstr == FALSE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "C", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_C_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### AD
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "A", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_AD_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "A", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_AD_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "A" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "A", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_AD_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### BE
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "B", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_BE_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "B", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_BE_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "B" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "B", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_BE_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### CE
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "C", time = "long", add_risk1 = FALSE, highlight_max = TRUE)
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_CE_long.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "C", time = "medium")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_CE_medium.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
grid(dat = stats %>%
       filter(HCR == "C" & BB == TRUE & TACconstr == TRUE &
                Ftrgt %in% round(seq(0, 1, 0.01), 2) & OM == "cod4"), 
     HCR = "C", time = "short")
ggsave(filename = "output/runs/cod4/1000_20/plots/grid/grid_CE_short.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")



### ------------------------------------------------------------------------ ###
### plot some runs ####
### ------------------------------------------------------------------------ ###
### function for plotting
plot_stk <- function(stats, OM_ = "cod4", HCR_ = "A", Ftrgt_ = 0.31,
                     Btrigger_ = 150000, TACconstr_ = FALSE, BB_ = FALSE,
                     probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                     path_out = "output/runs/cod4/1000_20/",
                     path_in = "input/cod4/1000_20/", file_in = "base_run.rds",
                     path_res = "output/runs/cod4/1000_20/plots/stock_plots/",
                     save_plot = TRUE,
                     get_history = TRUE, overwrite_catch_history = FALSE,
                     yr_start = 2018.5,
                     Blim = 107000, MSYBtrigger = 150000,
                     Flim = 0.54, Fmsy = 0.31,
                     plot_iter = FALSE,
                     width = 30, height = 20) {
  
  ### get scenario
  stats_i <- stats %>% filter(OM == OM_ & HCR == HCR_ & Ftrgt == Ftrgt_ & 
                                Btrigger == Btrigger_, TACconstr == TACconstr_ &
                                BB == BB_)
  ### load MSE results
  res_i <- readRDS(paste0(path_out, stats_i$file))
  
  ### load stock history
  if (isTRUE(get_history)) {
    hist_i <- readRDS(paste0(path_in, file_in))
    stk_i <- hist_i$om@stock
    stk_i[, dimnames(res_i@stock)$year] <- res_i@stock
    ### overwrite catch history
    if (isTRUE(overwrite_catch_history)) {
      catch_hist <- readRDS(paste0(path_in, "catch_n.rds"))
      catch.n(stk_i)[, dimnames(catch_hist)$year] <- catch_hist
      catch(stk_i) <- computeCatch(stk_i)
    }
  }
  
  ### plot
  if (!isTRUE(plot_iter)) {
    ### plot percentiles
    p <- plot(stk_i, probs = probs) +
      xlab("year") + geom_vline(xintercept = yr_start) +
        geom_hline(data = data.frame(qname = "SSB", data = Blim),
                   aes(yintercept = data), linetype = "dashed") +
        geom_hline(data = data.frame(qname = "SSB", data = MSYBtrigger),
                   aes(yintercept = data), linetype = "solid") +
        geom_hline(data = data.frame(qname = "F", data = Flim),
                   aes(yintercept = data), linetype = "dashed") +
        geom_hline(data = data.frame(qname = "F", data = Fmsy),
                   aes(yintercept = data), linetype = "solid") +
        theme_bw() + ylim(0, NA) + 
      facet_wrap(~ qname, ncol = 1, strip.position = "right", scales = "free_y",
                 labeller = as_labeller(c(
                   "Rec" = "Rec [1000]",
                   "SSB" = "SSB [t]",
                   "Catch" = "Catch [t]",
                   "F" = paste0("F (ages ", paste(range(stk_i)["minfbar"], 
                                      range(stk_i)["maxfbar"], sep = "-"),
                                ")"))))
  } else {
    ### plot iterations
    stk_df <- as.data.frame(FLQuants(`Rec [1000]` = rec(stk_i),
                                      `SSB [t]` = ssb(stk_i),
                                      `Catch [t]` = catch(stk_i),
                                      `F` = fbar(stk_i)))
    p <- ggplot(data = stk_df, 
           aes(x = year, y = data, group = iter)) +
      geom_line(alpha = 0.025) +
      facet_wrap(~ qname, ncol = 1, strip.position = "right",
                 scale = "free_y") +
      theme_bw() +
      ylim(c(0, NA)) + labs(y = "") + 
      geom_vline(xintercept = yr_start) +
      geom_hline(data = data.frame(qname = "SSB [t]", data = Blim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "SSB [t]", data = MSYBtrigger),
                 aes(yintercept = data), linetype = "solid") +
      geom_hline(data = data.frame(qname = "F", data = Flim),
                 aes(yintercept = data), linetype = "dashed") +
      geom_hline(data = data.frame(qname = "F", data = Fmsy),
                 aes(yintercept = data), linetype = "solid")
    
  }
  print(p)
  ### save plot
  if (isTRUE(save_plot)) {
    ggsave(filename = paste0(path_res, 
                             gsub(x = stats_i$file, pattern = ".rds", 
                                  replacement = ""),
                             ".png"),
           width = width, height = height, units = "cm", dpi = 300, 
           type = "cairo")
    
  }

}

### base OM
### current HCR
plot_stk(stats = stats, OM_ = "cod4", HCR_ = "A", Ftrgt_ = 0.31, 
         Btrigger_ = 150000, TACconstr_ = FALSE, BB_ = FALSE, 
         overwrite_catch_history = TRUE)
### A, B, C, AD, BE, CE
combs <- data.frame(name = c("F0", "A*", "A", "B", "C", "AD", "BE", "CE"),
                    HCR = c("F0", "A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(0, 150000, 170000, 160000, 170000, 190000, 130000,
                                 140000),
                    Ftrgt = c(0, 0.31, 0.38, 0.38, 0.38, 0.40, 0.36, 0.36))
combs <- rbind(cbind(combs, OM = "cod4"),
               cbind(combs, OM = "cod4_alt1"),
               cbind(combs, OM = "cod4_alt2"),
               cbind(combs, OM = "cod4_alt3"))
combs <- combs[-c(10, 18, 26), ]
lapply(X = split(combs, seq(nrow(combs))), FUN = function(i) {
  plot_stk(stats = stats, OM_ = i$OM, HCR_ = i$HCR, Ftrgt_ = i$Ftrgt, 
           Btrigger_ = i$Btrigger, TACconstr_ = i$TACconstr, BB_ = i$BB, 
           overwrite_catch_history = TRUE,
           path_in = paste0("input/", i$OM, "/1000_20/"))
})


### ------------------------------------------------------------------------ ###
### F=0 plot ####
### ------------------------------------------------------------------------ ###

stkF0 <- readRDS(file = "input/cod4/10000_100/data_F0.RData")$om@stock
stkF0_res <- readRDS("output/runs/cod4/F0_10000_100.rds")@stock

catch_n <- readRDS("input/cod4/10000_100/catch_n.rds")
catch.n(stkF0)[dimnames(catch_n)$age, dimnames(catch_n)$year] <- 
  catch_n
catch(stkF0) <- computeCatch(stkF0)

stkF0[, ac(2018:2118)] <- stkF0_res[, ac(2018:2118)]
plot(stkF0, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) +
  xlab("year") + geom_vline(xintercept = 2018.5) +
  geom_hline(data = data.frame(qname = "SSB", data = 107000),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "SSB", data = 150000),
             aes(yintercept = data), linetype = "solid") +
  geom_hline(data = data.frame(qname = "F", data = 0.54),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "F", data = 0.31),
             aes(yintercept = data), linetype = "solid") +
  theme_bw() #+
  # geom_blank(data = as.data.frame(FLQuants(`Rec` = rec(stkF0),
  #                                          `SSB` = ssb(stkF0),
  #                                          `Catch` = catch(stkF0),
  #                                          `F` = fbar(stkF0))), 
  #            aes(x = year, y = data, group = iter))
ggsave(filename = paste0("output/runs/cod4/1000_20/", 
                         "stk_F0_10000iters.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
### iterations
stkF0_df <- as.data.frame(FLQuants(`Rec [1000]` = rec(stkF0),
                                  `SSB [t]` = ssb(stkF0),
                                  `Catch [t]` = catch(stkF0),
                                  `F` = fbar(stkF0)))
ggplot(data = stkF0_df[stkF0_df$iter %in% 1:1000, ], 
       aes(x = year, y = data, group = iter)) +
  geom_line(alpha = 0.025) +
  facet_wrap(~ qname, ncol = 1, strip.position = "right",
             scale = "free_y") +
  theme_bw() +
  ylim(c(0, NA)) + labs(y = "") + 
  geom_vline(xintercept = 2018.5) +
  geom_hline(data = data.frame(qname = "SSB [t]", data = 107000),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "SSB [t]", data = 150000),
             aes(yintercept = data), linetype = "solid") +
  geom_hline(data = data.frame(qname = "F", data = 0.54),
             aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = data.frame(qname = "F", data = 0.31),
             aes(yintercept = data), linetype = "solid")
ggsave(filename = paste0("output/runs/cod4/1000_20/", 
                         "stk_F0_10000iters_iters.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### extrapolate 5% risk line ####
### ------------------------------------------------------------------------ ###

# df_risks <- data.frame(Btrigger = c(seq(from = 110000, to = 190000, 
#                                         by = 10000)),
#                        Ftrgt = c(0.355, 0.355, 0.365, 0.375, 0.375, 0.385,
#                                  0.395, 0.405, 0.425))
# 
# lm_risks <- lm(formula = Ftrgt ~ Btrigger, data = tail(df_risks, 5))
# plot(Ftrgt ~ Btrigger, data = df_risks, 
#      xlim = c(110000, 300000), ylim = c(0.3, 0.6))
# abline(lm_risks)


### ------------------------------------------------------------------------ ###
### summary plots: compare HCR options ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations
combs <- data.frame(name = c("F0", "A*", "A", "B", "C", "AD", "BE", "CE"),
                    OM = c("cod4"),
                    HCR = c("F0", "A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(0, 150000, 170000, 160000, 170000, 190000, 130000,
                                 140000),
                    Ftrgt = c(0, 0.31, 0.38, 0.38, 0.38, 0.40, 0.36, 0.36),
                    scenario = 0)
combs <- merge(combs, stats, all.x = TRUE)
combs2 <- gather(data = combs, key = "key", value = "value",
                catch_median_long, risk3_long, iav_long,
                ssb_median_long, recovery_proportion, recovery_time)
combs2$name <- factor(combs2$name, levels = c("F0", "A*", "A", "B", "C", "AD", "BE", "CE"))

ggplot(data = combs2, 
       mapping = aes(x = name, y = value, group = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ key, scales = "free_y") +
  theme_bw()

### load entire distribution for stats
stats_full <- function(data) {
  combs_full <- foreach(i = split(data, seq(nrow(data))), 
                      .packages = "FLCore", .combine = rbind) %dopar% {
                        
    stk_i <- readRDS(paste0("output/runs/cod4/1000_20/", i$file))
    MSYBtrigger <- 150000
    Blim <- ifelse(!i$OM == "cod4_alt2", 107000, 110000)
    res <- rbind(
    data.frame(name = i$name, scenario = i$scenario,
               key = "catch_long",
               value = c(window(catch(stk_i@stock), start = 2029))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "catch_medium",
               value = c(window(catch(stk_i@stock), start = 2019, end = 2023))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "catch_short",
               value = c(window(catch(stk_i@stock), start = 2024, end = 2028))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk1_long",
               value = mean(window(ssb(stk_i@stock), start = 2029) < Blim)),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk1_medium",
               value = mean(window(ssb(stk_i@stock), 
                                   start = 2024, end = 2028) < Blim)),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk1_short",
               value = mean(window(ssb(stk_i@stock), 
                                   start = 2019, end = 2023) < Blim)),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk3_long",
               value = max(iterMeans(window(ssb(stk_i@stock), 
                                            start = 2029) < Blim))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk3_medium",
               value = max(iterMeans(window(ssb(stk_i@stock), 
                                            start = 2024, end = 2028) < Blim))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "risk3_short",
               value = max(iterMeans(window(ssb(stk_i@stock), 
                                            start = 2019, end = 2023) < Blim))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "iav_long",
               value = c(iav(object = catch(window(stock(stk_i), start = 2028))))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "iav_medium",
               value = c(iav(object = catch(window(stock(stk_i), 
                                                   start = 2023, end = 2028))))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "iav_short",
               value = c(iav(object = catch(window(stock(stk_i), 
                                                   start = 2018, end = 2023))))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "ssb_long",
               value = c(window(ssb(stk_i@stock), start = 2029))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "ssb_medium",
               value = c(window(ssb(stk_i@stock), start = 2024, end = 2028))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "ssb_short",
               value = c(window(ssb(stk_i@stock), start = 2019, end = 2023))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "recovery_proportion",
               value = mean(apply(window(ssb(stk_i@stock), 
                                         start = 2019) >= MSYBtrigger, 6, max))),
    data.frame(name = i$name, scenario = i$scenario,
               key = "recovery_time",
               value = c(apply(window(ssb(stk_i@stock), 
                                      start = 2019)@.Data >= MSYBtrigger, 6,
                               function(x) {
                                 if (any(x)) {which(x)[1]} else {Inf}})))
  )
  if (i$HCR == "F0") {
    res$value[res$key %in% c("catch_long", "catch_medium", "catch_short",
                             "iav_long", "iav_medium", "iav_short")] <- 0
  }
  res <- merge(res, i[, c("name", "OM", "HCR", "BB", "TACconstr", "Btrigger",
                          "Ftrgt")])
  return(res)
  }
  combs_full$name <- factor(combs_full$name, 
                            levels = c("F0", "A*", "A", "B", "C", "AD", "BE", "CE"))
  return(combs_full)
}

### base OM
combs_base <- stats_full(data = combs)
ggplot(data = combs_base, 
       mapping = aes(x = name, y = value, group = name)) +
  #geom_bar(stat = "identity") +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

### get median for option A*
combs_base <- left_join(combs_base, 
                         combs_base %>%
  group_by(key, OM, name) %>%
  summarise(value_median = median(value)) %>%
  filter(name == "A*") %>%
    select(-name))
p_catch_long <- ggplot(data = combs_base[combs_base$key == "catch_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "long-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_catch_medium <- ggplot(data = combs_base[combs_base$key == "catch_medium", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "medium-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_catch_short <- ggplot(data = combs_base[combs_base$key == "catch_short", ], 
                         mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "short-term catch [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1))
p_risk1_long <- ggplot(data = combs_base[combs_base$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_blank(data = combs_base[combs_base$key == "risk3_long", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1")
p_risk1_medium <- ggplot(data = combs_base[combs_base$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk3_medium", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1")
p_risk1_short <- ggplot(data = combs_base[combs_base$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk3_short", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_base[combs_base$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_long", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3")
p_risk3_medium <- ggplot(data = combs_base[combs_base$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_medium", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3")
p_risk3_short <- ggplot(data = combs_base[combs_base$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_blank(data = combs_base[combs_base$key == "risk1_short", ]) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_base[combs_base$key == "iav_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_base[combs_base$key == "iav_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_base[combs_base$key == "iav_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1.2)) + 
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_base[combs_base$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1.2e+06)) +
  labs(x = "", y = "long-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2))
p_ssb_medium <- ggplot(data = combs_base[combs_base$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.9e+06)) +
  labs(x = "", y = "medium-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_ssb_short <- ggplot(data = combs_base[combs_base$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.65e+06)) +
  labs(x = "", y = "short-term SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_recovery_proportion <- 
  ggplot(data = combs_base[combs_base$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_base[combs_base$key == "recovery_time", ], 
         mapping = aes(x = name, y = value)) +
  geom_boxplot() + theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  labs(x = "", y = "recovery time [years]")

plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long, p_ssb_long,
          align = "hv")
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium, 
          p_ssb_medium,
          align = "hv")
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short, p_ssb_short,
          align = "hv")
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(p_recovery_proportion, p_recovery_time,
          align = "hv")
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_stats/", 
                         "summary_baseOM_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### summary plots: base OM, additional scenarios around maximum yield ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations
### add: 0.9 & 1.1 * Ftrgt
###      Fmsylower, Fmsyupper
combs <- data.frame(name = rep(c("A", "B", "C", "AD", "BE", "CE"), each = 5),
                    HCR = rep(c("A", "B", "C", "A", "B", "C"), each = 5),
                    BB = rep(c(rep(FALSE, 3), rep(TRUE, 3)), each = 5),
                    TACconstr = rep(c(rep(FALSE, 3), rep(TRUE, 3)), each = 5),
                    Btrigger = rep(c(170000, 160000, 170000, 190000, 130000,
                                     140000), each = 5),
                    Ftrgt = c(0.38 * c(0.9, 1, 1.1), 0.198, 0.46,
                              0.38 * c(0.9, 1, 1.1), 0.198, 0.46,
                              0.38 * c(0.9, 1, 1.1), 0.198, 0.46,
                              0.40 * c(0.9, 1, 1.1), 0.198, 0.46,
                              0.36 * c(0.9, 1, 1.1), 0.198, 0.46,
                              0.36 * c(0.9, 1, 1.1), 0.198, 0.46),
                    scenario = c("0.9*Ftrgt", "Ftrgt", "1.1*Ftrgt",
                                 "Fmsylower", "Fmsyupper"),
                    OM = "cod4")
combs <- merge(combs, stats)
combs_dat <- stats_full(data = combs)
combs_dat$scenario <- factor(combs_dat$scenario, 
                             levels = c("Fmsylower", "0.9*Ftrgt", "Ftrgt", 
                                        "1.1*Ftrgt", "Fmsyupper"))
ggplot(data = combs_dat, 
       mapping = aes(x = name, y = value, group = interaction(scenario, name), 
                     colour = scenario)) +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

p_catch_long <- ggplot(data = combs_dat[combs_dat$key == "catch_long", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term catch [t]")
p_catch_medium <- ggplot(data = combs_dat[combs_dat$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term catch [t]")
p_catch_short <- ggplot(data = combs_dat[combs_dat$key == "catch_short", ], 
                        mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term catch [t]")
p_risk1_long <- ggplot(data = combs_dat[combs_dat$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value, fill = scenario,
                                     colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1")
p_risk1_medium <- ggplot(data = combs_dat[combs_dat$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value, fill = scenario,
                                       colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1")
p_risk1_short <- ggplot(data = combs_dat[combs_dat$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value, fill = scenario,
                                      colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk3_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_dat[combs_dat$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value, fill = scenario,
                                     colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3")
p_risk3_medium <- ggplot(data = combs_dat[combs_dat$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value, fill = scenario,
                                       colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3")
p_risk3_short <- ggplot(data = combs_dat[combs_dat$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value, fill = scenario,
                                      colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_dat[combs_dat$key == "risk1_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_dat[combs_dat$key == "iav_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.3)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_dat[combs_dat$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.3)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_dat[combs_dat$key == "iav_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_dat[combs_dat$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term SSB [t]") +
  theme(legend.direction = "horizontal")
p_ssb_medium <- ggplot(data = combs_dat[combs_dat$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  coord_cartesian(ylim = c(0, 5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term SSB [t]") +
  theme(legend.direction = "horizontal")
p_ssb_short <- ggplot(data = combs_dat[combs_dat$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 4.5e+5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA))
p_recovery_proportion <- 
  ggplot(data = combs_dat[combs_dat$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value, fill = scenario, colour = scenario)) +
  geom_bar(show.legend = FALSE, stat = "identity",
           position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_dat[combs_dat$key == "recovery_time", ], 
         mapping = aes(x = name, y = value, colour = scenario)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery time [years]")

plot_grid(plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long,
                    p_ssb_long + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_long), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/baseOM_stats_combs/", 
                         "baseOM_combs_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### summary plots: compare alternative OMs ####
### ------------------------------------------------------------------------ ###

### alternative OMs
### select maximum yield combinations
combs_alt <- data.frame(name = c("F0", "A*", "A", "B", "C", "AD", "BE", "CE"),
                    HCR = c("F0", "A", "A", "B", "C", "A", "B", "C"),
                    BB = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    TACconstr = c(rep(FALSE, 5), TRUE, TRUE, TRUE),
                    Btrigger = c(0, 150000, 170000, 160000, 170000, 190000, 130000,
                                 140000),
                    Ftrgt = c(0, 0.31, 0.38, 0.38, 0.38, 0.40, 0.36, 0.36),
                    scenario = 0)
combs_alt <- rbind(cbind(combs_alt, OM = "cod4"),
                   cbind(combs_alt, OM = "cod4_alt1"),
                   cbind(combs_alt, OM = "cod4_alt2"),
                   cbind(combs_alt, OM = "cod4_alt3"))
combs_alt <- merge(combs_alt, stats)
combs_alt <- stats_full(data = combs_alt)
ggplot(data = combs_alt, 
       mapping = aes(x = name, y = value, group = interaction(OM, name), 
                     colour = OM)) +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)

p_catch_long <- ggplot(data = combs_alt[combs_alt$key == "catch_long", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term catch [t]") +
  coord_cartesian(ylim = c(0, 2e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA))
p_catch_medium <- ggplot(data = combs_alt[combs_alt$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term catch [t]") +
  coord_cartesian(ylim = c(0, 2e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_catch_short <- ggplot(data = combs_alt[combs_alt$key == "catch_short", ], 
                        mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term catch [t]") +
  coord_cartesian(ylim = c(0, 2.1e+05)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2),
                     limits = c(0, NA))
p_risk1_long <- ggplot(data = combs_alt[combs_alt$key == "risk1_long", ], 
                       mapping = aes(x = name, y = value, fill = OM, 
                                     colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_long", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 1") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk1_medium <- ggplot(data = combs_alt[combs_alt$key == "risk1_medium", ], 
                         mapping = aes(x = name, y = value, fill = OM,
                                       colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_medium", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 1") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk1_short <- ggplot(data = combs_alt[combs_alt$key == "risk1_short", ], 
                        mapping = aes(x = name, y = value, fill = OM,
                                      colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  geom_blank(data = combs_alt[combs_alt$key == "risk3_short", ]) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 1")
p_risk3_long <- ggplot(data = combs_alt[combs_alt$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value, fill = OM,
                                     colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term risk 3") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk3_medium <- ggplot(data = combs_alt[combs_alt$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value, fill = OM,
                                       colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term risk 3") +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1))
p_risk3_short <- ggplot(data = combs_alt[combs_alt$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value, fill = OM,
                                      colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity", 
           position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term risk 3")
p_iav_long <- ggplot(data = combs_alt[combs_alt$key == "iav_long", ], 
                     mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term inter-annual catch variability")
p_iav_medium <- ggplot(data = combs_alt[combs_alt$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term inter-annual catch variability")
p_iav_short <- ggplot(data = combs_alt[combs_alt$key == "iav_short", ], 
                      mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = FALSE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term inter-annual catch variability")
p_ssb_long <- ggplot(data = combs_alt[combs_alt$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "long-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  coord_cartesian(ylim = c(0, 2.3e+06))
p_ssb_medium <- ggplot(data = combs_alt[combs_alt$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) + 
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "medium-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  coord_cartesian(ylim = c(0, 1.8e+06))
p_ssb_short <- ggplot(data = combs_alt[combs_alt$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "short-term SSB [t]") +
  theme(legend.direction = "horizontal") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 1),
                     limits = c(0, NA)) +
  coord_cartesian(ylim = c(0, 8e+05))
p_recovery_proportion <- 
  ggplot(data = combs_alt[combs_alt$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value, fill = OM, colour = OM)) +
  geom_bar(show.legend = FALSE, stat = "identity",
           position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery proportion")
p_recovery_time <- 
  ggplot(data = combs_alt[combs_alt$key == "recovery_time", ], 
         mapping = aes(x = name, y = value, colour = OM)) +
  geom_boxplot(show.legend = TRUE, size = 0.25, outlier.size = 0.3,
               position = position_dodge2(preserve = "single")) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "recovery time [years]")

plot_grid(plot_grid(p_catch_long, p_risk1_long, p_risk3_long, p_iav_long,
                    p_ssb_long + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_long), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_long.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_medium, p_risk1_medium, p_risk3_medium, p_iav_medium,
                    p_ssb_medium + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_medium), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_medium.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_catch_short, p_risk1_short, p_risk3_short, p_iav_short,
                    p_ssb_short + theme(legend.position = "none"),
                    align = "hv"),
          get_legend(p_ssb_short), nrow = 2, rel_heights = c(1, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_short.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
plot_grid(plot_grid(p_recovery_proportion, 
                    p_recovery_time + theme(legend.position = "none"),
          align = "hv"),
          get_legend(p_recovery_time), ncol = 2, rel_widths = c(0.5, 0.1))
ggsave(filename = paste0("output/runs/cod4/1000_20/plots/altOMs_stats/", 
                         "summary_altOM_recovery.png"), 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


