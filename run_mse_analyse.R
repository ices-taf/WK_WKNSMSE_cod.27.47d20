### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)

### load additional functions
source("a4a_mse_WKNSMSE_funs.R")

### ------------------------------------------------------------------------ ###
### compare default HCR: full vs shortcut ####
### ------------------------------------------------------------------------ ###

### summary statistics
full_stats <- readRDS("../WK_WKNSMSE_cod.27.47d20/output/runs/cod4/1000_20/stats_full.rds")
short_stats <- readRDS("output/runs/cod4/1000_20/stats_cod4_HCR-A_Ftrgt-0.31_Btrigger-150000_TACconstr-FALSE_BB-FALSE.rds")
full_stats %>% 
  filter(OM == "cod4" & Ftrgt == 0.31 & Btrigger == 150000 & HCR == "A" &
           TACconstr == FALSE & BB == FALSE) %>%
  mutate(assessment = "SAM") %>%
  bind_rows(as.data.frame(short_stats) %>%
              mutate(assessment = "shortcut")) %>%
  select(catch_median_short, catch_median_medium, catch_median_long,
         ssb_median_short, ssb_median_medium, ssb_median_long,
         fbar_median_short, fbar_median_medium, fbar_median_long,
         risk1_short, risk1_medium, risk1_long,
         risk3_short, risk3_medium, risk3_long,
         iav_short, iav_medium, iav_long,
         assessment
         ) %>%
  pivot_longer(c(1:18)) %>%
  mutate(time = ifelse(grepl(x = name, pattern = "short"), "short", NA),
         time = ifelse(grepl(x = name, pattern = "medium"), "medium", time),
         time = ifelse(grepl(x = name, pattern = "long"), "long", time),
         stat = sapply(strsplit(x = name, split = "_"), "[[", 1)) %>%
  mutate(time = factor(time, levels = c("short", "medium", "long"))) %>%
  mutate(value = ifelse(stat %in% c("catch", "ssb"), value/1000, value)) %>%
  mutate(stat = factor(stat, levels = c("risk1", "risk3", "ssb", "catch", 
                                         "fbar", "iav"),
                        labels = c("Risk 1", "Risk 3", "SSB [1000t]", 
                                   "Catch [1000t]", "F (ages 2-4)", "ICV"))) %>%
  ggplot(aes(x = time, y = value, fill = assessment)) +
  geom_col(position = "dodge") +
  facet_wrap(~ stat, scales = "free_y", strip.position = "left") +
  theme_bw(base_size = 8) +
  ylim(c(0, NA)) +
  labs(y = "") +
  theme(strip.placement = "outside",
        strip.background = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave(filename = "output/plots/default_SAM_vs_shortcut_stats.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")

### time series
input_full <- readRDS("../WK_WKNSMSE_cod.27.47d20/input/cod4/1000_20/base_run.rds")
input_short <- readRDS("input/cod4/1000_20/base_run.rds")
input_short_catch_n <- readRDS("input/cod4/1000_20/catch_n.rds")
### historical period identical
res_full <- readRDS("../WK_WKNSMSE_cod.27.47d20/output/runs/cod4/1000_20/cod4_HCR-A_Ftrgt-0.31_Btrigger-150000_TACconstr-FALSE_BB-FALSE.rds")
res_short <- readRDS("output/runs/cod4/1000_20/MP_cod4_HCR-A_Ftrgt-0.31_Btrigger-150000_TACconstr-FALSE_BB-FALSE.rds")

### combine history and projection
stk_full <- input_short$om@stock
stk_full[, dimnames(res_full@stock)$year] <- res_full@stock
catch.n(stk_full)[, dimnames(input_short_catch_n)$year] <- input_short_catch_n
catch(stk_full) <- computeCatch(stk_full)
stk_short <- stk_full
stk_short[, dimnames(res_short@stock)$year] <- res_short@stock

probs <- c(0.05, 0.5, 0.95)
iters <- c(1:3)
qnts <- list(SAM = stk_full, shortcut = stk_short)
qnts <- lapply(qnts, function(x) {
  FLQuants(rec = rec(x), ssb = ssb(x), catch = catch(x), fbar = fbar(x))
})
qnts <- lapply(qnts, function(x) {
  lapply(x, function(y) {
    res <- y[,,,,, c(seq(length(probs)), iters)]
    res[,,,,, c(seq(length(probs)))] <- quantile(y, probs = probs)
    dimnames(res)$iter <- c(probs, iters)
    return(res)
  })
})
df <- lapply(names(qnts), function(x) cbind(as.data.frame(qnts[[x]]), 
                                            assessment = x))
df <- do.call(rbind, df)
df <- df %>% 
  select(year, iter, data, qname, assessment) %>%
  mutate(data = ifelse(qname %in% c("ssb", "catch"), data/1000, data)) %>%
  mutate(data = ifelse(qname %in% c("rec"), data/1e+6, data)) %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname, levels = c("catch", "rec", "fbar", "ssb"),
                        labels = c("Catch [1000t]", "Recruitment [millions]",
                                   "F (ages 2-4)", "SSB [1000t]")))
df %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `0.05`, ymax = `0.95`, fill = assessment),
              alpha = 0.1) +
  scale_fill_manual(values = c(SAM = "black", shortcut = "red")) +
  geom_line(aes(x = year, y = `0.5`, linetype = assessment, colour = assessment)) +
  geom_line(aes(x = year, y = `0.05`, linetype = assessment, colour = assessment), 
            show.legend = FALSE, alpha = 0.2) +
  geom_line(aes(x = year, y = `0.95`, linetype = assessment, colour = assessment), 
            show.legend = FALSE, alpha = 0.2) +
  scale_colour_manual(values = c(SAM = "black", shortcut = "red")) +
  geom_line(aes(x = year, y = `1`, linetype = assessment), 
            alpha = 0.1, size = 0.2, colour = scales::hue_pal()(3)[1]) +
  geom_line(aes(x = year, y = `2`, linetype = assessment), 
            alpha = 0.1, size = 0.2, colour = scales::hue_pal()(3)[2]) +
  geom_line(aes(x = year, y = `3`, linetype = assessment), 
            alpha = 0.1, size = 0.2, colour = scales::hue_pal()(3)[3]) +
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background = element_blank()) +
  labs(x = "year", y = "") +
  ylim(c(0, NA)) +
  geom_vline(xintercept = 2018.5, size = 0.2)
ggsave(filename = "output/plots/default_SAM_vs_shortcut_timeseries.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
df %>%
  filter(year >= 2010) %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `0.05`, ymax = `0.95`, fill = assessment),
              alpha = 0.1) +
  scale_fill_manual(values = c(SAM = "black", shortcut = "red")) +
  geom_line(aes(x = year, y = `0.5`, linetype = assessment, colour = assessment)) +
  geom_line(aes(x = year, y = `0.05`, linetype = assessment, colour = assessment), 
            show.legend = FALSE, alpha = 0.2) +
  geom_line(aes(x = year, y = `0.95`, linetype = assessment, colour = assessment), 
            show.legend = FALSE, alpha = 0.2) +
  scale_colour_manual(values = c(SAM = "black", shortcut = "red")) +
  geom_line(aes(x = year, y = `1`, linetype = assessment), 
            alpha = 0.2, size = 0.2, colour = scales::hue_pal()(3)[1]) +
  geom_line(aes(x = year, y = `2`, linetype = assessment), 
            alpha = 0.2, size = 0.2, colour = scales::hue_pal()(3)[2]) +
  geom_line(aes(x = year, y = `3`, linetype = assessment), 
            alpha = 0.2, size = 0.2, colour = scales::hue_pal()(3)[3]) +
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background = element_blank()) +
  labs(x = "year", y = "") +
  ylim(c(0, NA)) +
  geom_vline(xintercept = 2018.5, size = 0.2)
ggsave(filename = "output/plots/default_SAM_vs_shortcut_timeseries_zoom.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")


### ------------------------------------------------------------------------ ###
### SAM vs. shortcut - compare grid ####
### ------------------------------------------------------------------------ ###

### collate stats (on HPC)
# files <- list.files(path = "output/runs/cod4/1000_20/", pattern = "stats_cod",
#                     full.names = TRUE)
# stats_short <- lapply(files, readRDS)
# stats_short <- data.frame(do.call(bind_rows, stats_short))
# stats_short <- data.frame(lapply(stats_short, unlist))
# saveRDS(stats_short, file = "output/runs/cod4/1000_20/stats_combined.rds")
stats_short <- readRDS("output/runs/cod4/1000_20/stats_combined.rds")

stats_combined <- full_stats %>% 
  mutate(assessment = "OM1: SAM") %>%
  filter(Ftrgt %in% round(seq(0, 1, 0.01), 2)) %>%
  bind_rows(stats_short %>% 
              mutate(assessment = "OM1: shortcut") %>%
              filter(is.na(obs_sd))
  ) %>%
  filter(OM == "cod4" &
           HCR == "A" & TACconstr == FALSE & BB == FALSE) %>%
  select(Ftrgt, Btrigger, catch_median_long, risk3_long, assessment) %>%
  mutate(Ftrgt = round(Ftrgt, 2)) %>%
  mutate(risk_pos = ifelse(risk3_long <= 0.05, "below", "above"))
### find yield maximum
stats_combined_max <- stats_combined %>%
  group_by(assessment) %>%
  filter(risk3_long <= 0.05) %>%
  filter(catch_median_long == max(catch_median_long)) %>%
  ungroup() %>%
  mutate(assessment2 = assessment, assessment = NULL)
stats_combined_max
stats_combined_max_both <- stats_combined %>%
  filter((Ftrgt == 0.38 & Btrigger == 170000) | 
           (Ftrgt == 0.39 & Btrigger == 160000)) %>%
  select(Ftrgt, Btrigger)
# stats_combined <- stats_combined %>%
#   full_join(stats_combined_max) %>%
#   mutate(risk_pos = ifelse(!is.na(assessment2), "optimum", risk_pos),
#          assessment2 = NULL)
### Ftrgt steps where risk exceeds 5%
stats_combined_step <- stats_combined %>%
  group_by(assessment, Btrigger) %>%
  filter(risk3_long <= 0.05) %>%
  filter(risk3_long == max(risk3_long)) %>%
  ungroup() %>%
  arrange(Btrigger)
stats_combined_step <- stats_combined_step %>%
  mutate(Btrigger = Btrigger - 4999, Ftrgt = Ftrgt + 0.005) %>%
  bind_rows(
    stats_combined_step %>%
      mutate(Btrigger = Btrigger + 4999, Ftrgt = Ftrgt + 0.005)) %>%
  mutate(assessment2 = assessment, assessment = NULL)


ggplot() +
  geom_raster(data = stats_combined %>% 
                filter(risk3_long <= 0.05) %>%
                filter(catch_median_long >= 0.95 * max(catch_median_long)),
              aes(x = Btrigger, y = Ftrgt, fill = catch_median_long)) +
  scale_fill_gradient(paste0("yield maximum\narea [t]"), low = "red",
                      high = "green") +
  geom_text(data = stats_combined, 
            aes(x = Btrigger, y = Ftrgt, 
                label = round(catch_median_long), colour = risk3_long <= 0.05),
            size = 1.2, show.legend = FALSE) +
  scale_colour_manual("risk <= 0.05", 
                      values = c("FALSE" = "red", "TRUE" = "black",
                                 "OM1: SAM" = "black", "OM1: shortcut" = "blue")) +
  geom_line(data = stats_combined_step, 
            aes(x = Btrigger, y = Ftrgt, colour = assessment2),
            show.legend = FALSE, size = 0.3) +
  geom_tile(data = stats_combined_max, 
            aes(x = Btrigger, y = Ftrgt, colour = assessment2),
            width = 10000, height = 0.01,
            alpha = 0, colour = "black", size = 0.3) +
  geom_tile(data = data.frame(x = 170000, y = 0.38),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "black", linetype = 0) +
  geom_tile(data = data.frame(x = 160000, y = 0.39),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "blue", linetype = 0) +
  geom_text(data = stats_combined %>%
              filter((Btrigger == 170000 & Ftrgt == 0.38) |
                       (Btrigger == 160000 & Ftrgt == 0.39)),
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  scale_linetype_discrete("optimum") + 
  theme_bw() +
  facet_wrap(~ assessment) +
  scale_x_continuous(breaks = c(seq(from = 110000, to = 210000, by = 20000)),
                     labels = c(seq(from = 110000, to = 210000, by = 20000))/1000) +
  labs(x = expression(B[trigger]~"[1000t]"),
       y = expression(F[trgt]))
ggsave(filename = "output/plots/default_SAM_vs_shortcut_grid_A_both.png", 
       width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")


p <- ggplot() +
  geom_raster(data = stats_combined %>% 
                filter(risk3_long <= 0.05) %>%
                filter(catch_median_long >= 0.95 * max(catch_median_long)),
              aes(x = Btrigger, y = Ftrgt, fill = catch_median_long)) +
  scale_fill_gradient(paste0("yield maximum\narea [t]"), low = "red",
                      high = "green") +
  geom_text(data = stats_combined, 
            aes(x = Btrigger, y = Ftrgt, 
                label = round(catch_median_long), colour = risk3_long <= 0.05),
            size = 1.2, show.legend = FALSE) +
  scale_colour_manual("risk <= 0.05", 
                      values = c("FALSE" = "red", "TRUE" = "black",
                                 "OM1: SAM" = "black", "OM1: shortcut" = "blue")) +
  geom_line(data = stats_combined_step %>% 
              mutate(assessment = assessment2), 
            aes(x = Btrigger, y = Ftrgt, colour = assessment2),
            show.legend = FALSE, size = 0.3) +
  # geom_tile(data = stats_combined_max, 
  #           aes(x = Btrigger, y = Ftrgt, colour = assessment2),
  #           width = 10000, height = 0.01,
  #           alpha = 0, colour = "black", size = 0.3) +
  geom_tile(data = data.frame(x = 170000, y = 0.38, assessment = "OM1: SAM"),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "black", linetype = 0) +
  geom_tile(data = data.frame(x = 160000, y = 0.39, assessment = "OM1: shortcut"),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "blue", linetype = 0) +
  geom_text(data = stats_combined %>%
              filter((Btrigger == 170000 & Ftrgt == 0.38 & assessment == "OM1: SAM")),
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  geom_text(data = stats_combined %>%
              filter((Btrigger == 160000 & Ftrgt == 0.39 & assessment == "OM1: shortcut")),
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  scale_linetype_discrete("optimum") + 
  theme_bw() +
  facet_wrap(~ assessment) +
  scale_x_continuous(breaks = c(seq(from = 110000, to = 210000, by = 20000)),
                     labels = c(seq(from = 110000, to = 210000, by = 20000))/1000) +
  labs(x = expression(B[trigger]~"[1000t]"),
       y = expression(F[trgt])) +
  theme(legend.position = "left")
p
ggsave(filename = "output/plots/default_SAM_vs_shortcut_grid_A_1.png", 
       width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")
p <- p + geom_tile(data = data.frame(x = 160000, y = 0.39, assessment = "OM1: SAM"),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "blue", linetype = 0) +
  geom_text(data = stats_combined %>%
              filter((Btrigger == 160000 & Ftrgt == 0.39 & assessment == "OM1: SAM")),
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  geom_line(data = stats_combined_step %>% 
              filter(assessment2 == "OM1: shortcut") %>%
              mutate(assessment = "OM1: SAM"), 
            aes(x = Btrigger, y = Ftrgt, colour = assessment2),
            show.legend = FALSE, size = 0.3)
p
ggsave(filename = "output/plots/default_SAM_vs_shortcut_grid_A_2.png", 
       width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")
p <- p + geom_tile(data = data.frame(x = 170000, y = 0.38, assessment = "OM1: shortcut"),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "black", linetype = 0) +
  geom_text(data = stats_combined %>%
              filter((Btrigger == 170000 & Ftrgt == 0.38 & assessment == "OM1: shortcut")),
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  geom_line(data = stats_combined_step %>% 
              filter(assessment2 == "OM1: SAM") %>%
              mutate(assessment = "OM1: shortcut"), 
            aes(x = Btrigger, y = Ftrgt, colour = assessment2),
            show.legend = FALSE, size = 0.3)
p
ggsave(filename = "output/plots/default_SAM_vs_shortcut_grid_A_3.png", 
       width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")

### ------------------------------------------------------------------------ ###
### alternative OMs - grids ####
### ------------------------------------------------------------------------ ###
stats_short <- readRDS("output/runs/cod4/1000_20/stats_combined.rds")
stats_OMs <- stats_short %>% 
  filter(is.na(obs_sd)) %>%
  select(OM, Ftrgt, Btrigger, catch_median_long, risk3_long) %>%
  mutate(Ftrgt = round(Ftrgt, 2)) %>%
  mutate(risk_pos = ifelse(risk3_long <= 0.05, "below", "above")) %>%
  mutate(OM = factor(OM, 
                     levels = c("cod4", "cod4_alt1", "cod4_alt2", "cod4_alt3"),
                     labels = c("OM1 (baseline)", "OM2 (recruitment)",
                                "OM3 (year effects)", "OM4 (dens. dep. M)")))
### find yield maximum
stats_OM_max <- stats_OMs %>%
  group_by(OM) %>%
  filter(risk3_long <= 0.05) %>%
  filter(catch_median_long == max(catch_median_long)) %>%
  ungroup() %>%
  mutate(OM_ = OM)
stats_OM_max
### Ftrgt steps where risk exceeds 5%
stats_OM_step <- stats_OMs %>%
  group_by(OM, Btrigger) %>%
  filter(risk3_long <= 0.05) %>%
  filter(risk3_long == max(risk3_long)) %>%
  ungroup() %>%
  arrange(Btrigger)
stats_OM_step <- stats_OM_step %>%
  mutate(Btrigger = Btrigger - 4999, Ftrgt = Ftrgt + 0.005) %>%
  bind_rows(
    stats_OM_step %>%
      mutate(Btrigger = Btrigger + 4999, Ftrgt = Ftrgt + 0.005)) %>%
  mutate(OM_ = OM)

p <- ggplot() +
  geom_raster(data = stats_OMs %>% 
                group_by(OM) %>%
                filter(risk3_long <= 0.05) %>%
                filter(catch_median_long >= 0.95 * max(catch_median_long)) %>%
                mutate(catch_median_long = catch_median_long/max(catch_median_long)),
              aes(x = Btrigger, y = Ftrgt, fill = catch_median_long)) +
  scale_fill_gradient(paste0("yield maximum\narea (relative)"), low = "red",
                      high = "green") +
  geom_text(data = stats_OMs, 
            aes(x = Btrigger, y = Ftrgt, 
                label = round(catch_median_long), colour = risk3_long <= 0.05),
            size = 1.2, show.legend = FALSE) +
  scale_colour_manual("risk <= 0.05",
                      values = c("FALSE" = "red", "TRUE" = "black",
                                 "OM1 (baseline)" = "black", 
                                 "OM2 (recruitment)" = "blue",
                                 "OM3 (year effects)" = "darkgreen", 
                                 "OM4 (dens. dep. M)" = "brown")) +
  geom_line(data = stats_OM_step, 
            aes(x = Btrigger, y = Ftrgt, colour = OM_),
            show.legend = FALSE, size = 0.3) +
  new_scale_fill() +
  geom_tile(data = stats_OM_max, 
            aes(x = Btrigger, y = Ftrgt, fill = OM_),
            width = 10000, height = 0.01,
            size = 0.3) +
  scale_fill_manual("optimum", 
                    values = c("OM1 (baseline)" = "black", 
                               "OM2 (recruitment)" = "blue",
                               "OM3 (year effects)" = "darkgreen", 
                               "OM4 (dens. dep. M)" = "brown")) +
  geom_text(data = stats_OM_max,
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  theme_bw() +
  facet_wrap(~ OM) +
  scale_x_continuous(breaks = c(seq(from = 110000, to = 210000, by = 20000)),
                     labels = c(seq(from = 110000, to = 210000, by = 20000))/1000) +
  labs(x = expression(B[trigger]~"[1000t]"),
       y = expression(F[trgt])) +
  ylim(c(0.25, NA))
p
ggsave(filename = "output/plots/default_SAM_vs_shortcut_grid_OMs.png", 
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
### add alternative OMs steps to OM1
p + geom_line(data = stats_OM_step %>% 
                mutate(OM = "OM1 (baseline)"), 
              aes(x = Btrigger, y = Ftrgt, colour = OM_),
              show.legend = FALSE, size = 0.3)
ggsave(filename = "output/plots/default_SAM_vs_shortcut_grid_OMs_2.png", 
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")


### ------------------------------------------------------------------------ ###
### OM3: SAM vs. shortcut - compare grid ####
### ------------------------------------------------------------------------ ###

### shortcut
stats_short <- readRDS("output/runs/cod4/1000_20/stats_combined.rds")
### full MSE 
stats_full_OM3 <- readRDS("../WK_WKNSMSE_cod.27.47d20_mse2.0/output/runs/cod4/1000_20/stats_combined_alt2.rds")

stats_OM3_combined <- stats_full_OM3 %>% 
  mutate(assessment = "OM3: SAM") %>%
  bind_rows(stats_short %>% 
              mutate(assessment = "OM3: shortcut") %>%
              filter(is.na(obs_sd))
  ) %>%
  filter(OM == "cod4_alt2" &
           HCR == "A" & TACconstr == FALSE & BB == FALSE) %>%
  select(Ftrgt, Btrigger, catch_median_long, risk3_long, assessment) %>%
  mutate(Ftrgt = round(Ftrgt, 2)) %>%
  mutate(risk_pos = ifelse(risk3_long <= 0.05, "below", "above"))
### find yield maximum
stats_OM3_combined_max <- stats_OM3_combined %>%
  group_by(assessment) %>%
  filter(risk3_long <= 0.05) %>%
  filter(catch_median_long == max(catch_median_long)) %>%
  ungroup() %>%
  mutate(assessment2 = assessment, assessment = NULL)
stats_OM3_combined_max
stats_combined_max_both <- stats_OM3_combined %>%
  filter((Ftrgt == 0.35 & Btrigger == 180000) | 
           (Ftrgt == 0.37 & Btrigger == 170000)) %>%
  select(Ftrgt, Btrigger)
### Ftrgt steps where risk exceeds 5%
stats_OM3_combined_step <- stats_OM3_combined %>%
  group_by(assessment, Btrigger) %>%
  filter(risk3_long <= 0.05) %>%
  filter(risk3_long == max(risk3_long)) %>%
  ungroup() %>%
  arrange(Btrigger)
stats_OM3_combined_step <- stats_OM3_combined_step %>%
  mutate(Btrigger = Btrigger - 4999, Ftrgt = Ftrgt + 0.005) %>%
  bind_rows(
    stats_OM3_combined_step %>%
      mutate(Btrigger = Btrigger + 4999, Ftrgt = Ftrgt + 0.005)) %>%
  mutate(assessment2 = assessment, assessment = NULL) %>%
  mutate(Ftrgt = ifelse(assessment2 == "OM3: SAM" & Ftrgt <= 0.33, NA, Ftrgt))


p <- ggplot() +
  geom_raster(data = stats_OM3_combined %>% 
                filter(risk3_long <= 0.05) %>%
                filter(catch_median_long >= 0.95 * max(catch_median_long)),
              aes(x = Btrigger, y = Ftrgt, fill = catch_median_long)) +
  scale_fill_gradient(paste0("yield maximum\narea [t]"), low = "red",
                      high = "green") +
  geom_text(data = stats_OM3_combined, 
            aes(x = Btrigger, y = Ftrgt, 
                label = round(catch_median_long), colour = risk3_long <= 0.05),
            size = 1.2, show.legend = FALSE) +
  scale_colour_manual("risk <= 0.05", 
                      values = c("FALSE" = "red", "TRUE" = "black",
                                 "OM3: SAM" = "black", "OM3: shortcut" = "blue")) +
  geom_line(data = stats_OM3_combined_step %>% 
              mutate(assessment = assessment2), 
            aes(x = Btrigger, y = Ftrgt, colour = assessment2),
            show.legend = FALSE, size = 0.3) +
  geom_tile(data = data.frame(x = 180000, y = 0.35, assessment = "OM3: SAM"),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "black", linetype = 0) +
  geom_tile(data = data.frame(x = 170000, y = 0.37, assessment = "OM3: shortcut"),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "blue", linetype = 0) +
  geom_text(data = stats_OM3_combined %>%
              filter((Btrigger == 180000 & Ftrgt == 0.35 & assessment == "OM3: SAM")),
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  geom_text(data = stats_OM3_combined %>%
              filter((Btrigger == 170000 & Ftrgt == 0.37 & assessment == "OM3: shortcut")),
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  scale_linetype_discrete("optimum") + 
  theme_bw() +
  facet_wrap(~ assessment, nrow = 2) +
  scale_x_continuous(breaks = c(seq(from = 110000, to = 210000, by = 20000)),
                     labels = c(seq(from = 110000, to = 210000, by = 20000))/1000) +
  labs(x = expression(B[trigger]~"[1000t]"),
       y = expression(F[trgt])) +
  ylim(c(0.24, 0.42))
p
ggsave(filename = "output/plots/OM3_SAM_vs_shortcut_grid_A_1.png", 
       width = 10, height = 10, units = "cm", dpi = 600, type = "cairo")
p <- p + geom_tile(data = data.frame(x = 170000, y = 0.37, assessment = "OM3: SAM"),
            aes(x = x, y = y), width = 10000, height = 0.01,
            fill = "blue", linetype = 0) +
  geom_text(data = stats_OM3_combined %>%
              filter((Btrigger == 170000 & Ftrgt == 0.37 & assessment == "OM3: SAM")),
            aes(x = Btrigger, y = Ftrgt,
                label = round(catch_median_long)),
            size = 1.2, colour = "white", show.legend = FALSE) +
  geom_line(data = stats_OM3_combined_step %>% 
              filter(assessment2 == "OM3: shortcut") %>%
              mutate(assessment = "OM3: SAM"), 
            aes(x = Btrigger, y = Ftrgt, colour = assessment2),
            show.legend = FALSE, size = 0.3)
p
ggsave(filename = "output/plots/OM3_SAM_vs_shortcut_grid_A_2.png", 
       width = 10, height = 10, units = "cm", dpi = 600, type = "cairo")

### ------------------------------------------------------------------------ ###
### shortcut - compare assessment uncertainty ####
### ------------------------------------------------------------------------ ###

### results
stats_short <- readRDS("output/runs/cod4/1000_20/stats_combined.rds")
### default uncertainty values
obs_sd_def <- readRDS(paste0(path_data, "obs_sd.rds"))
obs_rho_def <- readRDS(paste0(path_data, "obs_rho.rds"))

stats_unc <- stats_short %>%
  filter(!is.na(obs_sd) & !is.na(obs_rho)) %>% 
  select(Ftrgt, Btrigger, catch_median_long, risk3_long, obs_sd, obs_rho) %>%
  mutate(Ftrgt = round(Ftrgt, 2)) %>%
  mutate(risk_pos = ifelse(risk3_long <= 0.05, "below", "above")) 
stats_unc <- bind_rows(
  stats_unc %>% 
    filter(round(obs_sd, 3) == round(obs_sd_def, 3)) %>% 
    mutate(scenario = "rho",
           uncertainty = obs_rho),
  stats_unc %>% 
    filter(round(obs_rho, 3) == round(obs_rho_def, 3)) %>% 
    mutate(scenario = "sd",
           uncertainty = obs_sd)
  ) %>%
  mutate(scenario = factor(scenario, levels = c("sd", "rho"),
                           labels = c("uncertainty~(italic(sd))",
                                      "autocorrelation~(italic(rho))")))
### find yield maximum
stats_unc_max <- stats_unc %>%
  group_by(obs_sd, obs_rho) %>%
  filter(risk3_long <= 0.05) %>%
  filter(catch_median_long == max(catch_median_long))
stats_unc_max
### Ftrgt steps where risk exceeds 5%
stats_unc_step <- stats_unc %>%
  group_by(uncertainty, obs_rho, obs_sd, Btrigger, scenario) %>%
  filter(risk3_long <= 0.05) %>%
  filter(risk3_long == max(risk3_long)) %>%
  ungroup() %>%
  arrange(Btrigger)
# stats_unc_step %>%
#   filter(scenario == "rho") %>%
#   filter(Btrigger == 150000)
stats_unc_step <- stats_unc_step %>%
  mutate(Btrigger = Btrigger - 4999, Ftrgt = Ftrgt + 0.005) %>%
  bind_rows(
    stats_unc_step %>%
      mutate(Btrigger = Btrigger + 4999, Ftrgt = Ftrgt + 0.005))
### plot 5% risk steps
p <- stats_unc_step %>%
  ggplot(aes(x = Btrigger, y = Ftrgt, colour = uncertainty, group = uncertainty)) +
  geom_line() +
  scale_colour_gradientn("uncertainty", 
                     colours = c(scales::hue_pal()(6)), # rainbow(n = 6, rev = TRUE)
                     values = c(0, 0.1, 0.2, 0.3, 0.5, 1)
                     ) +
  geom_line(data = stats_unc_step %>%
              filter(obs_rho == obs_rho_def & obs_sd == obs_sd_def),
            aes(linetype = "default"),
            colour = "black", show.legend = TRUE) +
  scale_linetype("") +
  theme_bw() +
  facet_wrap(~ scenario, labeller = "label_parsed") +
  scale_x_continuous(breaks = c(seq(from = 110000, to = 210000, by = 20000)),
                     labels = c(seq(from = 110000, to = 210000, by = 20000))/1000) +
  labs(x = expression(B[trigger]~"[1000t]"),
       y = expression(F[trgt])) +
  geom_tile(data = stats_unc_max,
              aes(x = Btrigger, y = Ftrgt), 
              width = 10000, height = 0.01,
              linetype = 0, alpha = 0) # invisble, for y-axis range
p
ggsave(filename = "output/plots/shortcut_A_uncertainty_risk_steps.png", 
       width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")
### add catch maxima
p + geom_tile(data = stats_unc_max,
              aes(x = Btrigger, y = Ftrgt, fill = uncertainty, group = uncertainty), 
              width = 10000, height = 0.01,
              linetype = 0) +
  scale_fill_gradientn("uncertainty", 
                       colours = c(scales::hue_pal()(6)),
                       values = c(0, 0.1, 0.2, 0.3, 0.5, 1)) +
  geom_line(data = stats_unc_step %>%
              filter(obs_rho == obs_rho_def & obs_sd == obs_sd_def),
            aes(linetype = "default"),
            colour = "black", show.legend = TRUE) +
  geom_tile(data = stats_unc_max %>%
              filter(obs_rho == obs_rho_def & obs_sd == obs_sd_def),
            aes(x = Btrigger, y = Ftrgt), 
            width = 10000, height = 0.01, linetype = 0, fill = "black") +
  geom_text(data = stats_unc_max,
            aes(x = Btrigger, y = Ftrgt, label = round(uncertainty, 2)), 
            colour = "white", size = 1.2, check_overlap = TRUE)
ggsave(filename = "output/plots/shortcut_A_uncertainty_risk_steps_2.png", 
       width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")

