### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)

### ------------------------------------------------------------------------ ###
### option A, first coarse grid optimisation ####
### ------------------------------------------------------------------------ ###

path_res <- "output/runs/cod4/1000_20/"
files_res <- data.frame(file = list.files(path_res, pattern = "*.rds"), 
                        stringsAsFactors = FALSE)

files_res$Ftrgt <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                     split = "_", fixed = TRUE),
                                            "[[", 2), gsub,
                                     pattern = "Ftrgt-", replacement = ""))
files_res$Btrigger <- as.numeric(lapply(lapply(strsplit(files_res$file, 
                                                     split = "_", fixed = TRUE),
                                            "[[", 3), gsub,
                                     pattern = "Btrigger-", replacement = ""))


res_list <- vector(mode = "list", length = nrow(files_res))
for (i in seq(nrow(files_res))) {
  res_list[[i]] <- readRDS(paste0(path_res, files_res$file[i]))
}


### last 10 years
### sum catch
files_res$catch_last10 <- sapply(res_list, function(x) {
  mean(window(catch(x@stock), start = 2029))
})
### risks
files_res$risk1_mean <- sapply(res_list, function(x) {
  mean(window(ssb(x@stock), start = 2019) < 107000)
})
files_res$risk1_last10 <- sapply(res_list, function(x) {
  mean(window(ssb(x@stock), start = 2029) < 107000)
})
files_res$risk3_last10 <- sapply(res_list, function(x) {
  max(iterMeans(window(ssb(x@stock), start = 2029) < 107000))
})
### inter-annual variation
files_res$iav_last10 <- sapply(res_list, function(x) {
  iav(object = catch(window(stock(x), start = 2028)), summary_per_iter = mean,
      summary = mean)
})

### catch
p1 <- ggplot(data = files_res, 
             aes(x = Btrigger, y = Ftrgt, fill = catch_last10)) +
  geom_raster() +
  scale_fill_gradient("mean catch\nlast 10 yrs", low = "red", high = "green") +
  geom_text(aes(label = round(catch_last10), colour = risk3_last10 <= 0.05),
            size = 2) +
  scale_colour_manual("risk <= 0.05", 
                      values = c("FALSE" = "red", "TRUE" = "black")) +
  theme_bw()
### risk
p2 <- ggplot(data = files_res, 
             aes(x = Btrigger, y = Ftrgt, fill = risk3_last10)) +
  geom_raster(alpha = 0.75) +
  geom_text(aes(label = round(risk3_last10, 3), colour = risk3_last10 <= 0.05),
            size = 2) +
  scale_colour_manual("risk <= 0.05", 
                      values = c("FALSE" = "red", "TRUE" = "black")) +
  scale_fill_gradient("risk 3\nlast 10 yrs", low = "green", high = "red") +
  theme_bw()
### iav
p3 <- ggplot(data = files_res, 
             aes(x = Btrigger, y = Ftrgt, fill = iav_last10)) +
  geom_raster() +
  geom_text(aes(label = round(iav_last10, 3), colour = risk3_last10 <= 0.05),
            size = 2) +
  scale_colour_manual("risk <= 0.05",
                      values = c("FALSE" = "red", "TRUE" = "black")) +
  scale_fill_gradient("mean inter-\nannual\nvariability\nof catch\nlast 10 yrs", 
                      low = "green", high = "red") +
  theme_bw()

plot_grid(p1, p2, p3, nrow = 2, ncol = 2, align = "v")
ggsave(filename = "output/runs/cod4/1000_20/coarse_grid2.png", 
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

