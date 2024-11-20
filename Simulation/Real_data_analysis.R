rm(list=ls())
library(matrixStats)
library(Rcpp)
library(RcppArmadillo)
library(changepoint.geo)
library(changepoint.np)
library(InspectChangepoint)
library(ecp)
library(dplyr)
library(tidyr)
library(RSpectra)
library(hdbcp)
library(ggplot2)
library(reshape2)
library(SP500R)
source("util_funcs.R")

# Parameters
## mxPBF
n_cores <- 1 # The number of cores for parallelization
nws <- c(25, 60, 100) # The set of window sizes

# Dataset
## Genomic Data
data("ACGH", package = "ecp")

medians <- colMedians(ACGH$data)
mad_values <- colMads(ACGH$data)

scaled_matrix <- sweep(ACGH$data, 2, medians, FUN = "-")
scaled_matrix <- sweep(scaled_matrix, 2, mad_values, FUN = "/")
n <- nrow(scaled_matrix)
p <- ncol(scaled_matrix)

## Stock Data
data("SP500_base", package = "SP500R")

filtered_data <- as.data.frame(SP500_base) %>%
  filter(Date >= as.Date("2015-07-01") & Date <= as.Date("2016-12-31"))
matrix_data <- filtered_data %>%
  select(Date, Name, Close) %>%
  pivot_wider(names_from = Name, values_from = Close)
matrix_data <- matrix_data %>%
  select(where(~ !any(is.na(.)))) %>%
  select(-Date) %>%
  as.matrix
date_column <- filtered_data %>% select(Date) %>% unique()

medians <- colMedians(matrix_data)
mad_values <- colMads(matrix_data)

scaled_matrix <- sweep(matrix_data, 2, medians, FUN = "-")
scaled_matrix <- sweep(scaled_matrix, 2, mad_values, FUN = "/")
n <- nrow(scaled_matrix)
p <- ncol(scaled_matrix)


# Apply methods
## mxPBF
res_mxPBF <- mxPBF_combined(given_data = scaled_matrix, nws, alps = seq(1,10,0.01), a0 =  0.01, b0 = 0.01, FPR_want = 0.05, n_sample =  300, n_cores = n_cores)
cp_mxPBF <- sort(c(res_mxPBF$Change_points_cov, res_mxPBF$Change_points_mean))

## E-Divisive
res_edv <- e.divisive(scaled_matrix, R = 499, min.size = nws[1], alpha = 1)
cp_edv <- res_edv$estimates[2:(length(res_edv$estimates) - 1)]

## Inspect
cp_inspect <- inspect(t(scaled_matrix))
cp_inspect_pruned <- inspect_prune(cp_inspect$changepoints, nws[1])

## GeomCP
cpt_distance <- cpt.np(distance.mapping(scaled_matrix), penalty = "CROPS", pen.value = c(5,200),
                       method="PELT", test.stat="empirical_distribution",
                       class=TRUE, minseglen=10,
                       nquantiles =4*log(length(distance.mapping(scaled_matrix))))
cpt_angle <- cpt.np(angle.mapping(scaled_matrix), penalty = "CROPS", pen.value = c(5,200),
                    method="PELT", test.stat="empirical_distribution",
                    class=TRUE, minseglen=10,
                    nquantiles =4*log(length(distance.mapping(scaled_matrix))))
### Diagnostic plots for the CROPS algorithm (Gene data)
plot(cpt_distance, diagnostic = TRUE)
title(main = "Diagnostic Plot for cpt_distance")
sort(cpt_distance@pen.value.full, decreasing = T)
cpt_distance@pen.value.full
arrows(46, 7, 46, 12, col = "black", lwd = 2, length = 0.05)
plot(cpt_angle, diagnostic = TRUE)
title(main = "Diagnostic Plot for cpt_angle")
cpt_angle@pen.value.full
arrows(42, 8, 42, 13, col = "black", lwd = 2, length = 0.05)

cp_Geomcp <- geom_CROPS(scaled_matrix, num_distance = 46, num_angle = 42, nws[1])

### Diagnostic plots for the CROPS algorithm (Stock data)
plot(cpt_distance, diagnostic = TRUE)
title(main = "Diagnostic Plot for cpt_distance")
cpt_distance@pen.value.full
arrows(13, 7, 13, 10, col = "black", lwd = 2, length = 0.05)
plot(cpt_angle, diagnostic = TRUE)
title(main = "Diagnostic Plot for cpt_angle")
cpt_angle@pen.value.full
arrows(5, 31, 5, 35, col = "black", lwd = 2, length = 0.05)

cp_Geomcp <- geom_CROPS(scaled_matrix, num_distance = 13, num_angle = 5, nws[1])

# Analyze Results
## The number of detected change points
length(cp_mxPBF)
length(cp_Geomcp)
length(cp_inspect_pruned)
length(cp_edv)

## Calculate Covering metric
mxPBF_partition <- cps_to_partition(cp_mxPBF, dim(scaled_matrix)[1])
edv_partition <- cps_to_partition(cp_edv, dim(scaled_matrix)[1])
geomcp_partition <- cps_to_partition(cp_Geomcp, dim(scaled_matrix)[1])
inspect_partition <- cps_to_partition(cp_inspect_pruned, dim(scaled_matrix)[1])
calculate_covering_metric(mxPBF_partition, edv_partition, dim(scaled_matrix)[1])
calculate_covering_metric(mxPBF_partition, geomcp_partition, dim(scaled_matrix)[1])
calculate_covering_metric(mxPBF_partition, inspect_partition, dim(scaled_matrix)[1])

## Calculate the proportion of cps within some points of any others
calc_prop_within_range(cp_Geomcp, cp_mxPBF, 15)
calc_prop_within_range(cp_edv, cp_mxPBF, 15)
calc_prop_within_range(cp_inspect_pruned, cp_mxPBF, 15)

calc_prop_within_range(cp_mxPBF, cp_Geomcp, 15)
calc_prop_within_range(cp_mxPBF, cp_edv, 15)
calc_prop_within_range(cp_mxPBF, cp_inspect_pruned, 15)

# Visualize
## Genomic Data
### Heatmap
scaled_melted <- melt(scaled_matrix)

plot <- ggplot(scaled_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_vline(xintercept = cp_mxPBF, linetype = "dotted", color = "black") +
  labs(title = "Detected Change Point (mxPBF)", x = "Positions", y = "Individuals", fill = "Intensity") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10)
  )
ggsave(
  filename = "real_gene.png",
  plot = plot,
  width = 10,
  height = 4,
  dpi = 600,
  device = "png"
)

## Stock Data
rownames(scaled_matrix) <- as.Date(unlist(date_column))
scaled_melted <- melt(scaled_matrix)
scaled_melted$Var1 <- as.factor(scaled_melted$Var1)

breaks_to_show <- unlist(date_column)[seq(1, length(unlist(date_column)), length.out = 5)]
labels_to_show <- as.character(format(as.Date(unlist(date_column))[seq(1, length(unlist(date_column)), length.out = 5)], "%Y-%m"))

plot_s <- ggplot(scaled_melted, aes(x = Var1, y = Var2, fill = value)) +  # Var1은 factor 형식
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_vline(xintercept = as.factor(as.numeric(date_column[cp_mxPBF, ])), linetype = "dotted", color = "black") +
  labs(title = "Detected Change Point (mxPBF)", x = "Dates", y = "Companies", fill = "Return") +
  scale_x_discrete(breaks = breaks_to_show, labels = labels_to_show) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(),
    legend.text = element_text(size = 10)
  )
ggsave(
  filename = "real_stock.png",
  plot = plot_s,
  width = 10,
  height = 4,
  dpi = 600,
  device = "png"
)

date_column[cp_mxPBF,]
