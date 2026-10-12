#### Longevity  ----
## Use framework from https://www.pnas.org/doi/10.1073/pnas.2218834120
## It seems we can assume sites age like humans. This can be tested and the hazard/saturation can be calculated

# load necessary functions
source(file.path('scripts', 'longevity_functions.R'))
source(file.path('scripts', 'sampling.R'))

# Load required packages
library(tidyverse)
library(survival)
library(flexsurv)  # This package handles multiple parametric survival models
library(MASS)      # For fitdistr

# # load MC sim results (or run the script)
targetdir <- file.path('data', 'derived_data')

all_sites_mc <- readRDS(file.path('iter_all.rds'))
sites <- read.csv(file.path(targetdir, 'sites_processed.csv'))
sites <- sites[,-1]

# calculate persistence 
mc_pers <- mc_persistence(pnum_clean)

## Q: do large sites have longer lifespans?
# extract large sites IDs
large_sites <- sites$Sitio[sites$Area >= 10 & sites$start <= 510]

ex_large <- unique(sites$Sitio[sites$Area >= 10])
ex_large <- setdiff(ex_large, large_sites)

## add size information
mc_pers$size[mc_pers$id %in% large_sites] <- 'Large'
mc_pers$size[mc_pers$id %in% ex_large] <- 'Failed Large'
mc_pers$size[is.na(mc_pers$size)] <- 'Small'

## add survival information 
mc_pers$survival <- ifelse(
  mc_pers$norm_pers == 1, 0, 1
  # 0 = alive, 1 = dead
)

# ## Plot persistence (or longevity)
# gg_long <- ggplot(mc_pers[mc_pers$norm_pers < 1,], aes(x = persistence)) +
#   geom_histogram(aes(y = after_stat(density)), 
#                  bins = 40, fill = "grey70", colour = "white") +
#   theme_minimal() +
#   labs(
#     x = "Persistence (years)",
#     y = "Density",
#     title = "Site persistence"
#   )
# ggsave(filename = 'site_persistence_MC.tiff', gg_long, dpi = 1200, width = 7, height = 7, units = 'in')

# # # Facet by size category
# ggplot(mc_pers, aes(x = persistence, fill = size)) +
#   geom_histogram(
#     bins = 40,
#     position = "identity",
#     alpha = 0.5,
#     colour = "white"
#   ) +
#   theme_minimal() +
#   labs(
#     x = "Persistence (years)",
#     y = "Density",
#     title = "Site persistence by size class"
#   ) +
#   geom_density(alpha = 0.3, adjust = 1.2, aes(color = size)) +
#   facet_wrap(~size)
 
# ggplot(mc_pers, aes(x = persistence, colour = size, fill = size)) +
#   geom_density(alpha = 0.3, adjust = 1.2) +
#   theme_minimal()

# visually, it seems that the best fitting function could be a saturating one (try the density plots)
# this aligns with the results of pnas.2218834120
# However, we can also perform the same tests they did to find the best fitting function

# check for the best fitting function for each iteration

# in 997 out of 1000 cases the best distribution is the generalised gamma dist.
# This makes sense as many distribution used for survival analysis (e.g., Gamma, exp, Weibull)
# are special cases of this one

# clus <- parallel::makeCluster(4, 'FORK')
# best_f <- parallel::parLapply(
#   cl = clus,
#   X = 1:max(mc_pers$iteration),
#   fun = function(i) {
#     times <- mc_pers$persistence[mc_pers$iteration == i]
#     results <- compare_functions(times)
#     return(results)
#   }
# )
# parallel::stopCluster(clus)
# 
# bests <- c()
# second <- c()
# for(i in seq_along(best_f)) {
#   bests <- c(bests, best_f[[i]]$Model[1])
#   second <- c(second, best_f[[i]]$Model[2])
# }
# table(bests)
# table(second)

clus <- parallel::makeCluster(4, 'FORK')
gengammas <- parallel::parLapply(
  cl = clus,
  X = 1:max(mc_pers$iteration),
  fun = function(i) {
    
    times <- mc_pers$persistence[mc_pers$iteration == i]
    survivals <- mc_pers$survival[mc_pers$iteration == i]
    
    results <- fit_functions(times, survivals, set_func = 'gengamma')
    
    return(results)
  }
)
parallel::stopCluster(clus)
# saveRDS(gengammas, file.path(targetdir, 'gengammas.RDS')) 

gengammas <- readRDS(file.path(targetdir, 'gengammas.RDS'))
# extract gengamma values
mus <-  c()
sigmas <- c()
Qs <- c()

for(i in 1:max(mc_pers$iteration)) {
  coefs <- gengammas[[i]]$gengamma$res[,'est']
  
  mus <- c(mus, coefs[1])
  sigmas <- c(sigmas, coefs[2])
  Qs <- c(Qs, coefs[3])
}

coefs <- data.frame(
  mu = mus,
  sigma = sigmas,
  Q = Qs
)

best_c <- colMeans(coefs)

# # check it
hist(rgengamma(
  n = 1000,
  mu = best_c[1],
  sigma = best_c[2],
  Q = best_c[2]
))
# for(i in sample(1:1000, 10)){
#   hist(rgengamma(
#     n = 100,
#     mu = coefs$mu[i],
#     sigma = coefs$sigma[i],
#     Q = coefs$Q[i]
#   ), breaks = seq(0, 8000, 200))
# }

ggplot(mc_pers[mc_pers$iteration == 123,], aes(x = persistence)) +
  geom_histogram(
    bins = 40,
    alpha = 0.4,
    position = "identity",
    show.legend = FALSE
  ) +
  theme_minimal() 

# calculate density of the gengamma function
gg_dens <- dgengamma(
  1:max(mc_pers$persistence),
  mu = best_c[1],
  sigma = best_c[2],
  Q = best_c[3]
)

gg_df <- data.frame(
  persistence = 1:max(mc_pers$persistence),
  density = gg_dens
)

## Plot the mean function and the data 
ggplot(mc_pers, aes(x = persistence) ) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    alpha = 0.35,
    position = "identity"
  ) +
  # geom_line(
  #   data = gg_df,
  #   aes(x = persistence, y = density), col = 'darkred',
  #   linewidth = 1.2
  # ) +
  theme_minimal()

# sites still occupied at the last phase
still_alive <- unique(mc_pers$id[mc_pers$survival == 0 & !is.na(mc_pers$survival)])
still_alive <- still_alive[!is.na(still_alive)]

# still_alive <- unique(mc_pers$id[mc_pers$survival == 0])
# still_alive <- still_alive[!is.na(still_alive)]
# 
# mc_alive <- mc_pers[mc_pers$id %in% still_alive, ]
# 
# mc_alive$new_persistence <- sample_remaining(mc_alive$persistence, mu = best_c[1], sigma = best_c[2],
#                  Q = best_c[3])
new_pers <- mc_pers

new_pers$persistence[new_pers$id %in% still_alive]  <- sample_remaining(
  new_pers$persistence[new_pers$id %in% still_alive], 
  mu = best_c[1], sigma = best_c[2], Q = best_c[3]
  )
# write.csv(new_pers, file.path(targetdir, 'new_pers.csv'))

## you can check the sampled distribution and see if it fits the original one
# hist(sample_remaining(
#   new_pers$persistence[new_pers$id %in% still_alive],
#   mu = best_c[1], sigma = best_c[2], Q = best_c[3]
# ), breaks = seq(0,200000, 50), xlab = '', xlim = c(0, 5000), ylim = c(0,1e-03),
# freq = FALSE)

# the histogram looks more like the function, this is because the sampling was done using a unique set of values

# Now replot persistence with the new data

## Plot persistence (or longevity)
# # you can also compare it with the ones previously plotted (without sampled times)
# hist(new_pers$persistence, xlim = c(0, 2500), breaks = seq(0,20000,10))
# hist(mc_pers$persistence, xlim = c(0, 2500), breaks = seq(0,20000,10))

# I sample uniformly over cumulative probability mass.
# The inverse CDF then maps that uniform probability into time,
# stretching and compressing intervals according to the shape of the density.
# Therefore, even though I sample probabilities uniformly, the resulting event 
# times follow the exponential (or generalised gamma) distribution.

#### Study persistence ----

new_summ <- long_persistence(new_pers)

# plot persistence over time

# hist(new_summ$mean)

ggplot(new_summ) +
  geom_point(aes(x = bcad_to_bp(start, bc_to_bp = FALSE), y = mean)) +
  geom_errorbar(
    aes(
      x = bcad_to_bp(start, bc_to_bp = FALSE),
      ymin = q025,
      ymax = q975    ),
    width = 0,
    alpha = 0.4
  ) 

ggplot(new_summ) +
  geom_point(aes(x = bcad_to_bp(start, bc_to_bp = FALSE), y = mean_n)) +
  geom_errorbar(
    aes(
      x = bcad_to_bp(start, bc_to_bp = FALSE),
      ymin = q025_n,
      ymax = q975_n
    ),
    width = 0,
    alpha = 0.4
  )


new_summ$size[new_summ$id %in% large_sites] <- 'Large'
new_summ$size[new_summ$id %in% ex_large] <- 'Failed Large'

new_summ$size[is.na(new_summ$size)] <- 'Small'

## plot 
ggplot(new_summ) +
  geom_point(aes(x = bcad_to_bp(start, bc_to_bp = FALSE), y = mean, colour = size), 
             alpha = 0.6) +
  geom_smooth(aes(x = bcad_to_bp(start, bc_to_bp = FALSE), y = mean, color = size, fill = size),
              method = "lm", se = TRUE) +
  # facet_wrap(~ size) +
  theme_minimal()

#### Time - persistence correlation ----
# one would expect a decreasing correlation with time, however, given the resampling of end dates this should have been accounted for


# Add transformed x as a new column
new_summ <- new_summ |>
  mutate(start_bp = bcad_to_bp(start, bc_to_bp = FALSE))

cor_results <- new_summ |>
  group_by(size) |>
  summarise(
    cor_value = cor(start_bp, mean, use = "complete.obs"),
    p_value = cor.test(start_bp, mean)$p.value,
    .groups = "drop"
  ) |>
  mutate(
    label = paste0("r = ", round(cor_value, 2), 
                   ", p = ", signif(p_value, 3))
  )

# Step 3: Plot
sizes <- unique(new_summ$size)
cols <- c('red3', 'skyblue3', 'green4')
names(cols) <- sizes

xlim <- range(new_summ$start_bp, na.rm = TRUE)
ylim <- range(new_summ$mean, na.rm = TRUE)

plot(NA, NA, xlim = xlim, ylim = ylim, 
     xlab = "Years BC/AD", ylab = "Persistence", 
     main = "")

for (s in sizes) {
  dat <- new_summ[new_summ$size == s, ]
  
  # Points
  points(dat$start_bp, dat$mean, col = cols[s], pch = 16, cex = 0.7)
  
  # Linear regression
  fit <- lm(mean ~ start_bp, data = dat)
  abline(fit, col = cols[s], lwd = 2)
  
  # Correlation
  cor_res <- cor.test(dat$start_bp, dat$mean)
  r_val <- round(cor_res$estimate, 2)
  p_val <- signif(cor_res$p.value, 3)
  
  # Annotate in plot (adjust positions to avoid overlap)
  x_pos <- max(dat$start_bp) * 0.5
  y_pos <- 3400 * (0.9 - 0.1*which(sizes == s))
  text(x_pos, y_pos, labels = paste0("r=", r_val, ", p=", p_val), 
       col = cols[s], adj = 0)
}

# Add legend
legend("topright", legend = sizes, col = cols, pch = 16, lwd = 2)

# summarise persistence by category
new_summ |>
  group_by(size) |>
  summarise(
    mean_persistence = mean(mean, na.rm = TRUE),
    median_persistence = median(mean, na.rm = TRUE),
    q975 = quantile(mean, probs = .975, na.rm = TRUE),
    q025 = quantile(mean, probs = .025, na.rm = TRUE),
    n = n()
  )

ggplot(new_summ, 
       aes(x = -start, y = mean,  color = size)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = size)) +
  ggpubr::stat_cor(method = "pearson") +
  labs(
    x = "Distance to nearest city ≥10ha (km)",
    y = "Total site duration (years)"
  ) +
  theme_minimal() 

# kruskal.test(mean ~ size, data = new_summ)

ggplot(new_summ, aes(x = size, y = mean, fill = size)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  ylab("Persistence") +
  theme_minimal()

#### Persistence by soil category ----
# Summary statistics by phase and zone

# load soil data 
YV_sf_joined <- st_read(
  dsn = file.path(targetdir, 'yautepec.gpkg'),
  layer = 'yautepec_soil'
) |>
  st_drop_geometry() |>
  dplyr::select(Sitio, GEO_ZONE, start_bp, end_bp) |>
  group_by(Sitio) |>
  ## some sites have multiple non-contemporary phases. Split them into different sites (add a letter to ID)
  # calculate the gap. if start - previous end > 0 there is a phase jump
  arrange(start_bp) |> # if BP use desc()
  mutate(gap = abs(lag(start_bp) - end_bp)) # with BP lag(start) - end 

YV_sf_joined <- YV_sf_joined |>
  group_by(Sitio) |>
  mutate(
    occupation_period = cumsum(coalesce(gap > 0, FALSE)),
    Sitio = if_else(
      occupation_period > 0,
      paste0(Sitio, letters[occupation_period + 1]),
      as.character(Sitio)
    ),
    geo = GEO_ZONE
  ) |>
  dplyr::select(-c(occupation_period, gap)) |>  # remove helper column if you don't need it
  ungroup()

soils <- YV_sf_joined |>
  filter(!is.na(GEO_ZONE)) |>
  group_by(Sitio) |>
  summarise(
    GEO  = names(which.max(table(GEO_ZONE))),
    n    = as.numeric(sort(table(GEO_ZONE), decreasing = TRUE)[1]),
    GEO2 = names(sort(table(GEO_ZONE), decreasing = TRUE))[2],
    n2   = as.numeric(sort(table(GEO_ZONE), decreasing = TRUE)[2]),
    .groups = 'drop'
  )

soils <- left_join(
  new_summ, soils, by = c('id' = 'Sitio')
)


## Account for uncertainty by use each individual draw fro the calculations
## for each MC iter we computed the mean persistence per soil type

# join soil and MC
soil_pers <- new_pers |>
  left_join(
    soils |> dplyr::select(id, GEO),
    by = "id"
  )
soil_pers <- soil_pers[!is.na(soil_pers$GEO), ]

# mean persistence per soil per iteration
soil_means <- soil_pers |>
  group_by(iteration, GEO) |>
  summarise(m = median(persistence), .groups = "drop")

# summarise by soil
soil_summary <- soil_means |>
  group_by(GEO) |>
  summarise(
    median = median(m),
    lower  = quantile(m, 0.025),
    upper  = quantile(m, 0.975)
  ) |>
  arrange(median)

## Plot credible intervals
ypos <- seq_len(nrow(soil_summary))
soil_summary$GEO
ylabels <- c('Bas.Foot.', 'Mount.', 'Alluvium', 'Hills-Lim.', 'C.-S.Foot', 'C.Foot')

write.csv(soil_summary, file.path(targetdir, 'soil_summary.csv'))
# tiff(file.path('figures', 'soil_summary_median.tiff'), width = 9, height = 5, units = 'in', res = 600)
plot(
  soil_summary$median, ypos,
  xlim = range(c(soil_summary$lower, soil_summary$upper)),
  pch = 19,
  yaxt = "n",
  ylab = "",
  xlab = "Median persistence"
)

segments(
  soil_summary$lower[-c(2,6)], ypos[-c(2,6)],
  soil_summary$upper[-c(2,6)], ypos[-c(2,6)],
  lwd = 2
)
points(
  soil_summary$median[c(2,6)], ypos[c(2,6)],
  col = 'firebrick', pch = 19
)
segments(
  soil_summary$lower[c(2,6)], ypos[c(2,6)],
  soil_summary$upper[c(2,6)], ypos[c(2,6)],
  col = 'firebrick', lwd = 2
)

axis(2, at = ypos, labels = ylabels, las = 1)
# dev.off()

## plot boxplot of it
gg_soil <- ggplot(
  soils[!is.na(soils$GEO),],
  aes(x = GEO, y = median, fill = GEO)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  ylab("Persistence") +
  xlab("Category") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )
# ggsave(gg_soil, filename = file.path('figures', 'soil_boxplot.tiff'), width = 10, height = 6, units = 'in', dpi = 600)

## pairwise probability matrix
soil_wide <- soil_means |>
  pivot_wider(names_from = GEO, values_from = m)

soils_vec <- colnames(soil_wide)[-1]
n <- length(soils_vec)

prob_matrix <- matrix(NA, n, n,
                      dimnames = list(soils_vec, soils_vec))

for(i in 1:n){
  for(j in 1:n){
    if(i != j){
      prob_matrix[i, j] <-
        mean(soil_wide[[soils_vec[i]]] >
               soil_wide[[soils_vec[j]]])
    }
  }
}

round(prob_matrix, 3)

# write.csv(round(prob_matrix, 3), file.path(targetdir, 'prob_soil.csv'))

## linear models

# set reference category
soil_pers <- soil_pers |>
  filter(!GEO %in% c("Mountains", "Conglomerate Foothills")) |>
  mutate(GEO = relevel(factor(GEO), ref = "Alluvium"))

# fit lm per MC iteration
lm_mc <- soil_pers |>
  group_by(iteration) |>
  nest() |>
  mutate(
    model = map(data, ~ lm(norm_pers ~ GEO, data = .x)),
    tidy  = map(model, broom::tidy)
  ) |>
  unnest(tidy)

# summarise coefficient distributions
lm_summary <- lm_mc |>
  group_by(term) |>
  summarise(
    median    = median(estimate),
    lower     = quantile(estimate, 0.025),
    upper     = quantile(estimate, 0.975),
    p_positive = mean(estimate > 0)
  )
write.csv(lm_summary, file.path(targetdir, 'lm_summaty.csv'))
write.csv(lm_summary, file.path(targetdir, 'lm_summaty_norm.csv'))

# summarise by iteration
apply(
  soil_wide[,-1],
  1,
  function(x) names(which.max(x))
) |> table()

## plot
# order soils by posterior median persistence
soil_order <- soil_summary$GEO[-c(2,6)]
prob_ord   <- prob_matrix[soil_order, soil_order]

image(
  1:nrow(prob_ord),
  1:ncol(prob_ord),
  t(prob_ord[nrow(prob_ord):1, ]),
  axes = FALSE,
  xlab = "",
  ylab = "",
  col = colorRampPalette(c("blue", "white", "red"))(100)
)

axis(1, at = 1:nrow(prob_ord), labels = soil_order, las = 2)
axis(2, at = 1:ncol(prob_ord),
     labels = rev(soil_order), las = 2)

title("Posterior P(Soil_i > Soil_j)")

## check for growth by env. zone

soils <- YV_sf_joined |>
  filter(!is.na(GEO_ZONE)) |>
  group_by(Sitio) |>
  summarise(
    GEO  = names(which.max(table(GEO_ZONE))),
    n    = as.numeric(sort(table(GEO_ZONE), decreasing = TRUE)[1]),
    GEO2 = names(sort(table(GEO_ZONE), decreasing = TRUE))[2],
    n2   = as.numeric(sort(table(GEO_ZONE), decreasing = TRUE)[2]),
    .groups = 'drop'
  )

## Account for uncertainty by use each individual draw fro the calculations
## for each MC iter we computed the mean persistence per soil type

# join soil and MC
soil_gr <- left_join(
  all_sites_mc, soils, by = c('id' = 'Sitio')
)
soil_gr <- soil_gr[!is.na(soil_gr$GEO), ]

# mean area growth per soil 
# since we work on areas - which are not sampled - we can just look at a single iteration (or the original data)
soil_gr_means <- soil_gr |>
  filter(iteration == 1) |>
  group_by(GEO, start) |>
  summarise(area = sum(area), 
            .groups = "drop") |>
  group_by(GEO) |>
  arrange(desc(start), .by_group = TRUE) |>
  mutate(
    prev_area = lag(area),
    prev_start = lag(start),
    
    # time difference (positive number of years)
    dt = prev_start - start,
    da = area - prev_area,
    # CAGR
    cagr = (area / prev_area)^(1 / dt) - 1
  ) |>
  ungroup()

# summarise by soil
soil_summary <- soil_gr_means |>
  filter(!is.na(cagr)) |>
  group_by(GEO) |>
  summarise(
    mean_growth = mean(cagr),
    median = median(cagr),
    min_growth  = quantile(cagr, 0.025),
    max_growth  = quantile(cagr, 0.975),
    sd_growth   = sd(cagr),
    .groups = "drop"
  )

ggplot(soil_summary, aes(y = reorder(GEO, max_growth))) +
  geom_point(size = 3, aes(x = median)) +
  geom_errorbarh(aes(xmin = min_growth, xmax = max_growth), height = 0.2) +
  geom_point(size = 2, aes(x = mean_growth), pch = 17, color = 'red') +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  xlab('Growth Rate') +
  ylab('Soil type')

# write.csv(soil_summary, file.path(targetdir, 'soil_growth.csv'))

#### persistence vs distance from large sites/rivers ----
# spatial data
library(sf)
sp_sites <- st_read(
  dsn = file.path(targetdir, 'yautepec.gpkg'),
  layer = 'sites_processed'
)
rivers <- st_read(
  dsn = file.path(targetdir, 'yautepec.gpkg'),
  layer = 'rivers'
)
# calculate distance from rivers
sp_sites$river <- as.numeric(st_distance(sp_sites, st_union(rivers)))

riv_pers <- left_join(new_pers, sp_sites, by = join_by('id' == 'Sitio'))

# mean_riv <- riv_pers |>
#   filter(iteration == 1) |>
#   mutate(bin = cut(river, breaks = 20)) |>
#   group_by(bin) |>
#   summarise(
#     mean_pers = mean(norm_pers),
#     mean_river = mean(river)
#   )
# AIC(lm(data = mean_riv, mean_pers ~ log(mean_river)))
# AIC(lm(data = mean_riv, mean_pers ~ mean_river))
# # simple form has lower AIC. non log is preferred. 

iter_cor <- riv_pers |>
  group_by(iteration) |>
  do({
    fit <- lm(norm_pers ~ river, data = .)
    data.frame(
      intercept = coef(fit)[1],
      slope = coef(fit)[2],
      r2 = summary(fit)$r.squared
    )
  }) |>
  ungroup()

# summary(iter_cor)

# Distribution of slopes
hist(iter_cor$slope, main = "Slope distribution", xlab = "Slope")
mean(iter_cor$slope)
quantile(iter_cor$slope, c(0.025, 0.975))  # 95% interval

# Same for intercepts and R²
mean(iter_cor$r2)
quantile(iter_cor$r2, c(0.025, 0.975))

# Get the range of river values from your data
min_river <- min(riv_pers$river)
max_river <- max(riv_pers$river)

# For each iteration, compute predicted persistence at min and max
iter_cor <- iter_cor |>
  mutate(
    pred_min = intercept + slope * min_river,
    pred_max = intercept + slope * max_river,
    pred_diff = pred_max - pred_min  # difference: largest minus smallest
  )

# Summarise the predicted difference
mean(iter_cor$pred_diff)
quantile(iter_cor$pred_diff, c(0.025, 0.975))

set.seed(1234)
ggplot(
  riv_pers |> filter(iteration %in% sample(1:1000, 9)),
  aes(x = river/1000, y = norm_pers)
) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  theme_minimal() +
  xlab('Distance from rivers [km]') +
  ylab('Persistence [0-1]') +
  facet_wrap(~iteration, nrow = 3) 

## Distance to cities ----
## longevity
phase_pers <- left_join(sp_sites, new_pers, by = join_by('Sitio' == 'id'))

# filter post 510 BP
phase_pers <- phase_pers |>
  filter(start.x >= 510)

phase_pers$East <- st_coordinates(phase_pers)[,1]
phase_pers$North <- st_coordinates(phase_pers)[,2]

new_phase_pers <-  phase_pers |>
  group_by(iteration) |>
  group_modify(~ dist_func(.x)) |>
  ungroup()

# check if it worked
# test <- new_phase_pers[new_phase_pers$iteration == 1 & new_phase_pers$start.x == 800,]
# test <- st_as_sf(test)
# library(tmap)
# tm_shape(test) +
#   tm_dots(size = .1) +
#   tm_shape(test[test$Sitio %in% c('596', '164b'),]) +
#   tm_dots(fill = 'red') +
#   tm_shape(test[test$Sitio == '148',]) +
#   tm_dots(fill = 'navy')
# 
# st_distance(test[test$Sitio == '596',], test[test$Sitio == '164b',])
# 
# test$dist[test$Sitio == '164b']

dist_iter <- new_phase_pers |>
  group_by(iteration) |>
  group_modify(~ {
    fit <- lm(new_pers ~ log(dist) , data = .)
       data.frame(
         intercept = coef(fit)[1],
         slope = coef(fit)[2],
         r2 = summary(fit)$r.squared
    )
    # s <- summary(fit)
    # data.frame(
    #   coef = coef(fit)[1],
    #   hr   = s$conf.int[1, "exp(coef)"],
    #   lo95 = s$conf.int[1, "lower .95"],
    #   hi95 = s$conf.int[1, "upper .95"]
    # )
  }) |>
  ungroup()

# Distribution of slopes
hist(dist_iter$slope, main = "Slope distribution", xlab = "Slope")
mean(dist_iter$slope)
quantile(dist_iter$slope, c(0.025, 0.975))  # 95% interval

# Same for intercepts and R²
mean(dist_iter$r2)
quantile(dist_iter$r2, c(0.025, 0.975))

# Get the range of river values from your data
min_dist <- min(new_phase_pers$dist[new_phase_pers$dist > 0], na.rm = TRUE)
max_dist <- max(new_phase_pers$dist, na.rm = TRUE)

# For each iteration, compute predicted persistence at min and max
iter_cdist <- dist_iter |>
  mutate(
    pred_min = intercept + slope * log(min_dist),
    pred_max = intercept + slope * log(max_dist),
    pred_diff = pred_max - pred_min  # difference: largest minus smallest
  )

# Summarise the predicted difference
mean(iter_cdist$pred_diff)
quantile(iter_cdist$pred_diff, c(0.025, 0.975))

set.seed(1234)
ggplot(
  new_phase_pers |> filter(iteration %in% sample(1:1000, 9) ),
  aes(x = dist/1000, y = new_pers)
) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  theme_minimal() +
  xlab('Distance from cities [km]') +
  ylab('Persistence [0-1]') +
  facet_wrap(~iteration, nrow = 3) 

## Try to fit a LM with all possible vars and then do an AIC

test <- new_phase_pers 

test <- left_join(
  test, soils, by = 'Sitio'
)

test <- test |>
  filter(!is.na(dist) & !is.na(GEO))

test_lm <- test |>
  group_by(start.x, iteration) |>
  group_modify(~ {
    fit <- stepAIC(
      lm(phases ~ East + North + dist + river + Area + GEO, data = .)
    )
    data.frame(
      intercept = coef(fit)[1],
      slope = coef(fit)[2],
      r2 = summary(fit)$r.squared,
      preds = as.character(fit$call$formula)[3]
    )
  }) |>
  ungroup()

best_model <- test_lm |>
  group_by(start.x, preds) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(start.x) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup()

model_freq <- test_lm |>
  group_by(start.x, preds) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(start.x) |>
  mutate(
    total = sum(n),
    sel_prob = n / total
  ) |>
  arrange(start.x, desc(sel_prob))

top_models <- model_freq |>
  filter(sel_prob >= 0.10)

param_summary <- test_lm |>
  inner_join(top_models, by = c("start.x", "preds")) |>
  group_by(start.x, preds, sel_prob) |>
  summarise(
    slope_med = median(slope),
    slope_lo  = quantile(slope, 0.025),
    slope_hi  = quantile(slope, 0.975),
    r2_med    = median(r2),
    .groups = "drop"
  )

# cor(test$Area, test$dist)

test_lm |>
  group_by(start.x) |>
  summarise(prob_neg = mean(slope < 0))

sp_test <- st_as_sf(test)

tm_shape(sp_test[sp_test$start.x == 1100 & sp_test$iteration == 1,]) +
  tm_dots(fill = 'new_n_pers') +
  tm_shape(sp_test[sp_test$start.x == 1100 & sp_test$iteration == 1 &
                     sp_test$Area > 10,]) +
  tm_dots(fill = 'red', size = .7)

large_summary <- test[test$iteration == 1,] |>
    group_by(start.x) |>
    summarise(
      n_total = n(),
      n_large = sum(Area >= 10),
      prop_large = n_large / n_total
    )
 
mean_dist_summary <- test[test$iteration == 1,] |>
  group_by(start.x) |>
  summarise(
    mean_dist = mean(dist, na.rm = TRUE),
    median_dist = median(dist, na.rm = TRUE)
  )

## now global relation, not by period. 
# what are the more likely predictors for persistence regardless of chronology
# include start.x as variable

tt <- stepAIC(
  lm(norm_pers ~ East + North + start.x + dist + river + Area + GEO, data = test[test$iteration == 1,])
)

pers_lm <- test |>
  group_by(iteration) |>
  group_modify(~ {
    fit <- stepAIC(
      lm(norm_pers ~ East + North + start.x + dist + river + Area + GEO, data = .)
    )
    # data.frame(
    #   intercept = coef(fit)[1],
    #   slope = coef(fit)[2],
    #   slope_start = coef(fit)["start.x"],
    #   r2 = summary(fit)$r.squared,
    #   preds = as.character(fit$call$formula)[3]
    # )
  }) |>
  ungroup()

model_freq <- pers_lm |>
  count(preds) |>
  mutate(
    total = sum(n),
    sel_prob = n / total
  ) |>
  arrange(desc(sel_prob))

models_keep <- model_freq |>
  filter(sel_prob >= 0.05) |>
  pull(preds)

library(purrr)

pers_coefs <- test |>
  group_by(iteration) |>
  group_modify(~ {
    
    fit <- stepAIC(
      lm(persistence ~ East + North + start.x + dist + river + Area + GEO, data = .),
      trace = FALSE
    )
    
    tibble(
      preds = as.character(fit$call$formula)[3],
      term = names(coef(fit)),
      estimate = coef(fit)
    )
  }) |>
  ungroup()

pers_coefs_filt <- pers_coefs |>
  filter(preds %in% models_keep)

coef_summary_filt <- pers_coefs_filt |>
  group_by(term) |>
  summarise(
    median = median(estimate, na.rm = TRUE),
    lo = quantile(estimate, 0.025, na.rm = TRUE),
    hi = quantile(estimate, 0.975, na.rm = TRUE),
    prop_pos = mean(estimate > 0, na.rm = TRUE),
    prop_neg = mean(estimate < 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(abs(median)))

best_pers_model <- pers_lm |>
  group_by(preds) |>
  summarise(n = n(), .groups = "drop") 

model_freq <- pers_lm |>
  group_by(preds) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(
    total = sum(n),
    sel_prob = n / total
  ) |>
  arrange(desc(sel_prob)) |>
  filter(sel_prob >= 0.10)

model_summary <- pers_lm |>
  inner_join(model_freq, by = c("preds")) |>
  group_by(preds, sel_prob) |>
  summarise(
    slope_med = median(slope),
    slope_lo  = quantile(slope, 0.025),
    slope_hi  = quantile(slope, 0.975),
    r2_med    = median(r2),
    .groups = "drop"
  )

# cor.test(test$start.x, test$norm_pers)

### Number of phases ----
## check if number of phases is a better predictor for site longevity

# get number of phases by site
pnum <- new_phase_pers
pnum$geom <- NULL

pnum <- pnum |>
  group_by(iteration, Sitio) |>
  mutate(
    phases = n(),
    prev_phases = row_number() - 1)

# join table 
pnum <- left_join(pnum, soils, by = 'Sitio')

saveRDS(pnum, file.path(targetdir, 'pnum.RDS'))
# new_phase_pers |> 
#   group_by(iteration, start.x) |>
#   group_modify(
#     ~ {
#       fit <- lm(phases ~ new_pers, data = .)
#       data.frame(
#         intercept = coef(fit)[1],
#         slope = coef(fit)[2],
#         r2 = summary(fit)$r.squared
#       )
#     }
#   )
# 
# 
# plot(new_phase_pers$phases[new_phase_pers$iteration == 1], new_phase_pers$persistence[new_phase_pers$iteration == 1])
# 
# plot(new_phase_pers$prev_phases[new_phase_pers$iteration == 1 & new_phase_pers$prev_phases == 0],
#      new_phase_pers$persistence[new_phase_pers$iteration == 1 & new_phase_pers$prev_phases == 0])

test <- pnum[pnum$iteration == sample(1:1000, 1),]

# Summarise to one row per site first
sites_summary <- test |>
  group_by(Sitio) |>
  summarise(
    total_phases = first(phases, na_rm = TRUE),
    East  = first(East, na_rm = TRUE),
    North = first(North, na_rm = TRUE),
    dist  = first(dist, na_rm = TRUE),
    river = first(river, na_rm = TRUE),
    Area  = first(Area, na_rm = TRUE),
    GEO   = first(GEO, na_rm = TRUE),
    persistence = first(persistence) 
  )

sites_summary <- sites_summary |>
  dplyr::select(total_phases, East, North, dist, river, Area, GEO, persistence) |>
  na.omit()  # remove NAs upfront so row count stays consistent

glm.nb(total_phases ~ East + North + dist + river + Area + GEO + persistence,
       data = sites_summary) |>
  stepAIC() |>
  summary()


# Build survival structure
sites_survival <- test |>
  group_by(Sitio) |>
  arrange(Sitio, prev_phases) |>
  mutate(
    survived = as.integer(!is.na(lead(prev_phases))),  # 1 if next phase exists
    was_active = as.integer(prev_phases > 0)           # your Q2 predictor
  )

sites_survival <- sites_survival[!is.na(sites_survival$GEO),]

# scale East/North coords
sites_survival$East_s <- (sites_survival$East - min(sites_survival$East))/(max(sites_survival$East) - min(sites_survival$East))
sites_survival$North_s <- (sites_survival$North - min(sites_survival$North))/(max(sites_survival$North) - min(sites_survival$North))

m_surv <- glm(survived ~ was_active + East + North + dist + river + Area + GEO,
              data = sites_survival,
              family = binomial) |>
  stepAIC()

summary(m_surv)
exp(coef(m_surv))
# was_active ≈ 0.73 — 27% lower odds of survival if previously active

m_surv2 <- glm(survived ~ prev_phases + East + dist + Area,
               data = sites_survival,
               family = binomial)
summary(m_surv2)
# If prev_phases is negative: each additional phase reduces survival probability

# Discrete-time hazard model (logistic regression with site random effect)
library(lme4)
glmer(survived ~ was_active + East + North + dist + river + Area + (1 | Sitio),
      data = sites_survival,
      family = binomial) |>
  summary()

#### Abs area variation ----

summary_site(all_sites_mc)



#### Aligning starts ----
