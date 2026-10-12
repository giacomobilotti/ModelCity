# Run the longevity analysis on data/examples/yautepec.csv.
#
# longevity.R reads Monte Carlo output and spatial layers that are not in this
# repository, and it leaves pnum_clean undefined. This runner builds that
# object with run_mc(), then calls mc_persistence() and the functions in
# longevity_functions.R. It writes the plots from the part of longevity.R that
# does not need yautepec.gpkg. Soil, river, and distance-to-city figures are
# omitted because those layers are not in the checkout.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg)) {
  script_file <- normalizePath(sub("^--file=", "", file_arg[[1]]))
  root <- normalizePath(file.path(dirname(script_file), "..", ".."))
} else {
  root <- normalizePath(getwd())
}

lib <- file.path(root, "R_libs")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(flexsurv)
  library(survival)
})

old_wd <- getwd()
setwd(root)
on.exit(setwd(old_wd), add = TRUE)

source(file.path("scripts", "r", "helpers_sampling.R"))
for (expr in parse(file.path("scripts", "r", "sampling.R"))) {
  if (!(is.call(expr) && identical(expr[[1]], as.name("source")))) eval(expr, globalenv())
}
source(file.path("scripts", "r", "longevity_functions.R"))

figdir <- file.path(root, "figures", "yautepec", "longevity")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

save_ggplot <- function(plot, name, width = 8, height = 6) {
  path <- file.path(figdir, paste0(name, ".png"))
  ggsave(path, plot, width = width, height = height, dpi = 150, bg = "white")
  message("Wrote ", path)
  invisible(path)
}

save_base_plot <- function(expr, name, width = 8, height = 6) {
  path <- file.path(figdir, paste0(name, ".png"))
  png(path, width = width, height = height, units = "in", res = 150)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  message("Wrote ", path)
  invisible(path)
}

sites <- read.csv(
  file.path(root, "data", "examples", "yautepec.csv"),
  stringsAsFactors = FALSE
)
sites$Sitio <- as.character(sites$Sitio)

iterations <- 1000L
set.seed(42)

message(
  "Running MC: ", length(unique(sites$Sitio)),
  " sites, ", iterations, " iterations"
)
mc <- run_mc(
  sites,
  iterations = iterations,
  site_id_col = "Sitio",
  area_col = "Area",
  start_method = "unif",
  end_method = "unif",
  peak_method = "unif",
  terminus = 510
)

# The object longevity.R calls pnum_clean. Terminus 431 BP is the Yautepec default.
mc_pers <- mc_persistence(mc, terminus = 431)

large_sites <- sites$Sitio[sites$Area >= 10 & sites$start <= 510]
ex_large <- setdiff(unique(sites$Sitio[sites$Area >= 10]), large_sites)

mc_pers$size <- NA_character_
mc_pers$size[mc_pers$id %in% large_sites] <- "Large"
mc_pers$size[mc_pers$id %in% ex_large] <- "Failed Large"
mc_pers$size[is.na(mc_pers$size)] <- "Small"
mc_pers$survival <- ifelse(mc_pers$norm_pers == 1, 0, 1)

message(
  "Persistence rows: ", nrow(mc_pers),
  "; censored draws: ", sum(mc_pers$survival == 0, na.rm = TRUE),
  "; observed ends: ", sum(mc_pers$survival == 1, na.rm = TRUE)
)

message("Fitting a generalised gamma to each iteration")
fit_log <- file.path(figdir, "gengamma-fit.log")
sink(fit_log)
gengammas <- vector("list", iterations)
for (i in seq_len(iterations)) {
  gengammas[[i]] <- fit_functions(
    mc_pers$persistence[mc_pers$iteration == i],
    mc_pers$survival[mc_pers$iteration == i],
    set_func = "gengamma"
  )
  if (i %% 100 == 0) message("Fitted ", i, "/", iterations)
}
sink()

failed <- which(vapply(gengammas, function(fit) is.null(fit$gengamma), logical(1)))
if (length(failed)) {
  stop("Generalised gamma failed for iterations: ", paste(failed, collapse = ", "))
}

coefs <- t(vapply(gengammas, function(fit) fit$gengamma$res[, "est"], numeric(3)))
colnames(coefs) <- c("mu", "sigma", "Q")
coefs <- as.data.frame(coefs)
best_c <- colMeans(coefs)
# sample_remaining() reads best_c from the global environment.
assign("best_c", best_c, envir = globalenv())

message(
  "Mean generalised-gamma coefficients: mu=", signif(best_c[["mu"]], 4),
  " sigma=", signif(best_c[["sigma"]], 4),
  " Q=", signif(best_c[["Q"]], 4)
)

max_persistence <- max(mc_pers$persistence)
density_grid <- data.frame(
  persistence = seq_len(max_persistence),
  density = dgengamma(
    seq_len(max_persistence),
    mu = best_c[["mu"]],
    sigma = best_c[["sigma"]],
    Q = best_c[["Q"]]
  )
)

save_base_plot(
  {
    hist(
      rgengamma(1000, mu = best_c[["mu"]], sigma = best_c[["sigma"]], Q = best_c[["Q"]]),
      breaks = 40,
      col = "grey80",
      border = "white",
      main = "Draws from the mean generalised gamma",
      xlab = "Persistence (years)"
    )
  },
  "gengamma-parameter-draw"
)

save_ggplot(
  ggplot(mc_pers[mc_pers$iteration == 123, ], aes(x = persistence)) +
    geom_histogram(bins = 40, alpha = 0.4, position = "identity", show.legend = FALSE) +
    theme_minimal() +
    labs(
      title = "Site persistence, iteration 123",
      x = "Persistence (years)",
      y = "Count"
    ),
  "persistence-iteration-123"
)

save_ggplot(
  ggplot(mc_pers, aes(x = persistence)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 40,
      alpha = 0.35,
      position = "identity"
    ) +
    theme_minimal() +
    labs(
      title = "Site persistence across draws",
      x = "Persistence (years)",
      y = "Density"
    ),
  "persistence-histogram"
)

save_ggplot(
  ggplot(mc_pers, aes(x = persistence)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 40,
      alpha = 0.35,
      position = "identity"
    ) +
    geom_line(
      data = density_grid,
      aes(x = persistence, y = density),
      colour = "darkred",
      linewidth = 1.2
    ) +
    theme_minimal() +
    labs(
      title = "Observed persistence and the mean generalised gamma",
      x = "Persistence (years)",
      y = "Density"
    ),
  "persistence-gengamma-overlay"
)

still_alive <- unique(mc_pers$id[mc_pers$survival == 0 & !is.na(mc_pers$survival)])
still_alive <- still_alive[!is.na(still_alive)]

new_pers <- mc_pers
imputed <- sample_remaining(
  new_pers$persistence[new_pers$id %in% still_alive],
  mu = best_c[["mu"]],
  sigma = best_c[["sigma"]],
  Q = best_c[["Q"]]
)
n_bad <- sum(!is.finite(imputed))
if (n_bad) {
  message("Dropping ", n_bad, " non-finite imputed lifespans")
  imputed[!is.finite(imputed)] <- NA_real_
}
new_pers$persistence[new_pers$id %in% still_alive] <- imputed

new_summ <- long_persistence(new_pers)
new_summ$size <- NA_character_
new_summ$size[new_summ$id %in% large_sites] <- "Large"
new_summ$size[new_summ$id %in% ex_large] <- "Failed Large"
new_summ$size[is.na(new_summ$size)] <- "Small"
new_summ$size <- factor(new_summ$size, levels = c("Large", "Failed Large", "Small"))
new_summ$start_bp <- bcad_to_bp(new_summ$start, bc_to_bp = FALSE)

write.csv(coefs, file.path(figdir, "gengamma-coefficients.csv"), row.names = FALSE)
write.csv(new_summ, file.path(figdir, "persistence-by-site.csv"), row.names = FALSE)
write.csv(
  new_summ |>
    group_by(size) |>
    summarise(
      mean_persistence = mean(mean, na.rm = TRUE),
      median_persistence = median(mean, na.rm = TRUE),
      q975 = quantile(mean, probs = 0.975, na.rm = TRUE),
      q025 = quantile(mean, probs = 0.025, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ),
  file.path(figdir, "persistence-by-size.csv"),
  row.names = FALSE
)

save_ggplot(
  ggplot(new_summ, aes(x = start_bp, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, alpha = 0.4) +
    theme_minimal() +
    labs(
      title = "Imputed persistence against foundation date",
      x = "Median foundation (BC negative, AD positive)",
      y = "Mean persistence (years)"
    ),
  "persistence-vs-foundation"
)

save_ggplot(
  ggplot(new_summ, aes(x = start_bp, y = mean_n)) +
    geom_point() +
    geom_errorbar(aes(ymin = q025_n, ymax = q975_n), width = 0, alpha = 0.4) +
    theme_minimal() +
    labs(
      title = "Normalised persistence against foundation date",
      subtitle = "Normalised values are from before lifespan imputation",
      x = "Median foundation (BC negative, AD positive)",
      y = "Mean share of available time"
    ),
  "norm-persistence-vs-foundation"
)

save_ggplot(
  ggplot(new_summ, aes(x = start_bp, y = mean, colour = size, fill = size)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_minimal() +
    labs(
      title = "Imputed persistence by size class",
      x = "Median foundation (BC negative, AD positive)",
      y = "Mean persistence (years)",
      colour = "Size",
      fill = "Size"
    ),
  "persistence-vs-foundation-by-size"
)

save_base_plot(
  {
    sizes <- levels(new_summ$size)
    cols <- c("Large" = "red3", "Failed Large" = "skyblue3", "Small" = "green4")
    xlim <- range(new_summ$start_bp, na.rm = TRUE)
    ylim <- range(new_summ$mean, na.rm = TRUE)
    plot(
      NA, NA, xlim = xlim, ylim = ylim,
      xlab = "Foundation (BC negative, AD positive)",
      ylab = "Mean persistence (years)",
      main = "Persistence and foundation date"
    )
    for (s in sizes) {
      dat <- new_summ[new_summ$size == s, ]
      points(dat$start_bp, dat$mean, col = cols[s], pch = 16, cex = 0.7)
      fit <- lm(mean ~ start_bp, data = dat)
      abline(fit, col = cols[s], lwd = 2)
      cor_res <- cor.test(dat$start_bp, dat$mean)
      y_pos <- ylim[2] - (ylim[2] - ylim[1]) * (0.06 + 0.08 * (match(s, sizes) - 1))
      text(
        xlim[1], y_pos,
        labels = paste0(s, ": r=", round(cor_res$estimate, 2), ", p=", signif(cor_res$p.value, 3)),
        col = cols[s], adj = 0
      )
    }
    legend("topright", legend = sizes, col = cols[sizes], pch = 16, lwd = 2)
  },
  "persistence-correlation-by-size"
)

save_ggplot(
  ggplot(new_summ, aes(x = size, y = mean, fill = size)) +
    geom_boxplot(alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4) +
    theme_minimal() +
    labs(
      title = "Imputed persistence by size class",
      x = NULL,
      y = "Mean persistence (years)",
      fill = "Size"
    ),
  "persistence-boxplot-by-size"
)

message("Longevity plots written to ", figdir)
