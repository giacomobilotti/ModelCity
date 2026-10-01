# Run chronological Monte Carlo sampling and write growth / persistence plots
# from data/examples/yautepec.csv

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
})

old_wd <- getwd()
setwd(root)
on.exit(setwd(old_wd), add = TRUE)

# sampling.R sources its helpers from a path of its own; load them from here
# and skip that call so the runner does not depend on the working layout.
source(file.path("scripts", "r", "helpers_sampling.R"))
for (expr in parse(file.path("scripts", "r", "sampling.R"))) {
  if (!(is.call(expr) && identical(expr[[1]], as.name("source")))) eval(expr, globalenv())
}

figdir <- file.path(root, "figures", "yautepec", "sampling")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

save_base_plot <- function(expr, name, width = 10, height = 7.5) {
  path <- file.path(figdir, paste0(name, ".png"))
  png(path, width = width, height = height, units = "in", res = 150)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  message("Wrote ", path)
  invisible(path)
}

# ---- data ----
sites <- read.csv(
  file.path(root, "data", "examples", "yautepec.csv"),
  stringsAsFactors = FALSE
)
sites$Sitio <- as.character(sites$Sitio)

# Example multi-phase sites from sampling.R comments, plus longest-lived sites
focus_ids <- c("160", "200")
phase_counts <- sort(table(sites$Sitio), decreasing = TRUE)
extra_ids <- names(phase_counts)[seq_len(min(6L, length(phase_counts)))]
focus_ids <- unique(c(focus_ids, extra_ids))

focus <- sites[sites$Sitio %in% focus_ids, ]
terminus <- 510L
iterations <- 500L
set.seed(42)

message(
  "Running MC: ", length(focus_ids), " sites, ",
  iterations, " iterations, terminus=", terminus
)
mc <- run_mc(
  focus,
  iterations = iterations,
  site_id_col = "Sitio",
  area_col = "Area",
  start_method = "unif",
  end_method = "unif",
  peak_method = "unif",
  terminus = terminus
)
by_site <- summary_site(mc)
by_period <- summary_period(mc)
pers <- persistence(mc, terminus = 431)

write.csv(by_site, file.path(figdir, "summary-by-site.csv"), row.names = FALSE)
write.csv(by_period, file.path(figdir, "summary-by-period.csv"), row.names = FALSE)
write.csv(pers, file.path(figdir, "persistence-by-site.csv"), row.names = FALSE)
message("Wrote summary CSVs to ", figdir)

for (sid in focus_ids) {
  site_mc <- mc[mc$id == sid, ]
  site_sum <- by_site[by_site$id == sid, ]
  # plot_site_cagr needs >= 2 summarised rows for CI ribbons
  if (nrow(site_sum) < 2) {
    message("Skipping site ", sid, " (need >=2 phase rows for CAGR plot; got ", nrow(site_sum), ")")
    next
  }

  save_base_plot(
    plot_site_cagr(site_sum),
    paste0("cagr-site-", sid)
  )
  save_base_plot(
    plot_spaghetti(site_mc, site_sum, iter_num = 50, type = "growth_cagr"),
    paste0("spaghetti-cagr-site-", sid)
  )
  save_base_plot(
    plot_spaghetti(site_mc, site_sum, iter_num = 50, type = "area"),
    paste0("spaghetti-area-site-", sid)
  )
}

save_base_plot({
  plot(
    by_period$med_year, by_period$median,
    type = "b", pch = 16, col = "steelblue", lwd = 2,
    xlab = "Year (BC/AD)", ylab = "CAGR",
    main = "Median CAGR by period (focus sites)",
    ylim = range(c(by_period$q025, by_period$q975), na.rm = TRUE)
  )
  polygon(
    c(by_period$med_year, rev(by_period$med_year)),
    c(by_period$q025, rev(by_period$q975)),
    col = "#4682B433", border = NA
  )
  lines(by_period$med_year, by_period$median, col = "steelblue", lwd = 2)
  points(by_period$med_year, by_period$median, pch = 16, col = "steelblue")
  abline(h = 0, lty = 2, col = "grey10")
}, "cagr-by-period")

save_base_plot({
  ord <- order(pers$median, decreasing = TRUE)
  p <- pers[ord, ]
  bp <- barplot(
    p$median,
    names.arg = p$id,
    las = 2,
    col = "steelblue",
    ylab = "Median persistence (years)",
    main = paste0("Site persistence to terminus 431 BP (n=", iterations, ")")
  )
  arrows(
    bp, p$q025, bp, p$q975,
    angle = 90, code = 3, length = 0.05
  )
}, "persistence-by-site", width = 11, height = 7)

demo <- sites[sites$Sitio %in% c("160", "200"), ]
message("Running exp/norm MC on sites 160 & 200...")
mc2 <- run_mc(
  demo,
  iterations = iterations,
  terminus = terminus,
  start_method = "exp",
  end_method = "exp",
  peak_method = "norm",
  perc = 5,
  peak_sd_perc = 0.1
)
sum2 <- summary_site(mc2)
for (sid in c("160", "200")) {
  site_mc <- mc2[mc2$id == sid, ]
  site_sum <- sum2[sum2$id == sid, ]
  if (nrow(site_sum) < 2) {
    message("Skipping exp/norm site ", sid, " (got ", nrow(site_sum), " rows)")
    next
  }
  save_base_plot(
    plot_site_cagr(site_sum),
    paste0("cagr-expnorm-site-", sid)
  )
  save_base_plot(
    plot_spaghetti(site_mc, site_sum, iter_num = 50, type = "growth_cagr"),
    paste0("spaghetti-expnorm-cagr-site-", sid)
  )
}

message("Done. Figures in ", figdir)
