# Generate Yautepec analysis plots from data/examples/yautepec.csv
# (notebooks/master-yautepec.qmd expects a GeoPackage that is not in this repo)

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
  library(tidyverse)
  library(ggplot2)
  library(ggsankey)
  library(rcarbon)
})

figdir <- file.path(root, "figures", "yautepec", "analysis")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 10, height = 7.5) {
  path <- file.path(figdir, paste0(name, ".png"))
  ggsave(path, p, width = width, height = height, dpi = 150, bg = "white")
  message("Wrote ", path)
  invisible(path)
}

# ---- data ----
sites <- read.csv(file.path(root, "data", "examples", "yautepec.csv"), stringsAsFactors = FALSE) |>
  select(Sitio, Area, start, end) |>
  mutate(
    Sitio = as.character(Sitio),
    start_bp = start,
    end_bp = end,
    Fase = paste0(start, "-", end)
  )

# ---- occupation / duration ----
occupation <- sites |>
  group_by(Sitio) |>
  mutate(duration = max(start, na.rm = TRUE) - min(end, na.rm = TRUE)) |>
  ungroup()
occupation <- occupation[!duplicated(occupation$Sitio), ]

png(file.path(figdir, "01-site-duration.png"), width = 10, height = 7.5, units = "in", res = 150)
d <- density(occupation$duration, bw = 150)
peaks_idx <- which(diff(sign(diff(d$y))) == -2) + 1
peaks_x <- d$x[peaks_idx]
peaks_y <- d$y[peaks_idx]
plot(d, main = "Average site duration", frame = FALSE, lwd = 2, col = "brown4",
     xlim = c(0, 3000), ylim = c(0, 1.5e-03), xlab = "Duration (years)", ylab = "Density")
points(peaks_x, peaks_y, pch = 19, col = "darkred")
text(peaks_x, peaks_y, labels = round(peaks_x), pos = 3, cex = 0.8, col = "darkred")
dev.off()
message("Wrote ", file.path(figdir, "01-site-duration.png"))
message("Median duration: ", median(occupation$duration))

# ---- foundations / abandonments ----
p_aban <- ggplot() +
  geom_density(data = occupation, aes(x = start, colour = "Site Foundations"),
               linewidth = 1, alpha = 0.25, bw = 150) +
  geom_density(data = occupation, aes(x = end, colour = "Site Abandonments"),
               linewidth = 1, alpha = 0.25, bw = 150) +
  scale_colour_manual(
    name = NULL,
    values = c("Site Foundations" = "lightblue", "Site Abandonments" = "orangered")
  ) +
  theme_minimal() +
  xlab("Years BP") + ylab("Density") +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_aban, "02-foundations-abandonments")

# ---- churn / diff ----
get_density <- function(x, from = 1, to = 3500, n = 3498) {
  dens <- density(x, from = from, to = to, n = n, na.rm = TRUE, bw = 150)
  data.frame(year = dens$x, density = dens$y)
}
df_found <- get_density(occupation$start)
df_aband <- get_density(occupation$end)
df_sum <- df_found |>
  mutate(sum = density + df_aband$density, diff = density - df_aband$density)

p_diff <- ggplot(df_sum, aes(x = year, y = sum)) +
  geom_ribbon(aes(ymin = 0, ymax = ifelse(diff > 0, diff, 0)),
              fill = "darkgreen", alpha = 0.2) +
  geom_ribbon(aes(ymin = ifelse(diff < 0, diff, 0), ymax = 0),
              fill = "firebrick", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1, colour = "lightblue3") +
  geom_line(linewidth = 1, colour = "grey55", aes(y = diff)) +
  labs(x = "Year", y = "Density") +
  geom_vline(xintercept = 430, lty = 2, col = "red") +
  scale_x_reverse(breaks = seq(3500, 0, by = -250),
                  labels = seq(-1500, 2000, by = 250)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_diff, "03-churn-diff")

# ---- total area over time ----
areas_sum <- sites |>
  group_by(start, end) |>
  summarise(Total_area = sum(Area, na.rm = TRUE), .groups = "drop")
total_sites <- sites |>
  group_by(start, end) |>
  summarise(Total_sites = n(), .groups = "drop")

p_area <- ggplot(areas_sum, aes(x = start, xend = end)) +
  geom_segment(aes(y = Total_area, yend = Total_area),
               linewidth = 3, alpha = 0.8, colour = "grey25") +
  geom_point(data = total_sites,
             aes(x = (start - (start - end) / 2), y = Total_sites * 4),
             fill = "lightblue2", size = 3, pch = 24) +
  scale_x_reverse(breaks = seq(4000, 250, by = -250),
                  labels = seq(-2000, 1750, by = 250)) +
  scale_y_continuous(
    name = "Total occupied area [ha]",
    sec.axis = sec_axis(~ . / 4, name = "Total number of sites")
  ) +
  labs(x = "Years BC/AD") +
  theme_minimal(base_size = 13) +
  geom_vline(xintercept = 430, lty = 2, col = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_area, "04-area-total")

# ---- rank-size / Zipf ----
rank_data <- sites |>
  group_by(Fase) |>
  arrange(desc(Area)) |>
  mutate(Rank = row_number(), ideal_freq = max(Area, na.rm = TRUE) / Rank) |>
  ungroup()

p_zipf <- ggplot(rank_data, aes(x = Rank, y = Area)) +
  geom_line(col = "steelblue") +
  geom_line(aes(y = ideal_freq), lty = 2, col = "grey40") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Rank (log scale)", y = "Site area (log scale)") +
  theme_minimal(base_size = 13) +
  facet_wrap(~Fase, scales = "free")
save_plot(p_zipf, "05-rank-size", width = 12, height = 10)

# ---- median / mean area ----
areas_median <- sites |>
  group_by(start, end) |>
  summarise(
    median = median(Area, na.rm = TRUE),
    mean = mean(Area, na.rm = TRUE),
    max = max(Area, na.rm = TRUE),
    q10 = quantile(Area, probs = .1, na.rm = TRUE),
    q90 = quantile(Area, probs = .9, na.rm = TRUE),
    .groups = "drop"
  )

p_med <- ggplot(areas_median) +
  geom_rect(aes(xmin = end, xmax = start, ymin = q10, ymax = q90),
            fill = "burlywood3", alpha = 0.3, colour = NA) +
  geom_segment(aes(x = start, xend = end, y = median, yend = median),
               colour = "brown4", linewidth = 1.2) +
  geom_segment(aes(x = start, xend = end, y = mean, yend = mean),
               colour = "sienna3", linewidth = 1, linetype = "dashed") +
  scale_x_reverse(breaks = seq(4000, 250, by = -250),
                  labels = seq(-2000, 1750, by = 250)) +
  scale_y_continuous(name = "Site area (ha)") +
  labs(x = "Years BC/AD", title = "Median (solid) / mean (dashed) site area") +
  theme_minimal(base_size = 13)
save_plot(p_med, "06-median-area")

p_max <- ggplot(areas_median, aes(x = (start - (start - end) / 2), y = max)) +
  geom_line() +
  geom_point(colour = "black", shape = 21, fill = "white", size = 2) +
  scale_x_reverse(breaks = seq(4000, 250, by = -250),
                  labels = seq(-2000, 1750, by = 250)) +
  scale_y_continuous(name = "Site area (ha)") +
  labs(x = "Years BC/AD", title = "Largest site per phase") +
  theme_minimal(base_size = 13)
save_plot(p_max, "07-largest-site")

# ---- urbanisation ----
urban_areas <- sites |>
  group_by(start_bp, end_bp) |>
  summarise(
    total = sum(Area, na.rm = TRUE),
    urban_40 = sum(Area[Area >= 40], na.rm = TRUE),
    urban_20 = sum(Area[Area >= 20], na.rm = TRUE),
    perc_40 = urban_40 / total * 100,
    perc_20 = urban_20 / total * 100,
    n = n(),
    .groups = "drop"
  )

p_urb <- ggplot(urban_areas) +
  geom_line(aes(x = (start_bp + end_bp) / 2, y = urban_20, colour = "≥20 ha")) +
  geom_line(aes(x = (start_bp + end_bp) / 2, y = urban_40, colour = "≥40 ha")) +
  scale_colour_manual(
    name = "Urban threshold",
    values = c("≥20 ha" = "#3182bd", "≥40 ha" = "#de2d26")
  ) +
  scale_x_reverse(
    breaks = seq(3500, 500, by = -250),
    labels = BPtoBCAD(seq(3500, 500, by = -250))
  ) +
  labs(x = "Years BC/AD", title = "Urbanisation trends") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
  geom_line(aes(x = (start_bp + end_bp) / 2, y = perc_20 * max(total) / 100),
            linetype = "dashed", colour = "#3182bd", linewidth = 1.2) +
  geom_line(aes(x = (start_bp + end_bp) / 2, y = perc_40 * max(total) / 100),
            linetype = "dashed", colour = "#de2d26", linewidth = 1.2) +
  scale_y_continuous(
    name = "Total occupied area [ha]",
    sec.axis = sec_axis(~ . / max(urban_areas$total) * 100,
                        name = "Percentage of total area (%)")
  )
save_plot(p_urb, "08-urbanisation")

# ---- trajectories grid ----
site_ids <- unique(sites$Sitio)
periods <- unique(sites[, c("start", "end")])
ext_mat <- merge(data.frame(ID = site_ids), periods, by = NULL)
sites_join <- sites[, c("Sitio", "start", "end", "Area")]
ext_df <- merge(
  ext_mat, sites_join,
  by.x = c("ID", "start", "end"),
  by.y = c("Sitio", "start", "end"),
  all.x = TRUE
)
ext_df$Area[is.na(ext_df$Area)] <- 0
ext_df$start_bp <- ext_df$start
ext_df$end_bp <- ext_df$end
ext_df$start_bc <- ext_df$start # BP values used with reverse calendar labels
ext_df <- ext_df[order(ext_df$ID, -ext_df$start), ]
ext_df <- within(ext_df, {
  prev_area <- ave(Area, ID, FUN = \(x) c(NA, head(x, -1)))
  growth_abs <- Area - prev_area
  growth_rate <- (Area - prev_area) / prev_area
})
ext_df$growth_rate[!is.finite(ext_df$growth_rate)] <- NA

teo <- ext_df[ext_df$ID %in% unique(ext_df$ID[ext_df$start == 1750 & ext_df$growth_rate >= .5]), ]
teo_large <- teo[teo$ID %in% unique(teo$ID[teo$start == 1750 & teo$Area > 20]), ]

late_large <- unique(ext_df$ID[ext_df$start < 600 & ext_df$Area > 20 & ext_df$start > 431])
early_large <- unique(ext_df$ID[ext_df$start >= 600 & ext_df$Area > 20 & ext_df$start > 431])

p_traj <- ggplot(teo[teo$start > 431, ], aes(x = start_bc, y = Area, group = ID)) +
  geom_line(alpha = 0.3) +
  geom_line(
    data = ext_df[ext_df$ID %in% early_large & ext_df$start > 431, ],
    colour = "brown4", linewidth = 1
  ) +
  geom_line(
    data = ext_df[ext_df$ID %in% late_large & ext_df$start > 431, ],
    colour = "blue4", linewidth = 1
  ) +
  scale_x_reverse(breaks = seq(3500, 500, by = -250),
                  labels = seq(-1500, 1500, by = 250)) +
  scale_y_continuous("Area [ha]") +
  theme_minimal(base_size = 13) +
  labs(x = "Years BC/AD", title = "Site trajectories (growth ≥50% at 1750 BP)")
save_plot(p_traj, "09-trajectories")

if (nrow(teo_large) > 0) {
  p_teo <- ggplot(teo[teo$start > 431, ], aes(x = start_bc, y = Area, group = ID)) +
    geom_line(alpha = 0.3) +
    geom_line(data = teo_large[teo_large$start > 431, ], colour = "brown4", linewidth = 1) +
    geom_line(
      data = teo_large[teo_large$ID %in% teo_large$ID[teo_large$end < 451 & teo_large$Area >= 20] &
                         teo_large$start > 431, ],
      colour = "blue4", linewidth = 1
    ) +
    scale_y_continuous("Area [ha]") +
    theme_minimal(base_size = 13) +
    xlab("Years BC/AD") +
    scale_x_continuous(breaks = seq(-1500, 1500, 250)) +
    labs(title = "Large Teotihuacan-period growers")
  # x was still BP; convert for continuous BC/AD scale
  teo2 <- teo |> mutate(year_bcad = BPtoBCAD(start))
  teo_l2 <- teo_large |> mutate(year_bcad = BPtoBCAD(start))
  persist <- unique(teo_l2$ID[teo_l2$end < 451 & teo_l2$Area >= 20])
  p_teo <- ggplot(teo2[teo2$start > 431, ], aes(x = year_bcad, y = Area, group = ID)) +
    geom_line(alpha = 0.3) +
    geom_line(data = teo_l2[teo_l2$start > 431, ], colour = "brown4", linewidth = 1) +
    geom_line(
      data = teo_l2[teo_l2$ID %in% persist & teo_l2$start > 431, ],
      colour = "blue4", linewidth = 1
    ) +
    scale_y_continuous("Area [ha]") +
    theme_minimal(base_size = 13) +
    xlab("Years BC/AD") +
    scale_x_continuous(breaks = seq(-1500, 1500, 250)) +
    labs(title = "Large Teotihuacan-period growers")
  save_plot(p_teo, "10-teo-large")
}

# ---- success scatter ----
p_succ <- ggplot(ext_df[ext_df$start_bp > 431, ], aes(x = prev_area, y = Area)) +
  geom_point(alpha = 0.5, colour = "grey70") +
  geom_point(
    data = ext_df[ext_df$ID %in% late_large & ext_df$start_bp > 431, ],
    aes(colour = as.factor(ID)), size = 2
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(
    x = "Previous phase area [ha]",
    y = "Current phase area [ha]",
    title = "Site area growth",
    colour = "Site ID"
  ) +
  theme_minimal(base_size = 13) +
  facet_wrap(~start_bp)
save_plot(p_succ, "11-success", width = 12, height = 10)

# ---- sankey for large late sites ----
ranked_sites <- ext_df |>
  group_by(start_bp, end_bp) |>
  arrange(desc(Area)) |>
  mutate(
    Rank = row_number(),
    Year = start_bp - (start_bp - end_bp) / 2
  ) |>
  ungroup()

df_sankey <- ranked_sites |>
  group_by(Year, Area) |>
  mutate(new_rank = min(Rank)) |>
  ungroup() |>
  group_by(ID) |>
  mutate(median_rank = median(new_rank, na.rm = TRUE)) |>
  ungroup() |>
  mutate(ID = forcats::fct_reorder(as.factor(ID), median_rank, .desc = TRUE))

large_sites <- unique(ext_df$ID[ext_df$Area > 9])
large_ids <- late_large
df_sankey$large <- ifelse(df_sankey$ID %in% large_ids, "Y", "N")

if (length(large_ids) > 0) {
  # one combined overview + up to 3 individual highlights
  overview <- ggplot(
    df_sankey[df_sankey$ID %in% large_sites & df_sankey$start_bp > 450, ],
    aes(x = Year, node = as.factor(ID), value = Area, fill = large)
  ) +
    geom_sankey_bump(space = 0, type = "alluvial", color = "transparent") +
    theme_sankey_bump(base_size = 14) +
    labs(x = "Year BP", y = "Area [ha]", title = "Sankey: sites > 9 ha") +
    scale_x_reverse() +
    scale_fill_manual(values = c("N" = "grey", "Y" = "lightblue")) +
    guides(fill = "none")
  save_plot(overview, "12-sankey-overview", width = 12, height = 8)

  for (i in seq_along(large_ids)[seq_len(min(3, length(large_ids)))]) {
    tmp_df <- df_sankey[!df_sankey$ID %in% large_ids[-i], ]
    p <- ggplot(
      tmp_df[tmp_df$ID %in% large_sites & tmp_df$start_bp > 450, ],
      aes(x = Year, node = as.factor(ID), value = Area, fill = large)
    ) +
      geom_sankey_bump(space = 0, type = "alluvial", color = "transparent") +
      theme_sankey_bump(base_size = 14) +
      labs(x = "Year BP", y = "Area [ha]", title = paste0("Site ID: ", large_ids[i])) +
      scale_x_reverse() +
      scale_fill_manual(values = c("N" = "grey", "Y" = "lightblue"), breaks = c("N", "Y")) +
      guides(fill = "none")
    save_plot(p, paste0("12-sankey-site-", large_ids[i]), width = 10, height = 7)
  }
}

# ---- summary table ----
write.csv(areas_median, file.path(figdir, "tbl-median-area.csv"), row.names = FALSE)
write.csv(urban_areas, file.path(figdir, "tbl-urban-areas.csv"), row.names = FALSE)

message("\nDone. Plots in: ", figdir)
print(list.files(figdir))
