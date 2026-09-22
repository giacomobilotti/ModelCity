# Generate Viabundus visualisations adapted from notebooks/analyses.qmd
# Missing upstream deps (linguistic_areas.gpkg, climate, conflict, 00_loading_data.R)
# are replaced with data/examples/viabundus.csv + approximate linguistic regions.

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
  library(sf)
  library(ggalluvial)
  library(patchwork)
  library(scales)
})

figdir <- file.path(root, "figures", "viabundus", "analysis")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 10, height = 7) {
  path <- file.path(figdir, paste0(name, ".png"))
  ggsave(path, p, width = width, height = height, dpi = 150, bg = "white")
  message("Wrote ", path)
  invisible(path)
}

# ---- load & clean ----
raw <- readr::read_csv(
  file.path(root, "data", "examples", "viabundus.csv"),
  show_col_types = FALSE
)

parse_xy <- function(g) {
  x <- as.numeric(sub("^c\\(([-0-9.]+),.*$", "\\1", g))
  y <- as.numeric(sub("^c\\([-0-9.]+,\\s*([-0-9.]+)\\)$", "\\1", g))
  cbind(x = x, y = y)
}
xy <- parse_xy(raw$geometry)

pop <- raw |>
  mutate(
    x = xy[, "x"],
    y = xy[, "y"],
    Inhabitants = as.numeric(Inhabitants),
    Year = as.integer(Year),
    Nodes_ID = as.integer(Nodes_ID)
  ) |>
  filter(!is.na(x), !is.na(y), !is.na(Inhabitants), Inhabitants > 0)

nodes_sf <- st_as_sf(pop, coords = c("x", "y"), crs = 3035, remove = FALSE)
nodes_ll <- st_transform(nodes_sf, 4326)
coords <- st_coordinates(nodes_ll)
nodes_sf$lon <- coords[, 1]
nodes_sf$lat <- coords[, 2]

# Approximate ethno-linguistic regions (stand-in for linguistic_areas.gpkg)
assign_region <- function(lon, lat) {
  dplyr::case_when(
    lon < 7.2 & lat >= 50.0 & lat < 54.5 ~ "dutch",
    lon >= 7.5 & lon < 16.5 & lat >= 54.0 ~ "danish",
    lon >= 14.5 & lat < 55.5 ~ "slavic",
    TRUE ~ "germanic"
  )
}
nodes_sf$region <- assign_region(nodes_sf$lon, nodes_sf$lat)

spatial_pop <- nodes_sf

message("Cities: ", n_distinct(spatial_pop$Name),
        " | rows: ", nrow(spatial_pop),
        " | regions: ", paste(names(table(st_drop_geometry(spatial_pop)$region)), collapse = ", "))
print(table(st_drop_geometry(spatial_pop) |> distinct(Name, region) |> pull(region)))

# ---- map: total inhabitants by year ----
europe <- tryCatch({
  rnaturalearth::ne_countries(scale = "medium", continent = "Europe", returnclass = "sf") |>
    st_transform(3035)
}, error = function(e) NULL)

years <- sort(unique(spatial_pop$Year))

if (!is.null(europe)) {
  bb <- st_bbox(st_buffer(spatial_pop, 1e5))
  europe_crop <- st_crop(europe, bb)

  map_list <- lapply(years, function(yr) {
    pts <- spatial_pop[spatial_pop$Year == yr, ]
    ggplot() +
      geom_sf(data = europe_crop, fill = "grey95", colour = "grey70", linewidth = 0.2) +
      geom_sf(
        data = pts,
        aes(size = Inhabitants),
        colour = "#1F4E79", alpha = 0.55, show.legend = "point"
      ) +
      scale_size_area(max_size = 8, name = "Inhabitants\n(×1000)") +
      labs(title = paste0("Viabundus settlements — ", yr)) +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "bottom")
  })
  p_maps <- wrap_plots(map_list, ncol = 3) +
    plot_annotation(title = "Settlement population over time (Viabundus)")
  save_plot(p_maps, "01-maps-by-year", width = 14, height = 10)
}

# region choropleth substitute: bubble map coloured by region for 1500
p_reg <- ggplot() +
  {if (!is.null(europe)) geom_sf(data = europe_crop, fill = "grey95", colour = "grey70", linewidth = 0.2)} +
  geom_sf(
    data = spatial_pop[spatial_pop$Year == 1500, ],
    aes(colour = region, size = Inhabitants), alpha = 0.7
  ) +
  scale_colour_manual(
    values = c(dutch = "#E31A1C", germanic = "#1F78B4", danish = "#33A02C", slavic = "#FF7F00")
  ) +
  scale_size_area(max_size = 10, name = "Inhabitants\n(×1000)") +
  labs(title = "Cities by approximate linguistic region (1500)",
       subtitle = "Regions approximated from geography (linguistic_areas.gpkg not in repo)",
       colour = "Region") +
  theme_void(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
save_plot(p_reg, "02-regions-1500", width = 10, height = 8)

# ---- aggregate population by region ----
agg_count <- spatial_pop |>
  st_drop_geometry() |>
  group_by(region, Year) |>
  summarise(pop = sum(Inhabitants * 1000, na.rm = TRUE), .groups = "drop")

p_agg <- ggplot(agg_count, aes(x = Year, y = pop / 1e6, colour = region)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c(dutch = "#E31A1C", germanic = "#1F78B4", danish = "#33A02C", slavic = "#FF7F00")
  ) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(title = "Total settlement population by region",
       y = "Population (millions)", x = "Year CE", colour = "Region") +
  theme_minimal(base_size = 13)
save_plot(p_agg, "03-regional-population")

# ---- relative city share (global) ----
total_pop_df <- spatial_pop |>
  st_drop_geometry() |>
  group_by(Year) |>
  summarise(Total_Settlement_Pop = sum(Inhabitants, na.rm = TRUE), .groups = "drop")

spatial_pop_df <- spatial_pop |>
  st_drop_geometry() |>
  left_join(total_pop_df, by = "Year") |>
  mutate(relative_pop = Inhabitants / Total_Settlement_Pop * 100)

# relative share by region
reg_pop_df <- spatial_pop_df |>
  group_by(region, Year) |>
  mutate(
    Total_Settlement_Pop = sum(Inhabitants, na.rm = TRUE),
    reg_relative_pop = Inhabitants / Total_Settlement_Pop * 100
  ) |>
  ungroup()

rel_plots <- lapply(sort(unique(reg_pop_df$region)), function(r) {
  ggplot(reg_pop_df[reg_pop_df$region == r, ],
         aes(x = Year, y = reg_relative_pop, group = Nodes_ID)) +
    geom_line(alpha = 0.25, colour = "grey50") +
    geom_point(alpha = 0.3, size = 0.8) +
    labs(title = paste0(r, " area"),
         y = "Share of regional pop (%)", x = "Year") +
    theme_minimal(base_size = 11)
})
save_plot(wrap_plots(rel_plots, ncol = 2), "04-city-relative-share", width = 12, height = 9)

# ---- ranks & top-hit trajectories ----
spatial_pop_clean <- spatial_pop_df |>
  filter(!is.na(relative_pop)) |>
  arrange(Nodes_ID, Year) |>
  group_by(Year) |>
  mutate(pop_rank = rank(-relative_pop, ties.method = "min")) |>
  ungroup()

top_hits <- spatial_pop_clean |>
  filter(pop_rank <= 10) |>
  distinct(Name)

p_rank <- ggplot(spatial_pop_clean, aes(x = Year, y = pop_rank, group = Name)) +
  geom_line(alpha = 0.15, colour = "grey60") +
  geom_line(
    data = filter(spatial_pop_clean, Name %in% top_hits$Name),
    aes(colour = Name), linewidth = 1
  ) +
  scale_y_reverse(breaks = seq(1, max(spatial_pop_clean$pop_rank), by = 20)) +
  labs(
    y = "Population rank (1 = largest)",
    title = "City–rank trajectories over time",
    subtitle = "Top-10 cities (any year) highlighted"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))
save_plot(p_rank, "05-rank-trajectories", width = 12, height = 8)

# ---- won / lost status ----
city_result <- spatial_pop_clean |>
  group_by(Nodes_ID) |>
  summarise(
    Name = first(Name),
    region = first(region),
    Last_Year = max(Year),
    Last_Share = relative_pop[Year == max(Year)][1],
    Max_Share = max(relative_pop),
    Max_Year = Year[which.max(relative_pop)],
    Second_Max_Share = sort(unique(relative_pop), decreasing = TRUE)[min(2, n_distinct(relative_pop))],
    .groups = "drop"
  ) |>
  mutate(
    Status = case_when(
      Max_Year == Last_Year & Max_Share >= Second_Max_Share * 1.10 ~ "Won",
      TRUE ~ "Lost"
    )
  )

p_status <- ggplot(city_result, aes(x = Second_Max_Share, y = Max_Share, colour = Status)) +
  geom_abline(slope = 1.10, intercept = 0, linetype = "dashed", colour = "darkgreen") +
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values = c(Won = "#2ECC71", Lost = "#E74C3C")) +
  labs(
    title = "City peak share vs second-best share",
    subtitle = "Won = peak at last year and ≥10% above second peak",
    x = "Second highest population share (%)",
    y = "Maximum population share (%)"
  ) +
  theme_minimal(base_size = 13)
save_plot(p_status, "06-won-lost")

# ---- regional ranks + CHURN ----
reg_pop_clean <- reg_pop_df |>
  group_by(region, Year) |>
  mutate(pop_rank = rank(-reg_relative_pop, ties.method = "min")) |>
  ungroup()

reg_ranks <- reg_pop_clean |>
  group_by(region, Nodes_ID) |>
  arrange(Year, .by_group = TRUE) |>
  mutate(
    rank_change = lag(pop_rank) - pop_rank,
    rank_normalised = rank_change / (Year - lag(Year))
  ) |>
  ungroup()

reg_change <- reg_ranks |>
  group_by(region, Year) |>
  summarise(
    total_change = sum(abs(rank_change), na.rm = TRUE),
    total_change_normal = sum(abs(rank_normalised), na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(region) |>
  mutate(scaled_change = total_change_normal / max(total_change_normal, na.rm = TRUE)) |>
  ungroup()

# global churn
pop_ranks <- spatial_pop_clean |>
  group_by(Nodes_ID) |>
  arrange(Year, .by_group = TRUE) |>
  mutate(
    rank_change = lag(pop_rank) - pop_rank,
    rank_normalised = rank_change / (Year - lag(Year))
  ) |>
  ungroup()

rank_change <- pop_ranks |>
  group_by(Year) |>
  summarise(
    total_change = sum(abs(rank_change), na.rm = TRUE),
    total_change_normal = sum(abs(rank_normalised), na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    region = "Global",
    scaled_change = total_change_normal / max(total_change_normal, na.rm = TRUE)
  )

ranks_collapsed <- bind_rows(reg_change, rank_change)

p_churn <- ggplot(ranks_collapsed, aes(x = Year, y = total_change_normal, colour = region)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c(dutch = "#E31A1C", germanic = "#1F78B4", danish = "#33A02C",
               slavic = "#FF7F00", Global = "grey40")
  ) +
  facet_wrap(~region, scales = "free_y") +
  labs(title = "Adapted CHURN: normalised rank change by region",
       y = "Normalised |Δrank| / year", x = "Year CE") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
save_plot(p_churn, "07-churn-by-region", width = 11, height = 7)

p_churn_scaled <- ggplot(ranks_collapsed, aes(x = Year, y = scaled_change, colour = region)) +
  geom_line(linewidth = 1.2, alpha = 0.9) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c(dutch = "#E31A1C", germanic = "#1F78B4", danish = "#33A02C",
               slavic = "#FF7F00", Global = "grey40")
  ) +
  labs(title = "Regional political/settlement change (scaled CHURN)",
       y = "Scaled change", x = "Year CE", colour = "Region") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")
save_plot(p_churn_scaled, "08-churn-scaled")

# ---- success / fail / stability ----
reg_ranks_classified <- reg_ranks |>
  group_by(region, Year) |>
  mutate(treshold = max(pop_rank, na.rm = TRUE) / 10) |>
  ungroup() |>
  mutate(
    Fail = ifelse(rank_change < -treshold & Year != min(Year), TRUE, FALSE),
    Success = ifelse(rank_change > treshold & Year != min(Year), TRUE, FALSE),
    Stability = ifelse(!Fail & !Success & Year != min(Year), TRUE, FALSE),
    Stability_Type = case_when(
      Stability & pop_rank <= treshold ~ "positive",
      Stability & pop_rank > treshold ~ "negative",
      TRUE ~ NA_character_
    )
  )

succ_plots <- lapply(sort(unique(reg_ranks_classified$region)), function(r) {
  reg_ranks_classified |>
    filter(region == r) |>
    group_by(Year) |>
    summarise(
      Success = mean(Success, na.rm = TRUE),
      Fail = mean(Fail, na.rm = TRUE),
      Stability_Pos = mean(Stability_Type == "positive", na.rm = TRUE),
      Stability_Neg = mean(Stability_Type == "negative", na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_longer(-Year, names_to = "Category", values_to = "Count") |>
    ggplot(aes(x = Year, y = Count, fill = Category)) +
    geom_area(alpha = 0.85) +
    scale_fill_manual(values = c(
      Success = "steelblue", Fail = "red3",
      Stability_Pos = "darkgreen", Stability_Neg = "grey60"
    )) +
    labs(title = paste0("Trajectories — ", r),
         y = "Share of cities", x = "Year", fill = "Trajectory") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
})
save_plot(wrap_plots(succ_plots, ncol = 2), "09-success-fail", width = 12, height = 9)

# alluvial of status by region
status_plots <- lapply(sort(unique(reg_ranks_classified$region)), function(r) {
  tmp <- reg_ranks_classified |>
    filter(region == r, Year != 1300) |>
    mutate(Status = case_when(
      Success ~ "Success",
      Fail ~ "Fail",
      Stability_Type == "positive" ~ "Stable+",
      Stability_Type == "negative" ~ "Stable−",
      TRUE ~ "Unclassified"
    ))
  ggplot(tmp, aes(x = Year, stratum = Status, alluvium = Name, y = 1, fill = Status)) +
    geom_flow(alpha = 0.65) +
    geom_stratum(width = 0.35) +
    scale_fill_manual(values = c(
      Success = "#2ECC71", Fail = "#E74C3C",
      `Stable+` = "#9ACD32", `Stable−` = "#3498DB",
      Unclassified = "grey80"
    )) +
    labs(title = paste0("Status flows — ", r), x = "Year", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
})
save_plot(wrap_plots(status_plots, ncol = 2), "10-status-alluvial", width = 12, height = 9)

# ---- Sankey/alluvial for top cities per region ----
reg_hits <- reg_pop_clean |>
  filter(pop_rank <= 7) |>
  distinct(region, Name)

make_city_alluvial <- function(city_name, region_name) {
  df_plot <- reg_pop_clean |>
    filter(region == region_name) |>
    mutate(
      Year_factor = factor(Year),
      Highlight = if_else(Name == city_name, city_name, "Other")
    )

  ggplot(df_plot,
         aes(x = Year_factor, stratum = reorder(Name, -pop_rank),
             alluvium = Name, y = Inhabitants, fill = Highlight)) +
    geom_alluvium(alpha = 0.75, decreasing = FALSE) +
    scale_fill_manual(values = setNames(c("lightgrey", "steelblue"), c("Other", city_name))) +
    labs(title = paste0(city_name, " (", region_name, ")"),
         x = "Year", y = "Inhabitants (×1000)") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}

# one highlight city per region (largest mean pop among top hits)
highlight_cities <- reg_pop_clean |>
  semi_join(reg_hits, by = c("region", "Name")) |>
  group_by(region, Name) |>
  summarise(mean_pop = mean(Inhabitants, na.rm = TRUE), .groups = "drop") |>
  group_by(region) |>
  slice_max(mean_pop, n = 2, with_ties = FALSE) |>
  ungroup()

sankey_plots <- pmap(
  list(highlight_cities$Name, highlight_cities$region),
  make_city_alluvial
)
save_plot(wrap_plots(sankey_plots, ncol = 2), "11-city-alluvial-highlights",
          width = 12, height = 14)

# ---- top cities bar by year ----
top10 <- spatial_pop_clean |>
  filter(pop_rank <= 10) |>
  mutate(Name = fct_reorder(Name, relative_pop, .fun = max, .desc = TRUE))

p_top <- ggplot(top10, aes(x = factor(Year), y = Inhabitants, fill = Name)) +
  geom_col(position = "dodge") +
  labs(title = "Top-10 cities by global population share, per year",
       x = "Year", y = "Inhabitants (×1000)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0)) +
  guides(fill = guide_legend(nrow = 4))
save_plot(p_top, "12-top10-bars", width = 12, height = 8)

# largest cities line chart
big <- spatial_pop_clean |>
  filter(Name %in% top_hits$Name)
p_big <- ggplot(big, aes(x = Year, y = Inhabitants * 1000, colour = Name)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Population of historically top-ranked cities",
       y = "Inhabitants", x = "Year CE", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))
save_plot(p_big, "13-top-city-population", width = 12, height = 7)

# zipf / rank-size per year
rank_size <- spatial_pop_clean |>
  group_by(Year) |>
  arrange(desc(Inhabitants), .by_group = TRUE) |>
  mutate(
    Rank = row_number(),
    ideal = max(Inhabitants) / Rank
  ) |>
  ungroup()

p_zipf <- ggplot(rank_size, aes(x = Rank, y = Inhabitants)) +
  geom_line(colour = "steelblue") +
  geom_line(aes(y = ideal), linetype = "dashed", colour = "grey40") +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~Year) +
  labs(title = "Rank–size (Zipf) by year",
       x = "Rank (log)", y = "Inhabitants ×1000 (log)") +
  theme_minimal(base_size = 12)
save_plot(p_zipf, "14-rank-size", width = 11, height = 7)

# tables
write_csv(city_result, file.path(figdir, "tbl-city-status.csv"))
write_csv(ranks_collapsed, file.path(figdir, "tbl-churn.csv"))
write_csv(agg_count, file.path(figdir, "tbl-regional-pop.csv"))

message("\nDone. Figures in: ", figdir)
print(list.files(figdir))
