# =============================================================================
# Puglia Rainfall Validation – Automated Multi-Domain Panel Maps
# Handles 5-Panel Stats (2x3) and 4-Panel Error Metrics (2x2) Layouts
# Products: PC (stations), CHIRPS v2, CHIRPS v3, ERA5, IMERG
#
# IMPROVEMENT:
#   Adaptive colour stretching per metric so maps with narrow spatial gradients
#   are still readable while keeping comparable panel scales.
# =============================================================================

# --- 1. Load Packages --------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(terra); library(gstat); library(ggplot2)
  library(ggspatial); library(haven); library(dplyr); library(RColorBrewer)
  library(scales); library(tidyterra); library(cowplot); library(grid); library(stringr)
})

# --- 2. Configuration & Paths ------------------------------------------------
base_shp    <- "G:/My Drive/Dottorato UNIBA/Puglia rainfall validation paper/Rawdata/stations/shapefiles"
data_path   <- "G:/My Drive/Dottorato UNIBA/Puglia rainfall validation paper/Temp/rain_allmetrics.dta"
output_path <- "G:/My Drive/Dottorato UNIBA/Puglia rainfall validation paper/Output/Interpolation maps"

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

target_crs <- 4326   # WGS 84 (Display)
proj_crs   <- 32633  # UTM 33N (Math/Metric)

# --- 3. Define Domains -------------------------------------------------------
# layout "5p" = 2x3 grid (PC + CHIRPS v2 + CHIRPS v3 + ERA5 + IMERG)
# layout "4p" = 2x2 grid (CHIRPS v2 + CHIRPS v3 + ERA5 + IMERG)
#
# stretch =
#   "robust"    -> clip colour limits to quantiles (better local contrast)
#   "symmetric" -> diverging scale centred on 0 with symmetric limits
#   "bounded01" -> fixed [0,100] (percent-like metrics)
#
# q_clip controls robust clipping percentiles.
map_configs <- list(
  # --- BLOCK A: RAINFALL STATISTICS (5 Panels: PC, C2, C3, E5, IMERG) ---
  list(prefix = "mean_season_2019",      title = "Mean Daily Rainfall",       unit = "mm/day", palette = "Blues",   layout = "5p", stretch = "robust", q_clip = c(0.02, 0.98)),
  list(prefix = "median_season_2019",    title = "Median Daily Rainfall",     unit = "mm/day", palette = "Blues",   layout = "5p", stretch = "robust", q_clip = c(0.02, 0.98)),
  list(prefix = "sd_season_2019",        title = "Std. Dev. Daily Rainfall",  unit = "mm",     palette = "Purples", layout = "5p", stretch = "robust", q_clip = c(0.02, 0.98)),
  list(prefix = "total_season_2019",     title = "Total Annual Rainfall",     unit = "mm",     palette = "YlGnBu",  layout = "5p", stretch = "robust", q_clip = c(0.02, 0.98)),
  list(prefix = "skew_season_2019",      title = "Skewness",                  unit = "-",      palette = "BrBG",    layout = "5p", stretch = "symmetric", q_clip = c(0.01, 0.99)),
  list(prefix = "norain_2019",           title = "Days with No Rain",         unit = "days",   palette = "OrRd",    layout = "5p", stretch = "robust", q_clip = c(0.02, 0.98)),
  list(prefix = "raindays_2019",         title = "Days with Rain",            unit = "days",   palette = "Greens",  layout = "5p", stretch = "robust", q_clip = c(0.02, 0.98)),
  list(prefix = "percent_raindays_2019", title = "Percentage of Rain Days",   unit = "%",      palette = "BuGn",    layout = "5p", stretch = "bounded01", q_clip = c(0.00, 1.00)),
  list(prefix = "dry_2019",              title = "Longest Dry Spell",         unit = "days",   palette = "YlOrRd",  layout = "5p", stretch = "robust", q_clip = c(0.02, 0.98)),

  # --- BLOCK B: VALIDATION METRICS (4 Panels: C2, C3, E5, IMERG) ---
  list(prefix = "r_rain",     title = "Pearson Correlation (r)", unit = "-",  palette = "Spectral", layout = "4p", stretch = "symmetric", q_clip = c(0.01, 0.99)),
  list(prefix = "rmse_rain",  title = "Root Mean Squared Error", unit = "mm", palette = "YlOrRd",   layout = "4p", stretch = "robust", q_clip = c(0.02, 0.98)),
  list(prefix = "cv_rain",    title = "Coefficient of Variation",unit = "-",  palette = "BuPu",     layout = "4p", stretch = "robust", q_clip = c(0.02, 0.98)),
  list(prefix = "smape_rain", title = "SMAPE",                   unit = "%",  palette = "Reds",     layout = "4p", stretch = "bounded01", q_clip = c(0.00, 1.00))
)

# Suffix & Label Definitions
suffixes_5p <- c("_p1", "_p2", "_p3", "_p4", "_p5")
labels_5p   <- c("(a) PC", "(b) CHIRPS v2", "(c) CHIRPS v3", "(d) ERA5", "(e) IMERG")

suffixes_4p <- c("_chirpsv2", "_chirpsv3", "_era5", "_imerg")
labels_4p   <- c("(a) CHIRPS v2", "(b) CHIRPS v3", "(c) ERA5", "(d) IMERG")

# --- 4. Import & Prep Spatial Data (Run Once) --------------------------------
message("Loading shapefiles and station data...")
zones <- st_read(file.path(base_shp, "zone_puglia_WGS84.shp"), quiet = TRUE) %>% st_transform(target_crs)
aoi   <- st_read(file.path(base_shp, "AOI_WGS84.shp"),         quiet = TRUE) %>% st_transform(target_crs)

zone_cols <- names(zones)
name_col  <- zone_cols[grepl("nome|name|zona", zone_cols, ignore.case = TRUE)][1]
if (is.na(name_col)) name_col <- names(zones)[1]
zones <- zones %>% rename(zone_name = all_of(name_col))

stations_raw <- read_dta(data_path) %>%
  mutate(across(where(is.labelled), as.character)) %>%
  mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# Grid Setup
aoi_utm   <- st_transform(aoi, crs = proj_crs)
bbox_utm  <- st_bbox(aoi_utm)
cell_size <- 5000
grid_sf   <- st_as_sf(
  expand.grid(
    x = seq(bbox_utm["xmin"], bbox_utm["xmax"], by = cell_size),
    y = seq(bbox_utm["ymin"], bbox_utm["ymax"], by = cell_size)
  ),
  coords = c("x", "y"), crs = proj_crs
)
grid_sf    <- grid_sf[st_intersects(grid_sf, aoi_utm, sparse = FALSE)[, 1], ]
r_template <- terra::rast(
  xmin = bbox_utm["xmin"], xmax = bbox_utm["xmax"],
  ymin = bbox_utm["ymin"], ymax = bbox_utm["ymax"],
  resolution = cell_size, crs = paste0("EPSG:", proj_crs)
)
aoi_vect <- terra::vect(aoi_utm)

bbox_aoi <- st_bbox(aoi)
pad      <- 0.05
map_lims <- c(
  x_min = bbox_aoi["xmin"] - (bbox_aoi["xmax"] - bbox_aoi["xmin"]) * pad,
  x_max = bbox_aoi["xmax"] + (bbox_aoi["xmax"] - bbox_aoi["xmin"]) * pad,
  y_min = bbox_aoi["ymin"] - (bbox_aoi["ymax"] - bbox_aoi["ymin"]) * pad,
  y_max = bbox_aoi["ymax"] + (bbox_aoi["ymax"] - bbox_aoi["ymin"]) * pad
)

# --- 5. Helper Functions -----------------------------------------------------
run_idw <- function(var_name, df) {
  st_sub <- df %>%
    select(station, latitude, longitude, val = all_of(var_name)) %>%
    filter(!is.na(val), !is.na(latitude), !is.na(longitude))

  if (nrow(st_sub) < 5) return(NULL)

  st_utm  <- st_as_sf(st_sub, coords = c("longitude", "latitude"), crs = target_crs) %>%
    st_transform(proj_crs)
  idw_obj <- gstat::idw(formula = val ~ 1, locations = st_utm, newdata = grid_sf,
                        idp = 2, debug.level = 0)
  v <- terra::vect(idw_obj)
  r <- terra::rasterize(v, r_template, field = "var1.pred", fun = "mean")
  r <- terra::mask(r, aoi_vect)
  terra::project(r, paste0("EPSG:", target_crs))
}

compute_scale <- function(values, cfg) {
  vals <- values[is.finite(values)]
  if (!length(vals)) return(list(limits = c(0, 1), breaks = waiver()))

  stretch <- cfg$stretch %||% "robust"
  q_clip  <- cfg$q_clip %||% c(0.02, 0.98)

  ql <- as.numeric(stats::quantile(vals, probs = q_clip[1], na.rm = TRUE))
  qh <- as.numeric(stats::quantile(vals, probs = q_clip[2], na.rm = TRUE))

  if (stretch == "symmetric") {
    edge <- max(abs(c(ql, qh, 0)), na.rm = TRUE)
    limits <- c(-edge, edge)
  } else if (stretch == "bounded01") {
    # for percentage-like variables, clamp to sensible metric bounds
    limits <- c(max(0, floor(min(vals, na.rm = TRUE))), min(100, ceiling(max(vals, na.rm = TRUE))))
  } else {
    limits <- c(ql, qh)
  }

  # hard constraints for known metrics
  if (cfg$prefix == "r_rain") limits <- c(max(-1, limits[1]), min(1, limits[2]))
  if (grepl("percent|smape", cfg$prefix, ignore.case = TRUE)) limits[1] <- max(0, limits[1])

  if (diff(limits) <= .Machine$double.eps) {
    eps <- abs(limits[1]) * 0.01 + 1e-6
    limits <- limits + c(-eps, eps)
  }

  # quantile-based breaks improve readability on skewed variables
  probs <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  breaks <- unique(as.numeric(stats::quantile(vals, probs = probs, na.rm = TRUE)))
  breaks <- breaks[breaks >= limits[1] & breaks <= limits[2]]

  list(limits = limits, breaks = breaks)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

create_panel <- function(raster_layer, subtitle, scale_info, pal_name, unit_label,
                         show_y = FALSE, show_x = FALSE, add_arrow = FALSE, add_scale = FALSE) {
  if (pal_name %in% row.names(brewer.pal.info)) {
    cols <- brewer.pal(min(9, brewer.pal.info[pal_name, "maxcolors"]), pal_name)
    if (pal_name == "Spectral") cols <- rev(cols)
    col_gradient <- colorRampPalette(cols)(255)
  } else {
    col_gradient <- colorRampPalette(c("white", "skyblue", "navy"))(255)
  }

  ggplot() +
    tidyterra::geom_spatraster(data = raster_layer) +
    scale_fill_gradientn(
      name = unit_label,
      colours = col_gradient,
      limits = scale_info$limits,
      breaks = if (length(scale_info$breaks) >= 2) scale_info$breaks else waiver(),
      labels = label_number(accuracy = 0.01),
      oob = scales::squish,
      na.value = NA
    ) +
    geom_sf(data = zones, fill = NA, color = "grey40", linewidth = 0.2) +
    geom_sf(data = aoi,   fill = NA, color = "#CC0000", linewidth = 0.8) +
    coord_sf(
      xlim = c(map_lims["x_min"], map_lims["x_max"]),
      ylim = c(map_lims["y_min"], map_lims["y_max"]),
      expand = FALSE
    ) +
    labs(title = subtitle, x = NULL, y = NULL) +
    theme_bw(base_family = "serif") +
    theme(
      panel.background = element_rect(fill = "#F4F8FB"),
      panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
      plot.title       = element_text(size = 10, face = "bold", margin = margin(b = 3)),
      axis.text        = element_text(size = 7, color = "grey30"),
      legend.position  = "none",
      axis.text.y = if (show_y) element_text() else element_blank(),
      axis.text.x = if (show_x) element_text() else element_blank()
    ) +
    {if (add_arrow) annotation_north_arrow(
      location = "tl", which_north = "true",
      style = north_arrow_fancy_orienteering(text_size = 6),
      pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"), height = unit(0.8, "cm")
    )} +
    {if (add_scale) annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.6)}
}

# --- 6. MAIN EXECUTION LOOP --------------------------------------------------
for (cfg in map_configs) {

  if (cfg$layout == "5p") {
    vars     <- paste0(cfg$prefix, suffixes_5p)
    labs     <- labels_5p
    msg_type <- "5-Panel Stats (2x3)"
  } else {
    vars     <- paste0(cfg$prefix, suffixes_4p)
    labs     <- labels_4p
    msg_type <- "4-Panel Metrics (2x2)"
  }

  message(paste0("Processing [", msg_type, "]: ", cfg$title))

  missing_vars <- vars[!vars %in% names(stations_raw)]
  if (length(missing_vars) > 0) {
    warning(paste("Missing columns:", paste(missing_vars, collapse = ", "), "- Skipping."))
    next
  }

  rasters <- lapply(vars, function(v) run_idw(v, stations_raw))
  names(rasters) <- vars

  all_values <- unlist(lapply(rasters, function(r) {
    if (!is.null(r)) terra::values(r, na.rm = TRUE) else NULL
  }))
  if (is.null(all_values) || length(all_values) == 0) {
    warning(paste("No valid raster values for", cfg$prefix, "- Skipping."))
    next
  }

  scale_info <- compute_scale(all_values, cfg)
  message(sprintf("  Adaptive scale: [%.3f, %.3f] %s", scale_info$limits[1], scale_info$limits[2], cfg$unit))

  if (cfg$layout == "5p") {
    p1 <- create_panel(rasters[[1]], labs[1], scale_info, cfg$palette, cfg$unit, TRUE,  FALSE, TRUE,  FALSE)
    p2 <- create_panel(rasters[[2]], labs[2], scale_info, cfg$palette, cfg$unit, FALSE, FALSE, FALSE, FALSE)
    p3 <- create_panel(rasters[[3]], labs[3], scale_info, cfg$palette, cfg$unit, FALSE, FALSE, FALSE, FALSE)
    p4 <- create_panel(rasters[[4]], labs[4], scale_info, cfg$palette, cfg$unit, TRUE,  TRUE,  FALSE, TRUE)
    p5 <- create_panel(rasters[[5]], labs[5], scale_info, cfg$palette, cfg$unit, FALSE, TRUE,  FALSE, FALSE)

    p_blank <- ggplot() + theme_void()

    row1 <- cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "hv", axis = "tblr")
    row2 <- cowplot::plot_grid(p4, p5, p_blank, nrow = 1, align = "hv", axis = "tblr")
    map_grid <- cowplot::plot_grid(row1, row2, ncol = 1, align = "v", axis = "lr")
    final_w <- 30; final_h <- 22
  } else {
    p1 <- create_panel(rasters[[1]], labs[1], scale_info, cfg$palette, cfg$unit, TRUE,  FALSE, TRUE,  FALSE)
    p2 <- create_panel(rasters[[2]], labs[2], scale_info, cfg$palette, cfg$unit, FALSE, FALSE, FALSE, FALSE)
    p3 <- create_panel(rasters[[3]], labs[3], scale_info, cfg$palette, cfg$unit, TRUE,  TRUE,  FALSE, TRUE)
    p4 <- create_panel(rasters[[4]], labs[4], scale_info, cfg$palette, cfg$unit, FALSE, TRUE,  FALSE, FALSE)

    map_grid <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv", axis = "tblr")
    final_w <- 24; final_h <- 20
  }

  p_leg <- create_panel(rasters[[1]], "Legend", scale_info, cfg$palette, cfg$unit) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.key.height = unit(1.5, "cm")
    )
  shared_legend <- cowplot::get_legend(p_leg)

  final_comp <- cowplot::plot_grid(map_grid, shared_legend, nrow = 1, rel_widths = c(1, 0.12))

  final_plot <- cowplot::plot_grid(
    ggdraw() + draw_label(cfg$title, fontface = "bold", fontfamily = "serif", size = 14, x = 0.02, hjust = 0),
    final_comp,
    ncol = 1, rel_heights = c(0.08, 1)
  )

  fname <- paste0("Panel_", cfg$prefix, ".png")
  ggsave(file.path(output_path, fname), final_plot,
         width = final_w, height = final_h, units = "cm", dpi = 300, bg = "white")
  message(paste("  Saved:", fname))
}

message("All panels complete.")
