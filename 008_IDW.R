#' @title Spatial Interpolation, Cross-Validation, and Mapping of Genetic Diversity (Ho)
#' @description Process genetic data, aggregate by SDM cells, optionally apply spatial thinning,
#' optimize IDW power, evaluate via LOOCV, and generate high-resolution maps across scenarios.
#' @author Chrystian Sosa (Modified via AI collaborator)
#' @date 2026-11-10

# 1. Load Required Libraries ----------------------------------------------
library(dplyr)
library(terra)
library(spatstat)
library(spatstat.explore)
library(sf)
library(readxl)
library(Metrics)
library(ggplot2)
library(ggpmisc)
library(tmap)

# 2. Define Function ------------------------------------------------------
run_idw_analysis <- function(
    sp_name,
    data_path,
    sdm_path,
    outdir,
    suffix = "_IDW",
    remove_outliers = FALSE,
    iqr_mult = 1.5,
    thinning = FALSE,
    thin_dist_km = 0
) {
  message(paste(">>> Starting IDW analysis for:", sp_name))
  
  # Adjust suffix automatically based on scenario parameters
  if (isTRUE(remove_outliers)) {
    suffix <- paste0(suffix, "_NOOUT")
  } else {
    suffix <- paste0(suffix, "_OUT")
  }
  
  if (isTRUE(thinning) && thin_dist_km > 0) {
    suffix <- paste0(suffix, sprintf("_THIN%dkm", thin_dist_km))
  } else {
    suffix <- paste0(suffix, "_NOTHIN")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Load SDM Raster
  sdm <- terra::rast(sdm_path) * 1
  sdm[sdm == 0] <- NA
  
  # Load, Filter, and Aggregate Data by SDM Raster Cell
  raw_data <- readxl::read_xlsx(data_path, sheet = "data") %>%
    dplyr::filter(sp == sp_name, !is.na(Ho), Ho > 0, !is.na(lat), !is.na(lon))
  
  coords_matrix <- as.matrix(raw_data[, c("lon", "lat")])
  raw_data$cell_id <- terra::cellFromXY(sdm, coords_matrix)
  raw_data <- raw_data %>% dplyr::filter(!is.na(cell_id))
  
  my_sf_object <- sf::st_as_sf(raw_data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  data_aggregated <- raw_data %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      Ho = mean(Ho, na.rm = TRUE),
      Ho_sd = sd(Ho, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  data_aggregated$Ho_sd[is.na(data_aggregated$Ho_sd)] <- 0
  
  cell_coords <- terra::xyFromCell(sdm, data_aggregated$cell_id)
  data_aggregated$lon <- cell_coords[, 1]
  data_aggregated$lat <- cell_coords[, 2]
  
  data_aggregated_sf <- sf::st_as_sf(data_aggregated, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  # ---> FILTRADO DE OUTLIERS POR IQR <---
  if (isTRUE(remove_outliers)) {
    q   <- stats::quantile(data_aggregated_sf$Ho, c(0.25, 0.75), na.rm = TRUE)
    iqr <- q[[2]] - q[[1]]
    
    if (iqr == 0) {
      warning("IQR = 0: outlier filter skipped.")
    } else {
      lo     <- q[[1]] - iqr_mult * iqr
      hi     <- q[[2]] + iqr_mult * iqr
      is_out <- data_aggregated_sf$Ho < lo | data_aggregated_sf$Ho > hi
      
      message(sprintf("IQR filter (k = %s): removing %d of %d cells (Ho outside [%.3f, %.3f])",
                      iqr_mult, sum(is_out), length(is_out), lo, hi))
      
      if (any(is_out)) {
        write.csv(sf::st_drop_geometry(data_aggregated_sf[is_out, ]),
                  file.path(outdir, paste0(sp_name, "_Ho_outliers_removed.csv")),
                  row.names = FALSE)
      }
      data_aggregated_sf <- data_aggregated_sf[!is_out, ]
      data_aggregated    <- data_aggregated[!is_out, ]
    }
  }
  
  # ---> SPATIAL THINNING (DISTANCIA MÍNIMA) <---
  if (isTRUE(thinning) && thin_dist_km > 0) {
    message(sprintf(">>> Applying spatial thinning (min distance: %g km)...", thin_dist_km))
    
    lon_center <- mean(data_aggregated$lon, na.rm = TRUE)
    lat_center <- mean(data_aggregated$lat, na.rm = TRUE)
    utm_zone   <- floor((lon_center + 180) / 6) + 1
    crs_utm    <- sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs", utm_zone)
    
    sf_utm <- sf::st_transform(data_aggregated_sf, crs = crs_utm)
    dist_matrix <- sf::st_distance(sf_utm)
    units(dist_matrix) <- "km"
    
    n_pts <- nrow(sf_utm)
    keep <- rep(TRUE, n_pts)
    
    for (i in seq_len(n_pts)) {
      if (keep[i]) {
        too_close <- which(as.numeric(dist_matrix[i, ]) < thin_dist_km & seq_len(n_pts) > i)
        keep[too_close] <- FALSE
      }
    }
    
    n_removed <- sum(!keep)
    message(sprintf("Thinning result: retained %d of %d cells (removed %d cells within %g km)",
                    sum(keep), n_pts, n_removed, thin_dist_km))
    
    data_aggregated_sf <- data_aggregated_sf[keep, ]
    data_aggregated    <- data_aggregated[keep, ]
  } else {
    message(">>> Spatial thinning bypassed: using all aggregated cells (conceptually reference).")
  }
  
  # Prepare Spatial Objects and Window
  points_vect <- terra::vect(data_aggregated_sf)
  coords <- terra::crds(points_vect)
  sdm_ext <- terra::ext(sdm)
  
  win <- spatstat.geom::owin(
    xrange = c(sdm_ext[1], sdm_ext[2]), 
    yrange = c(sdm_ext[3], sdm_ext[4])
  )
  
  ppp_obj <- spatstat.geom::ppp(x = coords[, 1], y = coords[, 2], window = win)
  spatstat.geom::marks(ppp_obj) <- data_aggregated$Ho
  
  # Cross-Validation to Optimize IDW Power
  powers <- seq(0.001, 10, 0.01)
  mse_result <- NULL
  
  for (power in powers) {
    CV_idw <- spatstat.explore::idw(ppp_obj, power = power, at = "points", leaveoneout = TRUE)
    mse_result <- c(mse_result, Metrics::mse(ppp_obj$marks, CV_idw))
  }
  
  optimal_power <- powers[which.min(mse_result)]
  
  # Evaluation metrics (LOOCV predictions)
  cv_preds <- spatstat.explore::idw(ppp_obj, power = optimal_power, 
                                    at = "points", leaveoneout = TRUE)
  rmse_val <- Metrics::rmse(ppp_obj$marks, as.numeric(cv_preds))
  mae_val  <- Metrics::mae(ppp_obj$marks, as.numeric(cv_preds))
  
  # Print key metrics directly to console
  message(sprintf("--- Results for %s ---", sp_name))
  message(sprintf("Optimal Power: %.4f", optimal_power))
  message(sprintf("LOOCV RMSE:    %.4f", rmse_val))
  
  # Performance Plot & Calculation of R2
  train_preds <- as.numeric(cv_preds)
  train_df <- data.frame(Observed = ppp_obj$marks, Predicted = train_preds)
  n_samples <- nrow(data_aggregated)
  
  r2_val <- summary(lm(Predicted ~ Observed, data = train_df))$r.squared
  
  file_to_save_plot <- file.path(outdir, paste0(sp_name, "_IDW_Performance_CELLS", suffix, ".png"))
  
  p_train <- ggplot(train_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "darkgreen") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "darkgreen") +
    labs(
      title = paste("A) IDW -", sp_name),
      subtitle = sprintf("n = %d cells | Opt. Power: %.3f | RMSE: %.4f", n_samples, optimal_power, rmse_val),
      x = "Observed Ho", y = "Predicted Ho"
    ) +
    theme_bw(13)
  
  ggsave(file_to_save_plot, p_train, width = 8, height = 5.5, dpi = 600)
  
  # Metrics Data Frame
  metrics_df <- data.frame(
    species = sp_name,
    n_cells = n_samples,
    optimal_power = optimal_power,
    loocv_rmse = rmse_val,
    loocv_mae = mae_val,
    r2 = r2_val,
    remove_outliers = remove_outliers,
    iqr_mult = if (remove_outliers) iqr_mult else NA,
    thinning = thinning,
    thin_dist_km = if (thinning) thin_dist_km else 0,
    date_processed = Sys.Date()
  )
  
  # Save metrics CSV
  csv_output_path <- file.path(outdir, paste0(sp_name, "_IDW_Metrics", suffix, ".csv"))
  write.csv(metrics_df, csv_output_path, row.names = FALSE)
  
  # Generate Final Raster Output
  idw_raster <- terra::rast(
    spatstat.explore::idw(ppp_obj, power = optimal_power, at = "pixels"),
    crs = "EPSG:4326"
  )
  
  idw_raster <- terra::resample(idw_raster, sdm)
  idw_raster <- terra::mask(idw_raster, sdm)
  names(idw_raster) <- "Predicted_Ho"
  
  # Save Raster as GeoTiff
  raster_output_path <- file.path(outdir, paste0(sp_name, "_Ho_raster_CELLS", suffix, ".tif"))
  writeRaster(idw_raster, raster_output_path, overwrite = TRUE)
  
  # Map Generation with tmap
  data("World", package = "tmap")
  
  file_to_save_pdf <- file.path(outdir, paste0(sp_name, "_Ho_Map_IDW_CELLS", suffix, ".pdf"))
  file_to_save_rdata <- file.path(outdir, paste0(sp_name, "_Ho_Workspace_CELLS", suffix, ".RData"))
  
  map_out <- tm_shape(idw_raster) + 
    tm_raster(
      col.legend = tm_legend(title = "Predicted Ho"),
      palette = RColorBrewer::brewer.pal(7, "YlGnBu")
    ) +
    tm_shape(World) + 
    tm_borders(col = "grey40") +
    tm_shape(my_sf_object) +
    tm_symbols(col = "grey20", size = 0.15, alpha = 0.6) +
    tm_shape(data_aggregated_sf) +
    tm_symbols(col = "red", size = 0.35, shape = 21, border.col = "black", border.lwd = 0.5) +
    tm_graticules(labels.size = 0.7) +
    tm_layout(inner.margins = 0, legend.outside = TRUE, legend.outside.position = "right")
  
  tmap::tmap_save(tm = map_out, filename = file_to_save_pdf, width = 25, height = 15, units = "cm", dpi = 600)
  
  # Save Workspace
  save.image(file_to_save_rdata)
  message(paste(">>> Successfully completed analysis for:", sp_name))
  
  return(metrics_df)
}

# ==============================================================================
# CONFIGURACIÓN GENERAL Y GRILLA DE ESCENARIOS IDW
# ==============================================================================
out_base_dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics"
data_path    <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx"
sdm_base     <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval"

species_info <- list(
  CA = list(
    name = "Crocodylus acutus",
    sdm  = file.path(sdm_base, "Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif")
  ),
  CM = list(
    name = "Crocodylus moreletii",
    sdm  = file.path(sdm_base, "Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif")
  )
)

# Grilla de combinaciones a testear:
# Incluye thin_dist_km = 0 para el escenario conceptualmente estándar (sin thinning).
idw_scenarios <- expand.grid(
  sp_code         = c("CA", "CM"),
  remove_outliers = c(TRUE, FALSE),
  thin_dist_km    = c(0, 3, 5, 10), # 0 = Sin thinning (Referencia conceptual)
  stringsAsFactors = FALSE
)

cat("Total de escenarios IDW a ejecutar:", nrow(idw_scenarios), "\n")

# ==============================================================================
# EJECUCIÓN AUTOMATIZADA
# ==============================================================================
all_idw_results <- list()

for (i in seq_len(nrow(idw_scenarios))) {
  row     <- idw_scenarios[i, ]
  sp_data <- species_info[[row$sp_code]]
  
  out_tag <- if (row$remove_outliers) "NOOUT" else "OUT"
  thin_tag <- if (row$thin_dist_km > 0) paste0("THIN", row$thin_dist_km, "KM") else "NOTHIN"
  
  folder_name <- paste("IDW", row$sp_code, thin_tag, out_tag, sep = "_")
  outdir_path <- file.path(out_base_dir, folder_name)
  
  cat(sprintf("\n------------------------------------------------------------\n"))
  cat(sprintf("[%d/%d] Ejecutando: %s\n", i, nrow(idw_scenarios), folder_name))
  cat(sprintf("------------------------------------------------------------\n"))
  
  all_idw_results[[folder_name]] <- run_idw_analysis(
    sp_name         = sp_data$name,
    data_path       = data_path,
    sdm_path        = sp_data$sdm,
    outdir          = outdir_path,
    remove_outliers = row$remove_outliers,
    iqr_mult        = 1.5,
    thinning        = (row$thin_dist_km > 0),
    thin_dist_km    = row$thin_dist_km
  )
  
  terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
}

# Consolidar todas las métricas en un único data.frame descriptivo
summary_idw_metrics <- bind_rows(all_idw_results)
write.csv(summary_idw_metrics, file.path(out_base_dir, "IDW_All_Scenarios_Summary_Metrics.csv"), row.names = FALSE)

cat("\n¡Procesamiento IDW finalizado exitosamente para todos los escenarios!\n")