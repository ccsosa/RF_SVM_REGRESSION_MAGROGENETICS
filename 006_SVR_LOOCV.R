# ==============================================================================
# DEPENDENCIES
# ==============================================================================
library(dplyr)
library(readxl)
library(sf)
library(terra)
library(caret)
library(doParallel)
library(kernlab)
library(yardstick)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(raster)
library(dismo)
library(usdm)
library(RColorBrewer)
library(tmap)
library(blockCV)
library(fastshap)   # Shapley values (>= 0.1.0: no plotting methods)
library(shapviz)    # SHAP visualisation (beeswarm additionally needs ggbeeswarm)

#' Spatial Prediction Workflow for Genetic Heterozygosity (Ho) using Support
#' Vector Regression (cell level; VIF- or correlation-based decorrelation)
#'
#' @description
#' Fits an SVR model (via \pkg{caret}) that relates observed heterozygosity (Ho)
#' to environmental predictors at the raster-cell level, evaluates it with
#' repeated k-fold CV, LOOCV or spatial block CV, explains it with SHAP values
#' (\pkg{fastshap} + \pkg{shapviz}), and projects Ho over the species' binary
#' SDM, masking non-analogue conditions with MESS.
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Load predictor rasters (must include layers named \code{lon} and
#'         \code{lat}) and genetic records for \code{sp_name}
#'         (sheet \code{"data"}; columns \code{sp}, \code{Ho}, \code{lon},
#'         \code{lat}). Records with missing or zero Ho are dropped.
#'   \item Assign each record to a raster cell (off-coverage points are snapped
#'         to the nearest valid cell within 5 km) and aggregate per cell:
#'         median Ho, MAD and number of samples. Optionally remove Ho outliers
#'         with the IQR rule.
#'   \item Extract predictor values per cell (coastal NAs are filled from a
#'         focal-mean raster) and, if more than 49 cells are available, thin
#'         cells spatially with \code{terra::thin()} at \code{thin_dist} m.
#'   \item Drop predictors dominated by a single value (\code{dom_max}), then
#'         decorrelate with VIF or Pearson correlation.
#'   \item Train the SVR (cell weights = \code{sqrt(n_samples_in_cell)}; inputs
#'         centred and scaled) over a tuning grid, minimising RMSE.
#'   \item Compute train (resubstitution) and cross-validated metrics.
#'   \item Compute SHAP values on the final model and export importance,
#'         beeswarm and dependence plots.
#'   \item Predict Ho over the SDM, mask extrapolation with MESS and export
#'         rasters, maps (PDF) and the R workspace.
#' }
#'
#' @param outdir \code{character}. Output directory path (created if missing).
#' @param sp_name \code{character}. Exact target species name (matched against
#'   column \code{sp}).
#' @param raster_dir \code{character}. Directory with predictor rasters
#'   (\code{.tif}); must contain layers named \code{lon} and \code{lat}.
#' @param data_path \code{character}. Excel file with genetic data (sheet
#'   \code{"data"}).
#' @param sdm_path \code{character}. Binary SDM raster (1 = suitable, 0 = not).
#' @param n_cores \code{integer}. Cores for model training. Default \code{6L}.
#' @param decorr_method \code{character}. \code{"vif"} (\code{usdm::vifstep})
#'   or \code{"cor"} (\code{caret::findCorrelation}). Default \code{"vif"}.
#' @param vif_th \code{numeric}. VIF threshold (only if
#'   \code{decorr_method = "vif"}). Default \code{10}.
#' @param cor_cutoff \code{numeric}. Absolute Pearson correlation threshold
#'   (only if \code{decorr_method = "cor"}). Default \code{0.5}.
#' @param dom_max \code{numeric}. Maximum fraction of cells allowed to share the
#'   same value of a predictor; predictors above it are dropped. Default
#'   \code{0.3}.
#' @param addLonLat \code{logical}. Retain \code{lon}/\code{lat} as predictors.
#'   Default \code{FALSE}.
#' @param use_loocv \code{logical}. Use leave-one-out CV (ignored if
#'   \code{use_blockcv = TRUE}). Default \code{FALSE}; if both are
#'   \code{FALSE}, 5-fold CV repeated 10 times is used.
#' @param use_blockcv \code{logical}. Use spatial block CV (k = 4, via
#'   \pkg{blockCV}); requires at least 49 cells, otherwise the scenario is
#'   discarded, \code{outdir} is deleted and \code{NULL} is returned. Default
#'   \code{FALSE}.
#' @param remove_outliers \code{logical}. Remove Ho outliers at the cell level
#'   (IQR rule); removed cells are written to a CSV. Default \code{FALSE}.
#' @param iqr_mult \code{numeric}. IQR multiplier for the outlier fences.
#'   Default \code{1.5}.
#' @param kernel_type \code{character}. \code{"svmRadial"} or
#'   \code{"svmLinear"}. Default \code{"svmRadial"}.
#' @param thin_dist \code{numeric}. Minimum spatial thinning distance in
#'   metres. Default \code{3000}.
#'
#' @return Invisibly returns \code{NULL}. Results are written to
#'   \code{outdir}: data used for modelling, predictor list, dominance and VIF
#'   tables, metrics and CV predictions (CSV), performance and SHAP plots (PNG),
#'   SHAP values (CSV), raw and MESS-masked Ho rasters (GeoTIFF), maps (PDF) and
#'   the R workspace (\code{.RData}).
#'
#' @export
svr_ho_pipeline <- function(outdir,
                            sp_name,
                            raster_dir,
                            data_path,
                            sdm_path,
                            n_cores         = 6L,
                            decorr_method   = c("vif", "cor"),
                            vif_th          = 10,
                            cor_cutoff      = 0.5,
                            dom_max         = 0.3,
                            addLonLat       = FALSE,
                            use_loocv       = FALSE,
                            use_blockcv     = FALSE,
                            remove_outliers = FALSE,
                            iqr_mult        = 1.5,
                            kernel_type     = "svmRadial",
                            thin_dist       = 3000) {
  
  set.seed(123)
  
  decorr_method <- match.arg(decorr_method)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  if (isTRUE(addLonLat)) message("Using lon and lat as predictors")
  
  # ------------------------------------------------------------------------------
  # 1. RASTER & GENETIC DATA LOADING
  # ------------------------------------------------------------------------------
  bios_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
  bios_names <- sub("\\.tif$", "", list.files(raster_dir, pattern = "\\.tif$"))
  
  bios <- terra::rast(bios_files)
  names(bios) <- bios_names
  
  if (!all(c("lon", "lat") %in% names(bios))) {
    stop("Error: 'lon' and/or 'lat' layers are missing from the raster directory.")
  }
  
  data <- readxl::read_xlsx(data_path, sheet = "data") %>%
    dplyr::filter(sp == sp_name, !is.na(Ho), Ho > 0, !is.na(lat), !is.na(lon))
  
  my_sf_object <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  points_vect  <- terra::vect(my_sf_object)
  
  # ------------------------------------------------------------------------------
  # 2. CELL ASSIGNMENT & SPATIAL AGGREGATION
  # ------------------------------------------------------------------------------
  cell_ids <- terra::cells(bios[[1]], points_vect)[, "cell"]
  na_cells <- which(is.na(cell_ids))
  
  if (length(na_cells) > 0) {
    message(sprintf("Assigning %d off-coverage points to nearest valid raster cell...", length(na_cells)))
    extracted_cells <- terra::extract(
      bios[[1]],
      points_vect[na_cells, ],
      buffer = 5000,
      fun = "first",
      cells = TRUE,
      exact = FALSE
    )
    cell_ids[na_cells] <- extracted_cells$cell
  }
  
  my_sf_object$cell_id  <- cell_ids
  coords_cells          <- terra::xyFromCell(bios[[1]], cell_ids)
  my_sf_object$cell_lon <- coords_cells[, 1]
  my_sf_object$cell_lat <- coords_cells[, 2]
  my_sf_object          <- my_sf_object[!is.na(my_sf_object$cell_id), ]
  
  data_aggregated_sf <- my_sf_object %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      n_samples_in_cell = dplyr::n(),
      Ho_mad = mad(Ho, na.rm = TRUE),
      Ho     = median(Ho, na.rm = TRUE),
      lon    = unique(cell_lon),
      lat    = unique(cell_lat),
      .groups = "drop"
    )
  
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
    }
  }
  
  message("Modelling with ", nrow(data_aggregated_sf), " aggregated grid cells.")
  
  bios_focal <- terra::focal(bios, w = 15, fun = mean, na.policy = "only", na.rm = TRUE)
  gc()
  
  unique_points_vect <- terra::vect(data_aggregated_sf)
  ext_direct         <- terra::extract(bios, unique_points_vect, "mean")
  
  if (any(!complete.cases(ext_direct))) {
    message("Coastal points with NA values detected. Extracting values from nearest valid pixel...")
    na_rows <- which(!complete.cases(ext_direct))
    for (i in na_rows) {
      ext_nearest <- terra::extract(bios_focal, unique_points_vect[i, ], nearest = TRUE)
      ext_direct[i, names(bios)] <- ext_nearest[, names(bios)]
    }
  }
  
  ext_direct_clean <- ext_direct %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::select(-ID)
  
  data_aggregated_sp <- cbind(data_aggregated_sf, ext_direct_clean)
  rm(ext_direct_clean); gc()
  
  if (nrow(data_aggregated_sp) > 49) {
    message(sprintf("More than 49 cells available, terra::thin at %d m", thin_dist))
    
    xy  <- terra::xyFromCell(bios[[1]], data_aggregated_sp$cell_id)    # one point per cell
    pts <- terra::vect(xy, type = "points", crs = "EPSG:4326")
    pts$cell_id <- data_aggregated_sp$cell_id
    
    kept <- terra::thin(pts, d = thin_dist)                            # d in metres
    
    data_aggregated_sp <- data_aggregated_sp[data_aggregated_sp$cell_id %in% kept$cell_id, ]
    
    D <- sf::st_distance(data_aggregated_sp)
    D <- matrix(as.numeric(D), nrow = nrow(D), ncol = ncol(D))
    diag(D) <- NA
    
    message(sprintf("Retained %d of %d cells | min. distance = %.1f km",
                    nrow(kept), nrow(pts),
                    min(D, na.rm = TRUE) / 1000)
    )
  }
  
  # ------------------------------------------------------------------------------
  # 3. FEATURE SELECTION & DECORRELATION
  # ------------------------------------------------------------------------------
  data_sel <- data_aggregated_sp %>%
    dplyr::select(Ho, n_samples_in_cell, dplyr::all_of(names(bios))) %>%
    sf::st_drop_geometry()
  
  message(sprintf("Retained %d cells | min. distance = %.1f km",
                  nrow(data_aggregated_sp),
                  min(terra::distance(terra::vect(data_aggregated_sp))) / 1000))
  
  valid_rows <- complete.cases(data_sel)
  data_sel   <- data_sel[valid_rows, ]
  n_samples  <- nrow(data_sel)
  
  bio_vars_only <- setdiff(names(bios), c("lon", "lat"))
  bio_df        <- as.data.frame(data_sel[, bio_vars_only])
  
  # Dominance screen: drop predictors where one value covers > dom_max of cells
  diag_var <- data.frame(
    var      = names(bio_df),
    dom_frac = sapply(bio_df, function(x) max(table(round(x, 6))) / length(x))
  )
  diag_var$drop <- diag_var$dom_frac > dom_max
  
  write.csv(diag_var, file.path(outdir, paste0(sp_name, "_dominance_screen.csv")),
            row.names = FALSE)
  message("Dropped by dominance (one value in more than ", dom_max * 100, "% of cells): ",
          paste(diag_var$var[diag_var$drop], collapse = ", "))
  
  bio_vars_only <- diag_var$var[!diag_var$drop]
  bio_df        <- bio_df[, bio_vars_only, drop = FALSE]
  
  if (decorr_method == "vif") {
    message(sprintf("Decorrelation: VIF stepwise (th = %s)", vif_th))
    vif_res       <- usdm::vifstep(bio_df, th = vif_th)
    selected_bios <- setdiff(bio_vars_only, vif_res@excluded)
    write.csv(
      as.data.frame(vif_res@results),
      file.path(outdir, paste0(sp_name, "_VIF_results.csv")),
      row.names = FALSE
    )
  } else {
    message(sprintf("Decorrelation: Pearson |r| cutoff = %s", cor_cutoff))
    cor_matrix    <- cor(bio_df)
    to_remove     <- caret::findCorrelation(cor_matrix, cutoff = cor_cutoff)
    selected_bios <- if (length(to_remove) > 0) bio_vars_only[-to_remove] else bio_vars_only
  }
  
  message(sprintf("Selected %d/%d predictors: %s",
                  length(selected_bios), length(bio_vars_only),
                  paste(selected_bios, collapse = ", ")))
  
  predictors_list <- if (isTRUE(addLonLat)) c("lon", "lat", selected_bios) else selected_bios
  
  data_sel_model <- as.data.frame(data_sel[, c("Ho", "n_samples_in_cell", predictors_list)])
  
  write.csv(data_sel_model,
            file.path(outdir, paste0(sp_name, "_data_used_for_model_CELLS.csv")),
            row.names = FALSE)
  write.csv(predictors_list,
            file.path(outdir, paste0(sp_name, "_predictor_CELLS.csv")),
            row.names = FALSE)
  
  # ------------------------------------------------------------------------------
  # 4. SVR MODEL TRAINING & VALIDATION SCHEME
  # ------------------------------------------------------------------------------
  weights_vec    <- sqrt(data_sel_model$n_samples_in_cell)
  data_sel_model <- data_sel_model %>% dplyr::select(-n_samples_in_cell)
  
  applied_blockcv <- FALSE
  
  if (isTRUE(use_blockcv)) {
    if (n_samples >= 49) {
      message(sprintf("Validation scheme: Spatial Block CV (blockCV) [n = %d >= 49]", n_samples))
      
      data_sf_model <- data_aggregated_sp[valid_rows, ]
      
      set.seed(123)
      sb <- blockCV::cv_spatial(x = data_sf_model, k = 4, selection = "random", progress = FALSE)
      
      train_index <- lapply(sb$folds_list, function(fold) fold[[1]])
      test_index  <- lapply(sb$folds_list, function(fold) fold[[2]])
      names(train_index) <- paste0("Fold", seq_along(train_index))
      names(test_index)  <- names(train_index)
      
      train_ctrl <- caret::trainControl(
        method          = "cv",
        number          = length(train_index),
        index           = train_index,
        indexOut        = test_index,
        savePredictions = "final",
        allowParallel   = TRUE
      )
      applied_blockcv <- TRUE
      
    } else {
      message(sprintf("[SCENARIO DISCARDED] 'use_blockcv = TRUE' requested, but n = %d (< 49).", n_samples))
      message("--> Removing empty/incomplete output folder and ending this scenario...\n")
      
      if (dir.exists(outdir)) {
        unlink(outdir, recursive = TRUE)
      }
      return(NULL)
    }
  }
  
  if (!applied_blockcv) {
    if (isTRUE(use_loocv)) {
      message(sprintf("Validation scheme: LOOCV [n = %d]", n_samples))
      train_ctrl <- caret::trainControl(method = "LOOCV", savePredictions = "final", allowParallel = TRUE)
    } else {
      message(sprintf("Validation scheme: Repeated 5-Fold CV (10 repeats) [n = %d]", n_samples))
      train_ctrl <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 10,
                                        savePredictions = "final", allowParallel = TRUE)
    }
  }
  
  val_scheme_str <- if (applied_blockcv) "BLOCKCV" else if (isTRUE(use_loocv)) "LOOCV" else "REPEATEDCV"
  
  # ------------------------------------------------------------------------------
  # SVM TRAINING
  # ------------------------------------------------------------------------------
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit({ try(stopCluster(cl), silent = TRUE); registerDoSEQ() }, add = TRUE)
  
  set.seed(123)
  
  if (kernel_type == "svmRadial") {
    svm_grid <- expand.grid(
      sigma = c(0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1),
      C     = c(0.05, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 500, 1000)
    )
  } else {
    svm_grid <- expand.grid(C = 10^seq(-3, 2, length.out = 30))
  }
  
  final_svr <- caret::train(
    Ho ~ .,
    data       = data_sel_model,
    method     = kernel_type,
    weights    = weights_vec,
    preProcess = c("center", "scale"),
    trControl  = train_ctrl,
    tuneGrid   = svm_grid,
    metric     = "RMSE"
  )
  
  stopCluster(cl)
  registerDoSEQ()
  
  # ------------------------------------------------------------------------------
  # 5. MODEL METRICS & PREDICTION CONSOLIDATION
  # ------------------------------------------------------------------------------
  train_pred_vals <- predict(final_svr, newdata = data_sel_model)
  train_df <- data.frame(
    rowIndex  = 1:nrow(data_sel_model),
    Observed  = data_sel_model$Ho,
    Predicted = train_pred_vals
  )
  
  train_rmse <- yardstick::rmse_vec(train_df$Observed, train_df$Predicted)
  train_mae  <- yardstick::mae_vec(train_df$Observed, train_df$Predicted)
  train_rsq  <- yardstick::rsq_vec(train_df$Observed, train_df$Predicted)
  
  cv_preds_raw <- dplyr::inner_join(
    final_svr$pred,
    final_svr$bestTune,
    by = names(final_svr$bestTune)
  )
  
  cv_df <- cv_preds_raw %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      Observed  = mean(obs,  na.rm = TRUE),
      Predicted = mean(pred, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(rowIndex)
  
  test_rmse <- yardstick::rmse_vec(cv_df$Observed, cv_df$Predicted)
  test_mae  <- yardstick::mae_vec(cv_df$Observed, cv_df$Predicted)
  test_rsq  <- yardstick::rsq_vec(cv_df$Observed, cv_df$Predicted)
  
  metrics_out <- data.frame(
    .metric           = c("rmse", "mae", "rsq", "rmse", "mae", "rsq"),
    mean              = c(test_rmse, test_mae, test_rsq, train_rmse, train_mae, train_rsq),
    dataset           = c("validation_test", "validation_test", "validation_test", "train", "train", "train"),
    validation_scheme = val_scheme_str,
    kernel_type       = kernel_type,
    add_lon_lat       = addLonLat,
    decorr_method     = decorr_method,
    thin_dist_m       = thin_dist,
    n_predictors      = length(predictors_list),
    remove_outliers   = remove_outliers
  )
  
  suffix <- paste0("_", kernel_type, "_", val_scheme_str, "_", toupper(decorr_method),
                   "_THIN", thin_dist, "m",
                   if (addLonLat) "_LONLAT" else "",
                   if (isTRUE(remove_outliers)) "_NOOUT" else "")
  
  write.csv(metrics_out,
            file.path(outdir, paste0(sp_name, "_SVM_Metrics_CELLS", suffix, ".csv")),
            row.names = FALSE)
  
  write.csv(cv_df,
            file.path(outdir, paste0(sp_name, "_CV_Predictions_Test_CELLS", suffix, ".csv")),
            row.names = FALSE)
  
  # ------------------------------------------------------------------------------
  # 6. EVALUATIVE PLOTS
  # ------------------------------------------------------------------------------
  p_train <- ggplot(train_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "darkgreen") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "darkgreen") +
    labs(title = "A) Training Fit",
         subtitle = paste("n =", n_samples, "cells (Resubstitution)"),
         x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  test_subtitle <- if (applied_blockcv) {
    "Spatial Block Cross-Validation (blockCV)"
  } else if (use_loocv) {
    "Leave-One-Out CV (LOOCV)"
  } else {
    "Repeated 5-Fold CV (Averaged)"
  }
  
  p_test <- ggplot(cv_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "blue") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "blue") +
    labs(title = "B) Cross-Validation (Test)",
         subtitle = test_subtitle,
         x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  ggsave(
    file.path(outdir, paste0(sp_name, "_SVM_Performance_CELLS", suffix, ".png")),
    p_train + p_test, width = 12, height = 5.5, dpi = 600
  )
  
  # ------------------------------------------------------------------------------
  # 6.1. EXPLAINABILITY WITH SHAP VALUES (fastshap + shapviz)
  # ------------------------------------------------------------------------------
  message("Calculating SHAP values for the SVR model...")
  
  # Prediction wrapper for fastshap
  pfun <- function(object, newdata) {
    as.numeric(predict(object, newdata = newdata))
  }
  
  # 'X' must contain only the predictors (no response 'Ho')
  X_shap <- data_sel_model[, predictors_list, drop = FALSE]
  
  # parallel = FALSE: the cluster is already closed at this point, and a
  # sequential foreach backend would emit one warning per call otherwise.
  set.seed(123)
  shap_values <- fastshap::explain(
    final_svr,
    X            = X_shap,
    pred_wrapper = pfun,
    nsim         = 50,        # Monte Carlo repetitions (increase for more precision)
    parallel     = FALSE
  )
  
  # Cell-level SHAP values
  shap_df <- as.data.frame(shap_values)
  shap_df$cell_row <- seq_len(nrow(shap_df))
  write.csv(
    shap_df,
    file.path(outdir, paste0(sp_name, "_SHAP_Values_CELLS", suffix, ".csv")),
    row.names = FALSE
  )
  
  # Mean |SHAP| per predictor
  shap_imp <- sort(colMeans(abs(as.matrix(shap_values))), decreasing = TRUE)
  write.csv(
    data.frame(predictor = names(shap_imp), mean_abs_shap = as.numeric(shap_imp)),
    file.path(outdir, paste0(sp_name, "_SHAP_Importance", suffix, ".csv")),
    row.names = FALSE
  )
  
  # shapviz object (baseline = mean prediction used by fastshap)
  # fastshap stores NA (not NULL) when adjust = FALSE; fall back to the mean
  # prediction over X, which is the reference value of the Shapley decomposition
  shap_baseline <- attr(shap_values, "baseline")
  if (is.null(shap_baseline) || length(shap_baseline) != 1L || !is.finite(shap_baseline)) {
    shap_baseline <- mean(pfun(final_svr, X_shap), na.rm = TRUE)
  }
  sv <- shapviz::shapviz(as.matrix(shap_values), X = X_shap, baseline = shap_baseline)
  
  # A) Global importance (bar)
  p_shap_imp <- shapviz::sv_importance(sv, kind = "bar", max_display = Inf) +
    labs(
      title    = "A) Global Importance (SHAP)",
      subtitle = paste("SVR Model |", val_scheme_str),
      x        = "Mean Absolute SHAP Value (Impact on Ho)"
    ) +
    theme_bw(13)
  
  ggsave(
    file.path(outdir, paste0(sp_name, "_SHAP_Importance", suffix, ".png")),
    p_shap_imp, width = 7, height = 6, dpi = 600
  )
  
  # B) Beeswarm (needs the 'ggbeeswarm' package)
  if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
    p_shap_bee <- shapviz::sv_importance(sv, kind = "beeswarm", max_display = Inf) +
      labs(
        title    = "B) SHAP Beeswarm",
        subtitle = paste("SVR Model |", val_scheme_str)
      ) +
      theme_bw(13)
    
    ggsave(
      file.path(outdir, paste0(sp_name, "_SHAP_Beeswarm", suffix, ".png")),
      p_shap_bee, width = 8, height = 6, dpi = 600
    )
  } else {
    message("Package 'ggbeeswarm' not installed: beeswarm plot skipped.")
  }
  
  # C) Dependence plot for the most important predictor
  top_var <- names(shap_imp)[1]
  p_dep <- shapviz::sv_dependence(sv, v = top_var) +
    labs(title = paste("C) SHAP Dependence:", top_var)) +
    theme_bw(13)
  
  ggsave(
    file.path(outdir, paste0(sp_name, "_SHAP_Dependence_", top_var, suffix, ".png")),
    p_dep, width = 6.5, height = 5, dpi = 600
  )
  
  # ------------------------------------------------------------------------------
  # 7. SPATIAL PREDICTION (RASTER)
  # ------------------------------------------------------------------------------
  sdm <- terra::rast(sdm_path) * 1
  sdm[sdm == 0] <- NA
  
  bios_selected <- bios[[predictors_list]]
  bios_selected <- terra::resample(bios_selected, sdm, method = "bilinear")
  bios_selected <- bios_selected * sdm
  
  valid_cells <- terra::cells(bios_selected[[1]])
  
  temp_dt <- as.data.frame(terra::extract(bios_selected, valid_cells))
  if ("ID" %in% colnames(temp_dt)) temp_dt$ID <- NULL
  temp_dt$cell_idx <- valid_cells
  
  valid_rows_pred <- complete.cases(temp_dt[, predictors_list, drop = FALSE])
  temp_dt         <- temp_dt[valid_rows_pred, ]
  temp_dt$prediction <- NA_real_
  
  chunk_size <- 50000
  n_rows     <- nrow(temp_dt)
  n_chunks   <- ceiling(n_rows / chunk_size)
  
  for (i in seq_len(n_chunks)) {
    s_idx <- ((i - 1) * chunk_size) + 1
    e_idx <- min(i * chunk_size, n_rows)
    newdata_chunk <- temp_dt[s_idx:e_idx, predictors_list, drop = FALSE]
    temp_dt$prediction[s_idx:e_idx] <- predict(final_svr, newdata = newdata_chunk)
  }
  
  pred_raster <- bios_selected[[1]] * NA
  names(pred_raster) <- "Predicted_Ho"
  pred_raster[temp_dt$cell_idx] <- temp_dt$prediction
  
  terra::writeRaster(
    pred_raster,
    file.path(outdir, paste0(sp_name, "_Ho_Map_SVM_CELLS_RAW", suffix, ".tif")),
    overwrite = TRUE
  )
  
  # ------------------------------------------------------------------------------
  # 8. MESS MASKING
  # ------------------------------------------------------------------------------
  bios_selected_stack <- raster::stack(bios_selected)
  ref_data <- as.data.frame(data_sel_model[, predictors_list, drop = FALSE])
  
  mess_rn <- terra::rast(dismo::mess(x = bios_selected_stack, v = ref_data))
  mess_rn <- terra::resample(mess_rn, sdm, method = "bilinear")
  
  mess_mask <- mess_rn
  mess_mask[mess_mask <= 0 | is.infinite(mess_mask)] <- NA
  mess_mask[!is.na(mess_mask)] <- 1
  
  pred_raster_no_interpolated <- pred_raster * mess_mask
  
  terra::writeRaster(
    pred_raster_no_interpolated,
    file.path(outdir, paste0(sp_name, "_Ho_Map_SVM_CELLS", suffix, ".tif")),
    overwrite = TRUE
  )
  
  # ------------------------------------------------------------------------------
  # 9. MAP EXPORT & WORKSPACE SAVING
  # ------------------------------------------------------------------------------
  data("World", package = "tmap")
  
  map_out <- tm_shape(pred_raster) +
    tm_raster(
      col.legend = tm_legend(title = "Ho"),
      palette = RColorBrewer::brewer.pal(7, "YlGnBu")
    ) +
    tm_shape(World) +
    tm_borders(col = "grey40") +
    tm_shape(my_sf_object) +
    tm_symbols(col = "grey20", size = 0.15, alpha = 0.6) +
    tm_shape(data_aggregated_sp) +
    tm_symbols(col = "red", size = 0.15, shape = 21, border.col = "black", border.lwd = 0.5) +
    tm_graticules(labels.size = 0.7) +
    tm_layout(inner.margins = 0, legend.outside = TRUE, legend.outside.position = "right")
  
  tmap::tmap_save(
    tm       = map_out,
    filename = file.path(outdir, paste0(sp_name, "_Ho_Map_SVM_CELLS_RAW", suffix, ".pdf")),
    width = 25, height = 15, units = "cm", dpi = 600
  )
  
  map_out <- tm_shape(pred_raster_no_interpolated) +
    tm_raster(
      col.legend = tm_legend(title = "Ho"),
      palette = RColorBrewer::brewer.pal(7, "YlGnBu")
    ) +
    tm_shape(World) +
    tm_borders(col = "grey40") +
    tm_shape(my_sf_object) +
    tm_symbols(col = "grey20", size = 0.15, alpha = 0.6) +
    tm_shape(data_aggregated_sp) +
    tm_symbols(col = "red", size = 0.15, shape = 21, border.col = "black", border.lwd = 0.5) +
    tm_graticules(labels.size = 0.7) +
    tm_layout(inner.margins = 0, legend.outside = TRUE, legend.outside.position = "right")
  
  tmap::tmap_save(
    tm       = map_out,
    filename = file.path(outdir, paste0(sp_name, "_Ho_Map_SVM_CELLS", suffix, ".pdf")),
    width = 25, height = 15, units = "cm", dpi = 600
  )
  
  save.image(file.path(outdir, paste0(sp_name, "_Ho_Workspace_CELLS", suffix, ".RData")))
  
  message(sprintf("Pipeline completed successfully for '%s' using [%s | %s | Thinning: %dm]!",
                  sp_name, val_scheme_str, toupper(decorr_method), thin_dist))
  invisible(NULL)
}

# ==============================================================================
# GENERAL PATHS AND PARAMETERS
# ==============================================================================
n_cores         <- 6L          # fixed number of cores; edit here if needed
decorr_method   <- "vif"
vif_th          <- 10
cor_cutoff      <- 0.5
dom_max         <- 0.3
remove_outliers <- TRUE
iqr_mult        <- 1.5
thin_dist       <- 5000

out_dir  <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics"
r_dir    <- "D:/test_layers_30s"

d_path   <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx"
sdm_base <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval"

sdm_ca <- file.path(sdm_base, "Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif")
sdm_cm <- file.path(sdm_base, "Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif")

species_info <- list(
  # CA = list(name = "Crocodylus acutus",     sdm = sdm_ca),
  CM = list(name = "Crocodylus moreletii",  sdm = sdm_cm)
)

cv_configs <- list(
  LOOCV   = list(use_loocv = TRUE,  use_blockcv = FALSE),
  BLOCKCV = list(use_loocv = FALSE, use_blockcv = TRUE)
)

# ==============================================================================
# SCENARIO GRID
# ==============================================================================
scenarios_grid <- expand.grid(
  sp_code          = c("CM"),   # c("CA", "CM")
  addLonLat        = c(TRUE, FALSE),
  cv_type          = c("LOOCV", "BLOCKCV"),
  thin_dist        = c(15000),
  stringsAsFactors = FALSE
)

scenarios_grid <- scenarios_grid[-c(5, 7), ]
cat("Total scenarios to run:", nrow(scenarios_grid), "\n")

# ==============================================================================
# AUTOMATED EXECUTION
# ==============================================================================
results <- list()

for (i in 1:nrow(scenarios_grid)) {
  row     <- scenarios_grid[i, ]
  sp_data <- species_info[[row$sp_code]]
  cv_data <- cv_configs[[row$cv_type]]
  
  ll_tag      <- if (isTRUE(row$addLonLat)) "WITH_LONLAT" else "NO_LONLAT"
  
  folder_name <- paste(row$sp_code, "CELLS", ll_tag, row$cv_type, toupper(decorr_method),
                       paste0("THIN", row$thin_dist, "M"),
                       if (remove_outliers) "NOOUT", sep = "_")
  
  outdir_path <- file.path(out_dir, folder_name)
  
  cat(sprintf("\n------------------------------------------------------------\n"))
  cat(sprintf("[%d/%d] Running: %s\n", i, nrow(scenarios_grid), folder_name))
  cat(sprintf("------------------------------------------------------------\n"))
  
  results[[folder_name]] <- svr_ho_pipeline(
    outdir          = outdir_path,
    sp_name         = sp_data$name,
    raster_dir      = r_dir,
    data_path       = d_path,
    sdm_path        = sp_data$sdm,
    n_cores         = n_cores,
    decorr_method   = decorr_method,
    vif_th          = vif_th,
    cor_cutoff      = cor_cutoff,
    dom_max         = dom_max,
    addLonLat       = row$addLonLat,
    use_blockcv     = cv_data$use_blockcv,
    use_loocv       = cv_data$use_loocv,
    remove_outliers = remove_outliers,
    iqr_mult        = iqr_mult,
    kernel_type     = "svmRadial",
    thin_dist       = row$thin_dist
  )
  
  terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
}

cat("\nAll scenarios finished successfully!\n")

# ==============================================================================
# CONSOLIDATED METRICS (TRAIN AND TEST) ACROSS ALL SCENARIOS
# ==============================================================================
cat("\nConsolidating metrics from all scenarios...\n")

metrics_files <- list.files(
  path       = out_dir,
  pattern    = "_SVM_Metrics_CELLS.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(metrics_files) > 0) {
  all_metrics <- do.call(rbind, lapply(metrics_files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    df$scenario_folder <- basename(dirname(f))
    return(df)
  }))
  
  consolidated_output_path <- file.path(out_dir, "Consolidated_SVM_Metrics_All_Scenarios.csv")
  write.csv(all_metrics, consolidated_output_path, row.names = FALSE)
  
  cat(sprintf("Consolidation successful! File saved to:\n%s\n", consolidated_output_path))
} else {
  warning("No metrics files found to consolidate.")
}