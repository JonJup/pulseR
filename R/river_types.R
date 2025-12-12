#' Classify River Types using Regional Gaussian Mixture Models
#'
#' This function trains Gaussian Mixture Models (GMM) on specified core regions and
#' projects these classifications onto a larger dataset. The final output weighs the
#' predicted cluster probabilities by the observation's membership to specific regions.
#'
#' @param all_data An \code{sf} object or data frame containing the complete dataset
#'   to be classified. Must contain the variables used in the core region models.
#' @param core_regions A list of \code{sf} objects or data frames. Each element represents
#'   a distinct geographic or logical region used to train the classification models.
#' @param n_river_types Integer or integer vector. The number of mixture components
#'   (clusters) to fit in \code{mclust::Mclust}. Default is \code{1:9}.
#' @param membership_id A matrix or data frame containing the fuzzy membership weights
#'   of each observation in \code{all_data} relative to the \code{core_regions}.
#'   Columns must correspond to the order of regions in \code{core_regions}.
#' @param membership_cols Character vector. Names of columns in \code{core_regions}
#'   related to prior class membership that should be excluded from training.
#'   Default is \code{paste0("X", 1:3)}.
#' @param non_value_cols Character vector. Names of identifier columns to be removed
#'   before model training (e.g., primary keys). Default is \code{"ID"}.
#'
#' @details
#' The function performs several preprocessing steps:
#' \itemize{
#'   \item Geometry columns are dropped using \code{sf::st_drop_geometry}.
#'   \item Constant variables within specific regions are identified and removed to avoid
#'     singular covariance matrices during GMM fitting.
#'   \item Specific logic is applied to geological variables (\code{area_sediment},
#'     \code{area_calcareous}, \code{area_siliceous}). If two are constant, the third
#'     is removed to prevent perfect collinearity (compositional data issues).
#' }
#'
#' The prediction for a specific river type \eqn{k} in region \eqn{r} is calculated as:
#' \eqn{P(Type_k | Region_r) \times W_{r}}
#' Where \eqn{W_r} is the weight provided in \code{membership_id}.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{\code{rivertypes}}{A matrix of weighted probabilities. Column names follow the format
#'   "\eqn{Region{i}RiverType{k}}".}
#'   \item{\code{models}}{A list of the fitted \code{mclust::Mclust} model objects for each region.}
#' }
#'
#' @importFrom sf st_drop_geometry
#' @importFrom mclust Mclust predict.Mclust
#'
#' @export
river_types <- function(all_data, 
                        core_regions, 
                        n_river_types = 1:9, 
                        membership_id, 
                        membership_cols = paste0("X", 1:3), 
                        non_value_cols = "ID") {
        
        n_regions <- length(core_regions)
        
        # --- 1. Preprocessing ---
        
        # Extract data for classification (drop geometry)
        data <- lapply(core_regions, sf::st_drop_geometry)
        data2 <- sf::st_drop_geometry(all_data)
        
        # Remove specific ID columns and membership columns from training data
        data <- lapply(data, function(x) {
                cols_to_remove <- c("core_region_id", membership_cols, non_value_cols)
                x[, !names(x) %in% cols_to_remove, drop = FALSE]
        })
        
        # Remove non-value columns from prediction data
        data2 <- data2[, !names(data2) %in% non_value_cols, drop = FALSE]
        
        # --- 2. Variable Selection & Collinearity Checks ---
        
        # Test if any values are constant within a region and handle compositional geology data
        for (i in 1:n_regions) {
                
                # Identify constant columns (variance = 0 / length unique = 1)
                is_constant <- apply(data[[i]], 2, function(x) length(unique(x))) == 1
                
                if (any(is_constant)) {
                        const_id <- which(is_constant)
                        const_names <- names(const_id)
                        
                        # Remove constant columns
                        data[[i]] <- data[[i]][, -const_id, drop = FALSE]
                        
                        # Special handling for geological compositional data to avoid singularity
                        # If 2 out of 3 geological variables are constant, remove the 3rd to avoid linear dependency
                        geo_vars <- c("area_sediment", "area_calcareous", "area_siliceous")
                        
                        if (all(c("area_sediment", "area_calcareous") %in% const_names) && "area_siliceous" %in% names(data[[i]])) {
                                data[[i]] <- data[[i]][, -which(names(data[[i]]) == "area_siliceous")]
                        }
                        if (all(c("area_siliceous", "area_calcareous") %in% const_names) && "area_sediment" %in% names(data[[i]])) {
                                data[[i]] <- data[[i]][, -which(names(data[[i]]) == "area_sediment")]
                        }
                        if (all(c("area_sediment", "area_siliceous") %in% const_names) && "area_calcareous" %in% names(data[[i]])) {
                                data[[i]] <- data[[i]][, -which(names(data[[i]]) == "area_calcareous")]
                        }
                }
        }
        
        # --- 3. Model Training ---
        
        # Fit Gaussian Mixture Models
        clusters <- lapply(data, mclust::Mclust, G = n_river_types)
        
        raw <- weighted <- list()
        
        # --- 4. Prediction & Weighting ---
        
        for (i in 1:n_regions) {
                
                iterMod <- clusters[[i]]
                
                # Extract only the variables used in this specific region's model
                model_vars <- rownames(iterMod$parameters$mean)
                model_data <- data2[, which(names(data2) %in% model_vars)]
                
                # Predict probabilities (z) using standard generic predict()
                # R will dispatch to predict.Mclust automatically because we imported it
                pred_rt <- predict(iterMod, newdata = model_data)
                
                # Mclust prediction returns a list, component 'z' holds the probabilities
                pred_rt <- pred_rt$z
                pred_rt <- round(pred_rt, 2)
                
                raw[[i]] <- pred_rt
                
                # Weight probabilities by regional membership
                weighted[[i]] <- pred_rt * membership_id[, i]
                colnames(weighted[[i]]) <- paste0("Region", i, "RiverType", 1:iterMod$G)
        }
        
        # Combine results
        out <- do.call(cbind, weighted)
        
        return(list(rivertypes = out, models = clusters))
}