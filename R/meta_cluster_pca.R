#' Calculate Bhattacharyya Distance between two Multivariate Normal Distributions
#'
#' Computes the Bhattacharyya distance, which measures the similarity of two probability
#' distributions. This implementation assumes multivariate normal distributions.
#'
#' @param mu1 Numeric vector. Mean vector of the first cluster.
#' @param sigma1 Matrix. Covariance matrix of the first cluster.
#' @param mu2 Numeric vector. Mean vector of the second cluster.
#' @param sigma2 Matrix. Covariance matrix of the second cluster.
#'
#' @return A single numeric value representing the distance.
#'
#' @export
bhattacharyya_dist <- function(mu1, sigma1, mu2, sigma2) {
        # Average covariance
        sigma_avg <- (sigma1 + sigma2) / 2
        
        # Term 1: Difference in means weighted by covariance
        diff <- mu1 - mu2
        term1 <- 1/8 * t(diff) %*% solve(sigma_avg) %*% diff
        
        # Term 2: Difference in covariance shapes
        det_avg <- det(sigma_avg)
        det1 <- det(sigma1)
        det2 <- det(sigma2)
        
        # Guard against negative determinants in numerical approximations
        term2 <- 0.5 * log(det_avg / sqrt(max(det1 * det2, 1e-10)))
        
        return(as.numeric(term1 + term2))
}

#' Visualize River Type Centroids in Environmental Space
#'
#' Extracts cluster centroids (means) from a list of regional Mclust models,
#' performs Principal Component Analysis (PCA) on the common environmental variables,
#' and visualizes the river types in 2D space.
#'
#' @param models A list of \code{mclust::Mclust} model objects. This is typically the
#'   \code{models} element returned by the \code{river_types} function.
#' @param plot Logical. If \code{TRUE} (default), returns a ggplot object. If \code{FALSE},
#'   returns the PCA data frame.
#'
#' @return A \code{ggplot2} object showing the centroids of river types, colored by region.
#'   If \code{plot = FALSE}, returns the data frame used for plotting.
#'
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point geom_text labs theme_minimal
#' @export
plot_river_centroids <- function(models, plot = TRUE) {
        
        # 1. Extract parameters from regional models
        all_types <- list()
        
        # Iterate through regions (models)
        for (i in seq_along(models)) {
                model <- models[[i]]
                region_name <- if (!is.null(names(models)[i])) names(models)[i] else paste0("Region", i)
                
                # Iterate through clusters (G)
                for (k in 1:model$G) {
                        type_id <- paste0(region_name, "_Type", k)
                        
                        # Handle cases where sigma might be stored differently (e.g., spherical models)
                        # Mclust stores sigma in $parameters$variance$sigma. If 3D array: [vars, vars, G]
                        sig <- if (length(dim(model$parameters$variance$sigma)) == 3) {
                                model$parameters$variance$sigma[,, k]
                        } else if (is.matrix(model$parameters$variance$sigma)) {
                                model$parameters$variance$sigma # If shared sigma
                        } else {
                                diag(model$parameters$variance$sigma[,k]) # If diagonal
                        }
                        
                        all_types[[type_id]] <- list(
                                mu = model$parameters$mean[, k],
                                sigma = sig,
                                region = region_name
                        )
                }
        }
        
        # 2. Variable Harmonization
        # Because 'river_types' drops constant variables per region, models may have
        # different sets of variables. We must intersect them for PCA.
        var_list <- lapply(all_types, function(x) names(x$mu))
        common_vars <- Reduce(intersect, var_list)
        
        if (length(common_vars) < 2) {
                stop("Fewer than 2 common variables found across regions. PCA cannot run.")
        }
        
        # Filter means to common variables and bind
        means_matrix <- do.call(rbind, lapply(all_types, function(x) x$mu[common_vars]))
        
        # 3. Run PCA
        # scale. = TRUE is critical as environmental vars likely have different units
        pca_res <- stats::prcomp(means_matrix, scale. = TRUE)
        
        # 4. Prepare Plot Data
        pca_df <- as.data.frame(pca_res$x)
        pca_df$TypeID <- rownames(pca_df)
        pca_df$Region <- sapply(all_types, function(x) x$region)
        
        if (!plot) {
                return(pca_df)
        }
        
        # 5. Plot
        p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = Region, label = TypeID)) +
                ggplot2::geom_point(size = 3) +
                ggplot2::geom_text(vjust = 1.5, size = 3, show.legend = FALSE) +
                ggplot2::labs(
                        title = "River Type Centroids in Environmental Space",
                        subtitle = paste("Based on", length(common_vars), "common variables across", length(models), "regions"),
                        x = "Principal Component 1",
                        y = "Principal Component 2"
                ) +
                ggplot2::theme_minimal()
        
        return(p)
}