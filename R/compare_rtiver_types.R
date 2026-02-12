#' Compare River Types Across Multiple Mixture Models
#'
#' Computes pairwise Bhattacharyya distances between all mixture components
#' (river types) across a set of fitted Gaussian mixture models. This produces
#' a symmetric distance matrix that quantifies the similarity of component
#' distributions, enabling cross-model comparison of derived typologies.
#'
#' @param models A list of fitted Gaussian mixture model objects. Each model
#'   must contain a `parameters` element with:
#'   \describe{
#'     \item{`mean`}{A named matrix of component means with variables as rows
#'       and components as columns.}
#'     \item{`variance$sigma`}{A three-dimensional array of covariance matrices
#'       with dimensions \eqn{p \times p \times K}, where \eqn{p} is the number
#'       of variables and \eqn{K} is the number of components.}
#'   }
#'   This structure is compatible with objects returned by
#'   \code{\link[mclust]{Mclust}}.
#'
#' @return A symmetric numeric matrix of pairwise Bhattacharyya distances
#'   between all component pairs across all models, with zeros on the diagonal.
#'   Rows and columns are labelled as \code{"M{i}_K{k}"}, where \code{i} is the
#'   model index and \code{k} is the component index within that model. Returns
#'   \code{NA_real_} with a warning if no variables are shared across all models.
#'
#' @details
#' The function first identifies the intersection of variable names (row names
#' of each model's mean matrix) shared across all models to ensure comparable
#' multivariate distributions. It then enumerates every model–component
#' combination as a distinct "type" and computes the Bhattacharyya distance
#' for each unique pair using their mean vectors and covariance matrices
#' restricted to the shared variable set.
#'
#' The Bhattacharyya distance between two multivariate Gaussian distributions
#' \eqn{\mathcal{N}(\mu_1, \Sigma_1)} and \eqn{\mathcal{N}(\mu_2, \Sigma_2)}
#' is defined as:
#' \deqn{D_B = \frac{1}{8} (\mu_1 - \mu_2)^T \Sigma^{-1} (\mu_1 - \mu_2) +
#'   \frac{1}{2} \ln \frac{|\Sigma|}{|\Sigma_1|^{1/2} |\Sigma_2|^{1/2}}}
#' where \eqn{\Sigma = \frac{\Sigma_1 + \Sigma_2}{2}}.
#'
#' @note The function currently references an undefined variable \code{all_vars}
#'   for subsetting means and covariance matrices. This should likely be
#'   \code{shared}, the variable computed from the intersection of model
#'   variable names.
#'
#' @seealso \code{\link{bhattacharyya_dist}} for the underlying distance
#'   calculation; \code{\link[mclust]{Mclust}} for fitting Gaussian mixture
#'   models.
#'
#' @examples
#' \dontrun{
#' library(mclust)
#'
#' # Fit mixture models to two river datasets
#' mod1 <- Mclust(river_data_alpine, G = 3)
#' mod2 <- Mclust(river_data_lowland, G = 4)
#'
#' # Compute pairwise distances between all 3 + 4 = 7 types
#' dist_mat <- compare_river_types(list(mod1, mod2))
#'
#' # Visualise as a heatmap
#' heatmap(dist_mat, symm = TRUE)
#'
#' # Hierarchical clustering of river types
#' hc <- hclust(as.dist(dist_mat))
#' plot(hc)
#' }
#'
#' @export


compare_river_types <- function(models) {
        
        # identify the global intersection of variables 
        shared <- Reduce(intersect, lapply(models, function(m) rownames(m$parameters$mean)))
        
        if (length(shared) == 0) {
                warning("No shared variables between all models")
                return(NA_real_)
        }
        
        # Build a list of all types (model + component combinations)
        types <- do.call(rbind, lapply(seq_along(models), function(i) {
                m <- models[[i]]
                K <- ncol(m$parameters$mean)
                data.frame(model = i, component = seq_len(K))
        }))

        # Give each type a label
        types$label <- paste0("M", types$model, "_K", types$component)
        
        n <- nrow(types)
        dist_mat <- matrix(NA, n, n, dimnames = list(types$label, types$label))
        
        for (i in 1:(n - 1)) {
                for (j in (i + 1):n) {
                        m1 <- models[[types$model[i]]]
                        m2 <- models[[types$model[j]]]
                        
                        mu1    <- m1$parameters$mean[shared, types$component[i]]
                        sigma1 <- m1$parameters$variance$sigma[shared, shared, types$component[i]]
                        
                        mu2    <- m2$parameters$mean[shared, types$component[j]]
                        sigma2 <- m2$parameters$variance$sigma[shared, shared, types$component[j]]
                        
                        d <- bhattacharyya_dist(mu1, sigma1, mu2, sigma2)
                        dist_mat[i, j] <- d
                        dist_mat[j, i] <- d
                }
        }
        
        diag(dist_mat) <- 0
}


