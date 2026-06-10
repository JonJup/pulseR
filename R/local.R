#' Classify local Types using Regional Gaussian Mixture Models
#'
#' Trains a Gaussian Mixture Model (GMM) on each of several user-defined
#' "core regions" and projects those classifications onto a larger target
#' dataset. For every observation the predicted local-type probabilities from
#' each regional model are weighted by that observation's (fuzzy or crisp)
#' membership in the corresponding region and concatenated into a single wide
#' probability matrix.
#'
#' @param all_data An \code{sf} object or \code{data.frame} containing the
#'   complete dataset to be classified. Must contain (at least) every predictor
#'   variable used by the regional models. Any \code{sf} geometry column is
#'   dropped automatically.
#' @param core_regions A non-empty list of \code{sf} objects or
#'   \code{data.frame}s. Each element holds the training observations for one
#'   region. All elements must contain the same predictor columns.
#' @param n_local_types Integer scalar or vector passed to the \code{G} argument
#'   of \code{\link[mclust]{Mclust}}; the number(s) of mixture components to
#'   consider. When a vector is supplied, \code{Mclust} selects the best value
#'   by BIC. Default \code{1:9}.
#' @param membership_id A numeric matrix or \code{data.frame} of regional
#'   membership weights with one row per observation in \code{all_data} and one
#'   column per element of \code{core_regions} (in the same order). Used to
#'   weight each regional model's predictions.
#' @param membership_cols Character vector naming prior-membership columns in
#'   \code{core_regions} that must be excluded from training. When \code{NULL}
#'   (default) it is set to \code{paste0("X", seq_len(ncol(membership_id)))}.
#' @param non_value_cols Character vector of identifier columns to drop before
#'   training (e.g. primary keys). Default \code{"ID"}. The column
#'   \code{"core_region_id"} is always dropped from the training data.
#' @param crisp Logical. When \code{TRUE} both the regional memberships and the
#'   per-region local-type probabilities are hardened to 0/1 (each observation
#'   assigned to its single most probable region and, within a region, its
#'   single most probable local type) before weighting. Default \code{FALSE}.
#' @param verbose Logical. When \code{TRUE} progress over regions is reported
#'   via \code{\link{message}}. Default \code{FALSE}.
#' @param modelNames Character vector indicating the models to be fitted in the 
#'   EM phase of clustering. Is passed directly to mclust::Mclust, see there for
#   for more details
#'
#' @details
#' Pre-processing applied before fitting:
#' \itemize{
#'   \item \code{sf} geometry is removed with \code{\link[sf]{st_drop_geometry}}
#'     (a no-op for plain data frames).
#'   \item Identifier and prior-membership columns are removed.
#'   \item Columns that are constant in \emph{any} region are removed from
#'     \emph{all} regions, since a zero-variance predictor produces a singular
#'     covariance matrix in \code{Mclust}.
#'   \item Compositional geology handling: the variables \code{area_sediment},
#'     \code{area_calcareous} and \code{area_siliceous} sum to a constant, so
#'     including all three is perfectly collinear. If any one is found constant,
#'     all three are dropped (see \strong{Note}).
#' }
#'
#' For region \eqn{r} and local type \eqn{k} the returned weighted probability
#' is \eqn{P(\mathrm{Type}_k \mid \mathrm{Region}_r)\, W_r}, where \eqn{W_r} is
#' the membership weight in column \eqn{r} of \code{membership_id}.
#'
#' @note The geology rule drops \emph{all three} compositional variables when
#'   \emph{any one} is constant. This is more aggressive than strictly required
#'   to break the simplex (sum-to-constant) constraint, for which dropping a
#'   single component would suffice. The behaviour is preserved from the
#'   original implementation; verify it matches your intended preprocessing.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{localtypes}}{Numeric matrix of membership-weighted local-type
#'     probabilities, one row per observation in \code{all_data}. Columns are
#'     named \code{"Region{i}localType{k}"}.}
#'   \item{\code{localtypes_raw}}{List (one element per region, named
#'     \code{"Region{i}"}) of the \emph{unweighted} predicted probability
#'     matrices.}
#'   \item{\code{models}}{List of the fitted \code{Mclust} objects, named by
#'     region.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' make_region <- function(n, shift) {
#'   data.frame(
#'     ID             = seq_len(n),
#'     core_region_id = 1,
#'     X1 = stats::runif(n), X2 = stats::runif(n), # prior-membership cols (dropped)
#'     temp     = stats::rnorm(n, shift),
#'     slope    = stats::rnorm(n, -shift),
#'     rainfall = stats::rnorm(n)
#'   )
#' }
#' regions  <- list(make_region(100, 0), make_region(100, 3))
#' all_data <- rbind(make_region(50, 0), make_region(50, 3))
#'
#' memb <- matrix(stats::runif(nrow(all_data) * 2), ncol = 2)
#' memb <- memb / rowSums(memb)                    # fuzzy weights summing to 1
#'
#' res <- get_local_types(all_data, regions, n_local_types = 1:4, membership_id = memb)
#' head(res$localtypes)
#' }
#'
#' @importFrom sf st_drop_geometry
#' @importFrom mclust Mclust mclustBIC
#' @importFrom stats predict complete.cases
#'
#' @export
get_local_types <- function(all_data,
                        core_regions,
                        n_local_types = 1:9,
                        membership_id,
                        membership_cols = NULL,
                        non_value_cols  = "ID",
                        crisp           = FALSE,
                        verbose         = FALSE,
                        modelNames = NULL) {
        
        ## ---- 0. Input validation --------------------------------------------------
        if (!is.data.frame(all_data)) {
                stop("`all_data` must be a data.frame or sf object.", call. = FALSE)
        }
        if (!is.list(core_regions) || length(core_regions) == 0L) {
                stop("`core_regions` must be a non-empty list of data.frames / sf objects.",
                     call. = FALSE)
        }
        if (!all(vapply(core_regions, is.data.frame, logical(1)))) {
                stop("Every element of `core_regions` must be a data.frame or sf object.",
                     call. = FALSE)
        }
        if (missing(membership_id)) {
                stop("`membership_id` is required.", call. = FALSE)
        }
        if (!is.logical(crisp) || length(crisp) != 1L || is.na(crisp)) {
                stop("`crisp` must be a single TRUE/FALSE value.", call. = FALSE)
        }
        if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
                stop("`verbose` must be a single TRUE/FALSE value.", call. = FALSE)
        }
        if (!is.numeric(n_local_types) || length(n_local_types) == 0L ||
            any(n_local_types < 1L) || any(n_local_types != as.integer(n_local_types))) {
                stop("`n_local_types` must be one or more positive integers.", call. = FALSE)
        }
        if (!is.null(non_value_cols) && !is.character(non_value_cols)) {
                stop("`non_value_cols` must be a character vector or NULL.", call. = FALSE)
        }
        if (!is.null(membership_cols) && !is.character(membership_cols)) {
                stop("`membership_cols` must be a character vector or NULL.", call. = FALSE)
        }
        
        n_regions <- length(core_regions)
        
        ## Coerce membership weights to a numeric matrix for predictable indexing,
        ## then check that its dimensions are consistent with the other inputs.
        membership_id <- as.matrix(membership_id)
        if (!is.numeric(membership_id)) {
                stop("`membership_id` must be numeric.", call. = FALSE)
        }
        if (ncol(membership_id) != n_regions) {
                stop(sprintf(
                        "`membership_id` has %d columns but there are %d core regions; they must match.",
                        ncol(membership_id), n_regions), call. = FALSE)
        }
        if (nrow(membership_id) != nrow(all_data)) {
                stop(sprintf(
                        "`membership_id` has %d rows but `all_data` has %d; they must match.",
                        nrow(membership_id), nrow(all_data)), call. = FALSE)
        }
        
        ## Default prior-membership column names are derived from the number of
        ## regions (X1, X2, ...). These are columns to *exclude from training*, not
        ## the weights themselves.
        if (is.null(membership_cols)) {
                membership_cols <- paste0("X", seq_len(ncol(membership_id)))
        }
        
        ## Local helper: harden each row of a numeric matrix to a 0/1 indicator of
        ## its single largest value (ties -> first column; an all-NA row -> all 0).
        ## Vectorised and shape-stable for any number of columns (incl. one).
        harden_rows <- function(m) {
                m   <- as.matrix(m)
                out <- matrix(0, nrow = nrow(m), ncol = ncol(m))
                j   <- max.col(replace(m, is.na(m), -Inf), ties.method = "first")
                out[cbind(seq_len(nrow(m)), j)] <- 1
                all_na <- rowSums(!is.na(m)) == 0L
                if (any(all_na)) out[all_na, ] <- 0
                out
        }
        
        ## ---- 1. Pre-processing ----------------------------------------------------
        
        ## Drop sf geometry (no-op on a plain data.frame) and standardise to data.frame.
        train_list   <- lapply(core_regions, function(x) as.data.frame(sf::st_drop_geometry(x)))
        predict_data <- as.data.frame(sf::st_drop_geometry(all_data))
        
        ## Columns that must never enter a model: the region key, the prior-membership
        ## columns and any user-supplied identifier columns.
        drop_cols  <- unique(c("core_region_id", membership_cols, non_value_cols))
        train_list <- lapply(train_list, function(x) x[, !names(x) %in% drop_cols, drop = FALSE])
        
        ## Prediction data only needs the identifier columns removed; the predictors
        ## actually used are selected by name per region further below.
        predict_data <- predict_data[, !names(predict_data) %in% non_value_cols, drop = FALSE]
        
        ## All regions must share the same predictor set (the global constant-column
        ## logic below assumes this).
        ref_cols <- names(train_list[[1L]])
        if (!all(vapply(train_list, function(x) setequal(names(x), ref_cols), logical(1)))) {
                stop("All core regions must contain the same predictor columns.", call. = FALSE)
        }
        
        ## ---- 1b. Optional crisp (hard) regional membership ------------------------
        if (crisp) membership_id <- harden_rows(membership_id)
        
        ## ---- 2. Drop constant & collinear compositional predictors ----------------
        
        ## A column constant within ANY region is dropped from ALL regions: a
        ## zero-variance predictor yields a singular covariance matrix in Mclust.
        ## NB: a column that is entirely NA also has a single unique value and is
        ## therefore also flagged here.
        const_names <- character(0)
        for (i in seq_len(n_regions)) {
                n_unique    <- vapply(train_list[[i]], function(x) length(unique(x)), integer(1))
                const_names <- union(const_names, names(n_unique)[n_unique <= 1L])
        }
        
        ## Remove constant variables by name
        if (length(const_names)) {
                train_list <- lapply(train_list, function(x) x[, !names(x) %in% const_names, drop = FALSE])
        }
        
        ## Predictors must be numeric and complete; otherwise Mclust fails cryptically.
        for (i in seq_len(n_regions)) {
                if (ncol(train_list[[i]]) == 0L) {
                        stop(sprintf("Region %d has no predictor columns left after pre-processing.", i),
                             call. = FALSE)
                }
                non_num <- !vapply(train_list[[i]], is.numeric, logical(1))
                if (any(non_num)) {
                        stop(sprintf("Region %d has non-numeric predictor column(s): %s.",
                                     i, paste(names(train_list[[i]])[non_num], collapse = ", ")),
                             call. = FALSE)
                }
                if (anyNA(train_list[[i]])) {
                        stop(sprintf("Region %d training data contains missing values; Mclust requires complete cases.", i),
                             call. = FALSE)
                }
        }
        
        ## ---- 3. Fit one Gaussian Mixture Model per region -------------------------
        ## (No scaling: Mclust models the full covariance structure, so standardising
        ## the predictors is unnecessary. The variables used to train each region are
        ## recorded explicitly rather than recovered from the fitted object, which is
        ## robust to the univariate case where `parameters$mean` has no rownames.)
        models     <- vector("list", n_regions)
        train_vars <- vector("list", n_regions)
        
        for (i in seq_len(n_regions)) {
                if (verbose) message(sprintf("Fitting GMM for region %d/%d ...", i, n_regions))
                
                train_vars[[i]] <- names(train_list[[i]])
                
                ## IMPORTANT: Mclust() rebuilds its arguments into a call to mclustBIC() and
                ## evaluates that call in the calling frame (this function's frame), not in
                ## mclust's own namespace. So mclustBIC must be reachable from here!
                fit <- mclust::Mclust(train_list[[i]], G = n_local_types, modelNames = modelNames)
                if (is.null(fit)) {
                        stop(sprintf("Mclust could not fit a model for region %d over G = %s.",
                                     i, paste(range(n_local_types), collapse = ":")), call. = FALSE)
                }
                models[[i]] <- fit
        }
        names(models) <- paste0("Region", seq_len(n_regions))
        
        ## ---- 4. Predict on `all_data` and weight by regional membership -----------
        raw      <- vector("list", n_regions)
        weighted <- vector("list", n_regions)
        
        for (i in seq_len(n_regions)) {
                vars_i <- train_vars[[i]]
                
                ## Every training variable must be present in the target data.
                missing_vars <- setdiff(vars_i, names(predict_data))
                if (length(missing_vars)) {
                        stop(sprintf("`all_data` is missing variable(s) required by region %d: %s.",
                                     i, paste(missing_vars, collapse = ", ")), call. = FALSE)
                }
                
                model_data <- predict_data[, vars_i, drop = FALSE]
                
                non_num <- !vapply(model_data, is.numeric, logical(1))
                if (any(non_num)) {
                        stop(sprintf("`all_data` column(s) used by region %d are not numeric: %s.",
                                     i, paste(vars_i[non_num], collapse = ", ")), call. = FALSE)
                }
                
                ## Mclust prediction needs complete, finite rows. Flag offenders in a
                ## vectorised way (avoids a slow row-wise apply on large datasets) and
                ## report how many and where.
                mat <- as.matrix(model_data)
                bad <- !stats::complete.cases(mat) | rowSums(!is.finite(mat)) > 0L
                if (any(bad)) {
                        offenders <- which(bad)
                        stop(sprintf(
                                "%d row(s) in `all_data` have missing/non-finite values for region %d predictors (%s). First offending row(s): %s.",
                                length(offenders), i, paste(vars_i, collapse = ", "),
                                paste(offenders[seq_len(min(6L, length(offenders)))], collapse = ", ")),
                             call. = FALSE)
                }
                
                ## predict.Mclust is dispatched through stats::predict. $z holds
                ##  the posterior cluster probabilities as an n-by-G matrix.
                pred <- stats::predict(models[[i]], newdata = model_data)$z
                ## guarantee matrix form (e.g. when G = 1)
                pred <- as.matrix(pred)   
                
                if (crisp) {
                        pred <- harden_rows(pred)
                }
                
                raw[[i]] <- pred
                
                ## Weight every local-type column by this region's membership weight. The
                ## length-n vector recycles down the columns, scaling each row by its weight.
                w <- pred * membership_id[, i]
                colnames(w) <- paste0("R", i, "\u00b7","T", seq_len(ncol(w)))
                weighted[[i]] <- w
        }
        names(raw) <- paste0("Region", seq_len(n_regions))
        
        ## ---- 5. Assemble output ---------------------------------------------------
        out <- do.call(cbind, weighted)
        
        list(
                localtypes     = out,
                localtypes_raw = raw,
                models         = models
        )
}
