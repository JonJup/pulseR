# -----------------------------------------------------------------------------*
# Helpers
# -----------------------------------------------------------------------------*

#' First non-NULL argument.
#' @keywords internal
#' @noRd
.first_non_null <- function(...) {
        for (x in list(...)) if (!is.null(x)) return(x)
        NULL
}

#' Extract the N x k region weight matrix from whatever the caller passed.
#'
#' Accepts a "pulse_regions" object, a "regions_tuning" object (whose fit lives
#' in $best_result), or a bare numeric matrix. `which = "typicality"` looks for
#' the possibilistic matrix under any of the names the contract has used.
#'
#' @keywords internal
#' @noRd
.region_weights <- function(regions, which = c("memberships", "typicality")) {
        which <- match.arg(which)
        
        if (is.matrix(regions) || is.data.frame(regions)) return(as.matrix(regions))
        if (!is.list(regions))
                stop("`regions` must be a 'pulse_regions' / 'regions_tuning' object or a matrix.",
                     call. = FALSE)
        
        ## FIX: tune_regions() returns the fit under $best_result, not at top level.
        src <- if (inherits(regions, "regions_tuning")) regions$best_result else regions
        if (is.null(src))
                stop("`regions` looks like a tuning object but carries no `$best_result`.",
                     call. = FALSE)
        
        W <- if (which == "memberships") {
                ## FIX: the old code read `$memberships_normalized`, which is not part of
                ## the pulse_regions contract. `$memberships` is already row-normalised.
                .first_non_null(src$memberships, src$fit$memberships)
        } else {
                ## The possibilistic matrix lives in $fit and is now called `cohesion`;
                ## accept every name it has had (typicality, affinity).
                .first_non_null(src$cohesion, src$typicality, src$affinity,
                                src$fit$cohesion, src$fit$typicality, src$fit$affinity)
        }
        
        if (is.null(W))
                stop(sprintf("Could not find a `%s` matrix on `regions`.%s", which,
                             if (which == "typicality")
                                     " Refit with a version of cluster_consensus() that returns `cohesion`."
                             else ""), call. = FALSE)
        W <- as.matrix(W)
        if (!is.numeric(W))
                stop("Region weights must be numeric.", call. = FALSE)
        if (anyNA(W))
                stop(sprintf(paste0(
                        "Region weights contain NA in %d row(s). These are units with ",
                        "pi_support == 0 (a one-unit piece in every ensemble partition); ",
                        "Level 2 handling for them is not yet defined."),
                        sum(rowSums(is.na(W)) > 0)), call. = FALSE)
        if (any(W < 0))
                stop("Region weights must be non-negative.", call. = FALSE)
        W
}

#' Weighted BIC of one me.weighted() fit.
#'
#' 2 * (weighted log-likelihood) - npar * log(n_eff), with n_eff = sum(weights).
#'
#' NOTE: this assumes `me.weighted()$loglik` is already the WEIGHTED
#' log-likelihood, which is the mclust convention. The previous version applied
#' a further `* mean(weights)` to the likelihood term but not to the penalty,
#' which changes the balance between them and therefore changes which G wins —
#' it is not a harmless rescaling. If you ever need the unweighted-scale
#' likelihood, rescale BOTH terms.
#'
#' @keywords internal
#' @noRd
.weighted_bic <- function(loglik, model_name, d, G, n_eff) {
        if (!is.finite(loglik) || !is.finite(n_eff) || n_eff <= 1) return(NA_real_)
        npar <- mclust::nMclustParams(modelName = model_name, d = d, G = G)
        2 * loglik - npar * log(n_eff)
}

#' Validate and resolve the candidate covariance models for d predictors.
#'
#' Multivariate names are mapped to their volume letter ("E"/"V") when only one
#' predictor is left, because mclust uses different names in one dimension.
#' Order is preserved: on an exact BIC tie the earlier candidate wins.
#'
#' @keywords internal
#' @noRd
.resolve_model_names <- function(model_name, d) {
        uni   <- c("E", "V")
        multi <- mclust::mclust.options("emModelNames")
        bad   <- setdiff(model_name, c(uni, multi))
        if (length(bad))
                stop(sprintf("Unknown mclust model name(s): %s. Valid: %s.",
                             paste(bad, collapse = ", "),
                             paste(c(multi, uni), collapse = ", ")), call. = FALSE)
        if (d == 1L) return(unique(substr(model_name, 1L, 1L)))
        m <- intersect(model_name, multi)
        if (!length(m))
                stop(sprintf("`model_name` holds only univariate models (%s) but %d predictors remain.",
                             paste(model_name, collapse = ", "), d), call. = FALSE)
        m
}

#' Column-wise variance without coercing the whole frame to one matrix.
#' @keywords internal
#' @noRd
.col_var <- function(X, rows) {
        vapply(seq_len(ncol(X)), function(j) stats::var(X[rows, j]), numeric(1))
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

#' Local Types within Regions via Weighted Gaussian Mixture Models
#'
#' Refines each region of a \code{"pulse_regions"} object into \emph{local
#' types} by fitting one membership-weighted Gaussian mixture model per region
#' and concatenating the weighted component probabilities into a single wide
#' matrix.
#'
#' @details
#' For region \eqn{r}, every unit carries a weight \eqn{w_{ir}} — its membership
#' in (or typicality for) that region. Units with \eqn{w_{ir} \ge}
#' \code{cutoff} form the region's \emph{core} and provide the model-based
#' hierarchical initialisation (\code{\link[mclust]{hc}} +
#' \code{\link[mclust]{hclass}}); every other unit starts from a flat
#' \eqn{1/G} row. The EM step (\code{\link[mclust]{me.weighted}}) then runs on
#' all units with \eqn{w_{ir} >} \code{weight_eps}, weighted by \eqn{w_{ir}},
#' so a unit informs a region's mixture in proportion to how much it belongs
#' to that region.
#'
#' The number of components \emph{and} the covariance model are chosen jointly
#' per region by weighted BIC, \eqn{2\ell - p \log n_{\mathrm{eff}}} with
#' \eqn{n_{\mathrm{eff}} = \sum_i w_{ir}} over the fitted units. A
#' (model, G) candidate is only fitted if its parameter count \eqn{p} is below
#' \eqn{n_{\mathrm{eff}}}; this is what lets small regions fall back to
#' shared, diagonal or spherical covariances instead of failing.
#'
#' The hierarchical initialisation uses \code{"VVV"} when a region's core has
#' more units than predictors and \code{"EII"} otherwise. It is computed once
#' per region and shared by every candidate model, so all candidates start
#' from the same partition.
#'
#' The returned local-type probability for region \eqn{r}, component \eqn{k} is
#' \eqn{z_{irk} w_{ir}}. Because each \eqn{z_i} row sums to one, the raw row
#' sums of the combined matrix equal \eqn{\sum_r w_{ir}} — exactly 1 when
#' \code{weights = "memberships"}, and something else when
#' \code{weights = "typicality"}, which is why \code{memberships} is returned
#' row-normalised and the unnormalised version is kept as \code{weighted}.
#'
#' Units with \eqn{w_{ir} \le} \code{weight_eps} are excluded from region
#' \eqn{r}'s fit and receive a zero in that region's columns. This is exact,
#' not an approximation: their contribution \eqn{z_{irk} w_{ir}} is negligible
#' by construction, and it is what keeps the fit affordable at continental
#' \eqn{N}.
#'
#' @param env A \code{data.frame} or \code{sf} object with one row per unit, in
#'   the same order as the rows of the region weight matrix. Geometry is
#'   dropped automatically.
#' @param regions A \code{"pulse_regions"} object from \code{\link{get_regions}},
#'   a \code{"regions_tuning"} object from \code{\link{tune_regions}}, or a bare
#'   N x k numeric weight matrix.
#' @param variables Character vector of predictor columns to use. \code{NULL}
#'   (default) uses every numeric column of \code{env} that is not listed in
#'   \code{non_value_cols}.
#' @param non_value_cols Character vector of identifier columns to exclude
#'   (e.g. \code{"ID"}). Default \code{NULL}.
#' @param n_local_types Integer vector of candidate component counts. Default
#'   \code{2:9}. Candidates larger than a region's core allow are dropped for
#'   that region only.
#' @param cutoff Numeric in \eqn{(0, 1]}. Weight at or above which a unit joins
#'   a region's core. Must exceed \code{weight_eps}. Default 0.1. Note that its
#'   meaning differs between the two \code{weights} settings: memberships shrink
#'   as \eqn{k} grows (a mean membership at \eqn{k = 20} is 0.05), whereas
#'   typicality is a co-association probability and does not.
#' @param weights One of \code{"memberships"} (default; row-normalised, rows sum
#'   to 1) or \code{"typicality"} (possibilistic co-association).
#' @param model_name Character vector of candidate \pkg{mclust} covariance
#'   models, compared by weighted BIC within each region. Default
#'   \code{c("VVV", "EEE", "VVI", "VII")}: full, shared full, diagonal,
#'   spherical. Pass \code{"VVV"} to reproduce the single-family behaviour of
#'   earlier versions. \code{"EEI"} and \code{"EII"} are valid but excluded
#'   by default: weighted EM for equal-volume models can need ~40x more
#'   iterations than \code{"VVV"}. With a single predictor the names are
#'   mapped to \code{"E"}/\code{"V"}.
#' @param hc_use Character scalar passed to \code{\link[mclust]{hc}} as
#'   \code{use}. Set explicitly (default \code{"VARS"}) because the
#'   \pkg{mclust} default has changed between versions.
#' @param drop_constant Logical. Drop predictors that are constant within any
#'   region's core from \emph{all} regions, since a zero-variance predictor
#'   gives a singular covariance. Default \code{TRUE}.
#' @param weight_eps Numeric. Units with weight at or below this are excluded
#'   from a region's fit. Default 1e-6.
#' @param prior Prior passed to \code{\link[mclust]{me.weighted}}. Default
#'   \code{mclust::priorControl()}; \code{NULL} disables regularisation.
#' @param keep_models Logical. Return the selected fits (with the large
#'   \code{$z} slot stripped). Default \code{TRUE}.
#' @param verbose Logical. Report progress. Default \code{FALSE}.
#'
#' @return An object of class \code{"pulse_local_types"}: a list with
#'   \describe{
#'     \item{\code{memberships}}{N x sum(G_r) matrix, rows sum to 1.}
#'     \item{\code{weighted}}{The same matrix before row-normalisation.}
#'     \item{\code{hard_classes}}{Integer vector, column index of the row maximum.}
#'     \item{\code{certainty}}{Row maximum of \code{memberships}.}
#'     \item{\code{certainty_raw}}{Row maximum of \code{weighted}.}
#'     \item{\code{entropy}}{Normalised row entropy of \code{memberships}.}
#'     \item{\code{key}}{Data frame mapping each column to its region and local type.}
#'     \item{\code{n_local_types}}{Selected G per region.}
#'     \item{\code{BIC}}{G x region matrix: the best weighted BIC over
#'       covariance models for each G.}
#'     \item{\code{BIC_by_model}}{G x model x region array of weighted BIC
#'       values; \code{NA} where a candidate was infeasible or failed.}
#'     \item{\code{model_names}}{Selected covariance model per region.}
#'     \item{\code{models}}{Selected fits, or \code{NULL}.}
#'     \item{\code{variables}, \code{cutoff}, \code{weights},
#'       \code{model_name}}{Settings used (\code{model_name} = candidates).}
#'   }
#'
#' @seealso \code{\link{get_regions}}, \code{\link[mclust]{me.weighted}}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200
#' env <- data.frame(ID = seq_len(n),
#'                   temp     = stats::rnorm(n),
#'                   slope    = stats::rnorm(n),
#'                   rainfall = stats::rnorm(n))
#' W <- matrix(stats::runif(n * 2), ncol = 2)
#' W <- W / rowSums(W)
#' res <- get_local_types(env, regions = W, non_value_cols = "ID",
#'                        n_local_types = 2:3)
#' res
#' }
#'
#' @importFrom mclust hc hclass unmap me.weighted priorControl nMclustParams mclust.options hcE hcV hcEII hcVII hcEEE hcVVV
#' @importFrom stats var
#' @export
get_local_types <- function(env,
                            regions,
                            variables      = NULL,
                            non_value_cols = NULL,
                            n_local_types  = 2:9,
                            cutoff         = 0.1,
                            weights        = c("memberships", "typicality"),
                            model_name     = c("VVV", "EEE", "VVI", "VII"),
                            hc_use         = "VARS",
                            drop_constant  = TRUE,
                            weight_eps     = 1e-6,
                            prior          = mclust::priorControl(),
                            keep_models    = TRUE,
                            verbose        = FALSE) {
        
        weights <- match.arg(weights)
        
        ## ---- 0. Input validation -------------------------------------------*
        if (!is.data.frame(env))
                stop("`env` must be a data.frame or sf object.", call. = FALSE)
        if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
                stop("`verbose` must be a single TRUE/FALSE value.", call. = FALSE)
        if (!is.numeric(n_local_types) || length(n_local_types) == 0L ||
            any(n_local_types < 1L) || any(n_local_types != as.integer(n_local_types)))
                stop("`n_local_types` must be one or more positive integers.", call. = FALSE)
        n_local_types <- sort(unique(as.integer(n_local_types)))
        if (!is.null(variables) && !is.character(variables))
                stop("`variables` must be a character vector or NULL.", call. = FALSE)
        ## FIX: `cutoff` was never range-checked, and core must be a subset of the
        ## fitted set or `match(core, act)` below produces NA row indices.
        if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) ||
            cutoff <= 0 || cutoff > 1)
                stop("`cutoff` must be a single number in (0, 1].", call. = FALSE)
        if (!is.numeric(weight_eps) || length(weight_eps) != 1L || weight_eps < 0)
                stop("`weight_eps` must be a single non-negative number.", call. = FALSE)
        if (cutoff <= weight_eps)
                stop("`cutoff` must be greater than `weight_eps`.", call. = FALSE)
        if (!is.character(model_name) || !length(model_name) || anyNA(model_name))
                stop("`model_name` must be a character vector of mclust model names.",
                     call. = FALSE)
        model_name <- unique(model_name)
        
        ## ---- 1. Region weights ---------------------------------------------*
        W         <- .region_weights(regions, weights)
        n_regions <- ncol(W)
        N         <- nrow(W)
        if (N != nrow(env))
                stop(sprintf("Region weights have %d rows but `env` has %d; they must match.",
                             N, nrow(env)), call. = FALSE)
        if (n_regions < 1L) stop("Region weight matrix has no columns.", call. = FALSE)
        
        ## ---- 2. Predictors --------------------------------------------------*
        ## FIX: `inherits()` rather than `"sf" %in% class()`.
        if (inherits(env, "sf")) env <- sf::st_drop_geometry(env)
        
        ## FIX: the old code built `drop_cols <- setdiff(names(env), variables)`, so
        ## the default `variables = NULL` dropped EVERY column and the function
        ## always died on its own "no predictor columns" error. Default now means
        ## "all numeric columns except the identifiers".
        if (is.null(variables)) {
                cand      <- setdiff(names(env), non_value_cols)
                is_num    <- vapply(env[cand], is.numeric, logical(1))
                variables <- cand[is_num]
                if (verbose)
                        message(sprintf("Using %d numeric predictor(s): %s",
                                        length(variables), paste(variables, collapse = ", ")))
        } else {
                ## FIX: a typo'd or absent variable used to be silently dropped.
                missing_v <- setdiff(variables, names(env))
                if (length(missing_v))
                        stop("`variables` not found in `env`: ",
                             paste(missing_v, collapse = ", "), call. = FALSE)
                variables <- setdiff(variables, non_value_cols)
        }
        if (!length(variables))
                stop("No predictor columns left after pre-processing.", call. = FALSE)
        
        Xdata <- env[, variables, drop = FALSE]
        
        ## FIX: `apply(Xdata, 2, is.numeric)` coerced the frame to one matrix first,
        ## so a single character column flagged every predictor and an all-numeric
        ## frame flagged none — the check could never do its job.
        non_num <- !vapply(Xdata, is.numeric, logical(1))
        if (any(non_num))
                stop(sprintf("Predictor(s) %s is/are non-numeric.",
                             paste(names(Xdata)[non_num], collapse = ", ")), call. = FALSE)
        if (anyNA(Xdata)) {
                bad <- names(Xdata)[vapply(Xdata, anyNA, logical(1))]
                stop(sprintf("Missing values in predictor(s) %s; mclust requires complete cases.",
                             paste(bad, collapse = ", ")), call. = FALSE)
        }
        ## OPT: one coercion to a numeric matrix; every downstream subset is then a
        ## cheap matrix slice rather than a data-frame subset.
        Xdata <- as.matrix(Xdata)
        storage.mode(Xdata) <- "double"
        
        ## ---- 3. Cores, active sets, constant predictors ---------------------*
        re   <- seq_len(n_regions)
        act  <- lapply(re, function(r) which(W[, r] > weight_eps))
        core <- lapply(re, function(r) which(W[, r] >= cutoff))
        
        ## FIX: the old empty-core check closed sprintf() before its argument, so it
        ## always failed with "too few arguments" instead of reporting the region.
        empty <- which(lengths(core) == 0L)
        if (length(empty))
                stop(sprintf("Region(s) %s have an empty core: no unit reaches cutoff = %g.",
                             paste(empty, collapse = ", "), cutoff), call. = FALSE)
        
        orphan <- which(rowSums(W > weight_eps) == 0L)
        if (length(orphan))
                warning(sprintf("%d unit(s) have no region weight above `weight_eps`; ",
                                length(orphan)),
                        "their local-type rows will be all zero.", call. = FALSE)
        
        if (drop_constant) {
                keep <- rep(TRUE, ncol(Xdata))
                for (r in re) keep <- keep & (.col_var(Xdata, core[[r]]) > 0) &
                                (.col_var(Xdata, act[[r]])  > 0)
                if (!all(keep)) {
                        if (verbose || sum(!keep) > 0)
                                message("Dropping predictor(s) constant within a region: ",
                                        paste(colnames(Xdata)[!keep], collapse = ", "))
                        Xdata <- Xdata[, keep, drop = FALSE]
                        variables <- colnames(Xdata)
                }
                if (!ncol(Xdata))
                        stop("Every predictor is constant within at least one region.",
                             call. = FALSE)
        }
        d      <- ncol(Xdata)
        models <- .resolve_model_names(model_name, d)
        
        ## ---- 4. Per-region fitting -----------------------------------------*
        ## OPT: the old code looped G on the outside and regions on the inside, so it
        ## held one full-N z matrix PER REGION simultaneously and refitted every
        ## region on all N rows for every candidate G. Here the region is the outer
        ## loop: the hc initialisation is computed once per region, the fit runs on
        ## the active rows only, and just the incumbent-best z is retained.
        BIC_all <- array(NA_real_,
                         dim = c(length(n_local_types), length(models), n_regions),
                         dimnames = list(G = as.character(n_local_types), model = models,
                                         region = paste0("Region", re)))
        best <- vector("list", n_regions)
        
        for (r in re) {
                a_r  <- act[[r]]
                c_r  <- core[[r]]
                w_r  <- W[a_r, r]
                X_r  <- Xdata[a_r, , drop = FALSE]
                ## core is a subset of act because cutoff > weight_eps (checked above).
                c_in_a <- match(c_r, a_r)
                n_eff  <- sum(w_r)
                
                ## Candidates a region can actually support.
                G_r <- n_local_types[n_local_types <= length(c_r) - 1L]
                if (!length(G_r))
                        stop(sprintf("Region %d supports no candidate in `n_local_types` (core has %d units).",
                                     r, length(c_r)), call. = FALSE)
                
                ## Initialisation: full-covariance hc needs more core units than
                ## predictors; below that, spherical hc is always defined. One hc
                ## per region, shared by every candidate model.
                hc_model <- if (length(c_r) > d) "VVV" else "EII"
                hcf <- tryCatch(mclust::hc(Xdata[c_r, , drop = FALSE],
                                           modelName = hc_model, use = hc_use),
                                error = function(e) e)
                if (inherits(hcf, "error") && hc_model != "EII") {
                        hc_model <- "EII"
                        hcf <- tryCatch(mclust::hc(Xdata[c_r, , drop = FALSE],
                                                   modelName = hc_model, use = hc_use),
                                        error = function(e) e)
                }
                if (inherits(hcf, "error"))
                        stop(sprintf("mclust::hc() failed for region %d: %s",
                                     r, conditionMessage(hcf)), call. = FALSE)
                cls <- mclust::hclass(hcf, G_r)          # length(c_r) x length(G_r)
                
                if (verbose)
                        message(sprintf("Region %d/%d: %d core, %d fitted, n_eff = %.1f, hc = %s, G in {%s}",
                                        r, n_regions, length(c_r), length(a_r), n_eff,
                                        hc_model, paste(G_r, collapse = ", ")))
                
                for (j in seq_along(G_r)) {
                        G  <- G_r[j]
                        z0 <- matrix(1 / G, nrow = length(a_r), ncol = G)
                        z0[c_in_a, ] <- mclust::unmap(cls[, j], groups = seq_len(G))
                        
                        for (m in models) {
                                ## Feasibility: BIC is meaningless with more parameters
                                ## than effective observations. This is the fallback
                                ## that lets small regions use simpler covariances.
                                npar <- mclust::nMclustParams(modelName = m, d = d, G = G)
                                if (npar >= n_eff) next
                                
                                fit <- tryCatch(
                                        mclust::me.weighted(data = X_r, modelName = m,
                                                            z = z0, weights = w_r, prior = prior),
                                        error = function(e) NULL)
                                if (is.null(fit) || is.null(fit$loglik) || !is.finite(fit$loglik))
                                        next
                                
                                bic <- .weighted_bic(fit$loglik, m, d, G, n_eff)
                                BIC_all[as.character(G), m, r] <- bic
                                
                                ## Keep only the incumbent's identity; its z is
                                ## recovered by one refit below.
                                if (is.finite(bic) &&
                                    (is.null(best[[r]]) || bic > best[[r]]$bic))
                                        best[[r]] <- list(bic = bic, G = G, j = j, model = m,
                                                          act = a_r, w = w_r)
                        }
                        if (verbose) {
                                b <- BIC_all[as.character(G), , r]
                                if (all(is.na(b))) {
                                        message(sprintf("    G = %d: no feasible model", G))
                                } else {
                                        k <- which.max(b)
                                        message(sprintf("    G = %d: best %s, BIC = %.1f (%d/%d models fitted)",
                                                        G, models[k], b[k], sum(!is.na(b)), length(b)))
                                }
                        }
                }
                
                if (is.null(best[[r]]))
                        stop(sprintf(paste0("Region %d: no (model, G) candidate could be fitted ",
                                            "(n_eff = %.1f, %d predictors, smallest candidate needs %d ",
                                            "parameters). Reduce `variables` or lower `cutoff`."),
                                     r, n_eff, d,
                                     min(vapply(models, function(m)
                                             mclust::nMclustParams(m, d, min(G_r)), numeric(1)))),
                             call. = FALSE)
                
                ## Refit the winner once to recover its z without storing every
                ## candidate's z: O(1) rather than O(|G| x |models|) z matrices.
                G  <- best[[r]]$G
                z0 <- matrix(1 / G, nrow = length(a_r), ncol = G)
                z0[c_in_a, ] <- mclust::unmap(cls[, best[[r]]$j], groups = seq_len(G))
                fit <- mclust::me.weighted(data = X_r, modelName = best[[r]]$model,
                                           z = z0, weights = w_r, prior = prior)
                best[[r]]$z     <- fit$z
                fit$z           <- NULL                  # stored model without the bulk
                best[[r]]$fit   <- fit
        }
        
        ## G x region summary (best over models), so existing BIC plots still work.
        BIC <- apply(BIC_all, c(1L, 3L), function(v)
                if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
        BIC <- matrix(BIC, nrow = length(n_local_types),
                      dimnames = list(as.character(n_local_types), paste0("Region", re)))
        
        ## ---- 5. Assemble ----------------------------------------------------*
        Gs     <- vapply(best, function(b) b$G, integer(1))
        offset <- c(0L, cumsum(Gs))
        total  <- offset[length(offset)]
        
        Zw <- matrix(0, nrow = N, ncol = total)
        for (r in re) {
                cols <- (offset[r] + 1L):offset[r + 1L]
                Zw[best[[r]]$act, cols] <- best[[r]]$z * best[[r]]$w
        }
        colnames(Zw) <- unlist(lapply(re, function(r)
                sprintf("local%d_%d", r, seq_len(Gs[r]))), use.names = FALSE)
        
        key <- data.frame(column     = seq_len(total),
                          region     = rep(re, times = Gs),
                          local_type = unlist(lapply(Gs, seq_len), use.names = FALSE),
                          name       = colnames(Zw),
                          stringsAsFactors = FALSE)
        
        rs <- rowSums(Zw)
        U  <- Zw / pmax(rs, .Machine$double.eps)
        
        ## FIX / OPT: max.col() once, instead of three separate apply(., 1, which.max)
        ## passes over an N x total matrix (and one whose result was discarded).
        hard <- max.col(U, ties.method = "first")
        idx  <- cbind(seq_len(N), hard)
        
        out <- list(
                memberships   = U,
                weighted      = Zw,
                hard_classes  = as.integer(hard),
                certainty     = U[idx],
                certainty_raw = Zw[idx],
                entropy       = .row_entropy(U),
                key           = key,
                n_local_types = stats::setNames(Gs, paste0("Region", re)),
                BIC           = BIC,
                BIC_by_model  = BIC_all,
                model_names   = stats::setNames(vapply(best, `[[`, character(1), "model"),
                                                paste0("Region", re)),
                models        = if (keep_models)
                        stats::setNames(lapply(best, `[[`, "fit"), paste0("Region", re))
                else NULL,
                variables     = variables,
                cutoff        = cutoff,
                weights       = weights,
                model_name    = models
        )
        structure(out, class = "pulse_local_types")
}


#' @rdname get_local_types
#' @param x An object of class \code{"pulse_local_types"}.
#' @param ... Unused.
#' @return \code{print.pulse_local_types} returns \code{x} invisibly.
#' @export
print.pulse_local_types <- function(x, ...) {
        cat(sprintf("<pulse_local_types>  %d units x %d local types across %d regions\n",
                    nrow(x$memberships), ncol(x$memberships), length(x$n_local_types)))
        cat(sprintf("  types per region: %s\n",
                    paste(x$n_local_types, collapse = ", ")))
        cat(sprintf("  weights = %s | cutoff = %g | %d predictor(s)\n",
                    x$weights, x$cutoff, length(x$variables)))
        sel <- if (is.null(x$model_names)) x$model_name else x$model_names
        tab <- table(factor(sel, levels = unique(sel)))
        cat(sprintf("  covariance model per region: %s\n",
                    paste(sprintf("%s x%d", names(tab), tab), collapse = ", ")))
        cat(sprintf("  mean certainty = %.3f | mean normalised entropy = %.3f\n",
                    mean(x$certainty), mean(x$entropy)))
        invisible(x)
}
#' Collapse local types across regions into meta-types
#'
#' Cuts the cross-region meta-clustering (see [plot_local_similarity()]) and
#' merges local types that fall in the same meta-cluster by summing their
#' membership columns. Meta-types are labelled sequentially (`meta.1`, `meta.2`,
#' ...) in increasing order of their meta-cluster id.
#'
#' Run this after inspecting [plot_local_similarity()] to commit to a particular
#' number of meta-types (`k`) or a dendrogram cut height (`h`).
#'
#' @param local A fitted local-type object (as returned by `get_local_types()`).
#'   Only used to (re)compute `similarity` when it is not supplied.
#' @param memberships The wide local-type membership matrix (catchments x types).
#'   Column names must resolve to the internal meta labels (`R<region>\u00b7T<type>`);
#'   columns named `Region<n>localType<g>` are relabelled to that form
#'   automatically.
#' @param k,h Cut the meta-clustering into `k` meta-types or at dendrogram height
#'   `h`. Supply exactly one of the two; supplying neither is an error.
#' @param similarity Optional precomputed `localtype_similarity` (from
#'   `.localtype_meta()`); recomputed from `local` when `NULL`.
#' @param distance Distance used only when `similarity` is recomputed: `"mean"`
#'   (component means), `"hellinger"` (mean and shape), or `"wasserstein"`.
#' @param sep Separator assumed between region and type in the membership column
#'   labels; used when building the mismatch error message.
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return A list with components:
#'   \describe{
#'     \item{memberships}{Catchments x meta-types matrix of summed memberships.
#'       Rows are *not* renormalised, so they need not sum to 1.}
#'     \item{hard}{Hard meta-type label per catchment (argmax of the normalised
#'       memberships).}
#'     \item{entropy}{Shannon entropy of each catchment's normalised membership
#'       vector; higher values flag more ambiguous assignments.}
#'     \item{row_sums}{Row sums of `memberships` before normalisation.}
#'     \item{mapping}{Data frame linking each original type (`region`, `type`,
#'       `label`) to its meta-cluster (`meta`) and label (`meta_label`).}
#'     \item{n_before,n_after}{Counts of local types before, and meta-types
#'       after, collapsing.}
#'     \item{k,h}{The cut arguments that were used.}
#'   }
#' @seealso [plot_local_similarity()], [plot_local_embedding()]
#' @export
collapse_local <- function(local, memberships, k = NULL, h = NULL, similarity = NULL,
                           distance = c("mean", "hellinger", "wasserstein"),
                           sep = ".", ...) {
        distance <- match.arg(distance)
        if (is.null(k) && is.null(h)) stop("supply k or h.")
        # Build the cross-region similarity / meta-clustering unless one was passed in
        if (is.null(similarity)) similarity <- .localtype_meta(local, distance = distance, sep = sep)
        
        # Cut the meta-dendrogram: grp is a named integer vector (names = type labels,
        # value = meta-cluster id)
        grp  <-
                if (!is.null(k)){
                        stats::cutree(similarity$hclust, k = k)
                } else {
                        stats::cutree(similarity$hclust, h = h)
                }
        info <- similarity$info
        info$meta <- grp[info$label]                   # attach each label's meta-cluster id
        
        # name each meta-type sequentially: meta.1, meta.2, ... in ascending id order
        ids <- sort(unique(info$meta))
        med <- paste("meta", seq_along(ids), sep = ".")
        names(med) <- ids
        info$meta_label <- med[as.character(info$meta)]
        
        # Align membership column names to the internal labels: e.g. "Region3localType2"
        # becomes "R3\u00b7T2". Columns already in that form are left untouched by sub().
        colnames(memberships) <-
                sub(pattern = "^Region([0-9]+)localType([0-9]+)$", replacement = "R\\1\u00b7T\\2", x = colnames(memberships))
        
        # NOTE: the labels use a middle dot ("\u00b7"); `sep` only formats the message below.
        if (!all(info$label %in% colnames(memberships)))
                stop("membership columns must match the meta labels (region", sep, "type). ",
                     "Confirm the get_local_types() column naming and pass a matching matrix or relabel.")
        
        # Restrict and order the membership matrix to the known labels, then map each
        # column to the meta-type it belongs to
        M   <- as.matrix(memberships[, info$label, drop = FALSE])
        key <- info$meta_label[match(colnames(M), info$label)]
        
        # Sum the member columns within each meta-type -> catchments x meta-types
        collapsed <- vapply(split(seq_len(ncol(M)), key),
                            function(j) rowSums(M[, j, drop = FALSE]), numeric(nrow(M)))
        
        collapsed <- collapsed[, unique(key), drop = FALSE]  # restore first-appearance order
        # rows are not guaranteed to sum to 1
        rs   <- rowSums(collapsed)
        # normalise only for hard label and entropy (the returned memberships stay raw)
        prob <- collapsed / ifelse(rs > 0, rs, 1)
        hard <- colnames(collapsed)[max.col(prob, ties.method = "first")]  # argmax meta-type
        ent  <- -rowSums(ifelse(prob > 0, prob * log(prob), 0))            # Shannon entropy
        
        list(memberships = collapsed, hard = hard, entropy = ent, row_sums = rs,
             mapping = info[order(info$meta), c("region", "type", "label", "meta", "meta_label")],
             n_before = nrow(info), n_after = ncol(collapsed), k = k, h = h)
}

#' Plot cross-region similarity of local types
#'
#' Compares the GMM components of every region's local types in a common,
#' standardised predictor space and draws a meta-clustered distance heatmap.
#' Off-diagonal blocks that mix regions flag local types that are candidates for
#' collapsing with [collapse_local()].
#'
#' @param local A fitted local-type object (uses `local$models`).
#' @param distance Component distance: `"mean"` (standardised component means,
#'   matches the fingerprint plot), `"hellinger"` (mean and shape), or
#'   `"wasserstein"`.
#' @param k,h Optional meta-clustering cut used *only* to colour the axis labels:
#'   `k` meta-types or dendrogram height `h`.
#' @param groups Optional grouping that overrides `k`/`h` for label colouring.
#'   Either a `collapse_local()` result (its `$mapping` is used) or a named
#'   vector mapping each label to a group.
#' @param triangle If `TRUE` (default) draw only the upper triangle; otherwise
#'   the full off-diagonal matrix. The diagonal is always dropped.
#' @param ... Passed to the internal `.localtype_meta()` (e.g. `scale`, `ridge`).
#'
#' @return A `ggplot` object: the cross-region similarity heatmap.
#' @seealso [collapse_local()], [plot_local_embedding()]
#' @export
plot_local_similarity <- function(local, distance = c("mean", "hellinger", "wasserstein"),
                                  k = NULL, h = NULL, groups = NULL, triangle = TRUE, ...) {
        .region_require("ggplot2")
        distance <- match.arg(distance)
        meta <- .localtype_meta(local, distance = distance, ...)
        ord  <- meta$hclust$order                      # dendrogram order, so blocks read along the diagonal
        labs <- meta$info$label[ord]
        n    <- length(labs)
        M    <- meta$matrix[ord, ord]
        
        ## reshape to long form; keep upper triangle, diagonal always dropped
        gi   <- expand.grid(i = seq_len(n), j = seq_len(n))
        gi$d <- M[cbind(gi$i, gi$j)]
        gi   <- if (triangle) gi[gi$j > gi$i, ] else gi[gi$i != gi$j, ]
        gi$a <- factor(labs[gi$i], levels = labs)
        gi$b <- factor(labs[gi$j], levels = labs)
        
        ## groups: a collapse_local() result, a named vector, or an internal cut
        grp <- if (is.list(groups) && !is.null(groups$mapping))
                stats::setNames(groups$mapping$meta_label, groups$mapping$label)[labs]
        else if (!is.null(groups)) groups[labs]
        else if (!is.null(k))      stats::cutree(meta$hclust, k = k)[labs]
        else if (!is.null(h))      stats::cutree(meta$hclust, h = h)[labs]
        else                       NULL
        
        ## colour each label by its group; merged types share a colour, singletons stay grey
        lab_col <- rep("grey35", n)
        if (!is.null(grp)) {
                tab   <- table(grp)
                multi <- names(tab)[tab >= 2]          # only groups with >= 2 members get a colour
                if (length(multi)) {
                        pal <- stats::setNames(grDevices::hcl.colors(length(multi), "Dark 3"), multi)
                        hit <- !is.na(grp) & grp %in% multi
                        lab_col[hit] <- pal[as.character(grp[hit])]
                }
        }
        
        ggplot2::ggplot(gi, ggplot2::aes(.data$a, .data$b, fill = .data$d)) +
                ggplot2::geom_tile() +
                ggplot2::scale_fill_gradient(low = "#08306B", high = "grey92", name = "distance") +
                ggplot2::scale_x_discrete(drop = FALSE) +
                ggplot2::scale_y_discrete(drop = FALSE) +
                ggplot2::coord_equal() +
                ggplot2::labs(title = "Cross-region local-type similarity",
                              subtitle = sprintf("%s distance between GMM components; shared label colour marks a collapse group", distance),
                              x = NULL, y = NULL) +
                theme_eco() +
                ggplot2::theme(panel.grid  = ggplot2::element_blank(),
                               axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                                                   colour = lab_col),
                               axis.text.y = ggplot2::element_text(colour = lab_col))
}

#' Per-component covariance matrices from an mclust fit
#'
#' Normalises the several shapes `mclust` uses to store variances (per-component
#' `sigma`, shared `Sigma`, or univariate `sigmasq`) into a list of `G` labelled
#' covariance matrices. Any `NULL` element in the result is the caller's signal
#' to fall back to a mean-only distance.
#'
#' @param fit A fitted `mclust` model for one region's local types.
#' @param vn Character vector of variable names used to label the matrices.
#' @return List of length `G` of d x d covariance matrices (entries may be NULL).
#' @noRd
.gmm_sigma <- function(fit, vn) {
        v <- fit$parameters$variance
        G <- if (is.null(dim(fit$parameters$mean))) {  # univariate: mean is a length-G vector
                length(fit$parameters$mean)
        }  else {
                ncol(fit$parameters$mean)
        }
        out <- vector("list", G)
        if (!is.null(v$sigma)) {                       # d x d x G (most multivariate models)
                for (g in seq_len(G)) { S <- v$sigma[, , g]; dimnames(S) <- list(vn, vn); out[[g]] <- S }
        } else if (!is.null(v$Sigma)) {                # shared covariance across components
                S <- v$Sigma; dimnames(S) <- list(vn, vn); for (g in seq_len(G)) out[[g]] <- S
        } else if (!is.null(v$sigmasq)) {              # univariate variance(s)
                s2 <- v$sigmasq; if (length(s2) == 1L) s2 <- rep(s2, G)
                for (g in seq_len(G)) out[[g]] <- matrix(s2[g], 1, 1, dimnames = list(vn, vn))
        }
        out                                            # any NULL element triggers the mean-only fallback
}

#' Symmetric matrix square root via eigendecomposition
#'
#' Returns the symmetric PSD square root of `S`, clamping tiny negative
#' eigenvalues to 0 for numerical safety.
#'
#' @param S A symmetric (covariance) matrix.
#' @return The symmetric square-root matrix of `S`.
#' @noRd
.mat_sqrt <- function(S) {
        e <- eigen(S, symmetric = TRUE); e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(S)) %*% t(e$vectors)
}

#' Distance between two GMM components
#'
#' Computes the `"mean"` (Euclidean between standardised means), `"wasserstein"`
#' (closed-form Gaussian W2), or `"hellinger"` (via the Bhattacharyya
#' coefficient) distance. Falls back to the mean distance whenever either
#' covariance is missing. A ridge proportional to the mean variance stabilises
#' the matrix inverse / square root.
#'
#' @param a,b Component descriptors with standardised mean `$m` and optional
#'   standardised covariance `$S`.
#' @param distance One of `"mean"`, `"hellinger"`, `"wasserstein"`.
#' @param ridge Relative ridge added to the covariance diagonals.
#' @return A scalar distance.
#' @noRd
.pair_dist <- function(a, b, distance, ridge) {
        dm <- a$m - b$m
        # mean distance, or fall back to it when either covariance is unavailable
        if (distance == "mean" || is.null(a$S) || is.null(b$S)) return(sqrt(sum(dm^2)))
        Sa <- a$S; Sb <- b$S
        rg <- ridge * mean(c(diag(Sa), diag(Sb)))      # ridge scaled to the typical variance
        Sa <- Sa + diag(rg, nrow(Sa)); Sb <- Sb + diag(rg, nrow(Sb))
        if (distance == "wasserstein") {               # closed-form 2-Wasserstein between Gaussians
                ra    <- .mat_sqrt(Sa)
                inner <- .mat_sqrt(ra %*% Sb %*% ra)
                return(sqrt(max(0, sum(dm^2) + sum(diag(Sa + Sb - 2 * inner)))))
        }
        Sbar <- (Sa + Sb) / 2                          # Hellinger via Bhattacharyya
        quad <- as.numeric(t(dm) %*% solve(Sbar, dm)) / 8                  # mean-separation term
        ld   <- as.numeric(determinant(Sbar, TRUE)$modulus) -             # covariance-mismatch term
                0.5 * (as.numeric(determinant(Sa, TRUE)$modulus) + as.numeric(determinant(Sb, TRUE)$modulus))
        sqrt(max(0, 1 - exp(-(quad + 0.5 * ld))))
}

#' Cross-region meta-clustering of local types
#'
#' Pools every region's GMM components into one standardised predictor space,
#' builds the pairwise component-distance matrix, and average-links them into a
#' meta-clustering. This is the engine behind both [plot_local_similarity()] and
#' [collapse_local()].
#'
#' @param lt A fitted local-type object (uses `lt$models`).
#' @param distance Component distance; see [plot_local_similarity()].
#' @param scale If `TRUE`, standardise predictors by their pooled SD so that no
#'   single variable dominates the distance.
#' @param sep Reserved separator for region/type labels (labels currently use a
#'   middle dot regardless of this value).
#' @param ridge Relative ridge forwarded to `.pair_dist()`.
#' @return A `localtype_similarity` object: the `dist` and full `matrix`, the
#'   `hclust` tree, per-label `info`, the common `vars`, and the `distance` used.
#' @noRd
.localtype_meta <- function(lt, distance = c("mean", "hellinger", "wasserstein"),
                            scale = TRUE, sep = ".", ridge = 1e-6) {
        distance <- match.arg(distance)
        
        comps  <- list()                               # one entry per (region, type) component
        pooled <- list()                               # per-region predictor tables, for SD scaling
        # Loop over Regions
        for (nm in names(lt$models)) {
                # Extract regional model
                fit <- lt$models[[nm]]
                # Means for each river type
                mu  <- fit$parameters$mean
                if (is.null(dim(mu))) {                # univariate fit: promote to a 1-row matrix
                        vn <- colnames(fit$data); if (is.null(vn)) vn <- "var"
                        mu <- matrix(mu, nrow = 1L, dimnames = list(vn[1], NULL))
                }
                # Variable names
                vn <- rownames(mu)
                # Add replacement V1...Vn if no rownames are given
                if (is.null(vn)) {
                        vn <- paste0("V", seq_len(nrow(mu)))
                        rownames(mu) <- vn
                }
                # run local helper which extracts covariance matrix
                sig <- .gmm_sigma(fit, vn)
                # register each component under a "R<region>\u00b7T<type>" label
                for (g in seq_len(ncol(mu))){
                        comps[[paste0(sub("Region","R", nm),"\u00b7", "T", g)]] <-
                                list(region = nm, type = g, mean = mu[, g], sigma = sig[[g]], vars = vn)
                }
                
                if (!is.null(fit$data)) pooled[[nm]] <- as.data.frame(fit$data)
        }
        if (length(comps) < 2L) stop("need at least two local types to compare.")
        
        vars <- Reduce(intersect, lapply(comps, `[[`, "vars"))   # common predictor space
        if (!length(vars)) stop("regional models share no common predictors.")
        
        # Per-variable scaling factor: pooled SD when scaling, else 1 (no scaling)
        sd_v <- stats::setNames(rep(1, length(vars)), vars)
        if (scale && length(pooled)) {
                pooled  <- do.call(rbind, lapply(pooled, function(d) d[, vars, drop = FALSE]))
                sd_v[]  <- vapply(vars, function(v) stats::sd(pooled[[v]], na.rm = TRUE), numeric(1))
                sd_v[!is.finite(sd_v) | sd_v <= 0] <- 1  # guard against zero / NA SDs
        }
        D <- diag(1 / sd_v, length(vars))              # diagonal standardising transform
        
        # Standardise each component's mean (and covariance) into the common space
        for (lab in names(comps)) {
                comps[[lab]]$m <- as.numeric(D %*% comps[[lab]]$mean[vars])
                S <- comps[[lab]]$sigma
                comps[[lab]]$S <- if (is.null(S)) NULL else D %*% S[vars, vars, drop = FALSE] %*% D
        }
        
        # Fill the symmetric pairwise distance matrix
        labs <- names(comps)
        P <- length(labs)
        M <- matrix(0, P, P, dimnames = list(labs, labs))
        for (a in seq_len(P - 1L)) for (b in (a + 1L):P){
                M[a, b] <- M[b, a] <- .pair_dist(comps[[a]],
                                                 comps[[b]],
                                                 distance,
                                                 ridge)
        }
        
        
        info <- data.frame(label  = labs,
                           region = vapply(comps, `[[`, "", "region"),
                           type   = vapply(comps, `[[`, 0L, "type"),
                           row.names = NULL, stringsAsFactors = FALSE)
        structure(list(dist = stats::as.dist(M), matrix = M,
                       hclust = stats::hclust(stats::as.dist(M), method = "average"),
                       info = info, vars = vars, distance = distance),
                  class = "localtype_similarity")
}

#' Plot local types in a reduced environmental space
#'
#' PCA-embeds all samples, then overlays original-type centroids and, optionally,
#' the collapsed-type centroids from [collapse_local()]. Useful for checking
#' whether a proposed collapse merges types that genuinely overlap in predictor
#' space.
#'
#' @param local Fitted local-type object; supplies the predictor set via
#'   `local$models`.
#' @param data All-sample predictor table, row-aligned to `memberships`.
#' @param memberships Wide local-type membership matrix; the argmax of each row
#'   gives that sample's original type.
#' @param collapse Optional `collapse_local()` result; colours samples by merged
#'   type and adds the collapsed centroids (diamonds).
#' @param dims 2 (static ggplot2) or 3 (interactive plotly).
#' @param n_max Maximum number of samples drawn per original type in the scatter.
#' @param scale If `TRUE`, scale predictors to unit variance before the PCA.
#' @param seed Optional RNG seed for the per-type sample thinning.
#' @param samples If `TRUE`, draw the sample points; if `FALSE` they are made
#'   fully transparent so only the centroids show (2D only).
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return For `dims = 2` a `ggplot` object; for `dims = 3` a `plotly` object.
#' @seealso [collapse_local()], [plot_local_similarity()]
#' @export
plot_local_embedding <- function(local, data, memberships, collapse = NULL,
                                 dims = 2, n_max = 100, scale = TRUE, seed = NULL, samples = TRUE, ...) {
        .region_require("ggplot2")
        if (!dims %in% c(2L, 3L)) stop("dims must be 2 or 3.")
        if (missing(data) || missing(memberships))
                stop("supply `data` (all-sample predictors) and `memberships` (wide river-type matrix).")
        memberships <- as.matrix(memberships)
        data <- as.data.frame(data)
        if (nrow(memberships) != nrow(data)) stop("`data` and `memberships` need the same rows.")
        
        # Predictor set: the model's variables, intersected with what `data` provides
        mvars <- rownames(local$models[[1]]$parameters$mean)
        if (is.null(mvars)) mvars <- names(Filter(is.numeric, data))
        mvars <- intersect(mvars, colnames(data))
        X <- as.matrix(data[, mvars, drop = FALSE])
        X <- X[, apply(X, 2, stats::sd) > 0, drop = FALSE]   # drop constant columns (PCA can't use them)
        
        # Original (pre-collapse) type per sample = argmax membership column
        orig <- colnames(memberships)[max.col(memberships, ties.method = "first")]
        coll <- NULL
        if (!is.null(collapse)) {
                # Map each original type to its collapsed (meta) type
                map <- stats::setNames(collapse$mapping$meta_label, collapse$mapping$label)
                if (!all(orig %in% names(map)))
                        warning("some membership columns are not in the collapse mapping; check the column naming.")
                coll <- unname(map[orig])
        }
        grp <- if (is.null(coll)) orig else coll       # colour samples by collapsed type if available
        
        ## region for each original-type label: from the mapping, else from per-region block sizes
        col_region <- if (!is.null(collapse)) {
                stats::setNames(collapse$mapping$region, collapse$mapping$label)
        } else {
                # rebuild region-per-column from the number of components in each region's model
                gper <- vapply(local$models, function(f) { mu <- f$parameters$mean
                if (is.null(dim(mu))) length(mu) else ncol(mu) }, integer(1))
                rcv  <- rep(names(local$models), gper)
                if (length(rcv) == ncol(memberships)) stats::setNames(rcv, colnames(memberships)) else NULL
        }
        have_region <- !is.null(col_region)
        if (!have_region) warning("could not resolve region per type; crosses drawn in one colour.")
        
        # compute PCA on the (optionally scaled) predictor matrix
        pc <- stats::prcomp(X, center = TRUE, scale. = scale)
        sc <- as.data.frame(pc$x[, seq_len(dims), drop = FALSE]); names(sc) <- paste0("PC", seq_len(dims))
        
        # helper: per-label centroid in PC space (ignores NA labels)
        cen <- function(lab) {
                ok <- !is.na(lab); stats::aggregate(sc[ok, , drop = FALSE], list(label = lab[ok]), mean)
        }
        cen_o <- cen(orig)                             # original-type centroids
        if (have_region) cen_o$region <- unname(col_region[cen_o$label])
        cen_c <- if (!is.null(coll)) cen(coll) else NULL   # collapsed-type centroids (optional)
        
        # Thin to at most n_max samples per original type so dense types don't dominate
        if (!is.null(seed)) set.seed(seed)
        pick <- unlist(lapply(split(seq_len(nrow(sc)), orig),
                              function(ix) if (length(ix) > n_max) sample(ix, n_max) else ix),
                       use.names = FALSE)
        pts <- cbind(sc[pick, , drop = FALSE], group = grp[pick])
        
        if (dims == 2L) {
                # `samples = FALSE` keeps the points but makes them invisible (centroids only)
                if (samples == TRUE){
                        pnt.alpha = 0.5
                } else {
                        pnt.alpha = 0
                }
                p <- ggplot2::ggplot(pts, ggplot2::aes(.data$PC1, .data$PC2)) +
                        ggplot2::geom_point(ggplot2::aes(fill = .data$group), shape = 21, stroke = 0,
                                            size = 1.6, alpha = pnt.alpha)
                # crosses = original centroids, coloured by region when known
                p <- p + if (have_region)
                        ggplot2::geom_point(data = cen_o, ggplot2::aes(colour = .data$region),
                                            shape = 4, size = 2.6, stroke = 1)
                else
                        ggplot2::geom_point(data = cen_o, shape = 4, size = 2.6, stroke = 1, colour = "grey20")
                p <- p +
                        ggplot2::geom_text(data = cen_o, ggplot2::aes(label = .data$label),
                                           size = 2.6, colour = "grey20", vjust = -0.8, check_overlap = TRUE) +
                        ggplot2::labs(title = "Local types in reduced environmental space",
                                      subtitle = sprintf("PCA of %d predictors; crosses are original centroids%s",
                                                         ncol(X),
                                                         if (!is.null(cen_c)) ", diamonds are collapsed centroids" else ""),
                                      fill = if (is.null(coll)) "type" else "collapsed type", colour = "region") +
                        theme_eco()
                # diamonds = collapsed centroids, overlaid on top
                if (!is.null(cen_c))
                        p <- p +
                        ggplot2::geom_point(data = cen_c, shape = 18, size = 4, colour = "black") +
                        ggplot2::geom_text(data = cen_c, ggplot2::aes(label = .data$label),
                                           size = 3, fontface = "bold", vjust = 1.6, check_overlap = TRUE)
                return(p)
        }
        
        ## 3D (interactive)
        .region_require("plotly")
        fig <- plotly::add_markers(plotly::plot_ly(), data = pts, x = ~PC1, y = ~PC2, z = ~PC3,
                                   color = ~group, marker = list(size = 3, opacity = 0.5), name = "samples")
        fig <- if (have_region)
                plotly::add_markers(fig, data = cen_o, x = ~PC1, y = ~PC2, z = ~PC3, color = ~region,
                                    text = ~label, marker = list(size = 5, symbol = "x"), name = "centroids")
        else
                plotly::add_markers(fig, data = cen_o, x = ~PC1, y = ~PC2, z = ~PC3, text = ~label,
                                    marker = list(size = 5, symbol = "x", color = "grey30"), name = "centroids")
        if (!is.null(cen_c))
                fig <- plotly::add_markers(fig, data = cen_c, x = ~PC1, y = ~PC2, z = ~PC3, text = ~label,
                                           marker = list(size = 8, symbol = "diamond", color = "black"),
                                           name = "collapsed centroids")
        plotly::layout(fig, title = "Local types in reduced environmental space",
                       scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"),
                                    zaxis = list(title = "PC3")))
}