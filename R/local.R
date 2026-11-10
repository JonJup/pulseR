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
#' Accepts a "pulse_regions" object, a "pulse_regions_tuning" object (whose fit lives
#' in $best_result), or a bare numeric matrix. `which = "cohesion"` reads
#' the fit's cohesion matrix.
#' The same checks (numeric, no NA, non-negative) run for every input type.
#'
#' @keywords internal
#' @noRd
.region_weights <- function(regions, which = c("memberships", "cohesion")) {
        which <- match.arg(which)
        
        bare <- is.matrix(regions) || is.data.frame(regions)
        if (bare) {
                W <- as.matrix(regions)
        } else {
                W <- .region_weights_from_object(regions, which)
        }
        
        if (!is.numeric(W))
                stop("Region weights must be numeric.", call. = FALSE)
        if (anyNA(W))
                stop(sprintf(paste0(
                        "Region weights contain NA in %d row(s); Level 2 handling for them ",
                        "is not defined. %s"),
                        sum(rowSums(is.na(W)) > 0),
                        if (bare) "Remove or fill those rows before calling get_local_types()."
                        else paste0("Current fits fill the rows of units with pi_support == 0 ",
                                    "from graph neighbours; refit the regions with get_regions().")),
                        call. = FALSE)
        if (any(W < 0))
                stop("Region weights must be non-negative.", call. = FALSE)
        W
}

#' Pull the region weight matrix out of a "pulse_regions" or "pulse_regions_tuning".
#' @keywords internal
#' @noRd
.region_weights_from_object <- function(regions, which) {
        if (!is.list(regions))
                stop("`regions` must be a 'pulse_regions' / 'pulse_regions_tuning' object or a matrix.",
                     call. = FALSE)
        src <- if (inherits(regions, "pulse_regions_tuning")) regions$best_result else regions
        if (is.null(src))
                stop("`regions` looks like a tuning object but carries no `$best_result`.",
                     call. = FALSE)
        W <- if (which == "memberships") .first_non_null(src$memberships, src$fit$memberships) else
                .first_non_null(src$cohesion, src$fit$cohesion)
        if (is.null(W))
                stop(sprintf("Could not find a `%s` matrix on `regions`.", which), call. = FALSE)
        as.matrix(W)
}

#' Level 1 hard region per unit: the regions object's own labels when present
#' (post contiguity repair), else the row argmax of the weight matrix.
#' @noRd
.region_hard <- function(regions, W) {
        src <- if (inherits(regions, "pulse_regions_tuning")) regions$best_result else regions
        h <- if (is.list(src) && !is.data.frame(src))
                .first_non_null(src$hard_labels, src$fit$hard_labels) else NULL
        if (is.null(h)) return(max.col(W, ties.method = "first"))
        h <- as.integer(h)
        if (length(h) != nrow(W) || anyNA(h) || any(h < 1L | h > ncol(W)))
                stop("`regions$hard_labels` must be integers in 1..ncol(memberships), ",
                     "one per unit.", call. = FALSE)
        h
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

#' mclust default prior with its scale taken from all units
#'
#' Mirrors \code{mclust::defaultPrior()} but computes \code{scale} from the
#' covariance \code{V} of all units, so a predictor that is constant within
#' one region still gets a positive prior variance there. A user-supplied
#' \code{scale} is kept.
#'
#' @param prior A \code{priorControl()} list for \code{"defaultPrior"}.
#' @param V Covariance matrix of all units (d x d).
#' @param model mclust model name.
#' @param G Number of components.
#' @return \code{prior} with \code{scale} set.
#' @keywords internal
#' @noRd
.global_scale_prior <- function(prior, V, model, G) {
        if (!is.null(prior$scale)) return(prior)
        d <- ncol(V)
        ## EXPR named: the first case label `E` would otherwise partially
        ## match switch()'s EXPR argument (R CMD check NOTE).
        prior$scale <- switch(EXPR = model,
                E = , V = V[1, 1] / G^2,
                EII = , VII = , EEI = , EVI = , VEI = , VVI = (1 / G)^(2 / d) * sum(diag(V)) / d,
                (1 / G)^(2 / d) * V)
        prior
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
#' in (or cohesion for) that region. Units with \eqn{w_{ir} \ge}
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
#' The hierarchical initialization uses \code{"VVV"} when a region's core has
#' more units than predictors and \code{"EII"} otherwise. It is computed once
#' per region and shared by every candidate model, so all candidates start
#' from the same partition.
#'
#' The membership of local type \eqn{k} of region \eqn{r} is \eqn{w_{ir} z_{irk}}:
#' Level 1 region membership × within-region type probability. Local types never
#' change region memberships: summed over the local types of region \eqn{r}, the
#' memberships return \eqn{w_{ir}}. Because each \eqn{z_i} row sums to one, the raw row
#' sums of the combined matrix equal \eqn{\sum_r w_{ir}} — exactly 1 when
#' \code{weights = "memberships"}, and something else when
#' \code{weights = "cohesion"}, which is why \code{memberships} is returned
#' row-normalized and the normalized version is kept as \code{weighted}.
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
#'   a \code{"pulse_regions_tuning"} object from \code{\link{tune_regions}}, or a bare
#'   N x k numeric weight matrix.
#' @param variables Character vector of predictor columns to use. \code{NULL}
#'   (default) uses every numeric column of \code{env} that is not listed in
#'   \code{non_value_cols}.
#' @param non_value_cols Character vector of identifier columns to exclude
#'   (e.g. \code{"ID"}). Default \code{NULL}.
#' @param n_local_types Candidate component counts. Either an integer vector
#'   (default \code{2:9}) shared by all regions, or a list with one integer
#'   vector per region, in column order of the region weights, e.g.
#'   \code{list(3, 2:5, 4)}. A length-1 element fixes that region's number of
#'   local types; its covariance model is still chosen by weighted BIC.
#'   Candidates larger than a region's core allows are dropped for that region
#'   (with a warning when given as a list); a region left with no candidate is
#'   an error.
#' @param cutoff Numeric in \eqn{(0, 1]}. Weight at or above which a unit joins
#'   a region's core. Must exceed \code{weight_eps}. Default 0.1. Note that its
#'   meaning differs between the two \code{weights} settings: memberships shrink
#'   as \eqn{k} grows (a mean membership at \eqn{k = 20} is 0.05), whereas
#'   cohesion does not depend on \eqn{k}.
#' @param weights One of \code{"memberships"} (default; composition memberships
#'   \eqn{\Pi}, rows sum to 1) or \code{"cohesion"} (the fit's
#'   \code{cohesion}: mean co-association with the members of each region,
#'   rows do not sum to 1). Cohesion shrinks as a region grows, so large
#'   homogeneous regions get low weights and small cores; prefer
#'   \code{"memberships"}.
#' @param model_name Character vector of candidate \pkg{mclust} covariance
#'   models, compared by weighted BIC within each region. Default
#'   \code{c("VVV", "EEE", "VVI", "VII")}: full, shared full, diagonal,
#'   spherical. Pass \code{"VVV"} to reproduce the single-family behaviour of
#'   earlier versions. \code{"EEI"} and \code{"EII"} are valid but excluded
#'   by default: weighted EM for equal-volume models can need ~40x more
#'   iterations than \code{"VVV"}. With a single predictor the names are
#'   mapped to \code{"E"}/\code{"V"}.
#' @param scale Logical. Standardise every predictor over all units before
#'   fitting (default \code{TRUE}), so results do not depend on predictor
#'   units. With \code{FALSE}, models whose fit depends on units
#'   (\code{"VVI"}, \code{"VII"} and the other diagonal or spherical
#'   families) are dropped from \code{model_name}.
#' @param hc_max_units Integer \eqn{\ge 2}. Largest number of core units
#'   passed to \code{\link[mclust]{hc}} per region. Larger cores are thinned
#'   to this many units, evenly spaced in unit order, because \code{hc()}
#'   needs quadratic memory; the other units start EM from a flat row.
#'   Default 2000, the \pkg{mclust} default \code{subset} size.
#' @param hc_use Character scalar passed to \code{\link[mclust]{hc}} as
#'   \code{use}. Set explicitly (default \code{"VARS"}) because the
#'   \pkg{mclust} default has changed between versions.
#' @param drop_constant Logical. How to handle a predictor that is constant
#'   within one region's core or active set, which would give that region a
#'   singular covariance. \code{FALSE} (default) keeps it everywhere and fits
#'   the affected regions with the \pkg{mclust} default prior whose
#'   \code{scale} is computed from the covariance of all units instead of the
#'   region's own units. \code{TRUE} drops it from \emph{all} regions (the
#'   behaviour of earlier versions), with a warning. Dropping is also used, with
#'   a warning, when \code{prior} is \code{NULL} or not the default prior.
#'   Predictors constant across all units are always dropped, with a warning.
#' @param weight_eps Numeric. Units with weight at or below this are excluded
#'   from a region's fit. Default 1e-6.
#' @param prior Prior passed to \code{\link[mclust]{me.weighted}}. Default
#'   \code{mclust::priorControl()}; \code{NULL} disables regularisation (and
#'   forces \code{drop_constant = TRUE} behaviour, see above).
#' @param keep_models Logical. Return the selected fits (with the large
#'   \code{$z} slot stripped). Default \code{TRUE}.
#' @param verbose Logical. Report progress. Default \code{FALSE}.
#'
#' @return An object of class \code{"pulse_local_types"}: a list with
#'   \describe{
#'     \item{\code{memberships}}{N x sum(G_r) matrix, rows sum to 1.}
#'     \item{\code{weighted}}{The same matrix before row-normalisation.}
#'     \item{\code{hard_labels}}{Integer vector: column of the most probable
#'       local type within the unit's Level 1 hard region
#'       (\code{region_labels}), so \code{key$region[hard_labels] == region_labels}.}
#'     \item{\code{region_labels}}{Level 1 hard region per unit: the regions
#'       object's \code{hard_labels}, or the row argmax of a bare weight matrix.}
#'     \item{\code{certainty}}{\code{memberships} at \code{hard_labels}.}
#'     \item{\code{certainty_raw}}{\code{weighted} at \code{hard_labels}.}
#'     \item{\code{entropy}}{Normalized row entropy of \code{memberships}.}
#'     \item{\code{key}}{Data frame mapping each column to its region and local type.}
#'     \item{\code{n_local_types}}{Selected G per region.}
#'     \item{\code{BIC}}{G x region matrix over all requested G: the best 
#'       weighted BIC over covariance models for each G; \code{NA} where G was 
#'       not a candidate for that region.}
#'     \item{\code{BIC_by_model}}{G x model x region array of weighted BIC
#'       values; \code{NA} where a candidate was in feasible or failed.}
#'     \item{\code{model_names}}{Selected covariance model per region.}
#'     \item{\code{models}}{Selected fits, or \code{NULL}.}
#'     \item{\code{constant}}{Logical predictor x region matrix: \code{TRUE}
#'       where the predictor is constant within the region and the region was
#'       fitted with the globally scaled prior.}
#'     \item{\code{variables}, \code{cutoff}, \code{weights},
#'       \code{model_name}}{Settings used (\code{model_name} = candidates).}
#'     \item{\code{scaling}}{List \code{center}, \code{scale}; \code{NULL} if \code{scale = FALSE}.}
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
#' @importFrom mclust hc hclass unmap me.weighted priorControl nMclustParams 
#'   mclust.options hcE hcV hcEII hcVII hcEEE hcVVV
#' @importFrom stats var
#' @export
get_local_types <- function(env,
                            regions,
                            variables      = NULL,
                            non_value_cols = NULL,
                            n_local_types  = 2:9,
                            cutoff         = 0.1,
                            weights        = c("memberships", "cohesion"),
                            model_name     = c("VVV", "EEE", "VVI", "VII"),
                            scale          = TRUE,
                            hc_use         = "VARS",
                            hc_max_units   = 2000L,
                            drop_constant  = FALSE,
                            weight_eps     = 1e-6,
                            prior          = mclust::priorControl(),
                            keep_models    = TRUE,
                            verbose        = FALSE) {
        
        weights <- match.arg(weights)
        user_models <- !missing(model_name)
        if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
                stop("`scale` must be a single TRUE/FALSE value.", call. = FALSE)
        if (!is.numeric(hc_max_units) || length(hc_max_units) != 1L ||
            is.na(hc_max_units) || hc_max_units < 2)
                stop("`hc_max_units` must be a single number >= 2.", call. = FALSE)
        hc_max_units <- as.integer(hc_max_units)
        
        
        ## ---- 0. Input validation -------------------------------------------*
        if (!is.data.frame(env))
                stop("`env` must be a data.frame or sf object.", call. = FALSE)
        if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
                stop("`verbose` must be a single TRUE/FALSE value.", call. = FALSE)
        ## n_local_types: one shared candidate vector, or a list with one
        ## candidate vector per region (length-1 element = fixed G). The list
        ## length is checked against the region count in step 1.
        .check_G <- function(g, what) {
                if (!is.numeric(g) || !length(g) || anyNA(g) ||
                    any(g < 1) || any(g != round(g)))
                        stop(sprintf("%s must be one or more positive integers.", what),
                             call. = FALSE)
                sort(unique(as.integer(g)))
        }
        G_per_region <- is.list(n_local_types)
        if (G_per_region) {
                if (!length(n_local_types))
                        stop("`n_local_types` is an empty list.", call. = FALSE)
                n_local_types <- lapply(seq_along(n_local_types), function(r)
                        .check_G(n_local_types[[r]], sprintf("`n_local_types[[%d]]`", r)))
        } else {
                n_local_types <- .check_G(n_local_types, "`n_local_types`")
        }
        
        
        if (!is.null(variables) && !is.character(variables))
                stop("`variables` must be a character vector or NULL.", call. = FALSE)
        ## The core must be a subset of the fitted set, or `match(core, act)`
        ## below produces NA row indices.
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
        G_list <- if (G_per_region) {
                if (length(n_local_types) != n_regions)
                        stop(sprintf(paste0("`n_local_types` is a list of length %d but there ",
                                            "are %d regions; supply one element per region, ",
                                            "in column order of the region weights."),
                                     length(n_local_types), n_regions), call. = FALSE)
                n_local_types
        } else rep(list(n_local_types), n_regions)
        G_all <- sort(unique(unlist(G_list)))
        ## ---- 2. Predictors --------------------------------------------------*
        if (inherits(env, "sf")) env <- sf::st_drop_geometry(env)
        
        ## Default: all numeric columns except the identifiers.
        if (is.null(variables)) {
                cand      <- setdiff(names(env), non_value_cols)
                is_num    <- vapply(env[cand], is.numeric, logical(1))
                variables <- cand[is_num]
                if (verbose)
                        message(sprintf("Using %d numeric predictor(s): %s",
                                        length(variables), paste(variables, collapse = ", ")))
        } else {
                missing_v <- setdiff(variables, names(env))
                if (length(missing_v))
                        stop("`variables` not found in `env`: ",
                             paste(missing_v, collapse = ", "), call. = FALSE)
                variables <- setdiff(variables, non_value_cols)
        }
        if (!length(variables))
                stop("No predictor columns left after pre-processing.", call. = FALSE)
        
        Xdata <- env[, variables, drop = FALSE]
        
        ## vapply() over columns: apply() would coerce the frame to one matrix.
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
        
        empty <- which(lengths(core) == 0L)
        if (length(empty))
                stop(sprintf("Region(s) %s have an empty core: no unit reaches cutoff = %g.",
                             paste(empty, collapse = ", "), cutoff), call. = FALSE)
        
        orphan <- which(rowSums(W > weight_eps) == 0L)
        if (length(orphan))
                warning(sprintf("%d unit(s) have no region weight above `weight_eps`; ",
                                length(orphan)),
                        "their local-type rows will be all zero.", call. = FALSE)
        
        ## A predictor constant across all units carries no information.
        glob_const <- !(.col_var(Xdata, seq_len(N)) > 0)
        if (any(glob_const)) {
                warning("Dropping predictor(s) constant across all units: ",
                        paste(colnames(Xdata)[glob_const], collapse = ", "), call. = FALSE)
                Xdata <- Xdata[, !glob_const, drop = FALSE]
                if (!ncol(Xdata))
                        stop("Every predictor is constant across all units.", call. = FALSE)
        }
        
        ## Predictor x region: constant within the region's core or active set.
        const <- matrix(vapply(re, function(r)
                !(.col_var(Xdata, core[[r]]) > 0) | !(.col_var(Xdata, act[[r]]) > 0),
                logical(ncol(Xdata))),
                nrow = ncol(Xdata), dimnames = list(colnames(Xdata), paste0("Region", re)))
        can_regularise <- !is.null(prior) && identical(prior$functionName, "defaultPrior")
        if (any(const)) {
                bad <- rownames(const)[rowSums(const) > 0]
                if (drop_constant || !can_regularise) {
                        warning(sprintf(paste0(
                                "Dropping predictor(s) %s from ALL regions: constant within ",
                                "region(s) %s%s."),
                                paste(bad, collapse = ", "),
                                paste(which(colSums(const) > 0), collapse = ", "),
                                if (drop_constant) "" else
                                        " and `prior` cannot regularise them (NULL or not the default prior)"),
                                call. = FALSE)
                        Xdata <- Xdata[, !(colnames(Xdata) %in% bad), drop = FALSE]
                        if (!ncol(Xdata))
                                stop("Every predictor is constant within at least one region.",
                                     call. = FALSE)
                        const <- const[colnames(Xdata), , drop = FALSE]
                } else {
                        message(sprintf(paste0(
                                "Predictor(s) %s constant within region(s) %s: those regions ",
                                "use a prior scaled by the covariance of all units."),
                                paste(bad, collapse = ", "),
                                paste(which(colSums(const) > 0), collapse = ", ")))
                }
        }
        variables <- colnames(Xdata)
        ## One coordinate system for all regions, so meta-type distances compare like with like.
        scaling <- NULL
        if (scale) {
                scaling <- list(center = colMeans(Xdata),
                                scale  = sqrt(.col_var(Xdata, seq_len(N))))
                Xdata <- sweep(sweep(Xdata, 2L, scaling$center, "-"), 2L, scaling$scale, "/")
        }
        V_all <- stats::cov(Xdata)
        d     <- ncol(Xdata)
        ## mclust's default prior shrinks diagonal/spherical variances toward one scalar.
        unit_dep <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
        if (!scale && d > 1L && any(model_name %in% unit_dep)) {
                drop_m <- intersect(model_name, unit_dep)
                if (user_models)
                        warning(sprintf("Model(s) %s depend on predictor units and need `scale = TRUE`; dropped.",
                                        paste(drop_m, collapse = ", ")), call. = FALSE)
                model_name <- setdiff(model_name, drop_m)
                if (!length(model_name))
                        stop("No `model_name` left with `scale = FALSE`; use \"VVV\" or \"EEE\".", call. = FALSE)
        }
        models <- .resolve_model_names(model_name, d)
        
        ## ---- 4. Per-region fitting -----------------------------------------*
        ## OPT: the old code looped G on the outside and regions on the inside, so it
        ## held one full-N z matrix PER REGION simultaneously and refitted every
        ## region on all N rows for every candidate G. Here the region is the outer
        ## loop: the hc initialisation is computed once per region, the fit runs on
        ## the active rows only, and just the incumbent-best z is retained.
        BIC_all <- array(NA_real_,
                         dim = c(length(G_all), length(models), n_regions),
                         dimnames = list(G = as.character(G_all), model = models,
                                         region = paste0("Region", re)))
        best <- vector("list", n_regions)
        
        for (r in re) {
                a_r  <- act[[r]]
                c_r  <- core[[r]]
                w_r  <- W[a_r, r]
                X_r  <- Xdata[a_r, , drop = FALSE]
                ## core is a subset of act because cutoff > weight_eps (checked above).
                ## hc() is quadratic in memory and roughly cubic in time, so it
                ## runs on at most `hc_max_units` core units, evenly spaced in
                ## unit order (deterministic, no RNG). Core units outside the
                ## subset start from a flat 1/G row, like non-core units.
                h_r <- if (length(c_r) > hc_max_units)
                        c_r[unique(round(seq(1, length(c_r), length.out = hc_max_units)))]
                else c_r
                h_in_a <- match(h_r, a_r)
                n_eff  <- sum(w_r)
                prior_r <- if (any(const[, r]))
                        function(m, G) .global_scale_prior(prior, V_all, m, G)
                else function(m, G) prior
                
                ## Candidates a region can actually support.
                G_req <- G_list[[r]]
                G_r   <- G_req[G_req <= length(c_r) - 1L]
                if (!length(G_r))
                        stop(sprintf(paste0("Region %d supports no requested G (%s): its core ",
                                            "has %d units, so G <= %d. Lower `cutoff` or ",
                                            "request fewer local types."),
                                     r, paste(G_req, collapse = ", "), length(c_r),
                                     length(c_r) - 1L), call. = FALSE)
                if (G_per_region && length(G_r) < length(G_req))
                        warning(sprintf("Region %d: dropped G = %s (core has %d units).",
                                        r, paste(setdiff(G_req, G_r), collapse = ", "),
                                        length(c_r)), call. = FALSE)
                
                ## Initialisation: full-covariance hc needs more core units than
                ## predictors; below that, spherical hc is always defined. One hc
                ## per region, shared by every candidate model.
                hc_model <- if (length(h_r) > d) "VVV" else "EII"
                hcf <- tryCatch(mclust::hc(Xdata[h_r, , drop = FALSE],
                                           modelName = hc_model, use = hc_use),
                                error = function(e) e)
                if (inherits(hcf, "error") && hc_model != "EII") {
                        hc_model <- "EII"
                        hcf <- tryCatch(mclust::hc(Xdata[h_r, , drop = FALSE],
                                                   modelName = hc_model, use = hc_use),
                                        error = function(e) e)
                }
                if (inherits(hcf, "error"))
                        stop(sprintf("mclust::hc() failed for region %d: %s",
                                     r, conditionMessage(hcf)), call. = FALSE)
                cls <- mclust::hclass(hcf, G_r)          # length(h_r) x length(G_r)
                
                if (verbose)
                        message(sprintf("Region %d/%d: %d core, %d fitted, n_eff = %.1f, hc = %s, G in {%s}",
                                        r, n_regions, length(c_r), length(a_r), n_eff,
                                        hc_model, paste(G_r, collapse = ", ")))
                
                for (j in seq_along(G_r)) {
                        G  <- G_r[j]
                        z0 <- matrix(1 / G, nrow = length(a_r), ncol = G)
                        z0[h_in_a, ] <- mclust::unmap(cls[, j], groups = seq_len(G))
                        
                        for (m in models) {
                                ## Feasibility: BIC is meaningless with more parameters
                                ## than effective observations. This is the fallback
                                ## that lets small regions use simpler covariances.
                                npar <- mclust::nMclustParams(modelName = m, d = d, G = G)
                                if (npar >= n_eff) next
                                
                                fit <- tryCatch(
                                        mclust::me.weighted(data = X_r, modelName = m,
                                                            z = z0, weights = w_r,
                                                            prior = prior_r(m, G)),
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
                z0[h_in_a, ] <- mclust::unmap(cls[, best[[r]]$j], groups = seq_len(G))
                fit <- mclust::me.weighted(data = X_r, modelName = best[[r]]$model,
                                           z = z0, weights = w_r,
                                           prior = prior_r(best[[r]]$model, G))
                best[[r]]$z     <- fit$z
                fit$z           <- NULL                  # stored model without the bulk
                best[[r]]$fit   <- fit
        }
        
        ## G x region summary (best over models), so existing BIC plots still work.
        BIC <- apply(BIC_all, c(1L, 3L), function(v)
                if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
        BIC <- matrix(BIC, nrow = length(G_all),
                      dimnames = list(as.character(G_all), paste0("Region", re)))
        
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
        
        ## Crisp local type = argmax of z inside the unit's Level 1 hard region.
        ## Within a region block Zw = z * w with w constant per row, so the
        ## block argmax equals argmax(z). A global max.col(U) could pick a type
        ## from a region the unit does not belong to.
        reg_hard <- .region_hard(regions, W)
        hard     <- integer(N)
        no_block <- integer(0)
        for (r in re) {
                i <- which(reg_hard == r)
                if (!length(i)) next
                cols <- (offset[r] + 1L):offset[r + 1L]
                blk  <- Zw[i, cols, drop = FALSE]
                hard[i]  <- offset[r] + max.col(blk, ties.method = "first")
                no_block <- c(no_block, i[rowSums(blk) == 0])
        }
        no_block <- setdiff(no_block, orphan)        # orphans already warned
        if (length(no_block))
                warning(sprintf(paste0("%d unit(s) have hard region weight <= `weight_eps`; ",
                                       "their local type is the region's first type with ",
                                       "certainty 0."), length(no_block)), call. = FALSE)
        idx <- cbind(seq_len(N), hard)
        
        out <- list(
                memberships   = U,
                weighted      = Zw,
                hard_labels  = as.integer(hard),
                region_labels = as.integer(reg_hard),
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
                constant      = const,
                variables     = variables,
                scaling = scaling,
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
