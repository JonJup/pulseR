# =============================================================================*
# Package utility helpers: neighbour lists, validation, messages, RNG, ARI
# =============================================================================*


#' Convert an spdep neighbour list to an edge list
#'
#' Flattens an \code{spdep} \code{nb} object into a two-column \code{data.frame}
#' of undirected edges (upper triangle only), using the neighbour list's own
#' 1-based region indices.
#'
#' @param nb An \code{spdep} \code{nb} object (e.g. from
#'   \code{\link[spdep]{poly2nb}}).
#'
#' @return A \code{data.frame} with integer columns \code{from} and \code{to}.
#'   Returns a 0-row (but correctly typed) data.frame when there are no edges.
#'
#' @details spdep encodes a region with no neighbours as the single value
#'   \code{0L}; those entries are dropped. Only pairs with \code{to > from} are
#'   kept, so each undirected edge appears once.
#'
#' @keywords internal
#' @noRd
nb_to_edgelist <- function(nb) {
        if (!inherits(nb, "nb"))
                stop("`nb` must be a neighbour list (an spdep 'nb' object).", call. = FALSE)
        
        n    <- length(nb)
        from <- rep.int(seq_len(n), lengths(nb))
        to   <- unlist(nb, use.names = FALSE)
        keep <- to != 0L & to > from        # drop "no neighbour" (0) and lower triangle
        
        data.frame(from = from[keep], to = to[keep])
}


#' Emit a progress message when verbose
#' @keywords internal
#' @noRd
.log <- function(verbose, fmt, ...) {
        if (isTRUE(verbose)) message(sprintf(fmt, ...))
}


#' Stop unless a package is installed
#' @keywords internal
#' @noRd
.need <- function(pkg) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
                stop(sprintf("Package '%s' is required for this operation. Install it with install.packages('%s').",
                             pkg, pkg),
                     call. = FALSE)
        }
}


#' Is a package installed? (A wrapper so tests can mock it.)
#' @keywords internal
#' @noRd
.has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)


#' Validate a logical(1) flag
#' @keywords internal
#' @noRd
.assert_flag <- function(x, nm) {
        if (!is.logical(x) || length(x) != 1L || is.na(x)) {
                stop(sprintf("`%s` must be a single TRUE/FALSE value.", nm), call. = FALSE)
        }
}


#' Validate a finite, non-negative numeric(1)
#' @keywords internal
#' @noRd
.assert_nonneg <- function(x, nm) {
        # is.finite() also rules out NA, NaN and Inf in one check.
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
                stop(sprintf("`%s` must be a single finite, non-negative number.", nm), call. = FALSE)
        }
}


#' Warn when working in a geographic (longitude/latitude) CRS
#'
#' Distance thresholds and centroid computations are unreliable in lon/lat;
#' this surfaces that before it silently corrupts results.
#'
#' @param crs A CRS or anything accepted by \code{\link[sf]{st_is_longlat}}.
#' @param verbose Logical; gates the (non-fatal) geographic-CRS message.
#'
#' @keywords internal
#' @noRd
.warn_if_longlat <- function(crs, verbose) {
        # Evaluate once: NA means the CRS is missing/unknown.
        longlat <- sf::st_is_longlat(crs)
        
        if (isTRUE(longlat)) {
                # Single-line string: the original spanned source lines, so the message
                # printed with literal newlines and the source indentation baked in.
                .log(verbose,
                     paste("CRS is geographic (lon/lat): distance thresholds are in degrees",
                           "and centroids are planar approximations. Consider projecting first."))
        } else if (is.na(longlat)) {
                warning("Input has no CRS; distance-based arguments are interpreted in raw ",
                        "coordinate units.", call. = FALSE)
        }
}

# `.data` is used inside ggplot2::aes() in plots.R; register it so
# R CMD check does not raise "no visible binding for global variable".
#' @importFrom utils globalVariables
utils::globalVariables(".data")

#' Build a throttled progress reporter for a loop of known length
#'
#' Returns a closure \code{f(i)} that prints an elapsed/ETA line at most once
#' per \code{report_every} fraction of progress.
#'
#' @param total Integer loop length.
#' @param stage_label Character label for the progress line.
#' @param report_every Numeric in (0, 1]; minimum progress increment between
#'   messages. Default 0.10.
#' @return A function of one argument \code{i} (the current iteration).
#' @keywords internal
#' @noRd
.make_progress <- function(total, stage_label, report_every = 0.10) {
        t0        <- proc.time()[["elapsed"]]
        last_frac <- -Inf
        function(i) {
                frac <- i / total
                if (frac < 1 && (frac - last_frac) < report_every) return(invisible(NULL))
                last_frac <<- frac
                elapsed <- proc.time()[["elapsed"]] - t0
                eta <- if (i > 0 && frac < 1) elapsed / frac * (1 - frac) else 0
                message(sprintf(
                        "  %s  %d/%d (%3.0f%%) | elapsed %s | ETA %s",
                        stage_label, i, total,
                        frac * 100, .format_duration(elapsed), .format_duration(eta)
                ))
        }
}
#' Format seconds into a human-readable duration
#'
#' @param seconds Numeric scalar. Elapsed seconds.
#' @return Character scalar such as \code{"3m 05s"}.
#' @keywords internal
#' @noRd
.format_duration <- function(seconds) {
        if (is.na(seconds) || !is.finite(seconds) || seconds < 0) return("--:--")
        seconds <- round(seconds)
        if (seconds < 60)   return(sprintf("%ds", seconds))
        if (seconds < 3600) return(sprintf("%dm %02ds", seconds %/% 60, seconds %% 60))
        sprintf("%dh %02dm %02ds", seconds %/% 3600, (seconds %% 3600) %/% 60, seconds %% 60)
}



# Partition comparison helpers ---------------------------------------------

#' Adjusted Rand Index (Hubert & Arabie 1985)
#'
#' Pure-base implementation. Permutation-invariant; ARI = 1 for identical
#' partitions (up to label permutation), expected ARI = 0 under random
#' assignment.
#'
#' @param a,b Integer or character label vectors of equal length.
#' @return Numeric scalar in approximately \code{[-0.5, 1]}; \code{NA_real_}
#'   if input length < 2.
#'
#' @details Used by the STAB indices (\code{.stability_consensus()},
#'   \code{.stability_pam_diss()}) and the meta-type crossing summary. Crisp ARI
#'   on argmax labels rather than a fuzzy-ARI generalisation because (i)
#'   stability of *structure* is what matters here, not stability of *softness*,
#'   and (ii) it has no dependencies. If a non-trivial fraction of objects have
#'   a top-two membership margin below ~0.1, consider substituting a
#'   Hüllermeier–Rifqi or Campello fuzzy Rand index (not implemented).
#'
#' @keywords internal
#' @noRd
.adjusted_rand_index <- function(a, b) {
        stopifnot(length(a) == length(b))
        if (length(a) < 2) return(NA_real_)
        tab <- table(a, b)
        n   <- sum(tab)
        if (n < 2) return(NA_real_)
        sum_nij <- sum(choose(tab, 2))
        ai <- rowSums(tab); bj <- colSums(tab)
        sum_ai <- sum(choose(ai, 2))
        sum_bj <- sum(choose(bj, 2))
        expected <- (sum_ai * sum_bj) / choose(n, 2)
        max_idx  <- 0.5 * (sum_ai + sum_bj)
        denom    <- max_idx - expected
        # Degenerate case: both partitions are singletons or a single block.
        if (denom == 0) return(1)
        (sum_nij - expected) / denom
}

#' Crisp Rousseeuw silhouette from a dissimilarity matrix
#'
#' Vectorised over the cluster axis; O(N k) after the N x N matrix is
#' materialised. Intended for a subsample of moderate size (a few thousand
#' rows), not the full graph.
#'
#' @param D Dissimilarity, either a \code{dist} object or a full N x N matrix.
#' @param labels Integer vector of cluster labels of length N.
#' @return Numeric vector of per-object silhouette widths.
#' @keywords internal
#' @noRd
.rousseeuw_silhouette <- function(D, labels) {
        if (inherits(D, "dist")) D <- as.matrix(D)
        N <- length(labels)
        clusters <- sort(unique(labels))
        K <- length(clusters)
        if (K < 2L) return(rep(0, N))
        cluster_sizes <- as.integer(table(factor(labels, levels = clusters)))
        
        # N x K matrix of summed distances from each point to each cluster.
        sums <- matrix(0, nrow = N, ncol = K)
        for (j in seq_len(K)) {
                in_j <- labels == clusters[j]
                sums[, j] <- rowSums(D[, in_j, drop = FALSE])
        }
        
        own_idx <- match(labels, clusters)
        s   <- numeric(N)
        eps <- .Machine$double.eps
        for (i in seq_len(N)) {
                o <- own_idx[i]
                if (cluster_sizes[o] <= 1L) { s[i] <- 0; next }
                means <- sums[i, ] / cluster_sizes
                # Own-cluster mean excludes self (self-distance is 0, so subtract nothing
                # from the sum but divide by size - 1).
                means[o] <- sums[i, o] / (cluster_sizes[o] - 1L)
                a_i <- means[o]
                b_i <- min(means[-o])
                s[i] <- (b_i - a_i) / max(a_i, b_i, eps)
        }
        s
}

#' Seed the RNG for the rest of the calling function
#'
#' Sets \code{seed} and registers an exit handler in the caller's frame that
#' restores the caller's RNG state, including "no state yet" (it then removes
#' the \code{.Random.seed} that \code{set.seed()} created). \code{NULL} does
#' nothing.
#' @keywords internal
#' @noRd
.local_seed <- function(seed, envir = parent.frame()) {
        if (is.null(seed)) return(invisible(NULL))
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                old     <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                restore <- function() assign(".Random.seed", old, envir = .GlobalEnv)
        } else {
                restore <- function()
                        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                                rm(".Random.seed", envir = .GlobalEnv)
        }
        do.call(on.exit, list(bquote(.(restore)()), add = TRUE), envir = envir)
        set.seed(seed)
        invisible(NULL)
}


#' Warn about columns that look like unit IDs
#'
#' A column whose values are all integers and all distinct is far more likely
#' an ID than an environmental variable. In a distance or a pseudo-F it
#' dominates every other column.
#'
#' @param X Numeric matrix.
#' @param how Text naming the argument that excludes columns.
#' @return \code{X}, invisibly.
#' @keywords internal
#' @noRd
.warn_id_like <- function(X, how) {
        if (nrow(X) < 10L) return(invisible(X))
        id_like <- vapply(seq_len(ncol(X)), function(j) {
                x <- X[, j]
                !anyNA(x) && all(x == round(x)) && !anyDuplicated(x)
        }, logical(1))
        if (any(id_like)) {
                nm <- colnames(X)
                if (is.null(nm)) nm <- paste0("column ", seq_len(ncol(X)))
                warning(sprintf(paste0("Column(s) %s look like unit IDs (integers, all distinct). ",
                                       "Exclude them via %s."),
                                paste(nm[id_like], collapse = ", "), how), call. = FALSE)
        }
        invisible(X)
}
## ARI on a mass-weighted contingency table. `wt` = soft unit count per
## element (local type), so a type holding 1% of units counts 1%, not 1/C.
## choose(x, 2) generalised to real x as x(x-1)/2.
.weighted_ari <- function(a, b, wt) {
        stopifnot(length(a) == length(b), length(a) == length(wt))
        tab <- tapply(wt, list(a, b), sum, default = 0)
        c2  <- function(x) sum(x * (x - 1) / 2)
        n <- sum(tab)
        if (n < 2) return(NA_real_)
        s_ij <- c2(tab); s_a <- c2(rowSums(tab)); s_b <- c2(colSums(tab))
        expected <- s_a * s_b / (n * (n - 1) / 2)
        denom <- 0.5 * (s_a + s_b) - expected
        if (abs(denom) < 1e-12) return(1)
        (s_ij - expected) / denom
}