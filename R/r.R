## R/ensemble_diagnostics.R  (new file)

#' Check whether cut placement follows the environment
#'
#' For a sweep of \code{size_influence}, cuts every tree as the ensemble does
#' and records where the cuts land: how heavy the cut edges are relative to
#' the other edges of the same tree, and how large the smaller piece is.
#' Answers whether the size factor leaves room for environmental contrast to
#' decide where a tree is cut, or whether the environment enters mainly
#' through the tree shape.
#'
#' @details
#' Each cut edge gets its within-tree weight percentile (mid-ranks, 0 to 1).
#' Two reference rows anchor the scale: \code{"random"} (uniform cuts, the
#' floor) and \code{size_influence = 0} (environment only, the ceiling).
#' \code{env_retained} is the mean percentile at each setting, rescaled so the
#' floor is 0 and the ceiling is 1: the share of the environment-only
#' preference for heavy edges that survives the size factor.
#'
#' All settings reuse the same exponential draws per tree (common random
#' numbers), so differences between rows are due to \code{size_influence},
#' not to sampling noise.
#'
#' Only cutting is simulated; pieces are not labelled. Cost is
#' \code{N_tree} x \code{length(size_influence)} partial sorts of \eqn{N - 1}
#' keys plus one subtree-size pass per tree: seconds on 10^4 catchments,
#' about a minute per 100 trees on 2.7 x 10^5.
#'
#' @param trees A \code{"spanning_trees"} object.
#' @param regions_per_partition Integer \eqn{R \ge 2}.
#' @param size_influence Numeric vector of settings to compare. 0 is always
#'   added (it is the ceiling reference).
#' @param weight_exponent \eqn{\gamma}. Default 1.
#' @param prior_strength \eqn{\lambda}, as in \code{\link{compute_ensemble_memberships}}.
#' @param N_tree Number of trees to use. Default \code{min(100, n_trees)}.
#' @param seed Optional integer; the caller's RNG state is restored.
#'
#' @return Object of class \code{"cut_placement"}: \code{summary} (one row per
#'   setting: \code{above_median} share of cuts on edges heavier than the tree
#'   median, \code{mean_pct}, \code{env_retained}, \code{singleton_share} share
#'   of cuts that split off one catchment, \code{median_small_piece}),
#'   \code{cuts} (one row per cut: setting, tree, percentile, small-piece
#'   size) and \code{config}.
#' @export
check_cut_placement <- function(trees, regions_per_partition,
                                size_influence  = c(0, 0.5, 1, 2, 4),
                                weight_exponent = 1, prior_strength = 0,
                                N_tree = NULL, seed = NULL) {
        stopifnot(inherits(trees, "spanning_trees"),
                  is.numeric(size_influence), all(size_influence >= 0),
                  length(regions_per_partition) == 1L, regions_per_partition >= 2,
                  length(weight_exponent) == 1L, weight_exponent >= 0,
                  length(prior_strength) == 1L, is.finite(prior_strength),
                  prior_strength >= 0)
        if (is.null(N_tree)) N_tree <- min(100L, trees$n_trees)
        stopifnot(N_tree >= 1, N_tree <= trees$n_trees)
        .local_seed(seed)

        si    <- sort(unique(c(0, size_influence)))
        lev   <- c("random", as.character(si))
        n_v   <- trees$n_vertices
        n_cut <- as.integer(regions_per_partition) - 1L
        wpe   <- trees$within_prior_edge
        out   <- vector("list", N_tree)

        for (j in seq_len(N_tree)) {
                tree     <- trees$trees[[j]]
                w        <- igraph::E(tree)$weight
                edge_ids <- igraph::E(tree)$eid
                n_e      <- length(w)
                if (n_cut >= n_e)
                        stop("regions_per_partition must be below the number of tree edges.",
                             call. = FALSE)
                n_c   <- .subtree_sizes(edge_ids, trees$edge_endpoints, n_v = n_v)
                small <- pmin(n_c, n_v - n_c)
                pct   <- (rank(w) - 0.5) / n_e
                above <- w > stats::median(w)
                le    <- log(stats::rexp(n_e))                  # shared draws
                wp    <- if (is.null(wpe)) NULL else wpe[edge_ids]

                pick <- function(key) {
                        cut <- which(key <= sort(key, partial = n_cut)[n_cut])
                        cut[seq_len(min(length(cut), n_cut))]
                }
                rows <- vector("list", length(si) + 1L)
                cut  <- pick(le)                                # uniform cuts
                rows[[1L]] <- data.frame(setting = "random", tree = j, pct = pct[cut],
                                         above = above[cut], small = small[cut])
                for (s in seq_along(si)) {
                        lp  <- .cut_log_propensity(w, n_c * (n_v - n_c) / n_v,
                                                   weight_exponent = weight_exponent,
                                                   size_influence  = si[s],
                                                   prior_strength  = prior_strength,
                                                   within_prior    = wp)$lp
                        cut <- pick(le - lp)
                        rows[[s + 1L]] <- data.frame(setting = lev[s + 1L], tree = j,
                                                     pct = pct[cut], above = above[cut],
                                                     small = small[cut])
                }
                out[[j]] <- do.call(rbind, rows)
        }
        cuts <- do.call(rbind, out)
        cuts$setting <- factor(cuts$setting, levels = lev)

        agg <- function(f) vapply(split(cuts, cuts$setting), f, numeric(1))
        mp  <- agg(function(d) mean(d$pct))
        rng <- mp[["0"]] - mp[["random"]]
        summary <- data.frame(
                setting            = lev,
                above_median       = agg(function(d) mean(d$above)),
                mean_pct           = mp,
                env_retained       = if (rng > 0) (mp - mp[["random"]]) / rng else NA_real_,
                singleton_share    = agg(function(d) mean(d$small == 1L)),
                median_small_piece = agg(function(d) stats::median(d$small)),
                row.names = NULL)
        structure(list(summary = summary, cuts = cuts,
                       config = list(regions_per_partition = as.integer(regions_per_partition),
                                     weight_exponent = weight_exponent,
                                     prior_strength = prior_strength,
                                     N_tree = as.integer(N_tree), n_vertices = n_v)),
                  class = "cut_placement")
}

#' @export
print.cut_placement <- function(x, ...) {
        cat(sprintf("<cut_placement>  N_tree = %d | R = %d | weight_exponent = %g | N = %d\n",
                    x$config$N_tree, x$config$regions_per_partition,
                    x$config$weight_exponent, x$config$n_vertices))
        s <- x$summary
        s[] <- lapply(s, function(v) if (is.numeric(v)) signif(v, 3) else v)
        print(s, row.names = FALSE)
        cat("env_retained: 0 = cuts ignore environment (random), 1 = as at size_influence = 0.\n")
        invisible(x)
}


#' Check whether the ensemble has enough trees (split-half convergence)
#'
#' Draws two disjoint sets of \code{m} ensemble columns, builds a consensus
#' from each with \code{\link{cluster_consensus}}, and compares the two hard
#' partitions by the adjusted Rand index (ARI). Repeated over a ladder of
#' sizes \code{m} up to half the ensemble, this traces how reproducible the
#' consensus is as the ensemble grows.
#'
#' @details
#' Read the curve, not one number. \code{N_tree} is large enough when the mean
#' split-half ARI stops rising: \code{converged = TRUE} when the gain over
#' the last doubling is below \code{tol}. The full ensemble is at least as
#' reproducible as either half, so a flat curve at \code{M / 2} is a
#' conservative pass for \code{M}.
#'
#' A plateau below 1 is not an ensemble-size problem. \code{search_floor} is
#' the ARI between two consensus runs on the same columns with different
#' random restarts. If split-half ARI is close to it, more trees will not
#' help; the remainder is salso search noise (raise \code{n_runs}) or real
#' ambiguity between near-equal partitions.
#'
#' Run it for the \code{final_regions} you intend to report; convergence at
#' one k does not imply it at another.
#'
#' Cost: \code{2 x B x length(sizes) + 2} consensus runs. With the defaults
#' and \code{M = 500}, 26 runs.
#'
#' @param ensemble An \code{"ensemble_memberships"} object or an integer
#'   label matrix (units x partitions).
#' @param k Integer. Number of consensus regions.
#' @param graph Optional \code{igraph}; if given, contiguity is repaired as in
#'   \code{cluster_consensus()} (consensus partition only).
#' @param sizes Integer vector of half-ensemble sizes. Default: \eqn{M/16},
#'   \eqn{M/8}, \eqn{M/4}, \eqn{M/2}, keeping sizes \eqn{\ge 10}.
#' @param B Integer. Split-half replicates per size. Default 3.
#' @param tol Numeric. ARI gain over the last doubling below which the curve
#'   counts as flat. Default 0.01.
#' @param large_n_threshold,sample_size,n_runs,n_cores Passed to
#'   \code{\link{cluster_consensus}}.
#' @param seed Optional integer; the caller's RNG state is restored.
#' @param verbose Logical.
#'
#' @return Object of class \code{"ensemble_size_check"}: \code{replicates}
#'   (one row per size and replicate), \code{summary} (per size: mean, min
#'   and max ARI, gain over the previous size), \code{search_floor},
#'   \code{converged}, \code{k}, \code{n_members}.
#' @export
check_ensemble_size <- function(ensemble, k, graph = NULL, sizes = NULL,
                                B = 3L, tol = 0.01,
                                large_n_threshold = 20000, sample_size = 5000L,
                                n_runs = 16L, n_cores = 1L,
                                seed = NULL, verbose = TRUE) {
        memb <- .ensemble_meta(ensemble)$memb_mat
        M    <- ncol(memb)
        stopifnot(length(k) == 1L, k >= 2, B >= 1, tol >= 0)
        if (M < 20L)
                stop("The ensemble needs at least 20 partitions for a split-half check.",
                     call. = FALSE)
        if (is.null(sizes)) {
                sizes <- unique(floor(M / c(16, 8, 4, 2)))
                sizes <- sizes[sizes >= 10L]
        }
        sizes <- sort(unique(as.integer(sizes)))
        if (any(sizes < 2L) || max(sizes) > M %/% 2L)
                stop(sprintf("`sizes` must lie in 2..%d (half the ensemble).", M %/% 2L),
                     call. = FALSE)
        .local_seed(seed)

        hard_of <- function(cols)
                cluster_consensus(memb[, cols, drop = FALSE], k, graph = graph,
                                  large_n_threshold = large_n_threshold,
                                  sample_size = sample_size, n_runs = n_runs,
                                  n_cores = n_cores, report_native = FALSE,
                                  verbose = FALSE)$hard_clusters

        reps <- vector("list", length(sizes) * B)
        i <- 0L
        for (m in sizes) for (b in seq_len(B)) {
                i <- i + 1L
                if (verbose) message(sprintf("Split-half: N_tree = %d, replicate %d / %d ...", m, b, B))
                idx <- sample.int(M, 2L * m)
                reps[[i]] <- data.frame(
                        N_tree = m, rep = b,
                        ARI = .adjusted_rand_index(hard_of(idx[seq_len(m)]),
                                                   hard_of(idx[m + seq_len(m)])))
        }
        reps <- do.call(rbind, reps)

        # Search-noise floor: same columns, two independent consensus runs.
        cols  <- sample.int(M, max(sizes))
        sfloor <- .adjusted_rand_index(hard_of(cols), hard_of(cols))

        s <- do.call(rbind, lapply(split(reps$ARI, reps$N_tree), function(a)
                data.frame(mean_ARI = mean(a), min_ARI = min(a), max_ARI = max(a))))
        s <- cbind(N_tree = sizes, s, row.names = NULL)
        s$gain <- c(NA_real_, diff(s$mean_ARI))
        converged <- nrow(s) >= 2L && is.finite(s$gain[nrow(s)]) && s$gain[nrow(s)] < tol

        structure(list(replicates = reps, summary = s, search_floor = sfloor,
                       converged = converged, tol = tol, k = as.integer(k),
                       n_members = M),
                  class = "ensemble_size_check")
}

#' @export
print.ensemble_size_check <- function(x, ...) {
        cat(sprintf("<ensemble_size_check>  k = %d | ensemble = %d partitions | tol = %g\n",
                    x$k, x$n_members, x$tol))
        s <- x$summary
        s[] <- lapply(s, function(v) if (is.numeric(v)) signif(v, 3) else v)
        print(s, row.names = FALSE)
        cat(sprintf("Search floor (same columns, two runs): ARI = %.3f\n", x$search_floor))
        cat(if (x$converged)
                sprintf("Flat over the last doubling: %d partitions are enough at k = %d.\n",
                        x$n_members, x$k)
            else "Still rising over the last doubling: add trees and rerun.\n")
        invisible(x)
}
