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
#' keys plus one subtree-size pass per tree: seconds on 10^4 units,
#' about a minute per 100 trees on 2.7 x 10^5.
#'
#' @param trees A \code{"pulse_trees"} object.
#' @param regions_per_partition Integer \eqn{R \ge 2}.
#' @param size_influence Numeric vector of settings to compare. 0 is always
#'   added (it is the ceiling reference).
#' @param weight_exponent \eqn{\gamma}. Default 1.
#' @param prior_strength \eqn{\lambda}, as in \code{\link{compute_ensemble}}.
#' @param N_tree Number of trees to use. Default \code{min(100, n_trees)}.
#' @param seed Optional integer; the caller's RNG state is restored.
#'
#' @return Object of class \code{"pulse_cut_placement"}: \code{summary} (one row per
#'   setting: \code{above_median} share of cuts on edges heavier than the tree
#'   median, \code{mean_pct}, \code{env_retained}, \code{singleton_share} share
#'   of cuts that split off one unit, \code{median_small_piece}),
#'   \code{cuts} (one row per cut: setting, tree, percentile, small-piece
#'   size) and \code{config}.
#' @export
check_cut_placement <- function(trees, regions_per_partition,
                                size_influence  = c(0, 0.5, 1, 2, 4),
                                weight_exponent = 1, prior_strength = 0,
                                N_tree = NULL, seed = NULL) {
        stopifnot(inherits(trees, "pulse_trees"),
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
                  class = "pulse_cut_placement")
}

#' @rdname check_cut_placement
#' @param x An object of class \code{"pulse_cut_placement"}.
#' @param ... Unused.
#' @return \code{print.pulse_cut_placement} returns \code{x} invisibly.
#' @export
print.pulse_cut_placement <- function(x, ...) {
        cat(sprintf("<pulse_cut_placement>  N_tree = %d | R = %d | weight_exponent = %g | N = %d\n",
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
#' from each with \code{\link{compute_consensus}}, and compares the two hard
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
#' @param ensemble An \code{"pulse_ensemble"} object or an integer
#'   label matrix (units x partitions).
#' @param k Integer. Number of consensus regions.
#' @param graph Optional \code{igraph}; if given, contiguity is repaired as in
#'   \code{compute_consensus()} (consensus partition only).
#' @param sizes Integer vector of half-ensemble sizes. Default: \eqn{M/16},
#'   \eqn{M/8}, \eqn{M/4}, \eqn{M/2}, keeping sizes \eqn{\ge 10}.
#' @param B Integer. Split-half replicates per size. Default 3.
#' @param tol Numeric. ARI gain over the last doubling below which the curve
#'   counts as flat. Default 0.01.
#' @param large_n_threshold,sample_size,n_runs,n_cores Passed to
#'   \code{\link{compute_consensus}}.
#' @param seed Optional integer; the caller's RNG state is restored.
#' @param verbose Logical.
#'
#' @return Object of class \code{"pulse_ensemble_size"}: \code{replicates}
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
                compute_consensus(memb[, cols, drop = FALSE], k, graph = graph,
                                  large_n_threshold = large_n_threshold,
                                  sample_size = sample_size, n_runs = n_runs,
                                  n_cores = n_cores, report_native = FALSE,
                                  verbose = FALSE)$hard_labels

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
                  class = "pulse_ensemble_size")
}

#' @rdname check_ensemble_size
#' @param x An object of class \code{"pulse_ensemble_size"}.
#' @param ... Unused.
#' @return \code{print.pulse_ensemble_size} returns \code{x} invisibly.
#' @export
print.pulse_ensemble_size <- function(x, ...) {
        cat(sprintf("<pulse_ensemble_size>  k = %d | ensemble = %d partitions | tol = %g\n",
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

#' Check how often edges are cut against their environmental distance
#'
#' For every graph edge, the share of ensemble partitions in which its two
#' endpoints fall in different regions (the cut frequency), set against the
#' environmental distance across the edge. If the ensemble follows the
#' environment, cut frequency rises with distance.
#'
#' @details
#' Edges are grouped into \code{n_bins} equal-count bins of distance. Per bin
#' the plot shows the mean cut frequency and its 90th percentile (cut
#' frequencies are skewed: most edges are never cut). The horizontal line is
#' the overall mean, the curve expected if cutting ignored distance.
#'
#' Two things push cut frequency up with distance: the tree shape (each tree
#' is the MST of perturbed weights, so it avoids heavy edges) and the cut rule
#' (\code{weight_exponent}, \code{size_influence}). To separate them, pass a
#' \code{reference} ensemble cut uniformly from the same trees:
#' \preformatted{
#' ref <- compute_ensemble(trees, N_tree, graph, R,
#'                                     weight_exponent = 0)
#' }
#' The gap between the curves is what the cut rule adds. γ = 0 ignores the 
#' weights; add size_influence = 0 for uniform random cuts.
#' Compare \code{lift} between the two, not \code{rho}: rank correlation
#' rewards the many single-unit cuts around environmental outliers.
#'
#' \code{rho} is the Spearman correlation between distance and cut frequency.
#' \code{lift} is the mean cut frequency in the top distance decile divided by
#' that in the bottom decile.
#'
#' Cost: one pass over the edges per partition; memory is two integer vectors
#' of edge length. Seconds on 10^4 units, minutes on 2.7 x 10^5 with
#' 500 trees.
#'
#' @param ensemble An \code{"pulse_ensemble"} object or an integer label
#'   matrix (units x partitions), rows in vertex order of \code{graph}.
#' @param graph The \code{igraph} the ensemble was built on.
#' @param distance Numeric, one value per edge of \code{graph}. Default
#'   \code{igraph::E(graph)$weight}. Edges with \code{NA} are left out: set
#'   bridge edges to \code{NA} to keep their penalty off the axis.
#' @param reference Optional second ensemble on the same graph (see Details).
#' @param n_bins Integer. Number of equal-count distance bins. Default 20.
#'
#' @return Object of class \code{"pulse_cut_frequency"} with \code{edges} (one row
#'   per edge: \code{distance}, \code{percentile}, \code{cut_freq},
#'   \code{ref_freq}), \code{bins} (per bin: \code{percentile},
#'   \code{distance}, \code{mean}, \code{q90}, \code{ref_mean}),
#'   \code{summary} (\code{rho}, \code{lift}, \code{mean}, \code{never_cut},
#'   and the same for the reference) and \code{M}. Plot it with \code{plot()}.
#' @export
check_cut_frequency <- function(ensemble, graph, distance = NULL,
                                reference = NULL, n_bins = 20L) {
        stopifnot(inherits(graph, "igraph"), length(n_bins) == 1L, n_bins >= 2)
        ep <- igraph::ends(graph, igraph::E(graph), names = FALSE)
        if (is.null(distance)) distance <- igraph::E(graph)$weight
        if (is.null(distance) || length(distance) != nrow(ep))
                stop("`distance` must have one value per edge of `graph` ",
                     "(default: the edge `weight` attribute).", call. = FALSE)
        
        freq <- function(ens) {
                memb <- .ensemble_meta(ens)$memb_mat
                if (nrow(memb) != igraph::vcount(graph))
                        stop(sprintf("`graph` has %d vertices but an ensemble has %d rows.",
                                     igraph::vcount(graph), nrow(memb)), call. = FALSE)
                f <- integer(nrow(ep))
                for (t in seq_len(ncol(memb))) {
                        lab <- memb[, t]
                        f <- f + (lab[ep[, 1]] != lab[ep[, 2]])
                }
                list(f = f / ncol(memb), M = ncol(memb))
        }
        main <- freq(ensemble)
        ref  <- if (is.null(reference)) NULL else freq(reference)
        
        keep <- !is.na(distance)
        d    <- distance[keep]
        pct  <- (rank(d) - 0.5) / length(d)
        cf   <- main$f[keep]
        rf   <- if (is.null(ref)) rep(NA_real_, length(d)) else ref$f[keep]
        bin  <- pmin(n_bins, floor(pct * n_bins) + 1L)
        
        by_bin <- function(v, fun) vapply(split(v, factor(bin, seq_len(n_bins))),
                                          function(z) if (length(z)) fun(z) else NA_real_,
                                          numeric(1))
        bins <- data.frame(percentile = (seq_len(n_bins) - 0.5) / n_bins,
                           distance   = by_bin(d, stats::median),
                           mean       = by_bin(cf, mean),
                           q90        = by_bin(cf, function(z) stats::quantile(z, 0.9, names = FALSE)),
                           ref_mean   = by_bin(rf, mean),
                           row.names  = NULL)
        
        stats_of <- function(f) {
                if (all(is.na(f))) return(c(rho = NA, lift = NA, mean = NA, never_cut = NA))
                lo <- mean(f[pct <= 0.1]); hi <- mean(f[pct > 0.9])
                c(rho       = suppressWarnings(stats::cor(d, f, method = "spearman")),
                  lift      = if (lo > 0) hi / lo else Inf,
                  mean      = mean(f),
                  never_cut = mean(f == 0))
        }
        structure(list(
                edges   = data.frame(edge = which(keep), distance = d, percentile = pct,
                                     cut_freq = cf, ref_freq = rf),
                bins    = bins,
                summary = list(main = stats_of(cf),
                               reference = if (is.null(ref)) NULL else stats_of(rf)),
                M = main$M, n_dropped = sum(!keep)),
                class = "pulse_cut_frequency")
}

#' @rdname check_cut_frequency
#' @param x An object of class \code{"pulse_cut_frequency"}.
#' @param ... Unused.
#' @return \code{print.pulse_cut_frequency} returns \code{x} invisibly.
#' @export
print.pulse_cut_frequency <- function(x, ...) {
        s <- x$summary$main
        cat(sprintf("<pulse_cut_frequency>  %d partitions | %d edges%s\n", x$M, nrow(x$edges),
                    if (x$n_dropped) sprintf(" (%d with NA distance left out)", x$n_dropped) else ""))
        line <- function(lab, s) cat(sprintf(
                "%-10s rho = %.3f | lift top/bottom decile = %.2f | mean = %.4f | never cut = %.1f%%\n",
                lab, s[["rho"]], s[["lift"]], s[["mean"]], 100 * s[["never_cut"]]))
        line("ensemble", s)
        if (!is.null(x$summary$reference)) line("reference", x$summary$reference)
        invisible(x)
}

#' @param x_scale \code{"percentile"} (default, rank of distance, robust to
#'   skew) or \code{"distance"} (raw values, bin medians).
#' @param n_points Integer. Edges drawn as background points. Default 3000.
#' @return \code{plot.pulse_cut_frequency} returns \code{x} invisibly.
#' @rdname check_cut_frequency
#' @export
plot.pulse_cut_frequency <- function(x, x_scale = c("percentile", "distance"),
                               n_points = 3000L, ...) {
        x_scale <- match.arg(x_scale)
        e  <- x$edges; b <- x$bins
        ex <- if (x_scale == "percentile") e$percentile else e$distance
        bx <- if (x_scale == "percentile") b$percentile else b$distance
        sub <- if (nrow(e) > n_points) sample.int(nrow(e), n_points) else seq_len(nrow(e))
        has_ref <- !all(is.na(b$ref_mean))
        col_ref <- "#0072B2"
        
        graphics::plot(ex[sub], e$cut_freq[sub], pch = 16, cex = 0.4,
                       col = grDevices::adjustcolor("grey40", 0.25),
                       ylim = c(0, max(b$q90, b$mean, b$ref_mean, na.rm = TRUE) * 1.15),
                       xlab = if (x_scale == "percentile") "Environmental distance (percentile)"
                       else "Environmental distance",
                       ylab = "Cut frequency (share of partitions)",
                       main = "Where the ensemble cuts")
        graphics::abline(h = x$summary$main[["mean"]], lty = 3)
        graphics::lines(bx, b$q90, lty = 2)
        graphics::lines(bx, b$mean, lwd = 2.5)
        if (has_ref) graphics::lines(bx, b$ref_mean, lwd = 2, col = col_ref)
        graphics::mtext(sprintf("Spearman rho %.2f | lift %.1f x%s",
                                x$summary$main[["rho"]], x$summary$main[["lift"]],
                                if (has_ref) sprintf(" (reference %.1f x)",
                                                     x$summary$reference[["lift"]]) else ""),
                        side = 3, line = 0.2, cex = 0.8)
        graphics::legend("topleft", bty = "n", cex = 0.8,
                         legend = c("mean per bin", "90th percentile per bin",
                                    if (has_ref) "reference (uniform cuts)",
                                    "overall mean (no relation)"),
                         col = c("black", "black", if (has_ref) col_ref, "black"),
                         lty = c(1, 2, if (has_ref) 1, 3),
                         lwd = c(2.5, 1, if (has_ref) 2, 1))
        invisible(x)
}
#' Check how much each new tree adds to the ensemble
#'
#' Treats every partition in the ensemble as a set of \emph{boundary edges}:
#' graph edges whose two endpoints get different labels. Because ensemble
#' regions are contiguous (pieces of a cut spanning tree), this set
#' determines the partition exactly, so nothing is lost by the change of
#' representation. Then asks two questions, tree by tree:
#' \describe{
#'   \item{Possibility space}{How many boundary edges does tree \eqn{t} use
#'     that no earlier tree used? Summed up, this is an accumulation curve,
#'     read like a species accumulation curve with edges as species and trees
#'     as samples.}
#'   \item{Difference}{How far is tree \eqn{t} from the most similar earlier
#'     tree (Jaccard distance between boundary-edge sets)?}
#' }
#'
#' @details
#' Trees are exchangeable draws, so their order is arbitrary. The smooth
#' curves are therefore the exact expectation over all orderings
#' (incidence-based rarefaction): with \eqn{f_e} the number of trees that use
#' edge \eqn{e} and \eqn{M} trees,
#' \deqn{E[S(m)] = \sum_e 1 - \binom{M - f_e}{m} / \binom{M}{m}.}
#' Points are the trees in their stored order.
#'
#' \code{coverage} is the Chao and Jost (2012) sample coverage for incidence
#' data. \code{1 - coverage} estimates the share of the boundary edges of the
#' \emph{next} tree that would be new: the direct answer to "is one more tree
#' worth it". \code{S_chao2} is the Chao2 lower bound on the number of
#' boundary edges the sampler can reach.
#'
#' Edge counts ignore distance: a boundary moved by one unit counts as
#' new edges. A slowly saturating curve with a small Jaccard distance means
#' the trees keep jittering boundaries locally; a large Jaccard distance means
#' they keep proposing different regionalisations.
#'
#' Cost: one pass over the edges per tree, plus a sparse
#' \eqn{M \times M} cross-product. Memory is about 12 bytes per boundary
#' edge per tree (e.g. 300 MB for 500 trees with 50,000 boundary edges each).
#'
#' @param ensemble An \code{"pulse_ensemble"} object or an integer label
#'   matrix (units x partitions), rows in vertex order of \code{graph}.
#' @param graph The \code{igraph} the ensemble was built on.
#'
#' @return Object of class \code{"pulse_tree_novelty"} with \code{per_tree} (one
#'   row per tree: \code{n_boundary}, \code{new_edges}, \code{share_new},
#'   \code{nn_jaccard}), \code{curve} (expected distinct edges and expected
#'   share new, by number of trees), \code{summary} (\code{S_obs},
#'   \code{S_chao2}, \code{coverage}, \code{Q1}, \code{Q2}) and \code{M}.
#'   Plot it with \code{plot()}.
#'
#' @references Chao A, Jost L (2012) Coverage-based rarefaction and
#'   extrapolation. \emph{Ecology} 93: 2533--2547.
#' @export
check_tree_novelty <- function(ensemble, graph) {
        memb <- .ensemble_meta(ensemble)$memb_mat
        stopifnot(inherits(graph, "igraph"))
        N <- nrow(memb); M <- ncol(memb)
        if (igraph::vcount(graph) != N)
                stop(sprintf("`graph` has %d vertices but the ensemble has %d rows.",
                             igraph::vcount(graph), N), call. = FALSE)
        if (M < 2L) stop("The ensemble needs at least 2 partitions.", call. = FALSE)
        
        ep <- igraph::ends(graph, igraph::E(graph), names = FALSE)
        bd <- lapply(seq_len(M), function(t) {
                lab <- memb[, t]
                which(lab[ep[, 1]] != lab[ep[, 2]])
        })
        n_t <- lengths(bd)
        e   <- unlist(bd, use.names = FALSE)
        tr  <- rep.int(seq_len(M), n_t)
        
        # Order of first use (trees in stored order)
        first     <- tr[!duplicated(e)]
        new_edges <- tabulate(first, M)
        
        # Incidence matrix over edges that are ever a boundary
        used <- sort(unique(e))
        B <- Matrix::sparseMatrix(i = match(e, used), j = tr, x = 1,
                                  dims = c(length(used), M))
        f <- Matrix::rowSums(B)
        
        # Nearest earlier tree, Jaccard distance on boundary-edge sets
        I  <- as.matrix(Matrix::crossprod(B))
        nn <- c(NA_real_, vapply(2:M, function(t) {
                s <- seq_len(t - 1L)
                u <- n_t[t] + n_t[s] - I[t, s]
                min(ifelse(u > 0, 1 - I[t, s] / u, 0))
        }, numeric(1)))
        
        # Exact expectation over tree orderings (incidence rarefaction)
        Q  <- tabulate(f, M)
        k  <- which(Q > 0)
        Es <- vapply(seq_len(M), function(m)
                sum(Q[k] * (1 - exp(lchoose(M - k, m) - lchoose(M, m)))), numeric(1))
        gain <- diff(c(0, Es))
        
        # Chao2 and sample coverage at M
        Q1 <- Q[1L]; Q2 <- if (M >= 2L) Q[2L] else 0
        S_obs  <- length(used)
        S_chao <- S_obs + (M - 1) / M *
                if (Q2 > 0) Q1^2 / (2 * Q2) else Q1 * (Q1 - 1) / 2
        cover  <- if (Q1 == 0) 1 else
                1 - Q1 / sum(f) * ((M - 1) * Q1 / ((M - 1) * Q1 + 2 * Q2))
        
        structure(list(
                per_tree = data.frame(tree = seq_len(M), n_boundary = n_t,
                                      new_edges = new_edges,
                                      share_new = ifelse(n_t > 0, new_edges / n_t, 0),
                                      nn_jaccard = nn),
                curve    = data.frame(N_tree = seq_len(M), expected_edges = Es,
                                      expected_share_new = gain / mean(n_t)),
                summary  = list(S_obs = S_obs, S_chao2 = S_chao, coverage = cover,
                                Q1 = Q1, Q2 = Q2),
                M = M),
                class = "pulse_tree_novelty")
}

#' @rdname check_tree_novelty
#' @param x An object of class \code{"pulse_tree_novelty"}.
#' @param ... Unused.
#' @return \code{print.pulse_tree_novelty} returns \code{x} invisibly.
#' @export
print.pulse_tree_novelty <- function(x, ...) {
        s <- x$summary
        cat(sprintf("<pulse_tree_novelty>  %d trees | %d boundary edges seen | Chao2 >= %.0f\n",
                    x$M, s$S_obs, s$S_chao2))
        cat(sprintf("Coverage %.3f: the next tree is expected to add %.1f%% new boundary edges.\n",
                    s$coverage, 100 * (1 - s$coverage)))
        cat(sprintf("Median Jaccard distance to nearest earlier tree (last 10%%): %.3f\n",
                    stats::median(utils::tail(x$per_tree$nn_jaccard,
                                              max(1L, x$M %/% 10L)), na.rm = TRUE)))
        invisible(x)
}

#' @param log_x Logical. Log-scale the tree axis. Default \code{FALSE}.
#' @return \code{plot.pulse_tree_novelty} returns \code{x} invisibly.
#' @rdname check_tree_novelty
#' @export
plot.pulse_tree_novelty <- function(x, log_x = FALSE, ...) {
        op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
        on.exit(graphics::par(op))
        s   <- x$summary; pt <- x$per_tree; cv <- x$curve
        lg  <- if (log_x) "x" else ""
        col_pt <- grDevices::adjustcolor("grey30", 0.5)
        col_nn <- "#D55E00"
        
        # A: possibility space
        graphics::plot(cv$N_tree, cv$expected_edges, type = "l", lwd = 2, log = lg,
                       ylim = c(0, max(s$S_chao2, s$S_obs)),
                       xlab = "Number of trees", ylab = "Distinct boundary edges",
                       main = "Possibility space explored")
        graphics::lines(pt$tree, cumsum(pt$new_edges), col = col_pt)
        graphics::abline(h = s$S_chao2, lty = 2)
        graphics::mtext(sprintf("coverage %.3f | dashed: Chao2 estimate", s$coverage),
                        side = 3, line = 0.2, cex = 0.8)
        
        # B: what each tree adds
        graphics::plot(pt$tree, pt$share_new, pch = 16, cex = 0.6, col = col_pt,
                       ylim = c(0, 1), log = lg, xlab = "Tree",
                       ylab = "Share of tree (0 to 1)",
                       main = "What each new tree adds")
        graphics::lines(cv$N_tree, pmin(cv$expected_share_new, 1), lwd = 2)
        ok <- is.finite(pt$nn_jaccard)
        graphics::points(pt$tree[ok], pt$nn_jaccard[ok], pch = 16, cex = 0.5,
                         col = grDevices::adjustcolor(col_nn, 0.3))
        sm <- stats::lowess(log(pt$tree[ok]), pt$nn_jaccard[ok], f = 0.3)
        graphics::lines(exp(sm$x), sm$y, col = col_nn, lwd = 2)
        graphics::legend("topright", bty = "n", cex = 0.8,
                         legend = c("new boundary edges (tree)", "expected",
                                    "distance to nearest earlier tree"),
                         col = c(col_pt, "black", col_nn), pch = c(16, NA, 16),
                         lty = c(NA, 1, 1), lwd = c(NA, 2, 2))
        invisible(x)
}