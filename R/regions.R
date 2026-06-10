# PULSE Regions: Fuzzy Contiguity-Constrained Clustering

# Modular pipeline implementation.
#
# Pipeline:
#   1. sample_spanning_trees()        -> sample random spanning trees
#   2. evaluate_nrst_stability()      -> pick ensemble size n.rst (knee)
#   3. compute_ensemble_memberships() -> SKATER edge removal -> region labels
#   4. cluster_consensus()            -> CLARA/PAM consensus + fuzzy U
#   5a. tune_regions()                -> hyper-parameter search (CVIs)
#   5b. get_regions()           -> one-call convenience wrapper
#
# Conventions:
#   * Public functions are exported and return small, S3-classed result objects.
#   * Internal helpers are prefixed "." and tagged @keywords internal.
#   * Any function that calls set.seed() restores the caller's RNG state on exit
#     (CRAN-friendly: never silently clobber the user's .Random.seed).



# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 1: Sample spanning trees
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Sample Random Spanning Trees from a Spatial Graph
#'
#' Draws \code{n} random spanning trees from a weighted, (ideally) connected
#' \code{igraph} object. The trees, together with the parent graph's edge
#' weights/endpoints and any resolved prior-typology information, are bundled
#' into a \code{"spanning_trees"} object that the rest of the pipeline consumes.
#'
#' @param graph An \code{igraph} object. Edges must carry a numeric
#'   \code{weight} attribute (dissimilarity / cost; finite and non-negative).
#'   Vertices may optionally carry a prior-typology attribute (see
#'   \code{prior_typology}).
#' @param n Integer scalar \eqn{\ge 1}. Number of spanning trees to sample.
#' @param seed Integer scalar or \code{NULL}. If supplied, used to seed the
#'   sampler; the caller's RNG state is restored on exit.
#' @param verbose Logical. Emit progress/diagnostic messages? Default
#'   \code{TRUE}.
#' @param prior_typology Optional. Either (a) a single character giving the
#'   name of a vertex attribute on \code{graph} that holds prior ecoregion
#'   labels (character / factor / integer); or (b) a vector of length
#'   \code{vcount(graph)} of such labels. \code{NA} indicates no prior at
#'   that vertex. When supplied, edges connecting two vertices with the same
#'   non-NA label are flagged as "within-prior" and become eligible for
#'   weight halving in the SKATER cutting step (see \code{prior_strength}
#'   in \code{\link{compute_ensemble_memberships}} and
#'   \code{\link{get_regions}}). Default \code{NULL} (no prior).
#'
#' @return An object of class \code{"spanning_trees"} — a list with:
#'   \describe{
#'     \item{\code{trees}}{List of sampled spanning-tree objects (igraph edge
#'       sequences or graphs, depending on the installed igraph version).}
#'     \item{\code{n_trees}}{Number of trees sampled.}
#'     \item{\code{n_vertices}}{Number of vertices in the original graph.}
#'     \item{\code{edge_weights}}{Numeric vector of all parent-graph edge
#'       weights, ordered as \code{E(graph)}.}
#'     \item{\code{edge_endpoints}}{Integer matrix (\code{|E| x 2}) of edge
#'       endpoint vertex indices, ordered as \code{E(graph)}.}
#'     \item{\code{seed}}{The seed used (or \code{NULL}).}
#'     \item{\code{prior_labels}}{Resolved prior labels (character vector
#'       length N), or \code{NULL}.}
#'     \item{\code{within_prior_edge}}{Logical vector length |E| flagging
#'       within-prior edges, or \code{NULL}.}
#'   }
#'
#' @seealso \code{\link{get_regions}} for the full one-call pipeline.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- sample_gnp(50, 0.2)
#' E(g)$weight <- runif(ecount(g))
#' trees <- sample_spanning_trees(g, n = 20, seed = 1)
#' print(trees)
#' }
#'
#' @importFrom igraph sample_spanning_tree E ends vcount vertex_attr vertex_attr_names edge_attr components
#' @export
sample_spanning_trees <- function(graph, n, seed = NULL, verbose = TRUE,
                                  prior_typology = NULL) {
        
        # --- Argument validation -------------------------------------------------
        stopifnot(
                inherits(graph, "igraph"),
                is.numeric(n), length(n) == 1L, is.finite(n), n >= 1
        )
        
        # Edges must have usable weights; fail early with a clear message rather than
        # producing NA-laden cut probabilities deep inside the sampler.
        if (is.null(igraph::edge_attr(graph, "weight"))) {
                stop("`graph` edges must carry a numeric `weight` attribute.", call. = FALSE)
        }
        edge_w <- igraph::E(graph)$weight
        if (anyNA(edge_w) || any(!is.finite(edge_w)) || any(edge_w < 0)) {
                stop("Edge `weight`s must be finite and non-negative.", call. = FALSE)
        }
        
        # A disconnected graph yields a spanning *forest* (n_v - n_components edges),
        # which breaks the "cut n_reg-1 edges of a tree" arithmetic downstream.
        comp <- igraph::components(graph)
        if (comp$no > 1L) {
                warning(sprintf(
                        paste0("`graph` has %d connected components: sample_spanning_tree() will ",
                               "return a spanning forest and the SKATER cut counts assume a single ",
                               "tree. Supply a connected graph for well-defined region counts."),
                        comp$no
                ), call. = FALSE)
        }
        
        # --- RNG handling: seed locally, restore the caller's state on exit ------
        if (!is.null(seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                set.seed(seed)
        }
        
        # --- Resolve prior typology (if supplied) --------------------------------
        if (!is.null(prior_typology) & prior_typology == TRUE){
                prior_typology = "prior"
                prior_labels      <- .resolve_prior_typology(graph, prior_typology)
                edge_endpoints    <- igraph::ends(graph, igraph::E(graph), names = FALSE)
                within_prior_edge <- .compute_within_prior_edge(prior_labels, edge_endpoints)
                
                if (verbose && !is.null(prior_labels)) {
                        n_cov    <- sum(!is.na(prior_labels))
                        n_within <- sum(within_prior_edge, na.rm = TRUE)
                        n_edges  <- length(within_prior_edge)
                        message(sprintf(
                                "  Prior typology: %d / %d vertices covered (%.1f%%), %d / %d edges within-prior (%.1f%%)",
                                n_cov, length(prior_labels), 100 * n_cov / length(prior_labels),
                                n_within, n_edges, 100 * n_within / max(n_edges, 1L)
                        ))
                }
        }

        
        # --- Sample the trees ----------------------------------------------------
        if (verbose) message(sprintf("Sampling %d random spanning trees ...", n))
        tick <- if (verbose) .make_progress(n, "Spanning trees") else function(i) NULL
        
        tree_list <- vector("list", n)
        for (j in seq_len(n)) {
                tree_list[[j]] <- igraph::sample_spanning_tree(graph)
                tick(j)
        }
        
        structure(
                list(
                        trees             = tree_list,
                        n_trees           = as.integer(n),
                        n_vertices        = igraph::vcount(graph),
                        edge_weights      = edge_w,
                        edge_endpoints    = edge_endpoints,
                        seed              = seed,
                        prior_labels      = prior_labels,
                        within_prior_edge = within_prior_edge
                ),
                class = "spanning_trees"
        )
}

#' @rdname sample_spanning_trees
#' @param x An object to print (e.g. of class \code{"spanning_trees"}).
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.spanning_trees} returns \code{x} invisibly.
#' @export
print.spanning_trees <- function(x, ...) {
        cat(sprintf(
                "<spanning_trees>  %d trees | %d vertices | seed = %s\n",
                x$n_trees, x$n_vertices,
                if (is.null(x$seed)) "NULL" else as.character(x$seed)
        ))
        if (!is.null(x$prior_labels)) {
                n_cov    <- sum(!is.na(x$prior_labels))
                n_within <- sum(x$within_prior_edge, na.rm = TRUE)
                cat(sprintf(
                        "  prior typology: %d/%d vertices covered | %d/%d edges within-prior\n",
                        n_cov, x$n_vertices, n_within, length(x$within_prior_edge)
                ))
        }
        invisible(x)
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 2: Evaluate n.rst stability
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Evaluate Coassociation Matrix Stability Across Ensemble Sizes
#'
#' Computes how much the consensus (Hamming) distance structure changes as more
#' spanning trees are added, and identifies a knee point where the marginal gain
#' plateaus. Used to choose the ensemble size \code{n.rst}.
#'
#' @param graph The parent \code{igraph} object the trees were sampled from.
#'   Required because membership re-computation needs the original edge set.
#' @param trees A \code{"spanning_trees"} object from
#'   \code{\link{sample_spanning_trees}}.
#' @param intermediate_regions Integer scalar \eqn{\ge 2}. Number of 
#'   intermediate regions used as anchor (typically the median of the candidate 
#'   range).
#' @param n.rst_candidates Integer vector (each \eqn{\ge 2}). Ensemble sizes to
#'   evaluate; all must be \code{<= trees$n_trees}.
#' @param stability_metric Character. One of \code{"mean_abs_diff"} (default),
#'   \code{"frobenius"}, or \code{"correlation"}.
#' @param subsample_n Integer. Number of nodes to subsample for stability
#'   comparisons, to avoid materialising N x N distance matrices on large
#'   graphs. Default 5000.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda_E} passed to
#'   the membership computation (see \code{\link{compute_ensemble_memberships}}).
#'   Default 0 (no prior weighting).
#' @param seed Integer or \code{NULL}. Seed for the subsample selection; the
#'   caller's RNG state is restored on exit.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return An object of class \code{"nrst_stability"} — a list with:
#'   \describe{
#'     \item{\code{best_nrst}}{Selected ensemble size (knee point).}
#'     \item{\code{stability}}{Data frame with columns \code{n.rst},
#'       \code{change}, and \code{is_knee}.}
#'     \item{\code{intermediate_regions}}{The anchor value used.}
#'     \item{\code{metric}}{The stability metric used.}
#'     \item{\code{prior_strength}}{The prior strength used.}
#'   }
#'
#' @examples
#' \dontrun{
#' stab <- evaluate_nrst_stability(g, trees, intermediate_regions = 20,
#'                                 n.rst_candidates = c(5, 10, 15, 20))
#' plot(stab)
#' }
#'
#' @importFrom parallelDist parallelDist
#' @importFrom stats cor
#' @export
evaluate_nrst_stability <- function(graph,
                                    trees,
                                    intermediate_regions,
                                    n.rst_candidates,
                                    stability_metric = c("mean_abs_diff",
                                                         "frobenius",
                                                         "correlation"),
                                    subsample_n    = 5000,
                                    prior_strength = 0,
                                    seed           = NULL,
                                    verbose        = TRUE) {
        
        # --- Argument validation -------------------------------------------------
        stopifnot(inherits(trees, "spanning_trees"))
        stability_metric <- match.arg(stability_metric)
        
        n.rst_candidates <- sort(unique(as.integer(n.rst_candidates)))
        stopifnot(all(n.rst_candidates >= 2),
                  max(n.rst_candidates) <= trees$n_trees)
        stopifnot(length(intermediate_regions) == 1L, intermediate_regions >= 2,
                  intermediate_regions <= trees$n_vertices)
        stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L,
                  prior_strength >= 0)
        
        # --- RNG handling --------------------------------------------------------
        if (!is.null(seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                set.seed(seed)
        }
        
        # --- Subsample selection (fixed across all candidate sizes) --------------
        n_v <- trees$n_vertices
        use_subsample <- n_v > subsample_n
        if (use_subsample) {
                sub_idx <- sort(sample.int(n_v, subsample_n))
                if (verbose) message(sprintf(
                        "  Subsampling %d / %d nodes for stability comparisons",
                        subsample_n, n_v
                ))
        } else {
                sub_idx <- seq_len(n_v)
        }
        
        if (verbose) message("  n.rst stability analysis (intermediate_regions = ",
                             intermediate_regions, ")")
        tick <- if (verbose) .make_progress(length(n.rst_candidates),
                                            "n.rst stability") else function(i) NULL
        
        # Distance over the subsample for a given ensemble size.
        compute_sub_dist <- function(n_use) {
                memb <- .compute_ensemble_memberships_internal(
                        graph = graph, trees = trees, n.rst = n_use,
                        intermediate_regions = intermediate_regions,
                        prior_strength = prior_strength
                )
                memb_sub <- memb[sub_idx, , drop = FALSE]
                parallelDist::parallelDist(memb_sub, method = "hamming")
        }
        
        # Change between two consecutive distance structures under the chosen metric.
        mat_change <- function(prev, curr) {
                v_prev <- as.numeric(prev)
                v_curr <- as.numeric(curr)
                switch(stability_metric,
                       mean_abs_diff = mean(abs(v_curr - v_prev)),
                       frobenius = {
                               sz <- attr(curr, "Size")
                               sqrt(sum((v_curr - v_prev)^2)) / sz
                       },
                       correlation = {
                               # cor() returns NA when either vector has zero variance (e.g. a
                               # degenerate single-region partition). Treat that as "undefined".
                               r <- suppressWarnings(stats::cor(v_prev, v_curr))
                               if (is.na(r)) NA_real_ else 1 - r
                       }
                )
        }
        
        prev_mat <- NULL
        changes  <- numeric(length(n.rst_candidates))
        
        for (i in seq_along(n.rst_candidates)) {
                curr_mat   <- compute_sub_dist(n.rst_candidates[i])
                changes[i] <- if (is.null(prev_mat)) NA_real_ else mat_change(prev_mat, curr_mat)
                prev_mat   <- curr_mat
                tick(i)
        }
        
        # Knee detection requires at least two finite change values.
        valid <- which(is.finite(changes))
        best  <- if (length(valid) >= 2L) {
                .find_knee(n.rst_candidates[valid], changes[valid])
        } else {
                max(n.rst_candidates)
        }
        
        stab_df <- data.frame(
                n.rst   = n.rst_candidates,
                change  = changes,
                is_knee = n.rst_candidates == best
        )
        
        if (verbose) {
                message(sprintf("  >> Knee point: n.rst = %d (change = %.6f)",
                                best, changes[n.rst_candidates == best]))
        }
        
        structure(
                list(
                        best_nrst      = best,
                        stability      = stab_df,
                        intermediate_regions = intermediate_regions,
                        metric         = stability_metric,
                        prior_strength = prior_strength
                ),
                class = "nrst_stability"
        )
}

#' @rdname evaluate_nrst_stability
#' @param x An object of class \code{"nrst_stability"}.
#' @param ... Further arguments. For \code{plot}, passed to \code{plot()};
#'   for \code{print}, unused.
#' @return \code{print.nrst_stability} returns \code{x} invisibly.
#' @export
print.nrst_stability <- function(x, ...) {
        cat(sprintf(
                "<nrst_stability>  best n.rst = %d | metric = %s | %d candidates evaluated\n",
                x$best_nrst, x$metric, nrow(x$stability)
        ))
        invisible(x)
}

#' @rdname evaluate_nrst_stability
#' @return \code{plot.nrst_stability} is called for its side effect (draws the
#'   stability curve) and returns \code{NULL} invisibly.
#' @importFrom graphics plot points legend
#' @export
plot.nrst_stability <- function(x, ...) {
        df <- x$stability[is.finite(x$stability$change), , drop = FALSE]
        if (nrow(df) == 0L) {
                warning("No finite change values to plot.", call. = FALSE)
                return(invisible(NULL))
        }
        plot(df$n.rst, df$change,
             type = "b", pch = 19,
             xlab = "n.rst", ylab = paste0("Change (", x$metric, ")"),
             main = "n.rst stability curve", ...)
        knee_row <- df[df$is_knee, , drop = FALSE]
        if (nrow(knee_row) > 0L) {
                points(knee_row$n.rst, knee_row$change, col = "red", pch = 17, cex = 2)
                legend("topright", legend = paste("knee =", knee_row$n.rst),
                       col = "red", pch = 17, bty = "n")
        }
        invisible(NULL)
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 3: Compute ensemble memberships (SKATER edge removal)
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Compute Ensemble Region Memberships via SKATER Edge Removal
#'
#' For each spanning tree, removes \code{intermediate_regions - 1} edges (sampled
#' without replacement, with probability proportional to edge weight) to produce
#' \code{intermediate_regions} connected components, then records each node's region
#' label. Returns an integer matrix (one column per tree).
#'
#' @param trees A \code{"spanning_trees"} object.
#' @param n.rst Integer scalar. Number of spanning trees to use
#'   (must be \code{<= trees$n_trees}).
#' @param graph The parent \code{igraph} object the trees were sampled from.
#'   Needed only to re-resolve tree edges on older igraph versions; always
#'   supplied for safety.
#' @param intermediate_regions Integer scalar \eqn{\ge 2}. Number of regions per tree;
#'   must not exceed the number of vertices.
#' @param prior_strength Non-negative finite numeric scalar \eqn{\lambda_E}.
#'   Halving rate for within-prior edge weights. \eqn{\lambda_E = 0} (default)
#'   disables prior weighting and reproduces the unconstrained algorithm.
#'   \eqn{\lambda_E = 1} halves within-prior weights, \eqn{\lambda_E = 2}
#'   quarters them, and so on. Has no effect if \code{trees} was constructed
#'   without a \code{prior_typology}. \code{Inf} is not accepted here (use the
#'   \code{prior_strength = Inf} short-circuit of \code{\link{get_regions}} /
#'   \code{\link{tune_regions}} instead).
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return An object of class \code{"ensemble_memberships"} — a list with:
#'   \describe{
#'     \item{\code{memberships}}{Integer matrix (\code{N x n.rst}). Each column
#'       is a region assignment from one spanning tree.}
#'     \item{\code{n.rst}}{Ensemble size used.}
#'     \item{\code{intermediate_regions}}{Number of regions per tree.}
#'     \item{\code{n_vertices}}{Number of vertices.}
#'     \item{\code{prior_strength}}{The prior strength used.}
#'   }
#'
#' @importFrom igraph E ends make_graph components
#' @export
compute_ensemble_memberships <- function(trees,
                                         n.rst,
                                         graph,
                                         intermediate_regions,
                                         prior_strength = 0,
                                         verbose = TRUE) {
        
        # --- Argument validation -------------------------------------------------
        stopifnot(inherits(trees, "spanning_trees"))
        stopifnot(length(n.rst) == 1L, n.rst >= 2, n.rst <= trees$n_trees)
        stopifnot(length(intermediate_regions) == 1L, intermediate_regions >= 2,
                  intermediate_regions <= trees$n_vertices)
        stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L,
                  prior_strength >= 0)
        if (is.infinite(prior_strength)) {
                stop("`prior_strength = Inf` is not supported here; use the exact-prior ",
                     "short-circuit in get_regions()/tune_regions().", call. = FALSE)
        }
        
        if (verbose) {
                msg <- sprintf(
                        "Computing ensemble memberships (n.rst = %d, intermediate_regions = %d",
                        n.rst, intermediate_regions
                )
                if (prior_strength > 0 && !is.null(trees$within_prior_edge)) {
                        msg <- paste0(msg, sprintf(", lambda_E = %g", prior_strength))
                }
                message(msg, ") ...")
        }
        
        t0 <- proc.time()[["elapsed"]]
        memb <- .compute_ensemble_memberships_internal(
                graph = graph, trees = trees, n.rst = n.rst,
                intermediate_regions = intermediate_regions, prior_strength = prior_strength
        )
        if (verbose) message(sprintf("  Done in %s",
                                     .format_duration(proc.time()[["elapsed"]] - t0)))
        
        structure(
                list(
                        memberships    = memb,
                        n.rst          = as.integer(n.rst),
                        intermediate_regions = as.integer(intermediate_regions),
                        n_vertices     = trees$n_vertices,
                        prior_strength = prior_strength
                ),
                class = "ensemble_memberships"
        )
}

#' @rdname compute_ensemble_memberships
#' @param x An object of class \code{"ensemble_memberships"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.ensemble_memberships} returns \code{x} invisibly.
#' @export
print.ensemble_memberships <- function(x, ...) {
        cat(sprintf(
                "<ensemble_memberships>  %d vertices x %d trees | intermediate_regions = %d\n",
                x$n_vertices, x$n.rst, x$intermediate_regions
        ))
        invisible(x)
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 4: Fuzzy consensus clustering
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Fuzzy Consensus Clustering from Ensemble Memberships
#'
#' Clusters nodes into \code{k} groups based on an ensemble membership matrix,
#' using CLARA (for large datasets) or PAM (for small ones), then converts
#' Hamming distances to medoids into fuzzy membership values.
#'
#' @param ensemble An \code{"ensemble_memberships"} object or a raw integer
#'   matrix (N x n_trees).
#' @param k Integer scalar with \eqn{2 \le k < N}. Number of final clusters.
#' @param fuzziness Numeric \eqn{> 1}. Fuzzifier \emph{m}. Default 2.
#' @param crisp Logical. If \code{TRUE}, returns one-hot (hard) memberships.
#' @param large_n_threshold Integer. Above this N, CLARA is used instead of PAM.
#'   Default 20000.
#' @param clara_samples Integer. Number of CLARA samples. Default 50.
#' @param clara_sampsize Integer or \code{NULL}. CLARA sample size; default
#'   \code{NULL} selects \code{min(N, 200 + 2k)}.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return An object of class \code{"fuzzy_clusters"} — a list with:
#'   \describe{
#'     \item{\code{memberships}}{Numeric matrix (N x k).}
#'     \item{\code{hard_clusters}}{Integer vector of hard cluster assignments.}
#'     \item{\code{medoid_idx}}{Integer vector of medoid row indices.}
#'     \item{\code{avg_sil_width}}{Hardened average silhouette width from
#'       CLARA/PAM. NOTE: computed on the argmax labels, i.e. a crisp index;
#'       do NOT use it as the primary k-selection criterion.}
#'     \item{\code{k}}{Number of clusters.}
#'     \item{\code{fuzziness}}{The fuzziness parameter used.}
#'     \item{\code{method}}{Either \code{"clara"} or \code{"pam"}.}
#'     \item{\code{cluster_result}}{The raw \code{clara} or \code{pam} object.}
#'     \item{\code{d_to_medoids}}{N x k matrix of Hamming distances from each
#'       node to each medoid. Reused by fuzzy CVI computations.}
#'     \item{\code{medoid_dist}}{k x k matrix of inter-medoid Hamming
#'       distances. Reused by XB*.}
#'   }
#'
#' @importFrom cluster clara pam
#' @importFrom parallelDist parallelDist
#' @export
cluster_consensus <- function(ensemble,
                              k,
                              fuzziness = 2,
                              crisp = FALSE,
                              large_n_threshold = 20000,
                              clara_samples = 50,
                              clara_sampsize = NULL,
                              verbose = TRUE) {
        
        # --- Coerce input to a membership matrix ---------------------------------
        if (inherits(ensemble, "ensemble_memberships")) {
                memb_mat <- ensemble$memberships
        } else if (is.matrix(ensemble)) {
                memb_mat <- ensemble
        } else {
                stop("`ensemble` must be an 'ensemble_memberships' object or an integer matrix.",
                     call. = FALSE)
        }
        
        n <- nrow(memb_mat)
        stopifnot(length(k) == 1L, k >= 2, k < n, fuzziness > 1)
        
        use_clara <- n > large_n_threshold
        
        if (verbose) message(sprintf(
                "Fuzzy consensus clustering (k = %d, N = %d, method = %s) ...",
                k, n, if (use_clara) "CLARA" else "PAM"
        ))
        
        t0 <- proc.time()[["elapsed"]]
        
        # --- Medoid clustering ---------------------------------------------------
        if (use_clara) {
                # CLARA operates directly on the data matrix (samples + PAM on each sample),
                # avoiding a full N x N dissimilarity matrix.
                sampsize <- if (is.null(clara_sampsize)) min(n, 200 + 2 * k) else clara_sampsize
                cl_res <- .clara_hamming(memb_mat, k = k,
                                         samples  = clara_samples,
                                         sampsize = sampsize,
                                         verbose  = verbose)
                medoid_idx  <- cl_res$i.med
                clustering  <- cl_res$clustering
                sil_width   <- cl_res$silinfo$avg.width
                method_used <- "clara"
        } else {
                # PAM on the full Hamming dissimilarity matrix.
                coassoc <- parallelDist::parallelDist(memb_mat, method = "hamming")
                pam_res <- cluster::pam(coassoc, k = k)
                medoid_idx  <- pam_res$id.med
                clustering  <- pam_res$clustering
                sil_width   <- pam_res$silinfo$avg.width
                method_used <- "pam"
                cl_res      <- pam_res
        }
        
        # --- Distances to medoids (always computed) ------------------------------
        # Needed both for fuzzy-membership derivation and for downstream CVIs (XB*,
        # SIL_F). Computing it unconditionally keeps the crisp = TRUE path and CVI
        # consumers well-defined.
        dists_to_med <- .hamming_to_refs(memb_mat, medoid_idx)
        
        # --- k x k inter-medoid distance matrix (separation term for XB*) --------
        if (length(medoid_idx) >= 2L) {
                medoid_dist <- as.matrix(
                        parallelDist::parallelDist(
                                memb_mat[medoid_idx, , drop = FALSE],
                                method = "hamming"
                        )
                )
        } else {
                medoid_dist <- matrix(0, 1L, 1L)
        }
        
        # --- Build memberships ---------------------------------------------------
        if (crisp) {
                # One-hot assignment, vectorised (avoid an N-length R loop).
                memberships <- matrix(0, nrow = n, ncol = k)
                memberships[cbind(seq_len(n), clustering)] <- 1
        } else {
                # Fuzzy c-means-style soft memberships from distances to medoids.
                m <- fuzziness
                memberships <- t(apply(dists_to_med, 1, function(d) {
                        # Exact zero distance => the node *is* a medoid: assign it fully.
                        if (any(d == 0)) {
                                u <- rep(0, length(d))
                                u[which.min(d)] <- 1
                                return(u)
                        }
                        d_pow <- d ^ (-2 / (m - 1))
                        d_pow / sum(d_pow)
                }))
        }
        
        if (verbose) message(sprintf("  Done in %s",
                                     .format_duration(proc.time()[["elapsed"]] - t0)))
        
        structure(
                list(
                        memberships    = memberships,
                        hard_clusters  = as.integer(clustering),
                        medoid_idx     = medoid_idx,
                        avg_sil_width  = sil_width,
                        k              = k,
                        fuzziness      = fuzziness,
                        method         = method_used,
                        cluster_result = cl_res,
                        d_to_medoids   = dists_to_med,
                        medoid_dist    = medoid_dist
                ),
                class = "fuzzy_clusters"
        )
}

#' @rdname cluster_consensus
#' @param x An object of class \code{"fuzzy_clusters"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.fuzzy_clusters} returns \code{x} invisibly.
#' @export
print.fuzzy_clusters <- function(x, ...) {
        cat(sprintf(
                "<fuzzy_clusters>  %d nodes x %d clusters | method = %s | avg silhouette = %.4f\n",
                nrow(x$memberships), x$k, x$method, x$avg_sil_width
        ))
        invisible(x)
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 5a: Tuning orchestrator
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Tune SKATER-CON Hyperparameters
#'
#' Searches over \code{intermediate_regions} and/or \code{final_regions} to find the
#' combination that optimises a clustering-quality metric, operating on a
#' pre-computed \code{"spanning_trees"} object.
#'
#' @param graph The parent \code{igraph} object the trees were sampled from.
#' @param trees A \code{"spanning_trees"} object from
#'   \code{\link{sample_spanning_trees}}.
#' @param n.rst Integer scalar ensemble size to use
#'   (\eqn{2 \le} \code{n.rst} \eqn{\le} \code{trees$n_trees}).
#' @param intermediate_regions Integer or integer vector (each \eqn{\ge 2}). Candidate
#'   number(s) of SKATER regions per tree. Default 50.
#' @param final_regions Integer or integer vector (each \eqn{\ge 2}). Candidate
#'   number(s) of final consensus clusters \emph{k}.
#' @param fuzziness Numeric \eqn{> 1}. Fuzzifier. Default 2.
#' @param crisp Logical. Produce hard memberships? Default \code{FALSE}.
#' @param strategy Character. Search strategy: \code{"grid"} (default; full
#'   factorial), \code{"sequential"} (tune \code{intermediate_regions} then
#'   \code{final_regions}), or \code{"iterative"} (alternate until convergence
#'   or \code{max_iter}).
#' @param metric Character. One of \code{"borda"} (default; rank-aggregate
#'   of \code{XB_star}, \code{SIL_F}, and \code{stability}),
#'   \code{"XB_star"} (Hamming-medoid Xie–Beni; minimise),
#'   \code{"SIL_F"} (Campello–Hruschka fuzzy silhouette; maximise),
#'   \code{"stability"} (bootstrap ARI; maximise),
#'   \code{"MPC"} / \code{"PE"} / \code{"PEN"} (membership-only diagnostics;
#'   not recommended as selectors), or the legacy
#'   \code{"silhouette"} / \code{"partition_coefficient"} /
#'   \code{"modified_partition_coefficient"} (kept for backward compatibility,
#'   but bias \emph{k} toward crisp partitions).
#' @param stability_B Integer. Bootstrap replicates for the \code{stability}
#'   sub-metric. Set to 0 to skip. Default 25.
#' @param silf_subsample_n Integer. Stratified subsample size for SIL_F when
#'   \emph{N} exceeds it. Default 5000.
#' @param silf_alpha Numeric. Campello–Hruschka \eqn{\alpha}. Default 1.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda_E}. Passed
#'   through to the membership computation. \code{Inf} short-circuits to an
#'   exact reproduction of \code{prior_typology} (which must then be supplied).
#'   Default 0.
#' @param prior_typology Optional. Vertex attribute name (single character) or
#'   length-\code{vcount(graph)} vector of prior labels. Used only by the
#'   \code{prior_strength = Inf} short-circuit; otherwise prior information is
#'   already carried by \code{trees}. Default \code{NULL}.
#' @param max_iter Integer. Maximum iterations for the \code{"iterative"}
#'   strategy. Default 10.
#' @param large_n_threshold Integer. CLARA/PAM cutoff, mirrors
#'   \code{\link{cluster_consensus}}. Default 20000.
#' @param seed Optional integer seed (controls stability subsampling); the
#'   caller's RNG state is restored on exit.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return An object of class \code{"regions_tuning"} — a list with
#'   \code{best_result} (a \code{"fuzzy_clusters"} object),
#'   \code{best_intermediate_regions}, \code{best_final_regions}, \code{best_score},
#'   \code{tuning_log} (a data frame exposing all CVI columns: PC, PE, MPC,
#'   PEN, XB_star, SIL_F, STAB, SILH_HARD, plus \code{score}), \code{metric},
#'   \code{strategy}, and \code{prior_strength}. When
#'   \code{prior_strength = Inf}, the exact-prior list from
#'   \code{.exact_prior_result} is returned instead.
#'
#' @export
tune_regions <- function(graph,
                         trees,
                         n.rst,
                         intermediate_regions = 50,
                         final_regions,
                         fuzziness  = 2,
                         crisp      = FALSE,
                         strategy   = c("grid", "sequential", "iterative"),
                         metric     = c("borda",
                                        "XB_star",
                                        "SIL_F",
                                        "stability",
                                        "MPC",
                                        "PE",
                                        "PEN",
                                        "silhouette",
                                        "partition_coefficient",
                                        "modified_partition_coefficient"),
                         stability_B      = 25L,
                         silf_subsample_n = 5000L,
                         silf_alpha       = 1,
                         prior_strength   = 0,
                         prior_typology   = NULL,
                         max_iter          = 10,
                         large_n_threshold = 20000,
                         seed              = NULL,
                         verbose           = TRUE) {
        
        stopifnot(inherits(trees, "spanning_trees"))
        strategy <- match.arg(strategy)
        metric   <- match.arg(metric)
        stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L,
                  prior_strength >= 0)
        
        # --- Exact prior reproduction short-circuit ------------------------------
        if (is.infinite(prior_strength)) {
                if (is.null(prior_typology)) {
                        stop("`prior_strength = Inf` requires `prior_typology`.", call. = FALSE)
                }
                prior_labels <- .resolve_prior_typology(graph, prior_typology)
                return(.exact_prior_result(prior_labels, verbose = verbose))
        }
        
        # Notify the user when they have opted into a known-biased legacy metric.
        legacy_biased <- c("silhouette", "partition_coefficient",
                           "modified_partition_coefficient")
        if (!crisp && metric %in% legacy_biased) {
                message("Note: '", metric, "' biases k-selection toward crisp ",
                        "partitions; 'borda' is recommended.")
        }
        
        # Crisp memberships make fuzzy-aware indices degenerate (U is one-hot, so
        # PC = 1, PE = 0, MPC = 1; SIL_F reduces to the crisp silhouette). Fall back
        # to the hardened silhouette in that case.
        crisp_degenerate <- c("partition_coefficient",
                              "modified_partition_coefficient",
                              "MPC", "PE", "PEN", "borda")
        if (crisp && metric %in% crisp_degenerate) {
                warning("Fuzzy / Borda metrics are degenerate with crisp = TRUE; ",
                        "falling back to 'silhouette'.", call. = FALSE)
                metric <- "silhouette"
        }
        
        intermediate_regions <- sort(unique(as.integer(intermediate_regions)))
        final_regions  <- sort(unique(as.integer(final_regions)))
        stopifnot(all(intermediate_regions >= 2), all(final_regions >= 2),
                  all(intermediate_regions <= trees$n_vertices))
        stopifnot(length(n.rst) == 1L, n.rst >= 2, n.rst <= trees$n_trees)
        
        # --- RNG handling --------------------------------------------------------
        if (!is.null(seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                set.seed(seed)
        }
        
        # --- Cache ensemble memberships keyed by intermediate_regions ------------------
        memb_cache <- new.env(parent = emptyenv())
        get_memb <- function(sr) {
                key <- as.character(sr)
                if (!exists(key, envir = memb_cache, inherits = FALSE)) {
                        ens <- compute_ensemble_memberships(
                                trees          = trees,
                                n.rst          = n.rst,
                                graph          = graph,
                                intermediate_regions = sr,
                                prior_strength = prior_strength,
                                verbose        = verbose
                        )
                        assign(key, ens$memberships, envir = memb_cache)
                }
                get(key, envir = memb_cache, inherits = FALSE)
        }
        
        # Canonical CVI column set (kept consistent so rbind() works across batches).
        cvi_cols <- c("PC", "PE", "MPC", "PEN",
                      "XB_star", "SIL_F", "STAB", "SILH_HARD")
        
        # --- Score one (sr, fr) pair -> one row of CVI columns -------------------
        score_pair <- function(sr, fr) {
                memb_mat <- get_memb(sr)
                res <- cluster_consensus(
                        memb_mat, k = fr, fuzziness = fuzziness, crisp = crisp,
                        large_n_threshold = large_n_threshold, verbose = FALSE
                )
                cvis <- .compute_fuzzy_cvis(
                        U                        = res$memberships,
                        d_to_medoids             = res$d_to_medoids,
                        medoid_dist              = res$medoid_dist,
                        memb_mat                 = memb_mat,
                        hard_clusters            = res$hard_clusters,
                        fuzziness                = res$fuzziness,
                        silhouette_hard          = res$avg_sil_width,
                        stability_B              = stability_B,
                        stability_subsample_frac = 0.8,
                        stability_seed           = seed,
                        silf_subsample_n         = silf_subsample_n,
                        silf_alpha               = silf_alpha,
                        large_n_threshold        = large_n_threshold,
                        verbose                  = verbose
                )
                if (verbose) {
                        message(sprintf(
                                "    sr=%d, fr=%d  | XB*=%.4f  SIL_F=%.4f  STAB=%.4f  MPC=%.4f",
                                sr, fr, cvis$XB_star, cvis$SIL_F, cvis$STAB, cvis$MPC
                        ))
                }
                data.frame(intermediate_regions = sr, final_regions = fr, cvis)
        }
        
        # --- tuning_log accumulator (CVI columns + derived score) ----------------
        tuning_log <- NULL
        
        rebuild_score_column <- function() {
                # Recompute the score column from the current sub-metric columns. Required
                # after every batch because Borda ranks are computed *globally* over all
                # rows evaluated so far.
                if (is.null(tuning_log) || nrow(tuning_log) == 0L) return(invisible(NULL))
                tuning_log$score <<- .derive_score_column(tuning_log, metric)
                invisible(NULL)
        }
        
        evaluate_pairs <- function(pairs_df, tick = NULL) {
                new_rows <- vector("list", nrow(pairs_df))
                for (i in seq_len(nrow(pairs_df))) {
                        sr <- pairs_df$intermediate_regions[i]
                        fr <- pairs_df$final_regions[i]
                        # Skip pairs already present (cache hit on the log itself).
                        if (!is.null(tuning_log)) {
                                hit <- which(tuning_log$intermediate_regions == sr &
                                                     tuning_log$final_regions  == fr)
                                if (length(hit) > 0L) {
                                        if (!is.null(tick)) tick(i)
                                        next
                                }
                        }
                        new_rows[[i]] <- score_pair(sr, fr)
                        if (!is.null(tick)) tick(i)
                }
                new_rows <- new_rows[!vapply(new_rows, is.null, logical(1))]
                if (length(new_rows) > 0L) {
                        new_df    <- do.call(rbind, new_rows)
                        col_order <- c("intermediate_regions", "final_regions", cvi_cols)
                        new_df    <- new_df[, col_order, drop = FALSE]
                        if (is.null(tuning_log)) {
                                tuning_log <<- new_df
                        } else {
                                # Drop the (stale) score column before rbind; it is rebuilt afterwards.
                                base       <- tuning_log[, col_order, drop = FALSE]
                                tuning_log <<- rbind(base, new_df)
                        }
                }
                rebuild_score_column()
                invisible(NULL)
        }
        
        # Best value of `vary_col` among rows where `fixed_col == fixed_val`.
        best_along <- function(fixed_col, fixed_val, vary_col) {
                idx <- which(tuning_log[[fixed_col]] == fixed_val)
                if (length(idx) == 0L) return(NA_integer_)
                sc <- tuning_log$score[idx]
                if (all(is.na(sc))) return(NA_integer_)
                tuning_log[[vary_col]][idx][which.max(sc)]
        }
        
        # --- Run the chosen strategy ---------------------------------------------
        if (verbose) message(sprintf("Tuning (%s strategy, metric = %s) ...",
                                     strategy, metric))
        
        if (strategy == "grid") {
                pairs <- expand.grid(intermediate_regions = intermediate_regions,
                                     final_regions  = final_regions,
                                     KEEP.OUT.ATTRS = FALSE)
                if (verbose) message(sprintf("  Grid: %d combinations", nrow(pairs)))
                tick <- if (verbose) .make_progress(nrow(pairs), "Grid search") else NULL
                evaluate_pairs(pairs, tick)
                
        } else if (strategy == "sequential") {
                fr_anchor <- final_regions[ceiling(length(final_regions) / 2)]
                if (length(intermediate_regions) > 1L) {
                        if (verbose) message("  Step 1: tuning intermediate_regions (fr = ", fr_anchor, ")")
                        tick <- if (verbose) .make_progress(length(intermediate_regions), "intermediate_regions") else NULL
                        evaluate_pairs(
                                data.frame(intermediate_regions = intermediate_regions, final_regions = fr_anchor),
                                tick
                        )
                        best_sr <- best_along("final_regions", fr_anchor, "intermediate_regions")
                        if (verbose) message("  >> Best intermediate_regions = ", best_sr)
                } else {
                        best_sr <- intermediate_regions
                }
                if (length(final_regions) > 1L) {
                        if (verbose) message("  Step 2: tuning final_regions (sr = ", best_sr, ")")
                        tick <- if (verbose) .make_progress(length(final_regions), "final_regions") else NULL
                        evaluate_pairs(
                                data.frame(intermediate_regions = best_sr, final_regions = final_regions),
                                tick
                        )
                }
                
        } else if (strategy == "iterative") {
                best_sr    <- intermediate_regions[ceiling(length(intermediate_regions) / 2)]
                best_fr    <- final_regions[ceiling(length(final_regions) / 2)]
                best_score <- -Inf
                iter       <- 0L
                repeat {
                        iter <- iter + 1L
                        if (verbose) message(sprintf("\n  --- Iteration %d/%d ---", iter, max_iter))
                        prev_score <- best_score
                        if (length(intermediate_regions) > 1L) {
                                if (verbose) message("  Tuning intermediate_regions (fr = ", best_fr, ")")
                                tick <- if (verbose) .make_progress(length(intermediate_regions),
                                                                    sprintf("Iter %d: sr", iter)) else NULL
                                evaluate_pairs(
                                        data.frame(intermediate_regions = intermediate_regions, final_regions = best_fr),
                                        tick
                                )
                                best_sr <- best_along("final_regions", best_fr, "intermediate_regions")
                        }
                        if (length(final_regions) > 1L) {
                                if (verbose) message("  Tuning final_regions (sr = ", best_sr, ")")
                                tick <- if (verbose) .make_progress(length(final_regions),
                                                                    sprintf("Iter %d: fr", iter)) else NULL
                                evaluate_pairs(
                                        data.frame(intermediate_regions = best_sr, final_regions = final_regions),
                                        tick
                                )
                                best_fr <- best_along("intermediate_regions", best_sr, "final_regions")
                        }
                        best_score <- suppressWarnings(max(tuning_log$score, na.rm = TRUE))
                        if (!is.finite(best_score) || best_score <= prev_score || iter >= max_iter) {
                                if (verbose) {
                                        if (is.finite(best_score) && best_score <= prev_score) message("  Converged.")
                                        else message("  Reached max_iter.")
                                }
                                break
                        }
                }
        }
        
        # --- Deduplicate and pick the best pair ----------------------------------
        tuning_log <- tuning_log[!duplicated(
                tuning_log[, c("intermediate_regions", "final_regions")]
        ), , drop = FALSE]
        rownames(tuning_log) <- NULL
        # One last rebuild in case dedup changed the row set Borda ranks against.
        tuning_log$score <- .derive_score_column(tuning_log, metric)
        
        if (all(is.na(tuning_log$score))) {
                stop("All hyper-parameter combinations produced an undefined score; ",
                     "check the metric and inputs.", call. = FALSE)
        }
        
        best_idx <- which.max(tuning_log$score)
        best_sr  <- tuning_log$intermediate_regions[best_idx]
        best_fr  <- tuning_log$final_regions[best_idx]
        best_sc  <- tuning_log$score[best_idx]
        
        if (verbose) {
                message(sprintf("\n>> Best: intermediate_regions = %d, final_regions = %d (%s score = %.4f)",
                                best_sr, best_fr, metric, best_sc))
                message("  Computing final memberships ...")
        }
        
        memb_mat <- get_memb(best_sr)
        best_res <- cluster_consensus(
                memb_mat, k = best_fr, fuzziness = fuzziness, crisp = crisp,
                large_n_threshold = large_n_threshold, verbose = verbose
        )
        
        structure(
                list(
                        best_result         = best_res,
                        best_intermediate_regions = best_sr,
                        best_final_regions  = best_fr,
                        best_score          = best_sc,
                        tuning_log          = tuning_log,
                        metric              = metric,
                        strategy            = strategy,
                        prior_strength      = prior_strength
                ),
                class = "regions_tuning"
        )
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 5b: Convenience wrapper
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Delineate Fuzzy Contiguity-Constrained Regions using SKATER (one-call wrapper)
#'
#' All-in-one convenience wrapper that samples spanning trees, optionally
#' chooses the ensemble size by stability analysis, optionally tunes
#' \code{intermediate_regions} / \code{final_regions}, and returns the final fuzzy
#' regionalisation. See \code{\link{tune_regions}} for the full menu of
#' \code{tuning_metric} options. Default is \code{"borda"}.
#'
#' @inheritParams sample_spanning_trees
#' @param intermediate_regions,final_regions,n.rst Integer or integer vectors of
#'   candidate values. If a vector of length > 1 is supplied for
#'   \code{intermediate_regions} or \code{final_regions}, tuning is triggered; if
#'   \code{n.rst} has length > 1, the ensemble size is chosen by stability.
#' @param fuzziness Numeric \eqn{> 1}. Default 2.
#' @param crisp Logical. Produce hard memberships? Default \code{FALSE}.
#' @param tuning_strategy Character. One of \code{"grid"} (default),
#'   \code{"sequential"}, \code{"iterative"}.
#' @param tuning_metric Character. See \code{\link{tune_regions}}. Default
#'   \code{"borda"}.
#' @param stability_B,silf_subsample_n,silf_alpha Passed through to
#'   \code{\link{tune_regions}}.
#' @param n.rst_stability_metric Character. One of \code{"mean_abs_diff"}
#'   (default), \code{"frobenius"}, \code{"correlation"}.
#' @param prior_typology Optional. Vertex attribute name (single character)
#'   or vector of length \code{vcount(graph)} giving a prior ecoregion label
#'   at each catchment. \code{NA} = no prior at that vertex. The user is
#'   expected to attach this to the polygons / graph themselves before
#'   calling. Default \code{NULL}.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda_E}. The
#'   parameter has a direct interpretation as the number of halvings applied
#'   to within-region edge weights: \eqn{\lambda_E = 0} (default) disables
#'   prior weighting, \eqn{\lambda_E = 1} halves the weight,
#'   \eqn{\lambda_E = 2} quarters it, and so on. \code{Inf} returns the prior
#'   typology exactly (requires \code{prior_typology}).
#' @param max_iter Integer. Maximum iterations for the \code{"iterative"}
#'   tuning strategy. Default 10.
#' @param verbose Logical. Default \code{TRUE}.
#' @param seed Optional integer seed; the caller's RNG state is restored on
#'   exit.
#'
#' @return A plain list with:
#'   \describe{
#'     \item{\code{memberships}}{Final N x k membership matrix.}
#'     \item{\code{hard_clusters}}{Integer vector of hard cluster labels.}
#'     \item{\code{best_intermediate_regions}, \code{best_final_regions},
#'       \code{best_n.rst}}{The selected hyper-parameters.}
#'     \item{\code{best_score}}{For a tuned run, the chosen \code{tuning_metric}
#'       score; for a single-configuration run, the hardened average silhouette
#'       width (these are NOT directly comparable).}
#'     \item{\code{tuning_log}}{The tuning data frame, or \code{NULL} if no
#'       tuning was performed.}
#'     \item{\code{n.rst_stability}}{The stability data frame, or \code{NULL}.}
#'     \item{\code{trees}}{The \code{"spanning_trees"} object used.}
#'     \item{\code{cluster_result}}{The raw CLARA/PAM object.}
#'     \item{\code{prior_strength}}{The prior strength used.}
#'   }
#'   When \code{prior_strength = Inf}, the exact-prior list (see
#'   \code{.exact_prior_result}) is returned instead.
#'
#' @seealso \code{\link{sample_spanning_trees}},
#'   \code{\link{evaluate_nrst_stability}},
#'   \code{\link{compute_ensemble_memberships}},
#'   \code{\link{cluster_consensus}}, \code{\link{tune_regions}}.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- sample_gnp(200, 0.05)
#' E(g)$weight <- runif(ecount(g))
#' res <- get_regions(g, intermediate_regions = c(10, 20), final_regions = c(3, 4),
#'                          n.rst = 30, seed = 1)
#' table(res$hard_clusters)
#' }
#'
#' @importFrom igraph vcount
#' @export
get_regions <- function(graph,
                              intermediate_regions = 50,
                              final_regions,
                              n.rst,
                              fuzziness  = 2,
                              crisp      = FALSE,
                              tuning_strategy = c("grid", "sequential", "iterative"),
                              tuning_metric   = c("borda", "XB_star", "SIL_F", "stability",
                                                  "MPC", "PE", "PEN",
                                                  "silhouette",
                                                  "partition_coefficient",
                                                  "modified_partition_coefficient"),
                              stability_B      = 25L,
                              silf_subsample_n = 5000L,
                              silf_alpha       = 1,
                              n.rst_stability_metric = c("mean_abs_diff",
                                                         "frobenius",
                                                         "correlation"),
                              prior_typology = NULL,
                              prior_strength = 0,
                              max_iter = 10,
                              verbose  = TRUE,
                              seed     = NULL) {
        
        tuning_strategy        <- match.arg(tuning_strategy)
        tuning_metric          <- match.arg(tuning_metric)
        n.rst_stability_metric <- match.arg(n.rst_stability_metric)
        
        stopifnot(inherits(graph, "igraph"))
        stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L,
                  prior_strength >= 0)
        
        # --- Exact prior reproduction short-circuit ------------------------------
        if (is.infinite(prior_strength)) {
                if (is.null(prior_typology)) {
                        stop("`prior_strength = Inf` requires `prior_typology`.", call. = FALSE)
                }
                prior_labels <- .resolve_prior_typology(graph, prior_typology)
                return(.exact_prior_result(prior_labels, verbose = verbose))
        }
        
        n.rst_values   <- sort(unique(as.integer(n.rst)))
        intermediate_regions <- sort(unique(as.integer(intermediate_regions)))
        final_regions  <- sort(unique(as.integer(final_regions)))
        
        # Friendly upper-bound check before doing any sampling work.
        if (any(intermediate_regions > igraph::vcount(graph))) {
                stop("`intermediate_regions` cannot exceed the number of vertices (",
                     igraph::vcount(graph), ").", call. = FALSE)
        }
        
        # --- Sample the tree ensemble (largest requested n.rst) ------------------
        trees <- sample_spanning_trees(graph, n = max(n.rst_values),
                                       seed = seed, verbose = verbose,
                                       prior_typology = prior_typology)
        
        # lambda_E > 0 without a prior is a silent no-op; warn the user.
        if (prior_strength > 0 && is.null(trees$prior_labels)) {
                warning("`prior_strength` > 0 but no `prior_typology` supplied; ",
                        "prior weighting has no effect.", call. = FALSE)
        }
        
        # --- Choose ensemble size n.rst via stability (if multiple candidates) ---
        n.rst_stability <- NULL
        if (length(n.rst_values) > 1L) {
                sr_anchor <- intermediate_regions[ceiling(length(intermediate_regions) / 2)]
                n.rst_stability <- evaluate_nrst_stability(
                        graph, trees, intermediate_regions = sr_anchor,
                        n.rst_candidates = n.rst_values,
                        stability_metric = n.rst_stability_metric,
                        prior_strength   = prior_strength,
                        seed             = seed,
                        verbose          = verbose
                )
                best_nrst <- n.rst_stability$best_nrst
        } else {
                best_nrst <- n.rst_values
        }
        
        needs_tuning <- length(intermediate_regions) > 1L || length(final_regions) > 1L
        
        if (needs_tuning) {
                tuning_result <- tune_regions(
                        graph          = graph,
                        trees          = trees,
                        n.rst          = best_nrst,
                        intermediate_regions = intermediate_regions,
                        final_regions  = final_regions,
                        fuzziness      = fuzziness, crisp = crisp,
                        strategy       = tuning_strategy,
                        metric         = tuning_metric,
                        stability_B      = stability_B,
                        silf_subsample_n = silf_subsample_n,
                        silf_alpha       = silf_alpha,
                        prior_strength   = prior_strength,
                        max_iter         = max_iter,
                        verbose          = verbose,
                        seed             = seed
                )
                list(
                        memberships         = tuning_result$best_result$memberships,
                        hard_clusters       = tuning_result$best_result$hard_clusters,
                        best_intermediate_regions = tuning_result$best_intermediate_regions,
                        best_final_regions  = tuning_result$best_final_regions,
                        best_n.rst          = best_nrst,
                        best_score          = tuning_result$best_score,
                        tuning_log          = tuning_result$tuning_log,
                        n.rst_stability     = if (!is.null(n.rst_stability)) n.rst_stability$stability else NULL,
                        trees               = trees,
                        cluster_result      = tuning_result$best_result$cluster_result,
                        prior_strength      = prior_strength
                )
        } else {
                ensemble <- compute_ensemble_memberships(
                        trees = trees, n.rst = best_nrst, graph = graph,
                        intermediate_regions = intermediate_regions,
                        prior_strength = prior_strength,
                        verbose = verbose
                )
                result <- cluster_consensus(
                        ensemble, k = final_regions, fuzziness = fuzziness,
                        crisp = crisp, large_n_threshold = 20000, verbose = verbose
                )
                list(
                        memberships         = result$memberships,
                        hard_clusters       = result$hard_clusters,
                        best_intermediate_regions = intermediate_regions,
                        best_final_regions  = final_regions,
                        best_n.rst          = best_nrst,
                        # NOTE: single-config mode has no tuning metric; report the hardened
                        # silhouette width so the slot is populated, but it is NOT comparable to
                        # a tuned `best_score`.
                        best_score          = result$avg_sil_width,
                        tuning_log          = NULL,
                        n.rst_stability     = if (!is.null(n.rst_stability)) n.rst_stability$stability else NULL,
                        trees               = trees,
                        cluster_result      = result$cluster_result,
                        prior_strength      = prior_strength
                )
        }
}


#
# HELPERS (internal, not exported)
#

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

#' Knee-point detection via maximum perpendicular distance
#'
#' Locates the point on a curve furthest from the straight line joining its
#' endpoints (after `[0, 1]` normalisation of both axes). Robust to degenerate
#' inputs (no spread, coincident endpoints) — falls back to the last x.
#'
#' @param x,y Numeric vectors of equal length.
#' @return The \code{x} value at the knee.
#' @keywords internal
#' @noRd
.find_knee <- function(x, y) {
        n <- length(x)
        if (n <= 2L) return(x[n])
        
        rx <- max(x) - min(x)
        ry <- max(y) - min(y)
        if (rx <= 0) return(x[n])                      # no spread in x
        
        x_norm <- (x - min(x)) / rx
        y_norm <- (y - min(y)) / (ry + .Machine$double.eps)
        
        p1 <- c(x_norm[1L], y_norm[1L])
        p2 <- c(x_norm[n],  y_norm[n])
        line_len <- sqrt(sum((p2 - p1)^2))
        if (line_len <= .Machine$double.eps) return(x[n])  # coincident endpoints
        
        distances <- vapply(seq_len(n), function(i) {
                abs((x_norm[i] - p1[1L]) * (p2[2L] - p1[2L]) -
                            (y_norm[i] - p1[2L]) * (p2[1L] - p1[1L])) / line_len
        }, numeric(1))
        
        x[which.max(distances)]
}

#' Internal: compute the N x n_trees ensemble membership matrix
#'
#' Shared by the public \code{compute_ensemble_memberships()} and the stability
#' evaluator. For each tree, removes \code{intermediate_regions - 1} edges (weighted
#' sampling without replacement) and labels the resulting connected components.
#'
#' @param graph Parent \code{igraph} (used only to re-resolve edges when the
#'   stored tree object is a subgraph rather than an edge sequence).
#' @param trees A \code{"spanning_trees"} object.
#' @param n.rst Number of trees to use.
#' @param intermediate_regions Number of regions per tree.
#' @param prior_strength Within-prior halving rate.
#' @return Integer matrix N x n.rst of region labels.
#' @keywords internal
#' @noRd
.compute_ensemble_memberships_internal <- function(graph,
                                                   trees,
                                                   n.rst,
                                                   intermediate_regions,
                                                   prior_strength = 0) {
        n_use   <- as.integer(n.rst)
        n_reg   <- as.integer(intermediate_regions)
        n_v     <- trees$n_vertices
        weights <- trees$edge_weights
        ep_all  <- trees$edge_endpoints
        
        if (n_reg > n_v) {
                stop(sprintf("intermediate_regions (%d) cannot exceed the number of vertices (%d).",
                             n_reg, n_v), call. = FALSE)
        }
        
        # One cut per extra region (a tree on n_v vertices has n_v - 1 edges).
        n_to_cut <- n_reg - 1L
        
        # Prior-typology weighting (no-op when no prior or prior_strength == 0).
        weights <- .apply_prior_weighting(weights, trees$within_prior_edge, prior_strength)
        
        memb_mat <- matrix(0L, nrow = n_v, ncol = n_use)
        
        for (j in seq_len(n_use)) {
                tree <- trees$trees[[j]]
                
                # Resolve parent-graph edge IDs spanned by this tree. Modern igraph's
                # sample_spanning_tree() returns an edge sequence (igraph.es) whose integer
                # IDs already index into E(graph) — and therefore into `weights`/`ep_all`,
                # which were captured from E(graph) in sample_spanning_trees(). Older igraph
                # may return a subgraph; handle both, and error loudly on anything else.
                if (inherits(tree, "igraph.es")) {
                        edge_ids <- as.integer(tree)
                } else if (inherits(tree, "igraph")) {
                        tree_ep  <- igraph::ends(tree, igraph::E(tree), names = FALSE)
                        edge_ids <- .match_tree_edges(tree_ep, ep_all)
                } else {
                        stop("Unexpected spanning-tree object of class '", class(tree)[1L],
                             "'; expected an igraph edge sequence or graph.", call. = FALSE)
                }
                
                if (anyNA(edge_ids)) {
                        stop("Failed to map a spanning-tree edge back to the parent graph.",
                             call. = FALSE)
                }
                if (n_to_cut >= length(edge_ids)) {
                        stop(sprintf(
                                paste0("Requested %d cuts but the spanning tree has only %d edges; ",
                                       "reduce intermediate_regions or supply a connected graph."),
                                n_to_cut, length(edge_ids)
                        ), call. = FALSE)
                }
                
                # Cut probability proportional to edge weight: high-dissimilarity edges are
                # the most likely to be removed (the SKATER principle). Guard against
                # degenerate weight vectors so sample.int() never errors.
                w <- weights[edge_ids]
                if (anyNA(w) || any(!is.finite(w)) || any(w < 0)) {
                        stop("Edge weights must be finite and non-negative.", call. = FALSE)
                }
                if (all(w == 0)) w <- NULL  # fall back to uniform sampling
                
                cut_idx    <- sample.int(length(edge_ids), n_to_cut, prob = w, replace = FALSE)
                kept_edges <- edge_ids[-cut_idx]
                
                # Rebuild the graph from surviving edges. n = n_v keeps singleton vertices
                # so components() counts every node.
                ep <- ep_all[kept_edges, , drop = FALSE]
                g  <- igraph::make_graph(as.vector(t(ep)), n = n_v, directed = FALSE)
                memb_mat[, j] <- igraph::components(g)$membership
        }
        
        memb_mat
}

#' Match spanning-tree edges to parent-graph edge indices
#'
#' Builds canonical undirected edge keys (smaller vertex first) and matches.
#'
#' @param tree_ep,parent_ep Two-column integer matrices of edge endpoints.
#' @return Integer vector of parent-edge indices (NA where unmatched).
#' @keywords internal
#' @noRd
.match_tree_edges <- function(tree_ep, parent_ep) {
        make_key <- function(ep) {
                paste(pmin(ep[, 1], ep[, 2]), pmax(ep[, 1], ep[, 2]), sep = "-")
        }
        match(make_key(tree_ep), make_key(parent_ep))
}

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
#' @details Used by \code{.stability_fuzzy()} to score how well two hardened
#'   partitions agree on the overlap of resampled objects. We use the crisp ARI
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
        if (length(a) < 2L) return(NA_real_)
        tab <- table(a, b)
        n   <- sum(tab)
        if (n < 2L) return(NA_real_)
        sum_nij <- sum(choose(tab, 2L))
        ai <- rowSums(tab); bj <- colSums(tab)
        sum_ai <- sum(choose(ai, 2L))
        sum_bj <- sum(choose(bj, 2L))
        expected <- (sum_ai * sum_bj) / choose(n, 2L)
        max_idx  <- 0.5 * (sum_ai + sum_bj)
        denom    <- max_idx - expected
        # Degenerate case: both partitions are singletons or a single block.
        if (denom == 0) return(1)
        (sum_nij - expected) / denom
}

#' Stratified subsample by hard-cluster labels
#'
#' Proportional allocation with a hard floor of 2 objects per cluster, so that
#' silhouette / ARI computations on the subsample are not pathological for small
#' clusters.
#'
#' @param hard_clusters Integer vector of length N.
#' @param n_target Approximate total subsample size.
#' @return Sorted integer vector of selected row indices.
#' @keywords internal
#' @noRd
.stratified_subsample <- function(hard_clusters, n_target) {
        N <- length(hard_clusters)
        if (n_target >= N) return(seq_len(N))
        tab   <- table(hard_clusters)
        props <- as.numeric(tab) / N
        n_per <- pmax(2L, as.integer(round(props * n_target)))
        names(n_per) <- names(tab)
        groups <- names(tab)
        idx <- unlist(lapply(groups, function(g) {
                in_g <- which(as.character(hard_clusters) == g)
                take <- min(length(in_g), n_per[[g]])
                if (length(in_g) <= take) in_g else sample(in_g, take)
        }), use.names = FALSE)
        sort(idx)
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

#' Fuzzy Silhouette (Campello & Hruschka 2006)
#'
#' Weights the per-object Rousseeuw silhouette by the fuzziness margin
#' \eqn{(u_{ip} - u_{iq})^\alpha}. \eqn{\alpha = 1} (default) is the canonical
#' fuzzy generalisation; \eqn{\alpha \to \infty} reduces to the crisp
#' silhouette; \eqn{\alpha \to 0} is an unweighted average of \eqn{s_i}.
#'
#' @param D N x N dissimilarity (or \code{dist}) over the same N rows as U.
#' @param U N x k membership matrix.
#' @param alpha Non-negative numeric. Default 1.
#' @return Numeric scalar; higher is better. \code{NA_real_} if \code{k < 2}.
#' @keywords internal
#' @noRd
.fuzzy_silhouette <- function(D, U, alpha = 1) {
        N <- nrow(U)       # number of objects
        k <- ncol(U)       # number of classes
        D_size <- if (inherits(D, "dist")) attr(D, "Size") else nrow(D)
        if (N != D_size) stop("Dimensions of D and U disagree.", call. = FALSE)
        if (k < 2L) return(NA_real_)
        
        # Hard labels = highest-probability class per object.
        labels <- max.col(U, ties.method = "first")
        s_i    <- .rousseeuw_silhouette(D, labels)
        
        # Per-object top-two memberships. apply(..., sort, decreasing = TRUE) returns
        # a k x N matrix (each column is one sorted row of U); rows 1 and 2 are the
        # per-object top-1 and top-2 memberships.
        sorted_U <- apply(U, 1L, sort, decreasing = TRUE)
        if (!is.matrix(sorted_U)) return(NA_real_)  # extra defense (k == 1)
        margins <- sorted_U[1L, ] - sorted_U[2L, ]
        
        w <- margins ^ alpha
        if (sum(w) <= 0) return(NA_real_)
        sum(w * s_i) / sum(w)
}

#' Bootstrap stability of a fuzzy partition (mean ARI vs. reference)
#'
#' On each of \code{B} resamples, draws a row subsample of size
#' \code{floor(N * subsample_frac)} WITHOUT replacement, re-runs CLARA/PAM on
#' the subsampled ensemble matrix to obtain a hardened partition of the
#' subsample, and computes the ARI against the reference labels restricted to
#' those rows. Returns the mean ARI.
#'
#' @details This is the Hennig-style "compare each replicate to the full-data
#'   reference" recipe (one ARI per replicate), not pairwise replicate-vs-
#'   replicate (B^2 ARIs). It is B-fold cheaper and is what
#'   \code{fpc::clusterboot} does by default. Sampling is without replacement so
#'   ARI's contingency table is well-defined.
#'
#' @param memb_mat Integer matrix N x n_trees (the ensemble matrix).
#' @param hard_clusters Reference hardened labels of length N (from the
#'   full-data fit).
#' @param k Number of clusters.
#' @param fuzziness Fuzzifier; unused for hardened-label stability but kept in
#'   the signature so future fuzzy-ARI extensions need no API break.
#' @param B Number of bootstrap replicates. \code{0} disables.
#' @param subsample_frac Fraction of rows to retain per replicate. Default 0.8.
#' @param large_n_threshold CLARA/PAM cutoff (mirrors
#'   \code{cluster_consensus}).
#' @param clara_samples,clara_sampsize CLARA controls.
#' @param seed Optional integer seed; the caller's RNG state is restored on
#'   exit.
#' @param verbose Logical.
#' @return Numeric scalar mean ARI in approximately \code{[-0.5, 1]}, or
#'   \code{NA_real_} if \code{B = 0} or the subsample is too small.
#' @keywords internal
#' @noRd
.stability_fuzzy <- function(memb_mat,
                             hard_clusters,
                             k,
                             fuzziness = 2,
                             B = 25L,
                             subsample_frac = 0.8,
                             large_n_threshold = 20000L,
                             clara_samples = 20L,
                             clara_sampsize = NULL,
                             seed = NULL,
                             verbose = FALSE) {
        if (B <= 0L) return(NA_real_)
        N     <- nrow(memb_mat)
        n_sub <- floor(N * subsample_frac)
        if (n_sub < k * 2L) {
                if (verbose) message("    [stability] subsample too small; returning NA")
                return(NA_real_)
        }
        
        # RNG handling: seed locally, restore the caller's state on exit.
        if (!is.null(seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                set.seed(seed)
        }
        
        aris <- numeric(B)
        for (b in seq_len(B)) {
                sub_idx  <- sort(sample.int(N, size = n_sub, replace = FALSE))
                memb_sub <- memb_mat[sub_idx, , drop = FALSE]
                ref_sub  <- hard_clusters[sub_idx]
                if (nrow(memb_sub) > large_n_threshold) {
                        sampsize <- if (is.null(clara_sampsize)) {
                                min(nrow(memb_sub), 200L + 2L * k)
                        } else clara_sampsize
                        fit_b <- cluster::clara(memb_sub, k = k,
                                                metric   = "manhattan",
                                                samples  = clara_samples,
                                                sampsize = sampsize,
                                                pamLike  = TRUE)
                        new_lab <- fit_b$clustering
                } else {
                        coassoc_b <- parallelDist::parallelDist(memb_sub, method = "hamming")
                        fit_b     <- cluster::pam(coassoc_b, k = k)
                        new_lab   <- fit_b$clustering
                }
                aris[b] <- .adjusted_rand_index(new_lab, ref_sub)
        }
        mean(aris, na.rm = TRUE)
}

#' Borda rank-aggregation across sub-metric columns
#'
#' For each column of \code{df}, ranks rows (1 = best per the direction); NAs
#' receive the worst rank within that column; fully-NA columns are dropped;
#' remaining ranks are summed row-wise. The optimal row MINIMISES the returned
#' Borda score.
#'
#' @param df Data frame of sub-metric values (one column per voter).
#' @param directions Character vector of length \code{ncol(df)}; each element
#'   \code{"min"} or \code{"max"}.
#' @return Numeric vector of length \code{nrow(df)}; all-NA if every voter is NA.
#'
#' @details Uses \code{ties.method = "average"} so genuine ties do not falsely
#'   separate candidates. Weighted Borda is intentionally not used: each voter
#'   is a different *aspect* of validity (compactness/separation, fuzziness-
#'   aware silhouette, stability), and over-weighting one would re-introduce the
#'   bias direction we are trying to balance.
#'
#' @keywords internal
#' @noRd
.borda_score <- function(df, directions) {
        stopifnot(ncol(df) == length(directions))
        ranks_list <- mapply(function(values, dir) {
                if (all(is.na(values))) return(rep(NA_real_, length(values)))
                if (dir == "min")      rank(values,  na.last = "keep", ties.method = "average")
                else if (dir == "max") rank(-values, na.last = "keep", ties.method = "average")
                else stop("direction must be 'min' or 'max'", call. = FALSE)
        }, df, directions, SIMPLIFY = FALSE)
        ranks_mat <- do.call(cbind, ranks_list)
        keep_col  <- !apply(ranks_mat, 2L, function(x) all(is.na(x)))
        if (!any(keep_col)) return(rep(NA_real_, nrow(df)))
        ranks_mat <- ranks_mat[, keep_col, drop = FALSE]
        # NAs within a kept column become the worst rank in that column.
        for (j in seq_len(ncol(ranks_mat))) {
                na_mask <- is.na(ranks_mat[, j])
                if (any(na_mask)) {
                        ranks_mat[na_mask, j] <- max(ranks_mat[, j], na.rm = TRUE) + 1
                }
        }
        rowSums(ranks_mat)
}

#' Derive the \code{score} column of \code{tuning_log} from sub-metric columns
#'
#' Arranged so that \code{which.max(score)} always selects the best row:
#' minimise-is-better indices are negated, and Borda rank-sums are negated. For
#' metrics whose primary value lives directly in the log (e.g. SIL_F),
#' \code{score} is a copy.
#'
#' @param df Data frame containing the CVI columns.
#' @param metric Character metric name.
#' @return Numeric vector of scores (higher = better).
#' @keywords internal
#' @noRd
.derive_score_column <- function(df, metric) {
        switch(metric,
               silhouette                     = df$SILH_HARD,
               partition_coefficient          = df$PC,
               modified_partition_coefficient = df$MPC,
               MPC                            = df$MPC,
               PE                             = -df$PE,
               PEN                            = -df$PEN,
               XB_star                        = -df$XB_star,
               SIL_F                          = df$SIL_F,
               stability                      = df$STAB,
               borda                          = -.borda_score(
                       df         = df[, c("XB_star", "SIL_F", "STAB"), drop = FALSE],
                       directions = c("min", "max", "max")
               ),
               stop("Unknown metric: ", metric, call. = FALSE)
        )
}

#' Compute fuzzy-aware Cluster Validity Indices on a final fuzzy partition
#'
#' Returns a 1-row data frame of values, NOT a single scalar. The caller derives
#' a \code{score} column according to the chosen \code{metric} (see
#' \code{.derive_score_column}).
#'
#' @param U N x k membership matrix (rows sum to 1).
#' @param d_to_medoids N x k matrix of distances from each object to each
#'   medoid (from \code{cluster_consensus}).
#' @param medoid_dist k x k matrix of inter-medoid distances.
#' @param memb_mat Integer matrix N x n_trees. Required if \code{SIL_F} or
#'   \code{STAB} are to be computed; may be \code{NULL} otherwise (those entries
#'   become \code{NA}).
#' @param hard_clusters Integer vector of length N (argmax of U). If \code{NULL},
#'   derived from \code{U}.
#' @param fuzziness Fuzzifier \emph{m}; used in the XB* numerator weighting.
#' @param silhouette_hard Numeric scalar; hardened average silhouette width from
#'   CLARA/PAM (passed through as \code{SILH_HARD}).
#' @param stability_B Integer; bootstrap replicates for STAB. \code{0} skips it.
#' @param stability_subsample_frac Fraction of rows per stability replicate.
#' @param stability_seed Optional seed.
#' @param silf_subsample_n Target subsample size for SIL_F when N exceeds it.
#' @param silf_alpha Campello–Hruschka fuzziness exponent (default 1).
#' @param large_n_threshold CLARA/PAM cutoff for stability re-fits.
#' @param clara_samples,clara_sampsize CLARA controls for stability re-fits.
#' @param verbose Logical.
#' @return One-row \code{data.frame} with columns PC, PE, MPC, PEN, XB_star,
#'   SIL_F, STAB, SILH_HARD.
#' @keywords internal
#' @noRd
.compute_fuzzy_cvis <- function(U,
                                d_to_medoids,
                                medoid_dist,
                                memb_mat                 = NULL,
                                hard_clusters            = NULL,
                                fuzziness                = 2,
                                silhouette_hard          = NA_real_,
                                stability_B              = 0L,
                                stability_subsample_frac = 0.8,
                                stability_seed           = NULL,
                                silf_subsample_n         = 5000L,
                                silf_alpha               = 1,
                                large_n_threshold        = 20000L,
                                clara_samples            = 20L,
                                clara_sampsize           = NULL,
                                verbose                  = FALSE) {
        
        N <- nrow(U); k <- ncol(U)
        eps <- sqrt(.Machine$double.eps)
        
        # --- PC / PE / MPC / PEN (membership-only diagnostics) -------------------
        PC  <- mean(rowSums(U^2))
        PE  <- -mean(rowSums(U * log(pmax(U, eps))))
        MPC <- if (k > 1L) 1 - (k / (k - 1)) * (1 - PC) else NA_real_
        PEN <- if (k > 1L) PE / log(k) else NA_real_
        
        # --- XB* (Hamming-medoid Xie-Beni) ---------------------------------------
        # Numerator: sum_{i,k} u_ik^m * d_ik^2 (compactness via medoid distances).
        # Denominator: N * min_{k != l} D_kl^2 (separation via inter-medoid dist).
        # medoid_dist's diagonal is 0; mask with lower.tri to get k != l pairs.
        compact <- sum((U ^ fuzziness) * (d_to_medoids ^ 2))
        if (k >= 2L) {
                off_diag <- medoid_dist[lower.tri(medoid_dist)]
                sep_min  <- min(off_diag ^ 2)
                XB_star  <- compact / (N * max(sep_min, eps))
        } else {
                XB_star <- NA_real_
        }
        
        # --- SIL_F (Campello-Hruschka fuzzy silhouette) --------------------------
        # Requires a dissimilarity matrix over the SAME rows U indexes. For large N
        # this is intractable, so subsample stratified by argmax(U) and recompute.
        SIL_F <- NA_real_
        if (!is.null(memb_mat) && k >= 2L) {
                if (is.null(hard_clusters)) hard_clusters <- max.col(U, ties.method = "first")
                if (N > silf_subsample_n) {
                        sub_idx <- .stratified_subsample(hard_clusters, silf_subsample_n)
                        if (verbose) message(sprintf(
                                "    [SIL_F] stratified subsample: %d / %d rows", length(sub_idx), N
                        ))
                } else {
                        sub_idx <- seq_len(N)
                }
                memb_sub <- memb_mat[sub_idx, , drop = FALSE]
                U_sub    <- U[sub_idx, , drop = FALSE]
                D_sub    <- parallelDist::parallelDist(memb_sub, method = "hamming")
                SIL_F <- tryCatch(
                        .fuzzy_silhouette(D_sub, U_sub, alpha = silf_alpha),
                        error = function(e) {
                                if (verbose) message("    [SIL_F] failed: ", conditionMessage(e))
                                NA_real_
                        }
                )
        }
        
        # --- STAB (bootstrap ARI) ------------------------------------------------
        STAB <- NA_real_
        if (!is.null(memb_mat) && stability_B > 0L) {
                if (is.null(hard_clusters)) hard_clusters <- max.col(U, ties.method = "first")
                STAB <- tryCatch(
                        .stability_fuzzy(memb_mat          = memb_mat,
                                         hard_clusters     = hard_clusters,
                                         k                 = k,
                                         fuzziness         = fuzziness,
                                         B                 = stability_B,
                                         subsample_frac    = stability_subsample_frac,
                                         large_n_threshold = large_n_threshold,
                                         clara_samples     = clara_samples,
                                         clara_sampsize    = clara_sampsize,
                                         seed              = stability_seed,
                                         verbose           = verbose),
                        error = function(e) {
                                if (verbose) message("    [STAB] failed: ", conditionMessage(e))
                                NA_real_
                        }
                )
        }
        
        data.frame(PC = PC, PE = PE, MPC = MPC, PEN = PEN,
                   XB_star = XB_star, SIL_F = SIL_F, STAB = STAB,
                   SILH_HARD = silhouette_hard)
}

#' Compute Hamming distances from all nodes to reference (medoid) nodes
#'
#' Returns an N x k matrix giving the fraction of ensemble columns in which each
#' node disagrees with each reference node.
#'
#' @param memb_mat Integer matrix (N x n_trees).
#' @param ref_idx Integer vector of length k (row indices of medoids).
#' @return Numeric matrix (N x k) of Hamming distances (fractions in `[0, 1]`).
#'
#' @details Accumulates per-tree disagreements into the N x k counter, iterating
#'   over the n_trees columns. This keeps peak memory at O(N k) instead of
#'   materialising an N x n_trees logical matrix for every medoid — important
#'   when N reaches 1e5-1e6.
#'
#' @keywords internal
#' @noRd
.hamming_to_refs <- function(memb_mat, ref_idx) {
        n_trees     <- ncol(memb_mat)
        N           <- nrow(memb_mat)
        k           <- length(ref_idx)
        medoid_rows <- memb_mat[ref_idx, , drop = FALSE]   # k x n_trees
        
        out <- matrix(0, nrow = N, ncol = k)
        for (t in seq_len(n_trees)) {
                # outer(col, medoid_col, "!=") is an N x k logical, coerced to numeric in +.
                out <- out + outer(memb_mat[, t], medoid_rows[, t], FUN = "!=")
        }
        out / n_trees
}

#' Resolve a \code{prior_typology} argument to a vertex-label vector
#'
#' Accepts either a single character (vertex-attribute name on \code{graph}) or
#' a vector of length \code{vcount(graph)}. Coerces to character; empty strings
#' become NA. Returns NULL if input is NULL.
#'
#' @param graph An \code{igraph} object.
#' @param prior_typology See \code{\link{sample_spanning_trees}}.
#' @return Character vector of length \code{vcount(graph)}, or \code{NULL}.
#' @keywords internal
#' @importFrom igraph vcount vertex_attr vertex_attr_names
#' @noRd
.resolve_prior_typology <- function(graph, prior_typology) {
        if (is.null(prior_typology)) return(NULL)
        
        n_v <- igraph::vcount(graph)
        
        if (is.character(prior_typology) && length(prior_typology) == 1L) {
                if (!prior_typology %in% igraph::vertex_attr_names(graph)) {
                        stop(sprintf("Vertex attribute '%s' not found on graph. ", prior_typology),
                             "Attach the typology via V(graph)$", prior_typology,
                             " <- ... before calling.", call. = FALSE)
                }
                labs <- igraph::vertex_attr(graph, prior_typology)
        } else if (length(prior_typology) == n_v) {
                labs <- prior_typology
        } else {
                stop("`prior_typology` must be either a single character giving the name of ",
                     "a vertex attribute on `graph`, or a vector of length vcount(graph) = ",
                     n_v, ".", call. = FALSE)
        }
        
        labs <- as.character(labs)
        labs[!nzchar(labs)] <- NA_character_
        
        # Sanity check: at least 2 distinct non-NA labels for weighting to do anything.
        uniq <- unique(labs[!is.na(labs)])
        if (length(uniq) < 2L) {
                warning("Prior typology contains fewer than 2 distinct non-NA labels; ",
                        "weighting will have no effect.", call. = FALSE)
        }
        labs
}

#' Flag edges whose endpoints share a non-NA prior label
#'
#' @param prior_labels Character vector of vertex labels, or \code{NULL}.
#' @param edge_endpoints Two-column integer matrix of edge endpoints.
#' @return Logical vector of length \code{nrow(edge_endpoints)}, or \code{NULL}.
#' @keywords internal
#' @noRd
.compute_within_prior_edge <- function(prior_labels, edge_endpoints) {
        if (is.null(prior_labels)) return(NULL)
        a <- prior_labels[edge_endpoints[, 1]]
        b <- prior_labels[edge_endpoints[, 2]]
        !is.na(a) & !is.na(b) & a == b
}

#' Apply the within-prior halving rule to edge weights
#'
#' \eqn{\tilde{w}(i,j) = w(i,j) / 2^{\lambda_E}} when both endpoints share a
#' non-NA prior label; \eqn{\tilde{w}(i,j) = w(i,j)} otherwise. Thus
#' \eqn{\lambda_E = 0} returns weights unchanged; \eqn{\lambda_E = 1} halves
#' within-prior weights; \eqn{\lambda_E = 2} quarters them, etc. Lowering these
#' weights makes the corresponding edges less likely to be cut, so prior regions
#' tend to stay intact.
#'
#' @param weights Numeric edge-weight vector.
#' @param within_prior_edge Logical vector (same length), or \code{NULL}.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda_E}.
#' @return Numeric vector of (possibly down-weighted) edge weights.
#' @keywords internal
#' @noRd
.apply_prior_weighting <- function(weights, within_prior_edge, prior_strength) {
        if (is.null(within_prior_edge) || prior_strength == 0) return(weights)
        weights / (2 ^ (prior_strength * as.numeric(within_prior_edge)))
}

#' Exact reproduction of a prior typology (\code{prior_strength = Inf})
#'
#' Returns the prior labels directly as a (degenerate) clustering result, with a
#' length-N \code{coverage} flag marking vertices that have no prior label.
#' Uncovered vertices receive \code{NA} hard labels and all-zero membership rows.
#'
#' @param prior_labels Character vector of length N (resolved labels).
#' @param verbose Logical.
#' @return A list mirroring the \code{get_regions()} return, plus \code{coverage}
#'   and \code{class_levels}.
#' @keywords internal
#' @noRd
.exact_prior_result <- function(prior_labels, verbose = TRUE) {
        N         <- length(prior_labels)
        coverage  <- !is.na(prior_labels)
        n_covered <- sum(coverage)
        n_dropped <- N - n_covered
        
        if (n_covered == 0L) {
                stop("`prior_strength = Inf` but no vertices have a prior label.",
                     call. = FALSE)
        }
        
        if (verbose) {
                message("prior_strength = Inf: returning prior typology directly.")
                if (n_dropped > 0L) {
                        message(sprintf(
                                "  %d / %d vertices have no prior label and are marked NA (coverage = FALSE).",
                                n_dropped, N
                        ))
                }
        }
        
        classes <- sort(unique(prior_labels[coverage]))
        K       <- length(classes)
        
        hard_clusters           <- rep(NA_integer_, N)
        hard_clusters[coverage] <- match(prior_labels[coverage], classes)
        
        # One-hot membership for covered vertices; all-zero rows for uncovered ones.
        memberships <- matrix(0, nrow = N, ncol = K)
        memberships[cbind(which(coverage), hard_clusters[coverage])] <- 1
        
        list(
                memberships         = memberships,
                hard_clusters       = hard_clusters,
                coverage            = coverage,        # length-N logical
                best_intermediate_regions = K,
                best_final_regions  = K,
                best_n.rst          = NA_integer_,
                best_score          = NA_real_,
                tuning_log          = NULL,
                n.rst_stability     = NULL,
                trees               = NULL,
                cluster_result      = NULL,
                prior_strength      = Inf,
                class_levels        = classes
        )
}