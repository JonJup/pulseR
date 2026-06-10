
# Fuzzy co-occurrence signatures and clustering
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Compute k-order reachability matrix
#'
#' Computes the k-order reachability matrix from a graph's adjacency matrix,
#' indicating which nodes can reach each other within k steps.
#'
#' @param graph An \code{igraph} graph object. Typically created with
#'   \code{polygon_to_network()}.
#' @param k Integer scalar \eqn{\ge 1}. The neighborhood order (number of steps).
#'
#' @return A sparse binary \code{CsparseMatrix} of dimension \eqn{n \times n},
#'   where entry \eqn{(i, j) = 1} if node \eqn{j} is reachable from node \eqn{i}
#'   within \code{k} steps (including self).
#'
#' @details The reachability pattern is the non-zero pattern of \eqn{(A + I)^k},
#'   where \eqn{A} is the (binary) adjacency matrix and \eqn{I} the identity.
#'   All stored non-zero entries are set to 1 to give a binary indicator.
#'
#'   The identity is added as a \emph{sparse} matrix (\code{Matrix::Diagonal});
#'   using base \code{diag(n)} would coerce the whole computation to a dense
#'   \eqn{n \times n} object. Note that fill-in grows with \code{k}, so very
#'   large \code{k} on a dense graph can still produce a large matrix.
#'
#' @importFrom igraph as_adjacency_matrix is_igraph
#' @importFrom Matrix Diagonal drop0
#' @importFrom methods as
#' @seealso \code{\link{compute_all_signatures}}
#' @export
get_reachability <- function(graph, k) {
        if (!igraph::is_igraph(graph))
                stop("`graph` must be an igraph object.", call. = FALSE)
        if (length(k) != 1L || is.na(k) || k < 1 || k != as.integer(k))
                stop("`k` must be a single positive integer.", call. = FALSE)
        k <- as.integer(k)
        
        # Binary adjacency (attr = NULL ignores any edge weights), kept sparse.
        A <- igraph::as_adjacency_matrix(graph, sparse = TRUE, attr = NULL)
        n <- nrow(A)
        
        # Sparse identity
        AI <- A + Matrix::Diagonal(n)
        
        # Reachability within k steps == non-zero pattern of (A + I)^k. Repeated
        # sparse products keep the memory footprint bounded by the actual fill-in.
        Rk <- AI
        if (k > 1L) {
                for (step in seq_len(k - 1L)) {
                        Rk <- Rk %*% AI
                }
        }
        
        # Binarize. Coerce to column-compressed form, drop any structural zeros
        # (so we never flip a true zero to one), then set all stored values to 1.
        Rk <- methods::as(Rk, "CsparseMatrix")
        Rk <- Matrix::drop0(Rk)
        Rk@x <- rep.int(1, length(Rk@x))
        Rk
}


#' Compute fuzzy co-occurrence signature for a focal polygon
#'
#' Computes a normalized fuzzy co-occurrence signature for a focal polygon
#' based on its k-order neighborhood, summarizing the pairwise class
#' co-occurrence structure.
#'
#' @param focal Integer. Index (1-based) of the focal polygon / graph vertex.
#' @param P Numeric matrix (\eqn{n \times c}) of fuzzy membership values, where
#'   \code{n} is the number of polygons and \code{c} the number of classes.
#'   Row \code{i} must correspond to graph vertex \code{i}.
#' @param graph An \code{igraph} graph object. Used only for shortest-path
#'   distances when \code{kernel != "none"}.
#' @param Rk Sparse reachability matrix as returned by
#'   \code{\link{get_reachability}}.
#' @param A Optional precomputed sparse adjacency matrix
#'   (\code{igraph::as_adjacency_matrix(graph, sparse = TRUE)}). Supplying it
#'   avoids rebuilding the adjacency on every call -- \code{\link{compute_all_signatures}}
#'   always passes it. If \code{NULL}, it is recomputed from \code{graph}.
#' @param kernel Decay function for distance weighting. One of:
#'   \describe{
#'     \item{"none"}{Uniform weights (no decay). Baseline; no shortest-path
#'       computation is performed.}
#'     \item{"gaussian"}{Weights decay as \eqn{e^{-d^2 / 2\sigma^2}}.}
#'     \item{"exponential"}{Weights decay as \eqn{e^{-d / \sigma}} (heavier tails;
#'       suited to river-network propagation such as dispersal or dilution).}
#'     \item{"linear"}{Weights decay as \eqn{\max(0, 1 - d / \sigma)} with a hard
#'       cutoff at \eqn{d = \sigma}; \code{sigma} is the effective radius.}
#'   }
#'   Tighter kernels (small \code{sigma}) give more locally distinct signatures
#'   and patchier clusters; wider kernels smooth spatial variation. Default
#'   "none".
#' @param sigma Bandwidth parameter; must be \eqn{> 0} when \code{kernel != "none"}
#'   (ignored otherwise). Default 1.
#' @param digits Integer or \code{NULL}. If non-\code{NULL}, the signature is
#'   rounded to \code{digits} decimal places (default \code{3}, matching the
#'   original behaviour). Rounding is lossy; pass \code{NULL} to disable it.
#' @param coords A numeric \eqn{n \times 2} matrix of centroids, rows aligned to P)
#' @return Numeric vector of length \eqn{c(c+1)/2} -- the upper triangle
#'   (including diagonal) of the normalised fuzzy co-occurrence matrix. Returns
#'   all \code{NA} if the neighborhood contains fewer than 2 polygons.
#'
#' @details The signature is computed as
#'   \eqn{C = P_w^{\top} A_{local} P_w}, where \eqn{A_{local}} is the adjacency
#'   submatrix on the k-order neighborhood and \eqn{P_w} is the local membership
#'   matrix with each row scaled by its kernel weight \eqn{w_i}. Because the
#'   weighted membership enters the quadratic form twice, the contribution of an
#'   adjacent pair \eqn{(a, b)} is scaled by \eqn{w_a w_b}: both endpoints are
#'   down-weighted by their distance to the focal polygon. This is intentional
#'   (pairs nearest the focal dominate the signature) -- change with care. The
#'   upper triangle of \eqn{C} is then normalised to a probability distribution.
#'
#'   Reachability of an undirected adjacency graph is symmetric, so the focal
#'   polygon's neighbors are read directly from the stored entries of column
#'   \code{focal} of \code{Rk} when it is a \code{CsparseMatrix} (fast); a dense
#'   row lookup is used as a fallback.
#'
#' @importFrom igraph as_adjacency_matrix distances
#' @importFrom Matrix crossprod
#' @importFrom methods is
#' @seealso \code{\link{get_reachability}}, \code{\link{compute_all_signatures}}
#' @export
compute_signature <- function(focal, P, graph, Rk, A = NULL, coords = NULL,
                              kernel = c("none", "gaussian", "exponential", "linear"),
                              sigma  = NULL,
                              digits = 3L) {
        kernel  <- match.arg(kernel)
        nc      <- ncol(P)
        sig_len <- nc * (nc + 1L) / 2L
        
        # ---- Neighborhood lookup (includes the focal itself, via the diagonal) ----
        if (methods::is(Rk, "CsparseMatrix")) {
                # Use the column-compressed structure directly; column `focal`'s stored
                # rows are exactly its k-order neighbors (reachability is symmetric).
                p <- Rk@p
                nb_idx <- if (p[focal + 1L] > p[focal]) {
                        Rk@i[(p[focal] + 1L):p[focal + 1L]] + 1L
                } else {
                        integer(0)
                }
        } else {
                nb_idx <- which(Rk[focal, ] > 0)
        }
        
        # A neighborhood of < 2 polygons carries no pairwise co-occurrence signal.
        if (length(nb_idx) < 2L) {
                return(rep(NA_real_, sig_len))
        }
        
        # ---- Adjacency: reuse the caller's matrix; only rebuild when standalone ----
        if (is.null(A)) {
                A <- igraph::as_adjacency_matrix(graph, sparse = TRUE, attr = NULL)
        }
        
        # ---- Distance weights -----------------------------------------------------
        # Spatial (Euclidean) distance between the focal centroid and each neighbor
        # centroid, replacing the previous graph shortest-path distance. `graph` is
        # now used only for the adjacency fallback above, not for distances.
        if (kernel == "none") {
                w <- rep.int(1, length(nb_idx))
        } else {
                if (is.null(coords))
                        stop("`coords` (centroid coordinates) are required when kernel != 'none'.",
                             call. = FALSE)
                
                # Centroid distances from the focal to each neighbour.
                diffs <- sweep(coords[nb_idx, , drop = FALSE], 2, coords[focal, ], "-")
                dists <- sqrt(rowSums(diffs^2))
                
                # Resolve bandwidth: median-heuristic default when sigma is NULL/"auto",
                # otherwise the supplied positive value.
                if (is.null(sigma) || identical(sigma, "auto")) {
                        sigma <- .default_sigma(dists)
                } else if (length(sigma) != 1L || !is.numeric(sigma) || is.na(sigma) || sigma <= 0) {
                        stop("`sigma` must be a single positive number, NULL, or \"auto\".",
                             call. = FALSE)
                }
                
                w <- switch(kernel,
                            gaussian    = exp(-dists^2 / (2 * sigma^2)),
                            exponential = exp(-dists / sigma),
                            linear      = pmax(0, 1 - dists / sigma))
        }
        
        # ---- Local membership / adjacency and weighting ---------------------------
        P_local <- P[nb_idx, , drop = FALSE]
        A_local <- A[nb_idx, nb_idx, drop = FALSE]
        P_w     <- sweep(P_local, 1, w, "*")              # scale each polygon's row by w_i
        
        # Fuzzy co-occurrence C = P_w^T A_local P_w (c x c, symmetric).
        C <- as.matrix(Matrix::crossprod(P_w, A_local %*% P_w))
        # Intratype neighborhoods are counted twice (from each direction).
        # Halving restores mass balance. 
        diag(C) <- diag(C)/2
        # Upper triangle (incl. diagonal) flattened to the signature vector.
        sig <- C[upper.tri(C, diag = TRUE)]
        
        # Normalise to a probability distribution (comparable under Jensen-Shannon).
        s <- sum(sig)
        if (s > 0) sig <- sig / s
        
        if (!is.null(digits)) sig <- round(sig, digits)
        as.numeric(sig)
}


#' Compute fuzzy co-occurrence signatures for all polygons
#'
#' Iterates over all polygons to compute their fuzzy co-occurrence signatures
#' based on a k-order spatial neighborhood.
#'
#' @param P Numeric matrix (\eqn{n \times c}) of fuzzy membership values. Rows
#'   must correspond, in order, to the vertices of \code{graph}.
#' @param graph An \code{igraph} graph object representing spatial adjacency.
#' @param k Integer scalar \eqn{\ge 1}. The neighborhood order.
#' @param verbose Logical. If \code{TRUE} (default), progress messages print.
#' @param kernel,sigma Distance-weighting kernel and bandwidth; see
#'   \code{\link{compute_signature}}. Defaults "none" / 1.
#' @param digits Rounding for signatures; see \code{\link{compute_signature}}.
#'   Default 3.
#' @param coords A numeric \eqn{n \times 2} matrix of centroids, rows aligned to P)
#'
#' @return A numeric matrix (\eqn{n \times sig\_len}) of signatures, where
#'   \eqn{sig\_len = c(c+1)/2}. Column names encode class pairs (e.g.
#'   \code{"c1:c2"}). Columns that are zero across all polygons are removed.
#'   Rows for isolated / near-isolated polygons contain \code{NA}.
#'
#' @details The adjacency matrix is computed once here and passed to
#'   \code{\link{compute_signature}} for every polygon (the original code rebuilt
#'   it on each call). The number of \code{P} rows is checked against the vertex
#'   count, since the method relies on row \code{i} of \code{P} corresponding to
#'   vertex \code{i} of \code{graph}.
#'
#' @importFrom igraph as_adjacency_matrix is_igraph vcount
#' @seealso \code{\link{compute_signature}}, \code{\link{get_reachability}}
#' @export
compute_all_signatures <- function(P, graph, k, verbose = TRUE, coords = NULL,
                                   kernel = c("none", "gaussian", "exponential", "linear"),
                                   sigma  = 1,
                                   digits = 3L) {
        kernel <- match.arg(kernel)
        
        # ---- Input validation -----------------------------------------------------
        if (!igraph::is_igraph(graph))
                stop("`graph` must be an igraph object.", call. = FALSE)
        P <- as.matrix(P)
        if (!is.numeric(P))
                stop("`P` must be a numeric membership matrix.", call. = FALSE)
        n  <- nrow(P)
        nc <- ncol(P)
        if (nc < 2L)
                stop("`P` must have at least two class columns.", call. = FALSE)
        if (n != igraph::vcount(graph))
                stop(sprintf(paste0("Row count of `P` (%d) must equal the number of graph ",
                                    "vertices (%d); rows of `P` must correspond, in order, ",
                                    "to graph vertices."),
                             n, igraph::vcount(graph)), call. = FALSE)
        
        sig_len <- nc * (nc + 1L) / 2L
        
        # Centroids are needed only when a distance kernel is in use.
        if (kernel != "none") {
                if (is.null(coords))
                        stop("`coords` (centroid coordinates) are required when kernel != 'none'.",
                             call. = FALSE)
                coords <- as.matrix(coords)
                if (!is.numeric(coords))
                        stop("`coords` must be numeric centroid coordinates.", call. = FALSE)
                if (nrow(coords) != n)
                        stop(sprintf("`coords` must have one row per polygon (%d); got %d.",
                                     n, nrow(coords)), call. = FALSE)
        }
        
        # Adjacency computed ONCE and reused for every polygon.
        A <- igraph::as_adjacency_matrix(graph, sparse = TRUE, attr = NULL)
        
        if (verbose) message("Computing ", k, "-order reachability for ", n, " polygons...")
        Rk <- get_reachability(graph, k)
        
        if (verbose) message("Computing signatures...")
        sigs <- matrix(NA_real_, nrow = n, ncol = sig_len)
        
        # Column names encode the (ordered) class pairs, e.g. "c1:c2".
        cnames <- colnames(P)
        if (is.null(cnames)) cnames <- paste0("c", seq_len(nc))
        pair_names <- character(sig_len)
        idx <- 1L
        for (i in seq_len(nc)) {
                for (j in i:nc) {
                        pair_names[idx] <- paste0(cnames[i], ":", cnames[j])
                        idx <- idx + 1L
                }
        }
        colnames(sigs) <- pair_names
        
        for (i in seq_len(n)) {
                sigs[i, ] <- compute_signature(i, P, graph, Rk, A = A, coords = coords,
                        kernel = kernel, sigma = sigma, digits = digits)
                if (verbose && i %% 1000L == 0L) message("  ", i, "/", n)
        }
        
        if (verbose) {
                n_na <- sum(is.na(sigs[, 1]))
                message("Done. ", n_na, " polygons with NA signatures (isolated/near-isolated).")
        }
        
        # Drop class-pair columns that are zero for every (non-NA) polygon.
        # na.rm = TRUE is essential: with even one isolated (NA) polygon the column
        # sums would otherwise be NA and the original `if (any(sum(x) == 0))` test
        # errored with "missing value where TRUE/FALSE needed".
        col_zero <- colSums(sigs, na.rm = TRUE) == 0
        if (any(col_zero)) {
                sigs <- sigs[, !col_zero, drop = FALSE]
        }
        
        sigs
}


#' Compute Jensen-Shannon distance matrix
#'
#' Computes pairwise Jensen-Shannon distances between signature vectors,
#' excluding rows with missing values.
#'
#' @param sigs Numeric matrix of signatures (\eqn{n \times sig\_len}), as
#'   returned by \code{\link{compute_all_signatures}}.
#'
#' @return A \code{dist} object of pairwise Jensen-Shannon distances (the square
#'   root of the divergence, which is a proper metric) over all complete cases.
#'
#' @details A tiny epsilon is added and each row renormalised to a distribution
#'   to avoid \eqn{\log(0)}. The divergence is clamped at 0 before the square
#'   root to absorb floating-point underflow. \strong{Note:} the result is over
#'   the complete-case subset only, and the underlying distance matrix is dense
#'   (\eqn{O(N^2)} memory); this stage is intended for landscape-scale \code{N},
#'   not the full polygon set.
#'
#' @importFrom philentropy distance
#' @importFrom stats complete.cases
#' @seealso \code{\link{compute_all_signatures}}, \code{\link{get_mosaic_types}}
#' @export
js_distance_fast <- function(sigs) {
        sigs  <- as.matrix(sigs)
        valid <- stats::complete.cases(sigs)
        if (sum(valid) < 2L)
                stop("Fewer than two complete signature rows; cannot compute a distance matrix.",
                     call. = FALSE)
        
        S <- sigs[valid, , drop = FALSE]
        
        # Add a tiny constant and renormalise so no probability is exactly 0
        # (avoids log(0)) while keeping each row a valid distribution.
        eps <- .Machine$double.eps
        S <- S + eps
        S <- S / rowSums(S)
        
        # philentropy expects rows = distributions and returns the divergence.
        d <- philentropy::distance(S, method = "jensen-shannon", unit = "log",
                                   as.dist.obj = TRUE, use.row.names = TRUE)
        
        # Clamp tiny negatives from rounding, then take sqrt (sqrt(JSD) is a metric).
        # `d[] <-` preserves the dist attributes (Size, Labels, class).
        d[] <- sqrt(pmax(as.numeric(d), 0))
        d
}


#' Compute fuzzy cluster memberships via PAM
#'
#' Performs PAM (Partitioning Around Medoids) clustering and derives fuzzy
#' membership values from distances to medoids using a fuzzifier. Also returns
#' the per-object distance-to-medoid matrix and the inter-medoid distance
#' matrix, which the fuzzy validity indices reuse.
#'
#' @param coassoc A \code{dist} object or a symmetric dissimilarity matrix.
#'   Anything that is not already a \code{dist} is coerced with
#'   \code{\link[stats]{as.dist}} so PAM never treats it as raw data.
#' @param k Integer \eqn{\ge 2}. Number of clusters.
#' @param fuzziness Numeric \eqn{> 1}. Fuzzifier \emph{m}. Default 2.
#'
#' @return A list with \code{memberships} (\eqn{n \times k}), \code{pam_result},
#'   \code{d_to_medoids} (\eqn{n \times k}), \code{medoid_dist} (\eqn{k \times k}),
#'   and \code{fuzziness}.
#'
#' @importFrom cluster pam
#' @importFrom stats as.dist
#' @export
compute_fuzzy <- function(coassoc, k, fuzziness = 2) {
        # Accept a dist OR a dissimilarity matrix; coerce so pam() always uses
        # dissimilarities (diss = TRUE) rather than interpreting a matrix as data.
        if (!inherits(coassoc, "dist")) coassoc <- stats::as.dist(coassoc)
        
        if (length(k) != 1L || is.na(k) || k < 2 || k != as.integer(k))
                stop("`k` must be a single integer >= 2.", call. = FALSE)
        if (length(fuzziness) != 1L || is.na(fuzziness) || fuzziness <= 1)
                stop("`fuzziness` (m) must be a single number > 1.", call. = FALSE)
        k <- as.integer(k)
        
        pam_res   <- cluster::pam(coassoc, k = k, diss = TRUE)
        dist_full <- as.matrix(coassoc)
        
        dists_to_med <- dist_full[, pam_res$id.med, drop = FALSE]
        medoid_dist  <- dist_full[pam_res$id.med, pam_res$id.med, drop = FALSE]
        
        m <- fuzziness
        memberships <- t(apply(dists_to_med, 1, function(d) {
                # Exactly on a medoid -> hard membership of 1 to the nearest medoid.
                if (any(d == 0)) {
                        u <- rep(0, length(d)); u[which.min(d)] <- 1; return(u)
                }
                # Standard FCM-style membership derived from distances to medoids.
                d_pow <- d ^ (-2 / (m - 1))
                d_pow / sum(d_pow)
        }))
        
        list(memberships  = memberships,
             pam_result   = pam_res,
             d_to_medoids = dists_to_med,
             medoid_dist  = medoid_dist,
             fuzziness    = m)
}


#' Cluster signatures with automatic selection of cluster number
#'
#' Clusters polygons by their fuzzy co-occurrence signatures using PAM with
#' Jensen-Shannon distance, selecting the number of clusters with a fuzzy-aware
#' rule: per candidate \code{k} it computes \code{XB_star} (min), \code{SIL_F}
#' (max) and a bootstrap stability ARI (max), then selects \code{k} by Borda
#' rank-consensus. The full validity profile is returned. With
#' \code{crisp = TRUE} the fuzzy indices degenerate, so selection falls back to
#' the crisp silhouette.
#'
#' @param sigs Signature matrix (\eqn{n \times sig\_len}).
#' @param n_clusters Candidate \code{k} values (all \eqn{\ge 2}). Default 2:10.
#' @param method Currently only "pam".
#' @param crisp Logical. Default \code{FALSE}.
#' @param metric Character: "borda" (default), "XB_star", "SIL_F", "stability",
#'   or "silhouette" (legacy crisp).
#' @param stability_B Bootstrap replicates for stability. Default 25; 0 skips.
#' @param silf_alpha Campello-Hruschka exponent. Default 1.
#' @param seed Optional seed. The caller's RNG state is restored on exit.
#'
#' @return A list with \code{clusters} (final fuzzy or crisp result),
#'   \code{best_k}, \code{validity} (data.frame), \code{best} (alias of
#'   \code{best_k}), and \code{valid} -- a logical vector mapping the clustered
#'   rows back to the rows of \code{sigs} (the complete cases).
#'
#' @importFrom cluster pam
#' @importFrom stats complete.cases
#' @seealso \code{\link{compute_fuzzy}}, \code{\link{js_distance_fast}}
#' @export
get_mosaic_types <- function(sigs,
                               n_clusters  = 2:10,
                               method      = c("pam"),
                               crisp       = FALSE,
                               metric      = c("borda", "XB_star", "SIL_F",
                                               "stability", "silhouette"),
                               stability_B = 25L,
                               silf_alpha  = 1,
                               seed        = NULL) {
        method <- match.arg(method)
        metric <- match.arg(metric)
        
        sigs  <- as.matrix(sigs)
        valid <- stats::complete.cases(sigs)
        if (sum(valid) < 2L)
                stop("Fewer than two complete signature rows to cluster.", call. = FALSE)
        S <- sigs[valid, , drop = FALSE]
        
        n_clusters <- sort(unique(as.integer(n_clusters)))
        if (any(is.na(n_clusters)) || any(n_clusters < 2L))
                stop("`n_clusters` must all be integers >= 2.", call. = FALSE)
        if (max(n_clusters) >= nrow(S))
                stop(sprintf("`n_clusters` must be < the number of clustered objects (%d).",
                             nrow(S)), call. = FALSE)
        
        # Reproducibility: restore the caller's RNG state on exit so the bootstrap
        # stability draws do not silently perturb downstream randomness.
        if (!is.null(seed) && exists(".Random.seed", envir = .GlobalEnv)) {
                old_seed <- get(".Random.seed", envir = .GlobalEnv)
                on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
        }
        
        message("Computing Jensen-Shannon distances for ", nrow(S), " polygons...")
        d <- js_distance_fast(S)          # S is already complete; inner filter is a no-op
        D_full <- as.matrix(d)
        
        if (crisp && metric != "silhouette") {
                message("Note: fuzzy metrics are degenerate with crisp = TRUE; using 'silhouette'.")
                metric <- "silhouette"
        }
        
        message("Evaluating cluster solutions (metric = ", metric, ") ...")
        rows <- vector("list", length(n_clusters))
        for (i in seq_along(n_clusters)) {
                nc <- n_clusters[i]
                if (crisp) {
                        pr <- cluster::pam(d, k = nc, diss = TRUE)
                        rows[[i]] <- data.frame(
                                k = nc, PC = 1, PE = 0, MPC = 1,
                                XB_star = NA_real_, SIL_F = NA_real_, STAB = NA_real_,
                                SILH_HARD = pr$silinfo$avg.width
                        )
                } else {
                        fz <- compute_fuzzy(d, nc)
                        rows[[i]] <- .mosaic_cvis(fz, D_full,
                                                     stability_B = stability_B,
                                                     silf_alpha  = silf_alpha,
                                                     seed        = seed)
                }
                message(sprintf("  k=%d | XB*=%.4f SIL_F=%.4f STAB=%.4f (hard sil=%.3f)",
                                nc, rows[[i]]$XB_star, rows[[i]]$SIL_F,
                                rows[[i]]$STAB, rows[[i]]$SILH_HARD))
        }
        validity <- do.call(rbind, rows)
        
        # Selection score (higher = better). Guard against an all-NA score column,
        # which would make which.max() return integer(0) and error downstream.
        validity$score <- .derive_score_column(validity, metric)
        if (all(is.na(validity$score)))
                stop("Selection metric '", metric, "' produced only NA scores; cannot pick k.",
                     call. = FALSE)
        best_idx <- which.max(validity$score)
        best_k   <- validity$k[best_idx]
        
        message("Optimal clusters: ", best_k, " (", metric,
                " score = ", round(validity$score[best_idx], 4), ")")
        
        final <- if (crisp) cluster::pam(d, k = best_k, diss = TRUE) else compute_fuzzy(d, best_k)
        
        list(clusters = final,
             best_k   = best_k,
             validity = validity,
             best     = best_k,
             valid    = valid)   # maps clustered rows back to the rows of `sigs`
}


#' Sweep over neighborhood orders and cluster numbers
#'
#' Evaluates each neighborhood order \code{k} and reports the (hard) silhouette
#' width of every candidate number of clusters, plus the selected best \code{k}
#' per order.
#'
#' @param P Numeric matrix (\eqn{n \times c}) of fuzzy membership values.
#' @param graph An \code{igraph} graph object (spatial adjacency). \strong{Note:}
#'   this must be an igraph object, not a raw adjacency matrix -- the original
#'   argument name \code{A} was misleading.
#' @param k_range Integer vector of neighborhood orders to evaluate (default
#'   \code{1:5}).
#' @param n_clusters Integer vector of candidate cluster counts (default
#'   \code{2:10}).
#' @param kernel,sigma Distance-weighting kernel and bandwidth passed to
#'   \code{\link{compute_all_signatures}}. Defaults "none" / 1.
#' @param coords A numeric \eqn{n \times 2} matrix of centroids, rows aligned to P)
#'
#' @return A \code{data.frame} with columns \code{k_order}, \code{n_clusters},
#'   \code{silhouette} (average hard silhouette width for that solution), and
#'   \code{best} (the selected number of clusters for that \code{k_order}). The
#'   overall best combination is also reported via a message.
#'
#' @details Iterates over \code{k_range}, computing signatures with
#'   \code{\link{compute_all_signatures}} and evaluating clustering with
#'   \code{\link{get_mosaic_types}}.
#'
#' @importFrom igraph is_igraph
#' @seealso \code{\link{compute_all_signatures}}, \code{\link{get_mosaic_types}}
#' @export
sweep_k_order <- function(P, graph, coords = NULL, k_range = 1:5, n_clusters = 2:10,
                          kernel = "none", sigma = 1) {
        if (!igraph::is_igraph(graph))
                stop("`graph` must be an igraph object.", call. = FALSE)
        
        k_range    <- sort(unique(as.integer(k_range)))
        n_clusters <- sort(unique(as.integer(n_clusters)))
        
        # Collect rows in a list (avoids the quadratic cost of rbind-in-a-loop).
        out <- vector("list", length(k_range))
        
        for (ki in seq_along(k_range)) {
                k <- k_range[ki]
                message("\n=== Neighborhood order k = ", k, " ===")
                sigs <- compute_all_signatures(P, graph, k, coords = coords, verbose = TRUE,
                                               kernel = kernel, sigma = sigma)
                cl   <- get_mosaic_types(sigs, n_clusters = n_clusters)
                
                # cl$validity has one row per candidate k (same order as the candidates);
                # SILH_HARD is the average hard silhouette for that solution. The original
                # code used `cl$validity[i]`, which indexes a *column* of the data.frame.
                out[[ki]] <- data.frame(
                        k_order    = k,
                        n_clusters = cl$validity$k,
                        silhouette = cl$validity$SILH_HARD,
                        best       = cl$best_k
                )
        }
        
        results <- do.call(rbind, out)
        
        best <- results[which.max(results$silhouette), ]
        message("\n=== Best overall: k_order=", best$k_order,
                ", n_clusters=", best$n_clusters,
                ", silhouette=", round(best$silhouette, 3), " ===")
        
        results
}


# .............................................................................
# Stability by subsetting the dissimilarity matrix and re-PAM-ing
# .............................................................................

#' Bootstrap stability for a dissimilarity-based PAM partition (mean ARI)
#'
#' On each of \code{B} resamples, subsets the full dissimilarity matrix to a
#' row/col subsample (without replacement) and re-runs PAM, then compares the
#' hardened labels to the reference labels restricted to that subsample via ARI.
#' Cheaper than recomputing JS distances because it slices the precomputed
#' matrix.
#'
#' @param D_full Full \eqn{N \times N} dissimilarity matrix.
#' @param hard_clusters Reference labels of length \eqn{N} (PAM on full data).
#' @param k Number of clusters.
#' @param B Bootstrap replicates. \eqn{\le 0} disables (returns \code{NA}).
#' @param subsample_frac Fraction of objects per replicate. Default 0.8.
#' @param seed Optional seed.
#'
#' @return Mean ARI scalar, or \code{NA_real_}.
#'
#' @importFrom cluster pam
#' @importFrom stats as.dist
#' @keywords internal
.stability_pam_diss <- function(D_full, hard_clusters, k,
                                B = 25L, subsample_frac = 0.8, seed = NULL) {
        if (B <= 0L) return(NA_real_)
        N <- nrow(D_full)
        n_sub <- floor(N * subsample_frac)
        if (n_sub < k * 2L) return(NA_real_)        # too few objects for a stable k-partition
        if (!is.null(seed)) set.seed(seed)          # RNG state restored by the caller
        aris <- numeric(B)
        for (b in seq_len(B)) {
                idx     <- sort(sample.int(N, n_sub, replace = FALSE))
                d_sub   <- stats::as.dist(D_full[idx, idx])
                lab     <- cluster::pam(d_sub, k = k, diss = TRUE)$clustering
                aris[b] <- .adjusted_rand_index(lab, hard_clusters[idx])
        }
        mean(aris, na.rm = TRUE)
}


# .............................................................................
# Fuzzy CVIs on a compute_fuzzy() result
# .............................................................................

#' Compute fuzzy-aware validity indices for one mosaic clustering
#'
#' @param fz Output of \code{\link{compute_fuzzy}}.
#' @param D_full Full \eqn{N \times N} dissimilarity matrix (JS distances over
#'   complete cases) -- the same one PAM was run on.
#' @param stability_B Bootstrap replicates for STAB. 0 to skip.
#' @param silf_alpha Campello-Hruschka exponent. Default 1.
#' @param silf_subsample_n If \eqn{N} exceeds this, SIL_F is computed on a
#'   stratified subsample. For this pipeline \eqn{N} is usually small enough to
#'   skip subsampling.
#' @param seed Optional seed (stability).
#'
#' @return One-row data.frame: \code{k}, \code{PC}, \code{PE}, \code{MPC},
#'   \code{XB_star}, \code{SIL_F}, \code{STAB}, \code{SILH_HARD}.
#'
#' @importFrom stats as.dist
#' @keywords internal
.mosaic_cvis <- function(fz, D_full,
                            stability_B = 25L, silf_alpha = 1,
                            silf_subsample_n = 8000L, seed = NULL) {
        U <- fz$memberships
        N <- nrow(U); k <- ncol(U)
        eps <- sqrt(.Machine$double.eps)
        
        # Partition coefficient, partition entropy, modified PC.
        PC  <- mean(rowSums(U^2))
        PE  <- -mean(rowSums(U * log(pmax(U, eps))))
        MPC <- if (k > 1L) 1 - (k / (k - 1)) * (1 - PC) else NA_real_
        
        # XB* in JS-distance / medoid space.
        compact <- sum((U ^ fz$fuzziness) * (fz$d_to_medoids ^ 2))
        off     <- fz$medoid_dist[lower.tri(fz$medoid_dist)]
        XB_star <- if (k >= 2L) compact / (N * max(min(off^2), eps)) else NA_real_
        
        # SIL_F directly on the full JS matrix (subsample only if very large).
        if (N > silf_subsample_n) {
                labels  <- max.col(U, ties.method = "first")
                sub_idx <- .stratified_subsample(labels, silf_subsample_n)
                D_sub   <- stats::as.dist(D_full[sub_idx, sub_idx])
                SIL_F   <- .fuzzy_silhouette(D_sub, U[sub_idx, , drop = FALSE], alpha = silf_alpha)
        } else {
                SIL_F <- .fuzzy_silhouette(stats::as.dist(D_full), U, alpha = silf_alpha)
        }
        
        STAB <- .stability_pam_diss(D_full, fz$pam_result$clustering, k,
                                    B = stability_B, seed = seed)
        
        data.frame(k = k, PC = PC, PE = PE, MPC = MPC,
                   XB_star = XB_star, SIL_F = SIL_F, STAB = STAB,
                   SILH_HARD = fz$pam_result$silinfo$avg.width)
}

#' Centroid coordinates for a set of polygons
#'
#' Convenience helper that extracts representative point coordinates (centroids
#' or points-on-surface) for a polygon layer, returning the n x 2 numeric matrix
#' that \code{\link{compute_all_signatures}} / \code{\link{compute_signature}}
#' expect as \code{coords}.
#'
#' @param polygons An \code{sf}/\code{sfc} polygon layer, or a \code{terra}
#'   \code{SpatVector} of polygons. Feature order is preserved, so the returned
#'   rows align with the rows of \code{P} and the vertices of the adjacency
#'   graph -- provided all three were built from this same layer, in this order.
#' @param type Representative point. \code{"centroid"} (default) uses the true
#'   centroid; \code{"point_on_surface"} guarantees a point inside the polygon,
#'   which is safer for concave or multipart catchments whose centroid can fall
#'   outside the geometry (e.g. in a river bend or between sub-basins).
#' @param warn_longlat Logical. If \code{TRUE} (default), warn when the layer is
#'   in geographic (lon/lat) coordinates, since the kernels use planar
#'   (Euclidean) distance. Project to a metric CRS (e.g. LAEA, EPSG:3035 for
#'   Europe) first.
#'
#' @return A numeric \eqn{n \times 2} matrix with columns \code{c("x", "y")},
#'   one row per input feature (order preserved). Row names are taken from the
#'   layer's feature IDs / row names when available.
#'
#' @seealso \code{\link{compute_all_signatures}}, \code{\link{compute_signature}}
#' @importFrom sf st_as_sf st_geometry st_centroid st_point_on_surface st_coordinates st_is_longlat
#' @export
polygon_to_centroids <- function(polygons,
                              type = c("centroid", "point_on_surface"),
                              warn_longlat = TRUE) {
        type <- match.arg(type)
        
        # Single code path: coerce terra SpatVector -> sf.
        if (inherits(polygons, "SpatVector")) {
                if (!requireNamespace("sf", quietly = TRUE))
                        stop("Package 'sf' is required.", call. = FALSE)
                polygons <- sf::st_as_sf(polygons)
        }
        if (!inherits(polygons, c("sf", "sfc")))
                stop("`polygons` must be an sf/sfc layer or a terra SpatVector.",
                     call. = FALSE)
        
        if (isTRUE(warn_longlat) && isTRUE(sf::st_is_longlat(polygons)))
                warning("`polygons` is in geographic (lon/lat) coordinates; the distance ",
                        "kernels assume planar distance. Project to a metric CRS ",
                        "(e.g. EPSG:3035) before computing signatures.", call. = FALSE)
        
        # st_geometry() avoids the "attribute variables assumed constant" warning.
        pts <- switch(type,
                      centroid         = sf::st_centroid(sf::st_geometry(polygons)),
                      point_on_surface = sf::st_point_on_surface(sf::st_geometry(polygons)))
        
        xy <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
        colnames(xy) <- c("x", "y")
        
        rn <- rownames(polygons)
        if (!is.null(rn)) rownames(xy) <- rn
        xy
}

#' Data-driven default bandwidth for the distance kernels (internal)
#'
#' Returns a characteristic neighbourhood distance to use as the kernel
#' bandwidth \code{sigma} when none is supplied: the median of the strictly
#' positive focal-to-neighbour centroid distances (median heuristic). The kernel
#' then "sees" structure at the typical local spacing. It is locally adaptive --
#' each focal polygon gets its own scale -- so denser regions get tighter
#' kernels.
#'
#' @param dists Numeric vector of centroid distances from the focal polygon to
#'   its k-order neighbours (derived from the coordinate matrix). The focal's
#'   own 0 and any non-finite entries are ignored.
#' @param fallback Positive bandwidth returned when no strictly positive
#'   distance exists (all neighbour centroids coincide with the focal). In that
#'   degenerate case every kernel weight is 1 regardless, so the value only
#'   needs to be > 0. Default 1.
#'
#' @return A single positive numeric bandwidth.
#'
#' @importFrom stats median
#' @keywords internal
.default_sigma <- function(dists, fallback = 1) {
        d <- dists[is.finite(dists) & dists > 0]
        if (length(d) == 0L) return(fallback)
        stats::median(d)
}
