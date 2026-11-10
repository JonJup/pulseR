# =============================================================================
# pulseR — LEVEL 3: mosaic types from fuzzy co-occurrence signatures
#
# INPUT      as_mosaic_input() adapts pulse_regions / pulse_local_types /
#            pulse_meta_types / bare matrices, and records which level and
#            which class labels the signatures are built from.
# SPEED      The signature of focal f is the sum of per-edge class-pair outer
#            products over edges inside its neighbourhood, so the whole matrix
#            is  crossprod(B, M)  with B[e, f] = w_f(u_e) w_f(v_e).
# SCALE      Distances are O(N^2), so clustering runs on a sample of
#            `sample_size` rows and assigns the rest by the fitted model
#            (DMM) or by distance to medoids (PAM).
# METRIC     sqrt-JSD saturates at log 2 when two signatures have disjoint
#            support; Hellinger, total variation, Bray-Curtis and Euclidean
#            are also available.
#
# Requires: Matrix, igraph, cluster, stats.
# =============================================================================


# -----------------------------------------------------------------------------*
# 1. Input adapters
# -----------------------------------------------------------------------------*

#' Membership matrix for the mosaic level
#'
#' Adapts the output of any earlier level into the \eqn{n \times c} class
#' membership matrix the signature machinery needs, and records what the columns
#' mean.
#'
#' @param x One of:
#'   \describe{
#'     \item{\code{"pulse_meta_types"}}{Cross-region meta-types. The usual input
#'       once merging is in play: meta-types recur across regions, so
#'       neighbourhoods in different regions can share classes.}
#'     \item{\code{"pulse_local_types"}}{Region-nested local types. Note that a
#'       local type determines its region, so signatures built on these encode
#'       region identity and the mosaic level will tend to recover the region
#'       map.}
#'     \item{\code{"pulse_regions"} / \code{"pulse_regions_tuning"}}{Level-1 regions.}
#'     \item{A numeric matrix}{Used as is.}
#'   }
#' @param use \code{"memberships"} (default; rows sum to one),
#'   \code{"cohesion"} (the fit's co-association cohesion, level 1 only) or
#'   \code{"weighted"} (the un-normalised \eqn{z \cdot w} of a local-type fit).
#' @param level \code{"as_is"} (default) or \code{"region"}. For local types,
#'   \code{"region"} sums each region's local-type columns back into one column.
#'   That aggregation is exact: local-type rows sum to one within a region, so
#'   the column sums reproduce the region weights. Not available for meta-types,
#'   which have no single parent region.
#' @param normalise Logical. Rescale rows to sum to one. The signature is
#'   bilinear in \eqn{P}, so rows that do not sum to one mean a polygon
#'   contributes proportionally less mass to its neighbours' signatures —
#'   defensible as an uncertainty weighting, but a modelling choice rather than
#'   a neutral one, and it interacts with the distance kernel. Default
#'   \code{FALSE} (leave the matrix as it came).
#' @param drop_empty Drop classes with no mass anywhere. They contribute only
#'   zero bins and inflate the signature dimension. Default \code{TRUE}.
#' @return An object of class \code{"pulse_mosaic_input"}: a list with \code{P},
#'   \code{labels}, \code{region_of} (parent region per column, \code{NA} where
#'   the class has none) and \code{level}.
#' @export
as_mosaic_input <- function(x,
                        use       = c("memberships", "cohesion", "weighted"),
                        level     = c("as_is", "region"),
                        normalise = FALSE,
                        drop_empty = TRUE) {
        use   <- match.arg(use)
        level <- match.arg(level)
        
        # --- meta-types ------------------------------------------------------*
        if (inherits(x, "pulse_meta_types")) {
                if (level == "region")
                        stop("`level = \"region\"` is undefined for meta-types: a meta-type ",
                             "spans several regions by construction. Aggregate the ",
                             "pulse_local_types object instead.", call. = FALSE)
                if (use != "memberships")
                        stop("Only `use = \"memberships\"` is defined for meta-types.",
                             call. = FALSE)
                ## region_of stays NA on purpose. Downstream code keys the smoothing
                ## kernel on it, and a fabricated parent region would reintroduce
                ## exactly the nesting the merge was meant to remove.
                return(.new_mosaic_input(x$memberships, colnames(x$memberships),
                                     rep(NA_integer_, ncol(x$memberships)), "meta",
                                     normalise, drop_empty))
        }
        
        # --- local types -----------------------------------------------------*
        if (inherits(x, "pulse_local_types")) {
                P <- switch(use,
                            memberships = x$memberships,
                            weighted    = x$weighted,
                            cohesion    = stop("`use = \"cohesion\"` is a level-1 quantity; ",
                                               "pass the pulse_regions object.", call. = FALSE))
                P   <- as.matrix(P)
                lb  <- x$key$name
                reg <- as.integer(x$key$region)
                
                if (level == "region") {
                        ## Numeric factor levels: split() on an integer vector would
                        ## group by as.character() and put region 10 before region 2.
                        lv  <- sort(unique(reg))
                        idx <- split(seq_len(ncol(P)), factor(reg, levels = lv))
                        P   <- vapply(idx, function(j) rowSums(P[, j, drop = FALSE]),
                                      numeric(nrow(P)))
                        return(.new_mosaic_input(P, paste0("R", lv), lv, "region",
                                             normalise, drop_empty))
                }
                return(.new_mosaic_input(P, lb, reg, "local", normalise, drop_empty))
        }
        
        # --- regions ---------------------------------------------------------*
        if (inherits(x, "pulse_regions") || inherits(x, "pulse_regions_tuning")) {
                src <- if (inherits(x, "pulse_regions_tuning")) x$best_result else x
                if (is.null(src))
                        stop("`pulse_regions_tuning` object carries no `$best_result`.", call. = FALSE)
                P <- switch(use,
                            memberships = .first_non_null(src$memberships, src$fit$memberships),
                            cohesion    = .first_non_null(src$cohesion, src$fit$cohesion),
                            weighted    = stop("`use = \"weighted\"` is a level-2 quantity.",
                                               call. = FALSE))
                if (is.null(P)) stop("Could not find a `", use, "` matrix on `x`.", call. = FALSE)
                P <- as.matrix(P)
                return(.new_mosaic_input(P, paste0("R", seq_len(ncol(P))),
                                     seq_len(ncol(P)), "region", normalise, drop_empty))
        }
        
        # --- bare matrix -----------------------------------------------------*
        if (is.matrix(x) || is.data.frame(x)) {
                P  <- as.matrix(x)
                lb <- colnames(P)
                if (is.null(lb)) lb <- paste0("c", seq_len(ncol(P)))
                return(.new_mosaic_input(P, lb, rep(NA_integer_, ncol(P)), "matrix",
                                     normalise, drop_empty))
        }
        
        stop("`x` must be a pulse_meta_types, pulse_local_types, pulse_regions, ",
             "pulse_regions_tuning, or a numeric matrix.", call. = FALSE)
}



# -----------------------------------------------------------------------------*
# 2. Reachability and neighbourhood weights
# -----------------------------------------------------------------------------*

#' k-order reachability matrix
#'
#' @param graph An \code{igraph} object.
#' @param k_order Integer \eqn{\ge 1}, neighbourhood order.
#' @return Sparse binary symmetric \code{CsparseMatrix}: the non-zero pattern of
#'   \eqn{(A + I)^k}.
#' @importFrom igraph as_adjacency_matrix is_igraph
#' @importFrom Matrix Diagonal drop0
#' @importFrom methods as
#' @export
get_reachability <- function(graph, k_order) {
        if (!igraph::is_igraph(graph)) stop("`graph` must be an igraph object.", call. = FALSE)
        if (length(k_order) != 1L || is.na(k_order) || k_order < 1 || k_order != as.integer(k_order))
                stop("`k_order` must be a single positive integer.", call. = FALSE)
        k <- as.integer(k_order)
        
        A  <- igraph::as_adjacency_matrix(graph, sparse = TRUE, attr = NULL)
        AI <- A + Matrix::Diagonal(nrow(A))
        Rk <- AI
        if (k > 1L) for (step in seq_len(k - 1L)) Rk <- Rk %*% AI
        
        Rk <- methods::as(Matrix::drop0(Rk), "CsparseMatrix")
        Rk@x <- rep.int(1, length(Rk@x))
        Rk
}

#' Kernel weights on the stored neighbourhood pattern
#'
#' Returns a sparse matrix with the pattern of \code{Rk} whose entry
#' \eqn{(u, f)} is the kernel weight of neighbour \eqn{u} for focal \eqn{f}.
#'
#' @param Rk Reachability matrix from \code{\link{get_reachability}}.
#' @param coords Numeric \eqn{n \times 2} centroid matrix.
#' @param kernel \code{"none"}, \code{"gaussian"}, \code{"exponential"} or
#'   \code{"linear"}.
#' @param sigma Bandwidth. A positive number (map units — the interpretable
#'   choice), \code{"global"} for the median of all stored positive
#'   focal-neighbour distances, or \code{"adaptive"} for a per-focal median.
#' @details
#' \strong{Prefer a fixed \code{sigma} in map units.} An adaptive per-focal
#' bandwidth makes the spatial extent of a "landscape" vary with unit
#' density, so two signatures can differ because the drainage network is denser
#' in one place, not because the configuration of types differs. That confounds
#' the very thing the mosaic level is supposed to measure. A fixed radius (the
#' 25--50 km scale from the work-package design is the natural starting point)
#' keeps the unit comparable across Europe. \code{"adaptive"} is kept because it
#' reproduces the previous default.
#' @return Sparse \code{CsparseMatrix}, same pattern as \code{Rk}.
#' @importFrom Matrix sparseMatrix
#' @export
compute_neighbourhood_weights <- function(Rk, coords,
                                  kernel = c("none", "gaussian", "exponential", "linear"),
                                  sigma  = "global") {
        kernel <- match.arg(kernel)
        n <- nrow(Rk)
        if (kernel == "none") return(Rk)
        if (is.null(coords)) stop("`coords` are required when kernel != \"none\".", call. = FALSE)
        coords <- as.matrix(coords)
        if (nrow(coords) != n) stop("`coords` must have one row per vertex.", call. = FALSE)
        
        ## Vectorised distances over the stored pattern (CSC: column = focal).
        row_i <- Rk@i + 1L
        col_j <- rep.int(seq_len(n), diff(Rk@p))
        dd    <- sqrt(rowSums((coords[row_i, , drop = FALSE] -
                                       coords[col_j, , drop = FALSE])^2))
        
        s <- if (identical(sigma, "global")) {
                pos <- dd[dd > 0 & is.finite(dd)]
                if (!length(pos)) 1 else stats::median(pos)
        } else if (identical(sigma, "adaptive") || identical(sigma, "auto")) {
                med <- vapply(split(dd, col_j), .default_sigma, numeric(1))
                med[as.character(col_j)]
        } else {
                if (length(sigma) != 1L || !is.numeric(sigma) || is.na(sigma) || sigma <= 0)
                        stop("`sigma` must be a positive number, \"global\" or \"adaptive\".",
                             call. = FALSE)
                sigma
        }
        
        w <- switch(kernel,
                    gaussian    = exp(-dd^2 / (2 * s^2)),
                    exponential = exp(-dd / s),
                    linear      = pmax(0, 1 - dd / s))
        
        W <- Rk
        W@x <- w
        Matrix::drop0(W)
}

#' @keywords internal
#' @noRd
.default_sigma <- function(dists, fallback = 1) {
        d <- dists[is.finite(dists) & dists > 0]
        if (!length(d)) return(fallback)
        stats::median(d)
}


# -----------------------------------------------------------------------------*
# 3. Signatures — vectorised over edges
# -----------------------------------------------------------------------------*

#' Fuzzy co-occurrence signatures for every polygon
#'
#' For focal \eqn{f} the signature is the upper triangle of
#' \eqn{C^{(f)} = P_w^\top A_{N(f)} P_w}, normalised to a probability
#' distribution. Rather than looping over focals, the identity
#' \deqn{C^{(f)} = \sum_{e \in N(f)} w_f(u_e)\,w_f(v_e)\, m_e}
#' is used, where \eqn{m_e} is the class-pair contribution of edge \eqn{e}. The
#' whole matrix is then one sparse product, which is orders of magnitude faster
#' than the per-focal implementation at continental \eqn{N}.
#'
#' @param mosaic_input Numeric \eqn{n \times c} membership matrix, rows in graph vertex
#'   order, or the list returned by \code{\link{as_mosaic_input}}.
#' @param graph An \code{igraph} adjacency graph.
#' @param k_order Neighbourhood order.
#' @param coords Centroid matrix; required when \code{kernel != "none"}.
#' @param kernel,sigma See \code{\link{compute_neighbourhood_weights}}.
#' @param drop_zero_cols Drop class pairs that never occur. Default \code{TRUE}.
#'   Note this makes signature dimension depend on the data, so matrices from
#'   two different \code{k_order} are not column-aligned; \code{\link{check_k_order}}
#'   handles that by clustering each \code{k} separately.
#' @param block Class-pair columns processed per sparse product. Lower it if
#'   memory is tight. Default 64.
#' @param digits Rounding, or \code{NULL} (default) for none. Rounding is lossy
#'   and was only ever cosmetic.
#' @param exclude_bridges Drop edges flagged \code{bridge == TRUE} (set by
#'   \code{\link{build_network}}) before building neighbourhoods, so
#'   neighbourhoods never cross a bridge. Default \code{TRUE}. A graph without
#'   a \code{bridge} edge attribute is used unchanged.
#' @param verbose Logical.
#' @return Numeric \eqn{n \times \binom{c+1}{2}} matrix with \code{"a:b"} column
#'   names and attribute \code{"exclude_bridges"}. Rows whose neighbourhood
#'   contains no edge are \code{NA}; with \code{exclude_bridges = TRUE} this
#'   includes units joined to the rest only by bridges.
#' @importFrom igraph as_adjacency_matrix ends E is_igraph vcount
#' @importFrom Matrix crossprod
#' @export
compute_all_signatures <- function(mosaic_input, graph, k_order,
                                   coords = NULL,
                                   kernel = c("none", "gaussian", "exponential", "linear"),
                                   sigma  = "global",
                                   drop_zero_cols = TRUE,
                                   block   = 64L,
                                   digits  = NULL,
                                   exclude_bridges = TRUE,
                                   verbose = TRUE) {
        kernel <- match.arg(kernel)
        P <- mosaic_input
        k <- k_order
        if (!is.logical(exclude_bridges) || length(exclude_bridges) != 1L || is.na(exclude_bridges))
                stop("`exclude_bridges` must be TRUE or FALSE.", call. = FALSE)
        
        lab <- NULL
        if (is.list(P) && !is.null(P$P)) { lab <- P$labels; P <- P$P }
        P <- as.matrix(P)
        if (!is.numeric(P)) stop("`mosaic_input` must be numeric.", call. = FALSE)
        if (!igraph::is_igraph(graph)) stop("`graph` must be an igraph object.", call. = FALSE)
        n <- nrow(P)
        nc <- ncol(P)
        if (nc < 2L) stop("`mosaic_input` must have at least two class columns.", call. = FALSE)
        if (n != igraph::vcount(graph))
                stop(sprintf("`mosaic_input` has %d rows but the graph has %d vertices.",
                             n, igraph::vcount(graph)), call. = FALSE)
        if (anyNA(P)) stop("`mosaic_input` contains NA.", call. = FALSE)
        if (is.null(lab)) lab <- colnames(P)
        if (is.null(lab)) lab <- paste0("c", seq_len(nc))
        
        if (exclude_bridges) graph <- .drop_bridges(graph, verbose)
        if (verbose) message(sprintf("Reachability (k = %d) for %d polygons ...", k, n))
        Rk <- get_reachability(graph, k)
        W  <- compute_neighbourhood_weights(Rk, coords, kernel, sigma)
        
        ep <- igraph::ends(graph, igraph::E(graph), names = FALSE)
        ep <- ep[ep[, 1] != ep[, 2], , drop = FALSE]      # self-loops carry no pair signal
        u <- ep[, 1]; v <- ep[, 2]
        ne <- nrow(ep)
        if (!ne) stop("`graph` has no edges.", call. = FALSE)
        
        ## B[e, f] = w_f(u_e) * w_f(v_e): non-zero only when BOTH endpoints are in
        ## N(f), which is exactly the "edge lies inside the neighborhood" test.
        if (verbose) message("Building edge-by-focal incidence ...")
        B <- W[u, , drop = FALSE] * W[v, , drop = FALSE]   # ne x n, sparse
        ## Kernel-weighted mean number of neighborhoods containing each edge.
        ## Every edge enters the DMM likelihood this many times.
        eB <- Matrix::rowSums(B)
        edge_mult <- if (any(eB > 0)) sum(eB) / sum(eB > 0) else NA_real_
        ## Class-pair index (upper triangle, diagonal included).
        pairs <- which(upper.tri(matrix(0, nc, nc), diag = TRUE), arr.ind = TRUE)
        pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
        sig_len <- nrow(pairs)
        pair_names <- paste0(lab[pairs[, 1]], ":", lab[pairs[, 2]])
        
        if (verbose) message(sprintf("Accumulating %d class pairs over %d edges ...",
                                     sig_len, ne))
        signatures <- matrix(0, n, sig_len, dimnames = list(NULL, pair_names))
        for (start in seq(1L, sig_len, by = block)) {
                cols <- start:min(sig_len, start + block - 1L)
                a <- pairs[cols, 1]; b <- pairs[cols, 2]
                ## Edge contribution. Off-diagonal: the two ways of assigning classes
                ## a and b to the endpoints. Diagonal: a single p_ua p_va, which is the
                ## "halve the diagonal" correction of the matrix form, done once here
                ## rather than patched afterwards.
                M <- P[u, a, drop = FALSE] * P[v, b, drop = FALSE]
                off <- a != b
                if (any(off))
                        M[, off] <- M[, off] + P[u, b[off], drop = FALSE] * P[v, a[off], drop = FALSE]
                signatures[, cols] <- as.matrix(Matrix::crossprod(B, M))
        }
        
        rs <- rowSums(signatures)
        empty <- !is.finite(rs) | rs <= 0
        signatures <- signatures / rs
        signatures[empty, ] <- NA_real_
        if (!is.null(digits)) signatures <- round(signatures, digits)
        
        if (verbose)
                message(sprintf("Done. %d polygon(s) with NA signatures (no edge in neighbourhood).",
                                sum(empty)))
        
        if (drop_zero_cols) {
                keep <- colSums(signatures, na.rm = TRUE) > 0
                if (any(!keep)) {
                        if (verbose) message(sprintf("Dropping %d never-occurring class pair(s).",
                                                     sum(!keep)))
                        signatures <- signatures[, keep, drop = FALSE]
                }
        }
        rs[empty] <- NA_real_
        structure(signatures, class = c("pulse_signatures", "matrix", "array"),
                  n_classes = nc, class_labels = lab, k_order = k,
                  kernel = kernel, sigma = sigma, exclude_bridges = exclude_bridges,
                  totals = rs, edge_multiplicity = edge_mult)
}

#' Signature of one focal polygon (reference implementation)
#'
#' Retained for testing and teaching. \code{\link{compute_all_signatures}} does
#' not call it.
#' @param focal Integer vertex index.
#' @param P Membership matrix.
#' @param A Sparse adjacency matrix.
#' @param W Weight matrix from \code{\link{compute_neighbourhood_weights}}.
#' @return Numeric signature vector, or all \code{NA}.
#' @keywords internal
#' @noRd
.compute_signature <- function(focal, P, A, W) {
        nb <- which(W[, focal] != 0)
        nc <- ncol(P); sig_len <- nc * (nc + 1) / 2
        if (length(nb) < 2L) return(rep(NA_real_, sig_len))
        Pw <- P[nb, , drop = FALSE] * as.numeric(W[nb, focal])
        C  <- as.matrix(Matrix::crossprod(Pw, A[nb, nb, drop = FALSE] %*% Pw))
        diag(C) <- diag(C) / 2
        sig <- C[upper.tri(C, diag = TRUE)]
        s <- sum(sig)
        if (s > 0) sig / s else rep(NA_real_, sig_len)
}


# -----------------------------------------------------------------------------*
# 4. Class-aware ground kernel and saturation
# -----------------------------------------------------------------------------*

#' Smooth signatures over a class-similarity ground metric
#'
#' Jensen-Shannon, Hellinger and every other bin-wise divergence treat class
#' pairs as exchangeable labels: "type 3 beside type 4" and "type 3 beside type
#' 40" are equally different from "type 3 beside type 3". Once the class count
#' grows, neighbourhood signatures have almost disjoint support and every
#' distance sits at the metric's ceiling — the partition is then noise. See
#' \code{\link{check_signature_saturation}}.
#'
#' Smoothing the signature over a similarity kernel on class pairs fixes this at
#' the source and costs one matrix product. It is a linearised stand-in for an
#' optimal-transport distance with the same ground cost, without the dependency
#' or the per-pair Sinkhorn iterations.
#'
#' @param signatures Signature matrix.
#' @param class_similarity \eqn{c \times c} similarity in \eqn{[0,1]} with unit
#'   diagonal. If \code{NULL}, built from \code{region_of}.
#' @param region_of Integer parent region per class (from
#'   \code{\link{as_mosaic_input}}). Classes in the same region get similarity
#'   \code{within_region}.
#' @param within_region Similarity between two local types of the same region.
#'   Default 0.5.
#' @param renormalise Renormalise rows to sum to one. Default \code{TRUE}.
#' @return Smoothed signature matrix, same dimensions.
#' @export
smooth_signatures <- function(signatures, class_similarity = NULL, region_of = NULL,
                              within_region = 0.5, renormalise = TRUE) {
        lab <- attr(signatures, "class_labels")
        pr  <- .parse_pair_names(colnames(signatures), lab)
        
        if (is.null(class_similarity)) {
                if (is.null(region_of))
                        stop("Supply `class_similarity` or `region_of`.", call. = FALSE)
                Sc <- outer(region_of, region_of, function(a, b) ifelse(a == b, within_region, 0))
                diag(Sc) <- 1
        } else {
                Sc <- as.matrix(class_similarity)
        }
        
        a <- pr$a; b <- pr$b
        ## K[(a,b),(a',b')] = mean of the two ways of matching the unordered pairs.
        K <- (Sc[a, a, drop = FALSE] * Sc[b, b, drop = FALSE] +
                      Sc[a, b, drop = FALSE] * Sc[b, a, drop = FALSE]) / 2
        out <- signatures %*% K
        if (renormalise) out <- out / rowSums(out)
        dimnames(out) <- dimnames(signatures)
        attributes(out) <- c(attributes(out), attributes(signatures)[
                setdiff(names(attributes(signatures)), names(attributes(out)))])
        out
}

#' @keywords internal
#' @noRd
.parse_pair_names <- function(cn, lab) {
        sp <- strsplit(cn, ":", fixed = TRUE)
        a  <- vapply(sp, `[`, character(1), 1L)
        b  <- vapply(sp, `[`, character(1), 2L)
        list(a = match(a, lab), b = match(b, lab), a_lab = a, b_lab = b)
}

#' Fraction of signature pairs at the metric ceiling
#'
#' The diagnostic for the failure mode above. Above roughly 0.3 the distance
#' matrix is mostly constant and any partition of it is arbitrary; smooth the
#' signatures, aggregate to a coarser level, or raise \code{k}.
#'
#' @param d A \code{dist} object.
#' @param method The metric used, for the ceiling. See
#'   \code{\link{compute_signature_distance}}.
#' @param thresh Fraction of the ceiling counted as saturated. Default 0.95.
#' @return Numeric scalar in \eqn{[0,1]}, or \code{NA} for unbounded metrics.
#' @export
check_signature_saturation <- function(d, method = "js", thresh = 0.95) {
        ceil <- switch(method,
                       js        = sqrt(log(2)),
                       hellinger = 1,
                       tv        = 1,
                       bray      = 1,
                       NA_real_)
        if (!is.finite(ceil)) return(NA_real_)
        mean(as.numeric(d) > thresh * ceil, na.rm = TRUE)
}


# -----------------------------------------------------------------------------*
# 5. Distances
# -----------------------------------------------------------------------------*

#' Pairwise distance between signatures
#'
#' @param signatures Signature matrix; complete cases only are used.
#' @param method One of
#'   \describe{
#'     \item{\code{"js"}}{Square root of the Jensen-Shannon divergence (nats).
#'       A proper metric, ceiling \eqn{\sqrt{\log 2}}. The default, and what the
#'       method was published against.}
#'     \item{\code{"hellinger"}}{\eqn{\|\sqrt p - \sqrt q\|_2/\sqrt 2}. A proper
#'       metric with ceiling 1, behaving almost identically to \code{"js"} —
#'       its real advantage is that it is Euclidean distance on
#'       \code{\link{embed_hellinger}(signatures)}, so k-means, CLARA, PCA and every
#'       other Euclidean tool applies without an \eqn{N \times N} matrix.}
#'     \item{\code{"tv"}}{Total variation, \eqn{\tfrac12\|p-q\|_1}.}
#'     \item{\code{"bray"}}{Bray-Curtis. On row-normalised signatures this
#'       equals total variation; offered because ecologists expect it.}
#'     \item{\code{"euclidean"}}{Plain \eqn{L_2}. Not a divergence; useful only
#'       as a null comparison.}
#'   }
#' @param signatures_ref Optional reference matrix. When supplied, returns the
#'   \eqn{N \times R} matrix of distances to those references instead of a full
#'   \code{dist} — this is what makes out-of-sample assignment cheap.
#' @return A \code{dist}, or a numeric matrix when \code{signatures_ref} is given.
#' @importFrom stats complete.cases dist
#' @export
compute_signature_distance <- function(signatures, method = c("js", "hellinger", "tv", "bray", "euclidean"),
                               signatures_ref = NULL) {
        method <- match.arg(method)
        S <- as.matrix(signatures)
        eps <- .Machine$double.eps
        
        if (method %in% c("hellinger", "euclidean")) {
                X <- if (method == "hellinger") embed_hellinger(S) else S
                if (is.null(signatures_ref)) return(stats::dist(X))
                R <- if (method == "hellinger") embed_hellinger(as.matrix(signatures_ref)) else as.matrix(signatures_ref)
                return(.cross_euclidean(X, R))
        }
        
        Sn <- (S + eps) / rowSums(S + eps)
        if (method %in% c("tv", "bray")) {
                if (is.null(signatures_ref)) return(stats::dist(Sn, method = "manhattan") / 2)
                Rn <- as.matrix(signatures_ref); Rn <- (Rn + eps) / rowSums(Rn + eps)
                return(.cross_l1(Sn, Rn) / 2)
        }
        
        ## Jensen-Shannon, computed directly (no philentropy dependency, and the
        ## previous call rebuilt an N x N object twice over).
        if (is.null(signatures_ref)) {
                # Built directly: no Euclidean dist as an attribute placeholder.
                return(structure(.js_pairwise(Sn), Size = nrow(Sn), Labels = rownames(Sn),
                                 Diag = FALSE, Upper = FALSE, method = "js",
                                 class = "dist"))
        }
        Rn <- as.matrix(signatures_ref); Rn <- (Rn + eps) / rowSums(Rn + eps)
        .js_cross(Sn, Rn)
}

#' Hellinger embedding of a signature matrix
#'
#' Euclidean distance on the returned matrix equals the Hellinger distance on
#' the input, so any Euclidean clustering method becomes available at
#' \eqn{O(N)} memory.
#' @param signatures Signature matrix.
#' @return Numeric matrix of the same shape.
#' @export
embed_hellinger <- function(signatures) {
        S <- as.matrix(signatures)
        S <- S / rowSums(S)
        sqrt(S) / sqrt(2)
}

#' @keywords internal
#' @noRd
.ent <- function(p) { z <- p * log(p); z[!(p > 0)] <- 0; -rowSums(z) }

#' @keywords internal
#' @noRd
.js_pairwise <- function(P) {
        n <- nrow(P)
        H <- .ent(P)
        out <- numeric(n * (n - 1) / 2)
        pos <- 1L
        for (i in seq_len(n - 1L)) {
                j  <- (i + 1L):n
                Mx <- (P[j, , drop = FALSE] + rep(P[i, ], each = length(j))) / 2
                jsd <- .ent(Mx) - (H[i] + H[j]) / 2
                out[pos:(pos + length(j) - 1L)] <- sqrt(pmax(jsd, 0))
                pos <- pos + length(j)
        }
        out
}

#' @keywords internal
#' @noRd
.js_cross <- function(P, R) {
        HP <- .ent(P); HR <- .ent(R)
        out <- matrix(0, nrow(P), nrow(R))
        for (r in seq_len(nrow(R))) {
                Mx <- (P + rep(R[r, ], each = nrow(P))) / 2
                out[, r] <- sqrt(pmax(.ent(Mx) - (HP + HR[r]) / 2, 0))
        }
        out
}

#' @keywords internal
#' @noRd
.cross_euclidean <- function(X, R) {
        d2 <- outer(rowSums(X^2), rowSums(R^2), "+") - 2 * tcrossprod(X, R)
        sqrt(pmax(d2, 0))
}

#' @keywords internal
#' @noRd
.cross_l1 <- function(P, R) {
        out <- matrix(0, nrow(P), nrow(R))
        for (r in seq_len(nrow(R))) out[, r] <- rowSums(abs(P - rep(R[r, ], each = nrow(P))))
        out
}



# -----------------------------------------------------------------------------*
# 6. Fuzzy memberships from a medoid partition
# -----------------------------------------------------------------------------*

#' Fuzzy memberships from distances to medoids
#'
#' @param d_to_medoids Numeric \eqn{N \times k} distance matrix.
#' @param fuzziness Numeric \eqn{> 1}.
#' @return Numeric \eqn{N \times k} membership matrix, rows sum to one.
#' @keywords internal
#' @noRd
.memberships_from_medoids <- function(d_to_medoids, fuzziness = 2) {
        D <- as.matrix(d_to_medoids)
        n <- nrow(D); k <- ncol(D)
        U <- matrix(0, n, k)
        dp  <- D ^ (-2 / (fuzziness - 1))
        bad <- (rowSums(D == 0) > 0L) | !is.finite(rowSums(dp))
        if (any(!bad)) {
                g <- dp[!bad, , drop = FALSE]
                U[!bad, ] <- g / rowSums(g)
        }
        if (any(bad)) {
                i <- which(bad)
                U[cbind(i, max.col(-D[i, , drop = FALSE], ties.method = "first"))] <- 1
        }
        U
}

#' Fuzzy PAM on a signature dissimilarity
#'
#' @param d A \code{dist}.
#' @param k Integer \eqn{\ge 2}.
#' @param fuzziness Numeric \eqn{> 1}. Default 2.
#' @param Dm Optional full matrix of \code{d}, reused across candidates
#'   instead of rebuilding it per call.
#' @return List with \code{memberships}, \code{pam_result}, \code{d_to_medoids},
#'   \code{medoid_dist}, \code{fuzziness}.
#' @importFrom cluster pam
#' @importFrom stats as.dist
#' @keywords internal
#' @noRd
.pam_memberships <- function(d, k, fuzziness = 2, Dm = NULL) {
        if (!inherits(d, "dist")) d <- stats::as.dist(d)
        if (length(k) != 1L || is.na(k) || k < 2 || k != as.integer(k))
                stop("`k` must be a single integer >= 2.", call. = FALSE)
        if (length(fuzziness) != 1L || is.na(fuzziness) || fuzziness <= 1)
                stop("`fuzziness` must be a single number > 1.", call. = FALSE)
        
        pam_res <- cluster::pam(d, k = as.integer(k), diss = TRUE)
        if (is.null(Dm)) Dm <- as.matrix(d)
        d_med <- Dm[, pam_res$id.med, drop = FALSE]
        list(memberships  = .memberships_from_medoids(d_med, fuzziness),
             pam_result   = pam_res,
             d_to_medoids = d_med,
             medoid_dist  = Dm[pam_res$id.med, pam_res$id.med, drop = FALSE],
             fuzziness    = fuzziness)
}


# -----------------------------------------------------------------------------*
# 7. Clustering
# -----------------------------------------------------------------------------*

#' Mosaic types from co-occurrence signatures
#'
#' Clusters polygons by their signatures and selects the number of mosaic types.
#'
#' @param signatures Signature matrix from \code{\link{compute_all_signatures}}. With
#'   \code{method = "dmm"} it must carry the \code{"totals"} attribute, which
#'   subsetting or \code{rbind()} drops.
#' @param n_mosaic_types Candidate counts, all \eqn{\ge 2}. Default 2:10.
#' @param method \code{"dmm"} (default): Dirichlet-multinomial mixture on the
#'   signatures converted back to edge pseudo-counts (Holmes et al. 2012).
#'   Memberships are posterior type probabilities ("which type, and how
#'   confidently"). \code{"pam"}: fuzzy PAM on a signature distance;
#'   memberships are distance-based fuzzy degrees that depend on
#'   \code{fuzziness}.
#' @param distance Metric; see \code{\link{compute_signature_distance}}. With
#'   \code{method = "dmm"} it is used only for saturation (estimated on at most
#'   2000 rows) and, when the metric needs them, silhouettes. The full
#'   \eqn{N \times N} matrix is built only in that last case.
#' @param crisp Logical. Hard memberships. With \code{"dmm"} only the returned
#'   memberships are hardened; selection still uses the soft fit.
#' @param metric Selection criterion: \code{"borda"}, \code{"CH_JS"},
#'   \code{"R2_JS"}, \code{"GD5"}, \code{"FS"}, \code{"stability"},
#'   \code{"silhouette_fuzzy"}, \code{"silhouette"}, \code{"ICL"} or
#'   \code{"BIC"}. Default \code{"borda"} for PAM, \code{"ICL"} for DMM.
#'   \code{"GD5"}, \code{"FS"} and \code{"stability"} need medoids (PAM only);
#'   \code{"ICL"} and \code{"BIC"} need a likelihood (DMM only).
#' @param fuzziness Numeric \eqn{> 1}. Default 2. PAM only.
#' @param stability_B PAM refits on 80\% subsamples of the same distance
#'   matrix for \code{STAB}. 0 skips. PAM only. \code{STAB} measures how
#'   sensitive PAM is to the sample, not whether new data would give the same
#'   types. Default 10.
#' @param sample_size Objects clustered exactly. Above this the routine samples,
#'   fits on the sample and assigns the remainder (PAM: nearest medoid; DMM:
#'   posterior under the fitted mixture). \code{Inf} forces the exact path.
#'   Default 5000.
#' @param dmm_control Named list overriding DMM defaults: \code{n_start} (5),
#'   \code{max_iter} (500), \code{tol} (1e-5, relative log-likelihood change, as in mclust), \code{A_init} (50),
#'   \code{alpha_min} (1e-4), \code{A_max} (1e6, precision cap: the multinomial
#'   limit), \code{n_scale} (1). \code{n_scale} multiplies the
#'   pseudo-counts; rerun at 0.5 and 2 to check that \code{k_sel} is not an
#'   artefact of the count scale.
#' @param seed Optional integer; the caller's RNG state is restored on exit.
#' @param verbose Logical.
#' @return A list of class \code{"pulse_mosaic_types"} with \code{clusters} (PAM fit
#'   or \code{pulse_dmm} object), \code{memberships}, \code{hard_labels},
#'   \code{d_to_medoids}, \code{medoid_signatures} (DMM: component mean
#'   signatures), \code{k_sel}, \code{validity}, \code{saturation},
#'   \code{distance}, \code{method}, \code{dmm_control}, \code{sampled} and
#'   \code{valid} (map back to rows of \code{signatures}).
#' @importFrom cluster pam
#' @importFrom stats complete.cases
#' @export
get_mosaic_types <- function(signatures,
                             n_mosaic_types  = 2:10,
                             method      = c("dmm","pam"),
                             distance    = c("js", "hellinger", "tv", "bray", "euclidean"),
                             crisp       = FALSE,
                             metric      = c("borda", "CH_JS", "R2_JS", "GD5", "FS",
                                             "stability", "silhouette_fuzzy", "silhouette",
                                             "ICL", "BIC"),
                             fuzziness   = 2,
                             stability_B = 10L,
                             sample_size = 5000,
                             dmm_control = list(),
                             seed        = NULL,
                             verbose     = TRUE) {
        method   <- match.arg(method)
        distance <- match.arg(distance)
        metric   <- if (missing(metric) && method == "dmm") "ICL" else match.arg(metric)
        
        if (method == "dmm" && metric %in% c("GD5", "FS", "stability"))
                stop(sprintf("metric = '%s' needs PAM medoids. With method = 'dmm' use ",
                             metric), "'ICL', 'BIC', 'CH_JS', 'R2_JS', 'silhouette_fuzzy', ",
                     "'silhouette' or 'borda'.", call. = FALSE)
        if (method == "pam" && metric %in% c("ICL", "BIC"))
                stop("metric = '", metric, "' needs a likelihood; use method = 'dmm'.",
                     call. = FALSE)
        
        ctl <- list(n_start = 5L, max_iter = 500L, tol = 1e-5,
                    A_init = 50, alpha_min = 1e-4, A_max = 1e6, n_scale = "auto")
        bad <- setdiff(names(dmm_control), names(ctl))
        if (length(bad))
                stop("Unknown `dmm_control` entries: ", paste(bad, collapse = ", "),
                     ". Allowed: ", paste(names(ctl), collapse = ", "), ".", call. = FALSE)
        ctl[names(dmm_control)] <- dmm_control
        
        totals <- attr(signatures, "totals")
        if (method == "dmm" && is.null(totals))
                stop("`signatures` has no \"totals\" attribute, which method = 'dmm' needs. ",
                     "Recompute with compute_all_signatures(); subsetting or rbind() drops it.",
                     call. = FALSE)
        if (method == "dmm") ctl$n_scale <- .resolve_n_scale(ctl$n_scale, signatures)
        signatures  <- as.matrix(signatures)
        valid <- stats::complete.cases(signatures)
        if (sum(valid) < 2L) stop("Fewer than two complete signature rows.", call. = FALSE)
        S <- signatures[valid, , drop = FALSE]
        N <- nrow(S)
        if (method == "dmm") tot <- totals[valid]
        
        n_mosaic_types <- sort(unique(as.integer(n_mosaic_types)))
        if (any(is.na(n_mosaic_types)) || any(n_mosaic_types < 2L))
                stop("`n_mosaic_types` must all be integers >= 2.", call. = FALSE)
        
        .local_seed(seed)
        
        sampled <- is.finite(sample_size) && N > sample_size
        idx <- if (sampled) sort(sample.int(N, sample_size)) else seq_len(N)
        if (max(n_mosaic_types) >= length(idx))
                stop(sprintf("`n_mosaic_types` must be < the number of clustered objects (%d).",
                             length(idx)), call. = FALSE)
        
        ## DMM needs no distances to fit. Build the N x N matrix only for the
        ## silhouette-based metrics; estimate saturation on <= 2000 rows.
        need_D <- method == "pam" || metric %in% c("borda", "silhouette", "silhouette_fuzzy")
        D_full <- NULL
        if (need_D) {
                if (verbose) message(sprintf("%s distances for %d of %d signatures ...",
                                             distance, length(idx), N))
                d      <- compute_signature_distance(S[idx, , drop = FALSE], distance)
                D_full <- as.matrix(d)
                sat    <- check_signature_saturation(d, distance)
        } else {
                isat <- if (length(idx) > 2000L) sort(sample(idx, 2000L)) else idx
                sat  <- check_signature_saturation(compute_signature_distance(S[isat, , drop = FALSE], distance),
                                             distance)
        }
        if (is.finite(sat) && sat > 0.3)
                warning(sprintf(paste0("%.0f%% of signature pairs sit at the %s ceiling. ",
                                       "Signatures barely overlap, so any partition of them ",
                                       "is close to arbitrary. Smooth with smooth_signatures(), ",
                                       "aggregate to region level, or raise k."),
                                100 * sat, distance), call. = FALSE)
        
        if (crisp && method == "pam" && metric != "silhouette") {
                if (verbose) message("Fuzzy metrics are degenerate with crisp = TRUE; using 'silhouette'.")
                metric <- "silhouette"
        }
        
        if (method == "dmm") {
                X <- .sig_to_counts(S[idx, , drop = FALSE], tot[idx], ctl$n_scale)
                n_u <- nrow(unique(sqrt(X / rowSums(X))))   # same check kmeans() uses
                bad <- n_mosaic_types > n_u
                if (all(bad)) stop("Only ", n_u, " distinct signature compositions; ",
                                   "every requested n_mosaic_types exceeds this.", call. = FALSE)
                if (any(bad)) {
                        warning("Dropping n_mosaic_types > ", n_u,
                                " (number of distinct signature compositions): ",
                                paste(n_mosaic_types[bad], collapse = ", "), call. = FALSE)
                        n_mosaic_types <- n_mosaic_types[!bad]
                }
        }
        
        if (verbose) message("Evaluating ", method, " solutions (metric = ", metric, ") ...")
        rows <- vector("list", length(n_mosaic_types))
        fits <- vector("list", length(n_mosaic_types))
        for (i in seq_along(n_mosaic_types)) {
                nc <- n_mosaic_types[i]
                if (method == "pam") {
                        fz <- .pam_memberships(d, nc, fuzziness, Dm = D_full)
                        if (crisp) {
                                U <- matrix(0, nrow(fz$memberships), nc)
                                U[cbind(seq_len(nrow(U)), fz$pam_result$clustering)] <- 1
                                fz$memberships <- U
                        }
                        fits[[i]] <- fz
                        rows[[i]] <- .mosaic_cvis(fz, signatures = S[idx, , drop = FALSE], D_full = D_full,
                                                  stability_B = stability_B, seed = seed)
                        if (verbose)
                                message(sprintf("  k=%2d | R2=%.3f CH=%8.1f GD5=%.3f STAB=%.3f SILF=%.3f",
                                                nc, rows[[i]]$R2_JS, rows[[i]]$CH_JS, rows[[i]]$GD5,
                                                rows[[i]]$STAB, rows[[i]]$SILH_fuzzy))
                } else {
                        fit <- .fit_dmm(X, nc, n_start = ctl$n_start, max_iter = ctl$max_iter,
                                        tol = ctl$tol, A_init = ctl$A_init,
                                        alpha_min = ctl$alpha_min, A_max = ctl$A_max)
                        if (!fit$converged)
                                warning(sprintf("DMM not converged at k = %d after %d iterations; ",
                                                nc, ctl$max_iter),
                                        "raise dmm_control$max_iter.", call. = FALSE)
                        fits[[i]] <- fit
                        rows[[i]] <- .mosaic_cvis_dmm(fit, signatures = S[idx, , drop = FALSE],
                                                      D_full = D_full)
                        if (verbose)
                                message(sprintf("  k=%2d | ICL=%.1f BIC=%.1f R2=%.3f empty=%d%s",
                                                nc, rows[[i]]$ICL, rows[[i]]$BIC, rows[[i]]$R2_JS,
                                                rows[[i]]$n_empty,
                                                if (fit$limit == "multinomial") " (multinomial limit)" else ""))
                }
        }
        validity <- do.call(rbind, rows)
        validity$SAT <- sat
        validity$score <- .mosaic_score(validity, metric)
        if (all(is.na(validity$score)))
                stop("Metric '", metric, "' produced only NA scores.", call. = FALSE)
        
        best_i <- which.max(validity$score)
        k_sel <- validity$k[best_i]
        if (verbose) message("Selected ", k_sel, " mosaic types (", metric,
                             " = ", round(validity$score[best_i], 4), ")")
        if (method == "dmm" && metric %in% c("ICL", "BIC") &&
            length(n_mosaic_types) > 1L && k_sel == max(n_mosaic_types))
                warning("Selected k is the largest candidate. Overlapping neighbourhoods count ",
                        "each edge many times, which inflates the DMM likelihood and pushes ",
                        metric, " upward. Extend `n_mosaic_types`, and rerun with ",
                        sprintf("dmm_control = list(n_scale = %.3g) to see whether the selected count moves.",
                                ctl$n_scale / 2),
                        call. = FALSE)
        
        final <- fits[[best_i]]
        if (method == "pam") {
                med_sig <- S[idx, , drop = FALSE][final$pam_result$id.med, , drop = FALSE]
                Dm <- if (sampled) {
                        if (verbose) message("Assigning the remaining ", N - length(idx), " polygons ...")
                        compute_signature_distance(S, distance, signatures_ref = med_sig)
                } else final$d_to_medoids
                U    <- .memberships_from_medoids(Dm, fuzziness)
                hard <- max.col(-Dm, ties.method = "first")
        } else {
                med_sig <- final$alpha / final$concentration     # E[p | type k]
                if (sampled) {
                        if (verbose) message("Assigning the remaining ", N - length(idx), " polygons ...")
                        U <- .score_dmm(final, .sig_to_counts(S, tot, ctl$n_scale))$U
                } else U <- final$U
                hard <- max.col(U, ties.method = "first")
                Dm   <- compute_signature_distance(S, distance, signatures_ref = med_sig)
        }
        if (crisp) {
                U <- matrix(0, nrow(U), ncol(U))
                U[cbind(seq_len(nrow(U)), hard)] <- 1
        }
        
        structure(list(clusters      = final,
                       memberships   = U,
                       hard_labels = hard,
                       d_to_medoids  = Dm,
                       medoid_signatures = med_sig,
                       n_mosaic_types = k_sel,
                       validity      = validity,
                       saturation    = sat,
                       distance      = distance,
                       method        = method,
                       dmm_control   = if (method == "dmm") ctl else NULL,
                       sampled       = sampled,
                       valid         = valid),
                  class = "pulse_mosaic_types")
}

#' @rdname get_mosaic_types
#' @param x An object of class \code{"pulse_mosaic_types"}.
#' @param ... Unused.
#' @return \code{print.pulse_mosaic_types} returns \code{x} invisibly.
#' @export
print.pulse_mosaic_types <- function(x, ...) {
        cat(sprintf("<pulse_mosaic_types>  %d units x %d mosaic types | %s | distance = %s%s\n",
                    nrow(x$memberships), x$n_mosaic_types, if (is.null(x$method)) "pam" else x$method,
                    x$distance, if (x$sampled) " (sampled)" else ""))
        cat(sprintf("  saturation = %.2f | mean certainty = %.3f\n",
                    x$saturation,
                    mean(x$memberships[cbind(seq_len(nrow(x$memberships)), x$hard_labels)])))
        invisible(x)
}


# -----------------------------------------------------------------------------*
# 8. Validity indices
# -----------------------------------------------------------------------------*

#' @keywords internal
#' @noRd
.ent_rows <- function(P) { P <- as.matrix(P); -rowSums(P * log(P + (P == 0))) }

#' @keywords internal
#' @noRd
.jsd_to_ref <- function(P, q) {
        Mx <- sweep(P, 2, q, "+") / 2
        .ent_rows(Mx) - 0.5 * .ent_rows(P) - 0.5 * .ent_rows(matrix(q, nrow = 1L))
}

#' Exact Jensen-Shannon decomposition of the signature mixture
#'
#' Shared by the PAM and DMM validity rows. Types with zero membership mass
#' (possible under DMM) are dropped before the decomposition.
#' @return List with \code{R2_JS} and \code{CH_JS}.
#' @keywords internal
#' @noRd
.js_decomp <- function(U, signatures) {
        N <- nrow(U); n_i <- colSums(U); keep <- n_i > 0
        U <- U[, keep, drop = FALSE]; n_i <- n_i[keep]; k <- ncol(U)
        Sbar <- crossprod(U, signatures) / n_i
        sbar <- colMeans(signatures)
        P_i  <- n_i / N
        Hc <- .ent_rows(Sbar); Hg <- .ent_rows(matrix(sbar, nrow = 1L))
        Ho <- mean(.ent_rows(signatures))
        Btw <- Hg - sum(P_i * Hc); Wth <- sum(P_i * Hc) - Ho
        if (Btw < -1e-8 || Wth < -1e-8) {
                warning(sprintf("Negative JS component at k = %d (Btw = %.3g, Wth = %.3g); ",
                                k, Btw, Wth), "returning NA for R2_JS/CH_JS.", call. = FALSE)
                return(list(R2_JS = NA_real_, CH_JS = NA_real_))
        }
        Btw <- max(Btw, 0); Wth <- max(Wth, 0)
        list(R2_JS = if (Btw + Wth > 0) Btw / (Btw + Wth) else NA_real_,
             CH_JS = if (Wth > 0 && k > 1L) (Btw / (k - 1)) / (Wth / (N - k)) else NA_real_)
}

#' Validity row for one DMM solution
#'
#' Same columns as \code{.mosaic_cvis()} so the two engines rbind and feed
#' \code{check_k_order()} unchanged. Medoid-based indices are \code{NA}, and
#' silhouettes are \code{NA} when \code{D_full} is \code{NULL}.
#' @keywords internal
#' @noRd
.mosaic_cvis_dmm <- function(fit, signatures, D_full = NULL) {
        U <- fit$U; k <- ncol(U)
        if (nrow(signatures) != nrow(U))
                stop(sprintf("`signatures` has %d rows but the fit has %d.", nrow(signatures), nrow(U)),
                     call. = FALSE)
        js  <- .js_decomp(U, signatures)
        lab <- max.col(U, ties.method = "first")
        has_D <- !is.null(D_full)
        data.frame(k = k, R2_JS = js$R2_JS, CH_JS = js$CH_JS, FS = NA_real_, GD5 = NA_real_,
                   STAB = NA_real_,
                   SILH_fuzzy = if (has_D) .fuzzy_silhouette(D = D_full, U = U) else NA_real_,
                   SILH_HARD  = if (has_D) mean(.rousseeuw_silhouette(D_full, lab)) else NA_real_,
                   loglik = fit$loglik, BIC = fit$BIC, ICL = fit$ICL,
                   n_empty = sum(tabulate(lab, k) == 0L),
                   multinomial_limit = identical(fit$limit, "multinomial"))
}

#' Validity indices for one mosaic solution
#'
#' @param fz Output of \code{.pam_memberships()}.
#' @param signatures Signature matrix on the SAME rows as \code{fz$memberships}.
#' @param D_full Full dissimilarity matrix PAM was run on.
#' @param stability_B PAM refits on 80\% subsamples.
#' @param seed Optional seed.
#' @return One-row data frame.
#' @keywords internal
#' @noRd
.mosaic_cvis <- function(fz, signatures, D_full, stability_B = 10L, seed = NULL) {
        U <- fz$memberships; N <- nrow(U); k <- ncol(U); m <- fz$fuzziness
        med <- fz$pam_result$id.med
        if (nrow(signatures) != N)
                stop(sprintf("`signatures` has %d rows but the partition has %d; pass the ",
                             "complete-case subset.", nrow(signatures), N), call. = FALSE)
        
        D  <- fz$d_to_medoids
        n_i <- colSums(U); s_w <- colSums(U^m)
        J_m <- sum(U^m * D^2)
        
        js <- .js_decomp(U, signatures)
        R2_JS <- js$R2_JS; CH_JS <- js$CH_JS
        sbar  <- colMeans(signatures)
        
        b_i <- .jsd_to_ref(signatures[med, , drop = FALSE], sbar)
        FS  <- J_m - sum(s_w * b_i)
        
        M <- crossprod(U, D)
        Delta3 <- 2 * diag(M) / n_i
        S5 <- (M + t(M)) / outer(n_i, n_i, "+")
        diag(S5) <- Inf
        GD5 <- if (max(Delta3) > 0) min(S5) / max(Delta3) else NA_real_
        
        STAB <- .stability_pam_diss(D_full, fz$pam_result$clustering, k,
                                    B = stability_B, seed = seed)
        SILF <- .fuzzy_silhouette(D = D_full, U = U)
        
        data.frame(k = k, R2_JS = R2_JS, CH_JS = CH_JS, FS = FS, GD5 = GD5,
                   STAB = STAB, SILH_fuzzy = SILF,
                   SILH_HARD = fz$pam_result$silinfo$avg.width)
}

#' @keywords internal
#' @noRd
.stability_pam_diss <- function(D_full, hard_labels, k, B = 10,
                                subsample_frac = 0.8, seed = NULL) {
        if (B <= 0L) return(NA_real_)
        N <- nrow(D_full); n_sub <- floor(N * subsample_frac)
        if (n_sub < k * 2L) return(NA_real_)
        .local_seed(seed)
        aris <- vapply(seq_len(B), function(b) {
                idx <- sort(sample.int(N, n_sub, replace = FALSE))
                lab <- cluster::pam(stats::as.dist(D_full[idx, idx]), k = k, diss = TRUE)$clustering
                .adjusted_rand_index(lab, hard_labels[idx])
        }, numeric(1))
        mean(aris, na.rm = TRUE)
}

#' @keywords internal
#' @noRd
.fuzzy_silhouette <- function(D, U, alpha = 1) {
        N <- nrow(U); k <- ncol(U)
        D_size <- if (inherits(D, "dist")) attr(D, "Size") else nrow(D)
        if (N != D_size) stop("Dimensions of D and U disagree.", call. = FALSE)
        if (k < 2L) return(NA_real_)
        labels <- max.col(U, ties.method = "first")
        s_i <- .rousseeuw_silhouette(D, labels)
        sorted <- apply(U, 1L, sort, decreasing = TRUE)
        if (!is.matrix(sorted)) return(NA_real_)
        w <- (sorted[1L, ] - sorted[2L, ]) ^ alpha
        if (sum(w) <= 0) return(NA_real_)
        sum(w * s_i) / sum(w)
}

#' @keywords internal
#' @noRd
.mosaic_score <- function(df, metric) {
        switch(metric,
               R2_JS            = df$R2_JS,
               CH_JS            = df$CH_JS,
               FS               = -df$FS,
               GD5              = df$GD5,
               stability        = df$STAB,
               silhouette       = df$SILH_HARD,
               silhouette_fuzzy = df$SILH_fuzzy,
               ICL              = -df$ICL,
               BIC              = -df$BIC,
               borda            = {
                       cols <- c("CH_JS", "STAB", "SILH_fuzzy")
                       cols <- cols[vapply(cols, function(cc) !all(is.na(df[[cc]])), logical(1))]
                       if (!length(cols)) return(rep(NA_real_, nrow(df)))
                       -.borda_score(df[, cols, drop = FALSE], rep("max", length(cols)))
               },
               stop("Unknown metric: ", metric, call. = FALSE))
}


# -----------------------------------------------------------------------------*
# 9. Sweep and diagnostics
# -----------------------------------------------------------------------------*

#' Sweep neighbourhood order and mosaic-type count
#'
#' @param mosaic_input Membership matrix or \code{\link{as_mosaic_input}} output.
#' @param graph An \code{igraph} adjacency graph.
#' @param coords Centroid matrix.
#' @param k_order Neighbourhood orders. Default 1:5.
#' @param n_mosaic_types Candidate mosaic-type counts. Default 2:10.
#' @param kernel,sigma,exclude_bridges Passed to
#'   \code{\link{compute_all_signatures}}. \code{exclude_bridges} defaults to
#'   \code{TRUE}.
#' @param verbose Logical. Report progress and the per-order summary.
#'   Passed on to \code{\link{compute_all_signatures}} and
#'   \code{\link{get_mosaic_types}}. Default \code{TRUE}.
#' @param ... Passed to \code{\link{get_mosaic_types}} (\code{distance},
#'   \code{metric}, \code{stability_B}, \code{sample_size}, \code{seed}).
#' @details \code{k} is bounded on both sides: too small and the signatures are
#'   noisy and saturate; too large and every neighbourhood converges on the
#'   global mixture, so the signature field goes flat and the clustering returns
#'   a few huge blobs. The returned \code{SAT} and \code{contrast} columns
#'   diagnose the two ends.
#'
#'   \code{score} is a Borda rank among the \code{n_mosaic_types} candidates of one
#'   \code{k_order}, so it selects \code{n_mosaic_types} within an order and must not be
#'   compared across orders. The function therefore does not pick an order.
#'   To compare orders, use the bounded columns at a fixed type count, e.g.
#'   \code{subset(res, n_mosaic_types == 6)[, c("k_order", "R2_JS", "SAT", "contrast")]}.
#'   \code{R2_JS} rises with \code{n_mosaic_types} and may rise with \code{k_order}
#'   as signatures smooth, so read it together with \code{SAT} and
#'   \code{contrast}, not alone.
#' @return A data frame with one row per (\code{k_order}, \code{n_mosaic_types}),
#'   with attribute \code{"by_order"}: one row per order at its own
#'   \code{n_mosaic_types}.
#' @export
check_k_order <- function(mosaic_input, graph, coords = NULL, k_order = 1:5, n_mosaic_types = 2:10,
                          kernel = "none", sigma = "global", exclude_bridges = TRUE,
                          verbose = TRUE, ...) {
        if (!igraph::is_igraph(graph)) stop("`graph` must be an igraph object.", call. = FALSE)
        .assert_flag(verbose, "verbose")
        k_order <- sort(unique(as.integer(k_order)))
        
        out <- vector("list", length(k_order))
        for (i in seq_along(k_order)) {
                k <- k_order[i]
                if (verbose) message("\n=== Neighbourhood order k = ", k, " ===")
                signatures <- compute_all_signatures(mosaic_input, graph, k, coords = coords,
                                               kernel = kernel, sigma = sigma,
                                               exclude_bridges = exclude_bridges,
                                               verbose = verbose)
                cl <- get_mosaic_types(signatures, n_mosaic_types = n_mosaic_types,
                                       verbose = verbose, ...)
                out[[i]] <- data.frame(k_order    = k,
                                       n_mosaic_types = cl$validity$k,
                                       R2_JS      = cl$validity$R2_JS,
                                       CH_JS      = cl$validity$CH_JS,
                                       STAB       = cl$validity$STAB,
                                       SILH_fuzzy = cl$validity$SILH_fuzzy,
                                       SAT        = cl$validity$SAT,
                                       score      = cl$validity$score,
                                       n_selected = cl$n_mosaic_types,
                                       contrast   = mean(compute_mosaic_contrast(signatures, graph),
                                                         na.rm = TRUE))
        }
        results <- do.call(rbind, out)
        
        # `score` is a Borda rank within one k_order and is not comparable across
        # orders. Summarise each order at its own n_mosaic_types; no cross-order pick.
        by_order <- results[results$n_mosaic_types == results$n_selected,
                            c("k_order", "n_selected", "R2_JS", "SAT", "contrast"), drop = FALSE]
        rownames(by_order) <- NULL
        if (verbose) message("\n=== Per order at its own n_mosaic_types (score is within-order only; ",
                "orders are not ranked) ===\n",
                paste(sprintf("  k_order %d: n_selected = %d, R2_JS = %.3f, SAT = %.2f, contrast = %.3f",
                              by_order$k_order, by_order$n_selected, by_order$R2_JS,
                              by_order$SAT, by_order$contrast), collapse = "\n"))
        hi <- by_order$k_order[is.finite(by_order$SAT) & by_order$SAT > 0.3]
        if (verbose && length(hi))
                message("Saturation > 0.3 at k_order = ", paste(hi, collapse = ", "),
                        ": signatures there carry little information.")
        attr(results, "by_order") <- by_order
        results
}

#' Signature contrast against the immediate neighbourhood
#'
#' Mean absolute difference between each polygon's signature and the mean of its
#' graph neighbours'. Near-zero everywhere means the large contiguous mosaics
#' are a real property of the landscape; sizeable contrasts that nevertheless
#' cluster into big blobs mean the overlapping neighbourhoods have smoothed real
#' structure away, and the kernel bandwidth should be tightened or the
#' neighbourhood order reduced.
#'
#' @param signatures Signature matrix.
#' @param graph The adjacency graph.
#' @return Numeric vector of length \eqn{n}.
#' @importFrom igraph as_adjacency_matrix
#' @export
compute_mosaic_contrast <- function(signatures, graph) {
        ## unclass(): Matrix's S4 %*% recurses without end on the S3
        ## "pulse_signatures" class (C stack overflow in check_k_order()).
        S <- unclass(as.matrix(signatures))
        A <- igraph::as_adjacency_matrix(graph, sparse = TRUE, attr = NULL)
        ok <- stats::complete.cases(S)
        S0 <- S; S0[!ok, ] <- 0
        deg <- as.numeric(A %*% ok)
        nb  <- as.matrix(A %*% S0) / pmax(deg, 1)
        out <- rowMeans(abs(S - nb))
        out[!ok | deg == 0] <- NA_real_
        out
}


#' Validate and package a mosaic input matrix
#' @keywords internal
#' @noRd
.new_mosaic_input <- function(P, labels, region_of, level, normalise, drop_empty) {
        P <- as.matrix(P)
        if (!is.numeric(P))  stop("Membership matrix must be numeric.", call. = FALSE)
        if (anyNA(P))        stop("Membership matrix contains NA.", call. = FALSE)
        if (any(P < 0))      stop("Membership matrix has negative entries.", call. = FALSE)
        if (nrow(P) < 2L)    stop("Membership matrix must have at least two rows.", call. = FALSE)
        
        if (drop_empty) {
                keep <- colSums(P) > 0
                if (any(!keep)) {
                        warning(sprintf("Dropping %d class(es) with no mass: %s.",
                                        sum(!keep), paste(labels[!keep], collapse = ", ")),
                                call. = FALSE)
                        P <- P[, keep, drop = FALSE]
                        labels <- labels[keep]
                        region_of <- region_of[keep]
                }
        }
        if (ncol(P) < 2L)
                stop("Membership matrix must have at least two non-empty class columns; ",
                     "a co-occurrence signature needs at least two classes to distinguish.",
                     call. = FALSE)
        
        rs <- rowSums(P)
        n_empty <- sum(rs <= 0)
        if (n_empty > 0L)
                warning(sprintf("%d row(s) have zero total membership; those polygons ",
                                n_empty), "contribute nothing and will get NA signatures.",
                        call. = FALSE)
        if (normalise) {
                P <- P / pmax(rs, .Machine$double.eps)
                rs <- rowSums(P)
        }
        
        if (length(labels) != ncol(P) || length(region_of) != ncol(P))
                stop("Internal: labels / region_of do not match the column count.",
                     call. = FALSE)
        if (anyDuplicated(labels))
                stop("Class labels must be unique; they become signature column names.",
                     call. = FALSE)
        colnames(P) <- labels
        
        structure(list(P = P, labels = labels,
                       region_of = as.integer(region_of), level = level,
                       row_sums = range(rs)),
                  class = "pulse_mosaic_input")
}

#' @rdname as_mosaic_input
#' @param ... Unused (print method).
#' @return \code{print.pulse_mosaic_input} returns \code{x} invisibly.
#' @export
print.pulse_mosaic_input <- function(x, ...) {
        cat(sprintf("<pulse_mosaic_input>  %d units x %d classes | level = %s\n",
                    nrow(x$P), ncol(x$P), x$level))
        cat(sprintf("  signature bins: %d | row sums in [%.3f, %.3f]\n",
                    ncol(x$P) * (ncol(x$P) + 1L) / 2L, x$row_sums[1], x$row_sums[2]))
        if (all(is.na(x$region_of))) {
                cat("  classes have no parent region (they cross regions, or none is known)\n")
        } else {
                nr <- length(unique(x$region_of))
                cat(sprintf("  classes nested in %d region(s)%s\n", nr,
                            if (ncol(x$P) == nr) "" else
                                    " -- class identity determines region, so mosaic types will track the region map"))
        }
        if (abs(diff(x$row_sums)) > 1e-6)
                cat("  note: rows do not all sum to the same total, so polygons contribute\n",
                    "  unequal mass to their neighbours' signatures (see `normalise`)\n")
        invisible(x)
}
