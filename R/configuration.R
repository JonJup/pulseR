# =============================================================================
# pulseR — LEVEL 3: label-invariant configuration features
#
# WHY        The co-occurrence signature is label-dependent: a unit inside a
#            type-1 patch looks "pure 1:1", inside a type-3 patch "pure 3:3".
#            Zones that differ only in arrangement (patchy vs mixed, same types
#            in the same shares) are then unions of several signature profiles,
#            and a mixture with one profile per component cannot recover them.
# WHAT       For each unit and graph lag d, the share of unit pairs at distance d
#            that carry the same type, compared with the share expected if the
#            same units were arranged at random (composition-adjusted, kappa-
#            style). Renaming the types leaves every feature unchanged.
# USE        Pass the result as `configuration` to get_mosaic_types(). The
#            signatures then carry composition (label-specific) and this matrix
#            carries arrangement (label-invariant).
# =============================================================================


#' Label-invariant configuration features for every polygon
#'
#' For each unit and graph lag \eqn{d}, measures how often units at distance
#' \eqn{d} share a type, relative to what the window's composition alone would
#' give. The result does not depend on the type labels, so a patch of type 1
#' and a patch of type 3 give the same value. Pass it to
#' \code{\link{get_mosaic_types}} as \code{configuration} to separate zones that
#' differ in spatial arrangement (patchy vs mixed) as well as in composition.
#'
#' @details
#' Memberships are first scaled to unit length, so the similarity of two units
#' is the cosine \eqn{\tilde p_u \cdot \tilde p_v}: 1 for two identical crisp
#' units, 0 for two different ones. Without this, diffuse memberships would
#' look like fine-grained mixing.
#'
#' For lag \eqn{d}, unit \eqn{u} contributes the summed similarity to all units
#' exactly \eqn{d} steps away, and the number of such units. Both are summed
#' over the focal's neighbourhood with the same weights as
#' \code{\link{compute_all_signatures}}, giving the observed same-type share
#' \eqn{q_d}. The share expected under random arrangement of the window's units
#' is \eqn{e = \lVert \bar p \rVert^2}, with \eqn{\bar p} the weighted mean of
#' the scaled memberships (for crisp input, \eqn{\sum_c p_c^2}) over a wider
#' reference window of order \code{k_reference}. The wider window matters: deep
#' inside a patch the focal window holds one type, so its own composition would
#' predict "all pairs alike" and hide the clumping. Against the reference
#' composition, a patch interior scores high and a random mix scores near 0.
#' The feature is
#' \deqn{\kappa_d = (q_d - e) / \max(1 - e, \code{kappa\_floor}).}
#' \eqn{\kappa_d = 0}: arrangement as random; \eqn{\kappa_d > 0}: clumped;
#' \eqn{\kappa_d < 0}: alternating. The floor keeps windows dominated by one
#' type, where arrangement is barely defined, from producing extreme values.
#'
#' \strong{Anisotropy.} With \code{anisotropy_lags}, the function also
#' measures whether same-type pairs line up along one direction, as in stripes
#' or bands. For each pair \eqn{(u, v)} at lag \eqn{d} it takes the bearing
#' \eqn{\theta_{uv}} between the two centroids and the doubled-angle vector
#' \eqn{(\cos 2\theta, \sin 2\theta)}, so that \eqn{\theta} and
#' \eqn{\theta + \pi} count as the same axis. Over the window,
#' \deqn{Z_d = \sum (s_{uv} - q_d)\, e^{2 i \theta_{uv}}, \qquad
#' a_d = |Z_d| / (n_d \max(1 - e, \code{kappa\_floor})),}
#' with \eqn{s_{uv}} the pair similarity, \eqn{q_d} the window's same-type
#' share at that lag and \eqn{n_d} the number of pairs. Subtracting \eqn{q_d}
#' removes the graph's own directional bias (on a lattice, more pairs in some
#' directions than others), so only direction-dependent similarity remains.
#' \eqn{a_d} is about 0 for isotropic patterns, whether patchy or random, and
#' large for bands. Taking the modulus makes it independent of the band
#' direction (rotation-invariant), and it stays label-invariant. Under random
#' arrangement \eqn{a_d} is small but positive, since it is a modulus. Use
#' lags up to about the band width: \code{1:2} for bands 1–2 units wide.
#'
#' \strong{Zone-scale smoothing.} With \code{smooth_k > 0}, each feature is
#' replaced by its weighted mean over the order-\code{smooth_k} neighbourhood
#' (same weights as above; \code{NA} values are skipped). Patch interiors and
#' patch edges of one zone then share one value, so a mixture component can
#' describe the zone instead of splitting it into "interior" and "edge" units.
#' The cost is blurrier zone boundaries, about \code{smooth_k} steps wide.
#' Choose \code{smooth_k} below half the expected zone width. The features stay
#' label-invariant.
#'
#' Several lags describe grain: \eqn{\kappa_d} stays high across lags in coarse
#' patches and falls quickly in fine ones.
#'
#' \strong{Scale.} Use a neighbourhood (\code{k_order}) at least as wide as the
#' largest lag, and wide enough to contain a few patches. Units deep inside a
#' patch otherwise see only that patch.
#'
#' @param mosaic_input Numeric \eqn{n \times c} membership matrix, or the list
#'   returned by \code{\link{as_mosaic_input}}.
#' @param graph An \code{igraph} adjacency graph, vertices in row order.
#' @param k_order Neighbourhood order of the window the features are summed
#'   over.
#' @param lags Integer graph lags. Default \code{1:3}.
#' @param coords,kernel,sigma Neighbourhood weighting, as in
#'   \code{\link{compute_all_signatures}}.
#' @param k_reference Order of the reference window for the expected share.
#'   Default \code{2 * k_order}. Make it wider than the largest patch. Cost
#'   grows with the window size (about \eqn{2 k^2} units per window on a grid).
#' @param kappa_floor Lower bound on the denominator \eqn{1 - e}. Default 0.05.
#' @param anisotropy_lags Integer lags for the anisotropy features, or
#'   \code{NULL} (default, none). Needs \code{coords}.
#' @param smooth_k Order of the zone-scale smoothing window. Default
#'   \code{0} (no smoothing).
#' @param exclude_bridges Drop edges flagged \code{bridge == TRUE} first.
#'   Default \code{TRUE}.
#' @param verbose Logical.
#' @return Numeric matrix, one row per polygon, with columns
#'   \code{"kappa_lag<d>"} and, with \code{anisotropy_lags},
#'   \code{"aniso_lag<d>"}; class \code{"pulse_configuration"}, and attributes
#'   \code{same_share} (observed \eqn{q_d}), \code{expected} (\eqn{e}),
#'   \code{lags} and \code{k_order}. Rows whose window has no pair at some lag
#'   are \code{NA}.
#' @importFrom igraph as_adjacency_matrix is_igraph vcount
#' @importFrom Matrix Diagonal crossprod drop0 rowSums colSums
#' @export
compute_configuration <- function(mosaic_input, graph, k_order,
                                  lags = 1:3,
                                  k_reference = 2L * k_order,
                                  coords = NULL,
                                  kernel = c("none", "gaussian", "exponential", "linear"),
                                  sigma  = "global",
                                  kappa_floor = 0.05,
                                  anisotropy_lags = NULL,
                                  smooth_k = 0L,
                                  exclude_bridges = TRUE,
                                  verbose = TRUE) {
        kernel <- match.arg(kernel)
        P <- mosaic_input
        if (is.list(P) && !is.null(P$P)) P <- P$P
        P <- as.matrix(P)
        if (!is.numeric(P)) stop("`mosaic_input` must be numeric.", call. = FALSE)
        if (anyNA(P)) stop("`mosaic_input` contains NA.", call. = FALSE)
        if (ncol(P) < 2L) stop("`mosaic_input` must have at least two class columns.", call. = FALSE)
        if (!igraph::is_igraph(graph)) stop("`graph` must be an igraph object.", call. = FALSE)
        n <- nrow(P)
        if (n != igraph::vcount(graph))
                stop(sprintf("`mosaic_input` has %d rows but the graph has %d vertices.",
                             n, igraph::vcount(graph)), call. = FALSE)
        if (!is.numeric(lags) || !length(lags) || anyNA(lags) || any(lags < 1) ||
            any(lags != as.integer(lags)) || anyDuplicated(lags))
                stop("`lags` must be distinct positive integers.", call. = FALSE)
        lags <- sort(as.integer(lags))
        if (length(kappa_floor) != 1L || !is.finite(kappa_floor) || kappa_floor <= 0 || kappa_floor > 1)
                stop("`kappa_floor` must be one number in (0, 1].", call. = FALSE)
        if (!is.logical(exclude_bridges) || length(exclude_bridges) != 1L || is.na(exclude_bridges))
                stop("`exclude_bridges` must be TRUE or FALSE.", call. = FALSE)
        if (length(k_reference) != 1L || is.na(k_reference) || k_reference < k_order ||
            k_reference != as.integer(k_reference))
                stop("`k_reference` must be one integer >= `k_order`.", call. = FALSE)
        if (length(smooth_k) != 1L || is.na(smooth_k) || smooth_k < 0 ||
            smooth_k != as.integer(smooth_k))
                stop("`smooth_k` must be one integer >= 0.", call. = FALSE)
        aniso <- length(anisotropy_lags) > 0L
        if (aniso) {
                al <- anisotropy_lags
                if (!is.numeric(al) || anyNA(al) || any(al < 1) || any(al != as.integer(al)) ||
                    anyDuplicated(al))
                        stop("`anisotropy_lags` must be distinct positive integers.", call. = FALSE)
                anisotropy_lags <- sort(as.integer(al))
                if (is.null(coords))
                        stop("`coords` are required for anisotropy features.", call. = FALSE)
                coords <- as.matrix(coords)
                if (!is.numeric(coords) || nrow(coords) != n || ncol(coords) < 2L ||
                    anyNA(coords[, 1:2]))
                        stop("`coords` must be a numeric matrix with one row per vertex and ",
                             "two columns, no NA.", call. = FALSE)
        }
        if (k_order < max(lags))
                warning("`k_order` (", k_order, ") is below the largest lag (", max(lags),
                        "). Each window then holds few pairs at that lag.", call. = FALSE)

        ## Unit-length rows: similarity = cosine, so fuzziness does not mimic mixing.
        nrm <- sqrt(rowSums(P^2))
        Pn  <- P / ifelse(nrm > 0, nrm, 1)

        if (exclude_bridges) graph <- .drop_bridges(graph, verbose)
        A  <- igraph::as_adjacency_matrix(graph, sparse = TRUE, attr = NULL)
        A  <- methods::as(Matrix::drop0(A), "CsparseMatrix")
        A@x <- rep.int(1, length(A@x))
        AI <- A + Matrix::Diagonal(n)

        d_max <- max(c(lags, anisotropy_lags))
        if (verbose) message(sprintf("Lag shells up to d = %d for %d polygons ...", d_max, n))
        same <- cnt <- matrix(0, n, length(lags))
        dir6 <- if (aniso) array(0, c(n, 6L, length(anisotropy_lags))) else NULL
        R_prev <- Matrix::Diagonal(n)                 # pattern of (A + I)^0
        R_cur  <- R_prev
        for (d in seq_len(d_max)) {
                R_cur <- R_cur %*% AI
                R_cur <- methods::as(Matrix::drop0(R_cur), "CsparseMatrix")
                R_cur@x <- rep.int(1, length(R_cur@x))
                j <- match(d, lags)
                if (!is.na(j)) {
                        Sh <- Matrix::drop0(R_cur - R_prev)   # exactly d steps away
                        same[, j] <- Matrix::rowSums((Sh %*% Pn) * Pn)
                        cnt[, j]  <- Matrix::rowSums(Sh)
                }
                ja <- if (aniso) match(d, anisotropy_lags) else NA_integer_
                if (!is.na(ja)) {
                        Sh <- Matrix::drop0(R_cur - R_prev)
                        dir6[, , ja] <- .direction_sums(Sh, Pn, coords)
                }
                R_prev <- R_cur
        }

        if (verbose) message(sprintf("Window sums (k = %d) ...", k_order))
        W  <- compute_neighbourhood_weights(get_reachability(graph, k_order), coords, kernel, sigma)
        ws <- Matrix::colSums(W)                      # column = focal
        S_w <- as.matrix(Matrix::crossprod(W, same))
        C_w <- as.matrix(Matrix::crossprod(W, cnt))
        Wr <- if (k_reference == k_order) W else
                compute_neighbourhood_weights(get_reachability(graph, k_reference), coords, kernel, sigma)
        wr <- Matrix::colSums(Wr)
        pbar <- as.matrix(Matrix::crossprod(Wr, Pn)) / ifelse(wr > 0, wr, NA_real_)
        e <- rowSums(pbar^2)

        q <- S_w / ifelse(C_w > 0, C_w, NA_real_)
        kap <- (q - e) / pmax(1 - e, kappa_floor)
        colnames(kap) <- paste0("kappa_lag", lags)
        colnames(q)   <- paste0("same_lag", lags)
        kap[!is.finite(kap)] <- NA_real_
        if (aniso) {
                an <- matrix(NA_real_, n, length(anisotropy_lags),
                             dimnames = list(NULL, paste0("aniso_lag", anisotropy_lags)))
                den_e <- pmax(1 - e, kappa_floor)
                for (ja in seq_along(anisotropy_lags)) {
                        D  <- as.matrix(Matrix::crossprod(W, dir6[, , ja]))
                        qd <- D[, 1] / ifelse(D[, 6] > 0, D[, 6], NA_real_)
                        zc <- D[, 2] - qd * D[, 4]
                        zs <- D[, 3] - qd * D[, 5]
                        an[, ja] <- sqrt(zc^2 + zs^2) / (D[, 6] * den_e)
                }
                an[!is.finite(an)] <- NA_real_
                kap <- cbind(kap, an)
        }
        if (smooth_k > 0L) {
                if (verbose) message(sprintf("Zone-scale smoothing (k = %d) ...", smooth_k))
                Ws <- compute_neighbourhood_weights(get_reachability(graph, smooth_k),
                                                    coords, kernel, sigma)
                ok <- !is.na(kap)
                K0 <- kap; K0[!ok] <- 0
                num <- as.matrix(Matrix::crossprod(Ws, K0))
                den <- as.matrix(Matrix::crossprod(Ws, ok * 1))
                kap[] <- num / ifelse(den > 0, den, NA_real_)
        }
        structure(kap, class = c("pulse_configuration", "matrix", "array"),
                  same_share = q, expected = e, lags = lags, k_order = k_order,
                  k_reference = k_reference, smooth_k = as.integer(smooth_k),
                  anisotropy_lags = if (aniso) anisotropy_lags else integer(0),
                  kappa_floor = kappa_floor)
}


## Per-unit sums over the pairs in shell Sh: similarity, similarity x cos 2theta,
## similarity x sin 2theta, cos 2theta, sin 2theta, count. n x 6. Pairs with
## coincident centroids carry no direction and are skipped. Chunked so memory
## stays bounded on large graphs.
.direction_sums <- function(Sh, Pn, coords, chunk = 1e6) {
        n  <- nrow(Pn)
        T  <- methods::as(Sh, "TsparseMatrix")
        i  <- T@i + 1L; j <- T@j + 1L
        out <- matrix(0, n, 6L)
        starts <- seq(1L, length(i), by = chunk)
        for (st in starts) {
                ix <- st:min(length(i), st + chunk - 1L)
                ii <- i[ix]; jj <- j[ix]
                dx <- coords[jj, 1] - coords[ii, 1]
                dy <- coords[jj, 2] - coords[ii, 2]
                ok <- (dx != 0 | dy != 0)
                if (!any(ok)) next
                ii <- ii[ok]; jj <- jj[ok]
                th <- 2 * atan2(dy[ok], dx[ok])
                sim <- rowSums(Pn[ii, , drop = FALSE] * Pn[jj, , drop = FALSE])
                c2 <- cos(th); s2 <- sin(th)
                agg <- rowsum(cbind(sim, sim * c2, sim * s2, c2, s2, 1), ii, reorder = FALSE)
                rr <- as.integer(rownames(agg))
                out[rr, ] <- out[rr, ] + agg
        }
        out
}
