# =============================================================================
# Spatial network construction + package utility helpers
# -----------------------------------------------------------------------------
# Hardened version. Key fixes over the original (see accompanying notes):
#   * Removed the `a:b` count-down footguns: `(off+1):(off+n)` and `1:n` both
#     misbehave when n == 0. Replaced with `seq_len()`-based indexing.
#   * Vectorised neighbour-list -> edge-list construction (no per-polygon
#     data.frame churn); both helpers now return a typed 0-row data.frame
#     instead of NULL when there are no edges.
#   * Dropped the unused `temp_dir` argument; exposed `snap` / `queen` instead
#     of hard-coding them (their units depend on the CRS).
#   * Added input validation, optional geometry repair, and a CRS check.
#   * Fixed the embedded-newline message string in .warn_if_longlat() and the
#     duplicate st_is_longlat() call; .assert_nonneg() now rejects Inf/NaN.
#   * Added roxygen to the previously undocumented functions.
# =============================================================================


#' Extract internal adjacency edges and boundary ("rim") polygons from one tile
#'
#' Reads a single spatial file (one tile of a larger, tiled polygon dataset),
#' assigns globally unique polygon IDs, computes the within-tile Queen/Rook
#' contiguity edges, and returns the subset of polygons touching the tile's
#' bounding box so that cross-tile adjacency can be "stitched" later.
#'
#' @param file_path Path to a polygon file readable by \code{\link[sf]{st_read}}.
#' @param global_offset Integer. IDs for this tile are assigned as
#'   \code{global_offset + seq_len(n)}, so the caller must pass a non-overlapping
#'   offset per tile (e.g. the running cumulative polygon count).
#' @param snap Numeric snapping tolerance passed to \code{\link[spdep]{poly2nb}}.
#'   \strong{Units follow the CRS} (metres for projected data, degrees for
#'   geographic) -- choose accordingly. Default \code{1e-6}.
#' @param queen Logical. \code{TRUE} (default) for Queen contiguity (shared
#'   edge \emph{or} vertex), \code{FALSE} for Rook (shared edge only).
#' @param make_valid Logical. If \code{TRUE} (default), repair geometries with
#'   \code{\link[sf]{st_make_valid}} before computing contiguity.
#' @param verbose Logical. Emit progress / CRS messages. Default \code{FALSE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{edges}{A \code{data.frame} of \code{from}/\code{to} \emph{global}
#'       IDs (upper triangle only; 0 rows if the tile has no internal edges).}
#'     \item{rim}{An \code{sf} object (global_id + geometry) of polygons touching
#'       the tile bounding box.}
#'     \item{bbox}{The tile \code{bbox}.}
#'     \item{count}{Number of polygons in the tile.}
#'     \item{file}{The input path.}
#'   }
#'
#' @details The "rim" set is a \emph{heuristic} for stitching: it keeps polygons
#'   whose geometry intersects the tile's (axis-aligned) bounding-box perimeter.
#'   This is exact for regular grid tiling, where shared tile edges coincide with
#'   bbox edges. For irregular partitions it can miss polygons that lie on a
#'   shared boundary without reaching the rectangular extent; in that case widen
#'   the catch (e.g. buffer the boundary) before stitching.
#'
#' @importFrom sf st_read st_bbox st_as_sfc st_cast st_intersects st_crs st_make_valid
#' @importFrom spdep poly2nb
#' @keywords internal
#' @noRd
process_file_internals <- function(file_path,
                                   global_offset,
                                   snap       = 1e-6,
                                   queen      = TRUE,
                                   make_valid = TRUE,
                                   verbose    = FALSE) {
        
        # ---- Validation ---------------------------------------------------------
        if (!is.character(file_path) || length(file_path) != 1L || !file.exists(file_path))
                stop("`file_path` must point to an existing file.", call. = FALSE)
        .assert_nonneg(global_offset, "global_offset")
        .assert_nonneg(snap, "snap")
        .assert_flag(queen, "queen")
        .assert_flag(make_valid, "make_valid")
        
        # ---- 1. Load ------------------------------------------------------------
        polys <- sf::st_read(file_path, quiet = TRUE)
        n <- nrow(polys)
        
        # Empty tile: return a well-formed, empty result rather than letting the
        # `:` operator or poly2nb() blow up downstream.
        if (n == 0L) {
                return(list(edges = data.frame(from = integer(0), to = integer(0)),
                            rim   = polys[, character(0)],
                            bbox  = sf::st_bbox(polys),
                            count = 0L,
                            file  = file_path))
        }
        
        .warn_if_longlat(sf::st_crs(polys), verbose)
        if (make_valid) polys <- sf::st_make_valid(polys)
        
        # ---- 2. Globally unique IDs --------------------------------------------
        # seq_len(n) is empty-safe; `(global_offset + 1):(global_offset + n)` counts
        # DOWN when n == 0 and is the classic source of silent off-by-one bugs.
        polys$global_id <- global_offset + seq_len(n)
        
        # ---- 3. Internal contiguity --------------------------------------------
        nb <- spdep::poly2nb(polys, queen = queen, snap = snap)
        
        # ---- 4. Vectorised edge extraction (upper triangle, global IDs) --------
        # spdep encodes "no neighbours" as the single value 0L, so lengths(nb) >= 1.
        from_local <- rep.int(seq_len(n), lengths(nb))
        to_local   <- unlist(nb, use.names = FALSE)
        keep       <- to_local != 0L & to_local > from_local   # one direction only
        edge_df <- data.frame(
                from = polys$global_id[from_local[keep]],
                to   = polys$global_id[to_local[keep]]
        )
        
        # ---- 5. Rim polygons (boundary candidates for stitching) ---------------
        bbox          <- sf::st_bbox(polys)
        bbox_boundary <- sf::st_cast(sf::st_as_sfc(bbox), "MULTILINESTRING")
        rim_hit       <- sf::st_intersects(polys, bbox_boundary, sparse = FALSE)[, 1]
        rim_polys     <- polys[rim_hit, "global_id"]            # sf keeps geometry (sticky)
        
        .log(verbose, "%s: %d polygons, %d internal edges, %d rim polygons.",
             basename(file_path), n, nrow(edge_df), nrow(rim_polys))
        
        # ---- 6. Return ----------------------------------------------------------
        list(edges = edge_df,
             rim   = rim_polys,
             bbox  = bbox,
             count = n,
             file  = file_path)
}


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

#' CLARA-style medoid clustering under Hamming distance
#'
#' A lightweight re-implementation of CLARA (Kaufman & Rousseeuw 1990) that
#' clusters the rows of an ensemble membership matrix under the \emph{Hamming}
#' distance, i.e. the fraction of spanning trees in which two nodes receive
#' different region labels. It exists because \code{cluster::clara()} offers no
#' Hamming metric. One-hot encoding then Manhattan would be a correct (it equals
#'  \eqn{2\times} Hamming) but memory-heavy alternative. Instead, we subsample
#' and reuse the \eqn{O(Nk)} distance-to-medoid routine.
#'
#' Algorithm: draw \code{samples} subsamples of \code{sampsize} rows; run PAM on
#' each subsample's (small) Hamming dissimilarity matrix; assign \emph{all} N
#' rows to the nearest resulting medoid by Hamming distance; keep the medoid set
#' with the lowest mean distance over all N rows. Following CLARA, the current
#' best medoids are forced into each subsequent subsample. Peak memory is
#' \eqn{O(N k + \text{sampsize}^2)} — no \eqn{N \times N} matrix is built.
#'
#' @param memb_mat Integer matrix (N x n_trees) of per-tree region labels.
#' @param k Integer \eqn{\ge 2}. Number of final clusters (medoids).
#' @param samples Integer \eqn{\ge 1}. Number of subsamples. Default 50. Cost
#'   scales linearly in this; 5–10 is often enough.
#' @param sampsize Integer or \code{NULL}. Rows per subsample. \code{NULL}
#'   selects \code{min(N, 40 + 2k)} (the \code{cluster::clara} default).
#' @param seed Integer or \code{NULL}. If supplied, seeds subsampling and the
#'   caller's RNG state is restored on exit.
#' @param verbose Logical. Progress + ETA. Default \code{TRUE}.
#'
#' @return A list with \code{i.med} (length-k integer vector of global medoid
#'   row indices), \code{clustering} (length-N integer hard labels, columns
#'   aligned to \code{i.med}), \code{cost} (mean Hamming distance to assigned
#'   medoid over all N rows), \code{silinfo} (silhouette info from PAM on the
#'   best subsample, matching \code{clara} semantics; \code{$avg.width} is a
#'   hardened, sample-based diagnostic), and \code{sample} (the winning
#'   subsample's global row indices).
#'
#' @importFrom cluster pam
#' @importFrom parallelDist parallelDist
#' @keywords internal
#' @noRd
.clara_hamming <- function(memb_mat,
                           k,
                           samples  = 50L,
                           sampsize = NULL,
                           seed     = NULL,
                           verbose  = TRUE) {
        
        N <- nrow(memb_mat)
        stopifnot(is.matrix(memb_mat), N >= 2L,
                  length(k) == 1L, k >= 2L, k < N,
                  length(samples) == 1L, samples >= 1L)
        
        if (is.null(sampsize)) sampsize <- min(N, 40L + 2L * k)
        sampsize <- as.integer(min(sampsize, N))
        if (sampsize <= k) {
                stop("`sampsize` (", sampsize, ") must exceed `k` (", k, ").", call. = FALSE)
        }
        
        # RNG: seed locally, restore the caller's state on exit (file-wide convention).
        if (!is.null(seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                set.seed(seed)
        }
        
        if (verbose) {
                message(sprintf(
                        "CLARA (Hamming): %d sample(s) of %d rows, N = %d, k = %d ...",
                        samples, sampsize, N, k
                ))
        }
        tick <- if (verbose) .make_progress(samples, "CLARA(Hamming)") else function(i) NULL
        
        best_cost       <- Inf
        best_medoids    <- NULL
        best_clustering <- NULL
        best_sample     <- NULL
        best_pam        <- NULL
        
        for (s in seq_len(samples)) {
                # Draw a subsample; force the current best medoids in (CLARA refinement).
                smp <- sample.int(N, sampsize, replace = FALSE)
                if (!is.null(best_medoids)) smp <- union(best_medoids, smp)
                smp <- sort(unique(smp))
                
                # PAM on the small (sampsize x sampsize) Hamming dissimilarity matrix.
                d_smp <- parallelDist::parallelDist(memb_mat[smp, , drop = FALSE],
                                                    method = "hamming")
                pr    <- cluster::pam(d_smp, k = k, diss = TRUE)
                med_global <- smp[pr$id.med]
                
                # Assign ALL N rows to the nearest medoid; cost = mean Hamming to that medoid.
                # .hamming_to_refs() accumulates an N x k matrix in O(N k) memory: no N x N.
                d_to_med <- .hamming_to_refs(memb_mat, med_global)        # N x k, in [0, 1]
                nearest  <- max.col(-d_to_med, ties.method = "first")     # argmin per row
                cost     <- mean(d_to_med[cbind(seq_len(N), nearest)])
                
                if (cost < best_cost) {
                        best_cost       <- cost
                        best_medoids    <- med_global
                        best_clustering <- nearest
                        best_sample     <- smp
                        best_pam        <- pr        # silhouette taken from this subsample's PAM
                }
                tick(s)
        }
        
        list(
                i.med      = best_medoids,
                clustering = as.integer(best_clustering),
                cost       = best_cost,
                silinfo    = best_pam$silinfo,
                sample     = best_sample
        )
}



# Register variables used via non-standard evaluation (ggplot2/dplyr) so
# R CMD check does not raise "no visible binding for global variable" NOTEs.
#' @importFrom utils globalVariables
utils::globalVariables(c(
        # plot_fuzzy_area
        "color", "tric_color", "Probability", "Dominant_Class", "Max_Prob", ".data",
        # plot_river_centroids
        "PC1", "PC2", "Region", "TypeID",
        # polygon_to_network
        "from", "connected", "to"
))