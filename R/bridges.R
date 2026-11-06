# =============================================================================
# Bridge edges: distances, penalty and exclusion (package-internal)
#
# Bridges are the edges connect_components() adds to join disconnected
# components (islands, topology gaps). They are flagged at graph creation
# (E(g)$bridge, E(g)$bridge_dist in metres), penalised by distance in
# add_edge_weight(), and can be excluded from Level 3 neighbourhoods.
#
# Distance definition differs by path (recorded as graph attribute
# `bridge_dist_method`):
#   * in-memory polygon_to_network(): "boundary" -- shortest polygon-to-polygon
#     distance between the two bridged units;
#   * chunked polygon_to_network():  "centroid" -- centroid-to-centroid
#     distance, because the polygons are no longer held in memory.
# Centroid distance is always >= boundary distance for the same pair, so the
# chunked path penalises the same bridge more.
# =============================================================================

#' Strip units from an sf distance, converting to metres
#'
#' @param d Output of \code{sf::st_distance()} (a \pkg{units} object) or a
#'   plain numeric vector. Plain numerics (no CRS) are returned unchanged.
#' @return Numeric vector in metres (or raw coordinate units without a CRS).
#' @keywords internal
#' @noRd
.to_metres <- function(d) {
        if (inherits(d, "units"))
                return(as.numeric(units::set_units(d, "m", mode = "standard")))
        as.numeric(d)
}

#' Metres per coordinate unit of a projected CRS
#'
#' @param crs An \code{sf::st_crs()} object.
#' @return Numeric scalar; 1 when the CRS is missing or has no linear unit.
#' @keywords internal
#' @noRd
.crs_metres_per_unit <- function(crs) {
        u <- tryCatch(crs$ud_unit, error = function(e) NULL)
        if (is.null(u)) return(1)
        .to_metres(u)
}

#' Element-wise distance between paired geometries, in metres
#'
#' @param geom_a,geom_b Equal-length \pkg{sf} geometry columns (sfc).
#' @return Numeric vector of length \code{length(geom_a)}.
#' @keywords internal
#' @noRd
.pair_dist_m <- function(geom_a, geom_b) {
        stopifnot(length(geom_a) == length(geom_b))
        if (!length(geom_a)) return(numeric(0))
        .to_metres(sf::st_distance(geom_a, geom_b, by_element = TRUE))
}

#' Centroid-to-centroid distance for vertex pairs, in metres
#'
#' Projected (or CRS-less) data use planar coordinates directly, which is fast
#' enough for every edge of a continental graph. Geographic data use
#' \code{sf::st_distance()} (great-circle metres).
#'
#' @param centroids \pkg{sf} point geometry column (sfc), one per vertex.
#' @param from,to Integer vertex indices.
#' @return Numeric vector of distances.
#' @keywords internal
#' @noRd
.centroid_pair_dist_m <- function(centroids, from, to) {
        if (!length(from)) return(numeric(0))
        crs <- sf::st_crs(centroids)
        if (isTRUE(sf::st_is_longlat(crs)))
                return(.pair_dist_m(centroids[from], centroids[to]))
        xy <- sf::st_coordinates(centroids)[, 1:2, drop = FALSE]
        sqrt(rowSums((xy[from, , drop = FALSE] - xy[to, , drop = FALSE])^2)) *
                .crs_metres_per_unit(crs)
}

#' Default bridge length scale (lambda) from contiguity edges
#'
#' Median, over all units with at least one contiguity neighbour, of the
#' centroid distance to the nearest such neighbour. A fixed rule, not a tuned
#' value: it follows the spatial resolution of the units.
#'
#' @param centroids \pkg{sf} point geometry column (sfc), one per vertex.
#' @param from,to Integer vertex indices of the contiguity (non-bridge) edges.
#' @return Numeric scalar in metres, or \code{NA_real_} without usable edges.
#' @keywords internal
#' @noRd
.default_bridge_scale <- function(centroids, from, to) {
        if (!length(from)) return(NA_real_)
        d   <- .centroid_pair_dist_m(centroids, from, to)
        v   <- c(from, to)
        dd  <- c(d, d)
        ok  <- is.finite(dd)
        v   <- v[ok]; dd <- dd[ok]
        if (!length(dd)) return(NA_real_)
        ord <- order(v, dd)                          # nearest neighbour first per vertex
        nn  <- dd[ord][!duplicated(v[ord])]
        nn  <- nn[nn > 0]
        if (!length(nn)) return(NA_real_)
        stats::median(nn)
}

#' Validate a bridge length scale
#' @keywords internal
#' @noRd
.assert_bridge_scale <- function(x, nm = "bridge_scale") {
        if (is.null(x)) return(invisible(NULL))
        if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0)
                stop(sprintf(paste0("`%s` must be NULL or a single positive number in metres ",
                                    "(Inf switches the bridge penalty off)."), nm),
                     call. = FALSE)
        invisible(NULL)
}

#' Distance-based penalty on bridge edge weights
#'
#' \deqn{w_b = w_b^{env} + s\,(1 - e^{-d_b / \lambda})}
#' with \eqn{s} the largest finite non-bridge weight. Additive, so a bridge
#' between environmentally identical units (\eqn{w^{env} = 0}) is still
#' penalised. \eqn{d \approx 0} (snapping gaps, corner contacts) adds about
#' nothing; \eqn{d \gg \lambda} makes the bridge as dissimilar as the most
#' dissimilar real neighbour pair. \code{scale = Inf} returns \code{w}
#' unchanged.
#'
#' @param w Numeric vector of environmental edge weights.
#' @param bridge Logical vector (or \code{NULL}), same length as \code{w}.
#' @param dist Numeric vector of bridge distances in metres (NA allowed on
#'   non-bridges).
#' @param scale Positive numeric scalar, the length scale lambda in metres.
#' @return Numeric vector of penalised weights.
#' @keywords internal
#' @noRd
.bridge_penalty <- function(w, bridge, dist, scale) {
        if (is.null(bridge)) return(w)
        bridge <- bridge %in% TRUE
        if (!any(bridge) || is.infinite(scale)) return(w)
        if (is.null(dist) || anyNA(dist[bridge]))
                stop("Bridge edges lack `bridge_dist`. Rebuild the graph with ",
                     "polygon_to_network() or connect_components().", call. = FALSE)
        ref <- w[!bridge & is.finite(w)]
        s   <- if (length(ref)) max(ref) else max(c(0, w[is.finite(w)]))
        w[bridge] <- w[bridge] + s * (1 - exp(-dist[bridge] / scale))
        w
}

#' Edge weights of a graph recomputed from an environment matrix
#'
#' Euclidean distance between the rows of \code{X} at each edge's end units,
#' plus the bridge penalty stored on \code{graph}. This is the weight function
#' of \code{\link{add_edge_weight}}, so \code{compute_null_baseline()} uses the
#' same function on permuted rows. Graphs without a stored
#' \code{bridge_scale} are treated as unpenalised.
#'
#' @param graph An \pkg{igraph} object.
#' @param X Numeric matrix, one row per vertex.
#' @return Numeric vector, one weight per edge in \code{E(graph)} order.
#' @keywords internal
#' @noRd
.edge_weights_from_env <- function(graph, X) {
        ep <- igraph::ends(graph, igraph::E(graph), names = FALSE)
        w  <- unname(sqrt(rowSums((X[ep[, 1], , drop = FALSE] -
                                   X[ep[, 2], , drop = FALSE])^2)))
        sc <- igraph::graph_attr(graph, "bridge_scale")
        if (is.null(sc) || is.na(sc)) sc <- Inf
        .bridge_penalty(w,
                        igraph::edge_attr(graph, "bridge"),
                        igraph::edge_attr(graph, "bridge_dist"),
                        sc)
}

#' Remove bridge edges from a graph (vertices kept)
#'
#' @param graph An \pkg{igraph} object.
#' @param verbose Logical.
#' @return \code{graph} without edges flagged \code{bridge == TRUE}.
#' @keywords internal
#' @noRd
.drop_bridges <- function(graph, verbose = FALSE) {
        br <- igraph::edge_attr(graph, "bridge")
        if (is.null(br)) {
                .log(verbose, paste("`graph` has no `bridge` edge attribute, so no edges",
                                    "are excluded. Rebuild it with polygon_to_network()",
                                    "and add_edge_weight() to flag bridges."))
                return(graph)
        }
        idx <- which(br %in% TRUE)
        if (!length(idx)) return(graph)
        .log(verbose, "Excluding %d bridge edge(s) from neighbourhoods.", length(idx))
        igraph::delete_edges(graph, idx)
}
