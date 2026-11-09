# =============================================================================
# connect_components()  --  greedily bridge disconnected components of a graph
#
# Companion to build_network(): it consumes the same vertex convention
# (integer ids 1..n in feature order, no name-based matching) and is called
# either with `polygons=` (in-memory path) or `centroids=` (chunked path).
#
# Relies on these package-internal utilities defined in utils.r:
#   .log(), .need(), .assert_flag(), .warn_if_longlat()
# =============================================================================

#' Connect disconnected components of a spatial graph
#'
#' Iteratively bridges disconnected components of an \pkg{igraph} object by
#' adding an edge between the geographically nearest pair of nodes across
#' components. At each step the smallest remaining component is linked to its
#' nearest neighbour, until the graph forms a single connected component.
#'
#' @param g An \pkg{igraph} object whose vertices correspond, \emph{in order},
#'   to the features in \code{polygons} / \code{centroids}. No name-based
#'   matching is performed; vertex \code{i} is feature \code{i}.
#' @param polygons Optional \pkg{sf} polygon object, one feature per vertex
#'   of \code{g}. Used to derive \code{centroids} when the latter is not
#'   supplied, and, whenever supplied, to measure bridge length as the
#'   shortest boundary-to-boundary distance.
#' @param centroids Optional \pkg{sf} point object or \pkg{sf} geometry column
#'   (\code{sfc}) giving one centroid per vertex of \code{g}. If \code{NULL},
#'   centroids are computed from \code{polygons} via \code{\link[sf]{st_centroid}}.
#'   Whatever is supplied is reduced to its geometry, so an \code{sf} object and
#'   its \code{sfc} are treated identically.
#' @param method Character; nearest-neighbour backend. One of:
#'   \describe{
#'     \item{\code{"rann"}}{kd-tree lookup via \code{\link[RANN]{nn2}}. Fastest,
#'       but operates on planar coordinates and ignores the CRS. If \pkg{RANN}
#'       is not installed, \code{"sf"} is used instead, with a message.}
#'     \item{\code{"nngeo"}}{\code{\link[nngeo]{st_nn}}; slower, but CRS-aware
#'       (geodesic distances for lon/lat data).}
#'     \item{\code{"sf"}}{Brute-force \code{\link[sf]{st_distance}}; no extra
#'       dependencies. Sub-samples to 500 points per side for large components
#'       (deterministically, evenly spaced), so the chosen bridge is only
#'       approximately nearest.}
#'   }
#' @param verbose Logical scalar; if \code{TRUE} (default) progress is reported
#'   via \code{message()}.
#'
#' @return The input graph \code{g} with additional bridge edges, guaranteed to
#'   form a single connected component. Every edge carries the logical
#'   attribute \code{bridge} (\code{TRUE} for added edges, \code{FALSE}
#'   otherwise) and the numeric attribute \code{bridge_dist} (bridge length in
#'   metres; \code{NA} on non-bridges). The graph attribute
#'   \code{bridge_dist_method} records how the length was measured and
#'   \code{bridge_search_method} which backend chose the pairs. If
#'   \code{g} is already connected only the attributes are added.
#'
#' @details
#' The greedy "smallest component first" strategy is cheap -- few query points
#' per iteration -- and tends to produce a small set of bridge edges, but is not
#' guaranteed to yield a Euclidean minimum spanning forest of the components.
#'
#' Either \code{polygons} or \code{centroids} must be supplied; if both are
#' given, \code{centroids} is used to choose the bridged pair. The number of
#' centroids (and of polygons, if given) must equal the number of vertices of
#' \code{g}, otherwise bridges would be added between the wrong nodes; a
#' mismatch is a hard error.
#'
#' @section Bridge length:
#' The bridged pair is always chosen by nearest centroid. Its recorded length
#' (\code{bridge_dist}) depends on what is supplied:
#' \describe{
#'   \item{\code{"boundary"}}{\code{polygons} supplied (in-memory
#'     \code{build_network()}): shortest distance between the two
#'     polygons. A topology gap or corner contact gets a length near 0.}
#'   \item{\code{"centroid"}}{only \code{centroids} supplied (chunked
#'     \code{build_network()}): centroid-to-centroid distance, which is
#'     always at least the boundary distance, so the same bridge receives a
#'     larger penalty in \code{\link{add_edge_weight}}.}
#' }
#' Lengths are in metres for projected and geographic CRSs (great-circle for
#' lon/lat), and in raw coordinate units when the data have no CRS.
#'
#' @section Coordinate reference system:
#' \code{method = "rann"} measures distance in raw coordinate space. For a
#' geographic (lon/lat) CRS those are degrees, so the "nearest" pair is only
#' approximate; use \code{"nngeo"} for geodesically correct distances or project
#' the data first.
#'
#' @seealso \code{\link[igraph]{components}}
#'
#' @examples
#' \dontrun{
#' library(sf); library(igraph)
#' g_connected <- connect_components(g, centroids = pts, method = "rann")
#' }
#'
#' @importFrom sf st_centroid st_geometry st_coordinates st_distance st_crs
#' @importFrom igraph is_igraph components vcount add_edges edge_attr set_graph_attr
#'
#' @export
connect_components <- function(g,
                               polygons  = NULL,
                               centroids = NULL,
                               method    = c("rann", "nngeo", "sf"),
                               verbose   = TRUE) {
        
        # ---- Argument validation ---------------------------------------------------
        if (!igraph::is_igraph(g)) {
                stop("`g` must be an igraph object.", call. = FALSE)
        }
        .assert_flag(verbose, "verbose")
        method <- match.arg(method)
        # RANN is in Suggests; locked-down machines and HPC nodes often lack it.
        if (method == "rann" && !.has_pkg("RANN")) {
                message("Package 'RANN' is not installed; using method = \"sf\" for bridges ",
                        "(approximate for components above 500 units). Install RANN for ",
                        "exact nearest pairs.")
                method <- "sf"
        }
        
        if (is.null(polygons) && is.null(centroids)) {
                stop("Either `polygons` or `centroids` must be supplied.", call. = FALSE)
        }
        
        # ---- Bridge attributes (set even when nothing is bridged) ------------------
        # Every edge carries `bridge`/`bridge_dist`, so downstream code never has to
        # distinguish "no bridges" from "graph built before bridges were flagged".
        if (is.null(igraph::edge_attr(g, "bridge")))
                igraph::E(g)$bridge <- rep(FALSE, igraph::ecount(g))
        if (is.null(igraph::edge_attr(g, "bridge_dist")))
                igraph::E(g)$bridge_dist <- rep(NA_real_, igraph::ecount(g))
        dist_method <- if (!is.null(polygons)) "boundary" else "centroid"
        
        # ---- Early exit ------------------------------------------------------------
        # Cheap connectivity test before any (potentially expensive) centroid work.
        comp <- igraph::components(g)
        if (comp$no <= 1L) {
                .log(verbose, "Graph is already connected.")
                return(g)
        }
        g <- igraph::set_graph_attr(g, "bridge_dist_method", dist_method)
        # Which backend chose the bridges: "sf" may pick other pairs than "rann".
        g <- igraph::set_graph_attr(g, "bridge_search_method", method)
        
        # ---- Derive / normalise centroids -----------------------------------------
        # `centroids` wins over `polygons` when both are given. Reduce to a geometry
        # column (sfc) so all feature indexing below is element-wise -- indexing an
        # `sf` object with a single `[idx]` selects COLUMNS, which would silently
        # bridge the wrong nodes for the "nngeo"/"sf" backends.
        if (is.null(centroids)) {
                if (!inherits(polygons, "sf")) {
                        stop("`polygons` must be an sf object when `centroids` is not supplied.",
                             call. = FALSE)
                }
                centroids <- suppressWarnings(sf::st_centroid(sf::st_geometry(polygons)))
        }
        centroids <- sf::st_geometry(centroids)   # idempotent for sfc input
        
        # ---- Alignment check -------------------------------------------------------
        if (length(centroids) != igraph::vcount(g)) {
                stop(sprintf(
                        "Centroid count (%d) does not match the number of graph vertices (%d).",
                        length(centroids), igraph::vcount(g)), call. = FALSE)
        }
        poly_geom <- NULL
        if (!is.null(polygons)) {
                poly_geom <- sf::st_geometry(polygons)
                if (length(poly_geom) != igraph::vcount(g)) {
                        stop(sprintf(
                                "Polygon count (%d) does not match the number of graph vertices (%d).",
                                length(poly_geom), igraph::vcount(g)), call. = FALSE)
                }
        }
        if (identical(method, "rann")) .warn_if_longlat(sf::st_crs(centroids), verbose)
        
        .log(verbose, "Connecting %d components...", comp$no)
        
        # Coordinate matrix extracted once; vertices/centroids never change (only
        # edges are added), so these indices stay valid throughout.
        coords <- sf::st_coordinates(centroids)
        
        # ---- Greedy bridging -------------------------------------------------------
        n_added  <- 0L
        max_iter <- comp$no   # exactly (comp$no - 1) bridges suffice; cap as a guard
        
        repeat {
                comp <- igraph::components(g)        # recompute: merging changes membership
                if (comp$no == 1L) break
                if (n_added >= max_iter) {
                        warning("Reached the iteration cap without full connection; ",
                                "returning a partially connected graph.", call. = FALSE)
                        break
                }
                
                n_added <- n_added + 1L
                
                # Search outward from the smallest component (fewest query points).
                comp_sizes  <- table(comp$membership)
                smallest_id <- as.integer(names(which.min(comp_sizes)))
                idx_small   <- which(comp$membership == smallest_id)
                idx_other   <- which(comp$membership != smallest_id)
                
                # Nearest pair across the boundary; from in smallest, to elsewhere -> the
                # endpoints are always in different components, so no self-loops or
                # parallel edges can be introduced.
                bridge <- find_nearest_pair(coords, idx_small, idx_other, method,
                                            centroids = centroids)
                
                # Bridge length in metres: boundary-to-boundary when polygons are
                # available, centroid-to-centroid otherwise (see "Bridge length").
                d_bridge <- if (!is.null(poly_geom)) {
                        .pair_dist_m(poly_geom[bridge$from], poly_geom[bridge$to])
                } else {
                        .centroid_pair_dist_m(centroids, bridge$from, bridge$to)
                }
                
                # Add by integer vertex id (graph carries no `name` attribute, so numeric
                # references are interpreted as ids -- matching build_network()).
                g <- igraph::add_edges(g, c(bridge$from, bridge$to),
                                       attr = list(bridge = TRUE, bridge_dist = d_bridge))
                
                if (verbose && (n_added %% 10L == 0L || n_added <= 5L || comp$no <= 3L)) {
                        .log(verbose, "  [%d] Linked component %d (%d nodes) -> %d remaining",
                             n_added, smallest_id, length(idx_small), comp$no - 1L)
                }
        }
        
        .log(verbose, "Done: added %d bridge edges (length by %s distance). Graph is connected.",
             n_added, dist_method)
        g
}


#' Nearest pair of points between two index sets
#'
#' Internal helper used by \code{\link{connect_components}} to locate the
#' closest pair of nodes across two disjoint index sets.
#'
#' @param coords Numeric matrix of coordinates (one row per node) from
#'   \code{\link[sf]{st_coordinates}}. Used by \code{method = "rann"}.
#' @param idx_a,idx_b Integer vectors of row indices into \code{coords} (and
#'   feature indices into \code{centroids}) defining the two sets to bridge.
#' @param method Character; one of \code{"rann"}, \code{"nngeo"}, \code{"sf"}.
#' @param centroids \pkg{sf} geometry column (\code{sfc}); required for
#'   \code{"nngeo"} and \code{"sf"}, ignored for \code{"rann"}.
#'
#' @return A list with elements \code{from} (drawn from \code{idx_a}) and
#'   \code{to} (drawn from \code{idx_b}) giving the nearest pair.
#'
#' @keywords internal
#' @noRd
find_nearest_pair <- function(coords, idx_a, idx_b, method, centroids = NULL) {
        
        switch(method,
               
               # --- Method 1: RANN (fastest; pure coordinate kd-tree lookup) ------------
               rann = {
                       .need("RANN")
                       # For each point in idx_a, find its nearest neighbour among idx_b.
                       nn <- RANN::nn2(
                               data  = coords[idx_b, , drop = FALSE],
                               query = coords[idx_a, , drop = FALSE],
                               k     = 1L
                       )
                       best <- which.min(nn$nn.dists)
                       list(from = idx_a[best], to = idx_b[nn$nn.idx[best, 1L]])
               },
               
               # --- Method 2: nngeo (sf-native, CRS/geodesic aware) ---------------------
               nngeo = {
                       .need("nngeo")
                       nn <- nngeo::st_nn(
                               centroids[idx_a], centroids[idx_b],
                               k = 1L, progress = FALSE, returnDist = TRUE
                       )
                       best <- which.min(unlist(nn$dist))
                       list(from = idx_a[best], to = idx_b[nn$nn[[best]]])
               },
               
               # --- Method 3: sf (no extra deps, slowest) -------------------------------
               sf = {
                       # Cap the cross-distance matrix size with a DETERMINISTIC, evenly spaced
                       # sub-sample so results are reproducible across runs (the previous
                       # random sample() made bridge choice non-reproducible).
                       max_pts <- 500L
                       thin <- function(idx) {
                               if (length(idx) > max_pts) {
                                       idx[round(seq(1L, length(idx), length.out = max_pts))]
                               } else {
                                       idx
                               }
                       }
                       sub_a <- thin(idx_a)
                       sub_b <- thin(idx_b)
                       
                       # st_distance returns a units matrix; strip units for plain comparison.
                       dmat <- matrix(
                               as.numeric(sf::st_distance(centroids[sub_a], centroids[sub_b])),
                               nrow = length(sub_a)
                       )
                       min_pos <- which(dmat == min(dmat), arr.ind = TRUE)[1L, ]
                       list(from = sub_a[min_pos[1L]], to = sub_b[min_pos[2L]])
               }
        )
}