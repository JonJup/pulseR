# ===========================================================================*
# build_network()  --  build a spatial contiguity network from polygons
#
# Public API (unchanged return shape) plus four new, fully-defaulted arguments
# (verbose, snap, queen, bbox_buffer). Hardened against the silent-corruption
# and CRS hazards documented in the package NEWS.
# ===========================================================================*

#' Create a Spatial Network from Polygon Geometries
#'
#' Constructs an \pkg{igraph} network object from polygon geometries based on
#' spatial contiguity. Polygons that share a boundary become connected nodes in
#' the network. The function supports two modes: an in-memory mode that takes a
#' single \pkg{sf} object, and a chunked mode that processes multiple polygon
#' files sequentially to handle datasets too large to fit in memory.
#'
#' @param polygons When \code{chunked = FALSE}, an \pkg{sf} object containing
#'   polygon geometries. When \code{chunked = TRUE}, a character vector of file
#'   paths to GeoParquet files containing polygon chunks. Each polygon becomes a
#'   node in the resulting network.
#' @param chunked Logical scalar. If \code{FALSE} (default) polygons are
#'   processed in memory as a single \pkg{sf} object. If \code{TRUE} polygons
#'   are read chunk by chunk from the file paths supplied in \code{polygons},
#'   with cross-chunk edges detected via bounding-box adjacency. Use the chunked
#'   path for datasets too large to hold in memory at once.
#' @param min_shared_length Numeric scalar \eqn{\ge 0}. Minimum length of shared
#'   boundary required for two polygons to be considered neighbours, expressed
#'   in the coordinate reference system's units (metres for a projected CRS,
#'   degrees for a geographic CRS). The default \code{0} keeps every touching
#'   pair. Any positive value triggers a per-candidate-edge boundary-length
#'   computation (see \emph{Details}), which is exact but considerably slower.
#' @param connect_islands Logical scalar. If \code{TRUE} (default), disconnected
#'   components in the network are joined by bridge edges using
#'   \code{\link{connect_components}}. Bridges are flagged and measured; see
#'   \emph{Bridge edges}.
#' @param control_paths Optional character vector of file paths, only used when
#'   \code{chunked = TRUE}; must be the same length as \code{polygons}. Each
#'   file should contain an identifier column (\code{ID} or \code{id}) listing
#'   the polygon identifiers to retain; polygons whose identifiers are not
#'   listed are filtered out before neighbour detection. Useful when
#'   environmental data are stored alongside polygons in separate files and only
#'   polygons with matching environmental records should enter the network.
#' @param return_combined_polygons Logical scalar. If \code{TRUE} and
#'   \code{chunked = TRUE}, the function additionally returns a single \pkg{sf}
#'   object combining all polygons processed across chunks (restricted to the
#'   columns common to every chunk). Ignored in in-memory mode. Default
#'   \code{FALSE}.
#' @param queen Logical scalar. Contiguity rule applied \emph{consistently} to
#'   within-chunk and cross-chunk edges. \code{FALSE} (default) uses Rook
#'   contiguity (polygons must share a boundary \emph{segment}); \code{TRUE} uses
#'   Queen contiguity (a single shared boundary \emph{point} is sufficient).
#' @param snap Numeric scalar \eqn{\ge 0}. Snapping tolerance passed to
#'   \code{spdep::poly2nb()} to absorb minor topological inconsistencies,
#'   expressed in CRS units. Default \code{1e-6}.
#' @param bbox_buffer Numeric scalar \eqn{\ge 0}, only used when
#'   \code{chunked = TRUE}. Tolerance added to each side of a chunk's bounding
#'   box when deciding whether two chunks are close enough to share edges,
#'   expressed in CRS units. Choose a value at least as large as \code{snap};
#'   too small a value risks missing genuine cross-chunk edges, too large a
#'   value merely costs extra (harmless) pairwise checks. Default \code{0.01}.
#' @param verbose Logical scalar. If \code{TRUE} (default) progress is reported
#'   via \code{message()} (stderr, suppressible with \code{suppressMessages()}).
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{graph}{An \pkg{igraph} object representing the spatial network.
#'       Nodes correspond to polygons (vertex ids \code{1..n} in input order)
#'       and edges represent spatial contiguity or bridges. The graph is
#'       undirected and simplified (no multi-edges or self-loops). Edge
#'       attributes \code{bridge} and \code{bridge_dist}, and graph attributes
#'       \code{bridge_scale_default} and (when bridges were added)
#'       \code{bridge_dist_method} are described under \emph{Bridge edges}.}
#'     \item{polygons}{Returned only when \code{chunked = FALSE}: the cleaned
#'       \pkg{sf} polygon object actually used to build the network (after
#'       validity repair and empty-geometry removal).}
#'     \item{edge_data}{A data frame with integer columns \code{from} and
#'       \code{to} giving connected polygon pairs as 1-based global indices,
#'       logical \code{bridge} and numeric \code{bridge_dist} (metres,
#'       \code{NA} for contiguity edges). Bridge rows come last.}
#'     \item{chunk_meta}{Returned only when \code{chunked = TRUE}: a list with
#'       one entry per chunk recording file paths, node count, global index
#'       range, and bounding box.}
#'     \item{combined_polygons}{Returned only when \code{chunked = TRUE} and
#'       \code{return_combined_polygons = TRUE}: a single \pkg{sf} object
#'       combining all polygons processed across chunks.}
#'   }
#'
#' @details
#' Neighbour detection within a chunk (and in memory) uses
#' \code{spdep::poly2nb()} with the requested \code{queen} rule and \code{snap}
#' tolerance. Cross-chunk neighbours use the same \code{poly2nb()} rule and
#' \code{snap}, run on the units of each chunk pair that lie near the other
#' chunk, so the two phases agree. The resulting graph is undirected, since
#' spatial contiguity is symmetric.
#'
#' In chunked mode the function operates in three phases:
#' \enumerate{
#'   \item \strong{Within-chunk edges.} Each chunk is read once through a single
#'     loader (\code{.load_chunk()}) that validates geometries, optionally
#'     filters by \code{control_paths}, and removes empties. Neighbours are
#'     detected locally and global index offsets are tracked so node ids remain
#'     unique across the full network.
#'   \item \strong{Cross-chunk edges.} Chunk pairs whose (buffered) bounding
#'     boxes overlap or touch are reloaded \emph{through the same loader} (so the
#'     row order matches phase 1 exactly) and the chunk-to-chunk predicate is
#'     evaluated. Pairs whose bounding boxes are clearly disjoint are skipped,
#'     which is what makes the chunked path scalable.
#'   \item \strong{Assembly.} Edges are concatenated, an \pkg{igraph} object is
#'     built directly from integer vertex ids and simplified, and -- when
#'     \code{connect_islands = TRUE} and more than one component remains --
#'     \code{connect_components()} is called using polygon centroids collected
#'     across chunks.
#' }
#'
#' If \code{min_shared_length > 0}, each candidate edge is checked by computing
#' the length of the shared boundary (\code{sf::st_intersection()} reduced to its
#' linear parts) and dropped when below the threshold. This is exact but scales
#' with the number of candidate edges, so it is off by default.
#'
#' If the input polygons contain an \code{ID} (or \code{id}) column, these
#' identifiers are stored as the \code{polygon_id} vertex attribute. In chunked
#' mode, chunks without such a column give \code{NA} ids.
#'
#' @section Bridge edges:
#' Edges added by \code{connect_islands = TRUE} are not contiguity. Every edge
#' carries \code{bridge} (\code{TRUE} for bridges) and \code{bridge_dist}
#' (bridge length in metres, \code{NA} otherwise). The length is measured
#' differently in the two modes, recorded in the graph attribute
#' \code{bridge_dist_method}:
#' \describe{
#'   \item{In memory (\code{"boundary"})}{shortest polygon-to-polygon
#'     distance.}
#'   \item{Chunked (\code{"centroid"})}{centroid-to-centroid distance, because
#'     the polygons are no longer in memory. It is always at least the
#'     boundary distance, so the two modes do not give identical bridge
#'     penalties for the same data.}
#' }
#' The graph attribute \code{bridge_scale_default} holds the default length
#' scale for the bridge penalty in \code{\link{add_edge_weight}}: the median,
#' over units with a contiguity neighbour, of the centroid distance to the
#' nearest such neighbour (metres). \code{\link{compute_all_signatures}}
#' excludes bridges from Level 3 neighbourhoods by default.
#'
#' @section Coordinate reference system:
#' All distance-like arguments (\code{min_shared_length}, \code{snap},
#' \code{bbox_buffer}) and centroid-based component joining are interpreted in
#' the CRS of the data. For a geographic (longitude/latitude) CRS these are
#' degrees and centroids are planar approximations; a warning is emitted. Project
#' to a suitable equal-area or local CRS for metric thresholds. In chunked mode
#' every chunk must share one CRS; a mismatch is a hard error.
#'
#' @note This function depends on a package-internal \code{connect_components()}
#'   helper (used when \code{connect_islands = TRUE}) which must accept either a
#'   \code{polygons} or a \code{centroids} argument alongside the graph.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # --- In-memory mode ---
#' polys   <- st_read("my_polygons.shp")
#' network <- build_network(polys)
#' plot(network$graph)
#'
#' # Queen contiguity, no auto-connection of islands
#' network <- build_network(polys, queen = TRUE, connect_islands = FALSE)
#'
#' # Only count neighbours sharing at least 50 m of boundary (projected CRS)
#' network <- build_network(polys, min_shared_length = 50)
#'
#' # --- Chunked mode across multiple GeoParquet files ---
#' files   <- list.files("data/chunks", pattern = "\\.parquet$", full.names = TRUE)
#' network <- build_network(files, chunked = TRUE,
#'                               return_combined_polygons = TRUE)
#'
#' # Chunked mode with paired environmental control files
#' env_files <- list.files("data/env", pattern = "\\.parquet$", full.names = TRUE)
#' network   <- build_network(files, chunked = TRUE, control_paths = env_files)
#' }
#'
#' @seealso \code{\link[spdep]{poly2nb}} for neighbour detection,
#'   \code{\link[igraph]{make_empty_graph}} for graph construction.
#'
#' @importFrom spdep poly2nb
#' @importFrom sf st_make_valid st_is_empty st_centroid st_geometry st_bbox st_intersects st_as_sfc st_length st_intersection st_collection_extract st_crs st_is_longlat
#' @importFrom igraph make_empty_graph add_edges set_vertex_attr vcount ecount components simplify ends set_graph_attr E
#'
#' @export
build_network <- function(polygons,
                               chunked                  = FALSE,
                               min_shared_length        = 0,
                               connect_islands          = TRUE,
                               control_paths            = NULL,
                               return_combined_polygons = FALSE,
                               queen                    = FALSE,
                               snap                     = 1e-6,
                               bbox_buffer              = 0.01,
                               verbose                  = TRUE) {
  
  # spdep, sf and igraph are Imports. The heavy I/O packages (sfarrow, arrow)
  # are only needed for chunked mode and are checked there (Suggests).
  
  # ---- Scalar argument validation ---- *-
  # Fail fast and informatively rather than letting a malformed value surface as
  # a cryptic error deep inside spdep / sf.
  .assert_flag(chunked,                  "chunked")
  .assert_flag(connect_islands,          "connect_islands")
  .assert_flag(return_combined_polygons, "return_combined_polygons")
  .assert_flag(queen,                    "queen")
  .assert_flag(verbose,                  "verbose")
  .assert_nonneg(min_shared_length, "min_shared_length")
  .assert_nonneg(snap,              "snap")
  .assert_nonneg(bbox_buffer,       "bbox_buffer")
  
  # ---- Mode-specific input validation ---------------------------------------
  if (chunked) {
    if (!is.character(polygons) || length(polygons) == 0L) {
      stop("When `chunked = TRUE`, `polygons` must be a non-empty character ",
           "vector of file paths.", call. = FALSE)
    }
    if (!all(file.exists(polygons))) {
      stop("These `polygons` files do not exist: ",
           paste(polygons[!file.exists(polygons)], collapse = ", "),
           call. = FALSE)
    }
    if (!is.null(control_paths)) {
      if (length(control_paths) != length(polygons)) {
        stop("`control_paths` must be the same length as `polygons`.",
             call. = FALSE)
      }
      if (!all(file.exists(control_paths))) {
        stop("These `control_paths` files do not exist: ",
             paste(control_paths[!file.exists(control_paths)], collapse = ", "),
             call. = FALSE)
      }
    }
  } else {
    if (!inherits(polygons, "sf")) {
      stop("When `chunked = FALSE`, `polygons` must be an sf object.",
           call. = FALSE)
    }
    if (nrow(polygons) == 0L) {
      stop("`polygons` contains no features.", call. = FALSE)
    }
    # These arguments only have meaning in chunked mode; warn instead of
    # silently ignoring so the caller learns their value had no effect.
    if (!is.null(control_paths)) {
      warning("`control_paths` is ignored when `chunked = FALSE`.",
              call. = FALSE)
    }
    if (return_combined_polygons) {
      warning("`return_combined_polygons` is ignored when `chunked = FALSE`.",
              call. = FALSE)
    }
  }
  
  # ---- Dispatch ---- *
  if (!chunked) {
    return(.network_in_memory(
      polygons          = polygons,
      min_shared_length = min_shared_length,
      connect_islands   = connect_islands,
      queen             = queen,
      snap              = snap,
      verbose           = verbose
    ))
  }
  
  .network_chunked(
    polygon_paths            = polygons,
    control_paths            = control_paths,
    connect_islands          = connect_islands,
    return_combined_polygons = return_combined_polygons,
    min_shared_length        = min_shared_length,
    queen                    = queen,
    snap                     = snap,
    bbox_buffer              = bbox_buffer,
    verbose                  = verbose
  )
}

# =============================================================================*
# connect_components()  --  greedily bridge disconnected components of a graph
#
# Companion to build_network(): it consumes the same vertex convention
# (integer ids 1..n in feature order, no name-based matching) and is called
# either with `polygons=` (in-memory path) or `centroids=` (chunked path).
# =============================================================================*

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
#' 
#' xy <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(pulse_grid)))
#' keep <- xy[,1] < 8000 | xy[,2] >11000
#' sub <- pulse_grid[keep,]
#' net <- build_network(sub, connect_islands = FALSE, verbose = FALSE)
#' 
#' g_connected <- connect_components(
#' net$graph, 
#' centroids = xy, 
#' method = "rann",
#' verbose = FALSE)
#' 
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
  
  # ---- Argument validation 
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
  
  # ---- Bridge attributes (set even when nothing is bridged) 
  # Every edge carries `bridge`/`bridge_dist`, so downstream code never has to
  # distinguish "no bridges" from "graph built before bridges were flagged".
  if (is.null(igraph::edge_attr(g, "bridge")))
    igraph::E(g)$bridge <- rep(FALSE, igraph::ecount(g))
  if (is.null(igraph::edge_attr(g, "bridge_dist")))
    igraph::E(g)$bridge_dist <- rep(NA_real_, igraph::ecount(g))
  dist_method <- if (!is.null(polygons)) "boundary" else "centroid"
  
  # ---- Early exit
  # Cheap connectivity test before any (potentially expensive) centroid work.
  comp <- igraph::components(g)
  if (comp$no <= 1L) {
    .log(verbose, "Graph is already connected.")
    return(g)
  }
  g <- igraph::set_graph_attr(g, "bridge_dist_method", dist_method)
  # Which backend chose the bridges: "sf" may pick other pairs than "rann".
  g <- igraph::set_graph_attr(g, "bridge_search_method", method)
  
  # ---- Derive / normalise centroids 
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
  
  # ---- Alignment check 
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
  
  # ---- Greedy bridging 
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
    bridge <- .find_nearest_pair(coords, idx_small, idx_other, method,
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




# ===========================================================================*
# Internal: in-memory implementation
# ===========================================================================*
.network_in_memory <- function(polygons, min_shared_length, connect_islands,
                               queen, snap, verbose) {
  
  .warn_if_longlat(sf::st_crs(polygons), verbose)
  
  # ---- Geometry hygiene before neighbour detection --------------------------
  # poly2nb is intolerant of invalid and empty geometries, so repair and drop
  # them first. Indices below are taken AFTER this step, so the returned
  # `polygons` object (also cleaned) stays aligned with the graph vertices.
  polygons   <- sf::st_make_valid(polygons)
  empty      <- sf::st_is_empty(polygons)
  if (any(empty)) {
    warning(sprintf("Dropping %d empty geometry/geometries.", sum(empty)),
            call. = FALSE)
    polygons <- polygons[!empty, ]
  }
  n <- nrow(polygons)
  if (n == 0L) stop("No valid, non-empty geometries remain.", call. = FALSE)
  
  .log(verbose, "Processing %d polygons with spdep::poly2nb()...", n)
  
  # ---- Contiguity neighbours -> from/to edgelist (local == global here) -----
  nb      <- suppressWarnings(spdep::poly2nb(polygons, queen = queen, snap = snap))
  edge_df <- .nb_to_edgelist(nb)
  
  # Optional shared-boundary-length filtering (exact but O(edges); see Details).
  if (nrow(edge_df) > 0 && min_shared_length > 0) {
    keep    <- .shared_lengths_pairs(sf::st_geometry(polygons)[edge_df$from],
                                     sf::st_geometry(polygons)[edge_df$to]) >=
      min_shared_length
    edge_df <- edge_df[keep, , drop = FALSE]
  }
  
  # ---- Assemble graph (isolated vertices preserved) -------------------------
  g <- .build_graph(edge_df, n_total = n)
  
  # Carry polygon identifiers through as a vertex attribute when available.
  id_col <- .id_col(polygons)
  if (!is.na(id_col)) {
    .check_unique_ids(polygons[[id_col]], id_col)
    g <- igraph::set_vertex_attr(g, "polygon_id", value = polygons[[id_col]])
  }
  
  comp <- igraph::components(g)$no
  .log(verbose, "Network created: %d nodes, %d edges, %d components.",
       igraph::vcount(g), igraph::ecount(g), comp)
  
  # Default bridge length scale from the contiguity edges, before bridging.
  centroids <- suppressWarnings(sf::st_centroid(sf::st_geometry(polygons)))
  g <- .init_bridge_attrs(g, centroids)
  
  if (comp != 1L && connect_islands) {
    g <- connect_components(g, polygons = polygons, centroids = centroids,
                            verbose = verbose)
  }
  
  igraph::E(g)$eid <- seq_len(ecount(g))
  
  list(
    graph     = g,
    polygons  = polygons,   # cleaned object, aligned with vertex ids
    edge_data = .edge_data_with_bridges(edge_df, g)
  )
}


# ===========================================================================*
# Internal: chunked, out-of-core implementation
# ===========================================================================*

#' @keywords internal
#' @noRd
.network_chunked <- function(polygon_paths, control_paths, connect_islands,
                             return_combined_polygons, min_shared_length,
                             queen, snap, bbox_buffer, verbose) {
  
  # Chunked mode pulls in the optional GeoParquet I/O stack.
  .need("sfarrow")
  .need("arrow")
  
  all_edges      <- list()   # list of from/to data frames (global indices)
  all_ids        <- NULL     # polygon_id values in global order
  offset         <- 0L       # running global-index offset across chunks
  chunk_meta     <- vector("list", length(polygon_paths))
  centroids_list <- vector("list", length(polygon_paths))
  combined_polys <- if (return_combined_polygons)
    vector("list", length(polygon_paths)) else NULL
  ref_crs        <- NULL     # CRS of the first non-empty chunk, for consistency
  saw_id         <- FALSE    # whether any chunk carried an id column
  
  # ---- Phase 1: within-chunk edges ---- *
  for (k in seq_along(polygon_paths)) {
    .log(verbose, "Processing chunk %d/%d", k, length(polygon_paths))
    
    polys_k <- .load_chunk(polygon_paths[k],
                           if (!is.null(control_paths)) control_paths[k] else NULL)
    n_k     <- nrow(polys_k)
    
    # Enforce a single CRS across the whole dataset; silent CRS mixing would
    # corrupt cross-chunk predicates and centroid distances.
    if (n_k > 0L) {
      crs_k <- sf::st_crs(polys_k)
      if (is.null(ref_crs)) {
        ref_crs <- crs_k
        .warn_if_longlat(ref_crs, verbose)
      } else if (crs_k != ref_crs) {
        stop(sprintf("Chunk %d has a different CRS than earlier chunks; ",
                     k),
             "all chunks must share one CRS.", call. = FALSE)
      }
    }
    
    # Global indices for this chunk. NB: seq() counts DOWN when n_k == 0, which
    # would silently invent two bogus indices for an empty chunk, so guard it.
    global_idx <- if (n_k > 0L) seq.int(offset + 1L, offset + n_k) else integer(0)
    
    if (return_combined_polygons) combined_polys[[k]] <- polys_k
    
    chunk_meta[[k]] <- list(
      path         = polygon_paths[k],
      control_path = if (!is.null(control_paths)) control_paths[k] else NULL,
      n            = n_k,
      global_idx   = global_idx,
      bbox         = if (n_k > 0L) sf::st_bbox(polys_k) else NULL
    )
    
    if (n_k > 0L) {
      # Within-chunk contiguity, mapped onto global indices.
      nb_k    <- suppressWarnings(spdep::poly2nb(polys_k, queen = queen, snap = snap))
      edges_k <- .nb_to_edgelist(nb_k)
      
      if (nrow(edges_k) > 0 && min_shared_length > 0) {
        keep    <- .shared_lengths_pairs(sf::st_geometry(polys_k)[edges_k$from],
                                         sf::st_geometry(polys_k)[edges_k$to]) >=
          min_shared_length
        edges_k <- edges_k[keep, , drop = FALSE]
      }
      if (nrow(edges_k) > 0) {
        edges_k$from <- edges_k$from + offset
        edges_k$to   <- edges_k$to   + offset
        all_edges[[length(all_edges) + 1L]] <- edges_k
      }
      
      # Polygon ids (in global order). Track whether ANY chunk had them so a
      # partial set does not get half-applied later.
      id_col <- .id_col(polys_k)
      if (!is.na(id_col)) {
        saw_id  <- TRUE
        all_ids <- c(all_ids, as.character(polys_k[[id_col]]))
      } else {
        all_ids <- c(all_ids, rep(NA_character_, n_k))
      }
      
      # Centroids for optional island-connection. Empty chunks contribute an
      # empty sfc so the per-chunk lists stay aligned with global indices.
      centroids_list[[k]] <-
        suppressWarnings(sf::st_centroid(sf::st_geometry(polys_k)))
    }
    
    offset <- offset + n_k
    # Release per-chunk objects promptly: out-of-core processing exists
    # precisely because we cannot afford to hold everything at once.
    if (!return_combined_polygons) rm(polys_k)
    gc(verbose = FALSE)
  }
  
  n_total <- offset
  
  # ---- Phase 2: cross-chunk edges ---- *-
  # Only inspect chunk pairs whose buffered bounding boxes are adjacent; the
  # rest cannot possibly share a boundary, which is what keeps this scalable.
  for (i in seq_along(chunk_meta)) {
    if (chunk_meta[[i]]$n == 0L)
      next
    
    # Cheap bbox screen first; collect adjacent partners j > i.
    j_adj <- integer(0)
    for (j in seq_along(chunk_meta)) {
      if (j <= i || chunk_meta[[j]]$n == 0L)
        next
      if (.bboxes_adjacent(chunk_meta[[i]]$bbox, chunk_meta[[j]]$bbox, buffer = bbox_buffer))
        j_adj <- c(j_adj, j)
    }
    if (length(j_adj) == 0L)
      next
    
    # Load chunk i ONCE, reuse across all its partners.
    polys_i <- .load_chunk(chunk_meta[[i]]$path, chunk_meta[[i]]$control_path)
    
    for (j in j_adj) {
      .log(verbose, "Checking cross-chunk edges: %d x %d", i, j)
      polys_j <- .load_chunk(chunk_meta[[j]]$path, chunk_meta[[j]]$control_path)
      
      # Same rule and snap tolerance as within a chunk (poly2nb), so a
      # near-touching pair on a chunk border is a neighbour, not a bridge.
      pairs <- .cross_chunk_pairs(polys_i, polys_j, queen = queen, snap = snap,
                                  buffer = bbox_buffer)
      
      if (!is.null(pairs) && nrow(pairs) > 0) {
        if (min_shared_length > 0) {
          keep  <- .shared_lengths_pairs(sf::st_geometry(polys_i)[pairs$row],
                                         sf::st_geometry(polys_j)[pairs$col]) >=
            min_shared_length
          pairs <- pairs[keep, , drop = FALSE]
        }
        if (nrow(pairs) > 0) {
          all_edges[[length(all_edges) + 1L]] <- data.frame(from = chunk_meta[[i]]$global_idx[pairs$row],
                                                            to   = chunk_meta[[j]]$global_idx[pairs$col])
        }
      }
      rm(polys_j, pairs)
    }
    rm(polys_i)
    gc(verbose = FALSE)
  }
  # ---- Phase 3: assemble graph ---- *
  all_edges <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, all_edges)
  edge_df   <- if (length(all_edges) > 0) {
    do.call(rbind, all_edges)
  } else {
    data.frame(from = integer(), to = integer())
  }
  
  g <- .build_graph(edge_df, n_total = n_total)
  
  # Attach ids if any chunk had an ID/id column. Chunks without one contribute
  # NA, so polygon_id can be partly NA. Duplicated non-NA ids only war
  if (saw_id && length(all_ids) == n_total && !all(is.na(all_ids))) {
    .check_unique_ids(all_ids, "ID/id")
    g <- igraph::set_vertex_attr(g, "polygon_id", value = all_ids)
  }
  
  comp <- igraph::components(g)$no
  .log(verbose, "Combined network: %d nodes, %d edges, %d components.",
       igraph::vcount(g), igraph::ecount(g), comp)
  
  # do.call(c, ...) concatenates the per-chunk sfc objects in global order;
  # NULL entries (chunks with zero nodes) are dropped harmlessly.
  all_centroids <- if (n_total > 0L)
    do.call(c, Filter(Negate(is.null), centroids_list)) else NULL
  g <- .init_bridge_attrs(g, all_centroids)
  
  if (connect_islands && comp > 1L) {
    # Polygons are not in memory here, so bridge length is centroid-to-centroid
    # (bridge_dist_method = "centroid"); the in-memory path uses boundaries.
    g <- connect_components(g, centroids = all_centroids, verbose = verbose)
  }
  
  igraph::E(g)$eid <- seq_len(ecount(g))
  
  out <- list(
    graph      = g,
    edge_data  = .edge_data_with_bridges(edge_df, g),
    chunk_meta = chunk_meta
  )
  
  if (return_combined_polygons) {
    .log(verbose, "Combining all chunked polygons into a single sf object...")
    # rbind.sf requires matching columns; restrict to the intersection so a
    # schema drift between chunks does not abort the whole run.
    present <- Filter(function(x) inherits(x, "sf") && nrow(x) > 0, combined_polys)
    if (length(present) > 0) {
      common <- Reduce(intersect, lapply(present, names))
      if (length(common) < max(vapply(present, length, integer(1)))) {
        warning("Chunks have differing columns; combined_polygons keeps only ",
                "the common ones: ", paste(common, collapse = ", "),
                call. = FALSE)
      }
      out$combined_polygons <- do.call(rbind, lapply(present, function(x) x[, common]))
    }
  }
  
  out
}


# ===========================================================================*
# Helpers (package-internal)
# ===========================================================================*

#' Identify the polygon identifier column, preferring "ID" then "id"
#' @return Column name, or \code{NA_character_} if neither is present.
#' @keywords internal
#' @noRd
.id_col <- function(x) {
  if ("ID" %in% names(x)) "ID" else if ("id" %in% names(x)) "id" else NA_character_
}

#' Read one polygon chunk, optionally control-filter, repair and drop empties
#'
#' This single loader is used by BOTH the within-chunk and the cross-chunk
#' phases so their row ordering is provably identical; the cross-chunk edge
#' mapping (\code{global_idx[row]}) relies on that invariant.
#'
#' @param path GeoParquet file with polygon geometries.
#' @param control_path Optional file with an identifier column listing the ids
#'   to retain.
#' @return An \pkg{sf} polygon object with empty geometries removed.
#' @keywords internal
#' @noRd
.load_chunk <- function(path, control_path = NULL) {
  polys <- sfarrow::st_read_parquet(path, quiet = TRUE)
  
  if (!is.null(control_path)) {
    poly_c   <- arrow::read_parquet(control_path)
    id_poly  <- .id_col(polys)
    id_ctrl  <- .id_col(poly_c)
    if (is.na(id_poly) || is.na(id_ctrl)) {
      stop("Control filtering requires an 'ID' (or 'id') column in both the ",
           "polygon file and its control file: ", path, call. = FALSE)
    }
    polys <- polys[polys[[id_poly]] %in% poly_c[[id_ctrl]], ]
  }
  
  polys <- sf::st_make_valid(polys)
  polys[!sf::st_is_empty(polys), ]
}

#' Build an undirected, simplified igraph from a global edge list
#'
#' Vertices are integer ids \code{1..n_total} (input order). Constructed without
#' name-based matching to avoid integer/character coercion pitfalls, then
#' simplified to drop any accidental duplicate edges or self-loops.
#'
#' @param edge_df Data frame with integer \code{from}/\code{to} columns.
#' @param n_total Total number of vertices.
#' @return An \pkg{igraph} object.
#' @keywords internal
#' @noRd
.build_graph <- function(edge_df, n_total) {
  g <- igraph::make_empty_graph(n = n_total, directed = FALSE)
  if (nrow(edge_df) > 0) {
    # add_edges() wants a flat vector v1,v2,v3,v4,... = from1,to1,from2,to2,...
    g <- igraph::add_edges(
      g, as.integer(c(rbind(edge_df$from, edge_df$to)))
    )
  }
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

#' Add bridge attributes and the default bridge length scale to a new graph
#'
#' Called before bridging, when every edge is a contiguity edge. Sets
#' \code{bridge = FALSE}, \code{bridge_dist = NA} on all edges and the graph
#' attribute \code{bridge_scale_default} (metres; \code{NA} without edges).
#'
#' @param g An \pkg{igraph} object built by \code{.build_graph()}.
#' @param centroids \pkg{sf} point geometry column, one per vertex, or NULL.
#' @return \code{g} with the attributes set.
#' @keywords internal
#' @noRd
.init_bridge_attrs <- function(g, centroids) {
  ne <- igraph::ecount(g)
  igraph::E(g)$bridge      <- rep(FALSE, ne)
  igraph::E(g)$bridge_dist <- rep(NA_real_, ne)
  lambda <- NA_real_
  if (ne > 0L && !is.null(centroids)) {
    ep     <- igraph::ends(g, igraph::E(g), names = FALSE)
    lambda <- .default_bridge_scale(centroids, ep[, 1], ep[, 2])
  }
  igraph::set_graph_attr(g, "bridge_scale_default", lambda)
}

#' Append bridge edges to the contiguity edge table
#'
#' @param edge_df Contiguity edge table (\code{from}, \code{to}).
#' @param g Final graph carrying the \code{bridge}/\code{bridge_dist} attributes.
#' @return \code{edge_df} with \code{bridge} and \code{bridge_dist} columns
#'   and one extra row per bridge edge.
#' @keywords internal
#' @noRd
.edge_data_with_bridges <- function(edge_df, g) {
  edge_df$bridge      <- rep(FALSE, nrow(edge_df))
  edge_df$bridge_dist <- rep(NA_real_, nrow(edge_df))
  br <- which(igraph::E(g)$bridge %in% TRUE)
  if (!length(br)) return(edge_df)
  ep <- igraph::ends(g, br, names = FALSE)
  rbind(edge_df,
        data.frame(from        = as.integer(pmin(ep[, 1], ep[, 2])),
                   to          = as.integer(pmax(ep[, 1], ep[, 2])),
                   bridge      = TRUE,
                   bridge_dist = igraph::E(g)$bridge_dist[br]))
}


#' Contiguity pairs between two chunks, by the within-chunk rule
#'
#' Runs \code{spdep::poly2nb()} with the same \code{queen} and \code{snap} as
#' the within-chunk phase on the units of each chunk that reach the other
#' chunk's bounding box (widened by \code{buffer + snap}), and keeps the pairs
#' that join the two chunks.
#'
#' @param polys_i,polys_j \pkg{sf} polygon chunks in one CRS.
#' @param queen,snap As in \code{build_network()}.
#' @param buffer Bounding-box buffer in CRS units.
#' @return Data frame with \code{row} (index into \code{polys_i}) and
#'   \code{col} (index into \code{polys_j}), or \code{NULL}.
#' @keywords internal
#' @noRd
.cross_chunk_pairs <- function(polys_i, polys_j, queen, snap, buffer) {
  near <- function(a, b) {
    bb <- sf::st_bbox(b)
    d  <- buffer + snap
    bb[c("xmin", "ymin")] <- bb[c("xmin", "ymin")] - d
    bb[c("xmax", "ymax")] <- bb[c("xmax", "ymax")] + d
    which(lengths(sf::st_intersects(sf::st_geometry(a), sf::st_as_sfc(bb))) > 0L)
  }
  ia <- near(polys_i, polys_j)
  ib <- near(polys_j, polys_i)
  if (!length(ia) || !length(ib)) return(NULL)
  both <- c(sf::st_geometry(polys_i)[ia], sf::st_geometry(polys_j)[ib])
  nb   <- suppressWarnings(spdep::poly2nb(both, queen = queen, snap = snap))
  el   <- .nb_to_edgelist(nb)                      # from < to
  na   <- length(ia)
  x    <- el$from <= na & el$to > na
  if (!any(x)) return(NULL)
  data.frame(row = ia[el$from[x]], col = ib[el$to[x] - na])
}

#' Do two bounding boxes touch or overlap (with a buffer)?
#'
#' @param bb1,bb2 Bounding boxes from \code{sf::st_bbox()}.
#' @param buffer Numeric tolerance added to each side before testing (CRS units).
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.bboxes_adjacent <- function(bb1, bb2, buffer = 0.01) {
  # TRUE unless the boxes are separated along x OR along y. Using `unname()`
  # keeps the result a bare logical scalar rather than a named one.
  unname(
    !(bb1["xmax"] + buffer < bb2["xmin"] |
        bb2["xmax"] + buffer < bb1["xmin"] |
        bb1["ymax"] + buffer < bb2["ymin"] |
        bb2["ymax"] + buffer < bb1["ymin"])
  )
}

#' Element-wise shared-boundary length for paired geometries
#'
#' For each \code{i}, returns the total length of the linear (1-D) parts of the
#' intersection of \code{geom_a[i]} and \code{geom_b[i]} -- i.e. the length of
#' the shared boundary. Point-only contacts and any areal slivers contribute 0.
#'
#' @param geom_a,geom_b Equal-length \pkg{sf} geometry columns (sfc).
#' @return Numeric vector (CRS units), with units stripped.
#' @keywords internal
#' @noRd
.shared_lengths_pairs <- function(geom_a, geom_b) {
  stopifnot(length(geom_a) == length(geom_b))
  vapply(seq_along(geom_a), function(i) {
    inter <- suppressWarnings(
      tryCatch(sf::st_intersection(geom_a[i], geom_b[i]),
               error = function(e) NULL)
    )
    if (is.null(inter) || length(inter) == 0) return(0)
    inter <- inter[!sf::st_is_empty(inter)]
    if (length(inter) == 0) return(0)
    # Reduce to linear components; on a plain LINESTRING this is a no-op, on a
    # GEOMETRYCOLLECTION it isolates the shared edge.
    lines <- suppressWarnings(
      tryCatch(sf::st_collection_extract(inter, "LINESTRING"),
               error = function(e) inter)
    )
    if (length(lines) == 0) return(0)
    as.numeric(sum(suppressWarnings(sf::st_length(lines))))
  }, numeric(1))
}
#' Warn when polygon identifiers are not unique
#' @param ids Identifier vector in vertex order.
#' @param id_col Column name, for the message.
#' @return \code{ids}, invisibly.
#' @keywords internal
#' @noRd
.check_unique_ids <- function(ids, id_col) {
  dup <- unique(ids[!is.na(ids) & duplicated(ids)])
  if (length(dup)) {
    warning(sprintf(paste0(
      "%d duplicated value(s) in `%s`, e.g. %s. Vertices are still built, ",
      "but crosswalk() and joins by `polygon_id` need one row per unit. ",
      "Deduplicate or dissolve the polygons first."),
      length(dup), id_col, paste(utils::head(dup, 3L), collapse = ", ")),
      call. = FALSE)
  }
  invisible(ids)
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
.find_nearest_pair <- function(coords, idx_a, idx_b, method, centroids = NULL) {
  
  switch(method,
         
         # --- Method 1: RANN (fastest; pure coordinate kd-tree lookup) 
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
         
         # --- Method 2: nngeo (sf-native, CRS/geodesic aware) 
         nngeo = {
           .need("nngeo")
           nn <- nngeo::st_nn(
             centroids[idx_a], centroids[idx_b],
             k = 1L, progress = FALSE, returnDist = TRUE
           )
           best <- which.min(unlist(nn$dist))
           list(from = idx_a[best], to = idx_b[nn$nn[[best]]])
         },
         
         # --- Method 3: sf (no extra deps, slowest) 
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

# =============================================================================*
# Bridge edges: distances, penalty and exclusion (package-internal)
#
# Bridges are the edges connect_components() adds to join disconnected
# components (islands, topology gaps). They are flagged at graph creation
# (E(g)$bridge, E(g)$bridge_dist in metres), penalised by distance in
# add_edge_weight(), and can be excluded from Level 3 neighbourhoods.
#
# Distance definition differs by path (recorded as graph attribute
# `bridge_dist_method`):
#   * in-memory build_network(): "boundary" -- shortest polygon-to-polygon
#     distance between the two bridged units;
#   * chunked build_network():  "centroid" -- centroid-to-centroid
#     distance, because the polygons are no longer held in memory.
# Centroid distance is always >= boundary distance for the same pair, so the
# chunked path penalises the same bridge more.
# =============================================================================*

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
#' enough for every edge of a graph with 10^5 or more units. Geographic data use
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
         "build_network() or connect_components().", call. = FALSE)
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
                        "are excluded. Rebuild it with build_network()",
                        "and add_edge_weight() to flag bridges."))
    return(graph)
  }
  idx <- which(br %in% TRUE)
  if (!length(idx)) return(graph)
  .log(verbose, "Excluding %d bridge edge(s) from neighbourhoods.", length(idx))
  igraph::delete_edges(graph, idx)
}

