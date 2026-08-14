# ===========================================================================*
# polygon_to_network()  --  build a spatial contiguity network from polygons
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
#'   components in the network are joined using the package-internal
#'   \code{connect_components()} helper.
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
#'       and edges represent spatial contiguity. The graph is undirected and
#'       simplified (no multi-edges or self-loops).}
#'     \item{polygons}{Returned only when \code{chunked = FALSE}: the cleaned
#'       \pkg{sf} polygon object actually used to build the network (after
#'       validity repair and empty-geometry removal).}
#'     \item{edge_data}{A data frame with integer columns \code{from} and
#'       \code{to} giving connected polygon pairs as 1-based global indices.}
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
#' tolerance. Cross-chunk neighbours use the matching \pkg{sf} predicate so the
#' two phases agree: \code{sf::st_touches()} for Queen contiguity, and
#' \code{sf::st_relate()} with the DE-9IM pattern \code{"F***1****"}
#' (shared 1-dimensional boundary) for Rook contiguity. The resulting graph is
#' undirected, since spatial contiguity is symmetric.
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
#' identifiers are stored as the \code{polygon_id} vertex attribute.
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
#' network <- polygon_to_network(polys)
#' plot(network$graph)
#'
#' # Queen contiguity, no auto-connection of islands
#' network <- polygon_to_network(polys, queen = TRUE, connect_islands = FALSE)
#'
#' # Only count neighbours sharing at least 50 m of boundary (projected CRS)
#' network <- polygon_to_network(polys, min_shared_length = 50)
#'
#' # --- Chunked mode across multiple GeoParquet files ---
#' files   <- list.files("data/chunks", pattern = "\\.parquet$", full.names = TRUE)
#' network <- polygon_to_network(files, chunked = TRUE,
#'                               return_combined_polygons = TRUE)
#'
#' # Chunked mode with paired environmental control files
#' env_files <- list.files("data/env", pattern = "\\.parquet$", full.names = TRUE)
#' network   <- polygon_to_network(files, chunked = TRUE, control_paths = env_files)
#' }
#'
#' @seealso \code{\link[spdep]{poly2nb}} for neighbour detection,
#'   \code{\link[igraph]{make_empty_graph}} for graph construction.
#'
#' @importFrom spdep poly2nb
#' @importFrom sf st_make_valid st_is_empty st_centroid st_geometry st_bbox st_touches st_relate st_length st_intersection st_collection_extract st_crs st_is_longlat
#' @importFrom sfarrow st_read_parquet
#' @importFrom arrow read_parquet
#' @importFrom igraph make_empty_graph add_edges set_vertex_attr vcount ecount components simplify
#'
#' @export
polygon_to_network <- function(polygons,
                               chunked                  = FALSE,
                               min_shared_length        = 0,
                               connect_islands          = TRUE,
                               control_paths            = NULL,
                               return_combined_polygons = FALSE,
                               queen                    = FALSE,
                               snap                     = 1e-6,
                               bbox_buffer              = 0.01,
                               verbose                  = TRUE) {

  # ---- Hard dependency checks ---- *
  # spdep is required for both modes; the heavy I/O packages (sfarrow, arrow)
  # are only needed for chunked mode and are checked there so they can live in
  # Suggests rather than Imports.
  .need("spdep")
  .need("sf")
  .need("igraph")

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
  edge_df <- nb_to_edgelist(nb)

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
    g <- igraph::set_vertex_attr(g, "polygon_id", value = polygons[[id_col]])
  }

  comp <- igraph::components(g)$no
  .log(verbose, "Network created: %d nodes, %d edges, %d components.",
       igraph::vcount(g), igraph::ecount(g), comp)

  if (comp != 1L && connect_islands) {
    g <- connect_components(g, polygons = polygons)
  }
  
  igraph::E(g)$eid <- seq_len(ecount(g))
  
  list(
    graph     = g,
    polygons  = polygons,   # cleaned object, aligned with vertex ids
    edge_data = edge_df
  )
}


# ===========================================================================*
# Internal: chunked, out-of-core implementation
# ===========================================================================*
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
      edges_k <- nb_to_edgelist(nb_k)

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
      if (bboxes_adjacent(chunk_meta[[i]]$bbox, chunk_meta[[j]]$bbox, buffer = bbox_buffer))
        j_adj <- c(j_adj, j)
    }
    if (length(j_adj) == 0L)
      next
    
    # Load chunk i ONCE, reuse across all its partners.
    polys_i <- .load_chunk(chunk_meta[[i]]$path, chunk_meta[[i]]$control_path)
    
    for (j in j_adj) {
      .log(verbose, "Checking cross-chunk edges: %d x %d", i, j)
      polys_j <- .load_chunk(chunk_meta[[j]]$path, chunk_meta[[j]]$control_path)
      
      sgbp <- if (queen) {
        sf::st_touches(polys_i, polys_j)
      } else {
        sf::st_relate(polys_i, polys_j, pattern = "F***1****")
      }
      
      sgbp   <- unclass(sgbp)
      counts <- lengths(sgbp)
      pairs  <- if (sum(counts) > 0L) {
        data.frame(row = rep.int(seq_along(sgbp), counts),
                   col = as.integer(unlist(sgbp, use.names = FALSE)))
      } else
        NULL
      
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
      rm(polys_j, sgbp, pairs)
    }
    rm(polys_i)
    gc(verbose = FALSE)
  }
  # for (i in seq_along(chunk_meta)) {
  #   if (chunk_meta[[i]]$n == 0L) next
  #   for (j in seq_along(chunk_meta)) {
  #     if (j <= i || chunk_meta[[j]]$n == 0L) next
  #     if (!bboxes_adjacent(chunk_meta[[i]]$bbox, chunk_meta[[j]]$bbox,
  #                          buffer = bbox_buffer)) next
  # 
  #     .log(verbose, "Checking cross-chunk edges: %d x %d", i, j)
  # 
  #     # Reload through the SAME loader as phase 1 to guarantee identical row
  #     # order, otherwise global_idx[row] would point at the wrong polygon.
  #     polys_i <- .load_chunk(chunk_meta[[i]]$path, chunk_meta[[i]]$control_path)
  #     polys_j <- .load_chunk(chunk_meta[[j]]$path, chunk_meta[[j]]$control_path)
  # 
  #     # Use the predicate that matches the requested contiguity rule so that
  #     # cross-chunk edges obey the same definition as within-chunk edges.
  #     sgbp <- if (queen) {
  #       sf::st_touches(polys_i, polys_j)
  #     } else {
  #       # DE-9IM "F***1****": interiors disjoint, boundaries meet in a line.
  #       sf::st_relate(polys_i, polys_j, pattern = "F***1****")
  #     }
  # 
  #     # Flatten the sparse predicate into local (row, col) candidate pairs.
  #     # pairs <- do.call(rbind, lapply(seq_along(sgbp), function(r) {
  #     #   cols <- sgbp[[r]]
  #     #   if (length(cols) == 0L) return(NULL)
  #     #   data.frame(row = r, col = cols)
  #     # }))
  #     sgbp   <- unclass(sgbp)          # plain list ops, no sgbp dispatch
  #     counts <- lengths(sgbp)
  #     pairs  <- if (sum(counts) > 0L) {
  #       data.frame(
  #         row = rep.int(seq_along(sgbp), counts),
  #         col = as.integer(unlist(sgbp, use.names = FALSE))
  #       )
  #     } else NULL
  # 
  #     if (!is.null(pairs) && nrow(pairs) > 0) {
  #       if (min_shared_length > 0) {
  #         keep  <- .shared_lengths_pairs(sf::st_geometry(polys_i)[pairs$row],
  #                                        sf::st_geometry(polys_j)[pairs$col]) >=
  #                  min_shared_length
  #         pairs <- pairs[keep, , drop = FALSE]
  #       }
  #       if (nrow(pairs) > 0) {
  #         all_edges[[length(all_edges) + 1L]] <- data.frame(
  #           from = chunk_meta[[i]]$global_idx[pairs$row],
  #           to   = chunk_meta[[j]]$global_idx[pairs$col]
  #         )
  #       }
  #     }
  # 
  #     rm(polys_i, polys_j, sgbp, pairs); gc(verbose = FALSE)
  #   }
  # }

  # ---- Phase 3: assemble graph ---- *
  all_edges <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, all_edges)
  edge_df   <- if (length(all_edges) > 0) {
    do.call(rbind, all_edges)
  } else {
    data.frame(from = integer(), to = integer())
  }

  g <- .build_graph(edge_df, n_total = n_total)

  # Only attach ids if every node has one (a complete vector); a partial set is
  # more misleading than none.
  if (saw_id && length(all_ids) == n_total && !all(is.na(all_ids))) {
    g <- igraph::set_vertex_attr(g, "polygon_id", value = all_ids)
  }

  comp <- igraph::components(g)$no
  .log(verbose, "Combined network: %d nodes, %d edges, %d components.",
       igraph::vcount(g), igraph::ecount(g), comp)

  if (connect_islands && comp > 1L) {
    # do.call(c, ...) concatenates the per-chunk sfc objects in global order;
    # NULL entries (chunks with zero nodes) are dropped harmlessly.
    all_centroids <- do.call(c, Filter(Negate(is.null), centroids_list))
    g <- connect_components(g, centroids = all_centroids)
  }

  igraph::E(g)$eid <- seq_len(ecount(g))
  
  out <- list(
    graph      = g,
    edge_data  = edge_df,
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

#' Convert an spdep neighbour list to an undirected edgelist
#'
#' @param nb An object of class \code{nb} from \code{spdep::poly2nb()}.
#' @return A data frame with integer columns \code{from} and \code{to}, with
#'   \code{from < to} so each undirected edge appears exactly once.
#' @keywords internal
#' @noRd
# nb_to_edgelist <- function(nb) {
#   edges <- lapply(seq_along(nb), function(i) {
#     neighbors <- nb[[i]]
#     # spdep encodes "no neighbours" as a single 0L.
#     if (length(neighbors) > 0 && neighbors[1] != 0) {
#       # Keep only neighbours with a higher index => each edge once, from < to.
#       valid_neighbors <- neighbors[neighbors > i]
#       if (length(valid_neighbors) > 0) {
#         return(data.frame(from = as.integer(i),
#                           to   = as.integer(valid_neighbors)))
#       }
#     }
#     NULL
#   })
#   edges <- Filter(Negate(is.null), edges)
#   if (length(edges) == 0) {
#     return(data.frame(from = integer(), to = integer()))
#   }
#   do.call(rbind, edges)
# }
nb_to_edgelist <- function(nb) {
  n <- length(nb)
  if (n == 0L) return(data.frame(from = integer(), to = integer()))
  
  # Keep only higher-indexed neighbours so each undirected edge appears once
  # (from < to). spdep encodes "no neighbours" as a single 0L.
  to_list <- lapply(seq_len(n), function(i) {
    nbr <- nb[[i]]
    if (length(nbr) && nbr[1L] != 0L) nbr[nbr > i] else integer(0)
  })
  
  counts <- lengths(to_list)
  if (sum(counts) == 0L) return(data.frame(from = integer(), to = integer()))
  
  data.frame(
    from = rep.int(seq_len(n), counts),
    to   = as.integer(unlist(to_list, use.names = FALSE))
  )
}
#' Do two bounding boxes touch or overlap (with a buffer)?
#'
#' @param bb1,bb2 Bounding boxes from \code{sf::st_bbox()}.
#' @param buffer Numeric tolerance added to each side before testing (CRS units).
#' @return Logical scalar.
#' @keywords internal
#' @noRd
bboxes_adjacent <- function(bb1, bb2, buffer = 0.01) {
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