#' Add Environmental Dissimilarity Weights to Graph Edges
#'
#' Calculates edge weights based on Euclidean distances between neighbouring
#' nodes using (optionally) scaled environmental variables, and returns a new
#' \pkg{igraph} object with weighted edges. Supports two modes: an in-memory
#' mode where environmental variables are taken directly from \code{g$polygons},
#' and a chunked mode where they are read from one or more parquet files and
#' stacked in input order.
#'
#' @param g A list containing:
#'   \describe{
#'     \item{graph}{An \code{igraph} object representing the network structure.
#'       Vertices should be named (the names are propagated to the output).}
#'     \item{polygons}{Required when \code{chunked = FALSE}. An \pkg{sf} object
#'       containing environmental variables for each node. Row order must
#'       correspond to graph vertex order (see \strong{Node ordering}).}
#'   }
#' @param chunked Logical scalar. If \code{FALSE} (default), environmental
#'   variables are taken from \code{g$polygons}. If \code{TRUE}, they are read
#'   from the parquet files supplied via \code{files} and concatenated in input
#'   order.
#' @param files Character vector of parquet file paths providing environmental
#'   variables. Required when \code{chunked = TRUE}, ignored otherwise. The
#'   concatenated row order must match the graph vertex order.
#' @param id_col Character vector of ID column names to exclude before
#'   weighting. Only used when \code{variables} is \code{NULL}. Default
#'   \code{NULL}.
#' @param variables Character vector of column names to use for environmental
#'   dissimilarity. If \code{NULL} (default), all non-geometry, non-\code{id_col}
#'   columns are used. If supplied, \code{id_col} is ignored. In chunked mode,
#'   supplying \code{variables} also restricts the columns read from disk,
#'   which can substantially reduce memory use on wide parquet files.
#' @param scale Logical scalar. Standardize variables (subtract mean, divide by
#'   SD) before computing distances. Default \code{TRUE}. Set to \code{FALSE}
#'   when variables are already on a comparable scale.
#'
#' @return An undirected \code{igraph} object whose edge \code{weight} attribute
#'   is the Euclidean distance between neighbouring nodes in (optionally scaled)
#'   environmental space. Vertex names from \code{g$graph} are preserved, all
#'   input vertices are retained (including isolated ones), and each undirected
#'   edge appears exactly once (no multi-edges).
#'
#' @section Node ordering:
#' This function assumes that the i-th row of the environmental data corresponds
#' to the i-th vertex of \code{g$graph} (i.e. the order returned by
#' \code{igraph::as_adjacency_matrix()}). This contract is \emph{not} verified
#' by content, only by count: the function errors if the number of environmental
#' rows differs from the number of vertices, but it cannot detect a same-length
#' \emph{re-ordering}. Ensure \code{g$polygons} / the stacked parquet rows are
#' built in the same order as the graph vertices.
#'
#' @details
#' Both modes share the following pipeline:
#' \enumerate{
#'   \item Extract the adjacency structure from \code{g$graph}.
#'   \item Convert it to an \pkg{spdep} neighbours list via
#'     \code{\link[spdep]{mat2listw}}.
#'   \item Acquire environmental variables (from \code{g$polygons} or by
#'     reading and stacking \code{files}).
#'   \item Coerce all variables to numeric and, if \code{scale = TRUE},
#'     standardize them. Zero-variance and non-numeric columns are handled
#'     defensively (see below).
#'   \item Compute neighbour Euclidean distances via
#'     \code{\link[spdep]{nbcosts}}.
#'   \item Build an undirected \pkg{igraph} object with these distances as
#'     edge weights, preserving vertex names and de-duplicating reciprocal
#'     edges.
#' }
#'
#' Variable selection precedence: if \code{variables} is supplied, it takes
#' precedence and \code{id_col} is ignored; otherwise all columns except those
#' in \code{id_col} (and common geometry-like columns such as \code{Shape},
#' \code{Shape_Length}, \code{Shape_Area}, \code{geom}, \code{geometry}) are
#' used.
#'
#' Data hygiene: columns that are entirely non-numeric after coercion trigger an
#' informative error (they are almost always a mis-classified ID/label column);
#' zero-variance columns are dropped with a warning (they cannot be scaled and
#' add nothing to a Euclidean distance); remaining missing values trigger a
#' warning because they propagate to \code{NA} edge weights, which will break a
#' downstream minimum spanning tree / SKATER step.
#'
#' Isolated vertices (degree 0) are retained in the output as isolated nodes but
#' contribute no edges; a warning is emitted when any are present.
#'
#' @seealso
#' \code{\link[spdep]{nbcosts}} for neighbour cost calculation,
#' \code{\link[igraph]{graph_from_data_frame}} for graph construction.
#'
#' @examples
#' \dontrun{
#' # --- In-memory ---
#' g_weighted        <- add_edge_weight(g)
#' g_weighted_subset <- add_edge_weight(g, variables = c("temp", "precip"))
#'
#' # --- Chunked across parquet files ---
#' files <- list.files("data/env", pattern = "\\.parquet$", full.names = TRUE)
#' g_weighted <- add_edge_weight(
#'   g,
#'   chunked = TRUE,
#'   files   = files,
#'   id_col  = c("ID", "country")
#' )
#' }
#'
#' @importFrom igraph as_adjacency_matrix graph_from_data_frame
#' @importFrom spdep mat2listw nbcosts
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom tidyselect all_of any_of
#' @importFrom arrow read_parquet
#' @importFrom data.table rbindlist
#' @importFrom stats sd
#'
#' @export
add_edge_weight <- function(g,
                            chunked   = FALSE,
                            files     = NULL,
                            id_col    = NULL,
                            variables = NULL,
                            scale     = TRUE) {
        
        # ---- Input validation -----------------------------------------------------
        # Validate logical scalars first so later branching is well defined.
        if (!is.logical(chunked) || length(chunked) != 1L || is.na(chunked)) {
                stop("`chunked` must be a single non-NA logical value.")
        }
        if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
                stop("`scale` must be a single non-NA logical value.")
        }
        
        # `g` must be a list carrying an igraph in `$graph`; otherwise the adjacency
        # extraction below would fail with an opaque internal error.
        if (!is.list(g)) {
                stop("`g` must be a list with a `graph` element (an igraph object).")
        }
        if (!inherits(g$graph, "igraph")) {
                stop("`g$graph` must be an igraph object.")
        }
        
        # Optional selection arguments must be character vectors if supplied.
        if (!is.null(variables) && !is.character(variables)) {
                stop("`variables` must be a character vector or NULL.")
        }
        if (!is.null(id_col) && !is.character(id_col)) {
                stop("`id_col` must be a character vector or NULL.")
        }
        
        # Mode-specific requirements.
        if (chunked) {
                if (!is.character(files) || length(files) == 0L) {
                        stop("`files` must be a non-empty character vector of parquet paths ",
                             "when `chunked = TRUE`.")
                }
                missing_files <- files[!file.exists(files)]
                if (length(missing_files)) {
                        stop("Parquet file(s) not found: ",
                             paste(missing_files, collapse = ", "), ".")
                }
        } else {
                if (is.null(g$polygons)) {
                        stop("`g$polygons` must be provided when `chunked = FALSE`.")
                }
                if (!is.null(files)) {
                        warning("`files` is ignored when `chunked = FALSE`.")
                }
        }
        
        # ---- Neighbour list from graph (shared by both modes) ---------------------
        # Sparse adjacency keeps memory bounded on large continental graphs.
        # `names = TRUE` carries the vertex names through, which we use both to map
        # neighbour indices back to IDs and to preserve all vertices in the output.
        adj_matrix <- igraph::as_adjacency_matrix(
                graph  = g$graph,
                type   = "both",
                names  = TRUE,
                sparse = TRUE
        )
        
        # Vertex names define both the row order of the environmental data and the
        # labelling of the output graph. Fall back to positional labels if unnamed.
        vertex_names <- rownames(adj_matrix)
        if (is.null(vertex_names)) {
                vertex_names <- as.character(seq_len(nrow(adj_matrix)))
        }
        n_vertices <- length(vertex_names)
        
        # style = "B" (binary) is irrelevant here -- we only consume the neighbour
        # structure ($neighbours), not the spatial weights -- but it avoids spdep
        # row-standardising and is the clearest declaration of intent.
        listw_object <- spdep::mat2listw(adj_matrix, style = "B")
        nb_object    <- listw_object$neighbours
        
        # ---- Acquire environmental variables --------------------------------------
        env_data <- if (chunked) {
                .env_from_parquet(files, variables = variables, id_col = id_col)
        } else {
                .env_from_polygons(g$polygons, variables = variables, id_col = id_col)
        }
        
        # Alignment guard: the i-th environmental row must be the i-th vertex. We can
        # only check the count, not the ordering -- see the "Node ordering" section.
        if (nrow(env_data) != n_vertices) {
                stop(sprintf(
                        paste0("Number of environmental rows (%d) does not match the number of ",
                               "graph vertices (%d). Rows must be in graph-vertex order."),
                        nrow(env_data), n_vertices
                ))
        }
        if (ncol(env_data) == 0L) {
                stop("No environmental variables left after column selection.")
        }
        
        # ---- Coerce to numeric, validate, optionally scale ------------------------
        # A shape-stable numeric matrix (unlike apply(..., as.numeric), which drops
        # to a vector for single-row inputs) is what scale() and nbcosts() expect.
        env_matrix <- .coerce_numeric_matrix(env_data)
        
        # Columns that are entirely NA after coercion were non-numeric to begin with
        # (e.g. a stray text/ID column). Report them by name rather than silently
        # producing NA distances.
        all_na_cols <- colnames(env_matrix)[colSums(!is.na(env_matrix)) == 0L]
        if (length(all_na_cols)) {
                stop("These columns are non-numeric (all-NA after coercion) and cannot be ",
                     "used as environmental variables: ",
                     paste(all_na_cols, collapse = ", "),
                     ". Exclude them via `id_col` or restrict `variables`.")
        }
        
        # Constant columns break standardization (SD = 0 -> NaN) and add nothing to a
        # Euclidean distance, so drop them with a warning.
        col_sd     <- apply(env_matrix, 2L, stats::sd, na.rm = TRUE)
        const_cols <- colnames(env_matrix)[is.na(col_sd) | col_sd == 0]
        if (length(const_cols)) {
                warning("Dropping zero-variance variable(s): ",
                        paste(const_cols, collapse = ", "), ".")
                env_matrix <- env_matrix[, !(colnames(env_matrix) %in% const_cols),
                                         drop = FALSE]
        }
        if (ncol(env_matrix) == 0L) {
                stop("No usable (numeric, non-constant) environmental variables remain.")
        }
        
        # Partial missingness propagates to NA edge weights, which silently breaks a
        # downstream MST/SKATER step -- warn loudly so the user imputes/removes first.
        if (anyNA(env_matrix)) {
                warning("Environmental data contain missing values; affected edges will ",
                        "receive NA weights. Impute or remove NAs before clustering.")
        }
        
        if (isTRUE(scale)) {
                # base::scale() is namespaced explicitly because the `scale` argument
                # shadows the function name.
                env_matrix <- base::scale(env_matrix)
        }
        
        # ---- Neighbour Euclidean costs --------------------------------------------
        # Identify isolated vertices: spdep encodes a no-neighbour region as a single
        # element equal to 0L. These are excluded from edge construction below but
        # retained as isolated vertices in the output.
        card     <- lengths(nb_object)
        sentinel <- card == 1L &
                vapply(nb_object, function(x) isTRUE(x[[1L]] == 0L), logical(1))
        isolated <- card == 0L | sentinel
        if (any(isolated)) {
                warning(sum(isolated), " vertex/vertices have no neighbours; they are ",
                        "retained as isolated nodes in the output but contribute no edges.")
        }
        
        netDist <- spdep::nbcosts(
                nb     = nb_object,
                data   = env_matrix,
                method = "euclidean"
        )
        
        # ---- Assemble weighted graph ----------------------------------------------
        # Build the edge list only from non-isolated nodes. For these, nb_object[[i]]
        # holds real (>= 1) neighbour indices and netDist[[i]] holds one cost per
        # neighbour in the same order, so concatenation stays aligned.
        kept     <- which(!isolated)
        from_idx <- rep.int(kept, lengths(nb_object[kept]))
        to_idx   <- unlist(nb_object[kept], use.names = FALSE)
        weight   <- unlist(netDist[kept],   use.names = FALSE)
        
        # Defensive alignment guard: if these lengths ever diverge, fail loudly
        # rather than silently pairing the wrong cost with the wrong edge.
        stopifnot(
                length(from_idx) == length(to_idx),
                length(to_idx)   == length(weight)
        )
        
        # A symmetric neighbour list stores each undirected edge twice (i->j and
        # j->i). Keep only the i < j copy so the output has one edge per pair (no
        # self-loops, no multi-edges).
        upper <- from_idx < to_idx
        
        edges_df <- data.frame(
                from   = vertex_names[from_idx[upper]],
                to     = vertex_names[to_idx[upper]],
                weight = weight[upper],
                stringsAsFactors = FALSE
        )
        
        # Passing an explicit vertex table preserves vertex names AND retains every
        # input vertex (including isolated ones, which never appear in `edges_df`).
        igraph::graph_from_data_frame(
                d        = edges_df,
                directed = FALSE,
                vertices = data.frame(name = vertex_names, stringsAsFactors = FALSE)
        )
}


# =============================================================================
# Internal helpers
# =============================================================================

#' Coerce a data frame of environmental variables to a numeric matrix
#'
#' Column-wise coercion that is shape-stable for single-row inputs and routes
#' factors through their labels (not their integer codes, as
#' \code{data.matrix()} would). Non-numeric columns become all-\code{NA},
#' which the caller detects and reports.
#'
#' @param df A data frame of environmental variables.
#' @return A numeric matrix with the same dimensions and column names as
#'   \code{df}.
#' @keywords internal
#' @noRd
.coerce_numeric_matrix <- function(df) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)
        coerced <- lapply(df, function(col) {
                if (is.numeric(col)) return(as.double(col))
                if (is.factor(col))  col <- as.character(col)
                suppressWarnings(as.numeric(col)) # genuine non-numerics -> NA (caught later)
        })
        matrix(
                unlist(coerced, use.names = FALSE),
                nrow     = nrow(df),
                ncol     = ncol(df),
                dimnames = list(NULL, names(df))
        )
}

#' Extract environmental data frame from an sf polygon object
#'
#' @param polygons An \pkg{sf} object.
#' @param variables Character vector of variable columns to keep, or
#'   \code{NULL} to keep everything except \code{id_col} and geometry-like
#'   columns.
#' @param id_col Character vector of ID columns to drop when \code{variables}
#'   is \code{NULL}.
#' @return A data frame of environmental variables (geometry dropped).
#' @keywords internal
#' @noRd
.env_from_polygons <- function(polygons, variables, id_col) {
        # Common ArcGIS-derived geometry/measurement columns to exclude alongside
        # the active sf geometry.
        geom_like <- c("Shape", "SHAPE", "Shape_Length", "Shape_Area",
                       "geom", "geometry")
        out <- if (!is.null(variables)) {
                dplyr::select(polygons, tidyselect::all_of(variables)) # errors if missing
        } else {
                dplyr::select(polygons, !tidyselect::any_of(c(id_col, geom_like)))
        }
        sf::st_drop_geometry(out)
}

#' Read and concatenate environmental data from parquet files
#'
#' @param files Character vector of parquet file paths.
#' @param variables Character vector of variable columns to keep, or
#'   \code{NULL} to keep everything except \code{id_col} and geometry-like
#'   columns.
#' @param id_col Character vector of ID columns to drop when \code{variables}
#'   is \code{NULL}.
#' @return A data frame of environmental variables, row-stacked across files in
#'   input order.
#' @keywords internal
#' @noRd
.env_from_parquet <- function(files, variables, id_col) {
        geom_like <- c("Shape", "SHAPE", "Shape_Length", "Shape_Area",
                       "geom", "geometry")
        
        chunks <- vector("list", length(files))
        for (i in seq_along(files)) {
                if (!is.null(variables)) {
                        # Read only the requested columns -- a large memory saving on wide,
                        # continental-scale parquet files. Errors clearly if a column is absent.
                        chunks[[i]] <- arrow::read_parquet(
                                files[[i]],
                                col_select = tidyselect::all_of(variables)
                        )
                } else {
                        d <- arrow::read_parquet(files[[i]])
                        chunks[[i]] <- dplyr::select(d, !tidyselect::any_of(c(id_col, geom_like)))
                }
        }
        
        # use.names = TRUE aligns columns by name across files; fill = FALSE (default)
        # makes a column-set mismatch an explicit error rather than silently padding
        # with NA.
        as.data.frame(data.table::rbindlist(chunks, use.names = TRUE))
}