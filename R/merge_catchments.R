#' Merge adjacent catchments based on environmental similarity
#'
#' Iteratively merges adjacent catchments whose environmental dissimilarity
#' falls below a specified threshold. Catchments are merged in increasing
#' order of dissimilarity within each round, with each catchment participating
#' in at most one merge per round.
#'
#' Environmental variables for a merged catchment are calculated as the
#' arithmetic mean of the corresponding values of the two catchments being
#' merged. Geometries are combined using [sf::st_union()]. Original catchment
#' IDs are retained in the .original_ids list column.
#'
#' The procedure stops when there are no eligible edges, only one catchment
#' remains, the requested target number of catchments is reached, or the
#' maximum number of rounds is exceeded.
#'
#' @param g An igraph graph object. The graph is rebuilt from x at the
#' beginning of each round, so its initial topology is not otherwise used.
#' @param x An sf object containing the catchment geometries and attributes.
#' @param id_col A character string giving the name of the column containing
#' unique catchment IDs.
#' @param env_cols A character vector giving the names of environmental
#' variables used to calculate edge weights and merged-catchment values.
#' @param threshold Numeric value. Only edges with a dissimilarity weight
#' strictly less than this value are eligible for merging.
#' @param max_original Maximum number of original catchments that may be
#' represented by a single merged catchment. Defaults to Inf.
#' @param target_n Optional target number of catchments. If supplied, merging
#' stops once nrow(x) <= target_n, and merges that would reduce the number
#' of catchments below target_n are not performed.
#' @param max_rounds Maximum number of merging rounds to perform. Defaults to
#' Inf.
#'
#' @return An sf object containing the resulting catchments. Merged
#' catchments receive generated IDs of the form "M<round>_<number>".
#' The output contains a .original_ids list column recording the original
#' catchment IDs represented by each catchment and a .merge_round column
#' indicating the round in which a catchment was created (0 for unchanged
#' original catchments). Merged catchments also contain .n_original,
#' giving the number of original catchments represented by the row.
#'
#' @details
#' Before merging, the function checks that x inherits from class sf,
#' that id_col and all env_cols exist, and that the original IDs are
#' unique. IDs are converted to character before processing.
#'
#' In each round, a spatial network is constructed from the current
#' catchments using [polygon_to_network()] and edge weights are calculated
#' using [add_edge_weight()]. Edges with weights below threshold are
#' considered in increasing order of dissimilarity. A catchment can be used
#' in at most one merge during a round.
#'
#' @seealso [polygon_to_network()], [add_edge_weight()], [sf::st_union()]
#'
#' @importFrom sf st_as_sf st_crs st_geometry st_sfc st_union
#' @importFrom dplyr bind_rows mutate tibble union
#' @importFrom igraph ecount E ends
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' result <- merge_catchments(
#' g = g,
#' x = catchments,
#' id_col = "catchment_id",
#' env_cols = c("temperature", "precipitation"),
#' threshold = 0.2,
#' max_original = 5,
#' target_n = 50,
#' max_rounds = 10
#' )
#' }
merge_catchments <- function(
                g,
                x,
                id_col,
                env_cols,
                threshold,
                max_original = Inf,
                target_n = NULL,
                max_rounds = Inf
) {
        
        # ------------------------------------------------------------*
        # Initial checks
        # ------------------------------------------------------------*
        
        stopifnot(inherits(x, "sf"))
        
        if (!id_col %in% names(x)) {
                stop("id_col not found in x")
        }
        
        if (!all(env_cols %in% names(x))) {
                stop("Some env_cols are not found in x")
        }
        
        if (anyDuplicated(x[[id_col]])) {
                stop("Original IDs must be unique.")
        }
        
        # Turn ID column into character to harmonize with outputs of this 
        # function
        x[[id_col]] <- as.character(x[[id_col]])
        geom_col <- attr(x, "sf_column")
        # Keep original IDs as a list column
        x <- dplyr::mutate(x,
                           .original_ids = as.list(.data[[id_col]])
        )
        
        round <- 0
        
        repeat {
                
                round <- round + 1
                
                if (nrow(x) <= 1) {
                        break
                }
                
                if (!is.null(target_n) && nrow(x) <= target_n) {
                        break
                }
                
                if (round > max_rounds) {
                        break
                }
                
                # Build current graph 
                g <- polygon_to_network(x)
                g <- add_edge_weight(g = g, variables = env_cols)
                
                
                
                if (igraph::ecount(g) == 0) {
                        break
                }
                
                # ----------------------------------------------------------
                # Find candidate edges below threshold
                # ----------------------------------------------------------
                
                eligible <- which(igraph::E(g)$weight < threshold)
                
                if (length(eligible) == 0) {
                        break
                }
                
                # Process edges in increasing dissimilarity order
                eligible <- eligible[
                        order(igraph::E(g)$weight[eligible])
                ]
                
                # Catchments already involved in a merge this round
                used <- logical(nrow(x))
                
                # Pairs to merge during this round
                merges <- list()
                
                for (e in eligible) {
                        
                        ends_e <- igraph::ends(g, e, names = FALSE)
                        
                        i <- as.integer(ends_e[1])
                        j <- as.integer(ends_e[2])
                        
                        # Already touched this round?
                        if (used[i] || used[j]) {
                                next
                        }
                        
                        # Respect maximum number of original catchments
                        n_orig_i <- length(x$.original_ids[[i]])
                        n_orig_j <- length(x$.original_ids[[j]])
                        
                        if ((n_orig_i + n_orig_j) > max_original) {
                                next
                        }
                        
                        # Respect target number of rows
                        #
                        # If target_n is supplied, don't perform a merge if it
                        # would take us below the target.
                        if (!is.null(target_n) &&
                            (nrow(x) - length(merges) - 1) < target_n) {
                                next
                        }
                        
                        # Accept this pair
                        used[i] <- TRUE
                        used[j] <- TRUE
                        
                        merges[[length(merges) + 1]] <- c(i, j)
                }
                
                # No merges possible in this round
                if (length(merges) == 0) {
                        break
                }
                
                # ----------------------------------------------------------
                # Construct the next set of catchments
                # ----------------------------------------------------------
                
                merged <- vector("list", length(merges))
                
                merged_indices <- integer(0)
                
                for (k in seq_along(merges)) {
                        
                        ij <- merges[[k]]
                        i <- ij[1]
                        j <- ij[2]
                        
                        # Environmental values = mean of the two current
                        # catchments (using vapply instead of map_dbl)
                        env_mean <- vapply(
                                env_cols,
                                function(v) {
                                        mean(
                                                c(x[[v]][i], x[[v]][j]),
                                                na.rm = TRUE
                                        )
                                },
                                FUN.VALUE = numeric(1)
                        )
                        
                        # Combine original IDs
                        original_ids <- c(
                                x$.original_ids[[i]],
                                x$.original_ids[[j]]
                        )
                        
                        # Geometry union
                        geom <- sf::st_union(
                                sf::st_geometry(x)[c(i, j)]
                        )[[1]]
                        
                        merged[[k]] <- list(
                                i = i,
                                j = j,
                                env = env_mean,
                                original_ids = original_ids,
                                geometry = geom
                        )
                        
                        merged_indices <- c(
                                merged_indices,
                                i,
                                j
                        )
                }
                
                # ----------------------------------------------------------
                # Rows that weren't merged
                # ----------------------------------------------------------
                
                keep <- setdiff(
                        seq_len(nrow(x)),
                        merged_indices
                )
                
                unchanged <- x[keep, , drop = FALSE]
                
                # ----------------------------------------------------------
                # Create merged rows
                # ----------------------------------------------------------
                
                # Using lapply + bind_rows instead of map_dfr
                new_rows_list <- lapply(
                        merged,
                        function(m) {
                                
                                out <- dplyr::tibble(
                                        .merge_round = round,
                                        .n_original = length(m$original_ids),
                                        .original_ids = list(m$original_ids)
                                )
                                
                                for (k in seq_along(env_cols)) {
                                        out[[env_cols[k]]] <- m$env[k]
                                }
                                
                                out
                        }
                )
                
                new_rows <- dplyr::bind_rows(new_rows_list)
                
                # Add geometry separately
                new_geom <- sf::st_sfc(
                        lapply(
                                merged,
                                function(m) m$geometry
                        ),
                        crs = sf::st_crs(x)
                )
                
                # new_rows <- sf::st_sf(
                #         new_rows,
                #         geometry = new_geom
                # )
                new_rows[[geom_col]] <- new_geom
                new_rows <- sf::st_sf(new_rows, sf_column_name = geom_col)
                # Give the newly merged catchments IDs
                new_rows[[id_col]] <- paste0(
                        "M",
                        round,
                        "_",
                        seq_len(nrow(new_rows))
                )
                
                
                if (!".merge_round" %in% names(unchanged)) unchanged$.merge_round <- 0L
                
                att_u <- sf::st_drop_geometry(unchanged)
                att_n <- sf::st_drop_geometry(new_rows)
                cols  <- dplyr::union(names(att_u), names(att_n))
                for (nm in setdiff(cols, names(att_u))) att_u[[nm]] <- NA
                for (nm in setdiff(cols, names(att_n))) att_n[[nm]] <- NA
                
                x <- dplyr::bind_rows(att_u[, cols], att_n[, cols])
                x[[geom_col]] <- c(sf::st_geometry(unchanged), sf::st_geometry(new_rows))
                x <- sf::st_sf(x, sf_column_name = geom_col)
                x <- sf::st_make_valid(x)

        }
        
        # ------------------------------------------------------------
        # Final result
        # ------------------------------------------------------------
        x$.n_original <- lengths(x$.original_ids)
        x
}
