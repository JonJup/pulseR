#' Determine Core Regions from Fuzzy Membership Degrees
#'
#' @description
#' Identifies *core regions* for each fuzzy class from per-polygon membership
#' degrees. A core region is a set of spatially contiguous polygons whose
#' membership degree for that class meets or exceeds a threshold (`cutoff`).
#'
#' @param x An `sf` object whose features are (multi)polygons carrying one
#'   column of fuzzy membership degrees per class.
#' @param membership_cols Character vector of column names holding the
#'   membership degrees. If `NULL` (default), every numeric column in `x`
#'   (excluding the geometry column) is treated as a membership column.
#' @param cutoff Single numeric in `[0, 1]`. Minimum membership degree for a
#'   polygon to be included in a class's core region. Default `0.5`.
#' @param queen Logical scalar. If `TRUE` (default) neighbours are defined by
#'   queen contiguity (polygons sharing an edge *or* a single vertex). If
#'   `FALSE`, rook contiguity is used (only polygons sharing an edge of
#'   positive length).
#'
#' @return
#' A named `list` with one element per membership column. Each element is an
#' `sf` object containing the polygons that belong to that class's core
#' region(s), with an additional **first** column `core_region_id` (integer)
#' labelling distinct contiguous regions. Classes with no qualifying polygons
#' yield a zero-row `sf` object with the same columns. The `core_region_id`
#' labels are arbitrary integers and are only comparable within a single class.
#'
#' @details
#' For each membership column the function:
#' \enumerate{
#'   \item selects polygons whose membership degree is `>= cutoff`
#'     (missing values are treated as *not* meeting the threshold);
#'   \item builds a sparse spatial neighbour list using the requested
#'     contiguity rule;
#'   \item derives an undirected graph and extracts its connected components,
#'     each of which is one contiguous core region;
#'   \item returns the selected polygons annotated with their region id.
#' }
#'
#' Contiguity is evaluated with DE-9IM relate patterns: `"F***T****"` for queen
#' (boundaries meet at a point *or* a line) and `"F***1****"` for rook
#' (boundaries meet along a line only). The neighbour list is kept sparse so the
#' routine scales to large polygon sets without materialising an
#' \eqn{n \times n} adjacency matrix. Contiguity is only meaningful for
#' polygonal geometries.
#'
#' @examples
#' \dontrun{
#' # 'fuzzy_data' is an sf object with one membership column per class.
#' core_regions <- get_core_regions(
#'   x = fuzzy_data,
#'   membership_cols = c("class_A", "class_B", "class_C"),
#'   cutoff = 0.6
#' )
#'
#' # Core regions for class_A (sf object, 'core_region_id' is the first column):
#' core_regions$class_A
#' }
#'
#' @importFrom sf st_drop_geometry st_relate
#' @importFrom igraph graph_from_adj_list components simplify
#' @export
get_core_regions <- function(x,
                             membership_cols = NULL,
                             cutoff = 0.5,
                             queen = TRUE) {
        
        ## ---- Argument validation -------------------------------------------------
        
        # `x` must be an sf object; we rely on its geometry and sf subsetting rules.
        if (!inherits(x, "sf")) {
                stop("'x' must be an sf object.", call. = FALSE)
        }
        
        # `cutoff` must be a single, non-missing number in the unit interval, because
        # membership degrees are fuzzy weights in [0, 1]. (An NA here would otherwise
        # crash the range test below, which combines NAs with `||`.)
        if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
                stop("'cutoff' must be a single non-missing numeric value.", call. = FALSE)
        }
        if (cutoff < 0 || cutoff > 1) {
                stop("'cutoff' must be between 0 and 1.", call. = FALSE)
        }
        
        # `queen` selects the contiguity rule and must be an unambiguous flag.
        if (!is.logical(queen) || length(queen) != 1L || is.na(queen)) {
                stop("'queen' must be a single logical value (TRUE or FALSE).",
                     call. = FALSE)
        }
        
        # `membership_cols`, if supplied, must be character.
        if (!is.null(membership_cols) && !is.character(membership_cols)) {
                stop("'membership_cols' must be NULL or a character vector.",
                     call. = FALSE)
        }
        
        # Warn before silently overwriting a pre-existing column of the same name.
        if ("core_region_id" %in% names(x)) {
                warning("'x' already contains a 'core_region_id' column; ",
                        "it will be overwritten in the output.", call. = FALSE)
        }
        
        # Resolve the active geometry column *by name* instead of assuming "geometry":
        # sf objects may use "geom", "SHAPE", etc.
        geom_col <- attr(x, "sf_column")
        
        ## ---- Resolve and check membership columns --------------------------------
        
        if (is.null(membership_cols)) {
                # Drop geometry first so the logical test aligns 1:1 with the remaining
                # column names. Indexing names(x) with this logical would misalign by one
                # column (the geometry) and select the wrong columns by recycling.
                attr_tbl <- sf::st_drop_geometry(x)
                is_num <- vapply(attr_tbl, is.numeric, logical(1))
                membership_cols <- names(attr_tbl)[is_num]
                
                if (length(membership_cols) == 0L) {
                        stop("No numeric columns found in 'x' to use as membership columns.",
                             call. = FALSE)
                }
                message("Using all numeric columns as membership columns: ",
                        paste(membership_cols, collapse = ", "))
        } else {
                # De-duplicate so repeated names cannot collapse list elements unexpectedly.
                membership_cols <- unique(membership_cols)
                
                # Every requested column must exist.
                missing_cols <- setdiff(membership_cols, names(x))
                if (length(missing_cols) > 0L) {
                        stop("The following membership columns are not in 'x': ",
                             paste(missing_cols, collapse = ", "), call. = FALSE)
                }
                
                # ... and must be numeric. Testing x[[cc]] directly also rejects the
                # geometry column (an sfc is not numeric).
                not_numeric <- membership_cols[
                        !vapply(membership_cols, function(cc) is.numeric(x[[cc]]), logical(1))
                ]
                if (length(not_numeric) > 0L) {
                        stop("The following membership columns are not numeric: ",
                             paste(not_numeric, collapse = ", "), call. = FALSE)
                }
        }
        
        # Sanity check: fuzzy membership degrees should lie in [0, 1]. Warn rather
        # than stop, so accidental inputs (counts, raw scores) are flagged without
        # blocking legitimate edge cases.
        finite_vals <- unlist(
                lapply(membership_cols, function(cc) {
                        v <- x[[cc]]
                        v[is.finite(v)]
                }),
                use.names = FALSE
        )
        if (length(finite_vals) > 0L &&
            (min(finite_vals) < 0 || max(finite_vals) > 1)) {
                warning("Some membership values fall outside [0, 1]; check that ",
                        "'membership_cols' hold fuzzy membership degrees.", call. = FALSE)
        }
        
        ## ---- Local helpers -------------------------------------------------------
        
        # DE-9IM relate pattern for the chosen contiguity rule.
        #   queen  "F***T****": interiors disjoint, boundaries meet at a point OR line
        #   rook   "F***1****": interiors disjoint, boundaries meet along a line only
        # (The leading "F" excludes self-matches, since a polygon's interior
        #  intersects itself, so no manual diagonal removal is required.)
        contiguity_pattern <- if (queen) "F***T****" else "F***1****"
        
        # Move `core_region_id` to the front, preserve all other attribute columns,
        # and keep the (possibly non-"geometry") geometry column last. Used for both
        # populated and empty results so every list element shares one column order.
        place_id_first <- function(obj) {
                middle <- setdiff(names(obj), c("core_region_id", geom_col))
                obj[, c("core_region_id", middle, geom_col)]
        }
        
        ## ---- Per-class core region detection -------------------------------------
        
        core_regions_list <- vector("list", length(membership_cols))
        names(core_regions_list) <- membership_cols
        
        for (class_col in membership_cols) {
                
                # Polygons meeting the threshold. Coercing NA comparisons to FALSE excludes
                # missing memberships; otherwise NA logical indices would inject all-NA
                # "phantom" rows when used to subset `x`.
                core_mask <- !is.na(x[[class_col]]) & x[[class_col]] >= cutoff
                
                if (!any(core_mask)) {
                        message("No polygons meet the cutoff for class '", class_col, "'.")
                        empty <- x[0, , drop = FALSE]
                        empty$core_region_id <- integer(0)
                        core_regions_list[[class_col]] <- place_id_first(empty)
                        next
                }
                
                core_polys <- x[core_mask, , drop = FALSE]
                
                if (nrow(core_polys) == 1L) {
                        # A lone qualifying polygon is trivially its own region; skip the (for
                        # n == 1, unnecessary) graph construction.
                        core_polys$core_region_id <- 1L
                } else {
                        nb <- sf::st_relate(core_polys,
                                            pattern = contiguity_pattern,
                                            sparse = TRUE)
                        
                        adj_graph <- igraph::simplify(
                                igraph::graph_from_adj_list(unclass(nb), mode = "all")
                        )
                        
                        comp <- igraph::components(adj_graph)
                        membership <- as.integer(comp$membership)
                        
                        # Keep ONLY the largest contiguous region. 
                        areas <- as.numeric(sf::st_area(core_polys))
                        largest_id <- as.integer(names(which.max(tapply(areas, membership, sum))))
                        core_polys  <- core_polys[membership == largest_id, , drop = FALSE]
                        core_polys$core_region_id <- 1L
                }
                
                core_regions_list[[class_col]] <- place_id_first(core_polys)
        }
        
        core_regions_list
}