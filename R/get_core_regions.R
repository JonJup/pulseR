#' Determine Core Regions from Fuzzy Membership Degrees
#'
#' @description
#' This function identifies core regions for each class based on fuzzy membership
#' degrees. Core regions are defined as contiguous groups of polygons where the
#' membership degree meets or exceeds a specified threshold.
#'
#' @param x An sf object containing polygons with fuzzy membership degrees.
#' @param membership_cols A character vector of column names containing the
#'   membership degrees for each class. If NULL, all numeric columns except
#'   geometry are used.
#' @param cutoff A numeric value between 0 and 1 specifying the minimum
#'   membership degree threshold for core region inclusion. Default is 0.5.
#' @param queen Logical. If TRUE (default
#', uses queen contiguity (polygons
 #'   sharing edges or vertices are considered neighbors). If FALSE, uses rook
 #'   contiguity (only shared edges).
 #'
 #' @return A named list where each element corresponds to a class (named by
 #'   the membership column). Each element is an sf object containing the core
 #'   region polygons for that class, with an additional column 'core_region_id'
 
 #'   identifying distinct contiguous regions.
 #'
 #' @details
 #' The function performs the following steps for each class:
 #' \enumerate{
 #'   \item Filters polygons where membership degree >= cutoff
 #'   \item Builds a spatial adjacency graph
 #'   \item Identifies connected components (contiguous regions)
 #'   \item Returns polygons with their core region assignment
 #' }
 #'
 #' @examples
 #' \dontrun
 #' # Assuming 'fuzzy_data' is an sf object with membership columns
 #' core_regions <- get_core_regions(
 #'   x = fuzzy_data,
 #'   membership_cols = c("class_A", "class_B", "class_C"),
 #'   cutoff = 0.6
 #' )
 #'
 #' # Access core regions for class_A
 #' core_regions$class_A
 #' }
 #'
 #' @importFrom sf st_touches st_drop_geometry
 #' @importFrom igraph graph_from_adjacency_matrix components
 #' @export
 get_core_regions <- function(x,
                              membership_cols = NULL,
                              cutoff = 0.5,
                              queen = TRUE) {
         
         # Input validation
         if (!inherits(x, "sf")) {
                 stop("Input 'x' must be an sf object.", call. = FALSE)
         }
         
         if (cutoff < 0 || cutoff > 1) {
                 stop("'cutoff' must be between 0 and 1.", call. = FALSE)
         }
         
         # Get membership columns if not specified
         if (is.null(membership_cols)) {
                 numeric_cols <- names(x)[sapply(sf::st_drop_geometry(x), is.numeric)]
                 membership_cols <- numeric_cols
                 message("Using all numeric columns as membership columns: ",
                         paste(membership_cols, collapse = ", "))
         }
         
         # Check that all membership columns exist
         missing_cols <- setdiff(membership_cols, names(x))
         if (length(missing_cols) > 0) {
                 stop("The following membership columns are not in 'x': ",
                      paste(missing_cols, collapse = ", "), call. = FALSE)
         }
         
         # Initialize result list
         core_regions_list <- vector("list", length(membership_cols))
         names(core_regions_list) <- membership_cols
         
         # Process each class
         for (class_col in membership_cols) {
                 
                 # Filter polygons meeting the cutoff threshold
                 core_mask <- x[[class_col]] >= cutoff
                 
                 if (sum(core_mask, na.rm = TRUE) == 0) {
                         message("No polygons meet the cutoff threshold for class '", class_col, "'.")
                         core_regions_list[[class_col]] <- x[0, ]  
                         core_regions_list[[class_col]]$core_region_id <- integer(0)
                         next
                 }
                 
                 core_polys <- x[core_mask, ]
                 
                 # Find contiguous regions using spatial adjacency
                 if (nrow(core_polys) == 1) {
                         # Single polygon is its own region
                         core_polys$core_region_id <- 1L
                 } else {
                         # Build adjacency matrix based on spatial relationships
                         # Queen: st_touches includes shared edges and vertices
                         # Rook: st_relate with pattern "F***T****" for shared edges only
                         if (queen) {
                                 adjacency <- sf::st_touches(core_polys, sparse = FALSE)
                         } else {
                                 # Rook contiguity: only shared edges (not just points)
                                 adjacency <- sf::st_relate(core_polys, core_polys,
                                                            pattern = "F***T****")
                         }
                         
                         # Create graph and find connected components
                         adj_graph <- igraph::graph_from_adjacency_matrix(
                                 adjacency,
                                 mode = "undirected",
                                 diag = FALSE
                         )
                         
                         components_result <- igraph::components(adj_graph)
                         core_polys$core_region_id <- as.integer(components_result$membership)
                 }
                 
                 # Reorder columns to put core_region_id near the front
                 col_order <- c("core_region_id",
                                setdiff(names(core_polys), c("core_region_id", "geometry")),
                                "geometry")
                 col_order <- col_order[col_order %in% names(core_polys)]
                 core_polys <- core_polys[, col_order]
                 
                 core_regions_list[[class_col]] <- core_polys
         }
         
         return(core_regions_list)
}