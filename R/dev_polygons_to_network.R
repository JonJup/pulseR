#' Create a Spatial Network from Polygon Geometries
#'
#' This function constructs an igraph network object from polygon geometries based on
#' spatial contiguity. Polygons that share boundaries become connected nodes in the network.
#' The function uses the spdep package for efficient neighbor detection and can optionally
#' connect disconnected components (islands) using spatial proximity.
#'
#' @param polygons An sf object containing polygon geometries. Each polygon will become
#'   a node in the resulting network. If the object contains an 'id' column, these values
#'   will be preserved as vertex attributes in the graph.
#' @param min_shared_length Numeric. Minimum length of shared boundary required for two
#'   polygons to be considered neighbors. Currently not implemented in the function body.
#'   Default is 0 (all touching polygons are connected regardless of shared boundary length).
#' @param connect_islands Logical. If TRUE, disconnected components in the network will be
#'   connected using the connect_components() helper function. This is useful for ensuring
#'   a fully connected network even when some polygons are spatially isolated. Default is TRUE.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{graph}{An igraph object representing the spatial network. Nodes correspond to
#'       polygons and edges represent spatial contiguity.}
#'     \item{polygons}{The original sf polygon object passed to the function.}
#'     \item{edge_data}{A data frame with columns 'from' and 'to' indicating connected
#'       polygon pairs (using 1-based indices).}
#'   }
#'
#' @details
#' The function uses Queen contiguity (polygons sharing any boundary point are neighbors)
#' implemented via spdep::poly2nb(). A small snap tolerance (1e-6) is applied to handle
#' minor topological inconsistencies in the polygon geometries.
#'
#' The resulting graph is undirected, as spatial contiguity is inherently symmetric.
#' If the input polygons contain an 'id' column, these identifiers are stored as the
#' 'polygon_id' vertex attribute in the graph.
#'
#' When connect_islands = TRUE and disconnected components are detected, the function
#' calls connect_components() to add edges between components (implementation not shown here).
#'
#' @note The min_shared_length parameter is defined but not currently used in filtering
#' neighbors. Future implementations may use this to require a minimum boundary length
#' before considering polygons as neighbors.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Create a simple polygon dataset
#' polys <- st_read("my_polygons.shp")
#'
#' # Create network with default settings
#' network <- dev_polygon_to_network(polys)
#'
#' # Access the graph
#' plot(network$graph)
#'
#' # Create network without connecting islands
#' network_islands <- dev_polygon_to_network(polys, connect_islands = FALSE)
#' }
#'
#' @seealso \code{\link[spdep]{poly2nb}} for neighbor detection,
#'   \code{\link[igraph]{graph_from_data_frame}} for graph construction
#'
#' @importFrom spdep poly2nb n.comp.nb nb2mat
#' @importFrom igraph graph_from_data_frame vcount ecount components V E
#' @importFrom dplyr mutate filter select
#' @importFrom tidyr pivot_longer
#' @importFrom igraph "E<-" "V<-"
#'
#' @export
dev_polygon_to_network <- function(polygons, min_shared_length = 0, connect_islands = TRUE) {
        
        # Dependency check: Ensure spdep package is available
        # spdep is required for spatial neighbor detection
        if (!requireNamespace("spdep", quietly = TRUE)) {
                stop("Please install spdep: install.packages('spdep')")
        }
        
        # Get the number of polygons to process
        n <- nrow(polygons)
        cat(sprintf("Processing %d polygons using spdep...\n", n))
        
        # Step 1: Find spatial neighbors using Queen contiguity
        # Queen contiguity: polygons sharing any boundary point (vertex or edge) are neighbors
        # snap parameter adds tolerance for minor geometric inconsistencies
        nb <-
                suppressWarnings(
                        spdep::poly2nb(polygons, queen = TRUE, snap = 1e-6)
                )
        
        # Check the number of disconnected components in the neighbor structure
        # This helps identify spatially isolated groups of polygons
        comp <- spdep::n.comp.nb(nb)
        
        # Step 2: Convert neighbor list to edge data frame DIRECTLY (Avoids Matrix)
        # spdep nb objects are lists. We can parse them directly.
        edge_list <- lapply(1:length(nb), function(i) {
                neighbors <- nb[[i]]
                # Check if neighbors exist and aren't 0 (spdep's code for no neighbor)
                if (length(neighbors) > 0 && neighbors[1] != 0) {
                        # Return edges only where 'from' < 'to' to make it undirected/unique
                        valid_neighbors <- neighbors[neighbors > i]
                        if (length(valid_neighbors) > 0) {
                                return(data.frame(from = i, to = valid_neighbors))
                        }
                }
                return(NULL)
        })
        
        # Combine into one lightweight data frame
        edge_df <- do.call(rbind, edge_list)
        
        # Step 3: Create igraph object from edge data frame
        # directed = FALSE because spatial contiguity is symmetric
        # vertices = 1:n ensures all polygons are included, even isolated ones
        g <- igraph::graph_from_data_frame(
                edge_df,
                directed = FALSE,
                vertices = data.frame(name = 1:n)
        )
        
        # Step 4: Add edge attributes if available
        # Check if shared boundary length was calculated
        # (Note: current implementation doesn't calculate this)
        if ("shared_length" %in% names(edge_df)) {
                E(g)$shared_length <- edge_df$shared_length
        }
        
        # Step 5: Add vertex attributes from original polygon data
        # Preserve polygon IDs if they exist in the input data
        if ("id" %in% names(polygons)) {
                V(g)$polygon_id <- polygons$id
        }
        
        # Report network statistics
        cat(sprintf(
                "Network created: %d nodes, %d edges, %d components\n",
                vcount(g),      # Number of vertices
                ecount(g),      # Number of edges
                components(g)$no  # Number of disconnected components
        ))
        
        # Step 6: Handle disconnected components (islands)
        # If network has multiple components and user wants them connected
        if (components(g)$no != 1 & connect_islands) {
                # Call helper function to add edges between components
                g <- connect_components(g, polygons)
        }
        
        # Step 7: Return results as a structured list
        # This allows access to multiple related objects
        return(list(
                graph = g,           # The igraph network object
                polygons = polygons, # Original spatial polygons
                edge_data = edge_df  # Data frame of edges
        ))
}
