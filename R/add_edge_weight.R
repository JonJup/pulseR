#' Add Environmental Dissimilarity Weights to Graph Edges
#'
#' Calculates edge weights based on Euclidean distances between neighboring 
#' nodes using scaled environmental variables, and returns a new igraph object 
#' with weighted edges.
#'
#' @param g A list containing:
#'   \describe{
#'     \item{graph}{An \code{igraph} object representing the network structure.}
#'     \item{polygons}{An \code{sf} object containing environmental variables 
#'       for each node. Node order must correspond to graph vertices.}
#'   }
#' @param variables Character vector of column names to use for calculating 
#'   environmental dissimilarity. If \code{NULL} (default), all non-geometry 
#'   columns in \code{g$polygons} are used.
#' @param id_col Character vector of column names of ID columns. If \code{NULL} (default), 
#'   it is assumed that none exist in \code{g$polygons}. These columns are removed prior 
#'   to weighting. This argument is commonly used in conjunction with \code{variables = NULL}.
#'   
#' @return An \code{igraph} object with edge weights representing Euclidean 
#'   distances between neighboring nodes in scaled environmental space.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the adjacency structure from the input graph
#'   \item Converts to a \code{spdep} neighbors object
#'   \item Scales (standardizes) the selected environmental variables
#'   \item Calculates Euclidean distances between neighbors using \code{\link[spdep]{nbcosts}}
#'   \item Constructs a new undirected graph with dissimilarity-weighted edges
#' }
#'
#' @seealso 
#' \code{\link[spdep]{nbcosts}} for neighbor cost calculation,
#' \code{\link[igraph]{graph_from_data_frame}} for graph construction
#'
#' @examples
#' \dontrun{
#' # Assuming g is a list with graph and polygons components
#' g_weighted <- add_edge_weight(g)
#' g_weighted_subset <- add_edge_weight(g, variables = c("temp", "precip"))
#' }
#'
#' @importFrom igraph as_adjacency_matrix graph_from_data_frame
#' @importFrom spdep mat2listw nbcosts
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom tidyselect all_of any_of  
#'
#' @export 
add_edge_weight <- function(g, id_col = NULL, variables = NULL) {
        
        adj_matrix <- igraph::as_adjacency_matrix(
                graph = g$graph, 
                type = "both", 
                names = TRUE, 
                sparse = TRUE
        )
        
        listw_object <- spdep::mat2listw(
                adj_matrix, 
                style = "B" 
        )
        
        nb_object <- listw_object$neighbours
        
        ## Prepare environmental variables 
        data_scaled <- g$polygons
        
        if (!is.null(variables)) {
                data_scaled <- dplyr::select(data_scaled, tidyselect::all_of(variables))
        } else {
                data_scaled <- dplyr::select(data_scaled, !tidyselect::any_of(id_col))
        }
        
        data_scaled <- sf::st_drop_geometry(data_scaled)
        data_scaled <- scale(data_scaled)
        
        netDist <- spdep::nbcosts(
                nb = nb_object, 
                data = data_scaled,
                method = "euclidean"
        )
        
        # Create edge data frame with weights
        edges_df <- data.frame(
                from = rep(1:length(nb_object), times = sapply(nb_object, length)),
                to = unlist(nb_object),
                weight = unlist(netDist)
        )
        
        # Create igraph object
        out <- igraph::graph_from_data_frame(
                d = edges_df,
                directed = FALSE,  
                vertices = NULL    
        )
        
        return(out)
}