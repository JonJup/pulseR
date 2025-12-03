#' Connect Disconnected Graph Components Using Spatial Proximity
#'
#' @description
#' Iteratively connects disconnected components of a graph by adding edges
#' between spatially closest nodes from different components. The function uses
#' polygon centroids to compute spatial distances and applies a minimum spanning
#' tree approach to ensure efficient connectivity.
#'
#' @param g An \code{igraph} object representing the graph with potentially
#'   disconnected components.
#' @param polygons An \code{sf} or \code{sfc} object containing spatial polygons
#'   corresponding to the nodes in the graph. The number of polygons must match
#'   the number of vertices in \code{g}.
#'
#' @return An \code{igraph} object with all components connected. The returned
#'   graph will have a single connected component.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Computing centroids of the input polygons
#'   \item Calculating the pairwise Euclidean distance matrix between all centroids
#'   \item Identifying disconnected components in the graph
#'   \item For each iteration while multiple components exist:
#'     \itemize{
#'       \item Computing minimum distances between all pairs of components
#'       \item Creating a component-level graph weighted by these distances
#'       \item Finding the minimum spanning tree (MST) of the component graph
#'       \item Adding an edge between the two nodes that represent the shortest
#'             connection between components
#'     }
#' }
#'
#' The algorithm continues until all components are connected into a single
#' connected component.
#'
#' @importFrom igraph components add_edges graph_from_adjacency_matrix mst
#'   as_edgelist edge.attributes
#' @importFrom sf st_centroid st_coordinates
#' @importFrom stats dist
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' library(sf)
#'
#' # Create a simple graph with two components
#' g <- make_empty_graph(n = 6) |>
#'   add_edges(c(1, 2, 2, 3)) |>
#'   add_edges(c(4, 5, 5, 6))
#'
#' # Create corresponding spatial polygons
#' polygons <- st_as_sf(data.frame(
#'   id = 1:6,
#'   geometry = st_sfc(
#'     st_point(c(0, 0)), st_point(c(1, 0)), st_point(c(2, 0)),
#'     st_point(c(10, 0)), st_point(c(11, 0)), st_point(c(12, 0))
#'   ) |> st_buffer(0.5)
#' ))
#'
#' # Connect the components
#' g_connected <- connect_components(g, polygons)
#' components(g_connected)$no  # Should be 1
#' }
#'
#' @seealso \code{\link[igraph]{components}}, \code{\link[igraph]{mst}}
#'
#' @export
connect_components <- function(g, polygons) {

        centroids <- suppressWarnings(st_centroid(polygons))
        coords <- st_coordinates(centroids)
        distMat <- dist(coords)
        distMat <- as.matrix(distMat)
        eval_id <- 1
        comps <- components(g)

        while(comps$no > 1) {

                #if (eval_id %% 10 == 0) cat(sprintf("Eval ID %d ...\n", eval_id))

                comps <- components(g)
                #cat(comps$no)
                if (comps$no == 1) return(g)
                comp_nodes <- split(seq_along(comps$membership), comps$membership)

                # function that gets min dist between component i and j
                min_dist <- function(a, b) min(distMat[a, b])

                # compute symmetric matrix (k × k)
                adj_matrix <- outer(
                        comp_nodes,
                        comp_nodes,
                        Vectorize(min_dist)
                )

                diag(adj_matrix) <- 0

                comp_graph <- graph_from_adjacency_matrix(
                        adj_matrix,
                        mode = "undirected",
                        weighted = TRUE
                )
                comp_mst <- mst(comp_graph)

                # Add edges from MST to original graph
                mst_edges <- as_edgelist(comp_mst)

                # Identify the combination with the shortest connection
                shory <- mst_edges[which.min(igraph::edge.attributes(comp_mst)$weight), ]

                # combine these two components
                shorty_id <- which(distMat == min(igraph::edge.attributes(comp_mst)$weight), arr.ind = TRUE)[1, ]
                g <- add_edges(g, shorty_id)

                comps <- components(g)
                eval_id <- eval_id + 1
        }

        return(g)
}
