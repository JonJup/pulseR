#' Fuzzy Contiguity-Constrained Clustering Using SKATER Algorithm
#'
#' Performs spatially constrained fuzzy clustering using a modified SKATER
#' (Spatial 'K'luster Analysis by Tree Edge Removal) algorithm. Multiple random
#' spanning trees generate an ensemble of regionalizations, which are then
#' combined via coassociation-based fuzzy clustering.
#'
#' @param graph An \code{igraph} object representing the spatial contiguity
#'   graph. Edges must have a \code{weight} attribute representing dissimilarity
#'   or cost between connected nodes.
#' @param skater_regions Integer. The number of regions to generate in each
#'   SKATER iteration. Default is 50. This determines the granularity of the
#'   initial partitioning before consensus clustering.
#' @param final_regions Integer. The number of final fuzzy clusters to produce
#'   via PAM clustering on the coassociation matrix.
#' @param n.rst Integer. The number of random spanning trees to generate. Each
#'   spanning tree produces an independent regionalization, which contributes
#'   to the coassociation matrix for ensemble consensus.
#' @param fuzzyness_parameter Numeric. The fuzziness exponent (m) for computing
#'   fuzzy memberships. Default is 2. Higher values produce fuzzier (more
#'   distributed) memberships; values approaching 1 produce crisper assignments.
#'
#' @return A numeric matrix of fuzzy membership values with dimensions
#'   \code{n_nodes x final_regions}. Each row sums to 1 and represents the
#'   degree of membership of each node to each of the \code{final_regions}
#'   clusters. Row order corresponds to node order in the input graph.
#'
#' @details
#' The function implements a fuzzy ensemble variant of the SKATER algorithm:
#' \enumerate{
#'   \item Generate \code{n.rst} random spanning trees from the input graph.
#'   \item For each spanning tree, iteratively remove the edge with the highest
#'     weight (dissimilarity) until \code{skater_regions} connected components
#'     remain.
#'   \item Compute a coassociation matrix using Hamming distance across all
#'     regionalization results, measuring how often pairs of nodes are assigned
#'     to the same region.
#'   \item Apply PAM (Partitioning Around Medoids) clustering to the
#'     coassociation matrix to identify \code{final_regions} representative
#'     partitions.
#'   \item Convert distances to medoids into fuzzy membership values using the
#'     specified fuzziness parameter.
#' }
#'
#' The fuzzy membership for node \eqn{i} to cluster \eqn{j} is calculated as:
#' \deqn{u_{ij} = \frac{d_{ij}^{-2/(m-1)}}{\sum_{k=1}^{K} d_{ik}^{-2/(m-1)}}}
#' where \eqn{d_{ij}} is the distance from node \eqn{i} to medoid \eqn{j}, and
#' \eqn{m} is the fuzziness parameter.
#'
#' @seealso
#' \code{\link[igraph]{sample_spanning_tree}} for spanning tree generation,
#' \code{\link[igraph]{components}} for connected component detection,
#' \code{\link[cluster]{pam}} for PAM clustering,
#' \code{\link[parallelDist]{parallelDist}} for distance computation,
#' \code{\link[spdep]{skater}} for the original SKATER implementation.
#'
#' @references
#' Assunção, R. M., Neves, M. C., Câmara, G., & Da Costa Freitas, C. (2006).
#' Efficient regionalization techniques for socio-economic geographical units
#' using minimum spanning trees. \emph{International Journal of Geographical
#' Information Science}, 20(7), 797-811.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#'
#' # Create example graph with weights
#' g <- make_lattice(c(10, 10))
#' E(g)$weight <- runif(ecount(g))
#'
#' # Generate fuzzy clustering with 5 final regions from 100 random trees
#' fuzzy_result <- skater_con(
#'   graph = g,
#'   skater_regions = 20,
#'   final_regions = 5,
#'   n.rst = 100,
#'   fuzzyness_parameter = 2
#' )
#'
#' # View fuzzy memberships for first 6 nodes
#' head(fuzzy_result)
#'
#' # Get hard assignments from fuzzy memberships
#' hard_clusters <- apply(fuzzy_result, 1, which.max)
#' }
#'
#' @importFrom igraph subgraph_from_edges sample_spanning_tree E delete_edges components
#' @importFrom purrr map2
#' @importFrom parallelDist parallelDist
#' @importFrom cluster pam
#'
#' @export
#'
skater_con <- function(graph, skater_regions = 50, final_regions, n.rst, fuzzyness_parameter = 2) {
        rst.l <- vector(mode = "list", length = n.rst)
        rst.l <-
                lapply(
                        seq_along(rst.l), 
                        function(x)
                                igraph::subgraph_from_edges(
                                        graph, 
                                        sample_spanning_tree(graph)
                                )
                )
        for (k in 1:(skater_regions - 1)) {
                high.cost.edge <- lapply(rst.l, function(x) which.max(E(x)$weight))
                high.cost.edge <- unlist(high.cost.edge)
                rst.l <- purrr::map2(rst.l, high.cost.edge, \(x, y) igraph::delete_edges(graph = x, edges = y))
        }
        mmbsp <- lapply(
                rst.l, 
                function(x)
                        igraph::components(x)$membership
        )
        mmbsp <- as.data.frame(mmbsp)
        names(mmbsp) <-
                paste("regionalization", 1:n.rst, sep = "_")
        
        coassociation_matrix <- parallelDist::parallelDist(x = as.matrix(mmbsp), method = "hamming")
        pam_result <- cluster::pam(coassociation_matrix, k = final_regions)
        medoid_indices <- pam_result$id.med
        
        # Step 2: Extract distances from each point to each medoid
        # This requires the full distance matrix
        dist_mat_full <- as.matrix(coassociation_matrix)
        distances_to_medoids <- dist_mat_full[, medoid_indices]
        
        # Step 3: Convert to fuzzy memberships
        # fuzziness parameter
        m <- fuzzyness_parameter  
        fuzzy_memberships <- t(apply(distances_to_medoids, 1, function(d) {
                if(any(d == 0)) {
                        # Handle case where point is a medoid
                        u <- rep(0, length(d))
                        u[which.min(d)] <- 1
                        return(u)
                }
                d_powered <- d^(-2/(m-1))
                d_powered / sum(d_powered)
        }))
        return(fuzzy_memberships)
}