
#' Exact reproduction of a prior typology (\code{prior_strength = Inf})
#'
#' Returns the prior labels directly as a (degenerate) clustering result, with a
#' length-N \code{coverage} flag marking vertices that have no prior label.
#' Uncovered vertices receive \code{NA} hard labels and all-zero membership rows.
#'
#' @param prior_labels Character vector of length N (resolved labels).
#' @param verbose Logical.
#' @return A list mirroring the \code{get_regions()} return, plus \code{coverage}
#'   and \code{class_levels}.
#' @keywords internal
#' @noRd
.exact_prior_result <- function(prior_labels, verbose = TRUE) {
        N         <- length(prior_labels)
        coverage  <- !is.na(prior_labels)
        n_covered <- sum(coverage)
        n_dropped <- N - n_covered
        
        if (n_covered == 0L) {
                stop("`prior_strength = Inf` but no vertices have a prior label.",
                     call. = FALSE)
        }
        
        if (verbose) {
                message("prior_strength = Inf: returning prior typology directly.")
                if (n_dropped > 0L) {
                        message(sprintf(
                                "  %d / %d vertices have no prior label and are marked NA (coverage = FALSE).",
                                n_dropped, N
                        ))
                }
        }
        
        classes <- sort(unique(prior_labels[coverage]))
        K       <- length(classes)
        
        hard_clusters           <- rep(NA_integer_, N)
        hard_clusters[coverage] <- match(prior_labels[coverage], classes)
        
        # One-hot membership for covered vertices; all-zero rows for uncovered ones.
        memberships <- matrix(0, nrow = N, ncol = K)
        memberships[cbind(which(coverage), hard_clusters[coverage])] <- 1
        
        list(
                memberships         = memberships,
                hard_clusters       = hard_clusters,
                coverage            = coverage,        # length-N logical
                best_regions_per_partition = K,
                best_final_regions  = K,
                best_N_tree          = NA_integer_,
                best_score          = NA_real_,
                tuning_log          = NULL,
                #N_tree_stability     = NULL,
                trees               = NULL,
                cluster_result      = NULL,
                prior_strength      = Inf,
                class_levels        = classes
        )
}



#' Resolve a \code{prior_typology} argument to a vertex-label vector
#'
#' Accepts either a single character (vertex-attribute name on \code{graph}) or
#' a vector of length \code{vcount(graph)}. Coerces to character; empty strings
#' become NA. Returns NULL if input is NULL.
#'
#' @param graph An \code{igraph} object.
#' @param prior_typology See \code{\link{sample_spanning_trees}}.
#' @return Character vector of length \code{vcount(graph)}, or \code{NULL}.
#' @keywords internal
#' @importFrom igraph vcount vertex_attr vertex_attr_names
#' @noRd
.resolve_prior_typology <- function(graph, prior_typology) {
        #if ((prior_typology)) return(NULL)
        
        n_v <- igraph::vcount(graph)
        
        if (is.character(prior_typology) && length(prior_typology) == 1L) {
                if (!prior_typology %in% igraph::vertex_attr_names(graph)) {
                        stop(sprintf("Vertex attribute '%s' not found on graph. ", prior_typology),
                             "Attach the typology via V(graph)$", prior_typology,
                             " <- ... before calling.", call. = FALSE)
                }
                labs <- igraph::vertex_attr(graph, prior_typology)
        } else if (length(prior_typology) == n_v) {
                labs <- prior_typology
        } else {
                stop("`prior_typology` must be either a single character giving the name of ",
                     "a vertex attribute on `graph`, or a vector of length vcount(graph) = ",
                     n_v, ".", call. = FALSE)
        }
        
        labs <- as.character(labs)
        labs[!nzchar(labs)] <- NA_character_
        
        # Sanity check: at least 2 distinct non-NA labels for weighting to do anything.
        uniq <- unique(labs[!is.na(labs)])
        if (length(uniq) < 2L) {
                warning("Prior typology contains fewer than 2 distinct non-NA labels; ",
                        "weighting will have no effect.", call. = FALSE)
        }
        labs
}

#' Flag edges whose endpoints share a non-NA prior label
#'
#' @param prior_labels Character vector of vertex labels, or \code{NULL}.
#' @param edge_endpoints Two-column integer matrix of edge endpoints.
#' @return Logical vector of length \code{nrow(edge_endpoints)}, or \code{NULL}.
#' @keywords internal
#' @noRd
.compute_within_prior_edge <- function(prior_labels, edge_endpoints) {
        if (is.null(prior_labels)) return(NULL)
        a <- prior_labels[edge_endpoints[, 1]]
        b <- prior_labels[edge_endpoints[, 2]]
        !is.na(a) & !is.na(b) & a == b
}

#' Apply the within-prior halving rule to edge weights
#'
#' \eqn{\tilde{w}(i,j) = w(i,j) / 2^{\lambda_E}} when both endpoints share a
#' non-NA prior label; \eqn{\tilde{w}(i,j) = w(i,j)} otherwise. Thus
#' \eqn{\lambda_E = 0} returns weights unchanged; \eqn{\lambda_E = 1} halves
#' within-prior weights; \eqn{\lambda_E = 2} quarters them, etc. Lowering these
#' weights makes the corresponding edges less likely to be cut, so prior regions
#' tend to stay intact.
#'
#' @param weights Numeric edge-weight vector.
#' @param within_prior_edge Logical vector (same length), or \code{NULL}.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda_E}.
#' @return Numeric vector of (possibly down-weighted) edge weights.
#' @keywords internal
#' @noRd
.apply_prior_weighting <- function(weights, within_prior_edge, prior_strength) {
        if (is.na(within_prior_edge) || prior_strength == 0) return(weights)
        weights <- lapply(weights, 
                          function(x) x / (2 ^ (prior_strength * as.numeric(within_prior_edge))
                          )
        )
}
