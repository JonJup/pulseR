# Contiguity check for tests: does every label form one connected component?
# Moved from R/region_helper.R (review D1); not used by package code.
.check_contiguity <- function(graph, lab) {
        n  <- igraph::vcount(graph)
        stopifnot(length(lab) == n)
        ep <- igraph::ends(graph, igraph::E(graph), names = FALSE)
        
        same <- !is.na(lab[ep[, 1]]) & !is.na(lab[ep[, 2]]) & lab[ep[, 1]] == lab[ep[, 2]]
        g2   <- igraph::make_graph(as.vector(t(ep[same, , drop = FALSE])),
                                   n = n, directed = FALSE)
        comp <- igraph::components(g2)$membership
        
        keep  <- !is.na(lab)
        frag  <- tapply(comp[keep], lab[keep], function(x) length(unique(x)))
        list(contiguous = all(frag == 1L),
             fragments  = frag[frag > 1L],
             n_orphans  = sum(keep) - sum(tapply(comp[keep], lab[keep],
                                                 function(x) max(tabulate(match(x, unique(x)))))))
}
