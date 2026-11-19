## Simulation tests check statistical properties over many seeds and take
## minutes. They never run on CRAN, and locally only when asked for:
##   Sys.setenv(PULSER_SLOW_TESTS = "true"); devtools::test(filter = "simulation")
skip_if_not_slow <- function() {
        skip_on_cran()
        if (!identical(tolower(Sys.getenv("PULSER_SLOW_TESTS")), "true"))
                skip("simulation test; set PULSER_SLOW_TESTS=true to run")
}

## n x n rook lattice whose edge weights are Euclidean distances on `env`,
## with `env` stored on the graph (as add_edge_weight(scale = FALSE) would).
.sim_graph <- function(n, env) {
        g  <- igraph::make_lattice(c(n, n))
        ep <- igraph::ends(g, igraph::E(g), names = FALSE)
        igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], , drop = FALSE] -
                                             env[ep[, 2], , drop = FALSE])^2))
        igraph::set_graph_attr(g, "env", env)
}
