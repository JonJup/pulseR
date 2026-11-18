test_that("default variable selection never uses true_* columns", {
        n <- nrow(pulse_grid)
        W <- cbind(as.numeric(seq_len(n) <= n / 2), as.numeric(seq_len(n) > n / 2))
        expect_message(
                loc <- suppressWarnings(get_local_types(pulse_grid, W, n_local_types = 1:2)),
                "true_region, true_local")
        expect_identical(loc$variables, c("temp", "prec"))
        
        net <- build_network(pulse_grid, verbose = FALSE)
        expect_message(gw <- add_edge_weight(net, verbose = FALSE), "true_region, true_local")
        expect_identical(colnames(igraph::graph_attr(gw, "env")), c("temp", "prec"))
        
        env <- sf::st_drop_geometry(pulse_grid)
        expect_message(r <- .resolve_env(env, "ID", gw), "true_region")
        expect_identical(colnames(r), c("temp", "prec"))
})

test_that("true_* columns named in `variables` are still used", {
        net <- build_network(pulse_grid, verbose = FALSE)
        gw  <- add_edge_weight(net, variables = c("temp", "true_region"), verbose = FALSE)
        expect_identical(colnames(igraph::graph_attr(gw, "env")), c("temp", "true_region"))
})