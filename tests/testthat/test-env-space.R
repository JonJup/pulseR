test_that("add_edge_weight() stores the matrix behind the weights", {
        x <- make_grid_sf()
        net <- suppressMessages(build_network(x, verbose = FALSE))
        g <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        E <- igraph::graph_attr(g, "env")
        expect_true(is.matrix(E))
        expect_identical(colnames(E), c("temp", "prec"))
        expect_equal(unname(colMeans(E)), c(0, 0), tolerance = 1e-12)
        expect_equal(unname(apply(E, 2, stats::sd)), c(1, 1))
        ep <- igraph::ends(g, igraph::E(g), names = FALSE)
        expect_equal(igraph::E(g)$weight,
                     unname(sqrt(rowSums((E[ep[, 1], ] - E[ep[, 2], ])^2))))
})

test_that("an ID/id column is excluded by default; other ID-like columns warn", {
        x <- make_grid_sf()
        net <- suppressMessages(build_network(x, verbose = FALSE))
        g0 <- expect_silent(add_edge_weight(net, verbose = FALSE))
        expect_identical(colnames(igraph::graph_attr(g0, "env")), c("temp", "prec"))
        expect_equal(igraph::E(g0)$weight,
                     igraph::E(add_edge_weight(net, id_col = "ID", verbose = FALSE))$weight)
        net2 <- net; names(net2$polygons)[names(net2$polygons) == "ID"] <- "unit_no"
        expect_warning(add_edge_weight(net2, verbose = FALSE), "look like unit IDs")
        g1 <- suppressWarnings(add_edge_weight(net, variables = c("ID", "temp"), verbose = FALSE))
        expect_identical(colnames(igraph::graph_attr(g1, "env")), c("ID", "temp"))
        X <- cbind(id = 1:20, v = stats::rnorm(20))
        expect_warning(.warn_id_like(X, "`env_id`"), "id")
        expect_silent(.warn_id_like(X[, "v", drop = FALSE], "`env_id`"))
})

test_that(".resolve_env() defaults to the stored matrix", {
        x <- make_grid_sf()
        net <- suppressMessages(build_network(x, verbose = FALSE))
        g <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        expect_identical(.resolve_env(NULL, NULL, g), igraph::graph_attr(g, "env"))
        env <- sf::st_drop_geometry(x)
        expect_warning(.resolve_env(env, NULL, g), "look like unit IDs")
        expect_silent(r <- .resolve_env(env, "ID", g))
        expect_identical(colnames(r), c("temp", "prec"))
        g0 <- igraph::delete_graph_attr(g, "env")
        expect_null(.resolve_env(NULL, NULL, g0))
})

test_that("an ID column does not reach CH_FIT under the default env", {
        x <- make_grid_sf()
        net <- suppressMessages(build_network(x, verbose = FALSE))
        g <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        tr <- sample_spanning_trees(g, n = 8, seed = 1, verbose = FALSE)
        run <- function(env) suppressWarnings(tune_regions(
                g, tr, N_tree = 8, regions_per_partition = 6, final_regions = 2:3,
                stability_B = 0, tuning_metric = "ch_fit", env = env, seed = 1,
                verbose = FALSE))$tuning_log$CH_FIT
        expect_equal(run(NULL), run(igraph::graph_attr(g, "env")))
        expect_false(isTRUE(all.equal(run(NULL), run(sf::st_drop_geometry(x)))))
})

test_that("compute_null_baseline() refuses an env that does not reproduce the weights", {
        x <- make_grid_sf()
        net <- suppressMessages(build_network(x, verbose = FALSE))
        g <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        raw <- sf::st_drop_geometry(x)[, c("temp", "prec")]
        expect_error(compute_null_baseline(g, raw, final_regions = 2, regions_per_partition = 6,
                                           N_tree = 8, B = 1, stability_B = 0, verbose = FALSE),
                     "do not reproduce")
        g0 <- igraph::delete_graph_attr(g, "env")
        expect_error(compute_null_baseline(g0, final_regions = 2, regions_per_partition = 6,
                                           N_tree = 8, B = 1, stability_B = 0, verbose = FALSE),
                     "`env` is required")
        nb <- compute_null_baseline(g, final_regions = 2, regions_per_partition = 6,
                                    N_tree = 8, B = 1, stability_B = 0, seed = 1,
                                    verbose = FALSE)
        expect_s3_class(nb, "pulse_null")
        expect_true(all(is.finite(nb$null$CH_FIT)))
})

test_that("deprecated CH_ENV / ch_env names map to CH_FIT / ch_fit with a warning", {
        expect_warning(v <- .renamed_cvi(c("CH_ENV", "STAB"), "voters"), "renamed to CH_FIT")
        expect_identical(v, c("CH_FIT", "STAB"))
        expect_warning(m <- .renamed_cvi("r2_env", "tuning_metric"), "r2_fit")
        expect_identical(m, "r2_fit")
        expect_silent(.renamed_cvi(c("CH_FIT", "STAB"), "voters"))
})
