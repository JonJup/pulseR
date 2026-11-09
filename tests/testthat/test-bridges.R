# 4 x 4 grid of 1 m squares (x 0..4, y 0..4) plus one detached square at
# x 6..7, y 1..2: boundary gap 2 m, centroid gap 3 m to (3.5, 1.5).
make_island_sf <- function() {
        grid   <- make_grid_sf(4L, 4L)
        island <- sf::st_sfc(sf::st_polygon(list(rbind(c(6, 1), c(7, 1), c(7, 2),
                                                        c(6, 2), c(6, 1)))),
                             crs = sf::st_crs(grid))
        rbind(grid, sf::st_sf(ID = 17L, temp = 20, prec = 0.2, geometry = island))
}

test_that("in-memory network flags one bridge with boundary length", {
        net <- build_network(make_island_sf(), verbose = FALSE)
        g   <- net$graph
        expect_true(igraph::is_connected(g))
        expect_identical(sum(igraph::E(g)$bridge), 1L)
        b <- which(igraph::E(g)$bridge)
        expect_equal(igraph::E(g)$bridge_dist[b], 2, tolerance = 1e-8)
        expect_true(all(is.na(igraph::E(g)$bridge_dist[-b])))
        expect_identical(igraph::graph_attr(g, "bridge_dist_method"), "boundary")
        expect_equal(igraph::graph_attr(g, "bridge_scale_default"), 1, tolerance = 1e-8)
        # edge_data matches the graph, bridge row last
        ed <- net$edge_data
        expect_equal(nrow(ed), igraph::ecount(g))
        expect_true(ed$bridge[nrow(ed)])
        expect_identical(sort(c(ed$from[nrow(ed)], ed$to[nrow(ed)])),
                         sort(as.integer(igraph::ends(g, b, names = FALSE))))
})

test_that("centroid-only bridging records centroid length", {
        x  <- make_island_sf()
        nb <- suppressWarnings(spdep::poly2nb(x, queen = FALSE))
        el <- nb_to_edgelist(nb)
        g0 <- igraph::make_graph(as.vector(t(as.matrix(el))), n = nrow(x), directed = FALSE)
        g  <- connect_components(g0, centroids = sf::st_centroid(sf::st_geometry(x)),
                                 verbose = FALSE)
        b  <- which(igraph::E(g)$bridge)
        expect_length(b, 1L)
        expect_equal(igraph::E(g)$bridge_dist[b], 3, tolerance = 1e-8)
        expect_identical(igraph::graph_attr(g, "bridge_dist_method"), "centroid")
})

test_that("a connected graph gets bridge attributes but no bridges", {
        net <- build_network(make_grid_sf(3L, 3L), verbose = FALSE)
        expect_false(any(igraph::E(net$graph)$bridge))
        expect_true("bridge" %in% names(net$edge_data))
})

test_that("add_edge_weight() carries the flag and applies the penalty", {
        net <- build_network(make_island_sf(), verbose = FALSE)
        g   <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        br  <- igraph::E(g)$bridge
        expect_identical(sum(br), 1L)
        w   <- igraph::E(g)$weight
        we  <- igraph::E(g)$weight_env
        s   <- max(we[!br])
        expect_equal(w[!br], we[!br])
        expect_equal(w[br], we[br] + s * (1 - exp(-2 / 1)))
        expect_equal(igraph::graph_attr(g, "bridge_scale"), 1, tolerance = 1e-8)
        expect_identical(igraph::graph_attr(g, "bridge_dist_method"), "boundary")
})

test_that("bridge_scale = Inf leaves weights unpenalised; bad values stop", {
        net <- build_network(make_island_sf(), verbose = FALSE)
        g   <- add_edge_weight(net, id_col = "ID", bridge_scale = Inf, verbose = FALSE)
        expect_equal(igraph::E(g)$weight, igraph::E(g)$weight_env)
        expect_error(add_edge_weight(net, id_col = "ID", bridge_scale = 0, verbose = FALSE),
                     "bridge_scale")
        expect_error(add_edge_weight(net, id_col = "ID", bridge_scale = c(1, 2), verbose = FALSE),
                     "bridge_scale")
})

test_that("a user-supplied bridge_scale overrides the default", {
        net <- build_network(make_island_sf(), verbose = FALSE)
        g1  <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        g4  <- add_edge_weight(net, id_col = "ID", bridge_scale = 4, verbose = FALSE)
        br  <- igraph::E(g4)$bridge
        expect_lt(igraph::E(g4)$weight[br], igraph::E(g1)$weight[br])
        expect_equal(igraph::graph_attr(g4, "bridge_scale"), 4)
})

test_that("the null weight function reproduces the penalised weights", {
        net <- build_network(make_island_sf(), verbose = FALSE)
        g   <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        X   <- igraph::graph_attr(g, "env")
        expect_equal(.edge_weights_from_env(g, X), igraph::E(g)$weight)
})

test_that("graphs without a bridge flag keep their old behaviour", {
        net <- build_network(make_grid_sf(3L, 3L), verbose = FALSE)
        net$graph <- igraph::delete_edge_attr(net$graph, "bridge")
        expect_warning(add_edge_weight(net, id_col = "ID", bridge_scale = 5, verbose = FALSE),
                       "no `bridge` edge attribute")
        g <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        expect_equal(igraph::E(g)$weight, igraph::E(g)$weight_env)
})

test_that("Level 3 excludes bridges from reachability and edges", {
        net <- build_network(make_island_sf(), verbose = FALSE)
        g   <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        n   <- igraph::vcount(g)
        P   <- cbind(a = rep(c(0.8, 0.2), length.out = n),
                     b = rep(c(0.2, 0.8), length.out = n))
        g_nob <- igraph::delete_edges(g, which(igraph::E(g)$bridge))
        for (k in 1:2) {
                s_ex  <- compute_all_signatures(P, g, k, verbose = FALSE)
                s_ref <- compute_all_signatures(P, g_nob, k, exclude_bridges = FALSE,
                                                verbose = FALSE)
                expect_equal(unclass(s_ex)[, ], unclass(s_ref)[, ])
                expect_true(all(is.na(s_ex[n, ])))          # island: no edge left
                expect_true(isTRUE(attr(s_ex, "exclude_bridges")))
        }
        s_in <- compute_all_signatures(P, g, 1, exclude_bridges = FALSE, verbose = FALSE)
        expect_false(anyNA(s_in[n, ]))
})

test_that("the permutation null accepts a penalised graph", {
        net <- build_network(make_island_sf(), verbose = FALSE)
        g   <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        nb  <- compute_null_baseline(g, final_regions = 2, regions_per_partition = 4,
                                     N_tree = 4, B = 1, stability_B = 0, seed = 1,
                                     verbose = FALSE)
        expect_s3_class(nb, "pulse_null")
})


