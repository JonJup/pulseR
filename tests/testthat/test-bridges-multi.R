## connect_components() with more than one island (round-5 review, section 4.2).
## 5 x 5 grid of 1 m squares (0..5) plus two detached squares, each 2 m from
## the grid and about 7.8 m from each other: both must bridge to the grid.

test_that("two islands get one bridge each, to the mainland", {
        grid <- make_grid_sf(5L, 5L)
        sq   <- function(x0, y0) sf::st_polygon(list(rbind(
                c(x0, y0), c(x0 + 1, y0), c(x0 + 1, y0 + 1), c(x0, y0 + 1), c(x0, y0))))
        isl  <- sf::st_sf(ID = 26:27, temp = c(20, 21), prec = c(0.2, 0.3),
                          geometry = sf::st_sfc(sq(7, 1), sq(0, 7), crs = sf::st_crs(grid)))
        net  <- build_network(rbind(grid, isl), verbose = FALSE)
        g    <- net$graph
        expect_true(igraph::is_connected(g))
        b <- which(igraph::E(g)$bridge)
        expect_length(b, 2L)
        expect_equal(igraph::E(g)$bridge_dist[b], c(2, 2), tolerance = 1e-8)
        eb <- igraph::ends(g, b, names = FALSE)
        expect_true(all(pmin(eb[, 1], eb[, 2]) <= 25L))  # one end on the grid
        expect_setequal(pmax(eb[, 1], eb[, 2]), 26:27)   # the other on each island
})
