test_that("chunked mode aligns rows by ID, not position", {
        skip_if_not_installed("arrow")
        g <- list(graph = igraph::set_vertex_attr(igraph::make_ring(6),
                                                  "polygon_id", value = 101:106))
        env <- data.frame(ID = 101:106, x = c(1, 5, 2, 8, 3, 9), y = c(6, 2, 5, 1, 4, 3))
        f1 <- tempfile(fileext = ".parquet"); f2 <- tempfile(fileext = ".parquet")
        arrow::write_parquet(env, f1)
        arrow::write_parquet(env[c(4, 1, 6, 2, 5, 3), ], f2)   # shuffled
        w1 <- igraph::E(add_edge_weight(g, chunked = TRUE, files = f1))$weight
        w2 <- igraph::E(add_edge_weight(g, chunked = TRUE, files = f2))$weight
        expect_equal(w2, w1)
        f3 <- tempfile(fileext = ".parquet")
        arrow::write_parquet(env[-2, ], f3)                     # missing vertex 102
        expect_error(add_edge_weight(g, chunked = TRUE, files = f3), "no environmental row")
})