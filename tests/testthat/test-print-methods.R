## Print methods: a short summary, never the full matrix.

.print_lines <- function(x) {
        out <- utils::capture.output(res <- withVisible(print(x)))
        expect_false(res$visible)
        expect_identical(res$value, x)
        out
}

test_that("print.pulse_regions_tuning summarises the search", {
        skip_if_not_installed("salso")
        x  <- make_grid_sf()
        g  <- add_edge_weight(suppressMessages(build_network(x, verbose = FALSE)),
                              id_col = "ID", verbose = FALSE)
        tr <- sample_spanning_trees(g, n = 8, seed = 1, verbose = FALSE)
        tu <- suppressWarnings(tune_regions(g, tr, N_tree = 8, regions_per_partition = 6,
                                            final_regions = 2:3, stability_B = 0,
                                            seed = 1, verbose = FALSE))
        out <- .print_lines(tu)
        expect_match(out[1], "^<pulse_regions_tuning>")
        expect_true(any(grepl("best: final_regions = ", out)))
        expect_lte(length(out), 10L)
        ex <- structure(list(prior_strength = Inf, best_final_regions = 3L,
                             best_result = list(memberships = diag(3))),
                        class = "pulse_regions_tuning")
        expect_match(.print_lines(ex), "exact prior")
})

test_that("print.pulse_signatures and print.pulse_configuration do not dump the matrix", {
        set.seed(1)
        g <- igraph::make_lattice(c(12, 12))
        P <- matrix(stats::runif(144 * 3), 144, 3)
        P <- P / rowSums(P)
        s <- compute_all_signatures(P, g, k_order = 2, verbose = FALSE)
        out <- .print_lines(s)
        expect_match(out[1], "^<pulse_signatures>  144 units x 6 class pairs")
        expect_lte(length(out), 5L)
        expect_false(inherits(head(s), "pulse_signatures"))

        cf <- compute_configuration(P, g, k_order = 3, lags = 1:2, verbose = FALSE)
        out <- .print_lines(cf)
        expect_match(out[1], "^<pulse_configuration>  144 units x 2 features")
        expect_true(any(grepl("kappa_lag1 = ", out)))
        expect_lte(length(out), 6L)
        expect_false(inherits(head(cf), "pulse_configuration"))
})
