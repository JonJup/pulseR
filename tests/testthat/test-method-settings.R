## tests/testthat/test-method-settings.R
## cut_rule, refine_boundary and min_size are arguments, recorded in the object.

.ms_graph <- function() {
        x   <- make_grid_sf(8L, 8L)
        add_edge_weight(build_network(x, verbose = FALSE), id_col = "ID", verbose = FALSE)
}

test_that("global options no longer change the cut rule", {
        g  <- .ms_graph()
        tr <- sample_spanning_trees(g, n = 6, seed = 1, verbose = FALSE)
        run <- function() { set.seed(3); compute_ensemble(tr, 6, g, 4, verbose = FALSE)$memberships }
        a <- run()
        op <- options(pulseR.cut_rule = "ssd", pulseR.min_size = 5L); on.exit(options(op))
        expect_identical(run(), a)
})

test_that("compute_ensemble() records the cut settings, also in tuner slices", {
        g  <- .ms_graph()
        tr <- sample_spanning_trees(g, n = 6, seed = 1, verbose = FALSE)
        e  <- compute_ensemble(tr, 6, g, 4, cut_rule = "ssd", min_size = 2L, verbose = FALSE)
        expect_identical(e$cut_rule, "ssd")
        expect_identical(e$min_size, 2L)
        expect_true(all(apply(e$memberships, 2, function(z) min(tabulate(z))) >= 2L))
        expect_error(compute_ensemble(tr, 6, g, 4, cut_rule = "nope", verbose = FALSE))
})

test_that("get_regions() stores the settings in config, tuned and untuned", {
        g <- .ms_graph()
        one <- suppressWarnings(get_regions(g, final_regions = 2, N_tree = 10,
                                            regions_per_partition = 4, 
                                            refine_boundary = TRUE, min_size = 2L,
                                            seed = 1, verbose = FALSE))
        expect_identical(one$config$cut_rule, "ssd")
        expect_true(one$config$refine_boundary)
        expect_identical(one$config$min_size, 2L)
        tun <- suppressWarnings(get_regions(g, final_regions = 2:3, N_tree = 10, cut_rule = "propensity",
                                            regions_per_partition = 4, stability_B = 0,
                                            refine_boundary = TRUE, seed = 1, verbose = FALSE))
        expect_identical(tun$config$cut_rule, "propensity")
        expect_true(tun$config$refine_boundary)
})

test_that("tuning scores the refined partition it returns", {
        g  <- .ms_graph()
        tr <- sample_spanning_trees(g, n = 10, seed = 1, verbose = FALSE)
        msgs <- character(0)
        withCallingHandlers(
                tu <- tune_regions(g, tr, N_tree = 10, regions_per_partition = 4,
                                   final_regions = 2:3, stability_B = 0,
                                   refine_boundary = TRUE, seed = 1, verbose = FALSE),
                warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
        expect_false(any(grepl("differ from the scored winner", msgs)))
        expect_true(tu$refine_boundary)
})

test_that("refine_boundary needs a graph", {
        memb <- cbind(c(1, 1, 2, 2), c(1, 2, 2, 2), c(1, 1, 1, 2))
        expect_error(compute_consensus(memb, 2, refine_boundary = TRUE, verbose = FALSE), "graph")
})

