test_that("compute_null_baseline re-raises errors that are not config errors", {
        fx  <- make_grid()                      # helper-grid.R
        obs <- get_regions(fx$weights, fx$graph, k = 3, seed = 1)
        local_mocked_bindings(
                .compute_ensemble_internal = function(...) stop("boom")
        )
        expect_error(
                compute_null_baseline(obs, fx$weights, fx$graph, B = 3, seed = 1),
                "boom"
        )
})