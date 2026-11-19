## B8 (round 5): get_regions() drew the cut race from different points of the
## RNG stream in its tuned and untuned branches, so refitting a tuned pick as a
## single configuration with the same seed gave a different ensemble.
## Needs the matching change in R/regions.R (untuned branch draws `ens_seed`).

.b8_graph <- function(n = 12) {
        set.seed(1)
        g  <- igraph::make_lattice(c(n, n))
        xy <- expand.grid(x = seq_len(n), y = seq_len(n))
        env <- cbind(e1 = 3 * cut(xy$x, 3, labels = FALSE) + stats::rnorm(n^2, sd = 0.3),
                     e2 = stats::rnorm(n^2, sd = 0.3))
        ep <- igraph::ends(g, igraph::E(g), names = FALSE)
        igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], ] - env[ep[, 2], ])^2)) + 1e-6
        igraph::set_graph_attr(g, "env", env)
}

test_that("B8: refitting the tuned pick with the same seed reproduces it", {
        skip_if_not_installed("salso")
        g <- .b8_graph()
        for (rule in c("propensity", "ssd")) {
                tu <- suppressWarnings(get_regions(
                        g, regions_per_partition = 9, final_regions = 2:4, N_tree = 20,
                        stability_B = 0, cut_rule = rule, seed = 11, verbose = FALSE))
                re <- suppressWarnings(get_regions(
                        g, regions_per_partition = 9, final_regions = tu$k, N_tree = 20,
                        cut_rule = rule, seed = 11, verbose = FALSE))
                expect_true(tu$tuned)
                expect_false(re$tuned)
                expect_identical(re$ensemble, tu$ensemble, info = rule)
                expect_equal(.adjusted_rand_index(re$hard_labels, tu$hard_labels), 1,
                             info = rule)
        }
})
