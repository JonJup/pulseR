## Simulation tests (round-5 review, section 4.3). Tolerances are counts over
## seeds ("k of n"), so a test fails when the method degrades, not by chance.
## Opt-in: see helper-slow.R. Default pipeline settings throughout.

test_that("null p-values are not anti-conservative on a field without signal", {
        skip_if_not_slow()
        skip_if_not_installed("salso")
        # iid noise: the observed map is exchangeable with every permutation,
        # so p is uniform on {1/20, ..., 1} and P(p <= 0.05) = 0.05.
        p <- vapply(1:20, function(s) {
                set.seed(s)
                g   <- .sim_graph(10, cbind(stats::rnorm(100), stats::rnorm(100)))
                obs <- suppressWarnings(get_regions(g, final_regions = 2, N_tree = 20,
                                                    stability_B = 0, seed = s,
                                                    verbose = FALSE))
                nb  <- suppressWarnings(compute_null_baseline(
                        g, final_regions = 2, N_tree = 20, B = 19, stability_B = 0,
                        observed = obs, seed = 1000 + s, verbose = FALSE))
                nb$summary$p[nb$summary$criterion == "CH_FIT"]
        }, numeric(1))
        expect_true(all(is.finite(p)))
        expect_lte(sum(p <= 0.05), 3L)     # expected 1 of 20; P(>= 4) is about 0.016
        expect_gt(mean(p), 0.3)            # uniform on {1/20, ..., 1} has mean 0.525
})

test_that("Level 1 recovers three planted regions in at least 9 of 10 seeds", {
        skip_if_not_slow()
        skip_if_not_installed("salso")
        n     <- 15
        xy    <- expand.grid(x = seq_len(n), y = seq_len(n))
        truth <- cut(xy$x, 3, labels = FALSE)            # three 5-column bands
        ari <- vapply(1:10, function(s) {
                set.seed(s)
                env <- cbind(truth + stats::rnorm(n^2, sd = 0.3), stats::rnorm(n^2, sd = 0.3))
                r   <- suppressWarnings(get_regions(.sim_graph(n, env), final_regions = 3,
                                                    N_tree = 30, seed = s, verbose = FALSE))
                .adjusted_rand_index(r$hard_labels, truth)
        }, numeric(1))
        expect_gte(sum(ari > 0.8), 9L)
})

test_that("DMM: ICL recovers the planted number of types in at least 8 of 10 seeds", {
        skip_if_not_slow()
        # Same generator as .dmm_sim_fixture() in test-mosaic-dmm.R, new seeds.
        k_hat <- vapply(1:10, function(s) {
                set.seed(s); K <- 3; C <- 20; N <- 600; A <- c(20, 50, 200)
                m <- matrix(stats::rgamma(K * C, 0.7), K); m <- m / rowSums(m)
                z <- sample(K, N, TRUE)
                X <- t(vapply(seq_len(N), function(i) {
                        p <- stats::rgamma(C, A[z[i]] * m[z[i], ]); p <- p / sum(p)
                        as.numeric(stats::rmultinom(1, 150, p))
                }, numeric(C)))
                icl <- vapply(2:5, function(k) .fit_dmm(X, k, n_start = 3, seed = s)$ICL,
                              numeric(1))
                (2:5)[which.min(icl)]
        }, numeric(1))
        expect_gte(sum(k_hat == 3), 8L)
})
