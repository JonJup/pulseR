## tests/testthat/test-ensemble-diagnostics.R

.toy_trees <- function(n_tree = 30, side = 20, prior = NULL) {
        set.seed(1)
        g  <- igraph::make_lattice(c(side, side))
        xy <- expand.grid(x = seq_len(side), y = seq_len(side))
        env <- cbind(cut(xy$x, 3, labels = FALSE) * 2 + stats::rnorm(side^2, sd = 0.5),
                     stats::rnorm(side^2))
        ep <- igraph::ends(g, igraph::E(g), names = FALSE)
        igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], ] - env[ep[, 2], ])^2))
        igraph::E(g)$eid    <- seq_len(igraph::ecount(g))
        list(g = g, tr = sample_spanning_trees(g, n = n_tree, seed = 1, verbose = FALSE,
                                               prior_typology = prior))
}

test_that("prior shift is lambda standard deviations of gamma * log(w) on any scale", {
        set.seed(3)
        w <- stats::rlnorm(200); s <- stats::runif(200, 1, 50)
        wp <- rep(c(TRUE, FALSE), 100)
        for (c in c(0.3, 1, 5)) for (gam in c(0.5, 2)) {
                a <- .cut_log_propensity(w^c, s, gam, prior_strength = 1.5, within_prior = wp)
                b <- .cut_log_propensity(w^c, s, gam, prior_strength = 0,   within_prior = wp)
                d <- a$lp - b$lp
                expect_equal((mean(d[wp]) - mean(d[!wp])) / stats::sd(gam * log(w^c)), -1.5)
        }
})

test_that("zero weight spread falls back to halvings, with one warning", {
        r <- .cut_log_propensity(rep(1, 50), stats::runif(50, 1, 10),
                                 prior_strength = 2, within_prior = rep(c(TRUE, FALSE), 25))
        expect_true(r$fallback)
        expect_equal(r$unit, log(2))
        x  <- .toy_trees(5, 15, prior = rep(1:3, each = 75))
        for (j in 1:5) igraph::E(x$tr$trees[[j]])$weight <- 1
        expect_warning(.compute_ensemble_internal(x$g, x$tr, 5, 4, prior_strength = 1),
                       "halvings")
})

test_that("lambda = 0 leaves cut selection unchanged", {
        x <- .toy_trees(10)
        set.seed(5); a <- .compute_ensemble_internal(x$g, x$tr, 10, 8)
        set.seed(5); b <- .compute_ensemble_internal(x$g, x$tr, 10, 8,
                                                                  prior_strength = 0)
        expect_identical(a, b)
})

test_that("check_cut_placement anchors random at ~0 and size_influence = 0 at 1", {
        x  <- .toy_trees(30)
        cp <- check_cut_placement(x$tr, regions_per_partition = 10,
                                  size_influence = c(1, 2), seed = 2)
        s  <- cp$summary
        expect_s3_class(cp, "pulse_cut_placement")
        expect_equal(s$setting, c("random", "0", "1", "2"))
        expect_equal(s$env_retained[s$setting == "random"], 0)
        expect_equal(s$env_retained[s$setting == "0"], 1)
        expect_gt(s$above_median[s$setting == "0"], s$above_median[s$setting == "random"])
        # More size weight -> fewer one-unit cuts.
        expect_lt(s$singleton_share[s$setting == "2"], s$singleton_share[s$setting == "0"])
})

test_that("check_ensemble_size returns a curve and a search floor", {
        skip_if_not_installed("salso")
        x   <- .toy_trees(80)
        ens <- compute_ensemble(x$tr, 80, x$g, 10, verbose = FALSE)
        chk <- check_ensemble_size(ens, k = 3, B = 2, seed = 1, verbose = FALSE)
        expect_s3_class(chk, "pulse_ensemble_size")
        expect_equal(chk$summary$N_tree, c(10L, 20L, 40L))
        expect_true(all(chk$replicates$ARI <= 1))
        expect_true(is.logical(chk$converged))
        expect_error(check_ensemble_size(ens, k = 3, sizes = 41), "half the ensemble")
})
