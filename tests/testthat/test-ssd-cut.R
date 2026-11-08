## tests/testthat/test-ssd-cut.R
## SSD cut rule: size tilt, prior penalty, min_size; refinement min_size guard.

.ssd_fixture <- function(n_tree = 20, side = 20, prior = TRUE) {
        set.seed(1)
        g   <- igraph::make_lattice(c(side, side))
        xy  <- expand.grid(x = seq_len(side), y = seq_len(side))
        env <- cbind(cut(xy$x, 3, labels = FALSE) * 2 + stats::rnorm(side^2, sd = 0.5),
                     stats::rnorm(side^2))
        ep  <- igraph::ends(g, igraph::E(g), names = FALSE)
        igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], ] - env[ep[, 2], ])^2))
        igraph::E(g)$eid    <- seq_len(igraph::ecount(g))
        g   <- igraph::set_graph_attr(g, "env", env)
        pr  <- if (prior) ifelse(xy$y <= side / 2, 1L, 2L) else NULL
        list(g = g, ep = ep, prior = pr,
             tr = sample_spanning_trees(g, n = n_tree, seed = 1, verbose = FALSE,
                                        prior_typology = pr))
}

.memb_ssd <- function(x, k, ...) {
        op <- options(pulseR.cut_rule = "ssd"); on.exit(options(op))
        .compute_ensemble_memberships_internal(x$g, x$tr, length(x$tr$trees), k, ...)
}

test_that(".forest_cut_gains returns Ward gain, h and smaller side", {
        ep <- cbind(1:3, 2:4)                                  # path 1-2-3-4
        fc <- .forest_cut_gains(ep, matrix(c(0, 0, 10, 10)))
        expect_equal(fc$h,     c(0.75, 1, 0.75))
        expect_equal(fc$gain,  c(0.75 * (20 / 3)^2, 100, 0.75 * (20 / 3)^2))
        expect_equal(fc$n_min, c(1, 2, 1))
})

test_that("defaults reproduce the plain Ward race draw for draw", {
        x <- .ssd_fixture()
        ref <- function(x, k, gam) {                           # pre-change branch
                X <- igraph::graph_attr(x$g, "env")
                vapply(x$tr$trees, function(tree) {
                        ep <- igraph::ends(tree, igraph::E(tree), names = FALSE)
                        for (s in seq_len(k - 1L)) {
                                g_e <- pmax(.forest_cut_gains(ep, X)$gain, .Machine$double.xmin)
                                ep  <- ep[-which.min(log(stats::rexp(nrow(ep))) - gam * log(g_e)), ,
                                          drop = FALSE]
                        }
                        as.integer(igraph::components(igraph::make_graph(
                                as.vector(t(ep)), n = nrow(X), directed = FALSE))$membership)
                }, integer(nrow(X)))
        }
        set.seed(7); a <- ref(x, 6, 4)
        set.seed(7); b <- .memb_ssd(x, 6, weight_exponent = 4, size_influence = 0,
                                    prior_strength = 0, min_size = 1)
        expect_identical(b, a)
})

test_that("min_size is respected in every partition", {
        x <- .ssd_fixture()
        set.seed(2)
        m <- .memb_ssd(x, 8, weight_exponent = 4, size_influence = 0, min_size = 25)
        expect_true(all(apply(m, 2, function(z) min(tabulate(z))) >= 25))
        expect_true(all(apply(m, 2, function(z) length(unique(z))) == 8))
})

test_that("min_size is validated", {
        x <- .ssd_fixture(n_tree = 2)
        expect_error(.memb_ssd(x, 8, min_size = 51), "exceeds the number of units")
        expect_error(.memb_ssd(x, 8, min_size = 2.5), "whole number")
        expect_error(.memb_ssd(x, 8, min_size = 0),   "whole number")
})

test_that("size_influence > 0 gives more balanced pieces", {
        x <- .ssd_fixture()
        set.seed(3); a <- .memb_ssd(x, 6, weight_exponent = 1, size_influence = 0)
        set.seed(3); b <- .memb_ssd(x, 6, weight_exponent = 1, size_influence = 3)
        small <- function(m) mean(apply(m, 2, function(z) min(tabulate(z))))
        expect_gt(small(b), small(a))
})

test_that("prior_strength moves cuts onto prior boundaries", {
        x <- .ssd_fixture()
        between <- x$prior[x$ep[, 1]] != x$prior[x$ep[, 2]]
        share <- function(m) mean(apply(m, 2, function(z) {
                cut <- z[x$ep[, 1]] != z[x$ep[, 2]]
                sum(cut & between) / sum(cut)
        }))
        set.seed(4); a <- .memb_ssd(x, 4, weight_exponent = 1, size_influence = 0, prior_strength = 0)
        set.seed(4); b <- .memb_ssd(x, 4, weight_exponent = 1, size_influence = 0, prior_strength = 3)
        expect_gt(share(b), share(a))
})

test_that(".refine_regions never shrinks a region below min_size", {
        x <- .ssd_fixture(n_tree = 5)
        set.seed(5)
        z <- .memb_ssd(x, 6, weight_exponent = 4, size_influence = 0, min_size = 20)[, 1]
        X <- igraph::graph_attr(x$g, "env")
        r <- .refine_regions(z, x$g, X, min_size = 20)
        expect_true(min(tabulate(r$z)) >= 20)
        expect_identical(.refine_regions(z, x$g, X)$z,
                         .refine_regions(z, x$g, X, min_size = 1L)$z)
})
