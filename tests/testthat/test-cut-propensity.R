test_that("zero weights are floored at the smallest positive weight", {
        set.seed(1)
        w <- stats::runif(500, 0.1, 2); s <- stats::runif(500, 1, 250)
        w0 <- w; w0[1:5] <- 0
        wf <- w0; wf[1:5] <- min(w0[w0 > 0])
        expect_equal(.cut_log_propensity(w0, s), .cut_log_propensity(wf, s))
        # 5 zeros among 500 edges used to raise alpha 6.8-fold
        a  <- .cut_log_propensity(w,  s)$alpha
        a0 <- .cut_log_propensity(w0, s)$alpha
        expect_lt(abs(a0 / a - 1), 0.1)
        expect_true(all(is.finite(.cut_log_propensity(w0, s)$lp)))
})

test_that("all-zero weights give uniform environmental propensity", {
        s <- c(1, 2, 3, 4)
        pr <- .cut_log_propensity(rep(0, 4), s)
        expect_true(all(is.finite(pr$lp)))
        expect_gt(pr$alpha, 0)
        expect_equal(pr$unit, log(2))
})

test_that("the ensemble runs on a graph with zero weights", {
        set.seed(2)
        g <- igraph::make_lattice(c(6, 6))
        w <- stats::runif(igraph::ecount(g)); w[1:4] <- 0
        igraph::E(g)$weight <- w
        tr <- sample_spanning_trees(g, n = 5, seed = 1, verbose = FALSE)
        em <- compute_ensemble(tr, N_tree = 5, graph = g,
                                           regions_per_partition = 4, verbose = FALSE)
        expect_true(all(apply(em$memberships, 2, function(l) length(unique(l))) == 4L))
})
