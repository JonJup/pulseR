# tests/testthat/test-moves-null.R
toy_graph <- function(n = 12, seed = 1) {
        set.seed(seed)
        g  <- igraph::make_lattice(c(n, n))
        xy <- expand.grid(x = seq_len(n), y = seq_len(n))
        env <- cbind(e1 = (xy$x > n / 3) + (xy$x > 2 * n / 3) + stats::rnorm(n^2, 0, 0.3),
                     e2 = stats::rnorm(n^2))
        ep <- igraph::ends(g, igraph::E(g), names = FALSE)
        igraph::E(g)$eid    <- seq_len(igraph::ecount(g))
        igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], ] - env[ep[, 2], ])^2)) + 1e-6
        list(g = g, env = env)
}

test_that("VI moves lower the loss monotonically and stop at a single-move optimum", {
        set.seed(2)
        base  <- rep(1:3, length.out = 300)
        memb  <- sapply(1:30, function(t) (base - 1L) * 4L + sample.int(4L, 300, TRUE))
        start <- base; flip <- runif(300) < 0.2
        start[flip] <- sample.int(3L, sum(flip), TRUE)
        a  <- 1.2
        l0 <- salso::VI(t(memb), start, a = a)
        mv <- .vi_moves(memb, start, k = 3L, a = a, max_rounds = 50L)
        expect_true(mv$converged)
        expect_equal(mv$blocked, 0L)
        expect_lt(salso::VI(t(memb), mv$labels, a = a), l0)
        expect_equal(sort(unique(mv$labels)), 1:3)
        gp <- .vi_move_gap(.ensemble_stats(memb, mv$labels, 3L, gvi_a = a), mv$labels)
        expect_true(all(gp$gap > -1e-9))
})

test_that("sampled path ends with no improving single move", {
        tg <- toy_graph()
        tr <- sample_spanning_trees(tg$g, n = 30, seed = 1, verbose = FALSE)
        en <- compute_ensemble_memberships(tr, N_tree = 30, graph = tg$g,
                                           regions_per_partition = 8, verbose = FALSE)
        set.seed(4)
        r <- suppressWarnings(cluster_consensus(en, 3, verbose = FALSE, report_native = FALSE,
                                                large_n_threshold = 100, sample_size = 60,
                                                refine_iter = 50L))
        expect_true(r$refine_info$converged)
        expect_true(all(r$loss_gap > -1e-9))
})

test_that("voters controls the Borda aggregate", {
        tg <- toy_graph()
        tr <- sample_spanning_trees(tg$g, n = 30, seed = 1, verbose = FALSE)
        expect_error(tune_regions(tg$g, tr, N_tree = 30, regions_per_partition = 8,
                                  final_regions = 2:3, env = tg$env, voters = "NOPE",
                                  verbose = FALSE), "voters")
        tu <- tune_regions(tg$g, tr, N_tree = 30, regions_per_partition = 8, final_regions = 2:4,
                           env = tg$env, stability_B = 0, voters = "CH_ENV",
                           verbose = FALSE, seed = 1)
        expect_equal(tu$best_final_regions,
                     tu$tuning_log$final_regions[which.max(tu$tuning_log$CH_ENV)])
})

test_that("null baseline returns one summary row per k and criterion", {
        tg <- toy_graph(n = 10)
        nb <- compute_null_baseline(tg$g, tg$env, final_regions = 2:3, regions_per_partition = 8,
                                    N_tree = 20, B = 2, stability_B = 0, seed = 1,
                                    verbose = FALSE, target_edge_turnover = 0.1)
        expect_s3_class(nb, "null_baseline")
        expect_equal(nrow(nb$summary), 2L * length(.CVI_DIRECTION))
        expect_equal(nrow(nb$null), 4L)
        expect_true(all(is.na(nb$summary$observed)))
})
