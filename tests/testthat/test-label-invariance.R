## Relabelling invariance (round-5 review, section 4.3). Region numbers, piece
## numbers, unit order and the row order of a tuning log are arbitrary, so no
## score or partition may depend on them.

test_that("CH_FIT and R2_FIT ignore region numbers and unit order", {
        set.seed(1)
        X <- matrix(stats::rnorm(120), 60, 2)
        z <- rep(1:3, each = 20)
        a <- .pseudo_F_env(X, z, 3)
        for (p in list(c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L)))
                expect_equal(.pseudo_F_env(X, p[z], 3), a)
        o <- sample.int(60)
        expect_equal(.pseudo_F_env(X[o, ], z[o], 3), a)
})

test_that("Pi, VI and loss ignore piece numbers; Pi columns follow region numbers", {
        set.seed(2)
        N <- 50; M <- 12; R <- 5; k <- 3
        memb <- matrix(sample.int(R, N * M, replace = TRUE), N, M)
        lab  <- sample(rep_len(seq_len(k), N))
        ref  <- .ensemble_stats(memb, lab, k, pi = TRUE, gvi_a = 1)

        # Renumber the pieces of every ensemble column independently.
        memb2 <- apply(memb, 2, function(v) sample.int(R)[v])
        st2   <- .ensemble_stats(memb2, lab, k, pi = TRUE, gvi_a = 1)
        expect_equal(st2$Pi, ref$Pi)
        expect_identical(st2$pi_support, ref$pi_support)
        expect_equal(st2$vi, ref$vi)
        expect_equal(st2$gvi$expected_loss, ref$gvi$expected_loss)

        # Renumber the consensus regions: Pi columns move, nothing else changes.
        p   <- c(3L, 1L, 2L)
        st3 <- .ensemble_stats(memb, p[lab], k, pi = TRUE, gvi_a = 1)
        expect_equal(st3$Pi[, p], ref$Pi)
        expect_equal(st3$vi, ref$vi)
        expect_equal(st3$gvi$expected_loss, ref$gvi$expected_loss)
        expect_equal(.vi_move_gap(st3, p[lab])$gap, .vi_move_gap(ref, lab)$gap)
})

test_that("the tuning winner does not depend on the row order of the log", {
        # Rows 1/2 and 4/6 tie on both voters; the axis rule must break them
        # the same way whatever order the rows arrive in.
        log <- data.frame(regions_per_partition = 30,
                          weight_exponent = c(1, 4, 4, 1, 4, 1),
                          N_tree = 50, final_regions = c(2, 2, 3, 3, 4, 4),
                          CH_FIT = c(120, 120, 90, 150, 60, 150),
                          STAB   = c(0.7, 0.7, 0.9, 0.6, 0.95, 0.6))
        score <- function(l) {
                b <- .derive_score_column(l, "borda")
                l$score <- as.numeric(b); l$score_tb <- attr(b, "tiebreak")
                l
        }
        axes <- c("final_regions", "weight_exponent")
        ref  <- score(log)
        best <- ref[.pick_best(ref), axes]
        for (s in 1:5) {
                set.seed(s)
                o  <- sample.int(nrow(log))
                l2 <- score(log[o, ])
                expect_equal(l2$score, ref$score[o])
                expect_equal(l2$score_tb, ref$score_tb[o])
                expect_equal(l2[.pick_best(l2), axes], best, ignore_attr = TRUE)
        }
})

test_that("compute_consensus() ignores piece numbers and unit order", {
        skip_if_not_installed("salso")
        set.seed(4)
        N <- 90; base <- rep(1:3, each = 30)
        memb <- sapply(1:20, function(t) {
                v <- (base - 1L) * 2L + sample.int(2L, N, TRUE)
                i <- sample.int(N, 6L)
                v[i] <- sample.int(6L, 6L, TRUE)             # some cross-block noise
                v
        })
        fit <- function(m) compute_consensus(m, 3, verbose = FALSE, report_native = FALSE)
        # map[j]: label in `y` of the region labelled j in `x`
        relabel <- function(x, y) tapply(y, x, function(v) v[1])
        a <- fit(memb)

        b   <- fit(apply(memb, 2, function(v) sample.int(6L)[v]))
        expect_equal(.adjusted_rand_index(a$hard_labels, b$hard_labels), 1)
        expect_equal(b$memberships[, relabel(a$hard_labels, b$hard_labels)], a$memberships)

        o  <- sample.int(N)
        cc <- fit(memb[o, ])
        expect_equal(.adjusted_rand_index(cc$hard_labels, a$hard_labels[o]), 1)
        expect_equal(cc$memberships[, relabel(a$hard_labels[o], cc$hard_labels)],
                     a$memberships[o, ])
})
