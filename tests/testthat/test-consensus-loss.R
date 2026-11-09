# tests/testthat/test-consensus-loss.R
make_memb <- function(N = 300, M = 40, R = 12, seed = 1) {
        set.seed(seed)
        base <- rep(1:3, length.out = N)[order(runif(N))]
        sapply(seq_len(M), function(t) {
                l <- (base - 1L) * 4L + sample.int(4L, N, replace = TRUE)
                flip <- runif(N) < 0.1
                l[flip] <- sample.int(R, sum(flip), replace = TRUE)
                l
        })
}

test_that("salso's a is kept and expected_loss is salso's generalised VI", {
        memb <- make_memb()
        r <- compute_consensus(memb, final_regions = 3, verbose = FALSE, report_native = FALSE)
        expect_true(is.finite(r$salso_a) && r$salso_a >= 0 && r$salso_a <= 2)
        ref <- salso::VI(t(memb), r$hard_labels, a = r$salso_a)
        expect_equal(r$expected_loss, ref, tolerance = 1e-10)
})

test_that("loss_gap equals the brute-force loss change of a single move", {
        memb <- make_memb()
        r <- compute_consensus(memb, final_regions = 3, verbose = FALSE, report_native = FALSE)
        N <- nrow(memb); ref <- r$expected_loss
        for (i in c(1L, 50L, which.min(r$loss_gap))) {
                h <- r$hard_labels; h[i] <- r$loss_alt[i]
                bf <- N * (salso::VI(t(memb), h, a = r$salso_a) - ref)
                expect_equal(r$loss_gap[i], bf, tolerance = 1e-8)
        }
})

test_that("label_step covers every unit and follows the sampled path", {
        memb <- make_memb(N = 600)
        r <- suppressWarnings(compute_consensus(memb, final_regions = 3, verbose = FALSE, report_native = FALSE,
                                                large_n_threshold = 400, sample_size = 200))
        expect_s3_class(r$label_step, "factor")
        expect_false(anyNA(r$label_step))
        expect_true(sum(r$label_step == "assign") > 0)
        expect_lte(sum(r$label_step %in% c("salso", "merge")), 200)
})
