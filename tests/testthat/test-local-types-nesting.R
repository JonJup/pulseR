test_that("crisp local types nest in the Level 1 hard region", {
        skip_if_not_installed("mclust")
        set.seed(1)
        n   <- 200
        env <- data.frame(ID = seq_len(n), a = rnorm(n), b = rnorm(n), c = rnorm(n))
        W   <- matrix(runif(n * 2), ncol = 2); W <- W / rowSums(W)
        h   <- max.col(W)
        flip <- 1:15
        h[flip] <- 3L - h[flip]          # mimic contiguity repair overriding argmax
        reg <- list(memberships = W, hard_labels = h)
        
        res <- get_local_types(env, regions = reg, non_value_cols = "ID",
                               n_local_types = 2:3)
        expect_identical(res$region_labels, as.integer(h))
        expect_identical(res$key$region[res$hard_labels], as.integer(h))
        expect_equal(res$certainty,
                     res$memberships[cbind(seq_len(n), res$hard_labels)])
})