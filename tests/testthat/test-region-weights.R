test_that(".region_weights() checks a bare matrix like a pulse_regions", {
        W <- cbind(c(0.9, 0.2, NA), c(0.1, 0.8, NA))
        expect_error(.region_weights(W), "NA in 1 row")
        expect_error(.region_weights(as.data.frame(W)), "NA in 1 row")
        expect_error(.region_weights(cbind(c(1, -0.1), c(0, 1.1))), "non-negative")
        expect_error(.region_weights(matrix(c("a", "b"), 1)), "numeric")
        ok <- W[1:2, ]
        expect_identical(.region_weights(ok), ok)
        reg <- structure(list(memberships = W), class = "pulse_regions")
        expect_error(.region_weights(reg), "pi_support")
})

test_that("get_local_types() stops on NA rows in a bare matrix", {
        set.seed(1)
        env <- data.frame(a = stats::rnorm(60), b = stats::rnorm(60))
        w <- stats::plogis(seq(-3, 3, length.out = 60)); W <- cbind(w, 1 - w)
        W[5, ] <- NA
        expect_error(get_local_types(env, W, n_local_types = 2), "NA in 1 row")
})
