test_that("a homogeneous region selects one local type", {
        set.seed(1)
        n   <- 200
        env <- data.frame(x1 = rnorm(2 * n), x2 = rnorm(2 * n))
        env$x1[(n + 1):(2 * n)] <- env$x1[(n + 1):(2 * n)] + 6   # region 2 differs
        W   <- cbind(rep(1:0, each = n), rep(0:1, each = n))
        loc <- get_local_types(env, W, n_local_types = 1:3)
        expect_equal(unname(loc$n_local_types), c(1L, 1L))
})