make_constant_case <- function(seed = 1) {
        set.seed(seed)
        n <- 300
        w1 <- c(rep(1, 200), rep(0, 100)); W <- cbind(w1, 1 - w1)
        geo <- c(stats::rbinom(200, 1, 0.4), rep(0, 100))      # constant in region 2
        env <- data.frame(a = c(stats::rnorm(100), stats::rnorm(100, 4), stats::rnorm(100)),
                          b = stats::rnorm(n), geo = geo)
        list(env = env, W = W)
}

test_that(".global_scale_prior() matches mclust::defaultPrior() on the same data", {
        set.seed(2)
        X <- matrix(stats::rnorm(200), 50, 4)
        V <- stats::cov(X)
        for (m in c("VVV", "EEE", "VVI", "VII")) {
                ref <- mclust::defaultPrior(X, G = 3, modelName = m)$scale
                got <- .global_scale_prior(mclust::priorControl(), V, m, 3)$scale
                expect_equal(got, ref, ignore_attr = TRUE)
        }
        expect_identical(.global_scale_prior(mclust::priorControl(scale = 1), V, "VVV", 3)$scale, 1)
})

test_that("a predictor constant in one region is kept and regularised by default", {
        cc <- make_constant_case()
        expect_message(
                f <- get_local_types(cc$env, cc$W, n_local_types = 2:3, model_name = "VVV"),
                "covariance of all units")
        expect_true("geo" %in% f$variables)
        expect_identical(unname(f$constant["geo", ]), c(FALSE, TRUE))
        expect_true(all(is.finite(f$BIC_by_model[, "VVV", 2][!is.na(f$BIC_by_model[, "VVV", 2])])))
        expect_false(all(is.na(f$BIC_by_model[, "VVV", 2])))
        expect_equal(f$n_local_types[["Region1"]], 2L)
})

test_that("drop_constant = TRUE and prior = NULL drop with a warning", {
        cc <- make_constant_case()
        expect_warning(f1 <- get_local_types(cc$env, cc$W, n_local_types = 2,
                                             drop_constant = TRUE), "ALL regions")
        expect_false("geo" %in% f1$variables)
        expect_warning(f2 <- get_local_types(cc$env, cc$W, n_local_types = 2, prior = NULL),
                       "cannot regularise")
        expect_false("geo" %in% f2$variables)
        env <- cc$env; env$k <- 1
        expect_warning(f3 <- suppressMessages(get_local_types(env, cc$W, n_local_types = 2)),
                       "constant across all units")
        expect_false("k" %in% f3$variables)
})
