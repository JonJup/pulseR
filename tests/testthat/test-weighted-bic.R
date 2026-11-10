## Pins the scale of the likelihood term in .weighted_bic().
## mclust's own `fit$loglik` changed scale in 6.1.x (now divided by mean(w)),
## so pulseR recomputes sum(w * log f) itself via .weighted_loglik().

wb_fixture <- function() {
        set.seed(1)
        X <- rbind(matrix(rnorm(200), ncol = 2), matrix(rnorm(200, 3), ncol = 2))
        list(X  = X,
             w  = runif(200, 0.2, 1),
             z0 = mclust::unmap(mclust::hclass(mclust::hc(X, "VVV"), 2)))
}

test_that(".weighted_loglik halves when weights halve, on any mclust version", {
        skip_if_not_installed("mclust")
        f  <- wb_fixture()
        ll <- function(w) .weighted_loglik(
                mclust::me.weighted(f$X, "VVV", f$z0, weights = w, prior = NULL),
                f$X, w)
        expect_equal(ll(f$w / 2) / ll(f$w), 0.5, tolerance = 1e-6)
})

test_that(".weighted_loglik equals sum(w * log f) at the fitted parameters", {
        skip_if_not_installed("mclust")
        f   <- wb_fixture()
        fit <- mclust::me.weighted(f$X, "VVV", f$z0, weights = f$w, prior = NULL)
        ld  <- mclust::dens(data = f$X, modelName = "VVV",
                            parameters = fit$parameters, logarithm = TRUE)
        expect_equal(.weighted_loglik(fit, f$X, f$w), sum(f$w * ld))
})

test_that(".weighted_bic = 2*loglik - npar*log(n_eff), NA on bad input", {
        npar <- mclust::nMclustParams("VVV", d = 2, G = 3)
        expect_equal(.weighted_bic(-100, "VVV", 2, 3, 50),
                     -200 - npar * log(50))
        expect_true(is.na(.weighted_bic(NA_real_, "VVV", 2, 3, 50)))
        expect_true(is.na(.weighted_bic(-100, "VVV", 2, 3, 1)))
})