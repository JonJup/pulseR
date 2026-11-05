## Level 2 covariance-model fallback (Section 9: small ecoregions).

.sim_two_regions <- function(n_big = 400, n_small = 8, d = 10, seed = 1) {
        set.seed(seed)
        n   <- n_big + n_small
        X   <- matrix(stats::rnorm(n * d), n, d,
                      dimnames = list(NULL, paste0("v", seq_len(d))))
        X[seq_len(n_big / 2), 1:3] <- X[seq_len(n_big / 2), 1:3] + 3   # two types in big region
        W   <- cbind(c(rep(0.95, n_big), rep(0.05, n_small)),
                     c(rep(0.05, n_big), rep(0.95, n_small)))
        list(env = data.frame(ID = seq_len(n), X), W = W, d = d)
}

test_that("a region with fewer core units than predictors no longer stops", {
        s   <- .sim_two_regions()
        res <- suppressWarnings(get_local_types(s$env, s$W, non_value_cols = "ID",
                                                n_local_types = 2:4))
        expect_s3_class(res, "pulse_local_types")
        expect_equal(unname(rowSums(res$memberships)), rep(1, nrow(s$env)), tolerance = 1e-8)
        ## the small region must have fallen back below full covariance
        expect_false(res$model_names[["Region2"]] %in% c("VVV", "EEE"))
})

test_that("selected fits never have more parameters than n_eff", {
        s   <- .sim_two_regions()
        res <- suppressWarnings(get_local_types(s$env, s$W, non_value_cols = "ID",
                                                n_local_types = 2:4))
        for (r in seq_along(res$model_names)) {
                npar  <- mclust::nMclustParams(res$model_names[[r]], s$d, res$n_local_types[[r]])
                n_eff <- sum(s$W[s$W[, r] > 1e-6, r])
                expect_lt(npar, n_eff)
        }
})

test_that("BIC keeps its G x region shape and equals the max over models", {
        s   <- .sim_two_regions()
        res <- suppressWarnings(get_local_types(s$env, s$W, non_value_cols = "ID",
                                                n_local_types = 2:4))
        expect_identical(dim(res$BIC), c(3L, 2L))
        expect_identical(dim(res$BIC_by_model)[c(1L, 3L)], c(3L, 2L))
        expect_equal(res$BIC[, 1], apply(res$BIC_by_model[, , 1], 1, max, na.rm = TRUE))
})

test_that("a single model name still works (one-family mode)", {
        s   <- .sim_two_regions(n_small = 60)
        res <- suppressWarnings(get_local_types(s$env, s$W, non_value_cols = "ID",
                                                n_local_types = 2:3, model_name = "VVI"))
        expect_true(all(res$model_names == "VVI"))
        expect_identical(dimnames(res$BIC_by_model)$model, "VVI")
})

test_that("bad model names fail early and clearly", {
        s <- .sim_two_regions()
        expect_error(get_local_types(s$env, s$W, non_value_cols = "ID", model_name = "XYZ"),
                     "Unknown mclust model name")
})

test_that("one predictor maps multivariate names to E/V", {
        s   <- .sim_two_regions(n_small = 60)
        res <- suppressWarnings(get_local_types(s$env, s$W, variables = "v1",
                                                n_local_types = 2:3))
        expect_true(all(res$model_names %in% c("E", "V")))
})
