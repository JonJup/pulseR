# A3 (round 5): compute_local_type_distance(method = "empirical") standardised
# `env` twice, the second time with full-length center/scale vectors, so
# distances were wrong whenever local$scaling was set, and recycled across the
# wrong columns when `variables` was a subset.
#
# Known answer without re-implementing the function: the empirical distance
# is invariant to an affine change of units. A model fitted in z-units with
# stored scaling must give the same distances as the same model expressed in
# raw units with scaling = NULL (the Jacobian is a constant factor per type
# and cancels in the column normalisation).

.a3_local <- function(mu, sigma, region, variables, scaling = NULL) {
        models <- lapply(sort(unique(region)), function(r) {
                j <- which(region == r)
                list(parameters = list(
                        mean     = do.call(cbind, mu[j]),
                        variance = list(sigma = array(unlist(sigma[j]),
                                                      c(length(variables), length(variables),
                                                        length(j))))))
        })
        structure(list(models    = models,
                       key       = data.frame(name   = paste0("t", seq_along(mu)),
                                              region = region),
                       variables = variables,
                       scaling   = scaling),
                  class = "pulse_local_types")
}

.a3_fixture <- function() {
        set.seed(11)
        v  <- c("a", "b", "c")
        ctr <- c(a = 100, b = 50, c = -20)
        scl <- c(a = 10, b = 5, c = 2)
        mu_z <- list(c(-1, 0, 0.5), c(1, 0.5, -0.5), c(0, -1, 1))
        S_z  <- lapply(1:3, function(i) {
                A <- matrix(stats::rnorm(9, sd = 0.4), 3); crossprod(A) + diag(0.3, 3)
        })
        reg  <- c(1L, 1L, 2L)
        Dm   <- diag(scl)
        z    <- matrix(stats::rnorm(3000 * 3), ncol = 3)
        env  <- as.data.frame(sweep(sweep(z, 2, scl, "*"), 2, ctr, "+"))
        names(env) <- v
        list(
                env = env,
                z   = .a3_local(mu_z, S_z, reg, v, list(center = ctr, scale = scl)),
                raw = .a3_local(lapply(mu_z, function(m) ctr + scl * m),
                                lapply(S_z, function(S) Dm %*% S %*% Dm), reg, v, NULL))
}

test_that("empirical distance is invariant to units (all variables)", {
        fx <- .a3_fixture()
        d_z   <- compute_local_type_distance(fx$z,   "empirical", env = fx$env, n_sample = 5000)
        d_raw <- compute_local_type_distance(fx$raw, "empirical", env = fx$env, n_sample = 5000)
        expect_equal(as.matrix(d_z), as.matrix(d_raw), tolerance = 1e-6)
        expect_gt(max(d_z), 0.1)                       # not trivially all zero
})

test_that("empirical distance is invariant to units on a variable subset", {
        fx <- .a3_fixture()
        vars  <- c("c", "a")
        d_z   <- expect_silent(compute_local_type_distance(
                fx$z, "empirical", env = fx$env, variables = vars, n_sample = 5000))
        d_raw <- compute_local_type_distance(
                fx$raw, "empirical", env = fx$env, variables = vars, n_sample = 5000)
        expect_equal(as.matrix(d_z), as.matrix(d_raw), tolerance = 1e-6)
        expect_identical(attr(d_z, "variables"), vars)
})

test_that("empirical distance ignores env column order and extra columns", {
        fx   <- .a3_fixture()
        vars <- c("b", "c")
        d1   <- compute_local_type_distance(fx$z, "empirical", env = fx$env,
                                            variables = vars, n_sample = 5000)
        env2 <- cbind(ID = seq_len(nrow(fx$env)), fx$env[, c("c", "a", "b")])
        d2   <- compute_local_type_distance(fx$z, "empirical", env = env2,
                                            variables = vars, n_sample = 5000)
        d3   <- compute_local_type_distance(fx$z, "empirical", env = as.matrix(env2),
                                            variables = vars, n_sample = 5000)
        expect_equal(as.matrix(d1), as.matrix(d2), tolerance = 1e-12)
        expect_equal(as.matrix(d1), as.matrix(d3), tolerance = 1e-12)
})

test_that("empirical distance names missing env columns", {
        fx <- .a3_fixture()
        expect_error(compute_local_type_distance(fx$z, "empirical",
                                                 env = fx$env[, c("a", "c")]),
                     "lacks variable\\(s\\).*b")
})
