# Every exported function that takes `seed` must leave the caller's RNG
# stream untouched, and give the same result for the same seed.

test_that(".local_seed() restores an existing stream and removes a new one", {
        f <- function(seed) { .local_seed(seed); stats::runif(1) }
        set.seed(5); before <- .Random.seed
        x1 <- f(7)
        expect_identical(.Random.seed, before)
        expect_identical(f(7), x1)
        
        had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (had) old <- get(".Random.seed", envir = .GlobalEnv)
        rm(".Random.seed", envir = .GlobalEnv)
        f(7)
        expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        if (had) assign(".Random.seed", old, envir = .GlobalEnv)
})

test_that(".fit_dmm() and .stability_pam_diss() restore the caller's RNG", {
        set.seed(3)
        X <- matrix(stats::rpois(60 * 4, 5) + 1, 60, 4)
        set.seed(99); before <- .Random.seed
        a <- .fit_dmm(X, 2L, n_start = 2L, max_iter = 20L, seed = 1)
        expect_identical(.Random.seed, before)
        b <- .fit_dmm(X, 2L, n_start = 2L, max_iter = 20L, seed = 1)
        expect_identical(a$cluster, b$cluster)
        
        D <- as.matrix(stats::dist(X))
        hc <- stats::cutree(stats::hclust(stats::as.dist(D)), 2)
        s1 <- .stability_pam_diss(D, hc, 2L, B = 3L, seed = 1)
        expect_identical(.Random.seed, before)
        expect_identical(.stability_pam_diss(D, hc, 2L, B = 3L, seed = 1), s1)
})
