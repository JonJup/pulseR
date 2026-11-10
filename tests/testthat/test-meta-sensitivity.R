test_that(".weighted_ari equals ARI under equal weights", {
        a <- c(1,1,2,2,3,3); b <- c(1,1,1,2,2,2)
        expect_equal(.weighted_ari(a, b, rep(1, 6)), .adjusted_rand_index(a, b))
})

test_that("sensitivity re-reads selection without refitting", {
        sel <- data.frame(n_meta_types = 2:5, mass_confined = c(0, .05, .15, .25))
        pick <- function(t) max(sel$n_meta_types[sel$mass_confined <= t + 1e-9])
        expect_equal(vapply(c(.1, .2, .3), pick, numeric(1)), c(3, 4, 5))
})