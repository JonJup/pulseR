test_that(".confidence_vec() defaults to Pi certainty, not cohesion", {
        U <- rbind(c(0.9, 0.1), c(0.6, 0.4), c(NA, NA), c(0.2, 0.8))
        A <- rbind(c(0.3, 0.1), c(0.2, 0.1), c(0.1, 0.1), c(0.1, 0.7))
        x <- structure(list(memberships = U, hard_clusters = c(1L, 1L, 1L, 2L),
                            fit = list(cohesion = A)), class = "pulse_regions")
        expect_equal(.confidence_vec(x), c(0.9, 0.6, NA, 0.8))
        expect_equal(.confidence_vec(x, "margin"), c(0.8, 0.2, NA, 0.6))
        expect_equal(.confidence_vec(x, "typicality"), c(0.3, 0.2, 0.1, 0.7))
        x$hard_clusters[2] <- NA
        expect_true(is.na(.confidence_vec(x)[2]))
        expect_identical(formals(plot_regions)$confidence,
                         quote(c("certainty", "margin", "typicality")))
})

test_that("joint local-type confidence counts Level 1 once", {
        loc <- list(certainty = c(0.8, 0.8), certainty_raw = c(0.4, 0.8))
        expect_equal(.local_confidence(loc, NULL, "joint"), c(0.4, 0.8))
})
