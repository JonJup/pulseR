## A11: smooth_signatures() must not return all-NA signatures silently.

.a11_sig <- function() {
        S <- rbind(c(0.50, 0.25, 0.25), c(0.10, 0.10, 0.80))
        colnames(S) <- c("m1:m1", "m1:m2", "m2:m2")
        structure(S, class_labels = c("m1", "m2"))
}

test_that("all-NA region_of (meta-type input) stops and points to class_similarity", {
        expect_error(smooth_signatures(.a11_sig(), region_of = c(NA_integer_, NA_integer_)),
                     "all NA.*class_similarity")
})

test_that("meta-type mosaic input gives the same error", {
        U  <- rbind(c(0.7, 0.3), c(0.2, 0.8), c(0.5, 0.5))
        colnames(U) <- c("m1", "m2")
        mi <- .new_mosaic_input(U, colnames(U), rep(NA_integer_, 2), "meta",
                                normalise = FALSE, drop_empty = TRUE)
        expect_error(smooth_signatures(.a11_sig(), region_of = mi$region_of), "all NA")
})

test_that("a partly NA region_of treats NA classes as similar only to themselves", {
        S   <- .a11_sig()
        got <- smooth_signatures(S, region_of = c(1L, NA_integer_))
        expect_false(anyNA(got))
        expect_equal(unclass(got)[, , drop = FALSE],
                     unclass(smooth_signatures(S, class_similarity = diag(2))),
                     ignore_attr = TRUE)
})

test_that("bad region_of length or class_similarity stops with a clear message", {
        S <- .a11_sig()
        expect_error(smooth_signatures(S, region_of = 1:3), "length 3 .* 2 classes")
        Sc <- matrix(c(1, NA, NA, 1), 2)
        expect_error(smooth_signatures(S, class_similarity = Sc), "2 x 2 matrix without NA")
        expect_error(smooth_signatures(S, class_similarity = diag(3)), "2 x 2")
})

test_that("regional smoothing is unchanged", {
        S  <- .a11_sig()
        Sc <- matrix(c(1, 0.5, 0.5, 1), 2)
        expect_equal(smooth_signatures(S, region_of = c(1L, 1L)),
                     smooth_signatures(S, class_similarity = Sc))
})
