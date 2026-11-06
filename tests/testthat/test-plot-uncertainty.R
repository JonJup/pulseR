mock_unc_levels <- function() {
        reg <- structure(list(memberships   = rbind(c(0.9, 0.1, 0), c(0.5, 0.3, 0.2),
                                                    c(NA, NA, NA), c(0.2, 0.7, 0.1)),
                              hard_clusters = c(1L, 2L, NA, 2L),   # row 2: repaired, hard != argmax
                              fit = list()),
                         class = "pulse_regions")
        ## Region 1 has two local types, region 2 has one.
        U_loc <- rbind(c(0.6, 0.4, 0), c(0.9, 0.1, 0), c(0, 0, 1), c(0, 0, 1))
        loc <- structure(list(memberships  = U_loc,
                              hard_classes = max.col(U_loc, ties.method = "first"),
                              certainty    = apply(U_loc, 1, max),
                              key = data.frame(column = 1:3, region = c(1L, 1L, 2L),
                                               local_type = c(1L, 2L, 1L))),
                         class = "pulse_local_types")
        U_mos <- rbind(c(0.7, 0.2, 0.1), c(0.4, 0.35, 0.25))
        mos <- structure(list(memberships = U_mos, hard_clusters = c(1L, 1L),
                              best_k = 3L, valid = c(TRUE, FALSE, TRUE, FALSE)),
                         class = "mosaic_types")
        list(reg = reg, loc = loc, mos = mos)
}

test_that(".top_two() matches a per-row sort, across row blocks and NA rows", {
        set.seed(1)
        U <- matrix(stats::runif(60), 20, 3); U[5, ] <- NA
        t2 <- .top_two(U, block = 9)                 # 3 rows per block
        ref <- t(apply(U, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
        expect_equal(t2$first, ref[, 1])
        expect_equal(t2$second, ref[, 2])
        expect_equal(.top_two(matrix(c(0.3, 1), 2))$second, c(0, 0))
})

test_that("per-level extraction uses the package's certainty and the region's k", {
        m <- mock_unc_levels()
        r <- .level_uncertainty(m$reg, "regions")
        expect_equal(r$certainty, c(0.9, 0.3, NA, 0.7))       # own-region Pi, not row max
        expect_equal(r$margin, c(0.8, 0.2, NA, 0.5))
        l <- .level_uncertainty(m$loc, "local")
        expect_equal(l$k, c(2L, 2L, 1L, 1L))
        expect_equal(l$margin, c(0.2, 0.8, 1, 1))
        expect_error(.level_uncertainty(m$loc, "regions"), "pulse_regions")
})

test_that("NA and single-candidate units are dropped and reported", {
        m  <- mock_unc_levels()
        df <- .uncertainty_df(list(regions = m$reg, local = m$loc), drop_trivial = TRUE)
        expect_equal(as.vector(table(df$level)), c(3L, 2L))
        expect_match(paste(attr(df, "notes"), collapse = " "), "single candidate")
        df2 <- .uncertainty_df(list(local = m$loc), drop_trivial = FALSE)
        expect_equal(nrow(df2), 4L)
        expect_error(.uncertainty_df(list(regions = NULL), TRUE), "at least one")
})

test_that("plot_certainty() and plot_margin() build for every option", {
        skip_if_not_installed("ggplot2")
        m <- mock_unc_levels()
        for (tp in c("histogram", "ecdf")) {
                p1 <- plot_certainty(regions = m$reg, local = m$loc, mosaic = m$mos, type = tp)
                p2 <- plot_certainty(regions = m$reg, mosaic = m$mos, rescale = TRUE, type = tp)
                p3 <- plot_margin(regions = m$reg, local = m$loc, mosaic = m$mos, type = tp)
                p4 <- plot_margin(mosaic = m$mos, threshold = NULL, type = tp)
                for (p in list(p1, p2, p3, p4)) {
                        expect_s3_class(p, "ggplot")
                        expect_silent(ggplot2::ggplot_build(p))
                }
        }
        expect_equal(plot_certainty(mosaic = m$mos, rescale = TRUE)$data$value,
                     c(0.55, 0.1))                    # (0.7 - 1/3) / (2/3), (0.4 - 1/3) / (2/3)
})

test_that("crisp memberships trigger a warning", {
        skip_if_not_installed("ggplot2")
        mos <- structure(list(memberships = diag(3), hard_clusters = 1:3, best_k = 3L,
                              valid = rep(TRUE, 3)), class = "mosaic_types")
        expect_warning(plot_certainty(mosaic = mos), "crisp")
})
