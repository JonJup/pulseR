cw_fixture <- function(h_p, h_e, ...) {
        n <- length(h_p)
        crosswalk_types(cw_regions(cw_onehot(h_p)), cw_ext(seq_len(n), h_e), "type",
                        units = seq_len(n), ...)
}
renders <- function(p) {
        expect_no_error(if (inherits(p, "patchwork")) print(p) else ggplot2::ggplot_build(p))
        invisible(p)
}

test_that("matrix order puts 1:1 relations on the diagonal", {
        .local_seed(1)
        h  <- sample(1:6, 300, TRUE)
        d  <- .cw_matrix_data(cw_fixture(h, paste0("E", h)), "crisp", "auto", 40L, 28L)$df
        expect_equal(nrow(d), 6L)
        expect_equal(as.integer(d$row), as.integer(d$col))
        expect_true(all(d$cond == 1))
})

test_that("matrix lumps small types into 'other' and keeps the total", {
        h_p <- rep(1:50, times = 50:1)
        cw  <- cw_fixture(h_p, ifelse(h_p %% 2 == 0, "even", "odd"))
        md  <- .cw_matrix_data(cw, "crisp", "auto", 40L, 28L)
        expect_equal(md$n_lumped[["pulse"]], 10L)
        expect_true("other (10 types)" %in% levels(md$df$row))
        expect_equal(sum(md$df$mass), 1)
        expect_equal(md$given, "pulse")
})

test_that("composition bars sum to one and respect top_n", {
        .local_seed(2)
        cw <- cw_fixture(sample(1:8, 400, TRUE), sample(letters[1:3], 400, TRUE))
        cd <- .cw_composition_data(cw, "crisp", "external", 3L, 40L)
        expect_equal(as.vector(tapply(cd$seg$share, cd$seg$bar, sum)), rep(1, 3))
        expect_true(all(tapply(!cd$seg$is_other, cd$seg$bar, sum) <= 3))
        expect_true(all(cd$seg$partner[cd$seg$is_other] == "other"))
})

test_that("map classes combine match and certainty; unmatched geometry is 'none'", {
        U  <- rbind(c(1, 0), c(1, 0), c(0.55, 0.45), c(0.55, 0.45))
        cw <- crosswalk_types(cw_regions(U), cw_ext(1:4, c("a", "b", "a", "b")),
                              "type", units = 1:4)
        g  <- make_grid_sf(5L, 1L)
        d  <- .cw_map_data(cw, g, NULL, 0.5)
        expect_equal(as.character(d$class),
                     c("agree_conf", "disagree_conf", "agree_unc", "disagree_unc", "none"))
})

test_that("certainty panel: chance rate is the independence expectation", {
        cw <- cw_fixture(rep(1:2, each = 50), rep(c("a", "b"), 50))
        d  <- .cw_certainty_data(cw)
        expect_equal(d$chance, 0.5)
        expect_equal(d$overall, 0.5)
})

test_that("plane places identity, refinement and independence correctly", {
        h_p <- rep(1:6, each = 20)
        same <- cw_fixture(h_p, paste0("E", h_p))
        refi <- cw_fixture(h_p, c("A", "B", "C")[ceiling(h_p / 2)])
        .local_seed(3)
        indp <- cw_fixture(sample(1:4, 4000, TRUE), sample(letters[1:4], 4000, TRUE))
        pts  <- .cw_plane_data(list(same = same, refines = refi, indep = indp), "crisp")$points
        expect_equal(unlist(pts[1, c("x", "y")], use.names = FALSE), c(0, 0))
        expect_equal(unlist(pts[2, c("x", "y")], use.names = FALSE), c(0, log(2) / log(6)))
        expect_gt(min(unlist(pts[3, c("x", "y")])), 0.99)
})

test_that("plane draws within-region satellites", {
        .local_seed(4)
        h  <- c(sample(1:2, 40, TRUE), sample(3:5, 60, TRUE))
        cw <- crosswalk_types(cw_local(cw_onehot(h, 5), c(1, 1, 2, 2, 2)),
                              cw_ext(1:100, sample(c("x", "y"), 100, TRUE)), "type",
                              units = 1:100, local_scope = "within_region")
        expect_equal(nrow(.cw_plane_data(cw, "crisp")$satellites), 2L)
})

test_that("every figure renders", {
        skip_if_not_installed("ggplot2")
        .local_seed(5)
        g  <- make_grid_sf(6L, 6L)
        h  <- ifelse(g$temp > stats::median(g$temp), 1L, 2L)
        h[1:4] <- 3L
        g$type <- ifelse(g$prec > stats::median(g$prec), "wet", "dry")
        cw <- crosswalk_types(cw_regions(cw_onehot(h)), g, "type", units = g$ID)
        for (w in c("matrix", "composition", "certainty", "plane"))
                renders(plot_crosswalk(cw, w))
        renders(plot_crosswalk(cw, "matrix", basis = "fuzzy", wording = "plain"))
        renders(plot_crosswalk(cw, "composition",
                               palette = c(R1 = "red", R2 = "blue", R3 = "green")))
        renders(plot_crosswalk(cw, "composition", side = "pulse"))
        renders(plot_crosswalk(cw, "composition",
                               palette = c(`1` = "red", `2` = "blue", `3` = "green")))
        renders(plot_crosswalk(cw, "map", geom = g, panel = FALSE))
        renders(plot_crosswalk(cw, "map", geom = g, dissolve = TRUE, panel = FALSE))
        renders(plot_crosswalk(cw, "map", geom = g, fill = "support", panel = FALSE))
        skip_if_not_installed("patchwork")
        renders(plot_crosswalk(cw, "map", geom = g))
})

test_that("min_coverage = 0 silences the partial-coverage warning", {
        obj <- cw_regions(cw_onehot(rep(1:2, 10)))
        expect_warning(crosswalk_types(obj, cw_ext(1:5, rep("a", 5)), "type", units = 1:20),
                       "min_coverage = 0")
        expect_no_warning(crosswalk_types(obj, cw_ext(1:5, rep("a", 5)), "type", units = 1:20,
                                          min_coverage = 0))
})
