## get_meta_types(): "gap" vs "crossing" on a known dendrogram
## (round-5 review, section 4.2).
##
## Six local types on a line, in three tight pairs. Average linkage merges
## at heights 0.1, 0.1, 0.1, 10, 16, so the largest jump (9.9) is the cut
## into the three pairs. Two pairs sit inside one region each, so only the
## two-group cut keeps every meta-type crossing regions.

.rule_fixture <- function() {
        pos <- c(0, 0.1, 10, 10.1, 21, 21.1)
        list(d = stats::dist(pos),
             local = structure(
                     list(memberships = diag(6),
                          key = data.frame(name = paste0("t", 1:6),
                                           region = c(1, 1, 2, 2, 1, 2))),
                     class = "pulse_local_types"))
}

test_that("gap rule cuts at the largest jump in merge height", {
        x <- .rule_fixture()
        m <- get_meta_types(x$local, n_meta_types = 2:5, distance = x$d, rule = "gap")
        expect_identical(m$n_meta_types, 3L)
        expect_equal(m$selection$height_gap, c(6, 9.9, 0, 0), tolerance = 1e-10)
        expect_equal(m$gap_ratio, 9.9 / 3, tolerance = 1e-10)      # median gap = 3
        expect_identical(m$hard_labels, rep(1:3, each = 2))
        expect_equal(unname(rowSums(m$memberships)), rep(1, 6))
        expect_match(m$rule, "^gap")
})

test_that("crossing rule takes the finest cut that still crosses regions", {
        x <- .rule_fixture()
        m <- get_meta_types(x$local, n_meta_types = 2:5, distance = x$d, rule = "crossing")
        expect_identical(m$n_meta_types, 2L)
        expect_equal(m$selection$mass_confined, c(0, 4 / 6, 1, 1))
        expect_identical(m$selection$passes, c(TRUE, FALSE, FALSE, FALSE))
        expect_true(attr(m$sensitivity, "stable"))
        expect_match(m$rule, "^crossing")
})

test_that("gap rule warns when the pick is the largest candidate", {
        x <- .rule_fixture()
        expect_warning(get_meta_types(x$local, n_meta_types = 2:3, distance = x$d, rule = "gap"),
                       "largest candidate")
})
