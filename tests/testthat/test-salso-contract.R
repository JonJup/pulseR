test_that("installed salso reports `a` for a k-targeted VI search", {
        set.seed(1)
        memb <- sapply(1:10, function(t) rep(1:3, each = 10) + 3L * sample.int(2L, 30, TRUE))
        est <- salso::salso(t(memb), loss = salso::VI(a = list(nClusters = 3)), nRuns = 2)
        a <- attr(est, "info")$a
        expect_length(a, 1L)
        expect_true(is.finite(a))
})

test_that(".salso_k() stops clearly when salso reports no `a`", {
        memb <- matrix(rep(1:2, each = 10), 20, 4)
        local_mocked_bindings(
                salso = function(x, ...) structure(rep(1:2, each = 10), info = list()),
                .package = "salso")
        expect_error(.salso_k(memb, 2L), "did not report the search parameter")
})
