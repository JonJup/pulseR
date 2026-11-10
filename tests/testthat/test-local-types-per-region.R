make_two_region_case <- function(seed = 1) {
        set.seed(seed)
        w1  <- rep(c(1, 0), each = 150)
        env <- data.frame(a = stats::rnorm(300), b = stats::rnorm(300))
        list(env = env, W = cbind(w1, 1 - w1))
}

test_that("a list n_local_types fixes G per region", {
        s <- make_two_region_case()
        f <- get_local_types(s$env, s$W, n_local_types = list(2, 4), model_name = "VVI")
        expect_identical(unname(f$n_local_types), c(2L, 4L))
        expect_identical(ncol(f$memberships), 6L)
        expect_identical(rownames(f$BIC), c("2", "4"))
        expect_true(is.na(f$BIC["4", "Region1"]))
})

test_that("list elements can mix fixed and candidate G", {
        s <- make_two_region_case()
        f <- get_local_types(s$env, s$W, n_local_types = list(3, 2:4), model_name = "VVI")
        expect_identical(f$n_local_types[["Region1"]], 3L)
        expect_true(f$n_local_types[["Region2"]] %in% 2:4)
})

test_that("bad per-region n_local_types errors or warns", {
        s <- make_two_region_case()
        expect_error(get_local_types(s$env, s$W, n_local_types = list(2)),
                     "length 1 but there are 2")
        expect_error(get_local_types(s$env, s$W, n_local_types = list(2, 200)),
                     "Region 2 supports no requested G")
        expect_error(get_local_types(s$env, s$W, n_local_types = list(2, 1.5)),
                     "n_local_types\\[\\[2\\]\\]")
        expect_warning(get_local_types(s$env, s$W, n_local_types = list(2, c(3, 200)),
                                       model_name = "VVI"),
                       "Region 2: dropped G = 200")
})