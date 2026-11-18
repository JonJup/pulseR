## A9: env rows are matched to the region weights by ID when `units` is given.

.a9_fixture <- function(n = 120) {
        set.seed(9)
        reg <- rep(1:2, each = n / 2)
        env <- data.frame(ID   = sprintf("u%03d", seq_len(n)),
                          temp = stats::rnorm(n, mean = c(0, 4)[reg]),
                          prec = stats::rnorm(n, mean = rep(c(0, 3), n / 2)))
        W <- cbind(ifelse(reg == 1, 0.9, 0.1), ifelse(reg == 1, 0.1, 0.9))
        list(env = env, W = W)
}

.a9_fit <- function(...) suppressMessages(
        get_local_types(..., n_local_types = 2, model_name = "VVV"))

test_that("a shuffled env with `units` gives the same local types as the ordered one", {
        f   <- .a9_fixture()
        ref <- .a9_fit(f$env, f$W, units = f$env$ID)
        set.seed(1)
        shuf <- f$env[sample(nrow(f$env)), ]
        got  <- .a9_fit(shuf, f$W, units = f$env$ID)
        expect_equal(got$memberships, ref$memberships)
        expect_identical(got$hard_labels, ref$hard_labels)
})

test_that("extra env rows are ignored and numeric IDs match", {
        f   <- .a9_fixture()
        f$env$ID <- seq_len(nrow(f$env))
        ref <- .a9_fit(f$env, f$W, units = f$env$ID)
        extra <- f$env[c(seq_len(nrow(f$env)), 1:5), ]
        extra$ID[seq(nrow(f$env) + 1, nrow(extra))] <- 1000 + 1:5
        got <- .a9_fit(extra[nrow(extra):1, ], f$W, units = f$env$ID)
        expect_equal(got$memberships, ref$memberships)
})

test_that("`units` accepts a graph with polygon_id", {
        f <- .a9_fixture()
        g <- igraph::make_ring(nrow(f$env))
        igraph::V(g)$polygon_id <- f$env$ID
        ref <- .a9_fit(f$env, f$W, units = f$env$ID)
        got <- .a9_fit(f$env[rev(seq_len(nrow(f$env))), ], f$W, units = g)
        expect_equal(got$memberships, ref$memberships)
})

test_that("ID problems stop with a message that names them", {
        f <- .a9_fixture()
        expect_error(.a9_fit(f$env[-1, ], f$W, units = f$env$ID), "no environmental row")
        dup <- f$env; dup$ID[2] <- dup$ID[1]
        expect_error(.a9_fit(dup, f$W, units = f$env$ID), "duplicated ID")
        expect_error(.a9_fit(f$env[, -1], f$W, units = f$env$ID), "no `ID`/`id` column")
        expect_error(.a9_fit(f$env, f$W, units = f$env$ID[-1]),
                     "the region weights")
})

test_that("without `units`, row order is assumed and a warning says so once", {
        f   <- .a9_fixture()
        old <- .pulse_state$local_row_order
        rm(list = "local_row_order", envir = .pulse_state)
        expect_warning(.a9_fit(f$env, f$W), "row position")
        expect_silent(suppressMessages(
                get_local_types(f$env, f$W, n_local_types = 2, model_name = "VVV")))
        expect_error(.a9_fit(f$env[-1, ], f$W), "must match")
        assign("local_row_order", old, envir = .pulse_state)
})
