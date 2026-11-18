test_that(".zsum_score orients, logs CH_FIT and handles NA", {
        df <- data.frame(CH_FIT = c(10, 100, 1000), EVI = c(3, 2, 1))
        z  <- .zsum_score(df, c("max", "min"))
        expect_equal(z, c(-1, 0, 1))                       # log10 spacing is even
        expect_lt(.zsum_score(df, c("max", "min"), log_cols = character(0))[2], 0)
        df$EVI[3] <- NA                                    # NA -> lowest z of the column
        z  <- .zsum_score(df, c("max", "min"))
        expect_equal(z[3], (1 - 1 / sqrt(2)) / 2)
        expect_true(all(is.na(.zsum_score(data.frame(a = c(NA, NA)), "max"))))
        expect_equal(.zsum_score(data.frame(a = c(1, 1, 1)), "max"), c(0, 0, 0))
})

test_that("z-sum breaks a Borda tie by the size of the differences", {
        # A wins CH_FIT by a lot, B wins STAB by a hair: Borda ties them.
        df <- data.frame(CH_FIT = c(200, 100, 50), STAB = c(0.80, 0.81, 0.20))
        b  <- .derive_score_column(df, "borda")
        expect_equal(as.numeric(b)[1], as.numeric(b)[2])
        log <- data.frame(regions_per_partition = 30, weight_exponent = 4,
                          N_tree = 50, final_regions = c(4, 3, 2),
                          score = as.numeric(b), score_tb = attr(b, "tiebreak"))
        expect_identical(.pick_best(log), 1L)              # old rule: row 2 (k = 3)
        expect_identical(.tie_resolution(log$score, log$score_tb), "zsum")
        log$score_tb <- NULL
        expect_identical(.pick_best(log), 2L)
})

test_that("z-sum never overrides a Borda rank difference", {
        df <- data.frame(CH_FIT = c(1000, 10, 9), STAB = c(0.10, 0.90, 0.89))
        b  <- .derive_score_column(df, "borda")
        log <- data.frame(regions_per_partition = 30, weight_exponent = 4,
                          N_tree = 50, final_regions = 2:4,
                          score = as.numeric(b), score_tb = attr(b, "tiebreak"))
        expect_identical(.pick_best(log), which.max(log$score))
        expect_identical(.tie_resolution(log$score, log$score_tb), "none")
})

test_that("two rows that split the voters fall through to the axis rule", {
        df <- data.frame(CH_FIT = c(200, 100), STAB = c(0.5, 0.9))
        b  <- .derive_score_column(df, "borda")
        log <- data.frame(regions_per_partition = 30, weight_exponent = 4,
                          N_tree = 50, final_regions = c(5, 3),
                          score = as.numeric(b), score_tb = attr(b, "tiebreak"))
        expect_identical(.pick_best(log), 2L)
        expect_identical(.tie_resolution(log$score, log$score_tb), "rule")
})

test_that("non-Borda metrics carry no tie-break", {
        df <- data.frame(CH_FIT = 1:3, STAB = 3:1)
        expect_null(attr(.derive_score_column(df, "ch_fit"), "tiebreak"))
})
