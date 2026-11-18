## A10: get_local_types() accepts exact-prior fits with uncovered units.

.a10_fixture <- function() {
        .local_seed(10)
        g <- igraph::make_lattice(c(10, 10))
        igraph::V(g)$polygon_id <- seq_len(100)
        lab <- rep(c("a", "b"), each = 50)
        lab[c(5, 60)] <- NA                                  # two uncovered units
        r <- get_regions(g, final_regions = 2, N_tree = 5, prior_typology = lab,
                         prior_strength = Inf, verbose = FALSE)
        env <- data.frame(ID = 1:100,
                          x  = stats::rnorm(100) + rep(c(0, 3), each = 50),
                          y  = stats::rnorm(100))
        list(r = r, env = env, uncovered = c(5L, 60L))
}

.a10_fit <- function(env, regions)
        suppressMessages(get_local_types(env, regions, units = env$ID,
                                         n_local_types = 2, model_name = "VVV"))

test_that("uncovered units of a prior_strength = Inf fit become NA orphans", {
        f <- .a10_fixture()
        expect_warning(lt <- .a10_fit(f$env, f$r), "2 of them have no Level 1 region")
        u <- f$uncovered
        expect_true(all(is.na(lt$hard_labels[u])))
        expect_true(all(is.na(lt$region_labels[u])))
        expect_true(all(is.na(lt$certainty[u])))
        expect_equal(unname(rowSums(lt$memberships[u, ])), c(0, 0))
        ## Covered units keep the nesting contract.
        ok <- setdiff(seq_len(100), u)
        expect_false(anyNA(lt$hard_labels[ok]))
        expect_identical(lt$key$region[lt$hard_labels[ok]], lt$region_labels[ok])
        expect_identical(lt$region_labels[ok], f$r$hard_labels[ok])
        expect_length(lt$certainty, 100L)
        expect_output(print(lt), "2 unit\\(s\\) without a Level 1 region")
})

test_that("the tuning wrapper of an exact-prior fit works the same way", {
        f  <- .a10_fixture()
        tr <- structure(list(best_result = list(memberships = f$r$memberships,
                                                hard_labels = f$r$hard_labels)),
                        class = "pulse_regions_tuning")
        lt <- suppressWarnings(.a10_fit(f$env, tr))
        expect_true(all(is.na(lt$hard_labels[f$uncovered])))
})

test_that("an NA region label on a unit with region weight is still an error", {
        f <- .a10_fixture()
        bad <- f$r
        bad$hard_labels[1] <- NA_integer_                    # unit 1 is covered
        expect_error(suppressWarnings(.a10_fit(f$env, bad)),
                     "NA hard label in `regions` but region weight")
})
