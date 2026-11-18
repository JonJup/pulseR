## Round-5 C2 (degenerate CH_FIT) and C7 (gap rule on non-monotone trees).

## ---- C2 ---------------------------------------------------------------
test_that("C2: .pseudo_F_env() gives NA, not Inf or NaN, for degenerate partitions", {
        set.seed(1)
        X <- matrix(stats::rnorm(40), 20, 2)
        one <- .pseudo_F_env(X, rep(1L, 20), k = 3)       # one non-empty region
        expect_true(is.na(one[["CH_FIT"]]))
        expect_equal(attr(one, "k_eff"), 1L)
        sing <- .pseudo_F_env(X, seq_len(20), k = 20)     # one unit per region
        expect_true(is.na(sing[["CH_FIT"]]))
        cst <- .pseudo_F_env(matrix(1, 20, 2), rep(1:2, 10), k = 2)
        expect_true(is.na(cst[["R2_FIT"]]))
        ok <- .pseudo_F_env(X, rep(1:2, each = 10), k = 2)
        expect_true(is.finite(ok[["CH_FIT"]]) && is.finite(ok[["R2_FIT"]]))
})

test_that("C2: an NA CH_FIT ranks last in Borda, not first", {
        df <- data.frame(CH_FIT = c(NA, 5, 3), STAB = c(1, 0.8, 0.7))
        s  <- .borda_score(df, c("max", "max"))
        expect_false(which.min(s) == 1L)
})

c2_toy <- function(n = 8, seed = 1) {
        set.seed(seed)
        g  <- igraph::make_lattice(c(n, n))
        xy <- expand.grid(x = seq_len(n), y = seq_len(n))
        env <- cbind(e1 = (xy$x > n / 2) + stats::rnorm(n^2, 0, 0.2),
                     e2 = stats::rnorm(n^2))
        ep <- igraph::ends(g, igraph::E(g), names = FALSE)
        igraph::E(g)$eid    <- seq_len(igraph::ecount(g))
        igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], ] - env[ep[, 2], ])^2)) + 1e-6
        igraph::graph_attr(g, "env") <- env
        list(g = g, env = env)
}

test_that("C2: tune_regions() logs a one-region consensus as failed and never picks it", {
        tg   <- c2_toy()
        orig <- .salso_k
        local_mocked_bindings(.salso_k = function(memb, k, ...) {
                r <- orig(memb, k, ...)
                if (k == 3L) r$labels[] <- 1L                 # collapse k = 3 to one region
                r
        })
        trees <- sample_spanning_trees(tg$g, n = 20, seed = 1, verbose = FALSE)
        msgs <- character(0)
        withCallingHandlers(
                tr <- tune_regions(tg$g, trees, N_tree = 20, regions_per_partition = 8,
                                   final_regions = 2:3, env = tg$env, stability_B = 0,
                                   seed = 1, cut_rule = "propensity",
                                   refine_boundary = FALSE, min_size = 1, verbose = FALSE),
                warning = function(w) { msgs <<- c(msgs, conditionMessage(w))
                                        invokeRestart("muffleWarning") })
        expect_true(any(grepl("1 of 2 configurations failed", msgs)))
        log <- tr$tuning_log
        bad <- log$final_regions == 3L
        expect_match(log$error[bad], "1 non-empty region")
        expect_true(is.na(log$score[bad]))
        expect_false(any(is.infinite(log$CH_FIT)))
        expect_equal(tr$best_config$final_regions, 2L)
})

test_that("C2: compute_validity() warns and returns NA CH_FIT for one region", {
        tg  <- c2_toy()
        res <- get_regions(tg$g, regions_per_partition = 8, final_regions = 2,
                           N_tree = 20, seed = 1, cut_rule = "propensity",
                           refine_boundary = FALSE, verbose = FALSE)
        res$fit$hard_labels[] <- 1L
        expect_warning(v <- compute_validity(res, env = tg$env, stability_B = 0),
                       "fewer than two non-empty regions")
        expect_true(is.na(v$CH_FIT))
})

## ---- C7 ---------------------------------------------------------------
c7_local <- function() {
        ## Five points in 2-D whose centroid-linkage tree has an inversion.
        P <- rbind(c(0, 0), c(2, 0), c(1, 1.7), c(10, 0), c(10, 3))
        d <- stats::dist(P)
        list(d = d, local = structure(
                list(memberships = diag(5),
                     key = data.frame(name = paste0("t", 1:5), region = c(1, 1, 2, 2, 1))),
                class = "pulse_local_types"))
}

test_that("C7: .n_inversions() detects centroid inversions, not monotone trees", {
        x <- c7_local()
        expect_gt(.n_inversions(stats::hclust(x$d^2, "centroid")), 0L)
        for (m in c("average", "complete", "single", "ward.D2", "mcquitty"))
                expect_identical(.n_inversions(stats::hclust(x$d, m)), 0L)
})

test_that("C7: the gap rule stops on a tree with inversions", {
        x <- c7_local()
        expect_error(get_meta_types(x$local, n_meta_types = 2:4, distance = x$d^2,
                                    linkage = "centroid", rule = "gap"),
                     "gap rule needs a monotone tree")
})

test_that("C7: other rules warn and set gap_ratio to NA", {
        x <- c7_local()
        expect_warning(m <- get_meta_types(x$local, n_meta_types = 2:4, distance = x$d^2,
                                           linkage = "centroid", rule = "crossing"),
                       "inversion")
        expect_true(is.na(m$gap_ratio))
        expect_warning(get_meta_types(x$local, n_meta_types = 2, distance = x$d^2,
                                      linkage = "centroid"), "inversion")
})

test_that("C7: a monotone linkage is unchanged and silent", {
        x <- c7_local()
        m <- get_meta_types(x$local, n_meta_types = 2:3, distance = x$d, linkage = "average")
        expect_true(is.finite(m$gap_ratio))
        expect_true(all(m$selection$height_gap >= 0))
})

test_that("C2: a one-region null replicate counts as failed, not as a score", {
        tg    <- c2_toy()
        orig  <- .salso_k
        calls <- 0L
        local_mocked_bindings(.salso_k = function(memb, k, ...) {
                calls <<- calls + 1L
                r <- orig(memb, k, ...)
                if (calls == 1L) r$labels[] <- 1L
                r
        })
        msgs <- character(0)
        withCallingHandlers(
                nb <- compute_null_baseline(tg$g, tg$env, regions_per_partition = 8,
                                            N_tree = 20, final_regions = 2, B = 5,
                                            stability_B = 0, seed = 1, verbose = FALSE,
                                            cut_rule = "propensity", refine_boundary = FALSE),
                warning = function(w) { msgs <<- c(msgs, conditionMessage(w))
                                        invokeRestart("muffleWarning") })
        expect_true(any(grepl("k = 2 (1 of 5)", msgs, fixed = TRUE)))
        expect_equal(sum(nb$null$failed), 1L)
        expect_match(nb$null$error[nb$null$failed], "1 non-empty region")
        expect_false(any(is.infinite(nb$null$CH_FIT)))
})
