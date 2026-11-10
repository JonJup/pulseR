# Regression tests for the round-3 review fixes.

two_region_env <- function(n = 300, seed = 2, p = 1) {
        .local_seed(seed)
        env <- data.frame(a = c(stats::rnorm(n / 2), stats::rnorm(n / 2, 4)))
        if (p > 1) env$b <- stats::rnorm(n)
        w1 <- stats::plogis(seq(-5, 5, length.out = n))
        list(env = env, W = cbind(w1, 1 - w1))
}

test_that("A4: get_meta_types() works with one predictor", {
        s   <- two_region_env()
        loc <- get_local_types(s$env, s$W, n_local_types = 2:3)
        m   <- get_meta_types(loc, n_meta_types = 2)
        expect_s3_class(m, "pulse_meta_types")
        expect_equal(unname(rowSums(m$memberships)), rep(1, nrow(s$env)))
        expect_equal(.mclust_sigma(list(sigmasq = c(2, 3)), 2, 1), matrix(3))
        expect_equal(.mclust_sigma(list(sigmasq = 2), 2, 1), matrix(2))
})

test_that("E1: hc runs on at most hc_max_units core units", {
        s <- two_region_env(p = 2)
        f <- get_local_types(s$env, s$W, n_local_types = 2, model_name = "VVV",
                             hc_max_units = 50)
        expect_equal(unname(f$n_local_types), c(2L, 2L))
        expect_equal(unname(rowSums(f$memberships)), rep(1, nrow(s$env)))
        expect_error(get_local_types(s$env, s$W, hc_max_units = 1), "hc_max_units")
})

test_that("B5: local certainty and margin are within the assigned region", {
        U   <- rbind(c(0.45, 0.15, 0.4), c(0.3, 0.2, 0.5))
        loc <- structure(list(memberships = U, weighted = U, hard_labels = c(1L, 3L),
                              certainty = apply(U, 1, max),
                              key = data.frame(column = 1:3, region = c(1L, 1L, 2L),
                                               local_type = c(1L, 2L, 1L))),
                         class = "pulse_local_types")
        l <- .level_uncertainty(loc, "local")
        expect_equal(l$certainty, c(0.75, 1))
        expect_equal(l$margin, c(0.5, 1))
        expect_equal(l$k, c(2L, 1L))
})

test_that("B10: ties are broken by one explicit rule", {
        log <- data.frame(regions_per_partition = c(20, 10, 10, 10),
                          weight_exponent = c(1, 2, 1, 1), N_tree = c(50, 50, 50, 100),
                          final_regions = c(3, 3, 3, 3), score = c(1, 1, 1, 0.5))
        expect_identical(.pick_best(log), 3L)
        log$final_regions[1] <- 2
        expect_identical(.pick_best(log), 1L)
        log$score[] <- NA
        expect_identical(.pick_best(log), NA_integer_)
})

test_that("B10: tune_regions() reports ties and the runner-up", {
        g <- igraph::make_lattice(c(6, 6))
        igraph::E(g)$weight <- stats::runif(igraph::ecount(g))
        tr  <- sample_spanning_trees(g, n = 8, seed = 1, verbose = FALSE)
        res <- tune_regions(g, tr, N_tree = 8, regions_per_partition = 6,
                            final_regions = 2:3, stability_B = 0,
                            tuning_metric = "expected_vi", seed = 1, verbose = FALSE)
        expect_gte(res$n_tied, 1L)
        expect_identical(res$tie_rule, .TIE_RULE)
        expect_identical(nrow(res$runner_up), 1L)
        expect_false(res$runner_up$final_regions == res$best_final_regions)
})

test_that("A8 + C3: an exact prior gives a tuning object the crosswalk accepts", {
        g <- igraph::make_lattice(c(6, 6))
        igraph::E(g)$weight <- stats::runif(igraph::ecount(g))
        tr  <- sample_spanning_trees(g, n = 4, seed = 1, verbose = FALSE)
        pri <- rep(c("a", "b", NA), each = 12)
        res <- tune_regions(g, tr, N_tree = 4, final_regions = 2, prior_strength = Inf,
                            prior_typology = pri, verbose = FALSE)
        expect_s3_class(res, "pulse_regions_tuning")
        expect_identical(res$best_result$hard_labels, rep(c(1L, 2L, NA), each = 12))
        ext <- data.frame(ID = 1:36, lab = rep(c("x", "y"), 18))
        cw  <- crosswalk_types(res, ext, label = "lab", units = 1:36, min_coverage = 0)
        expect_identical(cw$coverage$n_typed, 24L)
        expect_identical(cw$coverage$n_no_pulse_type, 12L)
})

test_that("C4: a zero MAD falls back to another scale and trees still vary", {
        expect_equal(.noise_scale(c(1, 2, 3, 4)), stats::mad(c(1, 2, 3, 4)))
        w <- c(rep(1, 10), 5)
        expect_warning(s <- .noise_scale(w), "MAD")
        expect_equal(s, stats::sd(w))
        expect_equal(suppressWarnings(.noise_scale(rep(2, 5))), 1)
        expect_silent(.noise_scale(w, warn = FALSE))
        g <- igraph::make_lattice(c(6, 6))
        igraph::E(g)$weight <- c(rep(1, igraph::ecount(g) - 3L), 2, 3, 4)
        expect_warning(tr <- sample_spanning_trees(g, n = 6, seed = 1, verbose = FALSE), "MAD")
        ids <- lapply(seq_len(tr$n_trees), function(j) sort(tr$edge_ids[, j]))
        expect_gt(length(unique(ids)), 1L)
        expect_gt(tr$noise_scale, 0)
})

test_that("C5: bare ensemble matrices are validated and 0-based labels relabelled", {
        M <- matrix(c(0, 0, 1, 1, 0, 1, 1, 1), 4)
        expect_identical(.ensemble_meta(M)$memb_mat,
                         matrix(c(1L, 1L, 2L, 2L, 1L, 2L, 2L, 2L), 4))
        expect_error(.ensemble_meta(matrix(c(1, NA, 2, 2), 2)), "whole-number")
        expect_error(.ensemble_meta(matrix(c(1.5, 1, 2, 2), 2)), "whole-number")
})

test_that("C7 + D3: duplicated N_tree does not tune; one ensemble constructor", {
        g <- igraph::make_lattice(c(6, 6))
        igraph::E(g)$weight <- stats::runif(igraph::ecount(g))
        res <- get_regions(g, regions_per_partition = 6, final_regions = 2,
                           N_tree = c(8, 8), seed = 1, verbose = FALSE)
        expect_false(res$tuned)
        tr <- sample_spanning_trees(g, n = 8, seed = 1, verbose = FALSE)
        e  <- compute_ensemble(tr, 8, g, 6, verbose = FALSE)
        ref <- .new_ensemble(e$memberships, 8, 6, 1, 0, tr$n_vertices,
                             cut_rule = "propensity", min_size = 1L)
        expect_identical(e, ref)
        expect_type(e$n_vertices, "integer")
        expect_type(.new_ensemble(e$memberships, 8, 6, 1, 0, 36)$n_vertices, "integer")
})

test_that("E2: the JS dist is built directly and the dense matrix is reusable", {
        .local_seed(1)
        S  <- matrix(stats::runif(40), 10); S <- S / rowSums(S)
        d  <- compute_signature_distance(S, "js")
        eps <- .Machine$double.eps
        Sn <- (S + eps) / rowSums(S + eps)
        expect_s3_class(d, "dist")
        expect_identical(attr(d, "Size"), 10L)
        D  <- as.matrix(d)
        expect_equal(D[2, 1], .js_pairwise(Sn[1:2, ]))
        expect_equal(D[10, 9], .js_pairwise(Sn[9:10, ]))
        f1 <- .pam_memberships(d, 3)
        f2 <- .pam_memberships(d, 3, Dm = D)
        expect_identical(f1$memberships, f2$memberships)
})

test_that("salso's cluster-cap warning becomes the `capped` flag", {
        local_mocked_bindings(salso = function(x, ...) {
                warning("The number of clusters equals the default maximum possible number of clusters.")
                x[1, ]
        }, .package = "salso")
        memb <- cbind(c(1L, 1L, 2L, 2L), c(1L, 2L, 2L, 2L))
        expect_no_warning(nat <- .salso_native(memb))
        expect_true(nat$capped)
        local_mocked_bindings(salso = function(x, ...) x[1, ], .package = "salso")
        expect_false(.salso_native(memb)$capped)
})
test_that("empirical distance samples n_sample rows, reproducibly, RNG untouched", {
        s   <- two_region_env(n = 300, p = 2)
        loc <- get_local_types(s$env, s$W, n_local_types = 2)
        full <- compute_local_type_distance(loc, "empirical", s$env, n_sample = 1e6)
        expect_equal(full, compute_local_type_distance(loc, "empirical", s$env))  # N < n_sample: no sampling
        set.seed(99); before <- .Random.seed
        a <- compute_local_type_distance(loc, "empirical", s$env, n_sample = 100, seed = 7)
        expect_identical(.Random.seed, before)
        expect_equal(a, compute_local_type_distance(loc, "empirical", s$env, n_sample = 100, seed = 7))
        expect_error(compute_local_type_distance(loc, "empirical", s$env, n_sample = 1), "n_sample")
})