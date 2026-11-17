# Planted 4-zone landscape on a lattice; small so it runs in a few seconds.
.dmm_fixture <- function(n = 20) {
        set.seed(1)
        g  <- igraph::make_lattice(c(n, n)); N <- igraph::vcount(g)
        xy <- as.matrix(expand.grid(1:n, 1:n))
        zone <- 1 + (xy[, 1] > n / 2) + 2 * (xy[, 2] > n / 2)
        mix <- rbind(c(.7, .2, .1, 0, 0), c(0, .1, .2, .7, 0),
                     c(.1, 0, 0, .2, .7), c(.3, .3, .4, 0, 0))
        lab <- vapply(zone, function(z) sample(5, 1, prob = mix[z, ]), 1L)
        P <- matrix(0.02, N, 5); P[cbind(seq_len(N), lab)] <- 0.92
        list(g = g, P = P, zone = zone,
             sig = compute_all_signatures(P, g, k_order = 2, verbose = FALSE))
}

test_that("compute_all_signatures() attaches totals", {
        fx <- .dmm_fixture()
        tot <- attr(fx$sig, "totals")
        expect_length(tot, nrow(fx$sig))
        expect_true(all(tot[stats::complete.cases(fx$sig)] > 0))
})

test_that("method = 'dmm' returns a valid pulse_mosaic_types object", {
        fx <- .dmm_fixture()
        m <- suppressWarnings(get_mosaic_types(fx$sig, 2:4, method = "dmm", seed = 1,
                                               verbose = FALSE,
                                               dmm_control = list(n_start = 2, n_scale = 1)))
        expect_s3_class(m, "pulse_mosaic_types")
        expect_s3_class(m$clusters, "pulse_dmm")
        expect_equal(m$method, "dmm")
        expect_equal(rowSums(m$memberships), rep(1, sum(m$valid)), tolerance = 1e-8)
        expect_true(all(c("ICL", "BIC", "n_empty") %in% names(m$validity)))
        expect_equal(dim(m$medoid_signatures), c(m$n_mosaic_types, ncol(fx$sig)))
        expect_gt(.adjusted_rand_index(m$hard_labels, fx$zone[m$valid]), 0.4)
})

test_that("sampled DMM path scores every valid row", {
        fx <- .dmm_fixture()
        m <- suppressWarnings(get_mosaic_types(fx$sig, 3, method = "dmm", sample_size = 150,
                                               seed = 1, verbose = FALSE,
                                               dmm_control = list(n_start = 1)))
        expect_true(m$sampled)
        expect_equal(nrow(m$memberships), sum(m$valid))
})

test_that("engine/metric mismatches and missing totals fail early", {
        fx <- .dmm_fixture()
        expect_error(get_mosaic_types(fx$sig, 2:3, method = "dmm", metric = "GD5",
                                      verbose = FALSE), "needs PAM medoids")
        expect_error(get_mosaic_types(fx$sig, 2:3, metric = "ICL", verbose = FALSE, method = "pam"),
                     "needs a likelihood")
        expect_error(get_mosaic_types(unclass(fx$sig)[, ], 2:3, method = "dmm",
                                      verbose = FALSE), "totals")
        expect_error(get_mosaic_types(fx$sig, 2:3, method = "dmm", verbose = FALSE,
                                      dmm_control = list(nstart = 2)), "Unknown")
})

test_that("check_k_order() runs for both engines (compute_mosaic_contrast stack overflow)", {
        fx <- .dmm_fixture(n = 12)
        expect_length(compute_mosaic_contrast(fx$sig, fx$g), nrow(fx$sig))
        sp <- suppressMessages(check_k_order(fx$P, fx$g, k_order = 1:2, n_mosaic_types = 2:3,
                                             stability_B = 0, verbose = FALSE, method = "pam"))
        sd <- suppressWarnings(suppressMessages(
                check_k_order(fx$P, fx$g, k_order = 1:2, n_mosaic_types = 2:3, method = "dmm",
                              verbose = FALSE, dmm_control = list(n_start = 1))))
        expect_equal(nrow(attr(sp, "by_order")), 2L)
        expect_equal(nrow(attr(sd, "by_order")), 2L)
})

test_that(".fit_dmm() recovers finite concentrations", {
        set.seed(4); K <- 3; C <- 20; N <- 900; A <- c(20, 50, 200)
        m <- matrix(stats::rgamma(K * C, 0.7), K); m <- m / rowSums(m)
        z <- sample(K, N, TRUE)
        X <- t(vapply(seq_len(N), function(i) {
                p <- stats::rgamma(C, A[z[i]] * m[z[i], ]); p <- p / sum(p)
                as.numeric(stats::rmultinom(1, 150, p))
        }, numeric(C)))
        f <- .fit_dmm(X, K, n_start = 2, seed = 1)
        expect_equal(f$limit, "dirichlet-multinomial")
        expect_true(f$converged)
        expect_lt(length(f$trace), 30)
        expect_gt(.adjusted_rand_index(f$cluster, z), 0.95)
        conc <- f$concentration[order(f$concentration)]
        expect_true(all(abs(log(conc / A)) < log(1.6)))      # within a factor of 1.6
})

test_that(".fit_dmm() returns the multinomial limit when there is no overdispersion", {
        set.seed(5); K <- 2; C <- 10; N <- 400
        th <- matrix(stats::rgamma(K * C, 1), K); th <- th / rowSums(th)
        z <- sample(K, N, TRUE)
        X <- t(vapply(z, function(k) as.numeric(stats::rmultinom(1, 40, th[k, ])), numeric(C)))
        f <- .fit_dmm(X, K, n_start = 1, seed = 1)
        expect_equal(f$limit, "multinomial")
        expect_equal(f$df, K * (C - 1) + (K - 1))
        expect_gt(.adjusted_rand_index(f$cluster, z), 0.95)
})

test_that("DMM with metric = 'ICL' skips the N x N distance matrix", {
        fx <- .dmm_fixture()
        m <- suppressWarnings(get_mosaic_types(fx$sig, 2:3, method = "dmm", seed = 1,
                                               verbose = FALSE, dmm_control = list(n_start = 1)))
        expect_true(all(is.na(m$validity$SILH_fuzzy)))
        expect_true(is.finite(m$saturation))
})
test_that("n_scale = 'auto' spends each edge once, whatever k_order", {
        fx <- .dmm_fixture()
        for (k in 1:3) {
                s <- compute_all_signatures(fx$P, fx$g, k_order = k, verbose = FALSE)
                X <- .sig_to_counts(s, n_scale = .resolve_n_scale("auto", s))
                expect_equal(sum(X), igraph::ecount(fx$g), tolerance = 1e-8)
        }
})
test_that(".fit_dmm() gives a clear error when K exceeds distinct compositions", {
        X <- rbind(c(1, 2), c(2, 4), c(3, 6), c(5, 1))   # 2 distinct compositions
        expect_error(.fit_dmm(X, 3L, seed = 1), "distinct row compositions")
        expect_silent(.fit_dmm(X, 2L, n_start = 1, max_iter = 20L, seed = 1))
})
## --- temperature ----------------------------------------------------------
.dmm_sim_fixture <- function(seed = 4, N = 600, n = 150) {
        set.seed(seed); K <- 3; C <- 20; A <- c(20, 50, 200)
        m <- matrix(stats::rgamma(K * C, 0.7), K); m <- m / rowSums(m)
        z <- sample(K, N, TRUE)
        X <- t(vapply(seq_len(N), function(i) {
                p <- stats::rgamma(C, A[z[i]] * m[z[i], ]); p <- p / sum(p)
                as.numeric(stats::rmultinom(1, n, p))
        }, numeric(C)))
        list(X = X, z = z, K = K)
}

test_that("temperature 0 returns pi; certainty falls monotonically with tau", {
        s <- .dmm_sim_fixture()
        f <- .fit_dmm(s$X, s$K, n_start = 3, seed = 1)
        expect_equal(f$tau, 1)
        U0 <- .score_dmm(f, s$X, tau = 0)$U
        expect_equal(U0, matrix(f$pi, nrow(s$X), s$K, byrow = TRUE))
        expect_equal(.score_dmm(f, s$X, tau = 0)$logpred, .score_dmm(f, s$X)$logpred)
        cert <- vapply(c(1, 0.5, 0.1, 0.02), function(t)
                mean(apply(.score_dmm(f, s$X, tau = t)$U, 1, max)), numeric(1))
        expect_true(all(diff(cert) < 0))
})

test_that("tau acts after fitting: same alpha, softer U, no collapse", {
        s <- .dmm_sim_fixture()
        f1 <- .fit_dmm(s$X, s$K, n_start = 3, seed = 1)
        ft <- .fit_dmm(s$X, s$K, n_start = 3, tau = 0.1, seed = 1)
        expect_equal(ft$tau, 0.1)
        expect_equal(ft$alpha, f1$alpha)
        expect_equal(ft$BIC - f1$BIC, -2 * (0.1 - 1) * f1$loglik)
        expect_gt(.adjusted_rand_index(ft$cluster, s$z), 0.9)
        expect_lt(mean(apply(ft$U, 1, max)), mean(apply(f1$U, 1, max)))
        expect_equal(.score_dmm(ft, s$X)$U, ft$U, tolerance = 1e-10)  # scoring reuses fit$tau
        expect_equal(.temper_dmm(f1, s$X, 0.1)$U, ft$U, tolerance = 1e-10)
})

test_that("calibration on a correctly specified DMM keeps tau near 1", {
        s <- .dmm_sim_fixture(N = 400, n = 30)
        f <- .fit_dmm(s$X, s$K, n_start = 3, seed = 1)
        cal <- .calibrate_dmm_tau(f, rowSums(s$X), B = 10, seed = 1)
        expect_true(cal$tau >= 0.5)
        expect_equal(nrow(cal$table), 25L)
        expect_lte(min(cal$table$ece), cal$table$ece[cal$table$tau == 1])
})

test_that("get_mosaic_types(): temperature 'auto' = 1/multiplicity, 'calibrate' records provenance", {
        fx <- .dmm_fixture()
        a <- suppressWarnings(get_mosaic_types(fx$sig, 4, method = "dmm", seed = 1, verbose = FALSE,
                                               dmm_control = list(n_start = 1)))
        expect_equal(a$dmm_control$temperature, 1 / attr(fx$sig, "edge_multiplicity"))
        expect_equal(a$clusters$tau, a$dmm_control$temperature)
        expect_null(a$calibration)
        expect_gt(.adjusted_rand_index(a$hard_labels, fx$zone[a$valid]), 0.4)
        cmb <- suppressWarnings(get_mosaic_types(fx$sig, 4, method = "dmm", seed = 1, verbose = FALSE,
                                                 dmm_control = list(n_start = 1, temperature = "calibrate",
                                                                    calibration_B = 5)))
        expect_equal(cmb$calibration$tau_select, a$dmm_control$temperature)
        expect_equal(cmb$clusters$tau, cmb$calibration$tau)
        expect_equal(rowSums(cmb$memberships), rep(1, sum(cmb$valid)), tolerance = 1e-8)
        expect_warning(get_mosaic_types(fx$sig, 4, method = "dmm", seed = 1, verbose = FALSE,
                                        dmm_control = list(n_start = 1, n_scale = 0.5)),
                       "corrected twice")
        expect_error(get_mosaic_types(fx$sig, 4, method = "dmm", verbose = FALSE,
                                      dmm_control = list(temperature = 2)), "\\(0, 1\\]")
})
