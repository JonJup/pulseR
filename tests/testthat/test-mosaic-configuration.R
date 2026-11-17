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


# Planted arrangement-only landscape: left half patchy (square blocks of one
# type), right half salt-and-pepper, three types in equal shares in both halves.
.arrangement_fixture <- function(n = 30, block = 5, seed = 1) {
        set.seed(seed)
        g  <- igraph::make_lattice(c(n, n)); N <- igraph::vcount(g)
        xy <- as.matrix(expand.grid(x = 1:n, y = 1:n))
        zone <- 1L + (xy[, 1] > n / 2)
        bx <- (xy[, 1] - 1) %/% block; by <- (xy[, 2] - 1) %/% block
        blocks <- unique(cbind(bx, by))
        btype <- sample(rep_len(1:3, nrow(blocks)))
        patchy <- btype[match(paste(bx, by), paste(blocks[, 1], blocks[, 2]))]
        mixed  <- sample(rep_len(1:3, N))
        lab <- ifelse(zone == 1L, patchy, mixed)
        P <- matrix(0, N, 3); P[cbind(seq_len(N), lab)] <- 1
        list(g = g, P = P, lab = lab, zone = zone)
}

test_that("compute_configuration() is invariant to relabelling the types", {
        fx <- .arrangement_fixture(n = 16, block = 4)
        a <- compute_configuration(fx$P, fx$g, k_order = 3, verbose = FALSE)
        b <- compute_configuration(fx$P[, c(3, 1, 2)], fx$g, k_order = 3, verbose = FALSE)
        expect_s3_class(a, "pulse_configuration")
        expect_equal(dim(a), c(nrow(fx$P), 3L))
        expect_equal(unclass(a)[, ], unclass(b)[, ], ignore_attr = TRUE)
})

test_that("kappa is high in patches and near zero in random mixing", {
        fx <- .arrangement_fixture()
        cf <- compute_configuration(fx$P, fx$g, k_order = 3, verbose = FALSE)
        k1 <- tapply(cf[, "kappa_lag1"], fx$zone, mean, na.rm = TRUE)
        expect_gt(k1[["1"]], 0.5)  # patch interiors count as clumped
        expect_lt(abs(k1[["2"]]), 0.15)
})

test_that("kappa is 0 for a single-type window and does not blow up", {
        fx <- .arrangement_fixture(n = 12, block = 12)        # one block = one type
        cf <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, verbose = FALSE)
        expect_true(all(is.finite(cf[fx$zone == 1L, ])))
        expect_true(all(abs(cf[fx$zone == 1L, ]) <= 1 / attr(cf, "kappa_floor")))
        expect_error(compute_configuration(fx$P, fx$g, 2, k_reference = 1, verbose = FALSE),
                     "k_reference")
})

test_that("configuration channel separates arrangement-only zones; signatures alone do not", {
        ## Patches (8 x 8) wider than the signature window (k = 2): the
        ## label-dependent regime in which signatures alone fail.
        fx  <- .arrangement_fixture(n = 32, block = 8)
        sig <- compute_all_signatures(fx$P, fx$g, k_order = 2, verbose = FALSE)
        cf  <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, k_reference = 8,
                                     verbose = FALSE)
        m0 <- suppressWarnings(get_mosaic_types(sig, 2, method = "dmm", seed = 1, verbose = FALSE,
                                                dmm_control = list(n_start = 3)))
        m1 <- suppressWarnings(get_mosaic_types(sig, 2, method = "dmm", seed = 1, verbose = FALSE,
                                                dmm_control = list(n_start = 3),
                                                configuration = cf))
        ari0 <- .adjusted_rand_index(m0$hard_labels, fx$zone[m0$valid])
        ari1 <- .adjusted_rand_index(m1$hard_labels, fx$zone[m1$valid])
        expect_gt(ari1, 0.8)
        expect_gt(ari1, ari0 + 0.3)
        expect_equal(rowSums(m1$memberships), rep(1, sum(m1$valid)), tolerance = 1e-8)
        expect_equal(m1$configuration$features, colnames(cf))
        expect_equal(dim(m1$configuration$means), c(2L, ncol(cf)))
})

test_that("free choice of k finds the two arrangement zones", {
        fx  <- .arrangement_fixture(n = 32, block = 8)
        sig <- compute_all_signatures(fx$P, fx$g, k_order = 2, verbose = FALSE)
        cf  <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, k_reference = 8,
                                     verbose = FALSE)
        m <- suppressWarnings(get_mosaic_types(sig, 2:5, method = "dmm", seed = 1, verbose = FALSE,
                                               dmm_control = list(n_start = 3),
                                               configuration = cf))
        expect_equal(m$n_mosaic_types, 2L)
        expect_gt(.adjusted_rand_index(m$hard_labels, fx$zone[m$valid]), 0.8)
})

test_that("composition zones are still recovered with the configuration channel on", {
        fx <- .dmm_fixture()
        cf <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, verbose = FALSE)
        m  <- suppressWarnings(get_mosaic_types(fx$sig, 2:5, method = "dmm", seed = 1,
                                                verbose = FALSE,
                                                dmm_control = list(n_start = 2),
                                                configuration = cf))
        expect_gt(.adjusted_rand_index(m$hard_labels, fx$zone[m$valid]), 0.4)
})

test_that("sampled path scores every valid row with configuration", {
        fx  <- .arrangement_fixture()
        sig <- compute_all_signatures(fx$P, fx$g, k_order = 3, verbose = FALSE)
        cf  <- compute_configuration(fx$P, fx$g, k_order = 3, verbose = FALSE)
        m <- suppressWarnings(get_mosaic_types(sig, 2, method = "dmm", sample_size = 300,
                                               seed = 1, verbose = FALSE,
                                               dmm_control = list(n_start = 1),
                                               configuration = cf))
        expect_true(m$sampled)
        expect_equal(nrow(m$memberships), sum(m$valid))
})

test_that("configuration input is validated", {
        fx  <- .arrangement_fixture(n = 12)
        sig <- compute_all_signatures(fx$P, fx$g, k_order = 2, verbose = FALSE)
        cf  <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, verbose = FALSE)
        expect_error(get_mosaic_types(sig, 2, method = "pam", verbose = FALSE,
                                      configuration = cf), "needs method = 'dmm'")
        expect_error(get_mosaic_types(sig, 2, verbose = FALSE, configuration = cf[-1, ]),
                     "one row per")
        expect_error(get_mosaic_types(sig, 2, verbose = FALSE, configuration = cf,
                                      dmm_control = list(temperature = "calibrate")),
                     "calibrate")
        expect_error(compute_configuration(fx$P, fx$g, 2, lags = c(0, 1), verbose = FALSE),
                     "lags")
})

test_that("smooth_k = 0 changes nothing; smoothing keeps invariance and cuts within-zone noise", {
        fx <- .arrangement_fixture(n = 32, block = 8)
        a  <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, k_reference = 8,
                                    verbose = FALSE)
        a0 <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, k_reference = 8,
                                    smooth_k = 0, verbose = FALSE)
        expect_equal(unclass(a)[, ], unclass(a0)[, ], ignore_attr = TRUE)
        s1 <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, k_reference = 8,
                                    smooth_k = 3, verbose = FALSE)
        s2 <- compute_configuration(fx$P[, c(2, 3, 1)], fx$g, k_order = 2, lags = 1:2,
                                    k_reference = 8, smooth_k = 3, verbose = FALSE)
        expect_equal(unclass(s1)[, ], unclass(s2)[, ], ignore_attr = TRUE)
        expect_equal(attr(s1, "smooth_k"), 3L)
        wsd <- function(x) mean(tapply(x[, 1], fx$zone, sd, na.rm = TRUE))
        expect_lt(wsd(s1), wsd(a))
        k1 <- tapply(s1[, 1], fx$zone, mean, na.rm = TRUE)
        expect_gt(k1[["1"]] - k1[["2"]], 0.4)        # zones still far apart
        expect_error(compute_configuration(fx$P, fx$g, 2, lags = 1:2, smooth_k = -1,
                                           verbose = FALSE), "smooth_k")
})

test_that("smoothed features still separate the arrangement zones", {
        fx  <- .arrangement_fixture(n = 32, block = 8)
        sig <- compute_all_signatures(fx$P, fx$g, k_order = 2, verbose = FALSE)
        cf  <- compute_configuration(fx$P, fx$g, k_order = 2, lags = 1:2, k_reference = 8,
                                     smooth_k = 3, verbose = FALSE)
        m <- suppressWarnings(get_mosaic_types(sig, 2:5, method = "dmm", seed = 1, verbose = FALSE,
                                               dmm_control = list(n_start = 3),
                                               configuration = cf))
        expect_gt(.adjusted_rand_index(m$hard_labels, fx$zone[m$valid]), 0.8)
})
