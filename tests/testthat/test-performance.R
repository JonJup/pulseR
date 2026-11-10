## tests/testthat/test-performance.R
## E1-E4 (round-4 review): faster SSD cuts, faster boundary refinement, cheaper
## null, single ensemble copy. Each speed-up must leave results unchanged.

.perf_graph <- function(side = 14, seed = 1) {
        set.seed(seed)
        g   <- igraph::make_lattice(c(side, side))
        xy  <- expand.grid(x = seq_len(side), y = seq_len(side))
        env <- cbind(e1 = (xy$x > side / 3) + (xy$x > 2 * side / 3) + stats::rnorm(side^2, 0, 0.3),
                     e2 = stats::rnorm(side^2))
        ep  <- igraph::ends(g, igraph::E(g), names = FALSE)
        igraph::E(g)$eid    <- seq_len(igraph::ecount(g))
        igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], ] - env[ep[, 2], ])^2)) + 1e-6
        g   <- igraph::set_graph_attr(g, "env", env)
        list(g = g, env = env, xy = xy)
}

# Gains by brute force: components and side means from scratch, per edge.
.brute_gains <- function(ep, X) {
        n  <- nrow(X)
        t(vapply(seq_len(nrow(ep)), function(e) {
                lab <- .labels_from_edges(ep[-e, , drop = FALSE], n)
                s1  <- lab == lab[ep[e, 1]]; s2 <- lab == lab[ep[e, 2]]
                n1  <- sum(s1); n2 <- sum(s2)
                d   <- colMeans(X[s1, , drop = FALSE]) - colMeans(X[s2, , drop = FALSE])
                h   <- n1 * n2 / (n1 + n2)
                c(h * sum(d^2), h, min(n1, n2))
        }, numeric(3))) |> matrix(ncol = 3)
}

test_that("E1: forest gains match brute force on forests with several components", {
        set.seed(11)
        for (r in 1:20) {
                n  <- sample(5:40, 1)
                ep <- igraph::as_edgelist(igraph::sample_tree(n), names = FALSE)
                ep <- ep[-sample.int(nrow(ep), sample(0:3, 1)), , drop = FALSE]
                X  <- cbind(stats::rnorm(n), stats::rnorm(n)) + 1e3      # offset: precision
                fc <- .forest_cut_gains(ep, X)
                bf <- .brute_gains(ep, X)
                expect_equal(fc$gain,  bf[, 1], tolerance = 1e-8)
                expect_equal(fc$h,     bf[, 2])
                expect_equal(fc$n_min, bf[, 3])
        }
})

test_that("E1: incremental cuts reproduce a from-scratch race draw for draw", {
        x  <- .perf_graph()
        tr <- sample_spanning_trees(x$g, n = 6, seed = 1, verbose = FALSE,
                                    prior_typology = ifelse(x$xy$y <= 7, 1L, 2L))
        X  <- x$env
        ref <- function(k, gam, si, ps, ms) {             # pre-E1 algorithm
                vapply(seq_len(tr$n_trees), function(j) {
                        ep <- tr$edge_endpoints[tr$edge_ids[, j], , drop = FALSE]
                        wp <- if (ps > 0) tr$within_prior_edge[tr$edge_ids[, j]] else NULL
                        for (s in seq_len(k - 1L)) {
                                fc <- .forest_cut_gains(ep, X)
                                lp <- gam * log(pmax(fc$gain, .Machine$double.xmin))
                                u  <- stats::sd(lp); ok <- is.finite(u) && u > 0
                                if (!ok) u <- log(2)
                                if (si > 0 && ok) {
                                        lh <- log(fc$h); sh <- stats::sd(lh)
                                        if (is.finite(sh) && sh > 0) lp <- lp + si * u / sh * lh
                                }
                                if (!is.null(wp) && any(wp)) lp <- lp - ps * u * wp
                                if (ms > 1L) lp[fc$n_min < ms] <- -Inf
                                d  <- which.min(log(stats::rexp(nrow(ep))) - lp)
                                ep <- ep[-d, , drop = FALSE]
                                if (!is.null(wp)) wp <- wp[-d]
                        }
                        .labels_from_edges(ep, nrow(X))
                }, integer(nrow(X)))
        }
        for (a in list(c(8, 1, 2, 0, 1), c(12, 4, 0, 1, 5), c(5, 0, 2, 2, 1))) {
                set.seed(9); r <- ref(a[1], a[2], a[3], a[4], a[5])
                set.seed(9); m <- suppressWarnings(.compute_ensemble_internal(
                        x$g, tr, tr$n_trees, a[1], prior_strength = a[4], weight_exponent = a[2],
                        size_influence = a[3], cut_rule = "ssd", min_size = a[5]))
                expect_identical(m, r)
        }
})

test_that("E2: boundary refinement makes the same moves as the whole-region check", {
        ref <- function(z, graph, X, max_pass = 20L, min_size = 1L) {   # pre-E2 algorithm
                z   <- match(z, sort(unique(z))); k <- max(z)
                nb  <- lapply(igraph::adjacent_vertices(graph, igraph::V(graph)), as.integer)
                n_r <- tabulate(z, k); S <- rowsum(X, z, reorder = TRUE); moves <- 0L
                for (pass in seq_len(max_pass)) {
                        moved <- 0L
                        bnd <- which(vapply(seq_along(z), function(i) any(z[nb[[i]]] != z[i]), logical(1)))
                        for (i in bnd) {
                                a <- z[i]; cand <- setdiff(unique(z[nb[[i]]]), a)
                                if (!length(cand) || n_r[a] <= min_size) next
                                out  <- n_r[a] / (n_r[a] - 1) * sum((X[i, ] - S[a, ] / n_r[a])^2)
                                gain <- vapply(cand, function(b)
                                        out - n_r[b] / (n_r[b] + 1) * sum((X[i, ] - S[b, ] / n_r[b])^2), numeric(1))
                                if (max(gain) <= 0) next
                                b <- cand[which.max(gain)]
                                same <- nb[[i]][z[nb[[i]]] == a]
                                if (length(same) > 1L) {
                                        ia <- setdiff(which(z == a), i)
                                        cm <- igraph::components(igraph::induced_subgraph(graph, ia))$membership
                                        if (length(unique(cm[match(same, ia)])) > 1L) next
                                }
                                z[i] <- b
                                S[a, ] <- S[a, ] - X[i, ]; S[b, ] <- S[b, ] + X[i, ]
                                n_r[a] <- n_r[a] - 1L; n_r[b] <- n_r[b] + 1L
                                moved <- moved + 1L
                        }
                        moves <- moves + moved
                        if (!moved) break
                }
                list(z = z, moves = moves, passes = pass)
        }
        for (sd in 1:3) {
                x   <- .perf_graph(side = 25, seed = sd)
                X   <- cbind(sin(x$xy$x / 4) + stats::rnorm(625, sd = 0.6),
                             cos(x$xy$y / 5) + stats::rnorm(625, sd = 0.6))
                ctr <- x$xy[sample(625, 8), ]
                z0  <- max.col(-sapply(1:8, function(r) (x$xy$x - ctr$x[r])^2 + (x$xy$y - ctr$y[r])^2),
                               ties.method = "first")
                for (ms in c(1L, 10L)) {
                        r <- .refine_regions(z0, x$g, X, min_size = ms)
                        expect_identical(r, ref(z0, x$g, X, min_size = ms))
                        expect_true(r$moves > 0L)
                }
        }
})

test_that("E2: the connectivity check is exact at cut vertices and around holes", {
        g  <- igraph::make_lattice(c(5, 5))
        nb <- lapply(igraph::adjacent_vertices(g, igraph::V(g)), as.integer)
        ok <- .region_connect_checker(nb, 25)
        z  <- rep(1L, 25); z[c(7, 8, 9, 12, 14, 17, 18, 19)] <- 2L       # ring around 13
        z[13] <- 3L
        expect_true(ok(8, c(7, 9), z))           # ring stays linked the long way round
        z2 <- z; z2[c(17, 18, 19)] <- 1L         # open ring: 7-8-9-14 / 12
        expect_false(ok(7, c(8, 12), z2))        # 7 joins 8 and 12: true cut vertex
        expect_true(ok(1, c(2, 6), rep(1L, 25)))
        expect_true(all(ok(13, integer(0), z), ok(13, 12L, z)))
})

test_that("E3: null calibrates eta once and computes STAB only where it is compared", {
        x   <- .perf_graph(side = 10)
        res <- get_regions(x$g, regions_per_partition = 8, final_regions = 2, N_tree = 12,
                           seed = 1, verbose = FALSE)
        n_cal <- 0L
        orig  <- .calibrate_eta_empirical
        local_mocked_bindings(.calibrate_eta_empirical = function(...) {
                n_cal <<- n_cal + 1L; orig(...)
        })
        nb <- compute_null_baseline(x$g, final_regions = 2:3, regions_per_partition = 8,
                                    N_tree = 12, B = 3, stability_B = 2, observed = res,
                                    seed = 1, verbose = FALSE)
        expect_equal(n_cal, 1L)
        expect_true(is.finite(nb$settings$eta))
        expect_identical(nb$settings$stability_k, 2L)
        expect_true(all(is.na(nb$null$STAB[nb$null$final_regions == 3L])))
        expect_true(all(is.finite(nb$null$STAB[nb$null$final_regions == 2L])))
        s3 <- nb$summary[nb$summary$criterion == "STAB" & nb$summary$final_regions == 3L, ]
        expect_true(is.na(s3$null_mean) && is.na(s3$p) && s3$n_null == 0L)
        nb0 <- compute_null_baseline(x$g, final_regions = 2:3, regions_per_partition = 8,
                                     N_tree = 12, B = 1, stability_B = 2, seed = 1, verbose = FALSE)
        expect_identical(nb0$settings$stability_k, 2:3)
        nb1 <- compute_null_baseline(x$g, final_regions = 2:3, regions_per_partition = 8,
                                     N_tree = 12, B = 1, stability_B = 2, stability_k = 3,
                                     seed = 1, verbose = FALSE)
        expect_identical(nb1$settings$stability_k, 3L)
})

test_that("E3: sample_spanning_trees() takes a fixed eta and validates it", {
        x  <- .perf_graph(side = 8)
        a  <- sample_spanning_trees(x$g, n = 5, seed = 3, verbose = FALSE)
        b  <- sample_spanning_trees(x$g, n = 5, seed = 3, verbose = FALSE, eta = a$eta)
        expect_identical(b$eta, a$eta)
        expect_identical(sample_spanning_trees(x$g, n = 5, seed = 3, verbose = FALSE, eta = 0)$eta, 0)
        expect_error(sample_spanning_trees(x$g, n = 5, verbose = FALSE, eta = -1), "eta")
        expect_error(sample_spanning_trees(x$g, n = 5, verbose = FALSE, eta = c(1, 2)), "eta")
})

test_that("E4: a pulse_regions object stores the ensemble once", {
        x   <- .perf_graph(side = 10)
        res <- get_regions(x$g, regions_per_partition = 8, final_regions = 2, N_tree = 12,
                           seed = 1, verbose = FALSE)
        expect_null(res$fit$ensemble_meta$memb_mat)
        expect_equal(dim(res$ensemble), c(100L, 12L))
        expect_true(is.finite(res$fit$ensemble_meta$chance_coassoc))
        tu  <- get_regions(x$g, regions_per_partition = 8, final_regions = 2:3, N_tree = 12,
                           env = x$env, stability_B = 0, seed = 1, verbose = FALSE)
        expect_null(tu$fit$ensemble_meta$memb_mat)
        expect_equal(dim(tu$ensemble), c(100L, 12L))
        expect_s3_class(compute_validity(res, env = x$env, stability_B = 0), "data.frame")
})
