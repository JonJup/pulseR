# Tests for suggest_regions() and its helpers (R/suggest_regions.R).

k_toy <- function(side = 10, R = 3, contrast = 3, seed = 1) {
  set.seed(seed)
  g  <- igraph::make_lattice(c(side, side))
  xy <- as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
  z  <- max.col(-as.matrix(stats::dist(rbind(xy[sample.int(nrow(xy), R), ], xy)))[-(1:R), 1:R])
  X  <- matrix(stats::rnorm(R * 3, sd = contrast), R, 3)[z, ] +
    matrix(stats::rnorm(length(z) * 3), ncol = 3)
  X  <- scale(X)
  g  <- igraph::set_graph_attr(g, "env", unname(X))
  ep <- igraph::ends(g, igraph::E(g), names = FALSE)
  igraph::E(g)$weight <- sqrt(rowSums((X[ep[, 1], ] - X[ep[, 2], ])^2))
  list(g = g, X = unname(X), z = z)
}
k_contig <- function(g, l) all(vapply(unique(l), function(r)
  igraph::components(igraph::induced_subgraph(g, which(l == r)))$no == 1L, logical(1)))

test_that("Ward path heights add up: W_1 = SST and W_k matches a direct count", {
  t <- k_toy()
  lab <- .superpixels(t$g, t$X, 20, 2)
  p <- .unit_path(lab, t$g, t$X)
  sst <- sum(scale(t$X, scale = FALSE)^2)
  expect_equal(p$W[1], sst, tolerance = 1e-10)
  direct <- vapply(seq_len(p$k0), function(k) {
    l <- .path_labels(p, k)
    sum((t$X - (rowsum(t$X, l) / tabulate(l))[l, ])^2)
  }, numeric(1))
  expect_equal(p$W, direct, tolerance = 1e-10)
  expect_equal(p$W[-p$k0] - p$W[-1], p$h[-1], tolerance = 1e-10)
})

test_that("path partitions are nested, have k regions and stay contiguous", {
  t <- k_toy()
  p <- .unit_path(.superpixels(t$g, t$X, 15, 2), t$g, t$X)
  for (k in c(1L, 2L, 5L, p$k0)) {
    l <- .path_labels(p, k)
    expect_equal(length(unique(l)), k)
    expect_true(k_contig(t$g, l))
  }
  fine <- .path_labels(p, 6L); coarse <- .path_labels(p, 3L)
  expect_true(all(tapply(coarse, fine, function(v) length(unique(v))) == 1L))
})

test_that("random contiguous partitions are contiguous with k regions", {
  t <- k_toy()
  adj <- lapply(igraph::adjacent_vertices(t$g, igraph::V(t$g)), as.integer)
  set.seed(2)
  for (i in 1:10) {
    l <- .random_contiguous(adj, 4L)
    expect_equal(length(unique(l)), 4L)
    expect_true(k_contig(t$g, l))
  }
})

test_that("spectral randomisation keeps means, covariance and Moran's I", {
  t <- k_toy()
  A <- .unit_adjacency(t$g); U <- .moran_basis(A)
  set.seed(3)
  Xs <- .spectral_randomise(t$X, U)
  n <- nrow(A)
  moran <- function(x) { x <- x - mean(x); n / sum(A) * sum(x * (A %*% x)) / sum(x^2) }
  expect_equal(colMeans(Xs), colMeans(t$X), tolerance = 1e-10)
  expect_equal(stats::cov(Xs), stats::cov(t$X), tolerance = 1e-10)
  expect_equal(apply(Xs, 2, moran), apply(t$X, 2, moran), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(Xs, t$X)))
})

test_that("slope heuristic finds a clean elbow", {
  W <- c(100, 60, 30, 10, 9.5, 9, 8.6, 8.2, 7.9, 7.6, 7.3, 7)
  expect_equal(.slope_heuristic(W, 1000)$k, 4L)
})

test_that("suggest_regions() validates its arguments", {
  t <- k_toy()
  expect_error(suggest_regions(t$g, k_range = c(2, 4, 5), B = 0), "contiguous")
  expect_error(suggest_regions(t$g, k_range = 2:4, k_base = 4, B = 0), "k_base")
  expect_warning(suggest_regions(t$g, k_range = 2:4,spectral_max = 50,
                               B = 1), "coarse")
  expect_warning(suggest_regions(t$g, k_range = 2:4, coarse = 30, B = 3, n_random = 10,
                                 n_boot = 10, seed = 1, verbose = FALSE), "smallest null")
})

test_that("coarse mode returns every rule inside k_range and rescoring is free", {
  t <- k_toy(side = 14)
  s <- suggest_regions(t$g, k_range = 2:6, coarse = 40, B = 19, n_random = 20,
                       n_boot = 20, seed = 1, verbose = FALSE)
  expect_s3_class(s, "pulse_k_suggestion")
  r <- s$rules[!is.na(s$rules)]
  expect_true(all(r >= 2L & r <= 6L))
  expect_equal(s$support$k, 2:6)
  expect_true(s$k %in% s$set)
  s2 <- .rescore_k(s, 3:5)
  expect_true(all(s2$rules[!is.na(s2$rules)] %in% 3:5))
  s_sp <- suggest_regions(t$g, k_range = 2:6, coarse = 40, B = 19,
                          n_random = 20, n_boot = 20, seed = 1, verbose = FALSE)
  expect_s3_class(s_sp, "pulse_k_suggestion")
  expect_output(print(s), "pulse_k_suggestion")
})

test_that("full pipeline runs end to end and is reproducible", {
  skip_on_cran()
  t <- k_toy(side = 12, R = 3)
  run <- function() suggest_regions(t$g, k_range = 2:5, N_tree = 12, B = 3,
                                    stability_B = 1, null_stability_B = 1,
                                    n_random = 10, n_boot = 10, min_size = 2,
                                    seed = 7, verbose = FALSE)
  a <- suppressWarnings(run()); b <- suppressWarnings(run())
  expect_identical(a$rules, b$rules)
  expect_equal(a$path$W, b$path$W)
  expect_equal(nrow(a$null$W), 3L)
  expect_true(all(is.finite(a$support$gap)))
  expect_equal(a$path$W[1], sum(scale(t$X, scale = FALSE)^2), tolerance = 1e-10)
})
