null_toy <- function(n = 10, seed = 1) {
  set.seed(seed)
  g  <- igraph::make_lattice(c(n, n))
  xy <- expand.grid(x = seq_len(n), y = seq_len(n))
  env <- cbind(e1 = (xy$x > n / 3) + (xy$x > 2 * n / 3) + stats::rnorm(n^2, 0, 0.3),
               e2 = stats::rnorm(n^2))
  ep <- igraph::ends(g, igraph::E(g), names = FALSE)
  igraph::E(g)$eid    <- seq_len(igraph::ecount(g))
  igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], ] - env[ep[, 2], ])^2)) + 1e-6
  list(g = g, env = env)
}

run_null <- function(tg, ...) {
  compute_null_baseline(tg$g, tg$env, regions_per_partition = 8, N_tree = 20,
                        stability_B = 0, seed = 1, verbose = FALSE,
                        cut_rule = "propensity", refine_boundary = FALSE, ...)
}

test_that("a non-config error in a null replicate aborts with context", {
  tg <- null_toy()
  local_mocked_bindings(.salso_k = function(...) stop("genuine bug"))
  expect_error(run_null(tg, final_regions = 2, B = 3),
               "Null replicate 1, k = 2: genuine bug")
})

test_that("too many config failures at one k stop the call early", {
  tg   <- null_toy()
  orig <- .salso_k
  local_mocked_bindings(.salso_k = function(memb, k, ...)
    if (k == 3L) .config_error("no consensus at k = 3") else orig(memb, k, ...))
  # 1 of 4 replicates = 25% > 20%: stops after replicate 1
  expect_error(run_null(tg, final_regions = 2:3, B = 4),
               "1 of 4 null replicates failed at k = 3 after replicate 1")
})

test_that("tolerated failures are reported and p counts them against observed", {
  tg  <- null_toy()
  res <- get_regions(tg$g, regions_per_partition = 8, final_regions = 2,
                     N_tree = 20, seed = 1, cut_rule = "propensity",
                     refine_boundary = FALSE, verbose = FALSE)
  orig  <- .salso_k
  calls <- 0L
  local_mocked_bindings(.salso_k = function(memb, k, ...) {
    calls <<- calls + 1L
    if (calls == 1L) .config_error("no consensus") else orig(memb, k, ...)
  })
  expect_warning(
    nb <- run_null(tg, final_regions = 2, B = 5, observed = res),   # 1 of 5 = 20%: allowed
    "k = 2 \\(1 of 5\\)")
  s <- nb$summary[nb$summary$criterion == "CH_FIT", ]
  expect_equal(s$n_fail, 1L)
  expect_equal(s$fail_share, 0.2)
  expect_equal(s$n_null, 4L)
  expect_true(is.finite(s$p) && is.finite(s$p_drop))
  expect_gte(s$p, s$p_drop)
  expect_equal(sum(nb$null$failed), 1L)
  expect_match(nb$null$error[nb$null$failed], "no consensus")
  expect_output(print(nb), "failed replicates: k = 2: 1")
})

test_that("without failures p equals p_drop", {
  tg  <- null_toy()
  res <- get_regions(tg$g, regions_per_partition = 8, final_regions = 2,
                     N_tree = 20, seed = 1, cut_rule = "propensity",
                     refine_boundary = FALSE, verbose = FALSE)
  nb <- run_null(tg, final_regions = 2, B = 3, observed = res)
  s  <- nb$summary[nb$summary$criterion == "CH_FIT", ]
  expect_equal(s$p, s$p_drop)
  expect_equal(s$n_fail, 0L)
  expect_output(print(nb), "failed replicates: none")
})

test_that("max_fail_share is validated", {
  tg <- null_toy()
  expect_error(run_null(tg, final_regions = 2, B = 2, max_fail_share = 1.5),
               "max_fail_share")
})
