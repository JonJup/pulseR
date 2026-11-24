# =========================================================================== *
# suggest_regions.R
#
# Opt-in, null-calibrated suggestion for the number of Level 1 regions k.
#
#   1. One ensemble at a fixed resolution (regions_per_partition = 3 * k_max),
#      one consensus at k_base = 2 * k_max regions, contiguity repair, optional
#      boundary refinement. Every candidate k is then read off ONE nested
#      path: contiguity-constrained Ward merges of the base regions down to 1.
#      (Removes the rpp = 3k confound; every k is scored on the same object.)
#   2. B null replicates run the same pipeline on permuted (or spectrally
#      randomised) env and give a null path each.
#   3. Rules, all from the stored paths:
#        gap        null-calibrated gap on log(W_k / W_1), Tibshirani 1-SE rule
#        merge      fixed-sequence test of the Ward merge cost against the null
#        slope      slope heuristic / dimension jump (no null)
#        composite  CH (and STAB) z-scored against random contiguous partitions
#        stab_norm  path stability z-scored against the null (opt-in, costly)
#   4. Output: a support table per k, a plausible set, one point per rule.
#
# Scales to large N with `coarse`: superpixels from one greedy SSD cut of the
# minimum spanning tree, then the path and the null on the superpixel graph.
# =========================================================================== *

#' Suggest the number of Level 1 regions
#'
#' Opt-in helper for choosing \code{final_regions}. It does not replace a
#' user-chosen \eqn{k}: it reports how well each candidate \eqn{k} is
#' supported, a plausible set, and one suggestion per rule. Fit the chosen
#' \eqn{k} with \code{\link{get_regions}} afterwards.
#'
#' @section Fixed resolution:
#' One ensemble is built at \code{regions_per_partition} (default
#' \code{3 * max(k_range)}) and summarized once at \code{k_base} regions
#' (default \code{2 * max(k_range)}). Smaller \eqn{k} come from merging
#' adjacent base regions with the smallest Ward cost
#' \eqn{n_r n_s / (n_r + n_s) ||\bar x_r - \bar x_s||^2}, so the partitions are
#' nested and \eqn{W_{k-1} = W_k + h_k}, with \eqn{W_k} the within-region sum
#' of squares of \code{graph_attr(graph, "env")} and \eqn{h_k} the cost of the
#' merge from \eqn{k} to \eqn{k - 1}. Path partitions are not boundary-refined
#' below \code{k_base}; they are for choosing \eqn{k}, not for mapping.
#'
#' @section Rules:
#' \describe{
#'   \item{\code{gap_global}}{\eqn{\mathrm{Gap}(k) = \bar{\ell}^*_k - \ell_k} with
#'     \eqn{\ell_k = \log(W_k / W_1)} and \eqn{\bar\ell^*_k} the mean over null
#'     replicates; \eqn{s_k = \mathrm{sd}^*_k \sqrt{1 + 1/B}}. The suggestion
#'     is \code{cluster::maxSE(method = "Tibs2001SEmax")} over
#'     \code{k_range}: the smallest \eqn{k} with
#'     \eqn{\mathrm{Gap}(k) \ge \mathrm{Gap}(k+1) - s_{k+1}}. The 1-SE set is
#'     \eqn{\{k : \mathrm{Gap}(k) \ge \max \mathrm{Gap} - s_{\arg\max}\}}. Under the
#'     permutation null \eqn{W_1} is the same in every replicate, so this is
#'     the usual gap statistic.}
#'   \item{\code{merge}}{\eqn{p_k = (1 + \#\{b : h^*_{k,b} \ge h_k\}) / (B + 1)}.
#'     Splits are tested in the fixed order \eqn{k = 2, 3, \dots}; the
#'     suggestion is the last \eqn{k} before the first \eqn{p_k > \alpha}.
#'     A fixed-sequence test controls the family-wise error at \eqn{\alpha}
#'     without adjustment. With \eqn{B = 19} the smallest \eqn{p} is 0.05.}
#'   \item{\code{slope}}{Dimension jump on \eqn{C_k = W_k / N} with penalty
#'     shape \eqn{(k/N)(5 + 2 \log(N/k))} (Lebarbier 2005) over the whole path
#'     \eqn{1, \dots,} \code{k_base}; \code{slope_ddse} uses
#'     \code{capushe::DDSE} when installed. Needs no null. The penalty shape
#'     is proven for 1-D segmentation only; on a graph it is a heuristic.}
#'   \item{\code{composite}}{Within-\eqn{k} calibration (Akhanli and Hennig
#'     2020): \eqn{\log} CH of the path partition and, when
#'     \code{stability_B > 0}, its split-half stability, each z-scored against
#'     \code{n_random} random contiguous \eqn{k}-partitions (graph Voronoi
#'     cells of random seeds); mean z, argmax.}
#'   \item{\code{stab_norm}}{Path stability z-scored against the null; the
#'     largest \eqn{k} with \eqn{z \ge} \code{qnorm(1 - alpha)}. Needs
#'     \code{null_stability_B > 0}.}
#' }
#'
#' @section Nulls model:
#' Each column is
#' expanded on the Moran eigenvectors of the graph and the coefficients get
#' random signs, the same signs for every column (Wagner and Dray 2015,
#' "singleton" variant). This keeps every column's Moran's I, variance and the
#' covariance between columns. It needs a dense eigendecomposition, so it is
#' allowed only up to \code{spectral_max} units, or on superpixels.
#'
#' @section Large N:
#' \code{coarse = 500} replaces the ensemble by superpixels: one near-greedy
#' SSD cut of the minimum spanning tree into \code{coarse} pieces. The path,
#' the null (superpixel means permuted, or spectral on the superpixel graph)
#' and all rules then run on the superpixel graph, which costs seconds at any
#' \eqn{N}. The superpixel null conditions on the over-segmentation, so it is
#' less conservative than the full null; compare both on a subset before
#' relying on it.
#'
#' @param graph Weighted \code{igraph} from \code{\link{add_edge_weight}}
#'   (or with \code{\link{set_env}}); the environment matrix must be stored.
#' @param k_range Contiguous integer range of candidate \eqn{k}, e.g.
#'   \code{2:12}. \eqn{k = 1} may be included (no regions).
#' @param N_tree,weight_exponent,size_influence,cut_rule,min_size,target_edge_turnover
#'   As in \code{\link{get_regions}}.
#' @param regions_per_partition Fixed ensemble resolution. Default
#'   \code{3 * max(k_range)}.
#' @param k_base Regions of the base consensus, the top of the path. Default
#'   \code{2 * max(k_range)}; must exceed \code{max(k_range)}.
#' @param B Null replicates. Default 19 (smallest p-value 0.05). 0 skips the
#'   null; \code{gap}, \code{merge} and \code{stab_norm} are then \code{NA}.
#' @param coarse \code{NULL} (default, full pipeline) or the number of
#'   superpixels.
#' @param rule Rule behind \code{$k} and \code{$set}. Default \code{"gap"}.
#' @param alpha Level for \code{merge} and \code{stab_norm}. Default 0.05.
#' @param stability_B Split-half repeats for the observed path stability
#'   (used by \code{composite} and \code{stab_norm}). Each costs two base
#'   consensus fits. Default 5; 0 skips. Ignored with \code{coarse}.
#' @param null_stability_B The same for every null replicate; needed by
#'   \code{stab_norm}. Default 0.
#' @param n_random Random contiguous partitions per \eqn{k} for
#'   \code{composite}. Default 100.
#' @param refine_base Logical. Boundary refinement of the base partition.
#'   Default \code{TRUE}.
#' @param trees Optional \code{"pulse_trees"} object, e.g. \code{$trees} of a
#'   \code{\link{get_regions}} fit, to reuse its trees.
#' @param n_boot Bootstrap resamples for the support columns. Default 200.
#' @param spectral_max Largest number of units (or superpixels) for
#'   \code{null = "spectral"}. Default 5000.
#' @param large_n_threshold,sample_size,n_cores As in \code{\link{get_regions}}.
#' @param seed Optional integer; the caller's RNG state is restored.
#' @param verbose Logical.
#'
#' @return An object of class \code{"pulse_k_suggestion"}:
#'   \describe{
#'     \item{\code{k}}{Suggested \eqn{k} of \code{rule}.}
#'     \item{\code{set}}{Plausible set of \code{rule}. gap: every \eqn{k}
#'       chosen in at least 5\% of bootstrap resamples of the null replicates,
#'       plus \code{k}; merge: \code{k} up to the largest \eqn{k} with
#'       \eqn{p_k \le \alpha}; slope: the \eqn{k} chosen for
#'       \eqn{\kappa \in [1.5, 2.5] \kappa_{min}}; composite: 1-SE set.}
#'     \item{\code{rules}}{Named integer vector: the suggestion of every rule
#'       (\code{gap}, \code{gap_first}, \code{gap_global}, \code{merge},
#'       \code{slope}, \code{slope_ddse}, \code{composite}, \code{stab_norm}).}
#'     \item{\code{support}}{One row per \eqn{k} in \code{k_range}: \code{W},
#'       \code{gap}, \code{gap_se}, \code{gap_support} (bootstrap share of null
#'       resamples that select \eqn{k}), \code{in_gap_set} (global 1-SE set),
#'       \code{merge_cost},
#'       \code{merge_z}, \code{merge_p}, \code{ch}, \code{stab},
#'       \code{z_ch}, \code{z_stab}, \code{composite}, \code{composite_se},
#'       \code{composite_support}, \code{stab_null_z}, \code{vi_next}
#'       (VI in bits between the path partitions at \eqn{k} and \eqn{k + 1}).}
#'     \item{\code{path}}{Base labels and merge sequence; partitions at any
#'       \eqn{k} with \code{pulseR:::.path_labels(x$path, k)}.}
#'     \item{\code{null}}{\code{W} and \code{h}: \eqn{B \times}
#'       \code{k_base} matrices of the null paths; \code{stab}.}
#'     \item{\code{random}}{Random-partition calibration draws.}
#'     \item{\code{settings}, \code{timing}}{Call settings; seconds per stage.}
#'   }
#'
#' @references
#' Tibshirani, Walther and Hastie (2001) JRSS-B 63:411-423.
#' Akhanli and Hennig (2020) Statistics and Computing 30:1523-1544.
#' Baudry, Maugis and Michel (2012) Statistics and Computing 22:455-470.
#' Lebarbier (2005) Signal Processing 85:717-736.
#' Wagner and Dray (2015) Methods in Ecology and Evolution 6:1169-1178.
#'
#' @examples
#' \dontrun{
#' s <- suggest_regions(gw, k_range = 2:12, N_tree = 200, B = 19, seed = 1)
#' s
#' fit <- get_regions(gw, final_regions = s$k, N_tree = 200, seed = 1)
#' }
#' @export
suggest_regions <- function(graph,
                            k_range,
                            N_tree                = 200,
                            regions_per_partition = NULL,
                            k_base                = NULL,
                            B                     = 19,
                            coarse                = NULL,
                            rule                  = c("gap_global", "merge", "slope",
                                                      "composite", "stab_norm"),
                            alpha                 = 0.05,
                            stability_B           = 5,
                            null_stability_B      = 0,
                            n_random              = 100,
                            weight_exponent       = 4,
                            size_influence        = 0,
                            cut_rule              = c("ssd", "propensity"),
                            refine_base           = TRUE,
                            min_size              = 5,
                            target_edge_turnover  = 0.1,
                            trees                 = NULL,
                            n_boot                = 200,
                            spectral_max          = 5000,
                            large_n_threshold     = 20000,
                            sample_size           = 5000,
                            n_cores               = 1,
                            seed                  = NULL,
                            verbose               = TRUE) {
  stopifnot(inherits(graph, "igraph"))
  null     <- match.arg(null)
  rule     <- match.arg(rule)
  cut_rule <- match.arg(cut_rule)
  X <- igraph::graph_attr(graph, "env")
  if (is.null(X))
    stop("suggest_regions() needs the environment matrix on `graph`: build it ",
         "with add_edge_weight(), or store one with set_env().", call. = FALSE)
  X <- as.matrix(X)
  if (anyNA(X)) stop("The environment matrix has missing values.", call. = FALSE)
  N <- nrow(X)
  k_range <- sort(unique(as.integer(k_range)))
  if (length(k_range) < 2L || any(diff(k_range) != 1L) || k_range[1] < 1L)
    stop("`k_range` must be a contiguous integer range such as 2:12.", call. = FALSE)
  k_max  <- max(k_range)
  k_base <- as.integer(if (is.null(k_base)) 2L * k_max else k_base)
  if (k_base <= k_max)
    stop("`k_base` (", k_base, ") must exceed max(k_range) (", k_max, ").", call. = FALSE)
  rpp <- as.integer(if (is.null(regions_per_partition)) 3L * k_max else regions_per_partition)
  if (is.null(coarse) && rpp < k_base)
    stop("`regions_per_partition` (", rpp, ") must be at least `k_base` (", k_base,
         ").", call. = FALSE)
  if (!is.null(coarse)) {
    coarse <- as.integer(coarse)
    if (length(coarse) != 1L || coarse <= k_base || coarse >= N)
      stop("`coarse` must be one integer between k_base and the number of units.",
           call. = FALSE)
  }
  n_nodes <- if (is.null(coarse)) N else coarse
  if (n_nodes > spectral_max){
    coarse = 500
    message(sprintf("The number of units exceeds the suggested maximum number of %s. 
                    Therefore, %s superpixels (argument coarse) are used.", spectral_max, coarse))
  }
    
  if (B > 0L && is.null(coarse)) .check_env_reproduces_weights(graph, X)
  if (B > 0L && (1 / (B + 1)) > alpha)
    warning(sprintf(paste0("With B = %d the smallest null p-value is %.3f > alpha = %g, ",
                           "so the merge rule cannot accept any split. Use B >= %d."),
                    B, 1 / (B + 1), alpha, ceiling(1 / alpha - 1)), call. = FALSE)
  .local_seed(seed)
  tm <- list()
  t0 <- proc.time()[["elapsed"]]
  tick <- function(name) {
    t1 <- proc.time()[["elapsed"]]; tm[[name]] <<- t1 - t0; t0 <<- t1
  }
  ens_args <- list(N_tree = N_tree, regions_per_partition = rpp,
                   weight_exponent = weight_exponent, size_influence = size_influence,
                   cut_rule = cut_rule, min_size = min_size)
  base_args <- list(k_base = k_base, refine = refine_base, min_size = min_size,
                    large_n_threshold = large_n_threshold, sample_size = sample_size,
                    n_cores = n_cores)

  # ---- observed path ------------------------------------------------------*
  stab_ks <- seq_len(k_max)
  if (is.null(coarse)) {
    if (is.null(trees)) {
      if (verbose) message(sprintf("Sampling %d spanning trees ...", N_tree))
      trees <- sample_spanning_trees(graph, n = N_tree, verbose = FALSE,
                                     target_edge_turnover = target_edge_turnover)
    }
    if (trees$n_trees < N_tree)
      stop("`trees` holds fewer than N_tree trees.", call. = FALSE)
    if (verbose) message(sprintf(
      "Ensemble at regions_per_partition = %d, base consensus at k_base = %d ...",
      rpp, k_base))
    memb <- do.call(.compute_ensemble_internal,
                    c(list(graph = graph, trees = trees), ens_args))
    base <- do.call(.base_partition, c(list(memb = memb, graph = graph, X = X), base_args))
    path <- .unit_path(base, graph, X)
    tick("observed")
    stab <- if (stability_B > 0L) {
      if (verbose) message(sprintf("Path stability (B = %d) ...", stability_B))
      .path_stability(memb, graph, X, stability_B, stab_ks, base_args)
    } else rep(NA_real_, k_max)
    tick("stability")
    rm(memb)
  } else {
    if (verbose) message(sprintf("Superpixels: %d pieces ...", coarse))
    sp   <- .superpixels(graph, X, coarse, min_size)
    path <- .unit_path(sp, graph, X)
    stab <- rep(NA_real_, k_max)
    tick("observed")
  }
  if (path$k0 <= k_max)
    stop(sprintf(paste0("The base partition has only %d regions after contiguity ",
                        "repair; raise k_base or regions_per_partition."), path$k0),
         call. = FALSE)

  # ---- null paths ---------------------------------------------------------*
  W0 <- matrix(NA_real_, B, path$k0)
  H0 <- W0
  S0 <- matrix(NA_real_, B, k_max)
  n_fail <- 0L
  if (B > 0L) {
    U <- if (is.null(coarse)){
      .moran_basis(.unit_adjacency(graph))
    } else {
      .moran_basis(.region_adjacency_matrix(path$region))
    }
    eta_null <- NULL
    for (b in seq_len(B)) {
      if (verbose) message(sprintf("Null replicate %d / %d ...", b, B))
      if (is.null(coarse)) {
        X_b <- .spectral_randomise(X, U)
        g_b <- graph
        igraph::E(g_b)$weight <- .edge_weights_from_env(graph, X_b)
        g_b <- igraph::set_graph_attr(g_b, "env", X_b)
        tr  <- sample_spanning_trees(g_b, n = N_tree, verbose = FALSE,
                                     target_edge_turnover = target_edge_turnover,
                                     eta = eta_null)
        if (is.null(eta_null)) eta_null <- tr$eta
        pb <- tryCatch({
          mb <- do.call(.compute_ensemble_internal,
                        c(list(graph = g_b, trees = tr), ens_args))
          bb <- do.call(.base_partition, c(list(memb = mb, graph = g_b, X = X_b), base_args))
          p  <- .unit_path(bb, g_b, X_b)
          s  <- if (null_stability_B > 0L)
            .path_stability(mb, g_b, X_b, null_stability_B, stab_ks, base_args) else NULL
          list(p = p, s = s)
        }, pulseR_config_error = function(e) NULL)
      } else {
        rg  <- path$region
        M   <- rg$S / rg$n
        M_b <- spectral_randomise(M, U)
        p   <- .ward_path(rg$pairs, rg$n, M_b * rg$n, rg$W_within)
        pb  <- list(p = p, s = NULL)
      }
      if (is.null(pb) || pb$p$k0 < k_max) { n_fail <- n_fail + 1L; next }
      kk <- seq_len(min(path$k0, pb$p$k0))
      W0[b, kk] <- pb$p$W[kk]
      H0[b, kk] <- pb$p$h[kk]
      if (!is.null(pb$s)) S0[b, ] <- pb$s
    }
    if (n_fail > 0L)
      warning(sprintf(paste0("%d of %d null replicates could not build a path to ",
                             "k_max = %d and were dropped."), n_fail, B, k_max),
              call. = FALSE)
  }
  tick("null")

  # ---- random contiguous partitions (composite) ---------------------------*
  rnd <- .random_calibration(graph, X, seq.int(max(2L, k_range[1]), k_max), n_random)
  tick("random")

  x <- structure(list(
    k = NA_integer_, set = integer(0), rules = NULL, support = NULL,
    path = path, stab = stab,
    null = list(W = W0, h = H0, stab = S0, type = null, n_fail = n_fail),
    random = rnd,
    settings = list(k_range = k_range, k_base = k_base, regions_per_partition = rpp,
                    N_tree = N_tree, B = B, null = null, coarse = coarse, rule = rule,
                    alpha = alpha, stability_B = stability_B,
                    null_stability_B = null_stability_B, n_random = n_random,
                    weight_exponent = weight_exponent, cut_rule = cut_rule,
                    refine_base = refine_base, min_size = min_size, n_boot = n_boot,
                    N = N),
    timing = unlist(tm)),
    class = "pulse_k_suggestion")
  .rescore_k(x, k_range)
}

#' Re-derive every rule for another candidate range, without refitting
#'
#' The stored paths and nulls cover \eqn{k = 1, \dots,} \code{k_base}, so any
#' contiguous \code{k_range} with \code{max(k_range) <} \code{k_base} (and
#' \code{max(k_range)} at most the original maximum, for the composite and
#' stability columns) can be scored again for free. Used by the benchmark to
#' score centred, offset and wide grids on one fit.
#' @keywords internal
#' @noRd
.rescore_k <- function(x, k_range) {
  st <- x$settings
  k_range <- sort(unique(as.integer(k_range)))
  stopifnot(length(k_range) >= 2L, all(diff(k_range) == 1L),
            max(k_range) < x$path$k0, max(k_range) <= max(st$k_range))
  r <- .k_rules(x, k_range)
  x$rules   <- r$rules
  x$support <- r$support
  x$sets    <- r$sets
  x$k       <- unname(r$rules[[st$rule]])
  x$set     <- r$sets[[st$rule]]
  x$settings$k_range_scored <- k_range
  x
}

# --------------------------------------------------------------------------- *
# Rules
# --------------------------------------------------------------------------- *

.k_rules <- function(x, k_range) {
  st   <- x$settings
  N    <- st$N
  W    <- x$path$W                       # length k0, W[k] = within SSD at k
  h    <- x$path$h                       # h[k] = cost of merge k -> k - 1 (h[1] = NA)
  K    <- length(k_range)
  clamp <- function(k) if (is.na(k)) NA_integer_ else
    as.integer(min(max(k, k_range[1]), k_range[K]))
  W0 <- x$null$W; H0 <- x$null$h
  ok <- if (nrow(W0)) stats::complete.cases(W0[, k_range, drop = FALSE]) else logical(0)
  W0 <- W0[ok, , drop = FALSE]; H0 <- H0[ok, , drop = FALSE]
  Bn <- nrow(W0)
  sup <- data.frame(k = k_range, W = W[k_range])

  # ---- gap ------------------------------------------------------------------*
  gap_pick <- function(Wn) {
    ln  <- log(Wn[, k_range, drop = FALSE] / Wn[, 1L])
    g   <- colMeans(ln) - log(W[k_range] / W[1L])
    s   <- apply(ln, 2, stats::sd) * sqrt(1 + 1 / nrow(Wn))
    s[!is.finite(s)] <- 0
    list(gap = g, se = s,
         tibs   = k_range[cluster::maxSE(g, s, method = "Tibs2001SEmax")],
         first  = k_range[cluster::maxSE(g, s, method = "firstSEmax")],
         global = k_range[cluster::maxSE(g, s, method = "globalSEmax")])
  }
  if (Bn >= 2L) {
    gp <- gap_pick(W0)
    sup$gap <- gp$gap; sup$gap_se <- gp$se
    best <- which.max(gp$gap)
    sup$in_gap_set <- gp$gap >= gp$gap[best] - gp$se[best]
    bt <- vapply(seq_len(st$n_boot), function(i)
      gap_pick(W0[sample.int(Bn, Bn, replace = TRUE), , drop = FALSE])$tibs, integer(1))
    sup$gap_support <- tabulate(match(bt, k_range), K) / st$n_boot
    k_gap <- gp$tibs; k_gap_first <- gp$first; k_gap_global <- gp$global
  } else {
    sup$gap <- sup$gap_se <- sup$gap_support <- NA_real_; sup$in_gap_set <- NA
    k_gap <- k_gap_first <- k_gap_global <- NA_integer_
  }

  # ---- merge-cost test ------------------------------------------------------*
  # Fixed sequence k = 2, 3, ...: accept k while p_k <= alpha.
  k_seq <- seq.int(2L, max(k_range))
  if (Bn >= 1L) {
    hp <- vapply(k_seq, function(k) (1 + sum(H0[, k] >= h[k])) / (Bn + 1), numeric(1))
    hz <- vapply(k_seq, function(k) {
      s <- stats::sd(H0[, k]); if (is.finite(s) && s > 0) (h[k] - mean(H0[, k])) / s else NA_real_
    }, numeric(1))
    acc  <- hp <= st$alpha
    stop_at <- which(!acc)[1L]
    k_merge <- if (is.na(stop_at)) max(k_seq) else k_seq[stop_at] - 1L
    sup$merge_p <- hp[match(k_range, k_seq)]
    sup$merge_z <- hz[match(k_range, k_seq)]
    sig  <- k_seq[acc]
    merge_set <- seq.int(k_merge, max(c(k_merge, sig)))
  } else {
    sup$merge_p <- sup$merge_z <- NA_real_
    k_merge <- NA_integer_; merge_set <- integer(0)
  }
  sup$merge_cost <- h[k_range]

  # ---- slope heuristic ------------------------------------------------------*
  sl <- .slope_heuristic(W, N)

  # ---- composite (within-k calibrated CH and STAB) ---------------------------*
  rnd <- x$random
  kc  <- intersect(k_range, rnd$k)
  ch_obs <- .ch_from_W(W, N)
  sup$ch   <- ch_obs[k_range]
  sup$stab <- x$stab[k_range]
  zc <- function(rc, ra) {
    zch <- vapply(kc, function(k) {
      v <- log(rnd$ch[rc, as.character(k)]); v <- v[is.finite(v)]
      if (length(v) > 1L && stats::sd(v) > 0) (log(ch_obs[k]) - mean(v)) / stats::sd(v) else NA_real_
    }, numeric(1))
    zst <- vapply(kc, function(k) {
      v <- rnd$ari[ra, as.character(k)]; v <- v[is.finite(v)]
      if (is.finite(x$stab[k]) && length(v) > 1L && stats::sd(v) > 0)
        (x$stab[k] - mean(v)) / stats::sd(v) else NA_real_
    }, numeric(1))
    A <- rowMeans(cbind(zch, zst), na.rm = TRUE); A[is.nan(A)] <- NA_real_
    list(zch = zch, zst = zst, A = A)
  }
  if (length(kc) >= 2L) {
    nr <- nrow(rnd$ch); na <- nrow(rnd$ari)
    z0 <- zc(seq_len(nr), seq_len(na))
    Ab <- vapply(seq_len(st$n_boot), function(i)
      zc(sample.int(nr, nr, replace = TRUE), sample.int(na, na, replace = TRUE))$A,
      numeric(length(kc)))
    Ab <- matrix(Ab, nrow = length(kc))
    pick <- apply(Ab, 2, function(a) if (all(is.na(a))) NA_integer_ else kc[which.max(a)])
    m <- match(k_range, kc)
    sup$z_ch <- z0$zch[m]; sup$z_stab <- z0$zst[m]; sup$composite <- z0$A[m]
    sup$composite_se <- apply(Ab, 1, stats::sd, na.rm = TRUE)[m]
    sup$composite_support <- (tabulate(match(pick, kc), length(kc)) / st$n_boot)[m]
    k_comp <- if (all(is.na(z0$A))) NA_integer_ else kc[which.max(z0$A)]
    jb <- which(kc == k_comp)
    comp_set <- kc[!is.na(z0$A) &
                     z0$A >= z0$A[jb] - apply(Ab, 1, stats::sd, na.rm = TRUE)[jb]]
  } else {
    sup$z_ch <- sup$z_stab <- sup$composite <- sup$composite_se <-
      sup$composite_support <- NA_real_
    k_comp <- NA_integer_; comp_set <- integer(0)
  }

  # ---- null-normalised stability -------------------------------------------*
  S0 <- x$null$stab
  S0 <- S0[stats::complete.cases(S0[, k_range[k_range > 1L], drop = FALSE]), , drop = FALSE]
  if (nrow(S0) >= 2L && any(is.finite(x$stab))) {
    zs <- vapply(k_range, function(k) {
      s <- stats::sd(S0[, k])
      if (k <= ncol(S0) && is.finite(s) && s > 0) (x$stab[k] - mean(S0[, k])) / s else NA_real_
    }, numeric(1))
    sup$stab_null_z <- zs
    pass <- k_range[!is.na(zs) & zs >= stats::qnorm(1 - st$alpha)]
    k_stab <- if (length(pass)) max(pass) else k_range[1]
    stab_set <- pass
  } else {
    sup$stab_null_z <- NA_real_; k_stab <- NA_integer_; stab_set <- integer(0)
  }

  # ---- change between adjacent path partitions -------------------------------*
  sup$vi_next <- vapply(k_range, function(k)
    if (k + 1L <= x$path$k0) .vi_bits(.path_labels(x$path, k), .path_labels(x$path, k + 1L))
    else NA_real_, numeric(1))

  rules <- c(gap = clamp(k_gap), gap_first = clamp(k_gap_first),
             gap_global = clamp(k_gap_global), merge = clamp(k_merge),
             slope = clamp(sl$k), slope_ddse = clamp(sl$k_ddse),
             composite = clamp(k_comp), stab_norm = clamp(k_stab))
  # Gap set: k with at least 5 % bootstrap support for the Tibshirani rule,
  # plus the rule's own pick. It reflects Monte Carlo uncertainty of the null,
  # not sampling uncertainty of the data; the global 1-SE set is in_gap_set.
  gap_set <- if (all(is.na(sup$gap_support))) integer(0) else
    sort(unique(c(k_range[sup$gap_support >= 0.05], clamp(k_gap))))
  sets <- list(gap = gap_set,
               merge = intersect(merge_set, k_range),
               slope = intersect(sl$set, k_range),
               composite = comp_set, stab_norm = stab_set)
  list(rules = rules, support = sup, sets = sets,
       unclamped = c(gap = k_gap, merge = k_merge, slope = sl$k, slope_ddse = sl$k_ddse))
}

# CH of the path partitions from the merge heights: SSB_k = W_1 - W_k.
.ch_from_W <- function(W, N) {
  k <- seq_along(W)
  ch <- ((W[1L] - W) / (k - 1)) / (W / (N - k))
  ch[1L] <- NA_real_
  ch
}

# Dimension jump with the Lebarbier penalty shape on the whole path.
.slope_heuristic <- function(W, N, n_kappa = 2000L) {
  k   <- seq_along(W)
  C   <- W / N
  pen <- (k / N) * (5 + 2 * log(N / k))
  out <- list(k = NA_integer_, k_ddse = NA_integer_, set = integer(0), kappa_min = NA_real_)
  if (length(k) < 4L) return(out)
  # kappa range: from 0 (largest model wins) to the value at which k = 1 wins.
  hi  <- max((C[1L] - C[-1L]) / (pen[-1L] - pen[1L]), na.rm = TRUE) * 4
  kap <- exp(seq(log(hi * 1e-6), log(hi), length.out = n_kappa))
  sel <- vapply(kap, function(a) which.min(C + a * pen), integer(1))
  drop <- -diff(sel)
  j <- which.max(drop)                       # largest jump in complexity
  kmin <- kap[j + 1L]
  out$kappa_min <- kmin
  out$k   <- which.min(C + 2 * kmin * pen)
  out$set <- sort(unique(vapply(c(1.5, 1.75, 2, 2.25, 2.5) * kmin,
                                function(a) which.min(C + a * pen), integer(1))))
  if (requireNamespace("capushe", quietly = TRUE)) {
    d <- data.frame(model = as.character(k), pen = pen, complexity = k, contrast = C)
    out$k_ddse <- tryCatch(as.integer(capushe::DDSE(d)@model),
                           error = function(e) NA_integer_,
                           warning = function(w) suppressWarnings(
                             tryCatch(as.integer(capushe::DDSE(d)@model),
                                      error = function(e) NA_integer_)))
  }
  out
}

# Variation of information in bits between two labelings.
.vi_bits <- function(a, b) {
  n  <- length(a)
  tj <- table(a, b) / n
  pa <- rowSums(tj); pb <- colSums(tj)
  H  <- function(p) -sum(p[p > 0] * log2(p[p > 0]))
  2 * H(as.numeric(tj)) - H(pa) - H(pb)
}

# --------------------------------------------------------------------------- *
# Paths
# --------------------------------------------------------------------------- *

#' Base consensus at k_base: salso, contiguity repair, optional refinement.
#' @keywords internal
#' @noRd
.base_partition <- function(memb, graph, X, k_base, refine = TRUE, min_size = 5L,
                            large_n_threshold = 20000, sample_size = 5000L,
                            n_cores = 1L) {
  fit  <- .consensus_partition(memb, k_base, large_n_threshold, sample_size,
                               refine_iter = 10L, n_runs = 16L, n_cores = n_cores,
                               verbose = FALSE)
  st   <- .ensemble_stats(memb, fit$labels, k_base, pi = TRUE)
  hard <- .enforce_contiguity(graph, fit$labels, st$Pi)
  hard <- match(hard, sort(unique(hard)))
  if (refine) {
    ms <- max(1L, min(as.integer(min_size), nrow(X) %/% max(hard)))
    hard <- .refine_regions(hard, graph, X, min_size = ms)$z
  }
  match(hard, sort(unique(hard)))
}

#' Region statistics of a unit labelling, then the Ward path.
#' @keywords internal
#' @noRd
.unit_path <- function(lab, graph, X) {
  lab <- match(lab, sort(unique(lab)))
  k0  <- max(lab)
  n   <- tabulate(lab, k0)
  S   <- rowsum(X, lab, reorder = TRUE)
  Wwithin <- sum((X - (S / n)[lab, , drop = FALSE])^2)
  ep  <- igraph::ends(graph, igraph::E(graph), names = FALSE)
  a <- lab[ep[, 1]]; b <- lab[ep[, 2]]
  cut <- a != b
  pairs <- unique(cbind(pmin(a[cut], b[cut]), pmax(a[cut], b[cut])))
  p <- .ward_path(pairs, n, S, Wwithin)
  p$base   <- lab
  p$region <- list(pairs = pairs, n = n, S = S, W_within = Wwithin)
  p
}

#' Contiguity-constrained Ward agglomeration of regions
#'
#' @param pairs Two-column matrix of adjacent region pairs (i < j).
#' @param n Region sizes. @param S Region column sums (k0 x p).
#' @param W_within Within-region SSD of the starting regions.
#' @return list(k0, merges = k0 - 1 rows (from, into, cost), W (length k0,
#'   W[k] at k regions), h (h[k] = cost of the merge k -> k - 1, h[1] = NA)).
#'   Costs are not forced to be monotone (constrained Ward can invert).
#' @keywords internal
#' @noRd
.ward_path <- function(pairs, n, S, W_within) {
  k0 <- length(n)
  n  <- as.numeric(n)
  S  <- as.matrix(S)
  A  <- matrix(FALSE, k0, k0)
  A[pairs] <- TRUE; A[pairs[, 2:1, drop = FALSE]] <- TRUE
  alive <- rep(TRUE, k0)
  cost_row <- function(i) {
    j <- which(A[i, ] & alive)
    if (!length(j)) return(list(j = integer(0), c = numeric(0)))
    mi <- S[i, ] / n[i]
    Mj <- S[j, , drop = FALSE] / n[j]
    d2 <- rowSums((Mj - rep(mi, each = length(j)))^2)
    list(j = j, c = n[i] * n[j] / (n[i] + n[j]) * d2)
  }
  Cm <- matrix(Inf, k0, k0)
  for (i in seq_len(k0)) { cr <- cost_row(i); Cm[i, cr$j] <- cr$c }
  merges <- matrix(NA_real_, k0 - 1L, 3L, dimnames = list(NULL, c("from", "into", "cost")))
  W <- numeric(k0); W[k0] <- W_within
  h <- rep(NA_real_, k0)
  for (s in seq_len(k0 - 1L)) {
    m <- which.min(Cm)
    if (!length(m) || !is.finite(Cm[m]))
      stop("Region graph is disconnected; cannot merge to one region.", call. = FALSE)
    i <- (m - 1L) %% k0 + 1L; j <- (m - 1L) %/% k0 + 1L
    a <- min(i, j); b <- max(i, j)              # b merges into a
    k_now <- k0 - s + 1L
    h[k_now] <- Cm[m]
    W[k_now - 1L] <- W[k_now] + Cm[m]
    merges[s, ] <- c(b, a, Cm[m])
    S[a, ] <- S[a, ] + S[b, ]; n[a] <- n[a] + n[b]
    A[a, ] <- A[a, ] | A[b, ]; A[, a] <- A[a, ]; A[a, a] <- FALSE
    alive[b] <- FALSE
    Cm[b, ] <- Inf; Cm[, b] <- Inf
    cr <- cost_row(a)
    Cm[a, ] <- Inf; Cm[, a] <- Inf
    Cm[a, cr$j] <- cr$c; Cm[cr$j, a] <- cr$c
  }
  list(k0 = k0, merges = merges, W = W, h = h)
}

#' Unit labels of the path partition with k regions (1..k).
#' @keywords internal
#' @noRd
.path_labels <- function(path, k) {
  stopifnot(k >= 1L, k <= path$k0)
  map <- seq_len(path$k0)
  n_m <- path$k0 - k
  if (n_m > 0L) for (s in seq_len(n_m)) {
    map[map == path$merges[s, "from"]] <- path$merges[s, "into"]
  }
  lab <- map[path$base]
  match(lab, sort(unique(lab)))
}

#' Path stability: split-half ensemble agreement at every k on one path each.
#' Failed halves count as ARI 0 (as STAB in tune_regions()).
#' @keywords internal
#' @noRd
.path_stability <- function(memb, graph, X, B, ks, base_args) {
  M <- ncol(memb); half <- M %/% 2L
  half_path <- function(cols) tryCatch({
    a <- base_args; a$refine <- FALSE
    bb <- do.call(.base_partition, c(list(memb = memb[, cols, drop = FALSE],
                                          graph = graph, X = X), a))
    .unit_path(bb, graph, X)
  }, pulseR_config_error = function(e) NULL)
  res <- matrix(0, B, length(ks))
  for (b in seq_len(B)) {
    perm <- sample.int(M)
    pa <- half_path(sort(perm[seq_len(half)]))
    pb <- half_path(sort(perm[half + seq_len(half)]))
    if (is.null(pa) || is.null(pb)) next
    res[b, ] <- vapply(ks, function(k)
      if (k <= min(pa$k0, pb$k0))
        .adjusted_rand_index(.path_labels(pa, k), .path_labels(pb, k)) else 0,
      numeric(1))
  }
  out <- colMeans(res)
  out[ks == 1L] <- NA_real_            # one region is trivially stable
  out
}

#' Superpixels: near-greedy SSD cut of the minimum spanning tree.
#' @keywords internal
#' @noRd
.superpixels <- function(graph, X, n_sp, min_size = 5L) {
  tr <- igraph::mst(graph, weights = igraph::E(graph)$weight)
  ep <- igraph::ends(tr, igraph::E(tr), names = FALSE)
  ms <- max(1L, min(as.integer(min_size), nrow(X) %/% (2L * n_sp)))
  .ssd_cut_tree(ep, X, n_sp - 1L, weight_exponent = 50, size_influence = 0,
                min_size = ms)$labels
}

# --------------------------------------------------------------------------- *
# Nulls
# --------------------------------------------------------------------------- *

.check_env_reproduces_weights <- function(graph, X) {
  r <- suppressWarnings(stats::cor(.edge_weights_from_env(graph, X),
                                   igraph::E(graph)$weight))
  if (!is.finite(r) || r < 0.999)
    stop(sprintf(paste0(
      "Euclidean distances on the stored env do not reproduce E(graph)$weight ",
      "(correlation %.4f, need > 0.999), so the null would not test ",
      "environmental signal. Use B = 0 (no null; slope and composite only) ",
      "or rebuild the weights with add_edge_weight()."), r), call. = FALSE)
  invisible(r)
}

.unit_adjacency <- function(graph) {
  A <- as.matrix(igraph::as_adjacency_matrix(graph, sparse = TRUE, attr = NULL))
  A[A > 0] <- 1
  diag(A) <- 0
  A
}

.region_adjacency_matrix <- function(region) {
  k <- length(region$n)
  A <- matrix(0, k, k)
  A[region$pairs] <- 1; A[region$pairs[, 2:1, drop = FALSE]] <- 1
  A
}

#' Moran eigenvector basis orthogonal to the constant vector.
#' @keywords internal
#' @noRd
.moran_basis <- function(A) {
  n  <- nrow(A)
  Cn <- diag(n) - 1 / n
  O  <- Cn %*% A %*% Cn
  # Push the constant vector to the top eigenvalue, then drop it, so every
  # remaining eigenvector is exactly orthogonal to 1.
  O  <- O + (max(rowSums(abs(A))) + 1) * 2 * matrix(1 / n, n, n)
  e  <- eigen((O + t(O)) / 2, symmetric = TRUE)
  e$vectors[, -1L, drop = FALSE]
}

#' Spectral randomisation: random signs on the Moran coefficients, the same
#' sign for every column. Keeps column means, variances, covariances and
#' Moran's I exactly.
#' @keywords internal
#' @noRd
.spectral_randomise <- function(X, U) {
  mu <- colMeans(X)
  Xc <- X - rep(mu, each = nrow(X))
  Cf <- crossprod(U, Xc)
  s  <- sample(c(-1, 1), nrow(Cf), replace = TRUE)
  U %*% (s * Cf) + rep(mu, each = nrow(X))
}

# --------------------------------------------------------------------------- *
# Random contiguous partitions (composite calibration)
# --------------------------------------------------------------------------- *

#' Graph Voronoi partition of k random seeds by synchronous multi-source BFS.
#' Every region is a subtree of the BFS forest, so it is contiguous.
#' @keywords internal
#' @noRd
.random_contiguous <- function(adj, k) {
  N <- length(adj)
  lab <- integer(N)
  seeds <- sample.int(N, k)
  lab[seeds] <- seq_len(k)
  front <- seeds
  while (length(front)) {
    front <- front[sample.int(length(front))]
    nb <- adj[front]
    to <- unlist(nb, use.names = FALSE)
    fl <- rep(lab[front], lengths(nb))
    keep <- lab[to] == 0L
    to <- to[keep]; fl <- fl[keep]
    first <- !duplicated(to)
    to <- to[first]; lab[to] <- fl[first]
    front <- to
  }
  lab[lab == 0L] <- 1L                       # unreachable units (disconnected graph)
  lab
}

.random_calibration <- function(graph, X, ks, n_random) {
  if (n_random < 2L || !length(ks))
    return(list(k = integer(0), ch = matrix(NA_real_, 0, 0), ari = matrix(NA_real_, 0, 0)))
  adj <- lapply(igraph::adjacent_vertices(graph, igraph::V(graph)), as.integer)
  ch  <- matrix(NA_real_, n_random, length(ks), dimnames = list(NULL, ks))
  ari <- matrix(NA_real_, n_random %/% 2L, length(ks), dimnames = list(NULL, ks))
  for (j in seq_along(ks)) {
    labs <- lapply(seq_len(n_random), function(i) .random_contiguous(adj, ks[j]))
    ch[, j] <- vapply(labs, function(l) .pseudo_F_env(X, l, ks[j])[["CH_FIT"]], numeric(1))
    for (i in seq_len(n_random %/% 2L))
      ari[i, j] <- .adjusted_rand_index(labs[[2L * i - 1L]], labs[[2L * i]])
  }
  list(k = as.integer(ks), ch = ch, ari = ari)
}

# --------------------------------------------------------------------------- *
# Print
# --------------------------------------------------------------------------- *

#' @rdname suggest_regions
#' @param x A \code{"pulse_k_suggestion"} object.
#' @param ... Unused.
#' @return \code{print.pulse_k_suggestion} returns \code{x} invisibly.
#' @export
print.pulse_k_suggestion <- function(x, ...) {
  st <- x$settings
  cat(sprintf("<pulse_k_suggestion>  k = %s (%s) | plausible set {%s}\n",
              format(x$k), st$rule, paste(x$set, collapse = ", ")))
  cat(sprintf("  k_range %d..%d | base %d regions | %s null, B = %d%s | %s\n",
              min(st$k_range_scored), max(st$k_range_scored), x$path$k0, st$null,
              st$B - x$null$n_fail,
              if (x$null$n_fail) sprintf(" (%d failed)", x$null$n_fail) else "",
              if (is.null(st$coarse)) sprintf("ensemble rpp = %d", st$regions_per_partition)
              else sprintf("%d superpixels", st$coarse)))
  r <- x$rules
  cat("  rules:", paste(sprintf("%s = %s", names(r), ifelse(is.na(r), "NA", r)),
                        collapse = " | "), "\n")
  s <- x$support[, intersect(c("k", "gap", "gap_se", "gap_support", "merge_p",
                               "composite", "composite_support", "vi_next"),
                             names(x$support))]
  s[] <- lapply(s, function(v) if (is.numeric(v)) signif(v, 3) else v)
  print(s, row.names = FALSE)
  cat("  A suggestion, not an estimate: k may be weakly identified. Fit the\n",
      " chosen k with get_regions(final_regions = ).\n")
  invisible(x)
}
