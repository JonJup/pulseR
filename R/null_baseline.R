#' Permutation null baseline for the Level 1 validity criteria
#'
#' Shuffles the rows of \code{env} among units, recomputes every edge weight
#' as the Euclidean distance between the shuffled rows of its two end units,
#' and reruns the Level 1 pipeline (spanning trees, ensemble, consensus,
#' validity criteria) at one ensemble configuration and every requested
#' \code{final_regions}. Contiguity is kept and environmental structure is
#' removed, so each replicate shows what the criteria report when only the
#' map drives the result: the "contiguity-only" floor.
#'
#' @details Edge weights are recomputed as Euclidean distance on \code{env}
#'   plus the bridge penalty stored on \code{graph} (\code{bridge},
#'   \code{bridge_dist}, \code{bridge_scale}), the same function
#'   \code{\link{add_edge_weight}} uses. This tests
#'   environmental signal only if \code{env} is the matrix the observed
#'   weights were computed from, so the function first recomputes the weights
#'   from the unshuffled \code{env} and stops unless they correlate with
#'   \code{E(graph)$weight} above 0.999. The default \code{env = NULL} uses
#'   the matrix \code{add_edge_weight()} stored on \code{graph}, which always
#'   passes. Shuffling leaves column means and standard deviations unchanged,
#'   so no rescaling is needed. A prior
#'   typology is not used in the null. With \code{B} replicates the smallest
#'   attainable p-value is \eqn{1 / (B + 1)}; \code{B = 19} is the minimum for
#'   a 5\% test.
#'
#'   The shuffle also removes spatial autocorrelation in \code{env}. Smooth
#'   but substantively meaningless gradients therefore count as signal here.
#'
#'   The tree noise multiplier \eqn{\eta} is calibrated to
#'   \code{target_edge_turnover} on the first replicate and reused for the
#'   others. Replicates are exchangeable permutations of the same matrix, so
#'   one calibration serves all of them; the noise scale itself is still each
#'   replicate's own weight MAD. The value is returned in \code{settings$eta}.
#'
#' @section Failed replicates:
#'   A replicate fails at a given k only when its consensus raises a
#'   configuration error (class \code{pulseR_config_error}: the null ensemble
#'   cannot support k regions). Any other error stops the call, with the
#'   replicate and k added to the message, because it signals a bug or bad
#'   input rather than a property of the null.
#'
#'   Failures are reported per k as \code{n_fail} and \code{fail_share}, and
#'   enter two p-values:
#'   \describe{
#'     \item{\code{p}}{Failures counted as at least as good as observed:
#'       \eqn{(1 + n_{\ge} + n_{fail}) / (n_{null} + n_{fail} + 1)}. Never
#'       smaller than the p-value with every failure resolved; use this one.}
#'     \item{\code{p_drop}}{Failures dropped:
#'       \eqn{(1 + n_{\ge}) / (n_{null} + 1)}, the test conditional on the
#'       null fitting at k.}
#'   }
#'   The two agree when no replicate fails. If the share of failed replicates
#'   at any k is certain to exceed \code{max_fail_share}, the call stops at
#'   that replicate instead of finishing: lower that k, raise
#'   \code{regions_per_partition}, or drop k from \code{final_regions}.
#'
#' @param graph,env,env_id As in \code{\link{tune_regions}}. \code{env}
#'   must reproduce the edge weights (see Details); \code{NULL} (default)
#'   uses the matrix stored on \code{graph}.
#' @param final_regions Integer vector of k values; all are evaluated on each
#'   null ensemble. Each k adds one consensus fit per replicate, and
#'   \eqn{2 \times} \code{stability_B} more if it is in \code{stability_k}.
#' @param regions_per_partition,N_tree,weight_exponent Single values: the
#'   ensemble configuration under test (e.g. \code{tune_regions()$best_config}).
#' @param size_influence,target_edge_turnover,cut_rule,refine_boundary,min_size
#'   As in \code{\link{get_regions}}; use the values of the observed run.
#' @param B Integer. Null replicates. Default 19.
#' @param stability_B Integer. Random splits for \code{STAB}; use the value of
#'   the observed run.
#' @param stability_k Integer vector: the k values at which \code{STAB} is
#'   computed for the null. STAB costs \eqn{2 \times} \code{stability_B}
#'   consensus fits per k and replicate, most of the run time. \code{NULL}
#'   (default): the k values with an observed \code{STAB} in \code{observed}
#'   (for a \code{"pulse_regions"} object, its k), or all
#'   \code{final_regions} when \code{observed} is \code{NULL}. At other k,
#'   \code{STAB} is \code{NA} in \code{null} and its summary row is
#'   \code{NA}.
#' @param observed Optional \code{"pulse_regions"} or
#'   \code{"pulse_regions_tuning"} object, or a data frame
#'   with \code{final_regions} and validity columns (e.g. a \code{tuning_log}).
#'   Rows are matched on the ensemble configuration when those columns exist.
#'   A \code{"pulse_regions"} object is rescored on \code{graph}, so its
#'   \code{STAB} halves get the same contiguity repair as the null replicate.
#' @param max_fail_share Number in \eqn{[0, 1]}. Largest share of the
#'   \code{B} replicates that may fail at any one k (see Failed replicates).
#'   Default 0.2.
#' @param large_n_threshold,sample_size,n_cores As in \code{\link{tune_regions}}.
#' @param seed Optional integer seed; the caller's RNG state is restored.
#' @param verbose Logical.
#'
#' @return An object of class \code{"pulse_null"}: \code{null} (one row per
#'   replicate and k, with \code{failed} and the configuration-error message
#'   in \code{error}), \code{summary} (one row per k and criterion: observed,
#'   null mean, null 95\% bound in the "better" direction, \code{excess} =
#'   observed minus null mean with the sign set so positive is better than the
#'   null, the empirical \code{p} and \code{p_drop}, \code{n_null},
#'   \code{n_fail}, \code{fail_share}), \code{config}, and \code{B}.
#'   \code{settings} holds \code{cut_rule}, \code{refine_boundary},
#'   \code{min_size}, \code{stability_k}, \code{max_fail_share} and the
#'   reused \code{eta}.
#' @export
compute_null_baseline <- function(
    graph,
    env = NULL,
    env_id = NULL,
    final_regions,
    regions_per_partition = NULL,
    N_tree,
    weight_exponent = 4,
    size_influence = 0,
    target_edge_turnover = 0.1,
    cut_rule = c("ssd", "propensity"),
    refine_boundary = TRUE,
    min_size = 5,
    B = 19,
    stability_B = 10,
    stability_k = NULL,
    observed = NULL,
    max_fail_share = 0.2,
    large_n_threshold = 20000,
    sample_size = 5000,
    n_cores = 1,
    seed = NULL,
    verbose = TRUE) {

  stopifnot(
    inherits(graph, "igraph"),
    length(regions_per_partition) <= 1L,
    length(N_tree) == 1L,
    length(weight_exponent) == 1L, B >= 1
  )
  if (!is.numeric(max_fail_share) || length(max_fail_share) != 1L ||
      is.na(max_fail_share) || max_fail_share < 0 || max_fail_share > 1)
    stop("`max_fail_share` must be one number in [0, 1].", call. = FALSE)
  cut_rule <- match.arg(cut_rule)
  env <- .resolve_env(env, env_id, graph, verbose = FALSE)
  if (is.null(env))
    stop("`env` is required: supply the matrix the edge weights were computed ",
         "from, or build `graph` with add_edge_weight(), which stores it.",
         call. = FALSE)
  if (anyNA(env))
    stop("`env` has missing values; the null needs complete rows.", call. = FALSE)
  final_regions <- sort(unique(as.integer(final_regions)))
  cvi_cols <- names(.CVI_DIRECTION)

  # The null is only a test of environmental signal if it uses the same
  # weight function on the same matrix as the observed graph, including the
  # distance penalty on bridge edges (bridge flag, length and scale are
  # geometry, so they are kept fixed across replicates).
  .null_weights <- function(X) .edge_weights_from_env(graph, X)
  r <- suppressWarnings(stats::cor(.null_weights(env), igraph::E(graph)$weight))
  if (!is.finite(r) || r < 0.999)
    stop(sprintf(paste0(
      "Euclidean distances on `env` do not reproduce E(graph)$weight ",
      "(correlation %.4f, need > 0.999), so the null would not test ",
      "environmental signal. Pass the matrix the weights were computed from ",
      "(scaled if add_edge_weight(scale = TRUE)), or env = NULL to use the one ",
      "stored on `graph`. Weights from another model (e.g. GDM) are not ",
      "Euclidean distances on any stored matrix, so this null does not apply ",
      "to them; see ?set_env."), r), call. = FALSE)
  rpp_tied <- is.null(regions_per_partition)
  rpp_of   <- function(k) if (rpp_tied) 3 * k else regions_per_partition
  cfg <-
    c(
      if (!rpp_tied)
        list(regions_per_partition = regions_per_partition),
      list(
        N_tree = N_tree,
        weight_exponent = weight_exponent
      ),
      if (rpp_tied) list(rpp_per_region = 3)
    )

  stability_k <- if (is.null(stability_k)) .null_stab_k(observed, cfg, final_regions)
  else intersect(final_regions, as.integer(stability_k))
  .local_seed(seed)

  na_cvis <- as.data.frame(as.list(stats::setNames(rep(NA_real_, length(cvi_cols)),
                                                   cvi_cols)))
  # Failed replicates per k so far. The call stops as soon as one k is
  # certain to exceed max_fail_share, instead of finishing all B replicates.
  n_fail_k <- stats::setNames(integer(length(final_regions)), final_regions)

  eta_null <- NULL
  out <- vector("list", B)
  for (b in seq_len(B)) {
    if (verbose) message(sprintf("Null replicate %d / %d ...", b, B))
    env_b <- env[sample.int(nrow(env)), , drop = FALSE]
    g_b   <- graph
    igraph::E(g_b)$weight <- .null_weights(env_b)
    ## The ssd rule and the boundary refinement read the environment
    ## from the graph, so it must be permuted there too.
    g_b   <- igraph::set_graph_attr(g_b, "env", env_b)
    tr  <- sample_spanning_trees(g_b, n = N_tree, verbose = FALSE,
                                 target_edge_turnover = target_edge_turnover,
                                 eta = eta_null)
    if (is.null(eta_null)) eta_null <- tr$eta
    build <- function(rpp) compute_ensemble(tr,
                                            N_tree                = N_tree,
                                            graph                 = g_b,
                                            regions_per_partition = rpp,
                                            weight_exponent       = weight_exponent,
                                            size_influence        = size_influence,
                                            cut_rule              = cut_rule,
                                            min_size              = min_size,
                                            verbose               = FALSE
    )
    ens_fixed <-
      if (rpp_tied) NULL else build(rpp_of(NA))
    stab_seed <- sample.int(.Machine$integer.max, 1L)
    out[[b]] <- do.call(rbind, lapply(final_regions, function(k) {
      ens <- if (rpp_tied) build(rpp_of(k)) else ens_fixed
      # Only a configuration error (the null ensemble cannot support k
      # regions) is a failed replicate. Anything else is a bug or bad
      # input and must not count in favour of the observed data.
      res <- tryCatch(
        compute_consensus(ens, k, graph = g_b,
                          refine_boundary = refine_boundary,
                          min_size = min_size,
                          large_n_threshold = large_n_threshold,
                          sample_size = sample_size, n_cores = n_cores,
                          report_native = FALSE, verbose = FALSE),
        pulseR_config_error = function(e) e,
        error = function(e)
          stop(sprintf("Null replicate %d, k = %d: %s", b, k,
                       conditionMessage(e)), call. = FALSE))
      failed <- inherits(res, "pulseR_config_error")
      cv <- if (failed) na_cvis
      else .consensus_cvis(res, ens$memberships, env = env_b,
                           stability_B = if (k %in% stability_k) stability_B else 0L,
                           graph = g_b,
                           stability_seed = stab_seed,
                           large_n_threshold = large_n_threshold,
                           sample_size = sample_size, n_cores = n_cores)
      cbind(data.frame(rep = b, final_regions = k, failed = failed,
                       error = if (failed) conditionMessage(res) else NA_character_,
                       stringsAsFactors = FALSE),
            cv[, cvi_cols, drop = FALSE])
    }))
    n_fail_k <- n_fail_k + as.integer(tapply(out[[b]]$failed,
                                             factor(out[[b]]$final_regions,
                                                    levels = final_regions), sum))
    over <- n_fail_k / B > max_fail_share
    if (any(over)) {
      k_bad <- final_regions[over][1L]
      msg   <- out[[b]]$error[out[[b]]$final_regions == k_bad & out[[b]]$failed]
      stop(sprintf(paste0(
        "%d of %d null replicates failed at k = %d after replicate %d, more ",
        "than max_fail_share = %g. A p-value from the remaining replicates ",
        "would be unreliable. Lower k, raise regions_per_partition, or drop ",
        "k = %d from final_regions. Last error: %s"),
        n_fail_k[over][1L], B, k_bad, b, max_fail_share, k_bad,
        if (length(msg)) msg[1L] else "(see earlier replicate)"), call. = FALSE)
    }
  }
  null <- do.call(rbind, out)
  if (any(n_fail_k > 0L))
    warning(sprintf(paste0(
      "Null replicates failed with a configuration error at %s. ",
      "See summary$fail_share; `p` counts failures as at least as good ",
      "as observed, `p_drop` excludes them."),
      paste(sprintf("k = %s (%d of %d)", names(n_fail_k)[n_fail_k > 0L],
                    n_fail_k[n_fail_k > 0L], B), collapse = ", ")),
      call. = FALSE)

  if (inherits(observed, "pulse_regions")) {
    cx  <- observed$config
    ax <- setdiff(names(cfg), "rpp_per_region")
    off <- ax[
      !vapply(
        ax,
        function(a)
          isTRUE(abs(cx[[a]] - cfg[[a]]) < 1e-9),
        logical(1)
      )
    ]
    if (rpp_tied && !isTRUE(cx$regions_per_partition == rpp_of(cx$final_regions)))
      off <- c(off, "regions_per_partition")
    if (length(off))
      warning("`observed` was fitted with different ", paste(off, collapse = ", "),
              " than the null; observed and null are not comparable.", call. = FALSE)
    if (!observed$k %in% final_regions)
      warning(sprintf("`observed` has k = %d, which is not in `final_regions`; ",
                      observed$k), "`observed` will be NA.", call. = FALSE)
    
    observed <- .validity_row(
      observed, 
      env, 
      graph = graph, 
      stability_B = stability_B, 
      seed = seed,
      large_n_threshold = large_n_threshold,
      sample_size = sample_size, 
      n_cores = n_cores
      )
  }


  structure(list(
    null = null,
    summary = .null_summary(
      null,
      observed,
      cfg,
      skip = list(
        STAB = setdiff(
          final_regions,
          stability_k
        )
      )
    ),
    config = cfg,
    B = as.integer(B),
    settings = list(
      cut_rule = cut_rule,
      refine_boundary = refine_boundary,
      min_size = min_size,
      stability_k = stability_k,
      max_fail_share = max_fail_share,
      eta = eta_null)
  ),
  class = "pulse_null"
  )
}

#' @keywords internal
#' @noRd
.null_summary <- function(null, observed, cfg, skip = list()) {
  observed <- .match_observed_config(observed, cfg)
  if (is.null(null$failed)) null$failed <- FALSE
  rows <- list()
  for (k in sort(unique(null$final_regions))) for (cr in names(.CVI_DIRECTION)) {
    at_k   <- null$final_regions == k
    n_fail <- sum(null$failed[at_k])
    x   <- null[[cr]][at_k & !null$failed]; x <- x[is.finite(x)]
    sgn <- if (.CVI_DIRECTION[[cr]] == "max") 1 else -1
    obs <- if (is.null(observed)) NA_real_ else {
      o <- observed[[cr]][observed$final_regions == k]
      if (length(o)) o[1] else NA_real_
    }
    if (k %in% skip[[cr]]) { x <- numeric(0); n_fail <- 0L }   # not computed for the null
    n_ge <- sum(sgn * x >= sgn * obs)
    rows[[length(rows) + 1L]] <- data.frame(
      final_regions = k, criterion = cr, observed = obs,
      null_mean = if (length(x)) mean(x) else NA_real_,
      null_95   = if (length(x)) unname(stats::quantile(x, if (sgn > 0) 0.95 else 0.05)) else NA_real_,
      excess    = if (length(x)) sgn * (obs - mean(x)) else NA_real_,
      # Failed replicates count as at least as good as observed: p can only
      # rise, never fall, compared with the p every failure resolved would give.
      p          = if (length(x) + n_fail > 0L && is.finite(obs))
        (1 + n_ge + n_fail) / (length(x) + n_fail + 1) else NA_real_,
      p_drop     = if (length(x) > 0L && is.finite(obs))
        (1 + n_ge) / (length(x) + 1) else NA_real_,
      n_null     = length(x),
      n_fail     = n_fail,
      fail_share = n_fail / max(sum(at_k), 1L))
  }
  do.call(rbind, rows)
}

# Observed validity rows for the null configuration. Accepts a tuning object
# or a data frame; rows are kept where the configuration columns match `cfg`.
#' @keywords internal
#' @noRd
.match_observed_config <- function(observed, cfg) {
  if (inherits(observed, "pulse_regions_tuning")) observed <- observed$tuning_log
  if (is.null(observed)) return(NULL)
  observed <- as.data.frame(observed)
  for (ax in setdiff(names(cfg), "rpp_per_region"))
    if (ax %in% names(observed))
      observed <- observed[abs(observed[[ax]] - cfg[[ax]]) < 1e-9, , drop = FALSE]
  if (!is.null(cfg$rpp_per_region) && "regions_per_partition" %in% names(observed))
    observed <- observed[observed$regions_per_partition ==
                           cfg$rpp_per_region * observed$final_regions, , drop = FALSE]
  observed
}

# k values at which the null needs STAB: those with an observed STAB, or all
# final_regions when nothing is observed.
#' @keywords internal
#' @noRd
.null_stab_k <- function(observed, cfg, final_regions) {
  if (is.null(observed)) return(final_regions)
  if (inherits(observed, "pulse_regions")) return(intersect(final_regions, observed$k))
  observed <- .match_observed_config(observed, cfg)
  if (!"STAB" %in% names(observed)) return(integer(0))
  intersect(final_regions, observed$final_regions[is.finite(observed$STAB)])
}

#' @rdname compute_null_baseline
#' @param x An object of class \code{"pulse_null"}.
#' @param ... Unused.
#' @return \code{print.pulse_null} returns \code{x} invisibly.
#' @export
print.pulse_null <- function(x, ...) {
  rpp <- if (is.null(x$config$rpp_per_region)){
    format(x$config$regions_per_partition)
  } else {
    paste0(x$config$rpp_per_region," x k")
  }
  cat(
    sprintf("<pulse_null>  %d replicates | regions_per_partition = %s | N_tree = %d | weight_exponent = %g\n",
            x$B,
            rpp,
            x$config$N_tree,
            x$config$weight_exponent
    )
  )
  if (!is.null(x$settings))
    cat(sprintf("  cut_rule = %s | refine_boundary = %s | min_size = %s\n",
                x$settings$cut_rule, x$settings$refine_boundary, x$settings$min_size))
  fs <- x$summary[x$summary$criterion == "CH_FIT", c("final_regions", "n_fail", "fail_share")]
  if (!is.null(fs$n_fail) && any(fs$n_fail > 0L, na.rm = TRUE)) {
    bad <- fs[fs$n_fail > 0L, , drop = FALSE]
    cat(sprintf("  failed replicates: %s (p counts them as >= observed; see p_drop)\n",
                paste(sprintf("k = %d: %d (%.0f%%)", bad$final_regions, bad$n_fail,
                              100 * bad$fail_share), collapse = ", ")))
  } else {
    cat("  failed replicates: none\n")
  }
  s <- x$summary[x$summary$criterion %in% c("CH_FIT", "STAB"), ]
  s[] <- lapply(s, function(v) if (is.numeric(v)) signif(v, 3) else v)
  print(s, row.names = FALSE)
  invisible(x)
}
