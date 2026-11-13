#' Delineate fuzzy contiguous regions (Level 1, one-call wrapper)
#'
#' Samples spanning trees, cuts them into an ensemble of contiguous
#' partitions, and summarises the ensemble as \code{final_regions} consensus
#' regions with composition memberships \eqn{\Pi}. If any of
#' \code{regions_per_partition}, \code{weight_exponent}, \code{N_tree} or
#' \code{final_regions} has more than one value, the configuration is chosen
#' by \code{\link{tune_regions}}; otherwise one ensemble and one consensus are
#' computed.
#'
#' @param graph A weighted \code{igraph} object, usually from
#'   \code{\link{add_edge_weight}}. Edges must carry a numeric \code{weight}
#'   attribute (finite, non-negative).
#' @param env Optional data frame or matrix, one row per vertex of
#'   \code{graph} in vertex order, for the validity indices \code{CH_FIT} and
#'   \code{R2_FIT} during tuning. \code{NULL} (default) uses the matrix
#'   \code{\link{add_edge_weight}} stored on \code{graph}, so the indices are
#'   computed in the space of the edge weights. A supplied \code{env} is used
#'   as given: non-numeric columns are dropped and nothing is scaled.
#' @param env_id Character vector of columns in \code{env} to exclude from the
#'   validity indices, for example numeric ID columns. Default \code{NULL}.
#' @param regions_per_partition Integer, or integer vector of candidates.
#'   Number of contiguous pieces each spanning tree is cut into. Default 50.
#' @param final_regions Integer, or integer vector of candidates. Number of
#'   consensus regions \eqn{k}.
#' @param N_tree Integer, or integer vector of candidates. Number of spanning
#'   trees, one partition each. Trees are sampled once at \code{max(N_tree)}.
#' @param weight_exponent Non-negative numeric scalar \eqn{\gamma}, an inverse
#'   temperature. It multiplies the whole log cut propensity. The size term and
#'   the prior shift are both expressed in \eqn{u = \mathrm{sd}(\gamma \log w)},
#'   so for \eqn{\gamma > 0}
#'   \deqn{\log \pi_e = \gamma \left[\log \tilde w_e
#'         + \texttt{size\_influence} \tfrac{\sigma}{\sigma_s} \log s_e
#'         - \lambda \sigma \, 1[e \in \text{prior}]\right],
#'         \quad \sigma = \mathrm{sd}(\log \tilde w).}
#'   \code{size_influence} and \code{prior_strength} decide which edges are
#'   preferred. \eqn{\gamma} decides how strongly that preference is followed.
#'   \eqn{\gamma \to 0} approaches uniform cuts, and \eqn{\gamma \to \infty}
#'   approaches deterministic removal of the top-scoring edges. Those are the
#'   heaviest edges only when \code{size_influence = 0} and
#'   \code{prior_strength = 0}. Raising \eqn{\gamma} does not shift emphasis
#'   from size to environment.
#'   At \eqn{\gamma = 0} cuts are uniform and \code{size_influence} is ignored.
#'   \code{prior_strength} falls back to one halving per unit.
#' @param size_influence Non-negative numeric. Spread of the size-balance term
#'   \eqn{\log(n_c (N - n_c) / N)} relative to the spread of the environmental
#'   term \eqn{\gamma \log w}, per tree. The default (0) ignores piece size.
#' @param prior_typology Optional. Vertex attribute name (single character)
#'   or vector of length \code{vcount(graph)} giving a prior ecoregion label
#'   at each vertex. \code{NA} means no prior at that vertex. Default
#'   \code{NULL}.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda}. Lowers the
#'   log cutting propensity of within-prior edges by \eqn{\lambda} standard
#'   deviations of the environmental term \eqn{\gamma \log w}, computed per
#'   tree: the same unit as \code{size_influence}. \eqn{\lambda = 0} (default)
#'   ignores the prior; \code{Inf} returns the prior itself as one-hot
#'   memberships. If that spread is zero (\code{weight_exponent = 0} or
#'   constant weights), one unit halves the propensity instead, with a
#'   warning.
#' @param cut_rule \code{"propensity"} or \code{"ssd"} (default): how each
#'   tree is cut. See \code{\link{compute_ensemble}}.
#' @param refine_boundary Logical. Contiguity-preserving boundary refinement
#'   after consensus. Default \code{TRUE}. See \code{\link{compute_consensus}}.
#' @param min_size Integer \eqn{\ge 1}. Smallest piece the \code{"ssd"} rule may
#'   cut off, and smallest region the boundary refinement may leave. Default 1.
#' @param target_edge_turnover Numeric in \eqn{[0, 1)}. Average share of the
#'   unperturbed minimum spanning tree's edges replaced in each sampled tree.
#'   Default 0.1.
#' @param tuning_strategy Character. One of \code{"grid"} (default),
#'   \code{"sequential"}, \code{"iterative"}. See \code{\link{tune_regions}}.
#' @param tuning_metric Character. One of \code{"borda"} (default),
#'   \code{"ch_fit"}, \code{"r2_fit"}, \code{"stability"},
#'   \code{"silhouette_ca"}, \code{"expected_vi"}. See
#'   \code{\link{tune_regions}}.
#' @param stability_B Integer. Random splits for the \code{STAB} index
#'   (split-half ensemble agreement, see \code{\link{tune_regions}}) during
#'   tuning. 0 skips it. Default 10.
#' @param large_n_threshold Integer. Above this number of vertices, salso
#'   runs on a sample of units and the rest are assigned and refined. Default
#'   20000. See \code{\link{compute_consensus}}.
#' @param sample_size As in \code{\link{tune_regions}}.
#' @param n_cores Integer. Cores passed to \code{salso::salso()}. Default 1.
#' @param trees Optional \code{"pulse_trees"} object from
#'   \code{\link{sample_spanning_trees}}. If supplied, no trees are sampled
#'   and \code{prior_typology} and \code{target_edge_turnover} are ignored.
#' @param max_iter Integer. Maximum passes for the \code{"iterative"}
#'   strategy. Default 10.
#' @param verbose Logical. Default \code{TRUE}.
#' @param seed Optional integer seed. One stream drives the whole call; the
#'   caller's RNG state is restored on exit.
#'
#' @return An object of class \code{"pulse_regions"}, a list with
#'   \describe{
#'     \item{\code{memberships}}{\eqn{N \times k} matrix \eqn{\Pi}. Rows sum
#'       to 1; filled from graph neighbours (see compute_consensus() Details);
#'        all-zero for uncovered units when
#'       \code{prior_strength = Inf}.}
#'     \item{\code{entropy}}{Row entropy of \code{memberships}.}
#'     \item{\code{hard_labels}}{Integer vector of consensus labels.}
#'     \item{\code{k}}{Number of regions.}
#'     \item{\code{config}}{One-row data frame: \code{regions_per_partition},
#'       \code{weight_exponent}, \code{N_tree}, \code{final_regions},
#'       \code{cut_rule}, \code{refine_boundary}, \code{min_size}.}
#'     \item{\code{tuned}}{Logical. \code{TRUE} if \code{tune_regions()} ran.}
#'     \item{\code{score}, \code{tuning_metric}, \code{tuning_strategy}}{Tuning result;
#'       \code{NA} when not tuned.}
#'     \item{\code{tuning_log}, \code{candidates},
#'       \code{n_ensembles_built}}{Search provenance; \code{NULL} or \code{NA}
#'       when not tuned.}
#'     \item{\code{prior_strength}}{\eqn{\lambda} as applied.}
#'     \item{\code{coverage}}{Logical per unit when
#'       \code{prior_strength = Inf}, else \code{NULL}.}
#'     \item{\code{fit}}{The \code{"pulse_consensus"} object from
#'       \code{\link{compute_consensus}}, with certainty, margin, cohesion,
#'       loss terms and diagnostics; \code{NULL} when
#'       \code{prior_strength = Inf}.}
#'     \item{\code{trees}}{The \code{"pulse_trees"} object used.}
#'     \item{\code{ensemble}}{The \eqn{N \times N_{tree}} label matrix.}
#'     \item{\code{p0}, \code{n_eff}}{Mean chance co-association of the
#'       ensemble and \eqn{1 / p_0}.}
#'   }
#'
#' @seealso \code{\link{sample_spanning_trees}},
#'   \code{\link{compute_ensemble}},
#'   \code{\link{compute_consensus}}, \code{\link{tune_regions}}.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_lattice(c(15, 15))
#' E(g)$weight <- runif(ecount(g))
#' res <- get_regions(g, regions_per_partition = c(10, 20), final_regions = c(3, 4),
#'                    N_tree = 30, seed = 1, cut_rule = "propensity", refine_boundary=FALSE)
#' table(res$hard_labels)
#' }
#'
#' @importFrom igraph vcount
#' @export
get_regions <- function(graph,
                        env = NULL,
                        env_id = NULL,
                        regions_per_partition = NULL,
                        final_regions,
                        N_tree,
                        weight_exponent     = 4,
                        size_influence      = 0,
                        cut_rule            = c("ssd", "propensity"),
                        refine_boundary     = TRUE,
                        min_size            = 5,
                        prior_typology = NULL,
                        prior_strength = 0,
                        target_edge_turnover = 0.1,
                        tuning_strategy = c("grid", "sequential", "iterative"),
                        tuning_metric   =c("borda", "ch_fit", "r2_fit", "stability", "silhouette_ca", "expected_vi"),
                        stability_B      = 10,
                        large_n_threshold = 20000,
                        sample_size = 5000,
                        n_cores = 1,
                        trees = NULL,
                        max_iter = 10,
                        verbose  = TRUE,
                        seed     = NULL) {
  
  
  # match args --------------------------------------------------------------*
  tuning_strategy        <- match.arg(tuning_strategy)
  tuning_metric          <- .renamed_cvi(tuning_metric, "tuning_metric")
  tuning_metric          <- match.arg(tuning_metric)
  cut_rule               <- match.arg(cut_rule)
  min_size               <- as.integer(min_size)
  
  # catch errors --------------------------------------------------------------*
  stopifnot(
    inherits(graph, "igraph"), 
    is.numeric(prior_strength) | is.infinite(prior_strength), 
    length(prior_strength) == 1L,
    prior_strength >= 0,
    is.numeric(regions_per_partition)|| is.null(regions_per_partition),
    is.numeric(N_tree),
    is.numeric(final_regions)
    
  )
  # One RNG stream for the whole call; the caller's state is restored on
  # exit. Downstream functions get seed = NULL so they do not reset it.
  .local_seed(seed)
  
  # Upper-bound check before doing any sampling work.
  if (any(regions_per_partition >= igraph::vcount(graph))) {
    stop("`regions_per_partition` must be smaller than the number of vertices (",
         igraph::vcount(graph), ").", call. = FALSE)
  }
  
  
  # --- Exact prior reproduction short-circuit ----------------------------*
  # Triggered by prior_strength = INF
  if (is.infinite(prior_strength)) {
    # Catch error 
    if (is.null(prior_typology)) {
      stop("`prior_strength = Inf` requires `prior_typology`.", call. = FALSE)
    }
    
    pr <- .exact_prior_result(.resolve_prior_typology(graph, prior_typology), verbose = verbose)
    K <- ncol(pr$memberships)
    return(.new_pulse_regions(
      fit           = NULL,
      memberships   = pr$memberships,
      hard_labels = pr$hard_labels,
      coverage      = pr$coverage,
      ensemble      = NULL,
      config        = list(regions_per_partition = K,
                           weight_exponent      = NA_real_,
                           N_tree                 = NA_integer_,
                           final_regions        = K,
                           cut_rule             = NA_character_,
                           refine_boundary      = NA,
                           min_size             = NA_integer_),
      prior_strength = Inf, tuned = FALSE))
  }
  
  # --- Exact tuneable parameters  ----------------------------*
  
  N_tree_values         <- sort(unique(as.integer(N_tree)))
  final_regions         <- sort(unique(as.integer(final_regions)))
  rpp_tied              <- is.null(regions_per_partition)
  regions_per_partition <- if(rpp_tied){
    if (length(final_regions) == 1) 3 * final_regions else NULL 
  } else {
    sort(unique(as.integer(regions_per_partition)))
  }
  
  if (any(regions_per_partition >= igraph::vcount(graph))){
    stop(
      "`regions_per_partition` (",
      paste0(regions_per_partition),
      ") must be smaller than the number of units (",
      paste(igraph::vcount(graph)),
      ").",
      call. = FALSE
    )
  }
  if (!rpp_tied && length(regions_per_partition) == 1 && regions_per_partition < max(final_regions)){
    stop(
      sprintf(
        paste0(
          "`regions_per_partition` (%d) must be as large as the largest",
          "`final_regions` (%d)." 
        ),
        regions_per_partition,
        max(final_regions)
      ),
      call. = FALSE
    )
  }
  
  # --- lambda_E > 0 needs prior labels; fail before any sampling ------*
  # Labels come from `trees` if supplied, otherwise from `prior_typology`.
  has_prior <- if (is.null(trees)) !is.null(prior_typology) else !all(is.na(trees$prior_labels))
  if (prior_strength > 0 && !has_prior) {
    stop("`prior_strength` > 0 but no prior labels found: supply ",
         "`prior_typology`, or `trees` sampled with one.", call. = FALSE)
  }
  
  
  
  # --- Sample the tree ensemble at the largest requested N_tree 
  if (is.null(trees))
  {
    # Nothing random runs before this, so these trees are identical to
    # sample_spanning_trees(..., seed = seed) called on its own.
    trees <- sample_spanning_trees(graph  = graph, 
                                   n               = max(N_tree_values),
                                   seed            = NULL, 
                                   verbose         = verbose,
                                   prior_typology  = prior_typology, 
                                   target_edge_turnover = target_edge_turnover,
                                   prior_strength       = prior_strength 
    )
    trees$seed <- seed
  }
  # Common random numbers for the stability metric across configurations,
  # drawn from the same stream so one `seed` fixes everything.
  stab_seed <- if (is.null(seed)) NULL else sample.int(.Machine$integer.max, 1L)
  
  # Candidate sets are de-duplicated first, so c(50, 50) is one value.
  needs_tuning <-
    length(regions_per_partition) > 1L || length(final_regions) > 1L ||
    length(N_tree_values) > 1L ||
    length(unique(weight_exponent)) > 1L
  
  if (needs_tuning) {
    if (verbose) message("Tuning Level 1 configurations ...")
    tr <- tune_regions(
      graph = graph, 
      env = env, 
      env_id = env_id, 
      large_n_threshold = large_n_threshold,
      n_cores = n_cores, 
      trees = trees, 
      N_tree = N_tree_values,
      regions_per_partition = regions_per_partition,
      final_regions = final_regions, 
      tuning_strategy = tuning_strategy,
      tuning_metric = tuning_metric, 
      stability_B = stability_B,
      prior_strength = prior_strength,
      weight_exponent = weight_exponent,
      max_iter = max_iter, 
      verbose = verbose, 
      seed = NULL,
      stability_seed = stab_seed,
      size_influence = size_influence,
      cut_rule = cut_rule,
      refine_boundary = refine_boundary,
      min_size = min_size,
      sample_size = sample_size
    )
    return(.new_pulse_regions(
      fit               = tr$best_result,
      config            = cbind(tr$best_config,
                                data.frame(cut_rule = cut_rule,
                                           refine_boundary = refine_boundary,
                                           min_size = min_size)),
      trees             = trees,
      prior_strength    = prior_strength,
      tuned             = TRUE,
      score             = tr$best_score,
      tuning_metric     = tr$tuning_metric,
      tuning_strategy   = tr$tuning_strategy,
      tuning_log        = tr$tuning_log,
      candidates        = tr$candidates,
      n_ensembles_built = tr$n_ensembles_built,
      ensemble          = tr$ensemble))
  } else {
    if (verbose) message("Single configuration; no tuning.")
    
    ensemble <- compute_ensemble(
      trees                 = trees, 
      N_tree                = N_tree_values, 
      graph                 = graph,
      regions_per_partition = regions_per_partition,
      weight_exponent       = weight_exponent,
      prior_strength        = prior_strength,
      size_influence        = size_influence,
      cut_rule              = cut_rule,
      min_size              = min_size, verbose = verbose)
    
    fit <- compute_consensus(
      ensemble, 
      final_regions = final_regions, 
      graph = graph,
      large_n_threshold = large_n_threshold,
      refine_boundary = refine_boundary,
      min_size = min_size,
      n_cores = n_cores, 
      verbose = verbose,
      sample_size = sample_size
    )
    .new_pulse_regions(
      fit      = fit,
      ensemble = ensemble$memberships,
      config   = list(regions_per_partition = regions_per_partition,
                      weight_exponent       = weight_exponent,
                      N_tree                = N_tree_values,
                      final_regions         = final_regions,
                      cut_rule              = cut_rule,
                      refine_boundary       = refine_boundary,
                      min_size              = min_size),
      trees          = trees,
      prior_strength = prior_strength,
      tuned          = FALSE)
    
  }
}

# -------------------------------------------------------------------------*
# -------------------------------------------------------------------------*
# -------------------------------------------------------------------------*



#' Sample spanning trees from a Spatial Graph
#'
#' Draws \code{n} spanning trees from a weighted, connected \code{igraph}
#' object; a disconnected graph is an error. The trees, together with the parent graph's edge
#' weights/endpoints and any resolved prior-typology information, are bundled
#' into a \code{"pulse_trees"} object that the rest of the pipeline consumes.
#'
#' @param graph A connected \code{igraph} object. Edges must carry a
#'   numeric \code{weight} attribute (dissimilarity; finite and non-negative).
#'   Any \code{eid} edge attribute is overwritten with the edge position.
#'   Vertices may optionally carry a prior-typology attribute (see
#'   \code{prior_typology}).
#' @param n Integer scalar \eqn{\ge 1}. Number of spanning trees to sample.
#' @param seed Integer scalar or \code{NULL}. If supplied, used to seed the
#'   sampler; the caller's RNG state is restored on exit.
#' @param verbose Logical. Emit progress/diagnostic messages? Default
#'   \code{TRUE}.
#' @param prior_typology Optional. Either (a) a single character giving the
#'   name of a vertex attribute on \code{graph} that holds prior ecoregion
#'   labels (character / factor / integer); or (b) a vector of length
#'   \code{vcount(graph)} of such labels. \code{NA} indicates no prior at
#'   that vertex. When supplied, edges connecting two vertices with the same
#'   non-NA label are flagged as "within-prior"; their cut propensity is
#'   lowered in the cutting step (see \code{prior_strength} in
#'   \code{\link{compute_ensemble}} and
#'   \code{\link{get_regions}}). Default \code{NULL} (no prior). Tree shape
#'   does not depend on the prior.
#' @param target_edge_turnover Numeric in \eqn{[0, 1)}. Average share of the
#'   unperturbed minimum spanning tree's edges replaced in each sampled tree.
#'   The noise scale is calibrated to reach it. 0 gives \code{n} copies of
#'   the minimum spanning tree. Default 0.1.
#' @param prior_strength Numeric. Stored in the result for
#'   \code{\link{compute_ensemble}}; not used for sampling.
#' @param eta \code{NULL} (default) or one non-negative number. The noise
#'   multiplier: noise sd = \code{eta} times the perturbation unit (see
#'   Details). \code{NULL} calibrates it to \code{target_edge_turnover};
#'   a number skips the calibration and \code{target_edge_turnover} is not
#'   used for it. \code{compute_null_baseline()} uses this to calibrate once
#'   for all replicates.
#'
#' @return An object of class \code{"pulse_trees"} — a list with:
#'   \describe{
#'    \item{\code{edge_ids}}{Integer matrix (\code{(N - 1) x n}) of parent-edge
#'       positions; column j holds tree j's edges in \code{igraph::mst()} order.
#'       Weights: \code{edge_weights[edge_ids[, j]]}; endpoints:
#'       \code{edge_endpoints[edge_ids[, j], ]}.}
#'     \item{\code{n_trees}}{Number of trees sampled.}
#'     \item{\code{n_vertices}}{Number of vertices in the original graph.}
#'     \item{\code{edge_weights}}{Numeric vector of all parent-graph edge
#'       weights, ordered as \code{E(graph)}.}
#'     \item{\code{edge_endpoints}}{Integer matrix (\code{|E| x 2}) of edge
#'       endpoint vertex indices, ordered as \code{E(graph)}.}
#'     \item{\code{seed}}{The seed used (or \code{NULL}).}
#'     \item{\code{prior_labels}}{Resolved prior labels (character vector
#'       length N), or \code{NULL}.}
#'     \item{\code{within_prior_edge}}{Logical vector length |E| flagging
#'       within-prior edges, or \code{NULL}.}
#'     \item{\code{prior_strength}}{As supplied.}
#'     \item{\code{noise_scale}, \code{eta}}{Base unit of the weight
#'       perturbation (MAD of the edge weights, see Details) and the calibrated
#'       multiplier; the noise standard deviation is their product.}
#'   }
#'
#' @details The perturbation unit is the median absolute deviation (MAD) of
#'   the edge weights. When more than half of the weights are equal, the MAD
#'   is 0 and would remove all perturbation; the unit then falls back to
#'   \code{IQR / 1.349}, then the standard deviation, then 1 (all weights
#'   equal), with a warning.
#'
#' @seealso \code{\link{get_regions}} for the full one-call pipeline.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- sample_gnp(50, 0.2)
#' E(g)$weight <- runif(ecount(g))
#' trees <- sample_spanning_trees(g, n = 20, seed = 1)
#' print(trees)
#' }
#'
#' @importFrom igraph E ends vcount vertex_attr vertex_attr_names edge_attr components
#' @importFrom stats mad rnorm
#' @export
sample_spanning_trees <- function(graph, n, seed = NULL, verbose = TRUE,
                                  prior_typology = NULL, target_edge_turnover = .1, 
                                  prior_strength = 0, eta = NULL) {
  
  # --- Argument validation ---- *--- *
  stopifnot(
    inherits(graph, "igraph"),
    is.numeric(n), length(n) == 1L, is.finite(n), n >= 1
  )
  
  
  # Edges must have usable weights; fail early with a clear message rather than
  # producing NA-laden cut probabilities deep inside the sampler.
  if (is.null(igraph::edge_attr(graph, "weight"))) {
    stop("`graph` edges must carry a numeric `weight` attribute.", call. = FALSE)
  }
  edge_w <- igraph::E(graph)$weight
  if (anyNA(edge_w) || any(!is.finite(edge_w)) || any(edge_w < 0)) {
    stop("Edge `weight`s must be finite and non-negative.", call. = FALSE)
  }
  
  # Trees are mapped back to parent edges through `eid`, which must be the
  # edge position. Set it here so any weighted igraph works, and so a
  # user-supplied, non-positional `eid` cannot mis-map edges.
  igraph::E(graph)$eid <- seq_len(igraph::ecount(graph))
  
  # Extract properties of weights later used for perturbation
  
  # A disconnected graph yields a spanning *forest* (n_v - n_components edges),
  # which breaks the "cut n_reg-1 edges of a tree" arithmetic downstream.
  comp <- igraph::components(graph)
  if (comp$no > 1L) {
    stop(sprintf(paste0(
      "`graph` has %d connected components (largest %d, smallest %d ",
      "vertices). Each tree would be a forest, every partition would get ",
      "%d extra pieces, and `regions_per_partition` would lose its meaning. ",
      "Connect the graph first, e.g. build_network(connect_islands = TRUE) ",
      "or connect_components()."),
      comp$no, max(comp$csize), min(comp$csize), comp$no - 1L), call. = FALSE)
  }
  
  # --- RNG handling: seed locally, restore the caller's state on exit ---- *-
  .local_seed(seed)
  edge_endpoints    <- igraph::ends(graph, igraph::E(graph), names = FALSE)
  # --- Resolve prior typology (if supplied) ---- *---- *---- *---- *---- *---- *--
  if (!is.null(prior_typology)){
    #prior_typology    <- "prior"
    prior_labels      <- .resolve_prior_typology(graph, 
                                                 prior_typology)
    within_prior_edge <- .compute_within_prior_edge(prior_labels, 
                                                    edge_endpoints)
    
    if (verbose && !is.null(prior_labels)) {
      n_cov    <- sum(!is.na(prior_labels))
      n_within <- sum(within_prior_edge, na.rm = TRUE)
      n_edges  <- length(within_prior_edge)
      message(sprintf(
        "  Prior typology: %d / %d vertices covered (%.1f%%), %d / %d edges within-prior (%.1f%%)",
        n_cov, length(prior_labels), 100 * n_cov / length(prior_labels),
        n_within, n_edges, 100 * n_within / max(n_edges, 1L)
      ))
    }
  } else {
    prior_labels      <- NULL
    within_prior_edge <- NULL
  }
  
  # Prior weighting is applied at the cutting step only
  # (compute_ensemble), never to tree shape.
  
  # --- Sample the trees ---- *---- *---- *---- *---- *---- *---- *---- *---- *---- *--
  if (verbose) message(sprintf("Sampling %d spanning trees ...", n))
  tick <- if (verbose) .make_progress(n, "Spanning trees") else function(i) NULL
  edge_mad <- .noise_scale(edge_w, warn = target_edge_turnover > 0)
  eta_emp <- if (is.null(eta)) {
    .calibrate_eta_empirical(graph  = graph, 
                             w      = edge_w, 
                             s      = edge_mad,
                             target = target_edge_turnover)
  } else {
    if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta < 0)
      stop("`eta` must be NULL or one finite number >= 0.", call. = FALSE)
    as.numeric(eta)
  }
  n_v      <- igraph::vcount(graph)
  edge_ids <- matrix(0L, nrow = n_v - 1L, ncol = n)
  for (j in seq_len(n)) {
    wp <- edge_w + eta_emp * edge_mad * rnorm(length(edge_w))
    # Keep the tree's own edge order. Cut indices and RNG draws in
    # compute_ensemble() are positional; sorting changes results.
    edge_ids[, j] <- igraph::E(igraph::mst(graph, weights = wp))$eid
    tick(j)
  }
  
  structure(
    list(
      edge_ids          = edge_ids,
      n_trees           = as.integer(n),
      n_vertices        = n_v,
      edge_weights      = edge_w,
      edge_endpoints    = edge_endpoints,
      seed              = seed,
      prior_labels      = prior_labels, 
      within_prior_edge = within_prior_edge,
      prior_strength    = prior_strength,
      noise_scale       = edge_mad,
      eta               = eta_emp
    ),
    class = "pulse_trees"
  )
}

# -------------------------------------------------------------------------*
# -------------------------------------------------------------------------*
# -------------------------------------------------------------------------*



#' Compute Ensemble Region Memberships via SKATER Edge Removal
#'
#' For each spanning tree, removes \code{regions_per_partition - 1} edges to
#' produce \code{regions_per_partition} connected pieces, then records each
#' vertex's piece label. Edges are drawn without replacement by an exponential
#' race on the log propensity
#' \eqn{\gamma \log w_e + \alpha \log s_e - \lambda u 1[e \text{ within prior}]},
#' where \eqn{s_e = n_c (N - n_c) / N} is the size balance of the cut.
#'
#' @param trees A \code{"pulse_trees"} object.
#' @param N_tree Integer scalar. Number of spanning trees to use
#'   (must be \code{<= trees$n_trees}).
#' @param graph The parent \code{igraph} object the trees were sampled from.
#'   Read only when \code{cut_rule = "ssd"}, which needs the scaled
#'   environment matrix \code{graph_attr(graph, "env")} stored by
#'   \code{\link{add_edge_weight}}. Tree edges and weights are taken from
#'   \code{trees}; \code{graph} is not checked against \code{trees}.
#' @param regions_per_partition Integer scalar \eqn{\ge 2}. Number of regions per tree;
#'   must not exceed the number of vertices.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda}. Lowers the
#'   log cutting propensity of within-prior edges by \eqn{\lambda} standard
#'   deviations of the environmental term \eqn{\gamma \log w}, computed per
#'   tree: the same unit as \code{size_influence}, so a value means the same
#'   on any graph. \eqn{\lambda = 0} (default) ignores the prior;
#'   \eqn{\lambda = 1} makes a within-prior edge as hard to cut as an edge
#'   one standard deviation lighter on the log scale. If that spread is zero
#'   (\code{weight_exponent = 0} or constant weights), one unit halves the
#'   propensity instead, with a warning.
#' @param verbose Logical. Default \code{TRUE}.
#' @param weight_exponent Non-negative numeric scalar \eqn{\gamma}, an inverse
#'   temperature. It multiplies the whole log cut propensity. The size term and
#'   the prior shift are both expressed in \eqn{u = \mathrm{sd}(\gamma \log w)},
#'   so for \eqn{\gamma > 0}
#'   \deqn{\log \pi_e = \gamma \left[\log \tilde w_e
#'         + \texttt{size\_influence} \tfrac{\sigma}{\sigma_s} \log s_e
#'         - \lambda \sigma \, 1[e \in \text{prior}]\right],
#'         \quad \sigma = \mathrm{sd}(\log \tilde w).}
#'   \code{size_influence} and \code{prior_strength} decide which edges are
#'   preferred. \eqn{\gamma} decides how strongly that preference is followed.
#'   \eqn{\gamma \to 0} approaches uniform cuts, and \eqn{\gamma \to \infty}
#'   approaches deterministic removal of the top-scoring edges. Those are the
#'   heaviest edges only when \code{size_influence = 0} and
#'   \code{prior_strength = 0}. Raising \eqn{\gamma} does not shift emphasis
#'   from size to environment.
#'   At \eqn{\gamma = 0} cuts are uniform and \code{size_influence} is ignored.
#'   \code{prior_strength} falls back to one halving per unit.
#' @param size_influence Non-negative numeric. Spread of the size term
#'   \eqn{\alpha \log s_e} relative to the spread \eqn{u} of the
#'   environmental term \eqn{\gamma \log w}, per tree. Default (0) ignores piece size.
#' @param cut_rule \code{"propensity"}: the
#'   \code{regions_per_partition - 1} cuts of a tree are drawn at once by the
#'   exponential race above. \code{"ssd"} (default): cuts are drawn one at a time, each
#'   by an exponential race on the Ward reduction in within-piece sum of
#'   squares of \code{graph_attr(graph, "env")}, recomputed after every cut;
#'   \code{weight_exponent}, \code{size_influence} and \code{prior_strength}
#'   act on that gain in the same units. Needs the environment matrix that
#'   \code{\link{add_edge_weight}} stores on the graph.
#' @param min_size Integer \eqn{\ge 1}. Smallest allowed piece. Used only by
#'   \code{cut_rule = "ssd"}, where cuts are sequential. Default 1.
#' @param seed Optional integer seed. One stream drives the whole call; the
#'   caller's RNG state is restored on exit. \code{NULL} (default) uses the
#'   current stream, so \code{set.seed()} before the call still works.
#'
#' @return An object of class \code{"pulse_ensemble"} — a list with:
#'   \describe{
#'     \item{\code{memberships}}{Integer matrix (\code{N x N_tree}). Each column
#'       is a region assignment from one spanning tree.}
#'     \item{\code{N_tree}}{Ensemble size used.}
#'     \item{\code{regions_per_partition}}{Number of regions per tree.}
#'     \item{\code{n_members}}{Number of columns.}
#'     \item{\code{weight_exponent}}{\eqn{\gamma} used.}
#'     \item{\code{n_vertices}}{Number of vertices.}
#'     \item{\code{prior_strength}}{The prior strength used.}
#'     \item{\code{chance_coassoc}}{Per-column chance co-association
#'       \eqn{p_{0,t} = \sum_r (n_r / N)^2}; its mean is the reference scale
#'       for cohesion.}
#'     \item{\code{cut_rule}, \code{min_size}}{The cut settings used.}
#'   }
#' @importFrom igraph E ends make_graph components
#' @export
compute_ensemble <- function(trees, 
                             N_tree, 
                             graph, 
                             regions_per_partition, 
                             prior_strength = 0, 
                             weight_exponent = 4, 
                             verbose = TRUE,
                             size_influence = 0,
                             cut_rule = c("ssd", "propensity"),
                             min_size = 5,
                             seed = NULL) {
  
  stopifnot(inherits(trees, "pulse_trees"))
  cut_rule <- match.arg(cut_rule)
  stopifnot(
    length(N_tree) == 1L,
    N_tree >= 2,
    N_tree <= trees$n_trees
  )
  stopifnot(
    length(regions_per_partition) == 1L,
    regions_per_partition >= 2,
    regions_per_partition < trees$n_vertices
  )
  stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L, prior_strength >= 0)
  stopifnot(is.numeric(weight_exponent), length(weight_exponent) == 1L, is.finite(weight_exponent), weight_exponent >= 0)
  if (is.infinite(prior_strength))  stop("`prior_strength = Inf` is not supported here; use the exact-prior ", "short-circuit in get_regions()/tune_regions().", call. = FALSE)
  .local_seed(seed)
  if (verbose) {
    msg <- sprintf(
      "Computing ensemble memberships (N_tree = %d, regions_per_partition = %d",
      N_tree, regions_per_partition
    )
    if (weight_exponent != 1)
      msg <- paste0(msg, sprintf(", weight_exponent = %g", weight_exponent))
    if (prior_strength > 0 && !is.null(trees$within_prior_edge))
      msg <- paste0(msg, sprintf(", lambda_E = %g", prior_strength))
    if (cut_rule != "propensity")
      msg <- paste0(msg, sprintf(", cut_rule = %s", cut_rule))
    message(msg, ") ...")
  }
  
  t0 <- proc.time()[["elapsed"]]
  memb <- .compute_ensemble_internal(
    graph = graph, 
    trees = trees, 
    N_tree = N_tree,
    regions_per_partition = regions_per_partition,
    prior_strength       = prior_strength,
    weight_exponent      = weight_exponent,
    size_influence       = size_influence,
    cut_rule             = cut_rule,
    min_size             = min_size
  )
  if (verbose) message(sprintf("  Done in %s", .format_duration(proc.time()[["elapsed"]] - t0)))
  
  # One constructor for every ensemble, so objects from tune_regions()
  # and from here have the same fields and types.
  .new_ensemble(memb,
                N_tree                = N_tree,
                regions_per_partition = regions_per_partition,
                weight_exponent       = weight_exponent,
                prior_strength        = prior_strength,
                n_vertices            = trees$n_vertices,
                cut_rule              = cut_rule,
                min_size              = min_size)
}

# -------------------------------------------------------------------------*
# -------------------------------------------------------------------------*
# -------------------------------------------------------------------------*


#' Log cutting propensity of the edges of one spanning tree
#'
#' \deqn{\log \pi_e = \gamma \log(w_e / \max w) + \alpha \log s_e
#'       - \lambda \, u \, 1[e \text{ within prior}]}
#' with \eqn{u = \mathrm{sd}(\gamma \log w)} and
#' \eqn{\alpha = \texttt{size\_influence} \cdot u / \mathrm{sd}(\log s)}.
#' Size and prior are both expressed in the unit \eqn{u}, the spread of the
#' environmental term on the log scale the exponential race competes on, so
#' \code{size_influence} and \code{prior_strength} are dimensionless and mean
#' the same on every graph. Zero weights are floored at the smallest positive
#' weight of the tree before taking logs. 
#' If \eqn{u = 0} (\eqn{\gamma = 0}, constant or all-zero weights) cuts are
#' uniform: \eqn{\alpha = 0}, and the prior falls back to \eqn{u = \log 2}, one
#' halving per unit of \eqn{\lambda}, and \code{fallback} is set.
#'
#' @param w Tree edge weights (non-negative, finite).
#' @param size_f Size factor \eqn{n_c (N - n_c) / N} per tree edge.
#' @param within_prior Logical per tree edge, or \code{NULL}.
#' @return List: \code{lp} (log-propensity, max 0), \code{unit} (\eqn{u}),
#'   \code{alpha}, \code{fallback} (logical).
#' @keywords internal
#' @noRd
.cut_log_propensity <- function(w, size_f, weight_exponent = 1,
                                size_influence = 2, prior_strength = 0,
                                within_prior = NULL) {
  # A zero weight (identical neighbours) is treated as the smallest
  # positive weight in the tree. log(0 + eps) = -36 would dominate
  # sd(lw) and change the meaning of size_influence and prior_strength.
  pos <- w > 0
  lw <- if (any(pos)) {
    w_n <- w / max(w)
    w_n[!pos] <- min(w_n[pos])
    weight_exponent * log(w_n)
  } else rep(0, length(w))
  ls   <- log(size_f)
  sd_w <- stats::sd(lw)
  sd_s <- stats::sd(ls)
  ok_w  <- is.finite(sd_w) && sd_w > 0
  ok_s  <- is.finite(sd_s) && sd_s > 0
  # no weight spread: one unit = one halving
  unit  <- if (ok_w) sd_w else log(2)
  
  alpha <- if (ok_w && ok_s) size_influence * unit / sd_s else 0   
  lp    <- lw + alpha * ls
  fallback <- FALSE
  if (prior_strength > 0 && !is.null(within_prior)) {
    within_prior[is.na(within_prior)] <- FALSE
    if (any(within_prior)) {
      fallback <- !ok_w
      lp <- lp - prior_strength * unit * within_prior
    }
  }
  list(lp = lp - max(lp), unit = unit, alpha = alpha, fallback = fallback)
}

#' Internal: compute the N x N_tree ensemble membership matrix
#'
#' Shared by the public \code{compute_ensemble()} and the stability
#' evaluator. For each tree, removes \code{regions_per_partition - 1} edges
#' (weighted sampling without replacement, Efraimidis-Spirakis) and labels the
#' resulting connected components. One cut-set, and one column, per tree.
#' Cut propensities come from \code{.cut_log_propensity()}.
#'
#' @param graph Parent \code{igraph}; read only for \code{cut_rule = "ssd"}
#'   (\code{graph_attr(graph, "env")}).
#' @param trees A \code{"pulse_trees"} object.
#' @param N_tree Number of trees to use.
#' @param regions_per_partition Number of regions per cut-set.
#' @param prior_strength \eqn{\lambda}, in units of sd(\eqn{\gamma} log w).
#' @param weight_exponent \eqn{\gamma}.
#' @param size_influence Size term spread relative to the environmental term.
#' @param cut_rule \code{"propensity"} or \code{"ssd"}; see \code{compute_ensemble()}.
#' @param min_size Smallest allowed piece size. Enforced only under
#'   \code{cut_rule = "ssd"}, where cuts are sequential.
#' @return Integer matrix N x N_tree of region labels.
#' @details Under \code{cut_rule = "ssd"} each of the
#'   \code{regions_per_partition - 1} cuts is drawn by an exponential race on
#'   \deqn{\log \pi_e = \gamma \log g_e + \alpha \log h_e - \lambda u \, 1[e \in \text{prior}],}
#'   recomputed on the current forest after every cut. \eqn{g_e = h_e D_e} is
#'   the Ward SSD reduction of cutting \eqn{e}, \eqn{h_e = n_1 n_2 / (n_1 + n_2)}
#'   and \eqn{D_e} the squared distance between the two side means of
#'   \code{graph_attr(graph, "env")}. Because \eqn{g_e} contains \eqn{h_e}, the
#'   Ward term already favors balanced cuts with exponent \eqn{\gamma};
#'   \code{size_influence} adds to that. \eqn{u = \mathrm{sd}(\gamma \log g)},
#'   \eqn{\alpha = \texttt{size\_influence} \cdot u / \mathrm{sd}(\log h)},
#'   with \eqn{\alpha = 0} when \eqn{u = 0} and the same \eqn{\log 2} prior 
#'   fallback as the propensity rule. Cuts leaving a side smaller than 
#'   \code{min_size} are excluded.
#'   \code{size_influence = 0}, \code{prior_strength = 0} and
#'   \code{min_size = 1} reproduce the plain Ward race draw for draw.
#' @keywords internal
#' @noRd
.compute_ensemble_internal <- function(graph,
                                       trees,
                                       N_tree,
                                       regions_per_partition,
                                       prior_strength  = 0,
                                       weight_exponent = 1,
                                       size_influence  = 2,
                                       cut_rule = "ssd",
                                       min_size = 1L) {
  n_use <- as.integer(N_tree)
  n_reg <- as.integer(regions_per_partition)
  n_v   <- trees$n_vertices
  wpe   <- trees$within_prior_edge
  if (n_reg > n_v) {
    stop(sprintf("regions_per_partition (%d) cannot exceed the number of vertices (%d).",
                 n_reg, n_v), call. = FALSE)
  }
  n_to_cut   <- n_reg - 1L
  n_fallback <- 0L
  min_size   <- .check_min_size(min_size, n_reg, n_v)
  
  ep_all   <- trees$edge_endpoints
  memb_mat <- matrix(0L, nrow = n_v, ncol = n_use)
  if (identical(cut_rule, "ssd")) {
    X <- igraph::graph_attr(graph, "env")
    if (is.null(X))
      stop("cut_rule 'ssd' needs graph_attr(graph, 'env').", call. = FALSE)
    if (anyNA(X))
      stop("cut_rule 'ssd' needs complete graph_attr(graph, 'env'); ",
           "impute or drop rows with NA first.", call. = FALSE)
  }
  for (j in seq_len(n_use)) {
    edge_ids <- .tree_edge_ids(trees, j)
    w        <- trees$edge_weights[edge_ids]
    n_e      <- length(w)
    if (anyNA(w) || any(is.infinite(w)) || any(w < 0))
      stop("Edge weights must be positive and finite numbers", call. = FALSE)
    if (anyNA(edge_ids))
      stop("Failed to map a spanning-tree edge back to the parent graph.",
           call. = FALSE)
    
    if (identical(cut_rule, "ssd")) {
      wp <- if (prior_strength > 0 && !is.null(wpe)) wpe[edge_ids] else NULL
      if (!is.null(wp)) wp[is.na(wp)] <- FALSE
      cut <- .ssd_cut_tree(ep_all[edge_ids, , drop = FALSE], X, n_to_cut,
                           weight_exponent = weight_exponent,
                           size_influence  = size_influence,
                           prior_strength  = prior_strength,
                           within_prior    = wp,
                           min_size        = min_size,
                           tree            = j)
      n_fallback <- n_fallback + cut$fallback
      memb_mat[, j] <- cut$labels
      next
    }
    
    if (n_to_cut >= n_e) {
      stop(sprintf(paste0(
        "Requested %d cuts but the spanning tree has only %d edges; ",
        "reduce regions_per_partition or supply a connected graph."),
        n_to_cut, n_e), call. = FALSE)
    }
    n_c    <- .subtree_sizes(edge_ids, ep_all, n_v = n_v)
    size_f <- n_c * (n_v - n_c) / n_v
    
    pr <- .cut_log_propensity(
      w, size_f,
      weight_exponent = weight_exponent,
      size_influence  = size_influence,
      prior_strength  = prior_strength,
      within_prior    = if (is.null(wpe)) NULL else wpe[edge_ids])
    n_fallback <- n_fallback + pr$fallback
    
    # Exponential race in log space: log(E / pi) = log E - log pi.
    # Same ordering, and the same RNG draws, as E / pi, without
    # underflow when gamma or the weight spread is large.
    key <- log(stats::rexp(n_e)) - pr$lp
    cut <- which(key <= sort(key, partial = n_to_cut)[n_to_cut])
    if (length(cut) > n_to_cut) cut <- cut[seq_len(n_to_cut)]
    memb_mat[, j] <- .labels_from_edges(
      ep_all[edge_ids[-cut], , drop = FALSE], n_v)
  }
  if (n_fallback > 0L)
    warning(sprintf(paste0(
      "%d of %d tree(s) have no spread in weight_exponent * log(weight) ",
      "(weight_exponent = 0 or constant weights). For these, prior_strength ",
      "was applied as halvings (1 unit = log 2), not in standard deviations."),
      n_fallback, n_use), call. = FALSE)
  memb_mat
}


.subtree_sizes <- function(edge_ids, ep_all, n_v) {
  ep  <- ep_all[edge_ids, , drop = FALSE]
  n_e <- nrow(ep)
  
  from <- c(ep[, 1], ep[, 2])
  to   <- c(ep[, 2], ep[, 1])
  eid  <- rep(seq_len(n_e), 2L)
  o    <- order(from)
  from_s <- from[o]; to_s <- to[o]; eid_s <- eid[o]
  ptr <- c(0L, cumsum(tabulate(from_s, n_v)))
  
  root   <- from[1L]
  parent <- integer(n_v)
  pedge  <- integer(n_v)
  bfs    <- integer(n_v)
  seen   <- logical(n_v)
  
  seen[root] <- TRUE
  bfs[1L]    <- root
  head <- 1L; tail <- 1L
  while (head <= tail) {
    v <- bfs[head]; head <- head + 1L
    rng <- (ptr[v] + 1L):ptr[v + 1L]
    if (ptr[v + 1L] > ptr[v]) {
      for (i in rng) {
        u <- to_s[i]
        if (!seen[u]) {
          seen[u]      <- TRUE
          parent[u]    <- v
          pedge[u]     <- eid_s[i]
          tail         <- tail + 1L
          bfs[tail]    <- u
        }
      }
    }
  }
  
  size <- rep(1L, n_v)
  for (i in n_v:2L) {
    v <- bfs[i]
    size[parent[v]] <- size[parent[v]] + size[v]
  }
  
  n_c <- integer(n_e)
  for (v in seq_len(n_v)) if (v != root) n_c[pedge[v]] <- size[v]
  n_c
}

# Rooted-forest state for Ward SSD cut gains.
#
# One DFS (igraph, compiled) gives each vertex a preorder position `pre` and
# its subtree size, so every subtree is one contiguous block of the preorder.
# Subtree sums of X then come from prefix sums, with no interpreted traversal.
# In a forest every edge is a DFS tree edge, so the parent of an edge is the
# endpoint with the smaller preorder position.
#
# ep: m x 2 integer edge matrix of a forest; X: n x p numeric matrix in vertex
# order, no NA. Returns a list of plain vectors (see .ssd_cut_tree()).
.ssd_forest_init <- function(ep, X) {
  X   <- as.matrix(X)
  n_v <- nrow(X); n_e <- nrow(ep); p <- ncol(X)
  # Gains depend on differences of means only. Centring keeps the
  # prefix sums small, so subtree sums lose no precision at large n.
  X   <- X - rep(colMeans(X), each = n_v)
  a   <- as.integer(ep[, 1]); b <- as.integer(ep[, 2])
  g   <- igraph::make_graph(as.vector(rbind(a, b)), n = n_v, directed = FALSE)
  d   <- igraph::dfs(g, root = 1L, mode = "all", unreachable = TRUE,
                     order = TRUE, order.out = TRUE)
  ord  <- as.integer(d$order)
  post <- integer(n_v); post[as.integer(d$order.out)] <- seq_len(n_v)
  pre  <- integer(n_v); pre[ord] <- seq_len(n_v)
  child  <- ifelse(pre[a] > pre[b], a, b)
  parent <- integer(n_v); parent[child] <- ifelse(pre[a] > pre[b], b, a)
  pedge  <- integer(n_v); pedge[child] <- seq_len(n_e)
  # Depth by pointer jumping: about log2(depth) vectorised steps.
  dep <- as.integer(parent > 0L); anc <- parent
  while (any(j <- anc > 0L)) {
    dep[j] <- dep[j] + dep[anc[j]]
    anc[j] <- anc[anc[j]]
  }
  # pre - 1 = ancestors + earlier vertices; post - 1 = earlier vertices
  # + descendants. Hence size = descendants + 1 = post - pre + depth + 1.
  size <- post - pre + dep + 1L
  tot  <- rep(1L, n_v)
  if (n_e) {
    rs <- rowsum(size[child], parent[child])
    ip <- as.integer(rownames(rs))
    tot[ip] <- tot[ip] + as.integer(rs[, 1L])
  }
  if (anyDuplicated(child) || !identical(size, tot))
    stop("Internal error: the cut edges do not form a forest, or ",
         "igraph::dfs() returned an inconsistent order.", call. = FALSE)
  Cx <- matrix(0, n_v + 1L, p)
  for (k in seq_len(p)) Cx[-1L, k] <- cumsum(X[ord, k])
  S  <- Cx[pre + size, , drop = FALSE] - Cx[pre, , drop = FALSE]
  roots <- ord[parent[ord] == 0L]                  # in preorder
  comp  <- roots[findInterval(pre, pre[roots])]    # current component root
  list(pre = pre, ord = ord, size0 = size, size = as.numeric(size), S = S,
       comp = comp, child = child, pedge = pedge)
}

# Ward gain, h = n1 n2 / (n1 + n2) and min(n1, n2) of cutting the edges whose
# lower endpoints are `v`, each within its current component.
.ssd_edge_gains <- function(v, comp, size, S) {
  r  <- comp[v]
  n1 <- size[v]; n2 <- size[r] - n1
  Sv <- S[v, , drop = FALSE]
  d  <- Sv / n1 - (S[r, , drop = FALSE] - Sv) / n2
  h  <- n1 * n2 / (n1 + n2)
  list(gain = h * rowSums(d^2), h = h, n_min = pmin(n1, n2))
}

# Ward SSD reduction of cutting each edge of a forest, within its own component.
# ep: m x 2 integer edge matrix; X: n x p matrix in vertex order.
# Returns list: gain = h * D, h = n1 n2 / (n1 + n2), n_min = min(n1, n2),
# all in the row order of ep.
.forest_cut_gains <- function(ep, X) {
  st <- .ssd_forest_init(ep, X)
  .ssd_edge_gains(st$child, st$comp, st$size, st$S)
}

# Sequential SSD cuts of one spanning tree.
#
# Same race as before, draw for draw: one exponential per remaining edge, in
# the row order of `ep` with cut rows removed. Only the gains change after a
# cut, and only inside the component that was split: sizes and sums of the
# cut vertex's ancestors drop by its subtree, its descendants move to a new
# component, and every edge of the two pieces is rescored. Other components
# keep their gains. Cost per cut is linear in the split component, with no
# interpreted traversal.
#
# Returns list(labels, fallback); labels are numbered like
# igraph::components() (by first vertex), as .labels_from_edges() does.
.ssd_cut_tree <- function(ep, X, n_to_cut, weight_exponent = 1,
                          size_influence = 2, prior_strength = 0,
                          within_prior = NULL, min_size = 1L, tree = 1L) {
  st    <- .ssd_forest_init(ep, X)
  pre   <- st$pre; ord <- st$ord; size0 <- st$size0; child <- st$child
  pedge <- st$pedge; size <- st$size; S <- st$S; comp <- st$comp
  rm(st)
  n_e   <- nrow(ep)
  alive <- rep(TRUE, n_e)
  fc    <- .ssd_edge_gains(child, comp, size, S)
  gain  <- fc$gain; h <- fc$h; n_min <- fc$n_min
  wp    <- within_prior
  fb    <- FALSE
  for (s in seq_len(n_to_cut)) {
    live <- which(alive)
    lp   <- weight_exponent * log(pmax(gain[live], .Machine$double.xmin))
    u    <- stats::sd(lp)                 # unit: sd(gamma log g)
    ok_u <- is.finite(u) && u > 0
    if (!ok_u) u <- log(2)
    if (size_influence > 0 && ok_u) {
      lh   <- log(h[live])
      sd_h <- stats::sd(lh)
      if (is.finite(sd_h) && sd_h > 0)
        lp <- lp + size_influence * u / sd_h * lh
    }
    if (!is.null(wp) && any(wp[live])) {
      if (!ok_u) fb <- TRUE
      lp <- lp - prior_strength * u * wp[live]
    }
    if (min_size > 1L) lp[n_min[live] < min_size] <- -Inf
    if (!any(is.finite(lp)))
      stop(sprintf(paste0(
        "Tree %d: no cut keeps both sides >= min_size (%d) ",
        "at cut %d of %d. Lower min_size or regions_per_partition."),
        tree, min_size, s, n_to_cut), call. = FALSE)
    e <- live[which.min(log(stats::rexp(length(live))) - lp)]
    alive[e] <- FALSE
    if (s == n_to_cut) break                  # no gains needed after the last cut
    # --- split the component of the cut edge ---------------------*
    v  <- child[e]; rt <- comp[v]; pv <- pre[v]
    M  <- ord[pre[rt] + seq_len(size0[rt]) - 1L]
    M  <- M[comp[M] == rt]                    # current component of rt
    pM <- pre[M]
    up <- M[pM < pv & pM + size0[M] > pv]     # ancestors of v
    dn <- M[pM >= pv & pM < pv + size0[v]]    # v and its descendants
    size[up] <- size[up] - size[v]
    S[up, ]  <- S[up, , drop = FALSE] - rep(S[v, ], each = length(up))
    comp[dn] <- v
    re <- pedge[M]; re <- re[re > 0L]; re <- re[alive[re]]
    if (length(re)) {
      fc <- .ssd_edge_gains(child[re], comp, size, S)
      gain[re] <- fc$gain; h[re] <- fc$h; n_min[re] <- fc$n_min
    }
  }
  list(labels = .labels_from_edges(ep[alive, , drop = FALSE], nrow(X)), fallback = fb)
}

# Validate min_size: one integer >= 1, and k regions of that size must fit.
.check_min_size <- function(min_size, k, n_v) {
  if (length(min_size) != 1L || is.na(min_size) || min_size < 1 ||
      min_size != round(min_size))
    stop("min_size must be one whole number >= 1.", call. = FALSE)
  min_size <- as.integer(min_size)
  if (min_size * k > n_v)
    stop(sprintf("min_size (%d) x regions (%d) exceeds the number of units (%d).",
                 min_size, k, n_v), call. = FALSE)
  min_size
}


#' Consensus Regions from Ensemble Memberships (salso + piece composition)
#'
#' Finds a hard consensus partition of the \eqn{N} units into \code{k} regions
#' by minimising expected variation of information over the ensemble with
#' \pkg{salso}, then derives graded memberships \eqn{\Pi} from the regional
#' composition of each unit's ensemble pieces.
#'
#' @details
#' \strong{Memberships \eqn{\Pi}.} For each ensemble partition \eqn{t}, take
#' the share of unit \eqn{i}'s piece-mates (the other units in its piece) that
#' belong to each consensus region, then average over partitions:
#' \deqn{\Pi_{ik} = \frac{1}{M_i} \sum_t
#'   \frac{\#\{\text{piece-mates of } i \text{ in region } k\}}
#'        {\#\{\text{piece-mates of } i\}},}
#' where \eqn{M_i} (returned as \code{pi_support}) counts the partitions in
#' which \eqn{i} has at least one piece-mate. Rows sum to 1 by construction,
#' there is no tuning parameter, and a boundary unit whose pieces fall half
#' into each of two regions gets 0.5 and 0.5 whatever the region sizes. Units
#' with \eqn{M_i = 0} are filled from graph neighbours (see compute_consensus()
#'  Details).
#'
#' \strong{Cohesion.} Let \eqn{c_{im}} be the fraction of ensemble partitions
#' in which \eqn{i} and \eqn{m} share a piece. The \emph{cohesion} of \eqn{i}
#' to region \eqn{j} is
#' \deqn{A_{ij} = \frac{1}{|C_j \setminus \{i\}|} \sum_{m \in C_j,\, m \ne i} c_{im}.}
#' Rows do not sum to one and keep the absolute level of co-association, so
#' cohesion answers "is there structure here at all?" rather than "which
#' region?". Because a piece holds about \eqn{N/R} units, cohesion shrinks
#' roughly as one over the region size; do not row-normalise it into
#' memberships (the pre-\eqn{\Pi} behaviour, which favoured small regions).
#'
#' \strong{Order of operations.} salso (plus merging and, above
#' \code{large_n_threshold}, cohesion-based refinement), then contiguity repair
#' guided by \eqn{\Pi} of the unrepaired labels, then \eqn{\Pi} and every
#' per-unit summary recomputed on the final, repaired labels. Repair uses
#' \eqn{\Pi} rather than cohesion because regrowth compares scores across
#' regions, where cohesion favours small regions.
#'
#' Everything is done through sparse indicator products and never materialises
#' the \eqn{N \times N} similarity matrix.
#'
#' @param ensemble An \code{"pulse_ensemble"} object or a bare integer
#'   matrix (\eqn{N \times M}).
#' @param final_regions Integer, number of consensus regions \eqn{k},
#'   \eqn{2 \le k < N}.
#' @param graph The parent \code{igraph}; needed for contiguity enforcement.
#' @param refine_boundary Logical. After contiguity repair, move boundary units
#'   to an adjacent region when that lowers the within-region sum of squares of
#'   \code{graph_attr(graph, "env")}, never disconnecting a region or shrinking
#'   it below \code{min_size}. Default \code{TRUE}. Not to be confused with
#'   \code{refine_iter}, which counts single-unit loss moves on the sampled path.
#' @param min_size Integer \eqn{\ge 1}. Smallest region size the boundary
#'   refinement may leave. Default 1.
#' @param large_n_threshold Above this \eqn{N}, use the sampled path.
#' @param sample_size Units passed to \pkg{salso} on the sampled path.
#' @param refine_iter Maximum rounds of single-unit moves on the VI loss
#'   (sampled path only). Each round starts from a full loss-gap pass.
#' @param enforce_contiguity Logical. Reattach fragments so every region is a
#'   single connected component of \code{graph}.
#' @param n_runs,n_cores Passed to \code{salso::salso()}.
#' @param report_native Logical. Also run salso without a target region count
#'   and return its partition as \code{salso_native}. Costs one extra salso
#'   run. Default \code{TRUE}; \code{tune_regions()} turns it off per
#'   configuration and on for the winner.
#' @param verbose Logical.
#'
#' @return A \code{"pulse_consensus"} object with
#'   \describe{
#'     \item{\code{memberships}}{\eqn{N \times k} matrix \eqn{\Pi}; rows sum to 1; rows of units with \code{pi_support == 0} are filled from graph neighbours (see Details).}
#'     \item{\code{hard_labels}}{Final labels (after contiguity repair).}
#'     \item{\code{certainty}}{\eqn{\Pi_{i, \mathrm{hard}_i}}.}
#'     \item{\code{margin}}{Best minus second-best entry of \eqn{\Pi}.}
#'     \item{\code{entropy}}{Row entropy of \eqn{\Pi}, normalized by \eqn{\log k}.}
#'     \item{\code{pi_support}}{\eqn{M_i}, partitions in which \eqn{i} had a piece-mate.}
#'     \item{\code{cohesion}}{Raw \eqn{A}; not a probability.}
#'     \item{\code{silhouette}, \code{avg_sil_width}}{Co-association silhouette, from cohesion.}
#'     \item{\code{expected_vi}, \code{region_coassoc}}{Mean VI(draw, consensus); \eqn{k \times k} mean co-association.}
#'     \item{\code{salso_a}}{The generalized-VI parameter \eqn{a} salso settled on (on the sample, if sampled).}
#'     \item{\code{expected_loss}}{Generalized VI of the final labels at \code{salso_a}, per unit, in bits.}
#'     \item{\code{loss_gap}, \code{loss_alt}}{Change in total loss (bits) if the unit alone moved to \code{loss_alt}; negative means a single move would lower the loss.}
#'     \item{\code{label_step}}{Factor: last step that changed the label (salso, merge, assign, refine, contiguity).}
#'     \item{\code{salso_info}, \code{refine_info}}{salso's own fit summary; on the sampled path, VI-move convergence and per-round log (else \code{NULL}).}
#'     \item{\code{salso_native}}{\code{NULL}, or a list: \code{k}, the
#'       region count salso picks under VI with \eqn{a = 1} and no target;
#'       \code{labels}, its partition (length N, \code{NA} for units outside
#'       the sample, no contiguity repair); \code{sampled}; \code{n_units};
#'       \code{capped}, \code{TRUE} when \code{k} reached the largest piece
#'       count among the ensemble partitions, the most salso searches, so
#'       \code{k} is a lower bound.
#'       This count can depend strongly on \code{regions_per_partition}; check
#'       with \code{\link{check_native_regions}} before reading it as a
#'       property of the landscape.}
#'   }
#'   plus \code{k}, \code{sampled}, \code{n_eff},
#'   \code{n_members}, \code{ensemble_meta}.
#' @importFrom Matrix sparseMatrix crossprod
#' @export
compute_consensus <- function(ensemble, 
                              final_regions, 
                              graph = NULL,
                              large_n_threshold  = 20000,
                              sample_size        = 5000L,
                              refine_iter        = 10L,
                              enforce_contiguity = !is.null(graph),
                              refine_boundary    = !is.null(graph),
                              min_size           = 5,
                              n_runs             = 16L,
                              n_cores            = 1L,
                              report_native      = TRUE,
                              verbose            = TRUE) {
  
  k    <- final_regions
  meta <- .ensemble_meta(ensemble)
  memb <- meta$memb_mat
  N <- nrow(memb); M <- ncol(memb)
  stopifnot(length(k) == 1L, k >= 2, k < N)
  if (enforce_contiguity && is.null(graph))
    stop("`enforce_contiguity = TRUE` requires `graph`.", call. = FALSE)
  if (!is.logical(refine_boundary) || length(refine_boundary) != 1L || is.na(refine_boundary))
    stop("`refine_boundary` must be TRUE or FALSE.", call. = FALSE)
  if (refine_boundary && is.null(graph))
    stop("`refine_boundary = TRUE` requires `graph`.", call. = FALSE)
  
  sampled <- N > large_n_threshold
  if (verbose) message(sprintf(
    "Consensus clustering (k = %d, N = %d, members = %d, %s) ...",
    k, N, M,
    if (sampled) sprintf("salso on %d-unit sample + cohesion assignment",
                         min(sample_size, N))
    else "salso on all units"))
  t0 <- proc.time()[["elapsed"]]
  
  # --- hard consensus partition ----------------------------------------*
  fit  <- .consensus_partition(memb, k, large_n_threshold, sample_size,
                               refine_iter, n_runs, n_cores, verbose)
  res <- .finish_consensus(meta, fit, k, graph = graph,
                           enforce_contiguity = enforce_contiguity,
                           large_n_threshold  = large_n_threshold,
                           sample_size        = sample_size,
                           n_runs = n_runs, n_cores = n_cores,
                           refine_boundary = refine_boundary, min_size = min_size,
                           report_native = report_native, verbose = verbose)
  if (verbose) message(sprintf("  Done in %s",
                               .format_duration(proc.time()[["elapsed"]] - t0)))
  res
}

#' Everything compute_consensus() does after the hard partition.
#'
#' Deterministic given (ensemble, fit), except the optional salso_native report,
#' which runs last so it cannot shift the labels. tune_regions() stores `fit` per
#' configuration and calls this again on the winner, so the returned object has
#' the partition that was scored, not a refit.
#' @keywords internal
#' @noRd
.finish_consensus <- function(meta, fit, k, graph = NULL,
                              enforce_contiguity = !is.null(graph),
                              large_n_threshold = 20000, sample_size = 5000,
                              n_runs = 16, n_cores = 1,
                              refine_boundary = FALSE, min_size = 1L,
                              report_native = TRUE, verbose = TRUE) {
  memb <- meta$memb_mat
  N <- nrow(memb); M <- ncol(memb)
  sampled <- N > large_n_threshold
  hard <- fit$labels
  step <- as.integer(fit$step)                 # stored fits keep step as raw
  a    <- fit$salso$a
  # --- stats on the consensus labels (Pi and loss terms included) ------*
  st  <- .ensemble_stats(memb, hard, k, pi = TRUE, gvi_a = a)
  coh <- .cohesion_from_stats(st, seq_len(N), hard, N)
  
  # --- contiguity (Pi-guided), then recompute on the repaired labels ----*
  # Regrowth ranks a freed unit's candidate regions against each other,
  # so it needs a size-neutral score: Pi, not cohesion.
  if (enforce_contiguity) {
    hard2 <- .enforce_contiguity(graph, hard, st$Pi)
    n_moved <- sum(hard2 != hard)
    if (n_moved > 0L) {
      if (verbose) message(sprintf("  contiguity: %d unit(s) reattached", n_moved))
      step[hard2 != hard] <- 5L                     # "contiguity"
      hard <- hard2
      st   <- .ensemble_stats(memb, hard, k, pi = TRUE, gvi_a = a)
      coh  <- .cohesion_from_stats(st, seq_len(N), hard, N)
    }
  }
  # --- boundary refinement: contiguity-preserving SSD descent ---------*
  # refine_boundary = TRUE. Needs contiguous input, so it runs after
  # contiguity repair; moves never disconnect a region.
  if (refine_boundary) {
    if (is.null(graph))
      stop("`refine_boundary = TRUE` requires `graph`.", call. = FALSE)
    X <- igraph::graph_attr(graph, "env")
    if (is.null(X))
      stop("`refine_boundary = TRUE` needs graph_attr(graph, 'env'); ",
           "build the graph with add_edge_weight().", call. = FALSE)
    if (anyNA(X))
      stop("`refine_boundary = TRUE` needs complete graph_attr(graph, 'env'); ",
           "impute or drop rows with NA first.", call. = FALSE)
    hard2   <- .refine_regions(hard, graph, X, min_size = min_size)$z
    n_moved <- sum(hard2 != hard)
    if (n_moved > 0L) {
      if (verbose) message(sprintf("  boundary refinement: %d unit(s) moved", n_moved))
      step[hard2 != hard] <- 6L                     # "boundary"
      hard <- hard2
      st   <- .ensemble_stats(memb, hard, k, pi = TRUE, gvi_a = a)
      coh  <- .cohesion_from_stats(st, seq_len(N), hard, N)
    }
  }
  # --- memberships and per-unit summaries, all from final labels -------*
  n_na <- sum(st$pi_support == 0L)
  Pi   <- .fill_unsupported_pi(st$Pi, hard, graph)
  sm   <- .membership_summaries(Pi, hard)
  if (n_na > 0L)
    warning(sprintf(paste0(
      "%d unit(s) were a one-unit piece in every partition (pi_support == 0). ",
      "Their memberships are the mean of their graph neighbours' rows, or ",
      "one-hot on the hard label without a graph. Raise N_tree or lower ",
      "regions_per_partition to reduce them."), n_na), call. = FALSE)
  
  # Co-association silhouette stays on cohesion (a co-association quantity).
  own_idx <- cbind(seq_len(N), hard)
  own     <- coh[own_idx]
  oth     <- coh; oth[own_idx] <- -Inf
  nxt     <- oth[cbind(seq_len(N), max.col(oth, ties.method = "first"))]
  sil     <- (own - nxt) / pmax(own, nxt, .Machine$double.eps)
  
  # k x k mean co-association; diagonal excludes self-pairs.
  n_k <- st$n_k
  B   <- st$P / (M * outer(n_k, n_k))
  dg  <- (diag(st$P) - M * n_k) / (M * n_k * pmax(n_k - 1, 1))
  dg[n_k < 2] <- NA_real_
  diag(B) <- dg
  
  # --- loss gaps on the final labels, at salso's a ----------------------*
  gap <- if (is.null(st$gvi)) list(gap = rep(NA_real_, N), alt = rep(NA_integer_, N))
  else .vi_move_gap(st, hard)
  
  # --- salso's own answer: no target k (VI, a = 1) ----------------------*
  # Run after everything above so it does not shift the RNG stream that
  # produced the k-region result.
  native <- NULL
  if (report_native) {
    idx_nat <- if (sampled) sort(sample.int(N, min(as.integer(sample_size), N)))
    else seq_len(N)
    native <- .salso_native(memb, idx_nat, n_runs, n_cores)
    if (verbose) message(sprintf(paste0(
      "  salso without a target k (VI, a = 1): %d regions%s; requested k = %d ",
      "(regions_per_partition = %s). See check_native_regions()."),
      native$k,
      paste0(if (native$sampled) sprintf(" on a %d-unit sample", native$n_units) else "",
             if (native$capped) " (at salso's cluster cap: a lower bound)" else ""),
      k, format(meta$regions_per_partition)))
  }
  
  structure(
    list(
      memberships    = Pi,
      hard_labels  = as.integer(hard),
      certainty      = sm$certainty,
      margin         = sm$margin,
      entropy        = sm$entropy,
      pi_support     = st$pi_support,
      cohesion       = coh,
      silhouette     = sil,
      avg_sil_width  = mean(sil, na.rm = TRUE),
      expected_vi    = mean(st$vi),
      salso_a        = a,
      expected_loss  = if (is.null(st$gvi)) NA_real_ else st$gvi$expected_loss,
      loss_gap       = gap$gap,
      loss_alt       = gap$alt,
      label_step     = factor(.LABEL_STEPS[step], levels = .LABEL_STEPS),
      salso_info     = fit$salso,
      refine_info    = fit$moves,
      refine_boundary = refine_boundary,
      region_coassoc = B,
      salso_native   = native,
      k              = as.integer(k),
      sampled        = sampled,
      n_eff          = 1 / meta$chance_coassoc,
      n_members      = M,
      ensemble_meta  = meta
    ),
    class = "pulse_consensus"
  )
}

#' salso with no target region count (VI, a = 1) on a subset of units.
#' @keywords internal
#' @noRd
.salso_native <- function(memb, idx = seq_len(nrow(memb)), n_runs = 16L, n_cores = 1L) {
  draws <- t(memb[idx, , drop = FALSE]); storage.mode(draws) <- "integer"
  # salso searches at most as many clusters as the largest draw has and
  # warns when its estimate reaches that cap. The count is then a lower
  # bound; record that as `capped` instead of passing the warning on.
  capped <- FALSE
  est <- withCallingHandlers(
    as.integer(salso::salso(draws, loss = salso::VI(),
                            nRuns = n_runs, nCores = n_cores)),
    warning = function(w) {
      if (grepl("maximum possible number of clusters", conditionMessage(w))) {
        capped <<- TRUE
        invokeRestart("muffleWarning")
      }
    })
  lab <- rep(NA_integer_, nrow(memb))
  lab[idx] <- match(est, sort(unique(est)))
  list(k = length(unique(est)), labels = lab,
       sampled = length(idx) < nrow(memb), n_units = length(idx),
       capped = capped)
}

#' Region count salso chooses on its own, across ensemble settings
#'
#' Cuts one set of trees at each \code{regions_per_partition} (and slices each
#' ensemble to each \code{N_tree}), then runs \pkg{salso} without a target
#' region count (VI loss, \eqn{a = 1}) on one fixed sample of units. Rows differ
#' only in the ensemble, so the table shows whether salso's count follows the
#' data or follows the ensemble settings.
#'
#' @inheritParams compute_ensemble
#' @param regions_per_partition Integer vector of piece counts to compare.
#' @param N_tree Integer vector of ensemble sizes, each \code{<= trees$n_trees}.
#'   One ensemble is built per \code{regions_per_partition} at
#'   \code{max(N_tree)} and sliced.
#' @param sample_size Units passed to salso; the same sample for every row.
#' @param cut_rule,min_size As in \code{\link{compute_ensemble}}.
#' @param n_runs,n_cores Passed to \code{salso::salso()}. Low \code{n_runs}
#'   keeps the sweep fast; confirm a final value with the default 16.
#' @param seed Optional integer; the caller's RNG state is restored on exit.
#' @details On a sample, small regions can be missed, so \code{native_regions}
#'   is a lower bound when \code{sample_size < N}. Compare rows with each other,
#'   not with the full-data count.
#' @return A data frame with \code{regions_per_partition}, \code{N_tree},
#'   \code{n_eff}, \code{native_regions}, \code{capped} (\code{TRUE} when
#'   salso hit its cluster cap, so \code{native_regions} is a lower bound),
#'   \code{ratio} (\code{native_regions / regions_per_partition}) and
#'   \code{seconds}.
#' @examples
#' g <- build_network(pulse_grid, verbose = FALSE)
#' g <- add_edge_weight(g, id_col = "ID", variables = c("temp", "prec"),
#' verbose = FALSE)
#' trees <- sample_spanning_trees(g, n = 200, seed = 1)
#' sw <- check_native_regions(trees, g, regions_per_partition = c(5,6),
#'                            N_tree = c(20,25), seed = 1)
#' plot(native_regions ~ regions_per_partition, data = sw, log = "xy")
#' 
#' @export
check_native_regions <- function(trees, graph,
                                 regions_per_partition = c(10L, 25L, 50L, 100L),
                                 N_tree              = trees$n_trees,
                                 weight_exponent     = 4,
                                 size_influence      = 0,
                                 prior_strength      = 0,
                                 cut_rule            = c("ssd", "propensity"),
                                 min_size            = 5,
                                 sample_size         = 2000L,
                                 n_runs              = 4L,
                                 n_cores             = 1L,
                                 seed                = NULL,
                                 verbose             = TRUE) {
  stopifnot(inherits(trees, "pulse_trees"))
  cut_rule <- match.arg(cut_rule)
  .local_seed(seed)
  R_vals <- sort(unique(as.integer(regions_per_partition)))
  N_vals <- sort(unique(as.integer(N_tree)))
  N      <- trees$n_vertices
  stopifnot(all(N_vals >= 2L), max(N_vals) <= trees$n_trees,
            all(R_vals >= 2L), max(R_vals) < N)
  idx <- if (sample_size < N) sort(sample.int(N, as.integer(sample_size))) else seq_len(N)
  
  rows <- vector("list", length(R_vals) * length(N_vals)); r <- 0L
  for (R in R_vals) {
    memb <- .compute_ensemble_internal(
      graph = graph, trees = trees, N_tree = max(N_vals),
      regions_per_partition = R, prior_strength = prior_strength,
      weight_exponent = weight_exponent,
      size_influence = size_influence,
      cut_rule = cut_rule, min_size = min_size)
    p0     <- .chance_coassoc_by_col(memb)
    memb_s <- memb[idx, , drop = FALSE]
    rm(memb)
    for (nt in N_vals) {
      cols <- seq_len(nt)
      t0   <- proc.time()[["elapsed"]]
      nat   <- .salso_native(memb_s[, cols, drop = FALSE],
                             n_runs = n_runs, n_cores = n_cores)
      k_nat <- nat$k
      sec  <- proc.time()[["elapsed"]] - t0
      r <- r + 1L
      rows[[r]] <- data.frame(regions_per_partition = R, N_tree = nt,
                              n_eff = 1 / mean(p0[cols]),
                              native_regions = k_nat,
                              capped = nat$capped,
                              ratio = k_nat / R, seconds = sec)
      if (verbose) message(sprintf(
        "  regions_per_partition = %d, N_tree = %d: %d regions (%.1f s)",
        R, nt, k_nat, sec))
    }
  }
  do.call(rbind, rows)
}

#' Coerce ensemble input to a meta list (fixes the bare-matrix bug of the
#' old compute_consensus, which read ensemble$memberships from a matrix).
#' @keywords internal
#' @noRd
.ensemble_meta <- function(ensemble) {
  if (inherits(ensemble, "pulse_ensemble")) {
    list(memb_mat              = ensemble$memberships,
         regions_per_partition = ensemble$regions_per_partition,
         N_tree                = ensemble$N_tree,
         weight_exponent       = ensemble$weight_exponent,
         prior_strength        = ensemble$prior_strength,
         chance_coassoc        = mean(ensemble$chance_coassoc))
  } else if (is.matrix(ensemble)) {
    if (!is.numeric(ensemble) || anyNA(ensemble) ||
        any(ensemble != round(ensemble)))
      stop("A bare `ensemble` matrix must hold whole-number labels without NA.",
           call. = FALSE)
    # Labels index sparse-matrix columns downstream, so they must
    # start at 1 (e.g. 0-based labels from other software).
    if (min(ensemble) < 1)
      ensemble <- matrix(apply(ensemble, 2L, function(v) match(v, sort(unique(v)))),
                         nrow = nrow(ensemble))
    storage.mode(ensemble) <- "integer"
    list(memb_mat              = ensemble,
         regions_per_partition = NA_integer_,
         N_tree                = NA_integer_,
         weight_exponent       = NA_real_,
         prior_strength        = NA_real_,
         chance_coassoc        = mean(.chance_coassoc_by_col(ensemble)))
  } else {
    stop("`ensemble` must be an 'pulse_ensemble' object or an integer matrix.",
         call. = FALSE)
  }
}

.xlogx <- function(p) { out <- p * log(p); out[!(p > 0)] <- 0; out }
.xlog2x <- function(p) { out <- p * log2(p); out[!(p > 0)] <- 0; out }

# Steps that can set a consensus label, in pipeline order. `label_step` records
# the LAST step that changed a unit's label.
.LABEL_STEPS <- c("salso", "merge", "assign", "refine", "contiguity", "boundary")

#' Fill Pi rows of units that never shared a piece
#'
#' A unit with pi_support == 0 gets the mean Pi row of its graph neighbours
#' that have one, in synchronous passes outward from supported units, so the
#' result does not depend on unit order. Units left without a filled
#' neighbour (no graph, or no supported unit in reach) get a one-hot row on
#' their hard label. pi_support stays 0 and flags every filled row.
#' @keywords internal
#' @noRd
.fill_unsupported_pi <- function(Pi, hard, graph = NULL) {
  miss <- which(rowSums(is.na(Pi)) > 0L)
  if (!length(miss)) return(Pi)
  if (!is.null(graph) && igraph::vcount(graph) == nrow(Pi)) {
    A    <- igraph::as_adjacency_matrix(graph, sparse = TRUE, attr = NULL)[miss, , drop = FALSE]
    done <- !(seq_len(nrow(Pi)) %in% miss)
    repeat {
      todo <- which(!done[miss])
      if (!length(todo)) break
      D  <- A[todo, , drop = FALSE] %*% Matrix::Diagonal(x = as.numeric(done))
      nd <- Matrix::rowSums(D)
      if (!any(nd > 0)) break
      P0 <- Pi; P0[!done, ] <- 0
      ok <- nd > 0
      Pi[miss[todo[ok]], ] <- as.matrix(D[ok, , drop = FALSE] %*% P0) / nd[ok]
      done[miss[todo[ok]]] <- TRUE
    }
  }
  left <- which(rowSums(is.na(Pi)) > 0L)
  if (length(left)) { Pi[left, ] <- 0; Pi[cbind(left, hard[left])] <- 1 }
  Pi
}

#' Co-labelling statistics of an ensemble against a reference labelling
#'
#' With A_t the N x R indicator of ensemble column t and Z the N x k indicator
#' of the reference labels (zero rows for non-reference units):
#'   S = sum_t A_t (A_t' Z)      N x k  co-label counts  (= M * PSM %*% Z)
#'   P = sum_t (A_t' Z)' (A_t' Z) k x k  co-label pair counts
#'   vi_t = VI(column t, reference) = 2 H(joint) - H(t) - H(ref)
#'
#' With `pi = TRUE` it also returns the composition matrix Pi (Section 6):
#'   Pi_ik = (1 / M_i) sum_t (T_t\\[p,k\\] - Z_ik) / (s_p - 1),  p = piece of i in t,
#' summed over the M_i columns in which i has at least one piece-mate
#' (s_p > 1). T_t = A_t' Z holds region counts per piece and s_p its row sums,
#' so Pi costs one extra sparse product per block: A_t %*% \\[D_w T_t | w | 1(s>1)\\]
#' with w_p = 1 / (s_p - 1). Rows sum to 1 by construction; rows with M_i = 0
#' (singleton piece in every column) are NA. Pi needs every unit labelled,
#' because piece-mate shares are only defined against a full reference.
#'
#' With `gvi_a` (a number in \\[0, 2\\]) it also returns `gvi`, the pieces of
#' salso's generalised VI at that a (log base 2, as salso):
#'   E\\[L\\] = ((2 - a) sum_k f(n_k) + (1/M) sum_t \\[a sum_p f(s_tp)
#'          - 2 sum_{p,k} f(n_tpk)\\]) / N,      f(x) = x log2 x,
#' plus Dp\\[i, s\\] = sum_t f(n_t,p(i),s + 1) - f(n_t,p(i),s) and the matching Dm
#' for removal, the per-draw terms of a one-unit move. Both are sparse in T_t
#' (f(1) - f(0) = 0), so they cost two more sparse products per block.
#'
#' Columns are processed in blocks so that the sparse indicator never exceeds
#' `budget` non-zeros (peak memory ~ 12 bytes * budget).
#'
#' @keywords internal
#' @noRd
.ensemble_stats <- function(memb, ref_lab, k, ref_idx = NULL, budget = 5e7,
                            pi = FALSE, gvi_a = NULL) {
  N <- nrow(memb); M <- ncol(memb)
  if (is.null(ref_idx)) ref_idx <- seq_len(N)
  stopifnot(length(ref_idx) == length(ref_lab), all(ref_lab >= 1L), all(ref_lab <= k))
  full <- length(ref_idx) == N
  if (pi && !full)
    stop("Pi needs a full reference labelling (ref_idx = all units).", call. = FALSE)
  do_gvi <- !is.null(gvi_a) && is.finite(gvi_a)
  if (do_gvi && !full)
    stop("Loss terms need a full reference labelling (ref_idx = all units).", call. = FALSE)
  R    <- as.integer(max(memb))                    # label stride
  Z    <- Matrix::sparseMatrix(i = ref_idx, j = ref_lab, x = 1, dims = c(N, k))
  n_k  <- tabulate(ref_lab, k)
  
  S  <- matrix(0, N, k)
  P  <- matrix(0, k, k)
  vi <- numeric(M)
  Hk <- -sum(.xlogx(n_k / N))
  if (pi) { Pn <- matrix(0, N, k + 2L) }           # [numerator | sum w | M_i]
  if (do_gvi) { Dp <- Dm <- matrix(0, N, k); sum_joint <- 0; sum_piece <- 0 }
  
  b <- max(1L, min(M, as.integer(budget %/% N)))
  for (start in seq(1L, M, by = b)) {
    cols <- start:min(M, start + b - 1L)
    m    <- length(cols)
    A <- Matrix::sparseMatrix(
      i = rep.int(seq_len(N), m),
      j = as.integer(memb[, cols]) + R * rep(seq_len(m) - 1L, each = N),
      x = 1, dims = c(N, R * m))
    Tm <- Matrix::crossprod(A, Z)                       # (R*m) x k
    S  <- S + as.matrix(A %*% Tm)
    P  <- P + as.matrix(Matrix::crossprod(Tm))
    if (full) {
      Td  <- as.matrix(Tm) / N
      grp <- rep(seq_len(m), each = R)
      Hj  <- -rowsum(rowSums(.xlogx(Td)), grp)[, 1]
      Ht  <- -rowsum(.xlogx(rowSums(Td)), grp)[, 1]
      vi[cols] <- 2 * Hj - Ht - Hk
    }
    if (pi) {
      s_p <- Matrix::rowSums(Tm)                  # piece sizes
      has <- s_p > 1                              # piece has mates
      w   <- ifelse(has, 1 / pmax(s_p - 1, 1), 0)
      Tw  <- cbind(Matrix::Diagonal(x = w) %*% Tm, w, as.numeric(has))
      Pn  <- Pn + as.matrix(A %*% Tw)
    }
    if (do_gvi) {
      Tc <- methods::as(methods::as(Tm, "CsparseMatrix"), "generalMatrix")
      x  <- Tc@x
      sum_joint <- sum_joint + sum(.xlog2x(x))
      sum_piece <- sum_piece + sum(.xlog2x(Matrix::rowSums(Tc)))
      Gp <- Tc; Gp@x <- .xlog2x(x + 1) - .xlog2x(x)
      Gm <- Tc; Gm@x <- .xlog2x(x - 1) - .xlog2x(x)
      Dp <- Dp + as.matrix(A %*% Gp)
      Dm <- Dm + as.matrix(A %*% Gm)
    }
  }
  out <- list(S = S, P = P, vi = if (full) vi else NULL, n_k = n_k, M = M)
  if (pi) {
    Mi  <- Pn[, k + 2L]
    num <- Pn[, seq_len(k), drop = FALSE]
    own <- cbind(seq_len(N), ref_lab[order(ref_idx)])
    num[own] <- num[own] - Pn[, k + 1L]               # remove self from piece counts
    Pi  <- num / pmax(Mi, 1)
    Pi[Mi == 0, ] <- NA_real_
    Pi[!is.na(Pi) & Pi < 0 & Pi > -1e-12] <- 0        # floating-point dust
    dimnames(Pi) <- NULL
    out$Pi         <- Pi
    out$pi_support <- as.integer(round(Mi))
  }
  if (do_gvi) {
    out$gvi <- list(
      a             = gvi_a,
      expected_loss = ((2 - gvi_a) * sum(.xlog2x(n_k)) +
                         (gvi_a * sum_piece - 2 * sum_joint) / M) / N,
      Dp = Dp, Dm = Dm)
  }
  out
}

#' Cohesion: mean co-association of every unit with each reference region,
#' self excluded. Not a probability: rows do not sum
#' to 1, and the value shrinks roughly as 1 / region size. Empty regions and
#' singleton self-cohesions are set to 0.
#' @keywords internal
#' @noRd
.cohesion_from_stats <- function(st, ref_idx, ref_lab, N) {
  S <- st$S; k <- ncol(S); M <- st$M
  self <- cbind(ref_idx, ref_lab)
  S[self] <- S[self] - M                       # a unit shares a label with itself in every column
  den <- matrix(st$n_k, N, k, byrow = TRUE)
  den[self] <- den[self] - 1
  A <- S / (M * pmax(den, 1))
  A[den == 0] <- 0
  A
}

#' Per-unit loss gap: the change in salso's generalised-VI loss (summed over
#' units, bits) if unit i alone moved from its region r to region s:
#'   delta_is = (2 - a) \\[f(n_s + 1) - f(n_s) + f(n_r - 1) - f(n_r)\\]
#'              - (2 / M) \\[Dp\\[i, s\\] + Dm\\[i, r\\]\\]
#' The draw-only term a * sum f(s_tp) does not depend on the labels and cancels.
#' gap_i = min over s != r. gap > 0: no single move improves the loss;
#' gap < 0: moving i to `alt` would lower it. The per-unit expected loss
#' changes by gap / N.
#' Needs `st` from .ensemble_stats(..., gvi_a = a) on the same labels.
#' @keywords internal
#' @noRd
.vi_move_gap <- function(st, hard) {
  stopifnot(!is.null(st$gvi))
  a <- st$gvi$a; M <- st$M; n_k <- st$n_k
  N <- length(hard); k <- length(n_k)
  own      <- cbind(seq_len(N), hard)
  join_k   <- .xlog2x(n_k + 1) - .xlog2x(n_k)              # entering region s
  leave_k  <- .xlog2x(n_k - 1) - .xlog2x(n_k)              # leaving region r
  delta <- (2 - a) * (matrix(join_k, N, k, byrow = TRUE) + leave_k[hard]) -
    (2 / M) * (st$gvi$Dp + st$gvi$Dm[own])
  delta[own] <- Inf
  alt <- max.col(-delta, ties.method = "first")
  list(gap = delta[cbind(seq_len(N), alt)], alt = alt)
}

#' Per-unit summaries of a membership matrix against the final labels.
#'
#' certainty = U\\[i, hard_i\\]; margin = best minus second-best entry (the
#' quantity the transition-zone map thresholds); entropy = normalized Shannon
#' entropy. NA rows of U give NA in all three.
#' @keywords internal
#' @noRd
.membership_summaries <- function(U, hard) {
  N  <- nrow(U); k <- ncol(U)
  ok <- !is.na(U[, 1L])
  certainty <- margin <- rep(NA_real_, N)
  certainty[ok] <- U[cbind(which(ok), hard[ok])]
  if (any(ok)) {
    Uo   <- U[ok, , drop = FALSE]
    top  <- max.col(Uo, ties.method = "first")
    rest <- Uo; rest[cbind(seq_len(nrow(Uo)), top)] <- -Inf
    margin[ok] <- Uo[cbind(seq_len(nrow(Uo)), top)] -
      rest[cbind(seq_len(nrow(Uo)), max.col(rest, ties.method = "first"))]
  }
  list(certainty = certainty, margin = margin, entropy = .row_entropy(U))
}

#' Greedily merge the most co-associated pair of regions until k remain.
#' Uses only the k0 x k0 pair-count matrix, so merges are O(k0^2) each.
#' @keywords internal
#' @noRd
.merge_to_k <- function(memb, lab, k) {
  k0 <- max(lab)
  st <- .ensemble_stats(memb, lab, k0)
  P  <- st$P; n <- st$n_k; M <- st$M
  groups <- as.list(seq_len(k0))
  while (length(n) > k) {
    B <- P / (M * outer(n, n)); diag(B) <- -Inf
    ij <- arrayInd(which.max(B), dim(B)); i <- min(ij); j <- max(ij)
    P[i, ] <- P[i, ] + P[j, ]
    P[, i] <- P[, i] + P[, j]
    P <- P[-j, -j, drop = FALSE]
    n[i] <- n[i] + n[j]; n <- n[-j]
    groups[[i]] <- c(groups[[i]], groups[[j]]); groups[[j]] <- NULL
  }
  new <- integer(length(lab))
  for (g in seq_along(groups)) new[lab %in% groups[[g]]] <- g
  new
}


# -----------------------------------------------------------------------------
# 3. Hard consensus partition
# -----------------------------------------------------------------------------

#' Signal a configuration-level failure
#'
#' Class `pulseR_config_error`: the inputs are valid but this configuration
#' cannot be fitted (e.g. the ensemble does not support k regions).
#' tune_regions() logs these as failed rows; every other error still aborts.
#' @keywords internal
#' @noRd
.config_error <- function(msg, parent = NULL) {
  stop(structure(class = c("pulseR_config_error", "error", "condition"),
                 list(message = msg, call = NULL, parent = parent)))
}

#' salso at exactly k regions (merging down if it overshoots).
#' Returns list(labels, step, salso): `step` codes index .LABEL_STEPS;
#' `salso` keeps the search parameter a that salso settled on and salso's own
#' expected loss for its unmerged estimate.
#' @keywords internal
#' @noRd
.salso_k <- function(memb, k, n_runs = 16L, n_cores = 1L) {
  draws <- t(memb); storage.mode(draws) <- "integer"
  est <- tryCatch(
    salso::salso(draws,
                 loss   = salso::VI(a = list(nClusters = k)),
                 nRuns  = n_runs, nCores = n_cores),
    error = function(e) .config_error(sprintf(paste0(
      "salso could not find a consensus partition with k = %d ",
      "regions (salso: %s). The ensemble likely cannot support ",
      "k = %d; lower k or raise regions_per_partition."),
      k, conditionMessage(e), k), parent = e))
  info <- attr(est, "info")
  raw  <- as.integer(est); rm(est)                    # est carries a copy of draws
  # salso's search for a solves a step function, and its final rerun at the
  # root can land on the low side (k - 1). Step a down until >= k regions.
  a_fit <- if (is.null(info$a)) NA_real_ else as.numeric(info$a)
  if (length(a_fit) != 1L || !is.finite(a_fit))
    stop(sprintf(paste0(
      "salso %s did not report the search parameter `a` ",
      "(attr(est, \"info\")$a). pulseR needs it for the loss terms and ",
      "the sampled path. Update salso: install.packages(\"salso\")."),
      as.character(utils::packageVersion("salso"))), call. = FALSE)
  n_retry <- 0L
  while (length(unique(raw)) < k && is.finite(a_fit) && a_fit > 0 && n_retry < 10L) {
    n_retry <- n_retry + 1L
    a_fit <- max(0, a_fit - 0.02)
    est  <- suppressWarnings(salso::salso(draws, loss = salso::VI(a = a_fit),
                                          nRuns = n_runs, nCores = n_cores))
    info <- attr(est, "info")
    raw  <- as.integer(est); rm(est)
  }
  raw  <- match(raw, sort(unique(raw)))
  k0   <- max(raw)
  if (k0 < k)
    .config_error(sprintf(paste0(
      "salso returned %d regions for the requested k = %d. The ensemble ",
      "does not support that many consensus regions; lower k or raise ",
      "regions_per_partition."), k0, k))
  lab  <- raw
  step <- rep(1L, length(raw))                        # "salso"
  if (k0 > k) {
    lab   <- .merge_to_k(memb, raw, k)
    n_raw <- tabulate(raw, k0)
    keep  <- vapply(seq_len(k), function(g) {
      r <- unique(raw[lab == g]); r[which.max(n_raw[r])]
    }, integer(1))
    step[!(raw %in% keep)] <- 2L                # absorbed by "merge"
  }
  list(labels = lab, step = step,
       salso  = list(a             = a_fit,
                     expected_loss = as.numeric(info$expectedLoss),
                     a_retries     = n_retry,
                     k_raw         = k0,
                     n_runs        = as.integer(info$nRuns)))
}

#' Single-unit moves on salso's generalised-VI loss
#'
#' Sequential local search: visit each candidate unit, move it to the region
#' with the most negative exact loss change (the delta of .vi_move_gap()), and
#' update the counts before the next unit. Every accepted move lowers the
#' loss, so the loss decreases monotonically. Moves that would empty a region
#' are not allowed, so k is preserved. Each round starts from the full gap
#' vector, so the search stops at a partition where no allowed single move
#' improves the loss (converged = TRUE) or after `max_rounds` rounds.
#' `blocked` counts units whose only improving move would empty their region:
#' a sign that the loss at this a prefers fewer than k regions.
#'
#' Counts are held as an (R*M) x k integer matrix T, row p + R (t - 1) for
#' piece p of draw t, the same layout as .ensemble_stats(). Move deltas use
#' lookup tables of f(c + 1) - f(c) and f(c - 1) - f(c), f(x) = x log2 x.
#' @keywords internal
#' @noRd
.vi_moves <- function(memb, hard, k, a, max_rounds = 10L, tol = 1e-9,
                      verbose = FALSE) {
  N <- nrow(memb); M <- ncol(memb); R <- as.integer(max(memb))
  if (!is.finite(a)) stop("`a` must be finite for VI moves.", call. = FALSE)
  off <- R * (seq_len(M) - 1L)
  Tn  <- matrix(0L, R * M, k)
  lab0 <- R * (hard - 1L)
  for (t in seq_len(M)) {
    rows <- off[t] + seq_len(R)
    Tn[rows, ] <- matrix(tabulate(memb[, t] + lab0, R * k), R, k)
  }
  n_k  <- tabulate(hard, k)
  cnt  <- 0:(N + 1L)
  F_in  <- .xlog2x(cnt + 1) - .xlog2x(cnt)        # index c + 1
  F_out <- c(0, .xlog2x(cnt[-1L] - 1) - .xlog2x(cnt[-1L]))
  log <- data.frame(round = integer(0), candidates = integer(0), moves = integer(0))
  converged <- FALSE
  n_blocked <- 0L
  
  for (rd in seq_len(max_rounds)) {
    gp   <- .vi_move_gap(.ensemble_stats(memb, hard, k, gvi_a = a), hard)
    neg  <- gp$gap < -tol
    n_blocked <- sum(neg & n_k[hard] <= 1L)       # move would empty a region
    cand <- which(neg & n_k[hard] > 1L)
    if (!length(cand)) { converged <- TRUE; break }
    cand <- cand[order(gp$gap[cand])]               # most negative first
    n_mv <- 0L
    for (i in cand) {
      r <- hard[i]
      if (n_k[r] <= 1L) next
      ri <- memb[i, ] + off
      Tc <- Tn[ri, , drop = FALSE]
      gain  <- colSums(matrix(F_in[Tc + 1L], M, k))
      rem   <- sum(F_out[Tc[, r] + 1L])
      delta <- (2 - a) * (F_in[n_k + 1L] + F_out[n_k[r] + 1L]) -
        (2 / M) * (gain + rem)
      delta[r] <- Inf
      s <- which.min(delta)
      if (delta[s] < -tol) {
        Tn[cbind(ri, r)] <- Tc[, r] - 1L
        Tn[cbind(ri, s)] <- Tc[, s] + 1L
        n_k[r] <- n_k[r] - 1L; n_k[s] <- n_k[s] + 1L
        hard[i] <- s
        n_mv <- n_mv + 1L
      }
    }
    log[nrow(log) + 1L, ] <- c(rd, length(cand), n_mv)
    if (verbose) message(sprintf("  VI moves, round %d: %d candidate(s), %d move(s)",
                                 rd, length(cand), n_mv))
    if (n_mv == 0L) break                             # only blocked moves remain
  }
  list(labels = hard, converged = converged, blocked = n_blocked, rounds = log)
}

#' Full-data salso below the threshold; sample + cohesion assignment above it,
#' then single-unit moves on the same loss salso minimised (VI at salso's a).
#' Returns the .salso_k() list, with `step` covering all N units and, on the
#' sampled path, `moves` (converged flag and per-round log).
#' @keywords internal
#' @noRd
.consensus_partition <- function(memb, k, large_n_threshold, sample_size,
                                 refine_iter, n_runs, n_cores, verbose) {
  N <- nrow(memb)
  if (N <= large_n_threshold)
    return(.salso_k(memb, k, n_runs, n_cores))
  
  n_s   <- min(as.integer(sample_size), N)
  idx   <- sort(sample.int(N, n_s))
  fit_s <- .salso_k(memb[idx, , drop = FALSE], k, n_runs, n_cores)
  lab_s <- fit_s$labels
  
  st   <- .ensemble_stats(memb, lab_s, k, ref_idx = idx)
  hard <- max.col(.cohesion_from_stats(st, idx, lab_s, N), ties.method = "first")
  hard[idx] <- lab_s
  step <- rep(3L, N)                                   # "assign"
  step[idx] <- fit_s$step
  
  mv <- .vi_moves(memb, hard, k, a = fit_s$salso$a,
                  max_rounds = refine_iter, verbose = verbose)
  step[mv$labels != hard] <- 4L                        # "refine"
  if (verbose && !mv$converged)
    message(sprintf(paste0("  VI moves stopped after %d round(s) with improving ",
                           "moves left; raise `refine_iter`."), nrow(mv$rounds)))
  if (mv$blocked > 0L)
    warning(sprintf(paste0("%d unit(s) sit alone in a region that the VI loss would ",
                           "rather dissolve: at a = %.3f the ensemble supports fewer ",
                           "than k = %d regions."), mv$blocked, fit_s$salso$a, k),
            call. = FALSE)
  list(labels = mv$labels, step = step, salso = fit_s$salso,
       moves = list(converged = mv$converged, blocked = mv$blocked, rounds = mv$rounds))
}

#' Keep, for each region, the connected component carrying the most score
#' mass; grow the freed units back onto adjacent regions in order of score.
#' `score` is any N x k matrix (compute_consensus() passes Pi of the
#' unrepaired labels). NA entries count as 0.
#' @keywords internal
#' @noRd
.enforce_contiguity <- function(graph, hard, score) {
  aff <- score; aff[is.na(aff)] <- 0
  n <- length(hard); k <- ncol(aff)
  ep  <- igraph::ends(graph, igraph::E(graph), names = FALSE)
  lab <- rep(NA_integer_, n)
  
  for (r in seq_len(k)) {
    in_r <- which(hard == r)
    if (!length(in_r)) next
    cmp  <- igraph::components(igraph::induced_subgraph(graph, in_r))$membership
    mass <- rowsum(aff[in_r, r], cmp)
    core <- as.integer(rownames(mass)[which.max(mass[, 1])])
    lab[in_r[cmp == core]] <- r
  }
  
  a <- ep[, 1]; b <- ep[, 2]
  repeat {
    free <- is.na(lab)
    if (!any(free)) break
    f1 <- free[a] & !free[b]
    f2 <- free[b] & !free[a]
    node <- c(a[f1], b[f2])
    if (!length(node)) { lab[free] <- hard[free]; break }   # isolated fragments
    reg  <- c(lab[b[f1]], lab[a[f2]])
    o    <- order(-aff[cbind(node, reg)])
    node <- node[o]; reg <- reg[o]
    keep <- !duplicated(node)
    lab[node[keep]] <- reg[keep]
  }
  lab
}


# ----- 
#' Construct an "pulse_ensemble" object
#'
#' Used both by compute_ensemble() and by the tuner when it slices
#' a cached ensemble down to a smaller N_tree, so that a sliced ensemble is
#' stamped with its own column count and chance co-association rather than
#' those of the larger matrix it came from.
#'
#' @keywords internal
#' @noRd
.new_ensemble <- function(memb, N_tree, regions_per_partition,
                          weight_exponent,
                          prior_strength, n_vertices,
                          p0_by_col = NULL,
                          cut_rule = NA_character_, min_size = NA_integer_) {
  if (is.null(p0_by_col)) p0_by_col <- .chance_coassoc_by_col(memb)
  structure(
    list(
      memberships           = memb,
      N_tree                  = as.integer(N_tree),
      regions_per_partition  = as.integer(regions_per_partition),
      n_members             = ncol(memb),
      weight_exponent       = as.numeric(weight_exponent),
      n_vertices            = as.integer(n_vertices),   # vcount() is double in igraph >= 2.0
      prior_strength        = prior_strength,
      chance_coassoc        = p0_by_col,
      cut_rule              = cut_rule,
      min_size              = as.integer(min_size)
    ),
    class = "pulse_ensemble"
  )
}

#' @rdname compute_ensemble
#' @param x An object of class \code{"pulse_ensemble"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.pulse_ensemble} returns \code{x} invisibly.
#' @export
print.pulse_ensemble <- function(x, ...) {
  cat(sprintf(
    "<pulse_ensemble>  %d vertices x %d trees | regions_per_partition = %d\n",
    x$n_vertices, x$N_tree, x$regions_per_partition
  ))
  invisible(x)
}


#' @rdname compute_consensus
#' @param x An object of class \code{"pulse_consensus"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.pulse_consensus} returns \code{x} invisibly.
#' @export
print.pulse_consensus <- function(x, ...) {
  cat(sprintf(
    "<pulse_consensus>  %d units x %d regions | avg silhouette = %.4f\n",
    nrow(x$memberships), x$k, x$avg_sil_width
  ))
  invisible(x)
}


#' @rdname sample_spanning_trees
#' @param x An object to print (e.g. of class \code{"pulse_trees"}).
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.pulse_trees} returns \code{x} invisibly.
#' @export
print.pulse_trees <- function(x, ...) {
  cat(sprintf(
    "<pulse_trees>  %d trees | %d vertices | seed = %s\n",
    x$n_trees, x$n_vertices,
    if (is.null(x$seed)) "NULL" else as.character(x$seed)
  ))
  if (!is.null(x$prior_labels)) {
    n_cov    <- sum(!is.na(x$prior_labels))
    n_within <- sum(x$within_prior_edge, na.rm = TRUE)
    cat(sprintf(
      "  prior typology: %d/%d vertices covered | %d/%d edges within-prior\n",
      n_cov, x$n_vertices, n_within, length(x$within_prior_edge)
    ))
  }
  invisible(x)
}



#' Environment for the validity indices and the null baseline
#'
#' \code{NULL} falls back to the matrix \code{add_edge_weight()} stored on the
#' graph, which is the space the edge weights were computed in. A supplied
#' \code{env} is used as given: \code{env_id} and non-numeric columns are
#' dropped, nothing is scaled.
#'
#' @return Numeric matrix with one row per vertex, or \code{NULL}.
#' @keywords internal
#' @noRd
.resolve_env <- function(env, env_id, graph, verbose = FALSE) {
  if (is.null(env)) {
    env <- igraph::graph_attr(graph, "env")
    if (is.null(env)) {
      if (verbose)
        message("Note: no `env` supplied and none stored on `graph`; ",
                "CH_FIT and R2_FIT are NA.")
      return(NULL)
    }
    if (verbose)
      message("Using the environment matrix stored on `graph` by ",
              "add_edge_weight() for CH_FIT and R2_FIT.")
    env <- as.matrix(env)
  } else {
    env <- as.data.frame(env)
    if (!is.null(env_id)) env <- env[, setdiff(names(env), env_id), drop = FALSE]
    env <- as.matrix(env[, vapply(env, is.numeric, logical(1)), drop = FALSE])
    .warn_id_like(env, "`env_id`")
  }
  if (nrow(env) != igraph::vcount(graph))
    stop("`env` must have one row per vertex of `graph`.", call. = FALSE)
  env
}

#' Canonical tuning axes, in coordinate-descent sweep order
#'
#' Order is deliberate: ensemble-shape parameters that change *which* partitions
#' the ensemble contains are swept before parameters that change *how many*, and
#' \code{final_regions} is swept last because it is conditionally cheap (it
#' reuses the cached ensemble) so the reported optimum has a \emph{k} that is
#' optimal for the selected ensemble.
#'
#' @keywords internal
#' @noRd
.TUNE_AXES <- c("regions_per_partition", "weight_exponent",
                "N_tree", "final_regions")

#' Ensemble-level axes (those that require recomputing the membership matrix)
#' @keywords internal
#' @noRd
.ENSEMBLE_AXES <- c("regions_per_partition", "weight_exponent")

#' Row of the best score, with an explicit tie-break
#'
#' Used by the final selection and by the coordinate-descent sweeps, so every
#' strategy breaks ties the same way.
#' @param log Data frame with the four axis columns.
#' @param score Score per row, higher is better.
#' @return Row index, or \code{NA_integer_} if every score is \code{NA}.
#' @keywords internal
#' @noRd
.pick_best <- function(log, score = log$score) {
  if (all(is.na(score))) return(NA_integer_)
  top <- which(score >= max(score, na.rm = TRUE) - .TIE_TOL)
  if (length(top) == 1L) return(top)
  top[order(log$final_regions[top], log$regions_per_partition[top],
            log$weight_exponent[top], -log$N_tree[top])[1L]]
}
.TIE_TOL  <- 1e-12
.TIE_RULE <- paste("fewest final_regions, then fewest regions_per_partition,",
                   "then smallest weight_exponent, then largest N_tree")



#' Canonical key for one hyper-parameter configuration
#'
#' Numeric axes are formatted at full double precision so that floating-point
#' candidates round-trip exactly through the log.
#'
#' @param df Data frame with the four axis columns.
#' @return Character vector of length \code{nrow(df)}.
#' @keywords internal
#' @noRd
.tune_key <- function(df) {
  paste(as.integer(df$regions_per_partition),
        sprintf("%.17g", as.numeric(df$weight_exponent)),
        as.integer(df$N_tree),
        as.integer(df$final_regions),
        sep = "|")
}

#' Tune SKATER-CON Hyperparameters
#'
#' Searches over up to four hyper-parameters -- \code{regions_per_partition},
#' \code{weight_exponent}, \code{N_tree} and \code{final_regions} -- to find the combination optimising a
#' clustering-validity criterion, operating on a pre-computed
#' \code{"pulse_trees"} object. Any argument supplied as a length-1 vector is
#' held fixed; any argument supplied as a vector of length > 1 becomes a search
#' axis.
#'
#' @section Cost structure:
#' Only \code{final_regions} is cheap to vary: it reuses the cached ensemble
#' membership matrix. \code{N_tree} is nearly free, because an ensemble of
#' \code{N_tree} trees is an exact column-subset of an ensemble of more trees
#' drawn from the same \code{trees} object; the tuner therefore builds each
#' ensemble once at \code{max(N_tree)} and slices columns. The remaining two
#' axes each require a fresh pass of
#' \code{\link{compute_ensemble}}, so the number of distinct
#' \code{(regions_per_partition, weight_exponent)} pairs
#' is the real cost driver, not the number of grid cells.
#'
#' Because ensembles are built at \code{max(N_tree)} regardless of the incumbent,
#' supplying a vector \code{N_tree} raises the cost of \emph{every} ensemble
#' build. This buys determinism (results do not depend on the order in which
#' configurations happen to be visited) and paired comparisons across
#' \code{N_tree} (all ensemble sizes share the same trees and the same cut
#' realisations). When \code{N_tree} is scalar there is no penalty.
#'
#' @section Caching:
#' A single-slot cache holds one ensemble matrix at a time, and evaluation is
#' ordered so that configurations sharing an ensemble are visited consecutively.
#' A keyed cache over a four-way grid would hold one N x N_tree integer
#' matrix per ensemble pair, which at continental N is tens of gigabytes.
#' Peak memory here is one full matrix plus at most one column slice.
#'
#' @section Search strategies:
#' \code{"grid"} evaluates the full factorial. \code{"sequential"} runs one
#' coordinate-descent pass over the active axes in the order
#' \code{regions_per_partition}, \code{weight_exponent}, \code{N_tree},
#' \code{final_regions}, starting from
#' the median candidate on each axis. \code{"iterative"} repeats that pass until
#' the incumbent configuration stops moving or \code{max_iter} is reached.
#'
#' @section Ties:
#' Several configurations can share the best score, most often under Borda,
#' whose scores are mean ranks. Ties are broken, in every search strategy, by
#' the fewest \code{final_regions}, then the fewest
#' \code{regions_per_partition}, then the smallest \code{weight_exponent},
#' then the largest \code{N_tree}. The number of tied rows and the runner-up
#' are returned.
#'
#' @section Caveat on the ensemble-size axis:
#' Most validity indices drift monotonically with \code{N_tree} (Monte Carlo
#' noise falls and co-association resolution rises as 1 / N_tree) rather than
#' showing an interior optimum, so a CVI-selected \code{N_tree} will typically
#' just be the largest candidate. Choose it by convergence, not by tuning.
#'
#' @section Validity indices:
#' \code{CH_FIT} and \code{R2_FIT}: pseudo-F and \eqn{R^2} of the hard labels
#' on \code{env}. By default \code{env} is the matrix the edge weights were
#' computed from, so these are internal fit indices: they measure how well
#' the regions summarise the space the trees were cut in, not agreement with
#' independent data. Every column adds to the sums of squares, so \eqn{p}
#' near-duplicate predictors carry about \eqn{p} times the weight of one
#' independent predictor; see \code{\link{add_edge_weight}}. The former
#' names \code{CH_ENV}, \code{R2_ENV}, \code{"ch_env"} and \code{"r2_env"}
#' are still accepted in \code{voters} and \code{tuning_metric}, with a
#' warning. \code{STAB}: split-half ensemble agreement. The ensemble
#' columns are split at random into two disjoint halves, the consensus is
#' computed at the same \code{final_regions} on each half, and the ARI
#' between the two partitions is averaged over \code{stability_B} splits.
#' The halves share no trees and no cuts, so \code{STAB} asks whether new
#' trees and cuts would give the same regions. Each half has half the
#' columns, so the value is conservative for the full ensemble, and it rises
#' with \code{N_tree}. \code{SIL_CA}, \code{CERT}, \code{EVI} and \code{PAC}
#' are computed from the ensemble itself and are reported, not used by the
#' default Borda vote.
#'
#' @param graph The parent \code{igraph} object the trees were sampled from.
#' @param env Optional data frame or matrix, one row per vertex of
#'   \code{graph}, for the validity indices \code{CH_FIT} and \code{R2_FIT}.
#'   \code{NULL} (default) uses the matrix \code{\link{add_edge_weight}}
#'   stored on \code{graph} (graph attribute \code{env}), so the indices are
#'   computed in the space of the edge weights. A supplied \code{env} is used
#'   as given: non-numeric columns are dropped and nothing is scaled. Columns
#'   that look like IDs trigger a warning.
#' @param env_id Character vector of columns in \code{env} to exclude, for
#'   example numeric ID columns. Default \code{NULL}.
#' @param trees A \code{"pulse_trees"} object from
#'   \code{\link{sample_spanning_trees}}.
#' @param N_tree Integer or integer vector (each \eqn{\ge 2} and
#'   \eqn{\le} \code{trees$n_trees}). Candidate ensemble size(s).
#' @param regions_per_partition Integer or integer vector (each \eqn{\ge 2}).
#'   Candidate number(s) of intermediate regions per cut-set. Default 50.
#' @param final_regions Integer or integer vector (each \eqn{\ge 2}). Candidate
#'   number(s) of final consensus clusters \emph{k}.
#' @param weight_exponent Numeric or numeric vector (each \eqn{\ge 0}).
#'   Candidate cut-probability exponent(s) \eqn{\gamma}. Default 1.
#' @param tuning_metric Character. Selection criterion, higher is better:
#'   \code{"borda"} (default; mean rank over \code{voters}),
#'   \code{"ch_fit"} (\code{CH_FIT}, pseudo-F of the hard labels on
#'   \code{env}), \code{"r2_fit"} (\code{R2_FIT}), \code{"stability"}
#'   (\code{STAB}), \code{"silhouette_ca"} (\code{SIL_CA}, co-association
#'   silhouette) or \code{"expected_vi"} (minus \code{EVI}, mean variation of
#'   information between ensemble columns and the consensus).
#' @param tuning_strategy Character. \code{"grid"} (default), \code{"sequential"}, or
#'   \code{"iterative"}. See Details.
#' @param stability_B Integer. Random splits for \code{STAB} (see
#'   Validity indices). Each split costs two consensus runs at half the
#'   ensemble size. 0 skips it. Default 10.
#' @param voters Character. Validity columns that vote in the Borda aggregate
#'   (\code{tuning_metric = "borda"}). Default \code{c("CH_FIT", "STAB")}. Check
#'   \code{STAB} against \code{\link{compute_null_baseline}} before relying on it.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda}. Lowers the
#'   log cutting propensity of within-prior edges by \eqn{\lambda} standard
#'   deviations of the environmental term \eqn{\gamma \log w}, computed per
#'   tree: the same unit as \code{size_influence}, so a value means the same
#'   on any graph. \eqn{\lambda = 0} (default) ignores the prior;
#'   \eqn{\lambda = 1} makes a within-prior edge as hard to cut as an edge
#'   one standard deviation lighter on the log scale. If that spread is zero
#'   (\code{weight_exponent = 0} or constant weights), one unit halves the
#'   propensity instead, with a warning.
#' @param prior_typology Optional. Vertex attribute name or label vector; used
#'   only by the \code{prior_strength = Inf} short-circuit.
#' @param max_iter Integer. Maximum passes for \code{"iterative"}. Default 10.
#' @param cut_rule,refine_boundary,min_size Cut rule, boundary refinement and
#'   smallest piece, as in \code{\link{get_regions}}. Fixed for the whole
#'   search; every configuration is scored on the refined partition when
#'   \code{refine_boundary = TRUE}, so the score belongs to the returned regions.
#' @param max_configs Integer or \code{Inf}. Guard against an accidentally
#'   enormous search: the call fails before any work if the number of distinct
#'   ensemble pairs exceeds this. Default 64.
#' @param size_influence Non-negative numeric, passed to
#'   \code{\link{compute_ensemble}}. Not a tuning axis. Default 2.
#' @param large_n_threshold Integer. Above this number of vertices, salso runs
#'   on \code{sample_size} units and the rest are assigned and refined.
#'   Default 20000.
#' @param sample_size Integer. Units passed to salso above
#'   \code{large_n_threshold}. Default 5000.
#' @param n_cores Integer. Cores passed to \code{salso::salso()}. Default 1.
#' @param seed Optional integer seed; the caller's RNG state is restored on
#'   exit.
#' @param stability_seed Optional integer. Seeds the subsampling inside the
#'   \code{STAB} index, identically for every configuration, so configurations
#'   are compared on the same column splits. Defaults to \code{seed}.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return An object of class \code{"pulse_regions_tuning"}: a list with
#'   \code{best_result} (a \code{"pulse_consensus"} object), \code{best_config}
#'   (one-row data frame of the winning values on all four axes), the individual
#'   slots \code{best_regions_per_partition}, \code{best_final_regions},
#'   \code{best_weight_exponent},
#'   \code{best_N_tree}, plus \code{best_score}, \code{tuning_log} (one row per
#'   configuration, with all four axis columns and all CVI columns),
#'   \code{candidates}, \code{active_axes}, \code{n_ensembles_built},
#'   \code{tuning_metric}, \code{voters}, \code{tuning_strategy}, \code{prior_strength} and
#'   \code{ensemble} (the winner's \eqn{N \times N_{tree}} label matrix).
#'   \code{n_tied} counts configurations sharing the best score,
#'   \code{runner_up} is the best other row of \code{tuning_log}, and
#'   \code{tie_rule} states how ties are broken (see Details).
#'   When \code{prior_strength = Inf}, the object has the same class and slots,
#'   with the prior as \code{best_result} (one-hot \code{memberships},
#'   \code{hard_labels}, \code{coverage}) and no \code{tuning_log}.
#'
#' @export
tune_regions <- function(graph, 
                         trees,
                         N_tree,
                         regions_per_partition = NULL,
                         final_regions,
                         weight_exponent     = 4,
                         tuning_metric     = c("borda", "ch_fit", "r2_fit", "stability", "silhouette_ca", "expected_vi"),
                         tuning_strategy   = c("grid", "sequential", "iterative"),
                         stability_B       = 10,
                         voters            = c("CH_FIT", "STAB"),
                         prior_strength    = 0,
                         prior_typology    = NULL,
                         max_iter          = 10,
                         size_influence    = 0,
                         cut_rule          = c("ssd", "propensity"),
                         refine_boundary   = TRUE,
                         min_size          = 5,
                         max_configs     = 64L,
                         large_n_threshold = 20000,
                         sample_size = 5000,
                         n_cores = 1,
                         env = NULL, 
                         env_id = NULL,
                         seed              = NULL,
                         stability_seed    = seed,
                         verbose           = TRUE) {
  
  # --- Checks -----------------------------------------------------*
  stopifnot(inherits(trees, "pulse_trees"))
  stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L,
            prior_strength >= 0)
  cut_rule <- match.arg(cut_rule)
  min_size <- as.integer(min_size)
  
  
  # Prepare env ------------------------------------------------------*
  env <- .resolve_env(env, env_id, graph, verbose)
  
  
  # --- Argument Matching  -----------------------------------------*
  tuning_strategy <- match.arg(tuning_strategy)
  tuning_metric   <- .renamed_cvi(tuning_metric, "tuning_metric")
  tuning_metric   <- match.arg(tuning_metric)
  voters   <- unique(.renamed_cvi(as.character(voters), "voters"))
  bad_v    <- setdiff(voters, names(.CVI_DIRECTION))
  if (length(bad_v) || !length(voters))
    stop("`voters` must be a non-empty subset of: ",
         paste(names(.CVI_DIRECTION), collapse = ", "), call. = FALSE)
  
  # --- Exact prior reproduction short-circuit --------------------------*
  if (is.infinite(prior_strength)) {
    if (is.null(prior_typology)) {
      stop("`prior_strength = Inf` requires `prior_typology`.", call. = FALSE)
    }
    prior_labels <- .resolve_prior_typology(graph, prior_typology)
    pr <- .exact_prior_result(prior_labels, verbose = verbose)
    K  <- ncol(pr$memberships)
    # Same class and slots as a tuned result, so get_local_types(),
    # as_mosaic_input() and crosswalk_types() accept it.
    best <- list(memberships = pr$memberships, hard_labels = pr$hard_labels,
                 coverage = pr$coverage, class_levels = pr$class_levels,
                 k = K)
    cfg <- data.frame(regions_per_partition = K, weight_exponent = NA_real_,
                      N_tree = NA_integer_, final_regions = K)
    return(structure(list(
      best_result = best, best_config = cfg,
      best_regions_per_partition = K, best_final_regions = K,
      best_weight_exponent = NA_real_, best_N_tree = NA_integer_,
      best_score = NA_real_, n_tied = NA_integer_, runner_up = NULL,
      tie_rule = NA_character_, tuning_log = NULL, candidates = NULL,
      active_axes = character(0), n_ensembles_built = 0L,
      tuning_metric = tuning_metric, voters = voters,
      tuning_strategy = tuning_strategy, prior_strength = Inf,
      cut_rule = NA_character_, refine_boundary = NA,
      min_size = NA_integer_, ensemble = NULL),
      class = "pulse_regions_tuning"))
  }
  
  # --- Normalize candidate sets ----------------------------------------*
  # weight_exponent stays double; the other four are counts.
  final_regions         <- sort(unique(as.integer(final_regions)))
  N_tree                <- sort(unique(as.integer(N_tree)))
  weight_exponent       <- sort(unique(as.numeric(weight_exponent)))
  rpp_tied <- is.null(regions_per_partition)
  regions_per_partition <- 
    if (rpp_tied){
      3 * final_regions
    } else {
      sort(unique(as.integer(regions_per_partition)))
    }
  
  
  # --- Checks -----------------------------------------------------*
  stopifnot(
    length(regions_per_partition) >= 1L,
    all(regions_per_partition >= 2),
    all(regions_per_partition < trees$n_vertices)
  )
  stopifnot(length(final_regions) >= 1L, all(final_regions >= 2))
  stopifnot(length(N_tree) >= 1L, all(N_tree >= 2), all(N_tree <= trees$n_trees))
  stopifnot(length(weight_exponent) >= 1L, all(is.finite(weight_exponent)),
            all(weight_exponent >= 0))
  stopifnot(is.numeric(max_configs), length(max_configs) == 1L,
            max_configs >= 1)
  
  cand <- list(
    regions_per_partition = regions_per_partition,
    weight_exponent      = weight_exponent,
    N_tree                 = N_tree,
    final_regions        = final_regions
  )
  active_axes <- .TUNE_AXES[vapply(cand[.TUNE_AXES], length, integer(1)) > 1L]
  if (rpp_tied) active_axes <- setdiff(active_axes, "regions_per_partition")      
  
  
  # Distinct ensemble pairs: the actual cost of the search. Tied: one ensemble per (final region, weight exponent)
  n_ens_triples <- length(regions_per_partition) * length(weight_exponent)
  if (tuning_strategy == "grid" && n_ens_triples > max_configs) {
    stop(sprintf(paste0(
      "Grid requires %d distinct ensemble computations (limit `max_configs` = %s). ",
      "Use tuning_strategy = \"sequential\" / \"iterative\", shrink the ",
      "regions_per_partition / weight_exponent ",
      "candidate sets, or raise `max_configs` deliberately."),
      n_ens_triples, format(max_configs)), call. = FALSE)
  }
  
  if (verbose && "N_tree" %in% active_axes) {
    message("Note: most validity indices drift monotonically with ensemble ",
            "size; every ensemble will be built at N_tree = ",
            max(N_tree), ".")
  }
  
  # --- RNG handling ----------------------------------------------------*
  .local_seed(seed)
  # One seed for every ensemble build (common random numbers for the cut
  # race): configurations differ only by their settings, and rebuilding the
  # winner's ensemble after a cache eviction reproduces it exactly.
  ens_seed <- sample.int(.Machine$integer.max, 1L)
  build_ensemble <- function(sr, we) {
    .local_seed(ens_seed)
    compute_ensemble(
      trees                 = trees,
      N_tree                = N_tree_build,
      graph                 = graph,
      regions_per_partition = sr,
      prior_strength        = prior_strength,
      weight_exponent       = we,
      size_influence        = size_influence,
      cut_rule              = cut_rule,
      min_size              = min_size,
      verbose               = verbose)
  }
  # Scored partitions, keyed like tuning_log rows (labels + raw step codes).
  fit_store <- new.env(parent = emptyenv())
  
  # --- Single-slot ensemble cache --------------------------------------*
  # Key is the ensemble pair only. N_tree is served by slicing columns of a
  # matrix built once at max(N_tree): column j comes from tree j, so the
  # first N_tree columns are exactly the ensemble N_tree trees would give.
  N_tree_build <- max(N_tree)
  ens_cache  <- new.env(parent = emptyenv())
  ens_cache$key   <- NULL
  ens_cache$mat   <- NULL
  ens_cache$p0c   <- NULL
  ens_cache$built <- 0L
  
  get_memb <- function(sr, we, ns) {
    
    key <- paste(sr, sprintf("%.17g", we), sep = "|")
    
    if (!identical(key, ens_cache$key)) {
      ens <- build_ensemble(sr, we)
      # Drop the previous matrix before holding two at once.
      ens_cache$mat   <- NULL
      ens_cache$mat   <- ens$memberships
      ens_cache$p0c   <- ens$chance_coassoc   # per-column vector here
      ens_cache$key   <- key
      ens_cache$built <- ens_cache$built + 1L
    }
    n_col <- ns
    mat <- if (n_col >= ncol(ens_cache$mat)) {
      ens_cache$mat
    } else {
      ens_cache$mat[, seq_len(n_col), drop = FALSE]
    }
    .new_ensemble(
      memb                 = mat,
      N_tree                 = ns,
      regions_per_partition = sr,
      weight_exponent      = we,
      prior_strength       = prior_strength,
      n_vertices           = nrow(mat),
      p0_by_col            = ens_cache$p0c[seq_len(n_col)],
      cut_rule             = cut_rule,
      min_size             = min_size
    )
  }
  
  # Canonical column sets (kept fixed so rbind() works across batches).
  axis_cols <- .TUNE_AXES
  diag_cols <- c("P0", "N_EFF")
  cvi_cols  <- c("CH_FIT", "R2_FIT", "STAB", "STAB_FAIL","SIL_CA", "CERT", "EVI", "PAC")
  col_order <- c(axis_cols, diag_cols, cvi_cols, "error")
  
  # A configuration that raised pulseR_config_error: NA diagnostics and CVIs.
  failed_row <- function(cfg, msg) {
    na <- as.data.frame(as.list(stats::setNames(
      rep(NA_real_, length(diag_cols) + length(cvi_cols)),
      c(diag_cols, cvi_cols))))
    cbind(cfg[, axis_cols, drop = FALSE], na, error = msg,
          stringsAsFactors = FALSE)
  }
  
  # Failed rows get score NA and stay out of the Borda ranking.
  score_rows <- function(log) {
    s  <- rep(NA_real_, nrow(log))
    ok <- is.na(log$error)
    if (any(ok)) s[ok] <- .derive_score_column(log[ok, , drop = FALSE],
                                               tuning_metric, voters)
    s
  }
  
  # --- Score one configuration -> one row of CVI columns ---------------*
  score_config <- function(cfg) {
    memb_mat <- get_memb(cfg$regions_per_partition, cfg$weight_exponent,
                         cfg$N_tree)
    meta <- .ensemble_meta(memb_mat)
    fit  <- .consensus_partition(meta$memb_mat, cfg$final_regions,
                                 large_n_threshold, sample_size,
                                 refine_iter = 10L, n_runs = 16L,
                                 n_cores = n_cores, verbose = FALSE)
    res  <- .finish_consensus(meta, fit, cfg$final_regions, graph = graph,
                              large_n_threshold = large_n_threshold,
                              sample_size = sample_size, n_cores = n_cores,
                              refine_boundary = refine_boundary,
                              min_size = min_size,
                              report_native = FALSE, verbose = FALSE)
    fit$step <- as.raw(fit$step)
    cvis <- .consensus_cvis(
      res, memb = memb_mat$memberships, env = env,
      stability_B = stability_B, graph = graph,
      refine_iter = 10L, n_runs = 16L,
      stability_seed = stability_seed, large_n_threshold = large_n_threshold,
      sample_size = sample_size, n_cores = n_cores, verbose = verbose)
    list(out = cbind(cfg[, axis_cols, drop = FALSE],
                     data.frame(P0    = mean(memb_mat$chance_coassoc),
                                N_EFF = 1 / mean(memb_mat$chance_coassoc)),
                     cvis, error = NA_character_, stringsAsFactors = FALSE),
         fit = fit)
  }
  # --- tuning_log accumulator ------------------------------------------*
  tuning_log <- NULL
  
  rebuild_score_column <- function() {
    # Borda ranks are computed globally over all rows evaluated so far, so
    # the score column must be rebuilt after every batch.
    if (is.null(tuning_log) || nrow(tuning_log) == 0L) return(invisible(NULL))
    tuning_log$score <<- score_rows(tuning_log)
    invisible(NULL)
  }
  
  evaluate_configs <- function(cfg_df, tick = NULL) {
    cfg_df <- cfg_df[, axis_cols, drop = FALSE]
    if (rpp_tied) cfg_df$regions_per_partition <- 3 * cfg_df$final_regions
    seen   <- if (is.null(tuning_log)) character(0) else .tune_key(tuning_log)
    keys   <- .tune_key(cfg_df)
    new_rows <- vector("list", nrow(cfg_df))
    for (i in seq_len(nrow(cfg_df))) {
      if (keys[i] %in% seen) {          # already evaluated
        if (!is.null(tick)) tick(i)
        next
      }
      sc <- tryCatch(score_config(cfg_df[i, , drop = FALSE]),
                     pulseR_config_error = function(e) e)
      if (inherits(sc, "pulseR_config_error")) {
        new_rows[[i]] <- failed_row(cfg_df[i, , drop = FALSE],
                                    conditionMessage(sc))
        if (verbose) message("  Failed: ", conditionMessage(sc))
      } else {
        new_rows[[i]] <- sc$out
        assign(keys[i], sc$fit, envir = fit_store)
      }
      seen <- c(seen, keys[i])
      if (!is.null(tick)) tick(i)
    }
    new_rows <- new_rows[!vapply(new_rows, is.null, logical(1))]
    if (length(new_rows) > 0L) {
      new_df <- do.call(rbind, new_rows)
      new_df <- new_df[, col_order, drop = FALSE]
      if (is.null(tuning_log)) {
        tuning_log <<- new_df
      } else {
        base       <- tuning_log[, col_order, drop = FALSE]
        tuning_log <<- rbind(base, new_df)
      }
    }
    rebuild_score_column()
    invisible(NULL)
  }
  
  # --- Coordinate-descent helpers --------------------------------------*
  # Incumbent starts at the median candidate on every axis.
  cur <- lapply(cand, function(v) v[ceiling(length(v) / 2)])
  if (rpp_tied) cur$regions_per_partition <- 3 * cur$final_regions
  cur_df <- function() as.data.frame(cur[axis_cols], stringsAsFactors = FALSE)
  
  # Configurations that vary `axis` and hold every other axis at the incumbent.
  axis_configs <- function(axis) {
    vals <- cand[[axis]]
    df   <- cur_df()[rep(1L, length(vals)), , drop = FALSE]
    df[[axis]] <- vals
    rownames(df) <- NULL
    df
  }
  
  # Best value on `axis` among logged rows that match the incumbent elsewhere.
  best_along <- function(axis) {
    others <- setdiff(axis_cols, c(axis, if (rpp_tied) "regions_per_partition"))
    keep   <- rep(TRUE, nrow(tuning_log))
    for (o in others) keep <- keep & (tuning_log[[o]] == cur[[o]])
    idx <- which(keep)
    if (length(idx) == 0L) return(cur[[axis]])
    sc <- tuning_log$score[idx]
    if (all(is.na(sc))) return(cur[[axis]])
    tuning_log[[axis]][idx][.pick_best(tuning_log[idx, , drop = FALSE], sc)]
  }
  
  sweep_axis <- function(axis, label) {
    cfgs <- axis_configs(axis)
    # Descending N_tree means the max-size ensemble is built on the first
    # visit and every smaller candidate is a free column slice.
    if (axis == "N_tree") cfgs <- cfgs[order(-cfgs$N_tree), , drop = FALSE]
    if (verbose) message(sprintf("  Sweeping %s (%d candidates)",
                                 axis, nrow(cfgs)))
    tick <- if (verbose) .make_progress(nrow(cfgs), label) else NULL
    evaluate_configs(cfgs, tick)
    cur[[axis]] <<- best_along(axis)
    if (rpp_tied) cur$regions_per_partition <<- 3 * cur$final_regions
    if (verbose) message(sprintf("  >> Best %s = %s", axis,
                                 format(cur[[axis]])))
    invisible(NULL)
  }
  
  # --- Run the chosen tuning_strategy -----------------------------------------*
  if (verbose) {
    message(sprintf("Tuning (%s tuning_strategy, tuning_metric = %s) ...", tuning_strategy, tuning_metric))
    message(sprintf("  Active axes: %s",
                    if (length(active_axes)) paste(active_axes, collapse = ", ")
                    else "none (single configuration)"))
  }
  
  if (length(active_axes) == 0L) {
    # Degenerate case: nothing to search, evaluate the one configuration.
    evaluate_configs(cur_df())
    
  } else if (tuning_strategy == "grid") {
    # expand.grid varies its first argument fastest; ordering it this way
    # keeps configurations that share an ensemble contiguous, and visits
    # N_tree descending so the max-size build happens first in each block.
    pairs <- 
      if (rpp_tied){
        p <- expand.grid(
          N_tree = rev(N_tree), 
          final_regions = final_regions, 
          weight_exponent = weight_exponent,
          KEEP.OUT.ATTRS = FALSE, stringsAsFactors=FALSE
        )
        p$regions_per_partition <- 3 * p$final_regions
        p
      } else {
        expand.grid(
          final_regions        = final_regions,
          N_tree                 = rev(N_tree),
          regions_per_partition = regions_per_partition,
          weight_exponent      = weight_exponent,
          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
        )
      }
    
    if (verbose) message(sprintf(
      "  Grid: %d configurations across %d ensemble computations",
      nrow(pairs), n_ens_triples))
    tick <- if (verbose) .make_progress(nrow(pairs), "Grid search") else NULL
    evaluate_configs(pairs, tick)
    
  } else if (tuning_strategy == "sequential") {
    evaluate_configs(cur_df())             # seed the incumbent
    for (axis in active_axes) sweep_axis(axis, axis)
    
  } else if (tuning_strategy == "iterative") {
    evaluate_configs(cur_df())             # seed the incumbent
    iter <- 0L
    repeat {
      iter     <- iter + 1L
      prev_cur <- cur
      if (verbose) message(sprintf("\n  --- Pass %d/%d ---", iter, max_iter))
      for (axis in active_axes) {
        sweep_axis(axis, sprintf("Pass %d: %s", iter, axis))
      }
      # Convergence on the incumbent configuration, not on the score.
      # Under Borda the score column is re-ranked globally after every
      # batch, so a score comparison across passes is not stable.
      if (identical(cur, prev_cur)) {
        if (verbose) message("  Converged (incumbent unchanged).")
        break
      }
      if (iter >= max_iter) {
        if (verbose) message("  Reached max_iter.")
        break
      }
    }
  }
  
  # --- Deduplicate and pick the best configuration ---------------------*
  tuning_log <- tuning_log[!duplicated(.tune_key(tuning_log)), , drop = FALSE]
  rownames(tuning_log) <- NULL
  tuning_log$score <- score_rows(tuning_log)
  
  failed <- !is.na(tuning_log$error)
  if (all(failed))
    stop(sprintf("All %d configurations failed. First error: %s",
                 nrow(tuning_log), tuning_log$error[1L]), call. = FALSE)
  if (any(failed))
    warning(sprintf(paste0(
      "%d of %d configurations failed and were scored NA ",
      "(see tuning_log$error). First error: %s"),
      sum(failed), nrow(tuning_log), tuning_log$error[failed][1L]),
      call. = FALSE)
  
  if (all(is.na(tuning_log$score))) {
    stop("All hyper-parameter combinations produced an undefined score; ",
         "check the tuning_metric and inputs.", call. = FALSE)
  }
  
  best_idx <- .pick_best(tuning_log)
  best_cfg <- tuning_log[best_idx, axis_cols, drop = FALSE]
  rownames(best_cfg) <- NULL
  best_sc  <- tuning_log$score[best_idx]
  n_tied   <- sum(tuning_log$score >= best_sc - .TIE_TOL, na.rm = TRUE)
  rest     <- tuning_log[-best_idx, , drop = FALSE]
  runner_up <- if (nrow(rest) && any(!is.na(rest$score))) {
    ru <- rest[.pick_best(rest), , drop = FALSE]; rownames(ru) <- NULL; ru
  } else NULL
  if (verbose && n_tied > 1L)
    message(sprintf("  %d configurations tie on the %s score; tie broken by: %s.",
                    n_tied, tuning_metric, .TIE_RULE))
  
  if (verbose) {
    message(sprintf(paste0(
      "\n>> Best: regions_per_partition = %d, weight_exponent = %g, ",
      "N_tree = %d, final_regions = %d ",
      "(%s score = %.4f)"),
      best_cfg$regions_per_partition, best_cfg$weight_exponent,
      best_cfg$N_tree, best_cfg$final_regions, tuning_metric, best_sc))
    message("  Finalising the scored winner (no refit) ...")
  }
  
  best_key <- .tune_key(best_cfg)
  best_fit <- get(best_key, envir = fit_store, inherits = FALSE)
  rm(list = ls(fit_store), envir = fit_store)   # free losers before any rebuild
  memb_mat <- get_memb(best_cfg$regions_per_partition,
                       best_cfg$weight_exponent,
                       best_cfg$N_tree)
  best_res <- .finish_consensus(
    .ensemble_meta(memb_mat), best_fit, best_cfg$final_regions, graph = graph,
    large_n_threshold = large_n_threshold,
    sample_size = sample_size, n_cores = n_cores,
    refine_boundary = refine_boundary, min_size = min_size,
    verbose = verbose)
  # SIL_CA and EVI depend only on (ensemble, partition): a mismatch means the
  # returned regions are not the ones that were scored.
  if (!isTRUE(all.equal(c(best_res$avg_sil_width, best_res$expected_vi),
                        c(tuning_log$SIL_CA[best_idx], tuning_log$EVI[best_idx])))) {
    warning("The returned regions differ from the scored winner. ",
            "Please report this as a bug.", call. = FALSE)
  }
  
  structure(
    list(
      best_result               = best_res,
      best_config               = best_cfg,
      best_regions_per_partition = best_cfg$regions_per_partition,
      best_final_regions        = best_cfg$final_regions,
      best_weight_exponent      = best_cfg$weight_exponent,
      best_N_tree                 = best_cfg$N_tree,
      best_score                = best_sc,
      n_tied                    = as.integer(n_tied),
      runner_up                 = runner_up,
      tie_rule                  = .TIE_RULE,
      tuning_log                = tuning_log,
      candidates                = cand,
      active_axes               = active_axes,
      n_ensembles_built         = ens_cache$built,
      tuning_metric                    = tuning_metric,
      voters                    = voters,
      tuning_strategy                  = tuning_strategy,
      prior_strength            = prior_strength,
      cut_rule                  = cut_rule,
      refine_boundary           = refine_boundary,
      min_size                  = min_size,
      ensemble                  = best_res$ensemble_meta$memb_mat
    ),
    class = "pulse_regions_tuning"
  )
}


#' Canonical field order for a "pulse_regions" object.
#'
#' Single source of truth for the return contract. Every exit path in
#' get_regions() produces exactly these names in this order, so a consumer can
#' rely on the shape without branching on how the object was made.
#' @keywords internal
#' @noRd
.PULSE_REGIONS_FIELDS <- c(
  "memberships", "entropy", "hard_labels", "k",
  "config", "tuned", "score", "tuning_metric", "tuning_strategy", "tuning_log",
  "candidates", "n_ensembles_built",
  "prior_strength", "coverage", "fit", "trees", "ensemble",
  "p0", "n_eff"
)


#' Row-wise Shannon entropy, normalised by log(k).
#'
#' Rows are renormalised before the calculation, so this is defined for a
#' possibilistic matrix as well, but for such a matrix it measures only the row
#' SHAPE and discards the overall level. Entropy of a cohesion matrix is
#' therefore not the same quantity as entropy of a membership matrix; the
#' contract stores the membership version.
#'
#' @param U Numeric N x k matrix with non-negative entries.
#' @return Numeric length-N vector within 0,1. All zeros when k == 1.
#' @keywords internal
#' @noRd
.row_entropy <- function(U) {
  eps <- .Machine$double.eps
  k   <- ncol(U)
  if (k < 2L) return(rep(0, nrow(U)))
  Un <- U / pmax(rowSums(U), eps)
  -rowSums(Un * log(pmax(Un, eps))) / log(k)
}

#' Coerce a configuration to the canonical one-row data frame.
#'
#' Accepts a named list (untuned branch) or a one-row data frame
#' (tune_regions()$best_config). Guarantees the same column names, order and
#' storage modes either way, so config rows from different runs rbind cleanly.
#'
#' @keywords internal
#' @noRd
## Method settings stored next to the tuning axes in every config.
.METHOD_COLS <- c("cut_rule", "refine_boundary", "min_size")

.as_config <- function(x) {
  if (is.list(x) && !is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  stopifnot(is.data.frame(x), nrow(x) == 1L)
  miss <- setdiff(.TUNE_AXES, names(x))
  if (length(miss))
    stop("configuration is missing axis column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  for (nm in .METHOD_COLS) if (!nm %in% names(x)) x[[nm]] <- NA
  x <- x[, c(.TUNE_AXES, .METHOD_COLS), drop = FALSE]
  for (nm in setdiff(.TUNE_AXES, "weight_exponent")) x[[nm]] <- as.integer(x[[nm]])
  x$weight_exponent <- as.numeric(x$weight_exponent)
  x$cut_rule        <- as.character(x$cut_rule)
  x$refine_boundary <- as.logical(x$refine_boundary)
  x$min_size        <- as.integer(x$min_size)
  rownames(x) <- NULL
  x
}

#' Construct a "pulse_regions" object
#'
#' The single exit point for get_regions(). Every branch (tuned, untuned, exact
#' prior) routes through here, so the returned shape cannot drift between them.
#' Derived quantities (entropy, k) are computed once, here, rather
#' than independently in each branch.
#'
#' @param fit A "pulse_consensus" object from compute_consensus(), or NULL for
#'   the exact-prior degenerate case, in which case \code{memberships} and
#'   \code{hard_labels} must be supplied directly.
#' @param config One-row data frame or named list giving the four tuning axes.
#' @param trees The "pulse_trees" object used, or NULL.
#' @param prior_strength Numeric lambda_E actually applied.
#' @param tuned Logical. TRUE if a search was run.
#' @param score,tuning_metric,tuning_strategy Tuning metadata. NA in the untuned branch: the
#'   hardened silhouette is deliberately NOT substituted, because it is not
#'   comparable to a tuned score and a populated-but-incomparable slot is worse
#'   than an empty one.
#' @param tuning_log,candidates,n_ensembles_built Search provenance, or NULL/NA.
#' @param coverage Length-N logical, exact-prior case only. NULL otherwise.
#' @param memberships,hard_labels Supplied directly only when \code{fit} is
#'   NULL. Ignored otherwise.
#' @param validate Logical. Run the shape assertions. TRUE except in tight
#'   internal loops.
#' @param ensemble N x N_tree label matrix, or NULL. Stored once: if it is
#'   identical to \code{fit$ensemble_meta$memb_mat}, that copy is dropped
#'   (saveRDS() would otherwise write the matrix twice); if NULL, it is taken
#'   from there.
#'
#' @return An object of class "pulse_regions": a list with exactly the fields in
#'   \code{.PULSE_REGIONS_FIELDS}.
#' @keywords internal
#' @noRd
.new_pulse_regions <- function(fit,
                               config,
                               trees             = NULL,
                               prior_strength    = 0,
                               tuned             = FALSE,
                               score             = NA_real_,
                               tuning_metric            = NA_character_,
                               tuning_strategy          = NA_character_,
                               tuning_log        = NULL,
                               candidates        = NULL,
                               n_ensembles_built = NA_integer_,
                               coverage          = NULL,
                               memberships       = NULL,
                               hard_labels     = NULL,
                               validate          = TRUE,
                               ensemble          = NULL) {
  
  degenerate <- is.null(fit)
  
  # --- Pull the primitives ---------------------------------------------*
  if (!degenerate) {
    need <- c("memberships", "hard_labels")
    miss <- setdiff(need, names(fit))
    if (length(miss))
      stop("compute_consensus() did not return: ",
           paste(miss, collapse = ", "),
           ". The pulse_regions contract cannot be met.", call. = FALSE)
    
    U             <- fit$memberships
    hard_labels <- fit$hard_labels
    mm <- fit$ensemble_meta$memb_mat
    if (!is.null(mm)) {
      if (is.null(ensemble)) ensemble <- mm
      if (identical(ensemble, mm)) fit$ensemble_meta$memb_mat <- NULL
    }
  } else {
    if (is.null(memberships) || is.null(hard_labels))
      stop("`fit = NULL` requires `memberships` and `hard_labels`.",
           call. = FALSE)
    U            <- memberships
  }
  .p0v <- if (degenerate) NA_real_ else {
    v <- tryCatch(fit$ensemble_meta$chance_coassoc,
                  error = function(e) NULL)
    if (is.numeric(v) && length(v) > 0L) as.numeric(v) else NA_real_
  }
  # --- Derived quantities, computed exactly once ------------------------*
  k <- ncol(U)
  
  # Entropy is of the MEMBERSHIPS, always. See .row_entropy() on why the
  # cohesion version is a different quantity and is not stored here.
  entropy <- .row_entropy(U)
  
  config <- .as_config(config)
  
  # --- Shape assertions -------------------------------------------------*
  if (isTRUE(validate)) {
    N <- nrow(U)
    stopifnot(is.matrix(U), is.numeric(U), N >= 1L, k >= 1L)
    stopifnot(length(hard_labels) == N)
    stopifnot(length(entropy) == N)
    
    # Rows sum to 1 up to floating point. Uncovered vertices in the
    # exact-prior case are all-zero rows and are exempt, and so are
    # NA rows (units with pi_support == 0, flagged by the fit).
    rs <- rowSums(U)
    chk <- if (!is.null(coverage)) coverage else rep(TRUE, N)
    chk <- chk & !is.na(rs)
    if (any(abs(rs[chk] - 1) > 5e-3))
      warning("membership rows do not sum to 1 (max deviation ",
              signif(max(abs(rs[chk] - 1)), 3), ").", call. = FALSE)
    if (!is.null(coverage)) stopifnot(length(coverage) == N)
    if (!degenerate && config$final_regions != k)
      warning("config$final_regions (", config$final_regions,
              ") disagrees with ncol(memberships) (", k, ").",
              call. = FALSE)
  }
  
  p0        = mean(.p0v, na.rm = TRUE)
  n_eff <- if (is.finite(p0) && p0 > 0) 1 / p0 else NA_real_ 
  
  out <- list(
    memberships       = U,
    entropy           = entropy,
    hard_labels     = as.integer(hard_labels),
    k                 = as.integer(k),
    config            = config,
    tuned             = isTRUE(tuned),
    score             = as.numeric(score),
    tuning_metric            = as.character(tuning_metric),
    tuning_strategy          = as.character(tuning_strategy),
    tuning_log        = tuning_log,
    candidates        = candidates,
    n_ensembles_built = as.integer(n_ensembles_built),
    prior_strength    = prior_strength,
    coverage          = coverage,
    fit               = fit,
    trees             = trees,
    ensemble          = ensemble,
    p0                = p0,
    n_eff             = n_eff
  )
  
  # Enforce the canonical order, and fail loudly if a field was added here
  # without being added to the contract.
  stopifnot(setequal(names(out), .PULSE_REGIONS_FIELDS))
  structure(out[.PULSE_REGIONS_FIELDS], class = "pulse_regions")
}


#' @rdname get_regions
#' @param x An object of class \code{"pulse_regions"}.
#' @param ... Unused.
#' @return \code{print.pulse_regions} returns \code{x} invisibly.
#' @export
print.pulse_regions <- function(x, ...) {
  cat(sprintf("<pulse_regions>  %d units x %d regions | %s\n",
              nrow(x$memberships), x$k,
              if (x$tuned) sprintf("tuned (%s, %s)", x$tuning_metric, x$tuning_strategy)
              else "single configuration"))
  cfg <- x$config
  cat(sprintf("  ir = %d | gamma = %g | N_tree = %d\n",
              cfg$regions_per_partition, cfg$weight_exponent, cfg$N_tree))
  cat(sprintf("  cut_rule = %s | refine_boundary = %s | min_size = %s\n",
              cfg$cut_rule, cfg$refine_boundary, cfg$min_size))
  cat(sprintf("  mean normalised entropy = %.3f | lambda_E = %s\n",
              mean(x$entropy, na.rm = TRUE), format(x$prior_strength)))
  if (!is.null(x$coverage))
    cat(sprintf("  prior coverage: %d/%d vertices\n",
                sum(x$coverage), length(x$coverage)))
  invisible(x)
}

#' Robust scale of the edge weights, the base unit of the tree perturbation
#'
#' MAD, falling back to IQR / 1.349, then sd, then 1 when a scale is 0, so a
#' graph where most neighbour pairs are identical still gets diverse trees.
#' @keywords internal
#' @noRd
.noise_scale <- function(w, warn = TRUE) {
  s <- stats::mad(w)
  if (is.finite(s) && s > 0) return(s)
  cand <- c(`IQR / 1.349` = stats::IQR(w) / 1.349, sd = stats::sd(w), `1` = 1)
  use  <- names(cand)[which(is.finite(cand) & cand > 0)[1L]]
  if (warn)
    warning(sprintf(paste0(
      "The MAD of the edge weights is 0 (more than half are equal), which ",
      "would make every tree the same minimum spanning tree. Using %s ",
      "(= %.4g) as the perturbation unit instead."), use, cand[[use]]),
      call. = FALSE)
  cand[[use]]
}

.turnover <- function(t0_ids, t1_ids) 1 - length(intersect(t0_ids, t1_ids)) / length(t0_ids)
.calibrate_eta_empirical <- function(graph,
                                     w,
                                     s,
                                     target,
                                     n_probe = 20,
                                     tol = 0.02,
                                     max_iter = 12,
                                     warn = TRUE) {
  if (target == 0) return(0)
  
  # Cyclomatic number = edges that are not in any spanning tree.
  # 0 means the graph is a tree or forest, so it has one spanning tree.
  n_v     <- igraph::vcount(graph)
  n_comp  <- igraph::components(graph)$no
  n_tedge <- n_v - n_comp
  n_cyc   <- igraph::ecount(graph) - n_tedge
  if (n_cyc == 0L) {
    if (warn) warning(paste0(
      "The graph is a tree (no cycles), so it has exactly one spanning ",
      "tree: every sampled tree is identical and the ensemble varies ",
      "only through its cuts. target_edge_turnover is ignored (eta = 0). ",
      "Add adjacency (e.g. bridges or a wider contiguity rule) if tree ",
      "diversity is needed."), call. = FALSE)
    return(0)
  }
  # Each swapped tree edge needs a non-tree edge to replace it.
  max_turn <- min(1, n_cyc / n_tedge)
  
  t0 <- igraph::E(igraph::mst(graph, weights = w))$eid
  measure <- function(eta) {
    mean(replicate(n_probe, {
      # Not clamped at 0: the MST depends only on the weight
      # order, and the sampler does not clamp either.
      wp <- w + eta * s * stats::rnorm(length(w))
      t1 <- igraph::E(igraph::mst(graph, weights = wp))$eid
      .turnover(t0, t1)
    }))
  }
  
  lo <- 1e-4
  hi <- 1000
  best_eta <- NA_real_
  best_tm  <- NA_real_
  for (i in seq_len(max_iter)) {
    mid <- sqrt(lo * hi)
    tm  <- measure(mid)
    if (is.na(best_tm) || abs(tm - target) < abs(best_tm - target)) {
      best_eta <- mid
      best_tm  <- tm
    }
    if (abs(tm - target) < tol) return(mid)
    if (tm < target) lo <- mid else hi <- mid
  }
  
  if (warn) {
    why <- if (target > max_turn)
      sprintf(paste0("The graph has only %d edge(s) outside a spanning tree, ",
                     "so turnover cannot exceed %.3f."), n_cyc, max_turn)
    else
      "Try a lower target_edge_turnover, or raise n_probe / max_iter."
    warning(sprintf(paste0(
      "eta calibration did not reach target_edge_turnover = %.3f ",
      "(closest: %.3f at eta = %.3g, after %d steps). %s"),
      target, best_tm, best_eta, max_iter, why), call. = FALSE)
  }
  best_eta
}


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
#'   but ecologically meaningless gradients therefore count as signal here.
#'
#'   The tree noise multiplier \eqn{\eta} is calibrated to
#'   \code{target_edge_turnover} on the first replicate and reused for the
#'   others. Replicates are exchangeable permutations of the same matrix, so
#'   one calibration serves all of them; the noise scale itself is still each
#'   replicate's own weight MAD. The value is returned in \code{settings$eta}.
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
#' @param large_n_threshold,sample_size,n_cores As in \code{\link{tune_regions}}.
#' @param seed Optional integer seed; the caller's RNG state is restored.
#' @param verbose Logical.
#'
#' @return An object of class \code{"pulse_null"}: \code{null} (one row per
#'   replicate and k), \code{summary} (one row per k and criterion: observed,
#'   null mean, null 95\% bound in the "better" direction, \code{excess} =
#'   observed minus null mean with the sign set so positive is better than the
#'   null, and the empirical \code{p}), \code{config}, and \code{B}.
#'   \code{settings} holds \code{cut_rule}, \code{refine_boundary},
#'   \code{min_size}, \code{stability_k} and the reused \code{eta}.
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
      "stored on `graph`."), r), call. = FALSE)
  # cfg  <- list(
  #   regions_per_partition = as.integer(regions_per_partition),
  #   N_tree = as.integer(N_tree),
  #   weight_exponent = as.numeric(weight_exponent)
  # )
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
    # ens <- compute_ensemble(tr,
    #                         N_tree                = N_tree, 
    #                         graph                 = g_b,
    #                         regions_per_partition = regions_per_partition,
    #                         weight_exponent       = weight_exponent,
    #                         size_influence        = size_influence,
    #                         cut_rule              = cut_rule, 
    #                         min_size              = min_size,
    #                         verbose               = FALSE
    #                       )
    ens_fixed <- 
      if (rpp_tied) NULL else build(rpp_of(NA))
    stab_seed <- sample.int(.Machine$integer.max, 1L)
    out[[b]] <- do.call(rbind, lapply(final_regions, function(k) {
      ens <- if (rpp_tied) build(rpp_of(k)) else ens_fixed
      res <- tryCatch(
        compute_consensus(ens, k, graph = g_b,
                          refine_boundary = refine_boundary,
                          min_size = min_size,
                          large_n_threshold = large_n_threshold,
                          sample_size = sample_size, n_cores = n_cores,
                          report_native = FALSE, verbose = FALSE),
        error = function(e) NULL)
      cv <- if (is.null(res))
        as.data.frame(as.list(stats::setNames(rep(NA_real_, length(cvi_cols)), cvi_cols)))
      else .consensus_cvis(res, ens$memberships, env = env_b,
                           stability_B = if (k %in% stability_k) stability_B else 0L,
                           graph = g_b,
                           stability_seed = stab_seed,
                           large_n_threshold = large_n_threshold,
                           sample_size = sample_size, n_cores = n_cores)
      cbind(data.frame(rep = b, final_regions = k, failed = is.null(res)), cv[, cvi_cols, drop = FALSE])
    }))
  }
  null <- do.call(rbind, out)
  
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
    observed <- .validity_row(observed, env, stability_B = stability_B, seed = seed,
                              large_n_threshold = large_n_threshold,
                              sample_size = sample_size, n_cores = n_cores)
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
    rows[[length(rows) + 1L]] <- data.frame(
      final_regions = k, criterion = cr, observed = obs,
      null_mean = if (length(x)) mean(x) else NA_real_,
      null_95   = if (length(x)) unname(stats::quantile(x, if (sgn > 0) 0.95 else 0.05)) else NA_real_,
      excess    = if (length(x)) sgn * (obs - mean(x)) else NA_real_,
      p          = if (length(x) + n_fail > 0L && is.finite(obs))
        (1 + sum(sgn * x >= sgn * obs)) / (length(x) + n_fail + 1) else NA_real_,
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
  s <- x$summary[x$summary$criterion %in% c("CH_FIT", "STAB"), ]
  s[] <- lapply(s, function(v) if (is.numeric(v)) signif(v, 3) else v)
  print(s, row.names = FALSE)
  invisible(x)
}
# Contiguity-preserving SSD descent on a hard partition. Moves a boundary
# unit to an adjacent region if total within-region SSD drops and the old
# region stays connected and keeps at least min_size units. Same objective
# as the "ssd" cut rule.
.refine_regions <- function(z, graph, X, max_pass = 20L, min_size = 1L) {
  z   <- match(z, sort(unique(z))); k <- max(z)
  min_size <- .check_min_size(min_size, k, length(z))
  nb  <- lapply(igraph::adjacent_vertices(graph, igraph::V(graph)), as.integer)
  connected <- .region_connect_checker(nb, length(z))
  el  <- igraph::as_edgelist(graph, names = FALSE)
  n_r <- tabulate(z, k)
  S   <- rowsum(X, z, reorder = TRUE)
  moves <- 0L
  for (pass in seq_len(max_pass)) {
    moved <- 0L
    cut <- z[el[, 1]] != z[el[, 2]]
    bnd <- sort(unique(c(el[cut, 1], el[cut, 2])))
    for (i in bnd) {
      a <- z[i]; cand <- unique(z[nb[[i]]]); cand <- cand[cand != a]
      if (!length(cand) || n_r[a] <= min_size) next
      out  <- n_r[a] / (n_r[a] - 1) * sum((X[i, ] - S[a, ] / n_r[a])^2)
      gain <- vapply(cand, function(b)
        out - n_r[b] / (n_r[b] + 1) * sum((X[i, ] - S[b, ] / n_r[b])^2), numeric(1))
      if (max(gain) <= 0) next
      b    <- cand[which.max(gain)]
      same <- nb[[i]][z[nb[[i]]] == a]
      if (!connected(i, same, z)) next
      z[i] <- b
      S[a, ] <- S[a, ] - X[i, ]; S[b, ] <- S[b, ] + X[i, ]
      n_r[a] <- n_r[a] - 1L;     n_r[b] <- n_r[b] + 1L
      moved <- moved + 1L
    }
    moves <- moves + moved
    if (!moved) break
  }
  list(z = z, moves = moves, passes = pass)
}

# Returns check(i, same, z): does region z[i] stay connected when unit i
# leaves it? `same` are i's neighbours in that region. Exact.
#
# Breadth-first search inside the region, started from all of `same` at once,
# one label per start. Labels merge when their searches meet. The answer is
# TRUE as soon as one label is left, FALSE as soon as a label runs out of
# units to visit while others remain. A unit whose neighbours are linked
# around it (most moves in planar adjacency) is settled after two or three
# steps; at a true cut vertex the search stops when the smaller piece is
# exhausted. The whole region is searched only if i closes a long loop.
# `mark` is shared across calls and reset after each one, so a call costs
# what it visits, not N.
.region_connect_checker <- function(nb, n) {
  mark <- integer(n)
  function(i, same, z) {
    same <- unique(same[same != i])
    L    <- length(same)
    if (L < 2L) return(TRUE)
    a    <- z[i]
    up   <- seq_len(L)                      # union-find over labels
    find <- function(x) {
      while (any(p <- up[x] != x)) x[p] <- up[x[p]]
      x
    }
    seen <- same
    on.exit(mark[seen] <<- 0L)
    mark[same] <<- seq_len(L)
    front <- same
    repeat {
      nbr  <- nb[front]
      to   <- unlist(nbr, use.names = FALSE)
      fl   <- rep(find(mark[front]), lengths(nbr))
      keep <- to != i & z[to] == a
      to   <- to[keep]; fl <- fl[keep]
      m    <- mark[to]
      old  <- m > 0L
      tn   <- to[!old]; ln <- fl[!old]
      first <- ln[match(tn, tn)]          # label of first arrival
      pa <- c(fl[old], ln)
      pb <- c(find(m[old]), first)
      d  <- pa != pb
      if (any(d)) {
        pr <- unique(cbind(pa[d], pb[d]))
        for (r in seq_len(nrow(pr))) {
          x <- find(pr[r, 1L]); y <- find(pr[r, 2L])
          if (x != y) up[max(x, y)] <- min(x, y)
        }
      }
      roots <- unique(find(seq_len(L)))
      if (length(roots) == 1L) return(TRUE)
      keep  <- !duplicated(tn)
      front <- tn[keep]
      if (!length(front)) return(FALSE)
      mark[front] <<- first[keep]
      seen <- c(seen, front)
      if (!all(roots %in% find(mark[front]))) return(FALSE)
    }
  }
}

#' Internal: parent-edge positions of tree j
#'
#' Column j of \code{trees$edge_ids}, in the order \code{igraph::mst()}
#' returned the edges. Stops with a re-run hint on objects saved before
#' the edge-id format.
#' @keywords internal
#' @noRd
.tree_edge_ids <- function(trees, j) {
  ids <- trees$edge_ids
  if (is.null(ids))
    stop(paste0("This `pulse_trees` object uses the old format (`$trees`, ",
                "no `$edge_ids`). Re-run sample_spanning_trees()."),
         call. = FALSE)
  if (nrow(ids) != trees$n_vertices - 1L)
    stop("Spanning tree does not span the parent vertex set", call. = FALSE)
  ids[, j]
}

#' Internal: component labels of the forest given by an edge-endpoint matrix
#'
#' Labels follow igraph::components(): numbered by first vertex index, so
#' they equal the labels of the same forest built any other way.
#' @keywords internal
#' @noRd
.labels_from_edges <- function(ep, n_v) {
  as.integer(igraph::components(igraph::make_graph(
    as.vector(t(ep)), n = n_v, directed = FALSE))$membership)
}

