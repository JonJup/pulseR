#' Borda rank-aggregation across sub-metric columns
#'
#' For each column of \code{df}, ranks rows (1 = best per the direction); NAs
#' receive the worst rank within that column; fully-NA columns are dropped;
#' remaining ranks are summed row-wise. The optimal row MINIMISES the returned
#' Borda score.
#'
#' @param df Data frame of sub-metric values (one column per voter).
#' @param directions Character vector of length \code{ncol(df)}; each element
#'   \code{"min"} or \code{"max"}.
#' @return Numeric vector of length \code{nrow(df)}; all-NA if every voter is NA.
#'
#' @details Uses \code{ties.method = "average"} so genuine ties do not falsely
#'   separate candidates. Weighted Borda is intentionally not used: each voter
#'   is a different *aspect* of validity (compactness/separation, fuzziness-
#'   aware silhouette, stability), and over-weighting one would re-introduce the
#'   bias direction we are trying to balance.
#'
#' @keywords internal
#' @noRd
.borda_score <- function(df, directions) {
        stopifnot(ncol(df) == length(directions))
        if (!is.null(names(directions))) {
                stopifnot(all(names(directions) %in% names(df)))
                df <- df[, names(directions), drop = FALSE]
        }
        n <- nrow(df)
        
        ranks_list <- Map(function(values, dir) {
                if (all(is.na(values))) return(NULL)
                r <- switch(dir,
                            min = rank(values,  na.last = "keep", ties.method = "average"),
                            max = rank(-values, na.last = "keep", ties.method = "average"),
                            stop("direction must be 'min' or 'max'", call. = FALSE))
                na_mask <- is.na(r)
                if (any(na_mask)) {
                        # NAs share the average of the unused worst ranks, so every
                        # column contributes the same total regardless of coverage.
                        m <- sum(!na_mask)
                        r[na_mask] <- mean(seq.int(m + 1L, n))
                }
                r
        }, df, directions)
        
        ranks_list <- Filter(Negate(is.null), ranks_list)
        if (!length(ranks_list)) return(rep(NA_real_, n))
        rowMeans(do.call(cbind, ranks_list))
}



#' Chance co-association per ensemble column
#'
#' For each column, the probability that two units drawn uniformly at random
#' share an intermediate region. Averaged over columns this is the floor
#' against which co-association to a medoid should be read; its reciprocal is
#' the effective number of intermediate regions (inverse Simpson), which is
#' the same n_eff diagnostic the size-biased cut fix was tuned against.
#'
#' Returned per column so that an ensemble sliced to n_st' < n_st can take
#' the mean over its own columns without an O(N M) rescan.
#'
#' @param memb_mat Integer matrix N x M of intermediate region labels.
#' @return Numeric vector of length M, each element in (0, 1].
#' @keywords internal
#' @noRd
.chance_coassoc_by_col <- function(memb_mat) {
        N <- nrow(memb_mat)
        vapply(seq_len(ncol(memb_mat)), function(t) {
                n <- tabulate(memb_mat[, t])
                sum((n / N)^2)
        }, numeric(1))
}


.check_contiguity <- function(graph, lab) {
        n  <- igraph::vcount(graph)
        stopifnot(length(lab) == n)
        ep <- igraph::ends(graph, igraph::E(graph), names = FALSE)
        
        same <- !is.na(lab[ep[, 1]]) & !is.na(lab[ep[, 2]]) & lab[ep[, 1]] == lab[ep[, 2]]
        g2   <- igraph::make_graph(as.vector(t(ep[same, , drop = FALSE])),
                                   n = n, directed = FALSE)
        comp <- igraph::components(g2)$membership
        
        keep  <- !is.na(lab)
        frag  <- tapply(comp[keep], lab[keep], function(x) length(unique(x)))
        list(contiguous = all(frag == 1L),
             fragments  = frag[frag > 1L],
             n_orphans  = sum(keep) - sum(tapply(comp[keep], lab[keep],
                                                 function(x) max(tabulate(match(x, unique(x)))))))
}

#' Calinski-Harabasz index computed from a distance matrix
#'
#' Ratio of between-group to within-group dispersion, corrected for the number
#' of groups. Needs only a distance matrix and a membership matrix; no
#' coordinates and no centroids are required, so the same estimator works on
#' consensus Hamming distances (level 1) and Jensen-Shannon distances (level 3).
#'
#' @param D Distance matrix or \code{dist} object over N objects.
#' @param U N x k membership matrix.
#' @param fuzziness Fuzzifier m. Ignored when \code{crisp = TRUE}.
#' @param square Logical. If TRUE, \code{D} is squared before use. Set this so
#'   that the quantity entering the formula is the SQUARED distance in the
#'   sense of an ordinary geometric layout. See details in the notes below.
#' @param crisp Logical. Use one-hot argmax weights instead of u^m.
#' @return Numeric scalar; higher is better. NA_real_ if undefined.
#' @keywords internal
#' @noRd
.ch_distance <- function(D, U, fuzziness = 2, square = T, crisp = FALSE){
        
        if(inherits(D, "dist")) D <- as.matrix(D)
        N <- nrow(U) 
        k <- ncol(U)
        if (N != nrow(D)) stop("Dimensions of D and N do not agree", call. = FALSE)
        if (k < 2 || N < k) return(NA_real_)
        D2 <- if(square) D*D else D
        if (crisp){
                Wt <- matrix(0, N, k)
                Wt[cbind(seq_len(N), max.col(U, ties.method = "first"))] <- 1
        } else {
                Wt <- U ^ fuzziness
        }
        
        # Group size-like notion. Excatlcy the group size for crisp. 
        mass <- colSums(Wt)
        if (any(!is.finite(mass)) || any(mass<=0)) return(NA_real_)
        
        # Within group distance. Wt provides fuzzy notion of within-groupies. 
        # The 0.5 corrects for counting every relationship twice (once from each side). 
        within_k <- 0.5 * colSums(Wt * (D2 %*% Wt)) / mass
        W <- sum(within_k)
        
        # Total spread 
        Tot <- 0.5 * sum(D2) / N
        # Between groups 
        s   <- rowSums(D2)
        h   <- rowSums(Wt)
        M   <- sum(mass)
        B   <- sum(h * s) / N - W - M * Tot / N
        if (!is.finite(B) || !is.finite(W) || B <= 0 || W <= 0) return(NA_real_)
        (B/(k - 1))/(W/(N-k))
}


#' Compute Ensemble Region Memberships via SKATER Edge Removal
#'
#' For each spanning tree, removes \code{intermediate_regions - 1} edges (sampled
#' without replacement, with probability proportional to edge weight) to produce
#' \code{intermediate_regions} connected components, then records each node's region
#' label. Returns an integer matrix (one column per tree).
#'
#' @param trees A \code{"spanning_trees"} object.
#' @param n_st Integer scalar. Number of spanning trees to use
#'   (must be \code{<= trees$n_trees}).
#' @param graph The parent \code{igraph} object the trees were sampled from.
#'   Needed only to re-resolve tree edges on older igraph versions; always
#'   supplied for safety.
#' @param intermediate_regions Integer scalar \eqn{\ge 2}. Number of regions per tree;
#'   must not exceed the number of vertices.
#' @param prior_strength Non-negative finite numeric scalar \eqn{\lambda_E}.
#'   Halving rate for within-prior edge weights. \eqn{\lambda_E = 0} (default)
#'   disables prior weighting and reproduces the unconstrained algorithm.
#'   \eqn{\lambda_E = 1} halves within-prior weights, \eqn{\lambda_E = 2}
#'   quarters them, and so on. Has no effect if \code{trees} was constructed
#'   without a \code{prior_typology}. \code{Inf} is not accepted here (use the
#'   \code{prior_strength = Inf} short-circuit of \code{\link{get_regions}} /
#'   \code{\link{tune_regions}} instead).
#' @param verbose Logical. Default \code{TRUE}.
#' @param partitions_per_tree Integer scalar \eqn{\ge 1}. Number of independent
#'   cut-sets applied to each spanning tree. Each cut-set is one regionalisation
#'   and contributes one ensemble column, so the returned matrix has
#'   \code{n_st * partitions_per_tree} columns. Default 1 (one partition per
#'   tree, the original behavior). Partitions from the same tree share its
#'   topology and are more correlated than partitions from distinct trees.
#' @param weight_exponent Non-negative numeric scalar \eqn{\gamma}. The cut
#'   probability of an edge is proportional to its (prior-adjusted) weight raised
#'   to \eqn{\gamma}. \eqn{\gamma = 0} cuts edges uniformly at random,
#'   \eqn{\gamma = 1} (default) cuts with probability proportional to weight, and
#'   larger \eqn{\gamma} concentrates cuts on the heaviest edges (approaching
#'   deterministic heaviest-edge removal as \eqn{\gamma \to \infty}).
#'
#' @return An object of class \code{"ensemble_memberships"} — a list with:
#'   \describe{
#'     \item{\code{memberships}}{Integer matrix (\code{N x n_st}). Each column
#'       is a region assignment from one spanning tree.}
#'     \item{\code{n_st}}{Ensemble size used.}
#'     \item{\code{intermediate_regions}}{Number of regions per tree.}
#'     \item{\code{n_vertices}}{Number of vertices.}
#'     \item{\code{prior_strength}}{The prior strength used.}
#'     \item{\code{chance_coassoc}}{Mean chance co-association \eqn{\bar p_0}; the reference scale for typicality.}
#'   }
#'@details When both a prior and \eqn{\gamma \ne 1} are in use, weights are first
#'   exponentiated and the weighted by the prior. Thus the halvings are relative
#'   to the exponentiated weights. 
#' @importFrom igraph E ends make_graph components
#' @export
compute_ensemble_memberships <- function(trees, 
                                         n_st, 
                                         graph, 
                                         intermediate_regions, 
                                         prior_strength = 0, 
                                         partitions_per_tree = 1,
                                         weight_exponent = 1, 
                                         verbose = TRUE,
                                         size_influence = 2) {
        
        stopifnot(inherits(trees, "spanning_trees"))
        stopifnot(length(n_st) == 1L, n_st >= 2, n_st <= trees$n_trees)
        stopifnot(length(intermediate_regions) == 1L, intermediate_regions >= 2, intermediate_regions <= trees$n_vertices)
        stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L, prior_strength >= 0)
        stopifnot(length(partitions_per_tree) == 1L, is.finite(partitions_per_tree), partitions_per_tree >= 1)
        stopifnot(is.numeric(weight_exponent), length(weight_exponent) == 1L, is.finite(weight_exponent), weight_exponent >= 0)
        if (is.infinite(prior_strength))  stop("`prior_strength = Inf` is not supported here; use the exact-prior ", "short-circuit in get_regions()/tune_regions().", call. = FALSE)
        partitions_per_tree <- as.integer(partitions_per_tree)
        
        if (verbose) {
                msg <- sprintf(
                        "Computing ensemble memberships (n_st = %d, intermediate_regions = %d",
                        n_st, intermediate_regions
                )
                if (partitions_per_tree > 1L)
                        msg <- paste0(msg, sprintf(", partitions_per_tree = %d", partitions_per_tree))
                if (weight_exponent != 1)
                        msg <- paste0(msg, sprintf(", weight_exponent = %g", weight_exponent))
                if (prior_strength > 0 && !is.null(trees$within_prior_edge))
                        msg <- paste0(msg, sprintf(", lambda_E = %g", prior_strength))
                message(msg, ") ...")
        }
        
        t0 <- proc.time()[["elapsed"]]
        memb <- .compute_ensemble_memberships_internal(
                graph = graph, 
                trees = trees, 
                n_st = n_st,
                intermediate_regions = intermediate_regions,
                prior_strength       = prior_strength,
                partitions_per_tree  = partitions_per_tree,
                weight_exponent      = weight_exponent,
                size_influence       = size_influence
        )
        if (verbose) message(sprintf("  Done in %s", .format_duration(proc.time()[["elapsed"]] - t0)))
        
        structure(
                list(
                        memberships         = memb,
                        n_st                = as.integer(n_st),
                        intermediate_regions = as.integer(intermediate_regions),
                        partitions_per_tree = partitions_per_tree,
                        n_members           = ncol(memb),
                        weight_exponent     = weight_exponent,
                        n_vertices          = trees$n_vertices,
                        prior_strength      = prior_strength,
                        chance_coassoc      = .chance_coassoc_by_col(memb)
                ),
                class = "ensemble_memberships"
        )
}


#' Compute fuzzy-aware Cluster Validity Indices on a final fuzzy partition
#'
#' Returns a 1-row data frame of values, NOT a single scalar. The caller derives
#' a \code{score} column according to the chosen \code{metric} (see
#' \code{.derive_score_column}).
#'
#' @param res results of cluster_consensus()
#' @param env Environmental data passed through to \code{.compute_fuzzy_cvis()}.
#'  The same data used to weight the edges in \code{\link{add_edge_weight}}.
#' @param env_id Name of any ID columns in `env`, which should not be used in 
#'   the computation of cluster validity metrics. If `NULL` (default) no columns
#'   are removed from `env`.
#' @param d_to_medoids N x k matrix of distances from each object to each
#'   medoid (from \code{cluster_consensus}).
#' @param medoid_dist k x k matrix of inter-medoid distances.
#' @param memb_mat Integer matrix N x n_trees. Required if \code{SIL_F} or
#'   \code{STAB} are to be computed; may be \code{NULL} otherwise (those entries
#'   become \code{NA}).
#' @param fuzziness Fuzzifier \emph{m}; used in the XB* numerator weighting.
#' @param silhouette_hard Numeric scalar; hardened average silhouette width from
#'   CLARA/PAM (passed through as \code{SILH_HARD}).
#' @param stability_B Integer; bootstrap replicates for STAB. \code{0} skips it.
#' @param stability_subsample_frac Fraction of rows per stability replicate.
#' @param stability_seed Optional seed.
#' @param large_n_threshold CLARA/PAM cutoff for stability re-fits.
#' @param clara_samples,clara_sampsize CLARA controls for stability re-fits.
#' @param verbose Logical.
#' @return One-row \code{data.frame} with columns PC, PE, MPC, PEN, XB_star,
#'   SIL_F, STAB, SILH_HARD.
#' @keywords internal
#' @noRd
.compute_fuzzy_cvis <- function(res,
                                # env,
                                # env_id = NULL,
                                memb_mat = NULL,
                                
                                stability_B              = 0L,
                                stability_subsample_frac = 0.8,
                                stability_seed           = NULL,
                                silf_alpha               = 2,
                                large_n_threshold        = 20000L,
                                clara_samples            = 20L,
                                clara_sampsize           = NULL,
                                verbose                  = FALSE) {
        
        
        if (is.null(memb_mat)) stop("No membership matrix provided", call. = FALSE)
        fuzziness = res$fuzziness
        # U = Membership matrix 
        U <- res$memberships
        D <- parallelDist::parallelDist(memb_mat, method = "hamming")
        SIL_F <- .fuzzy_silhouette(D, U, alpha = silf_alpha)
        # Hamming already plays the role of the squared distance, so square = FALSE.
        CH_D <- .ch_distance(D, U, fuzziness = fuzziness, square = FALSE)
        # --- STAB (bootstrap ARI) ---- *
        STAB <- NA_real_
        if (stability_B > 0) {
                hard_clusters <- res$hard_clusters
                
                STAB <- tryCatch(
                        .stability_fuzzy(memb_mat = memb_mat,
                                         hard_clusters = hard_clusters, 
                                         k = res$k, 
                                         fuzziness = res$fuzziness, 
                                         B = stability_B, 
                                         subsample_frac = stability_subsample_frac, 
                                         large_n_threshold = large_n_threshold,
                                         clara_samples = clara_samples, 
                                         clara_sampsize = clara_sampsize,
                                         seed = stability_seed, 
                                         verbose  = verbose),
                        error = function(e) {
                                if (verbose) message("    [STAB] failed: ", conditionMessage(e))
                                NA_real_
                        }
                )
        }
        
        
        data.frame(SIL_F = SIL_F, 
                   CH_D = CH_D,
                   STAB = STAB)


    
}





#' Typicality of Each Unit Relative to Each Region
#'
#' Reports, for every unit and every region, the fraction of ensemble
#' partitions in which that unit shared an intermediate region with that
#' region's medoid.
#'
#' @param x A fitted object: a \code{"fuzzy_clusters"} object from
#'   \code{\link{cluster_consensus}}, or a \code{"pulse_regions"} object from
#'   \code{\link{get_regions}}.
#' @param ... Passed to methods.
#' @return Numeric matrix, N units by k regions.
#' @export
compute_typicality <- function(x, ...) UseMethod("compute_typicality")


#' @rdname compute_typicality
#'
#' @details
#' Raw typicality \eqn{s_{ij} = 1 - d_{ij}} is the proportion of ensemble
#' partitions in which unit \eqn{i} and the medoid of region \eqn{j} fall in
#' the same intermediate region. Its resolution is \eqn{1 / n\_members}, and
#' its scale depends on \code{intermediate_regions}: each partition applies
#' \code{intermediate_regions - 1} cuts, and every cut is a further
#' opportunity to separate the pair, so raw values shrink as that setting
#' rises. Raw values are therefore comparable within a fit but not across
#' fits at different resolutions.
#'
#' \code{scale = "rescaled"} re-expresses a fit at a common resolution via
#' \eqn{s^{(R_{ref}-1)/(R-1)}}. This is exact if cuts separate a given pair
#' independently, which is an approximation that has not been validated on
#' real catchment graphs; treat it as a sensitivity tool. With \code{R_ref}
#' left at \code{NULL} it is the identity and returns raw values.
#'
#' \code{scale = "chance"} returns \eqn{1 - \log s / \log \bar p_0}, where
#' \eqn{\bar p_0} is the chance co-association of two randomly chosen units.
#' It is 0 at chance, 1 for a unit never separated from the medoid, and
#' negative for units systematically cut away from the region they were
#' assigned to.
#'
#' To reduce the matrix to one value per unit for its own assigned region:
#' \preformatted{
#'   tp <- compute_typicality(fit)
#'   tp[cbind(seq_len(nrow(tp)), fit$hard_clusters)]
#' }
#'
#' @param scale One of \code{"raw"} (default), \code{"rescaled"},
#'   \code{"chance"}.
#' @param R_ref Reference \code{intermediate_regions} for
#'   \code{scale = "rescaled"}. \code{NULL} (default) uses the fit's own
#'   value, i.e. no rescaling.
#' @export
compute_typicality.fuzzy_clusters <- function(x,
                                              scale = c("raw", "rescaled", "chance"),
                                              R_ref = NULL,
                                              ...) {
        scale <- match.arg(scale)
        
        meta <- x$ensemble_meta
        M    <- x$n_members
        if (is.null(meta) || is.null(M)) {
                stop("This `fuzzy_clusters` object carries no `ensemble_meta`, so the ",
                     "resolution its distances were measured at is unknown. Refit with ",
                     "the current version of `cluster_consensus()`.", call. = FALSE)
        }
        
        s <- 1 - x$d_to_medoids
        
        if (scale == "raw") return(s)
        
        if (scale == "rescaled") {
                R <- meta$intermediate_regions
                if (is.null(R) || is.na(R) || R < 2L) {
                        stop("`intermediate_regions` is unavailable on this object; ",
                             "rescaling is undefined.", call. = FALSE)
                }
                if (is.null(R_ref)) R_ref <- R
                stopifnot(is.numeric(R_ref), length(R_ref) == 1L,
                          is.finite(R_ref), R_ref >= 2)
                return(s ^ ((R_ref - 1) / (R - 1)))
        }
        
        # scale == "chance".
        # Half-count smoothing: s is a count out of M, so a unit never co-assigned
        # with a medoid has s = 0 and log(0) = -Inf. Nudging by half a count keeps
        # the logarithm finite while leaving all other values essentially unmoved.
        p0 <- meta$chance_coassoc
        s  <- (s * M + 0.5) / (M + 1)
        1 - log(s) / log(p0)
}

#' @rdname compute_typicality
#' @export
compute_typicality.pulse_regions <- function(x, ...) {
        if (is.null(x$fit)) {
                stop("This `pulse_regions` object has no `$fit`; it was produced by an ",
                     "older version of `get_regions()`.", call. = FALSE)
        }
        compute_typicality(x$fit, ...)
}


# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 4: Fuzzy consensus clustering
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Fuzzy Consensus Clustering from Ensemble Memberships
#'
#' Clusters nodes into \code{k} groups based on an ensemble membership matrix,
#' using CLARA (for large datasets) or PAM (for small ones), then converts
#' Hamming distances to medoids into fuzzy membership values.
#'
#' @param ensemble An \code{"ensemble_memberships"} object or a raw integer
#'   matrix (N x n_trees).
#' @param k Integer scalar with \eqn{2 \le k < N}. Number of final clusters.
#' @param graph The parent \code{igraph} object the trees were sampled from.  
#' @param fuzziness Numeric \eqn{> 1}. Fuzzifier \emph{m}. Default 2.
#' @param crisp Logical. If \code{TRUE}, returns one-hot (hard) memberships.
#' @param large_n_threshold Integer. Above this N, CLARA is used instead of PAM.
#'   Default 20000.
#' @param clara_samples Integer. Number of CLARA samples. Default 50.
#' @param clara_sampsize Integer or \code{NULL}. CLARA sample size; default
#'   \code{NULL} selects \code{min(N, 200 + 2k)}.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return An object of class \code{"fuzzy_clusters"} — a list with:
#'   \describe{
#'     \item{\code{memberships}}{Numeric matrix (N x k).}
#'     \item{\code{hard_clusters}}{Integer vector of hard cluster assignments.}
#'     \item{\code{medoid_idx}}{Integer vector of medoid row indices.}
#'     \item{\code{avg_sil_width}}{Hardened average silhouette width from
#'       CLARA/PAM. NOTE: computed on the argmax labels, i.e. a crisp index;
#'       do NOT use it as the primary k-selection criterion.}
#'     \item{\code{k}}{Number of clusters.}
#'     \item{\code{fuzziness}}{The fuzziness parameter used.}
#'     \item{\code{method}}{Either \code{"clara"} or \code{"pam"}.}
#'     \item{\code{cluster_result}}{The raw \code{clara} or \code{pam} object.}
#'     \item{\code{d_to_medoids}}{N x k matrix of Hamming distances from each
#'       node to each medoid. Reused by fuzzy CVI computations.}
#'     \item{\code{medoid_dist}}{k x k matrix of inter-medoid Hamming
#'       distances. Reused by XB*.}
#'   }
#'
#' @importFrom cluster clara pam
#' @importFrom parallelDist parallelDist
#' @export
cluster_consensus <- function(ensemble, k, graph, fuzziness = 2, crisp = FALSE,
                              large_n_threshold = 20000, clara_samples = 50,
                              clara_sampsize = NULL, verbose = TRUE) {
        
        # --- Coerce input to a membership matrix ---- *
        if (inherits(ensemble, "ensemble_memberships")) {
                # memb_mat <- ensemble$memberships
                ens_meta <- list(
                        memb_mat = ensemble$memberships,
                        intermediate_regions = ensemble$intermediate_regions,
                        n.rst                = ensemble$n.rst,
                        partitions_per_tree  = ensemble$partitions_per_tree,
                        weight_exponent      = ensemble$weight_exponent,
                        prior_strength       = ensemble$prior_strength,
                        chance_coassoc       = ensemble$chance_coassoc
                )
        } else if (is.matrix(ensemble)) {
                # Bare-matrix path: no provenance available, so p0 is computed here
                # rather than left NA. Everything else is unknowable
                # memb_mat <- ensemble
                ens_meta <- list(
                        memb_mat = ensemble$memberships,
                        intermediate_regions = NA_integer_,
                        n.rst                = NA_integer_,
                        partitions_per_tree  = NA_integer_,
                        weight_exponent      = NA_real_,
                        prior_strength       = NA_real_,
                        chance_coassoc       = .chance_coassoc_by_col(ensemble)
                )
        } else {
                stop("`ensemble` must be an 'ensemble_memberships' object or an integer matrix.", call. = FALSE)
        }
        
        # --- Coerce input to a membership matrix ---- *
        n <- nrow(ens_meta$memb_mat)
        stopifnot(length(k) == 1L, k >= 2, k < n, fuzziness > 1)
        
        use_clara <- n > large_n_threshold
        
        if (verbose) message(sprintf(
                "Fuzzy consensus clustering (k = %d, N = %d, method = %s) ...",
                k, n, if (use_clara) "CLARA" else "PAM"
        ))
        
        t0 <- proc.time()[["elapsed"]]
        
        # --- Medoid clustering w .clara_hamming ---- *
        if (use_clara) {
                # CLARA operates directly on the data matrix (samples + PAM on each sample),
                # avoiding a full N x N dissimilarity matrix.
                sampsize <- if (is.null(clara_sampsize)) min(n, 200 + 2 * k) else clara_sampsize
                cl_res <- .clara_hamming(ens_meta$memb_mat,
                                         k = k,
                                         samples  = clara_samples,
                                         sampsize = sampsize,
                                         verbose  = verbose)
                medoid_idx  <- cl_res$i.med
                clustering  <- cl_res$clustering
                sil_width   <- cl_res$silinfo$avg.width
                n_eff = "clara"
                method_used <- "clara"
        } else {
                # PAM on the full Hamming dissimilarity matrix.
                coassoc <- parallelDist::parallelDist(ens_meta$memb_mat, method = "hamming")
                n_eff <- 1/(1-mean(coassoc,na.rm =T))
                pam_res <- cluster::pam(coassoc, k = k)
                medoid_idx  <- pam_res$id.med
                clustering  <- pam_res$clustering
                sil_width   <- pam_res$silinfo$avg.width
                method_used <- "pam"
                cl_res      <- pam_res
        }
        
        # --- Distances to medoids (always computed) ---- *---- *---- *---- *---- *---- *
        # Needed both for fuzzy-membership derivation and for downstream CVIs (XB*,
        # SIL_F). Computing it unconditionally keeps the crisp = TRUE path and CVI
        # consumers well-defined.
        dists_to_med <- .hamming_to_refs(ens_meta$memb_mat, medoid_idx)
        G <- .region_step_distance(graph, lab = clustering, k = k)

        # --- k x k inter-medoid distance matrix ---- *
        if (length(medoid_idx) >= 2L) {
                medoid_dist <- as.matrix(
                        parallelDist::parallelDist(
                            ens_meta$memb_mat[medoid_idx, , drop = FALSE],
                                method = "hamming"
                        )
                )
        } else {
                medoid_dist <- matrix(0, 1L, 1L)
        }
        
        # --- Build memberships ---- *---- *---- *---- *---- *---- *---- *---- *---- *---- *-
        if (crisp) {
                # One-hot assignment, vectorised (avoid an N-length R loop).
                memberships <- matrix(0, nrow = n, ncol = k)
                memberships[cbind(seq_len(n), clustering)] <- 1
        } else {
                contiguity_bool <- FALSE
                h_var = 20
                #print("Checking contiguity")
                while (!contiguity_bool){
                        #print(paste("h = ", h_var))
                        memberships <- .fuzzy_memberships(
                                d_to_medoids = dists_to_med,
                                G = G,
                                fuzziness = fuzziness,
                                h = h_var
                        )
                        memberships <- round(memberships,3)
                        clustering <- apply(memberships,1,which.max)
                        conti <- .check_contiguity(graph, clustering)
                        contiguity_bool <- conti$contiguous
                        if (!contiguity_bool & h_var != 1){
                                h_var <- h_var-1
                        } else if (h_var == 1){
                                contiguity_bool <- TRUE
                        }
                }
        }
        if (verbose) message(sprintf("  Done in %s",
                                     .format_duration(proc.time()[["elapsed"]] - t0)))

        structure(
                list(
                        memberships    = memberships,
                        hard_clusters  = as.integer(clustering),
                        medoid_idx     = medoid_idx,
                        avg_sil_width  = sil_width,
                        k              = k,
                        fuzziness      = fuzziness,
                        method         = method_used,
                        cluster_result = cl_res,
                        d_to_medoids   = dists_to_med,
                        medoid_dist    = medoid_dist,
                        attenuation    = h_var,
                        n_eff          = n_eff, 
                        n_members      = ncol(ens_meta$memb_mat),
                        ensemble_meta  = ens_meta
                ),
                class = "fuzzy_clusters"
        )
}

#' Construct an "ensemble_memberships" object
#'
#' Used both by compute_ensemble_memberships() and by the tuner when it slices
#' a cached ensemble down to a smaller n_st, so that a sliced ensemble is
#' stamped with its own column count and chance co-association rather than
#' those of the larger matrix it came from.
#'
#' @keywords internal
#' @noRd
.new_ensemble_memberships <- function(memb, n_st, intermediate_regions,
                                      partitions_per_tree, weight_exponent,
                                      prior_strength, n_vertices,
                                      p0_by_col = NULL) {
        if (is.null(p0_by_col)) p0_by_col <- .chance_coassoc_by_col(memb)
        structure(
                list(
                        memberships           = memb,
                        n_st                  = as.integer(n_st),
                        intermediate_regions  = as.integer(intermediate_regions),
                        partitions_per_tree   = as.integer(partitions_per_tree),
                        n_members             = ncol(memb),
                        weight_exponent       = as.numeric(weight_exponent),
                        n_vertices            = n_vertices,
                        prior_strength        = prior_strength,
                        chance_coassoc        = mean(p0_by_col),
                        chance_coassoc_by_col = p0_by_col
                ),
                class = "ensemble_memberships"
        )
}




# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 5a: Tuning orchestrator (five-axis)
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾





# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Step 5b: Convenience wrapper
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

#' Delineate Fuzzy Contiguity-Constrained Regions using SKATER (one-call wrapper)
#'
#' All-in-one convenience wrapper that samples spanning trees, optionally
#' chooses the ensemble size by stability analysis, optionally tunes
#' \code{intermediate_regions} / \code{final_regions}, and returns the final fuzzy
#' regionalisation. See \code{\link{tune_regions}} for the full menu of
#' \code{tuning_metric} options. Default is \code{"borda"}.
#'
#' @inheritParams sample_spanning_trees
#' @param intermediate_regions,final_regions,n_st Integer or integer vectors of
#'   candidate values. If a vector of length > 1 is supplied for
#'   \code{intermediate_regions} or \code{final_regions}, tuning is triggered; if
#'   \code{n_st} has length > 1, the ensemble size is chosen by stability.
#' @param fuzziness Numeric \eqn{> 1}. Default 2.
#' @param trees A \code{"spanning_trees"} object from
#'   \code{\link{sample_spanning_trees}}.
#' @param crisp Logical. Produce hard memberships? Default \code{FALSE}.
#' @param tuning_strategy Character. One of \code{"grid"} (default),
#'   \code{"sequential"}, \code{"iterative"}.
#' @param tuning_metric Character. See \code{\link{tune_regions}}. Default
#'   \code{"borda"}.
#' @param prior_typology Optional. Vertex attribute name (single character)
#'   or vector of length \code{vcount(graph)} giving a prior ecoregion label
#'   at each catchment. \code{NA} = no prior at that vertex. The user is
#'   expected to attach this to the polygons / graph themselves before
#'   calling. Default \code{NULL}.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda_E}. The
#'   parameter has a direct interpretation as the number of halvings applied
#'   to within-region edge weights: \eqn{\lambda_E = 0} (default) disables
#'   prior weighting, \eqn{\lambda_E = 1} halves the weight,
#'   \eqn{\lambda_E = 2} quarters it, and so on. \code{Inf} returns the prior
#'   typology exactly (requires \code{prior_typology}).
#' @param partitions_per_tree Integer scalar \eqn{\ge 1}. Number of independent
#'   cut-sets applied to each spanning tree. Each cut-set is one regionalisation
#'   and contributes one ensemble column, so the returned matrix has
#'   \code{n_st * partitions_per_tree} columns. Default 1 (one partition per
#'   tree, the original behaviour). Partitions from the same tree share its
#'   topology and are more correlated than partitions from distinct trees.
#' @param weight_exponent Non-negative numeric scalar \eqn{\gamma}. The cut
#'   probability of an edge is proportional to its (prior-adjusted) weight raised
#'   to \eqn{\gamma}. \eqn{\gamma = 0} cuts edges uniformly at random,
#'   \eqn{\gamma = 1} (default) cuts with probability proportional to weight, and
#'   larger \eqn{\gamma} concentrates cuts on the heaviest edges (approaching
#'   deterministic heaviest-edge removal as \eqn{\gamma \to \infty}).
#' @param max_iter Integer. Maximum iterations for the \code{"iterative"}
#'   tuning strategy. Default 10.
#' @param stability_B Integer. Only used if tuning is necesary. Passed to 
#'   \code{\link{tune_regions}}. Bootstrap replicates for the stability
#'   sub-metric. Set to 0 to skip. Default 25.
#' @param verbose Logical. Default \code{TRUE}.
#' @param seed Optional integer seed; the caller's RNG state is restored on
#'   exit.
#'
#' @return A plain list with:
#'   \describe{
#'     \item{\code{memberships}}{Final N x k membership matrix.}
#'     \item{\code{hard_clusters}}{Integer vector of hard cluster labels.}
#'     \item{\code{best_intermediate_regions}, \code{best_final_regions},
#'       \code{best_n_st}}{The selected hyper-parameters.}
#'     \item{\code{best_score}}{For a tuned run, the chosen \code{tuning_metric}
#'       score; for a single-configuration run, the hardened average silhouette
#'       width (these are NOT directly comparable).}
#'     \item{\code{tuning_log}}{The tuning data frame, or \code{NULL} if no
#'       tuning was performed.}
#'     \item{\code{n_st_stability}}{The stability data frame, or \code{NULL}.}
#'     \item{\code{trees}}{The \code{"spanning_trees"} object used.}
#'     \item{\code{cluster_result}}{The raw CLARA/PAM object.}
#'     \item{\code{prior_strength}}{The prior strength used.}
#'   }
#'   When \code{prior_strength = Inf}, the exact-prior list (see
#'   \code{.exact_prior_result}) is returned instead.
#'
#' @seealso \code{\link{sample_spanning_trees}},
#'   \code{\link{compute_ensemble_memberships}},
#'   \code{\link{cluster_consensus}}, \code{\link{tune_regions}}.
#'
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- sample_gnp(200, 0.05)
#' E(g)$weight <- runif(ecount(g))
#' res <- get_regions(g, intermediate_regions = c(10, 20), final_regions = c(3, 4),
#'                          n_st = 30, seed = 1)
#' table(res$hard_clusters)
#' }
#'
#' @importFrom igraph vcount
#' @export
get_regions <- function(graph,
                        intermediate_regions = 50,
                        final_regions,
                        n_st,
                        fuzziness  = 2,
                        trees = NULL,
                        crisp      = FALSE,
                        tuning_strategy = c("grid", "sequential", "iterative"),
                        tuning_metric   = c(
                                "borda"
                                # "fuzzy_calinski_harabasz",
                                # "Fukuyama_Sugeno",
                                # "dunn5",
                                # "fuzzy_hyper_volume",
                                # "stability"
                        ),
                        stability_B      = 25L,
                        prior_typology = NULL,
                        prior_strength = 0,
                        partitions_per_tree = 1L,
                        weight_exponent     = 1L,
                        size_influence      = 2,
                        max_iter = 10,
                        verbose  = TRUE,
                        seed     = NULL) {
        
        
        tuning_strategy        <- match.arg(tuning_strategy)
        tuning_metric          <- match.arg(tuning_metric)
        
        stopifnot(inherits(graph, "igraph"))
        stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L,
                  prior_strength >= 0)
        
        # --- Exact prior reproduction short-circuit ---- *---- *---- *---- *---- *---- *
        if (is.infinite(prior_strength)) {
                if (is.null(prior_typology)) {
                        stop("`prior_strength = Inf` requires `prior_typology`.", call. = FALSE)
                }
                pr <- .exact_prior_result(.resolve_prior_typology(graph, prior_typology), verbose = TRUE)
                K <- ncol(pr$memberships)
                return(.new_pulse_regions(
                    fit           = NULL,
                    memberships   = pr$memberships,
                    hard_clusters = pr$hard_clusters,
                    coverage      = pr$coverage,
                    ensemble      = NULL,
                    config        = list(intermediate_regions = K,
                                         weight_exponent      = NA_real_,
                                         partitions_per_tree  = NA_integer_,
                                         n_st                 = NA_integer_,
                                         final_regions        = K),
                    prior_strength = Inf, tuned = FALSE))
        }
        
        n_st_values   <- sort(unique(as.integer(n_st)))
        intermediate_regions <- sort(unique(as.integer(intermediate_regions)))
        final_regions  <- sort(unique(as.integer(final_regions)))
        
        # Upper-bound check before doing any sampling work.
        if (any(intermediate_regions > igraph::vcount(graph))) {
                stop("`intermediate_regions` cannot exceed the number of vertices (", igraph::vcount(graph), ").", call. = FALSE)
        }
        
        # --- Sample the tree ensemble (largest requested n_st) 
        if (is.null(trees)){
                trees <- sample_spanning_trees(graph, n = max(n_st_values), seed = seed, verbose = verbose, prior_typology = prior_typology)   
        }

        # lambda_E > 0 without a prior is a silent no-op; warn the user.
        if (prior_strength > 0 && is.null(trees$prior_labels)) {
                warning("`prior_strength` > 0 but no `prior_typology` supplied; ",
                        "prior weighting has no effect.", call. = FALSE)
        }
        
        needs_tuning <-
                length(intermediate_regions) > 1L || length(final_regions) > 1L ||
                length(n_st) > 1L || length(partitions_per_tree) > 1L ||
                length(weight_exponent) > 1L
        
        if (needs_tuning) {
                print("Start Tuning")
            tr <- tune_regions(
                graph = graph, trees = trees, n_st = n_st_values,
                intermediate_regions = intermediate_regions,
                final_regions = final_regions, fuzziness = fuzziness,
                crisp = crisp, strategy = tuning_strategy,
                metric = tuning_metric, stability_B = stability_B,
                prior_strength = prior_strength,
                partitions_per_tree = partitions_per_tree,
                weight_exponent = weight_exponent,
                max_iter = max_iter, verbose = verbose, seed = seed, 
                size_influence       = size_influence
            )
            return(.new_pulse_regions(
                fit               = tr$best_result,
                config            = tr$best_config,
                trees             = trees,
                prior_strength    = prior_strength,
                tuned             = TRUE,
                score             = tr$best_score,
                metric            = tr$metric,
                strategy          = tr$strategy,
                tuning_log        = tr$tuning_log,
                candidates        = tr$candidates,
                n_ensembles_built = tr$n_ensembles_built,
                ensemble          = tr$ensemble))
        } else {
            message("Single configuration; no tuning.")
            ensemble <- compute_ensemble_memberships(
                trees = trees, n_st = n_st_values, graph = graph,
                intermediate_regions = intermediate_regions,
                partitions_per_tree  = partitions_per_tree,
                weight_exponent      = weight_exponent,
                prior_strength       = prior_strength, 
                size_influence       = size_influence,
                verbose = verbose)
            
            fit <- cluster_consensus(
                ensemble, k = final_regions, graph = graph,
                fuzziness = fuzziness, crisp = crisp,
                large_n_threshold = 20000, verbose = verbose)
            
            .new_pulse_regions(
                fit    = fit,
                ensemble = fit$ensemble_meta$memb_mat,
                config   = list(intermediate_regions = intermediate_regions,
                              weight_exponent      = weight_exponent,
                              partitions_per_tree  = partitions_per_tree,
                              n_st                 = n_st_values,
                              final_regions        = final_regions),
                trees          = trees,
                prior_strength = prior_strength,
                tuned          = FALSE)          # score/metric/strategy stay NA
        }
}


.enforce_contiguity <- function(graph, hard_clusters, d_to_medoids, medoid_idx) {
        n   <- length(hard_clusters)
        ep  <- igraph::ends(graph, igraph::E(graph), names = FALSE)
        lab <- rep(NA_integer_, n)
        
        for (r in seq_along(medoid_idx)) {
                in_r <- which(hard_clusters == r)
                cmp  <- igraph::components(igraph::induced_subgraph(graph, in_r))$membership
                own  <- match(medoid_idx[r], in_r)
                if (is.na(own)) next
                lab[in_r[cmp == cmp[own]]] <- r
        }
        
        a <- ep[, 1]; b <- ep[, 2]
        repeat {
                free <- is.na(lab)
                if (!any(free)) break
                f1 <- free[a] & !free[b]
                f2 <- free[b] & !free[a]
                node <- c(a[f1], b[f2])
                if (!length(node)) { lab[free] <- hard_clusters[free]; break }
                reg <- c(lab[b[f1]], lab[a[f2]])
                o   <- order(d_to_medoids[cbind(node, reg)])
                node <- node[o]; reg <- reg[o]
                keep <- !duplicated(node)
                lab[node[keep]] <- reg[keep]
        }
        lab
}







#' Knee-point detection via maximum perpendicular distance
#'
#' Locates the point on a curve furthest from the straight line joining its
#' endpoints (after `[0, 1]` normalisation of both axes). Robust to degenerate
#' inputs (no spread, coincident endpoints) — falls back to the last x.
#'
#' @param x,y Numeric vectors of equal length.
#' @return The \code{x} value at the knee.
#' @keywords internal
#' @noRd
.find_knee <- function(x, y) {
        n <- length(x)
        if (n <= 2L) return(x[n])
        
        rx <- max(x) - min(x)
        ry <- max(y) - min(y)
        if (rx <= 0) return(x[n])                      # no spread in x
        
        x_norm <- (x - min(x)) / rx
        y_norm <- (y - min(y)) / (ry + .Machine$double.eps)
        
        p1 <- c(x_norm[1L], y_norm[1L])
        p2 <- c(x_norm[n],  y_norm[n])
        line_len <- sqrt(sum((p2 - p1)^2))
        if (line_len <= .Machine$double.eps) return(x[n])  # coincident endpoints
        
        distances <- vapply(seq_len(n), function(i) {
                abs((x_norm[i] - p1[1L]) * (p2[2L] - p1[2L]) -
                            (y_norm[i] - p1[2L]) * (p2[1L] - p1[1L])) / line_len
        }, numeric(1))
        
        x[which.max(distances)]
}

#' Internal: compute the N x (n_st * partitions_per_tree) ensemble membership matrix
#'
#' Shared by the public \code{compute_ensemble_memberships()} and the stability
#' evaluator. For each tree, removes \code{intermediate_regions - 1} edges
#' (weighted sampling without replacement) and labels the resulting connected
#' components. Each tree is cut \code{partitions_per_tree} times independently;
#' every cut-set contributes one column.
#'
#' @param graph Parent \code{igraph} (used only to re-resolve edges when the
#'   stored tree object is a subgraph rather than an edge sequence).
#' @param trees A \code{"spanning_trees"} object.
#' @param n_st Number of trees to use.
#' @param intermediate_regions Number of regions per cut-set.
#' @param prior_strength Within-prior halving rate.
#' @param partitions_per_tree Integer >= 1. Independent cut-sets per tree.
#' @param weight_exponent Non-negative numeric. Cut probability is proportional
#'   to (prior-adjusted edge weight)^weight_exponent. 1 = proportional (default).
#' @return Integer matrix N x (n_st * partitions_per_tree) of region labels.
#' @keywords internal
#' @noRd
.compute_ensemble_memberships_internal <- function(graph,
                                                   trees,
                                                   n_st,
                                                   intermediate_regions,
                                                   prior_strength      = 0,
                                                   partitions_per_tree = 1L,
                                                   weight_exponent     = 1,
                                                   size_influence = 2) {
        n_use   <- as.integer(n_st)
        n_reg   <- as.integer(intermediate_regions)
        n_part  <- as.integer(partitions_per_tree)
        n_v     <- trees$n_vertices
        #weights <- lapply(trees$trees, function(g) igraph::E(g)$weight)
        #ep_all  <- trees$edge_endpoints
        wpe    <- trees$within_prior_edge
        if (n_reg > n_v) {
                stop(sprintf("intermediate_regions (%d) cannot exceed the number of vertices (%d).",
                             n_reg, n_v), call. = FALSE)
        }
        
        # One cut per extra region (a tree on n_v vertices has n_v - 1 edges).
        n_to_cut <- n_reg - 1L
        
        memb_mat <- matrix(0, nrow = n_v, ncol = n_use * n_part)
        col <- 0
        for (j in seq_len(n_use)) {
                tree <- trees$trees[[j]]
                if (igraph::vcount(tree) != n_v){
                        stop("Spanning tree does not span the parent vertex set", call. = FALSE)
                }
                w   <- igraph::E(tree)$weight      
                n_e <- length(w)   
                edge_ids <- igraph::E(tree)$eid
                n_c    <- .subtree_sizes(edge_ids, trees$edge_endpoints, n_v = n_v)
                size_f <- n_c * (n_v - n_c) / n_v
                if (any(is.na(w)) || any(is.infinite(w)) || any(w < 0)){
                        stop("Edge weights must be positive and finite numbers", .call = FALSE)
                }
                if (n_to_cut > n_v){
                        stop(
                                sprintf(
                                        "Requested %d cuts but the tree onlt has %d edges",
                                        n_to_cut, n_e
                                ) , .call = FALSE
                                )
                }        
                edge_ids <- igraph::E(tree)$eid
                if (anyNA(edge_ids)) {
                        stop("Failed to map a spanning-tree edge back to the parent graph.",
                             call. = FALSE)
                }
                if (n_to_cut >= length(edge_ids)) {
                        stop(sprintf(
                                paste0("Requested %d cuts but the spanning tree has only %d edges; ",
                                       "reduce intermediate_regions or supply a connected graph."),
                                n_to_cut, length(edge_ids)
                        ), call. = FALSE)
                }
                

                ## Normalize within the tree (heaviest edge -> 1) before exponentiating.
                # sample.int() (further below) rescales what you pass to prop internally,
                # so dividing by a constant leaves the sampling distribution 
                # unchanged and prevents over- or underflow errors when weight_exponent 
                # or weights are large/small. Wrapping in if clause keep legacy applications intact.
                w_n <- w / max(w)
                
                ## size_influence is the ratio of the size term's spread to the
                ## environmental term's spread, both on the log scale the exponential
                ## race competes on. Solving for the exponent per tree makes the
                ## argument mean the same thing whatever the graph size or the weight
                ## distribution: 0 = size-blind, 1 = balance and environment contribute
                ## equally, 2 = balance dominates two to one.
                lw <- weight_exponent * log(w_n)
                ls <- log(size_f)
                sd_w <- stats::sd(lw)
                sd_s <- stats::sd(ls)
                alpha <- if (sd_s > 0 && sd_w > 0) size_influence * sd_w / sd_s else 0
                w <- exp(lw + alpha * ls)

                # if (weight_exponent != 1){
                #         w = size_f * ((w / max(w))^weight_exponent)
                # } else {
                #         w <- size_f * (w)
                # }
                        
                if (prior_strength != 0 && !is.na(wpe) && !is.null(wpe))
                        #w <- .apply_prior_weighting(w, trees$within_prior_edge, prior_strength)
                        w <- w / 2^(prior_strength * wpe[igraph::E(tree)$eid])
                # partitions_per_tree independent cut-sets of the SAME tree; each is an
                # independent draw contributing one ensemble column.
                for (p in seq_len(n_part)) {
                        col        <- col + 1
                        # fix 29.07.26 exponential race
                        key <- stats::rexp(n_e) / w  
                        cut <- which(key <= sort(key, partial = n_to_cut)[n_to_cut])
                        if (length(cut) > n_to_cut) cut <- cut[seq_len(n_to_cut)]
                        memb_mat[, col] <- as.integer(
                                igraph::components(igraph::delete_edges(tree, cut))$membership)
                        # old slower version 
                        # cut_idx    <- sample.int(length(edge_ids), n_to_cut,
                        #                          prob = w, replace = FALSE)
                        # kept_edges <- edge_ids[-cut_idx]
                        # ep <- ep_all[kept_edges, , drop = FALSE]
                        # g  <- igraph::make_graph(as.vector(t(ep)), n = n_v, directed = FALSE)
                        # memb_mat[, col] <- igraph::components(g)$membership
                }
        }
        # A matrix with n_st columns and one row for each polygon
        memb_mat
}



# ## --- Prior Weighting ---- *
# # Prior-typology weighting (no-op when no prior or prior_strength == 0). This
# # acts on the exponentiated dissimilarities (just above).
# # Prior_strength is a fixed number of halvings of the cut probability
# # regardless of weight_exponent.
# if (prior_strength != 0)
#         weights <- .apply_prior_weighting(weights, trees$within_prior_edge, prior_strength)
# 
# 
# total_cols <- n_use * n_part

#' Match spanning-tree edges to parent-graph edge indices
#'
#' Builds canonical undirected edge keys (smaller vertex first) and matches.
#'
#' @param tree_ep,parent_ep Two-column integer matrices of edge endpoints.
#' @return Integer vector of parent-edge indices (NA where unmatched).
#' @keywords internal
#' @noRd
.match_tree_edges <- function(tree_ep, parent_ep) {
        make_key <- function(ep) {
                paste(pmin(ep[, 1], ep[, 2]), pmax(ep[, 1], ep[, 2]), sep = "-")
        }
        match(make_key(tree_ep), make_key(parent_ep))
}

#' Adjusted Rand Index (Hubert & Arabie 1985)
#'
#' Pure-base implementation. Permutation-invariant; ARI = 1 for identical
#' partitions (up to label permutation), expected ARI = 0 under random
#' assignment.
#'
#' @param a,b Integer or character label vectors of equal length.
#' @return Numeric scalar in approximately \code{[-0.5, 1]}; \code{NA_real_}
#'   if input length < 2.
#'
#' @details Used by \code{.stability_fuzzy()} to score how well two hardened
#'   partitions agree on the overlap of resampled objects. We use the crisp ARI
#'   on argmax labels rather than a fuzzy-ARI generalisation because (i)
#'   stability of *structure* is what matters here, not stability of *softness*,
#'   and (ii) it has no dependencies. If a non-trivial fraction of objects have
#'   a top-two membership margin below ~0.1, consider substituting a
#'   Hüllermeier–Rifqi or Campello fuzzy Rand index (not implemented).
#'
#' @keywords internal
#' @noRd
.adjusted_rand_index <- function(a, b) {
        stopifnot(length(a) == length(b))
        if (length(a) < 2) return(NA_real_)
        tab <- table(a, b)
        n   <- sum(tab)
        if (n < 2) return(NA_real_)
        sum_nij <- sum(choose(tab, 2))
        ai <- rowSums(tab); bj <- colSums(tab)
        sum_ai <- sum(choose(ai, 2))
        sum_bj <- sum(choose(bj, 2))
        expected <- (sum_ai * sum_bj) / choose(n, 2)
        max_idx  <- 0.5 * (sum_ai + sum_bj)
        denom    <- max_idx - expected
        # Degenerate case: both partitions are singletons or a single block.
        if (denom == 0) return(1)
        (sum_nij - expected) / denom
}

#' Stratified subsample by hard-cluster labels
#'
#' Proportional allocation with a hard floor of 2 objects per cluster, so that
#' silhouette / ARI computations on the subsample are not pathological for small
#' clusters.
#'
#' @param hard_clusters Integer vector of length N.
#' @param n_target Approximate total subsample size.
#' @return Sorted integer vector of selected row indices.
#' @keywords internal
#' @noRd
.stratified_subsample <- function(hard_clusters, n_target) {
        N <- length(hard_clusters)
        if (n_target >= N) return(seq_len(N))
        tab   <- table(hard_clusters)
        props <- as.numeric(tab) / N
        n_per <- pmax(2L, as.integer(round(props * n_target)))
        names(n_per) <- names(tab)
        groups <- names(tab)
        idx <- unlist(lapply(groups, function(g) {
                in_g <- which(as.character(hard_clusters) == g)
                take <- min(length(in_g), n_per[[g]])
                if (length(in_g) <= take) in_g else sample(in_g, take)
        }), use.names = FALSE)
        sort(idx)
}

#' Crisp Rousseeuw silhouette from a dissimilarity matrix
#'
#' Vectorised over the cluster axis; O(N k) after the N x N matrix is
#' materialised. Intended for a subsample of moderate size (a few thousand
#' rows), not the full graph.
#'
#' @param D Dissimilarity, either a \code{dist} object or a full N x N matrix.
#' @param labels Integer vector of cluster labels of length N.
#' @return Numeric vector of per-object silhouette widths.
#' @keywords internal
#' @noRd
.rousseeuw_silhouette <- function(D, labels) {
        if (inherits(D, "dist")) D <- as.matrix(D)
        N <- length(labels)
        clusters <- sort(unique(labels))
        K <- length(clusters)
        if (K < 2L) return(rep(0, N))
        cluster_sizes <- as.integer(table(factor(labels, levels = clusters)))
        
        # N x K matrix of summed distances from each point to each cluster.
        sums <- matrix(0, nrow = N, ncol = K)
        for (j in seq_len(K)) {
                in_j <- labels == clusters[j]
                sums[, j] <- rowSums(D[, in_j, drop = FALSE])
        }
        
        own_idx <- match(labels, clusters)
        s   <- numeric(N)
        eps <- .Machine$double.eps
        for (i in seq_len(N)) {
                o <- own_idx[i]
                if (cluster_sizes[o] <= 1L) { s[i] <- 0; next }
                means <- sums[i, ] / cluster_sizes
                # Own-cluster mean excludes self (self-distance is 0, so subtract nothing
                # from the sum but divide by size - 1).
                means[o] <- sums[i, o] / (cluster_sizes[o] - 1L)
                a_i <- means[o]
                b_i <- min(means[-o])
                s[i] <- (b_i - a_i) / max(a_i, b_i, eps)
        }
        s
}
#' Sample spanning trees from a Spatial Graph
#'
#' Draws \code{n} spanning trees from a weighted, (ideally) connected
#' \code{igraph} object. The trees, together with the parent graph's edge
#' weights/endpoints and any resolved prior-typology information, are bundled
#' into a \code{"spanning_trees"} object that the rest of the pipeline consumes.
#'
#' @param graph An \code{igraph} object. Edges must carry a numeric
#'   \code{weight} attribute (dissimilarity / cost; finite and non-negative).
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
#'   non-NA label are flagged as "within-prior" and become eligible for
#'   weight halving in the SKATER cutting step (see \code{prior_strength}
#'   in \code{\link{compute_ensemble_memberships}} and
#'   \code{\link{get_regions}}). Default \code{NULL} (no prior).
#'@param target_eta numeric. How many edges of the MST should be altered on 
#'   average in each iteration.
#'   
#' @return An object of class \code{"spanning_trees"} — a list with:
#'   \describe{
#'     \item{\code{trees}}{List of sampled spanning-tree objects (igraph edge
#'       sequences or graphs, depending on the installed igraph version).}
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
#'   }
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
#' @importFrom igraph sample_spanning_tree E ends vcount vertex_attr vertex_attr_names edge_attr components
#' @importFrom stats mad rnorm
#' @export
sample_spanning_trees <- function(graph, n, seed = NULL, verbose = TRUE,
                                  prior_typology = FALSE, target_eta = .1) {
        
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
        
        # Extract properties of weights later used for perturbation
        
        # A disconnected graph yields a spanning *forest* (n_v - n_components edges),
        # which breaks the "cut n_reg-1 edges of a tree" arithmetic downstream.
        comp <- igraph::components(graph)
        if (comp$no > 1L) {
                warning(sprintf(
                        paste0("`graph` has %d connected components: sample_spanning_tree() will ",
                               "return a spanning forest and the intermediate region cut counts assume a single ",
                               "tree. Supply a connected graph for well-defined region counts."),
                        comp$no
                ), call. = FALSE)
        }
        
        # --- RNG handling: seed locally, restore the caller's state on exit ---- *-
        if (!is.null(seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                set.seed(seed)
        }
        edge_endpoints    <- igraph::ends(graph, igraph::E(graph), names = FALSE)
        # --- Resolve prior typology (if supplied) ---- *---- *---- *---- *---- *---- *--
        if (prior_typology){
                prior_typology    <- "prior"
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
        }
        
        
        # --- Sample the trees ---- *---- *---- *---- *---- *---- *---- *---- *---- *---- *--
        if (verbose) message(sprintf("Sampling %d spanning trees ...", n))
        tick <- if (verbose) .make_progress(n, "Spanning trees") else function(i) NULL
        
        .turnover <- function(t0_ids, t1_ids) 1 - length(intersect(t0_ids, t1_ids)) / length(t0_ids)
        .calibrate_eta_empirical <- function(graph,
                                             w,
                                             s,
                                             target,
                                             n_probe = 20,
                                             tol = 0.02,
                                             max_iter = 12) {
                # if target is 0 hard code an eta of zero and exit 
                if (target == 0){
                        return(0)
                }
                
                t0 <- igraph::E(igraph::mst(graph, weights = w))$eid
                measure <- function(eta) {
                        mean(replicate(n_probe, {
                                wp <- w + eta * s * stats::rnorm(length(w))
                                t1 <- igraph::E(igraph::mst(graph, weights = pmax(wp, 0)))$eid
                                .turnover(t0, t1)
                        }))
                }
                lo <- 1e-4
                hi <- 1000
                
                for (i in seq_len(max_iter)) {
                        mid <- sqrt(lo * hi)
                        tm <- measure(mid)
                        if (abs(tm - target) < tol)
                                return(mid)
                        if (tm < target)
                                lo <- mid
                        else
                                hi <- mid
                }
        }
        
        eta_emp <- .calibrate_eta_empirical(graph, 
                                            w = edge_w, 
                                            s = stats::mad(edge_w),
                                            target = target_eta)
        tree_list <- vector("list", n)
        
        for (j in seq_len(n)) {
                # Perturb weights 
                wp <- edge_w + eta_emp * mad(edge_w) * rnorm(length(edge_w))
                tree_list[[j]] <- igraph::mst(graph, weights = wp)
                tick(j)
        }
        
        structure(
                list(
                        trees             = tree_list,
                        n_trees           = as.integer(n),
                        n_vertices        = igraph::vcount(graph),
                        edge_weights      = edge_w,
                        edge_endpoints    = edge_endpoints,
                        seed              = seed,
                        prior_labels      = ifelse(exists("prior_labels"), prior_labels, NA),
                        within_prior_edge = ifelse(exists("within_prior_edge"), within_prior_edge, NA)
                ),
                class = "spanning_trees"
        )
}



#' Bootstrap stability of a fuzzy partition (mean ARI vs. reference)
#'
#' On each of \code{B} resamples, draws a row subsample of size
#' \code{floor(N * subsample_frac)} WITHOUT replacement, re-runs CLARA/PAM on
#' the subsampled ensemble matrix to obtain a hardened partition of the
#' subsample, and computes the ARI against the reference labels restricted to
#' those rows. Returns the mean ARI.
#'
#' @details This is the Hennig-style "compare each replicate to the full-data
#'   reference" recipe (one ARI per replicate), not pairwise replicate-vs-
#'   replicate (B^2 ARIs). It is B-fold cheaper and is what
#'   \code{fpc::clusterboot} does by default. Sampling is without replacement so
#'   ARI's contingency table is well-defined.
#'
#' @param memb_mat Integer matrix N x n_trees (the ensemble matrix).
#' @param hard_clusters Reference hardened labels of length N (from the
#'   full-data fit).
#' @param k Number of clusters.
#' @param fuzziness Fuzzifier; unused for hardened-label stability but kept in
#'   the signature so future fuzzy-ARI extensions need no API break.
#' @param B Number of bootstrap replicates. \code{0} disables.
#' @param subsample_frac Fraction of rows to retain per replicate. Default 0.8.
#' @param large_n_threshold CLARA/PAM cutoff (mirrors
#'   \code{cluster_consensus}).
#' @param clara_samples,clara_sampsize CLARA controls.
#' @param seed Optional integer seed; the caller's RNG state is restored on
#'   exit.
#' @param verbose Logical.
#' @return Numeric scalar mean ARI in approximately \code{[-0.5, 1]}, or
#'   \code{NA_real_} if \code{B = 0} or the subsample is too small.
#' @keywords internal
#' @noRd
.stability_fuzzy <- function(memb_mat,
                             hard_clusters,
                             k,
                             fuzziness = 2,
                             B = 25L,
                             subsample_frac = 0.8,
                             large_n_threshold = 5,
                             clara_samples = 20L,
                             clara_sampsize = NULL,
                             seed = NULL,
                             verbose = FALSE) {
        if (B <= 0L) return(NA_real_)
        N     <- nrow(memb_mat)
        n_sub <- floor(N * subsample_frac)
        if (n_sub < k * 2L) {
                if (verbose) message("    [stability] subsample too small; returning NA")
                return(NA_real_)
        }
        
        # RNG handling: seed locally, restore the caller's state on exit.
        if (!is.null(seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                set.seed(seed)
        }
        
        aris <- numeric(B)
        for (b in seq_len(B)) {
                sub_idx  <- sort(sample.int(N, size = n_sub, replace = FALSE))
                memb_sub <- memb_mat[sub_idx, , drop = FALSE]
                ref_sub  <- hard_clusters[sub_idx]
                if (nrow(memb_sub) > large_n_threshold) {
                        
                        sampsize <- 
                                if (is.null(clara_sampsize)) {
                                min(nrow(memb_sub), 200L + 2L * k)
                                } else clara_sampsize
                        
                        # CRUD: exchanged with code below 08.07.2026
                        # fit_b <- cluster::clara(memb_sub, k = k,
                        #                         metric   = "manhattan",
                        #                         samples  = clara_samples,
                        #                         sampsize = sampsize,
                        #                         pamLike  = TRUE)
                        
                        fit_b <- .clara_hamming(memb_sub, 
                                                k = k,
                                                 samples  = clara_samples,
                                                 sampsize = sampsize)
                        
                        new_lab <- fit_b$clustering
                } else {
                        coassoc_b <- parallelDist::parallelDist(memb_sub, method = "hamming")
                        fit_b     <- cluster::pam(coassoc_b, k = k)
                        new_lab   <- fit_b$clustering
                }
                aris[b] <- .adjusted_rand_index(new_lab, ref_sub)
        }
        mean(aris, na.rm = TRUE)
}

#' Derive the \code{score} column of \code{tuning_log} from sub-metric columns
#'
#' Arranged so that \code{which.max(score)} always selects the best row:
#' minimise-is-better indices are negated, and Borda rank-sums are negated. For
#' metrics whose primary value lives directly in the log (e.g. SIL_F),
#' \code{score} is a copy.
#'
#' @param df Data frame containing the CVI columns.
#' @param metric Character metric name.
#' @return Numeric vector of scores (higher = better).
#' @keywords internal
#' @noRd
.derive_score_column <- function(df, metric) {
        
        switch(metric,
               #R2_JS     = df$R2_JS,
               CH_JS     = df$CH_D,
               #FS        = -df$FS,
              # GD5       = df$GD5,
               stability = df$STAB,
              # silhouette = df$SILH_HARD,
               silhouette_fuzzy = df$SIL_F,
               borda     = -.borda_score(df[, c("CH_D","STAB", "SIL_F")],
                                         c("max", "max", "max")),
               stop("Unknown metric: ", metric, call. = FALSE))
}

#' Compute Hamming distances from all nodes to reference (medoid) nodes
#'
#' Returns an N x k matrix giving the fraction of ensemble columns in which each
#' node disagrees with each reference node.
#'
#' @param memb_mat Integer matrix (N x n_trees).
#' @param ref_idx Integer vector of length k (row indices of medoids).
#' @return Numeric matrix (N x k) of Hamming distances (fractions in `[0, 1]`).
#'
#' @details Accumulates per-tree disagreements into the N x k counter, iterating
#'   over the n_trees columns. This keeps peak memory at O(N k) instead of
#'   materialising an N x n_trees logical matrix for every medoid — important
#'   when N reaches 1e5-1e6.
#'
#' @keywords internal
#' @noRd
.hamming_to_refs <- function(memb_mat, ref_idx) {
        n_trees     <- ncol(memb_mat)
        N           <- nrow(memb_mat)
        k           <- length(ref_idx)
        medoid_rows <- memb_mat[ref_idx, , drop = FALSE]   # k x n_trees
        
        out <- matrix(0, nrow = N, ncol = k)
        for (t in seq_len(n_trees)) {
                # outer(col, medoid_col, "!=") is an N x k logical, coerced to numeric in +.
                out <- out + outer(memb_mat[, t], medoid_rows[, t], FUN = "!=")
        }
        out / n_trees
}

#' Resolve a \code{prior_typology} argument to a vertex-label vector
#'
#' Accepts either a single character (vertex-attribute name on \code{graph}) or
#' a vector of length \code{vcount(graph)}. Coerces to character; empty strings
#' become NA. Returns NULL if input is NULL.
#'
#' @param graph An \code{igraph} object.
#' @param prior_typology See \code{\link{sample_spanning_trees}}.
#' @return Character vector of length \code{vcount(graph)}, or \code{NULL}.
#' @keywords internal
#' @importFrom igraph vcount vertex_attr vertex_attr_names
#' @noRd
.resolve_prior_typology <- function(graph, prior_typology) {
        if (is.null(prior_typology)) return(NULL)
        
        n_v <- igraph::vcount(graph)
        
        if (is.character(prior_typology) && length(prior_typology) == 1L) {
                if (!prior_typology %in% igraph::vertex_attr_names(graph)) {
                        stop(sprintf("Vertex attribute '%s' not found on graph. ", prior_typology),
                             "Attach the typology via V(graph)$", prior_typology,
                             " <- ... before calling.", call. = FALSE)
                }
                labs <- igraph::vertex_attr(graph, prior_typology)
        } else if (length(prior_typology) == n_v) {
                labs <- prior_typology
        } else {
                stop("`prior_typology` must be either a single character giving the name of ",
                     "a vertex attribute on `graph`, or a vector of length vcount(graph) = ",
                     n_v, ".", call. = FALSE)
        }
        
        labs <- as.character(labs)
        labs[!nzchar(labs)] <- NA_character_
        
        # Sanity check: at least 2 distinct non-NA labels for weighting to do anything.
        uniq <- unique(labs[!is.na(labs)])
        if (length(uniq) < 2L) {
                warning("Prior typology contains fewer than 2 distinct non-NA labels; ",
                        "weighting will have no effect.", call. = FALSE)
        }
        labs
}

#' Flag edges whose endpoints share a non-NA prior label
#'
#' @param prior_labels Character vector of vertex labels, or \code{NULL}.
#' @param edge_endpoints Two-column integer matrix of edge endpoints.
#' @return Logical vector of length \code{nrow(edge_endpoints)}, or \code{NULL}.
#' @keywords internal
#' @noRd
.compute_within_prior_edge <- function(prior_labels, edge_endpoints) {
        if (is.null(prior_labels)) return(NULL)
        a <- prior_labels[edge_endpoints[, 1]]
        b <- prior_labels[edge_endpoints[, 2]]
        !is.na(a) & !is.na(b) & a == b
}

#' Apply the within-prior halving rule to edge weights
#'
#' \eqn{\tilde{w}(i,j) = w(i,j) / 2^{\lambda_E}} when both endpoints share a
#' non-NA prior label; \eqn{\tilde{w}(i,j) = w(i,j)} otherwise. Thus
#' \eqn{\lambda_E = 0} returns weights unchanged; \eqn{\lambda_E = 1} halves
#' within-prior weights; \eqn{\lambda_E = 2} quarters them, etc. Lowering these
#' weights makes the corresponding edges less likely to be cut, so prior regions
#' tend to stay intact.
#'
#' @param weights Numeric edge-weight vector.
#' @param within_prior_edge Logical vector (same length), or \code{NULL}.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda_E}.
#' @return Numeric vector of (possibly down-weighted) edge weights.
#' @keywords internal
#' @noRd
.apply_prior_weighting <- function(weights, within_prior_edge, prior_strength) {
        if (is.na(within_prior_edge) || prior_strength == 0) return(weights)
        weights <- lapply(weights, 
                          function(x) x / (2 ^ (prior_strength * as.numeric(within_prior_edge))
                                  )
        )
}

#' Exact reproduction of a prior typology (\code{prior_strength = Inf})
#'
#' Returns the prior labels directly as a (degenerate) clustering result, with a
#' length-N \code{coverage} flag marking vertices that have no prior label.
#' Uncovered vertices receive \code{NA} hard labels and all-zero membership rows.
#'
#' @param prior_labels Character vector of length N (resolved labels).
#' @param verbose Logical.
#' @return A list mirroring the \code{get_regions()} return, plus \code{coverage}
#'   and \code{class_levels}.
#' @keywords internal
#' @noRd
.exact_prior_result <- function(prior_labels, verbose = TRUE) {
        N         <- length(prior_labels)
        coverage  <- !is.na(prior_labels)
        n_covered <- sum(coverage)
        n_dropped <- N - n_covered
        
        if (n_covered == 0L) {
                stop("`prior_strength = Inf` but no vertices have a prior label.",
                     call. = FALSE)
        }
        
        if (verbose) {
                message("prior_strength = Inf: returning prior typology directly.")
                if (n_dropped > 0L) {
                        message(sprintf(
                                "  %d / %d vertices have no prior label and are marked NA (coverage = FALSE).",
                                n_dropped, N
                        ))
                }
        }
        
        classes <- sort(unique(prior_labels[coverage]))
        K       <- length(classes)
        
        hard_clusters           <- rep(NA_integer_, N)
        hard_clusters[coverage] <- match(prior_labels[coverage], classes)
        
        # One-hot membership for covered vertices; all-zero rows for uncovered ones.
        memberships <- matrix(0, nrow = N, ncol = K)
        memberships[cbind(which(coverage), hard_clusters[coverage])] <- 1
        
        list(
                memberships         = memberships,
                hard_clusters       = hard_clusters,
                coverage            = coverage,        # length-N logical
                best_intermediate_regions = K,
                best_final_regions  = K,
                best_n_st          = NA_integer_,
                best_score          = NA_real_,
                tuning_log          = NULL,
                #n_st_stability     = NULL,
                trees               = NULL,
                cluster_result      = NULL,
                prior_strength      = Inf,
                class_levels        = classes
        )
}

## TODO ???
.region_step_distance <- function(graph, lab, k) {
        n <- igraph::vcount(graph)
        out <- matrix(Inf, n, k)
        for (j in seq_len(k)) {
                src <- which(lab == j)
                if (!length(src)) next
                g2 <- igraph::add_edges(igraph::add_vertices(graph, 1),
                                        as.vector(rbind(n + 1L, src)))
                d <- igraph::distances(g2, v = n + 1L, weights = NA, mode = "all")
                out[, j] <- d[1, seq_len(n)] - 1
        }
        out
}

.fuzzy_memberships <- function(d_to_medoids,
                               G = G,
                               fuzziness = fuzziness,
                               h = 20
) {
        n <- nrow(d_to_medoids)
        k <- ncol(d_to_medoids)
        stopifnot(fuzziness > 1, h >= 0)
        attenuate <- !is.null(G) && is.finite(h)
        if (attenuate && (nrow(G) != n || ncol(G) != k))
                stop("G must have the same dimensions as d_to_medoids.", call. = FALSE)
        
        # fuzzyness exponent
        p    <- -2 / (fuzziness - 1)
        # Inf where d == 0
        logu <- p * log(d_to_medoids)
        
        if (attenuate) {
                # Inf where a region is unreachable within h steps
                pen <- G / h                                
                logu[!is.infinite(logu)] <- 
                        logu[!is.infinite(logu)] - pen[!is.infinite(logu)]
        }
        
        # Membership matrix
        U <- matrix(0, n, k)
        
        # Rows containing a zero distance: the node coincides with one 
        # or more medoids. Among those, take the nearest in graph steps;
        #  ties split evenly.
        zr <- which(apply(is.infinite(logu) & logu > 0, 1L, any))
        if (length(zr)) {
                for (i in zr) {
                        cand <- which(is.infinite(logu[i, ]) & logu[i, ] > 0)
                        if (attenuate &&
                            length(cand) > 1L)
                                cand <- cand[G[i, cand] == min(G[i, cand])]
                        U[i, cand] <- 1 / length(cand)
                }
        }
        
        ok <- setdiff(seq_len(n), zr)
        if (length(ok)) {
                L  <- logu[ok, , drop = FALSE]
                mx <- do.call(pmax, c(as.data.frame(L), na.rm = TRUE))
                E  <- exp(L - mx)                            # rowmax is now exactly 0
                E[!is.finite(E)] <- 0
                s  <- rowSums(E)
                bad <- !is.finite(s) |
                        s <= 0                # every region unreachable
                if (any(bad)) {
                        E[bad, ] <- exp(p * log(d_to_medoids[ok[bad], , drop = FALSE]) -
                                                max(p * log(d_to_medoids[ok[bad], , drop = FALSE])))
                        s[bad] <- rowSums(E[bad, , drop = FALSE])
                }
                U[ok, ] <- E / s
        }
        U
}


#' @rdname compute_ensemble_memberships
#' @param x An object of class \code{"ensemble_memberships"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.ensemble_memberships} returns \code{x} invisibly.
#' @export
print.ensemble_memberships <- function(x, ...) {
        ppt <- if (is.null(x$partitions_per_tree)) 1L else x$partitions_per_tree
        if (ppt > 1L) {
                cat(sprintf(
                        "<ensemble_memberships>  %d vertices x %d members (%d trees x %d partitions) | intermediate_regions = %d\n",
                        x$n_vertices, ncol(x$memberships), x$n_st, ppt, x$intermediate_regions
                ))
        } else {
                cat(sprintf(
                        "<ensemble_memberships>  %d vertices x %d trees | intermediate_regions = %d\n",
                        x$n_vertices, x$n_st, x$intermediate_regions
                ))
        }
        invisible(x)
}


#' @rdname cluster_consensus
#' @param x An object of class \code{"fuzzy_clusters"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.fuzzy_clusters} returns \code{x} invisibly.
#' @export
print.fuzzy_clusters <- function(x, ...) {
        cat(sprintf(
                "<fuzzy_clusters>  %d nodes x %d clusters | method = %s | avg silhouette = %.4f\n",
                nrow(x$memberships), x$k, x$method, x$avg_sil_width
        ))
        invisible(x)
}


#' @export
print.pulse_regions <- function(x, ...) {
        cat(sprintf("<pulse_regions>  %d units x %d regions\n",
                    nrow(x$memberships), ncol(x$memberships)))
        m <- x$fit$ensemble_meta
        if (!is.null(m)) {
                cat(sprintf(
                        "  ensemble: %d members | intermediate_regions = %s | chance floor = %.4f (n_eff = %.1f)\n",
                        x$fit$n_members, format(m$intermediate_regions),
                        m$chance_coassoc, 1 / m$chance_coassoc))
        }
        cat(sprintf("  mean typicality (assigned region) = %.3f | mean entropy = %.3f\n",
                    mean(x$typicality[cbind(seq_len(nrow(x$typicality)),
                                            x$hard_clusters)]),
                    mean(x$entropy)))
        invisible(x)
}


#' @rdname tune_regions
#' @param x An object of class \code{"regions_tuning"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.regions_tuning} returns \code{x} invisibly.
#' @export
print.regions_tuning <- function(x, ...) {
        cat(sprintf("<regions_tuning>  %d configurations | %d ensembles built | metric = %s | strategy = %s\n",
                    nrow(x$tuning_log),
                    if (is.null(x$n_ensembles_built)) NA_integer_ else x$n_ensembles_built,
                    x$metric, x$strategy))
        cat(sprintf("  axes searched: %s\n",
                    if (length(x$active_axes)) paste(x$active_axes, collapse = ", ") else "none"))
        cat(sprintf("  best: ir = %d | gamma = %g | ppt = %d | n_st = %d | k = %d  (score = %.4f)\n",
                    x$best_intermediate_regions, x$best_weight_exponent,
                    x$best_partitions_per_tree, x$best_n_st,
                    x$best_final_regions, x$best_score))
        invisible(x)
}


#' @rdname sample_spanning_trees
#' @param x An object to print (e.g. of class \code{"spanning_trees"}).
#' @param ... Further arguments passed to or from other methods (unused).
#' @return \code{print.spanning_trees} returns \code{x} invisibly.
#' @export
print.spanning_trees <- function(x, ...) {
        cat(sprintf(
                "<spanning_trees>  %d trees | %d vertices | seed = %s\n",
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
.TUNE_AXES <- c("intermediate_regions", "weight_exponent",
                "partitions_per_tree", "n_st", "final_regions")

#' Ensemble-level axes (those that require recomputing the membership matrix)
#' @keywords internal
#' @noRd
.ENSEMBLE_AXES <- c("intermediate_regions", "weight_exponent", "partitions_per_tree")

#' Canonical key for one hyper-parameter configuration
#'
#' Numeric axes are formatted at full double precision so that floating-point
#' candidates round-trip exactly through the log.
#'
#' @param df Data frame with the five axis columns.
#' @return Character vector of length \code{nrow(df)}.
#' @keywords internal
#' @noRd
.tune_key <- function(df) {
        paste(as.integer(df$intermediate_regions),
              sprintf("%.17g", as.numeric(df$weight_exponent)),
              as.integer(df$partitions_per_tree),
              as.integer(df$n_st),
              as.integer(df$final_regions),
              sep = "|")
}

#' Tune SKATER-CON Hyperparameters
#'
#' Searches over up to five hyper-parameters -- \code{intermediate_regions},
#' \code{weight_exponent}, \code{partitions_per_tree}, \code{n_st} and
#' \code{final_regions} -- to find the combination optimising a
#' clustering-validity criterion, operating on a pre-computed
#' \code{"spanning_trees"} object. Any argument supplied as a length-1 vector is
#' held fixed; any argument supplied as a vector of length > 1 becomes a search
#' axis.
#'
#' @section Cost structure:
#' Only \code{final_regions} is cheap to vary: it reuses the cached ensemble
#' membership matrix. \code{n_st} is nearly free, because an ensemble of
#' \code{n_st} trees is an exact column-subset of an ensemble of more trees
#' drawn from the same \code{trees} object; the tuner therefore builds each
#' ensemble once at \code{max(n_st)} and slices columns. The remaining three
#' axes each require a fresh pass of
#' \code{\link{compute_ensemble_memberships}}, so the number of distinct
#' \code{(intermediate_regions, weight_exponent, partitions_per_tree)} triples
#' is the real cost driver, not the number of grid cells.
#'
#' Because ensembles are built at \code{max(n_st)} regardless of the incumbent,
#' supplying a vector \code{n_st} raises the cost of \emph{every} ensemble
#' build. This buys determinism (results do not depend on the order in which
#' configurations happen to be visited) and paired comparisons across
#' \code{n_st} (all ensemble sizes share the same trees and the same cut
#' realisations). When \code{n_st} is scalar there is no penalty.
#'
#' @section Caching:
#' A single-slot cache holds one ensemble matrix at a time, and evaluation is
#' ordered so that configurations sharing an ensemble are visited consecutively.
#' A keyed cache over a five-way grid would hold one N x (n_st * ppt) integer
#' matrix per ensemble triple, which at continental N is tens of gigabytes.
#' Peak memory here is one full matrix plus at most one column slice.
#'
#' @section Search strategies:
#' \code{"grid"} evaluates the full factorial. \code{"sequential"} runs one
#' coordinate-descent pass over the active axes in the order
#' \code{intermediate_regions}, \code{weight_exponent},
#' \code{partitions_per_tree}, \code{n_st}, \code{final_regions}, starting from
#' the median candidate on each axis. \code{"iterative"} repeats that pass until
#' the incumbent configuration stops moving or \code{max_iter} is reached.
#'
#' @section Caveats on the ensemble-size axes:
#' \code{n_st} and \code{partitions_per_tree} both act on the ensemble by adding
#' columns (the matrix has \code{n_st * partitions_per_tree} of them), so they
#' are partially confounded: what most validity indices respond to is the total
#' column count, through the resolution of the Hamming consensus distance
#' (quantized in units of 1 / n_columns) and through Monte Carlo noise. Most
#' indices therefore drift monotonically in ensemble size rather than exhibiting
#' an interior optimum, and a CVI-selected \code{n_st} will typically just be
#' the largest candidate. Tuning \code{partitions_per_tree} at fixed
#' total columns (trading it against \code{n_st}) is the better-posed version of
#' that question: it asks how ensemble diversity should be split between tree
#' topologies and cut-sets.
#'
#' @param graph The parent \code{igraph} object the trees were sampled from.
#' @param env Environmental data passed through to \code{.compute_fuzzy_cvis()}.
#'  The same data used to weight the edges in \code{\link{add_edge_weight}}.
#' @param env_id Name of any ID columns in `env`, which should not be used in 
#'   the computation of cluster validity metrics. If `NULL` (default) no columns
#'   are removed from `env`.
#' @param trees A \code{"spanning_trees"} object from
#'   \code{\link{sample_spanning_trees}}.
#' @param n_st Integer or integer vector (each \eqn{\ge 2} and
#'   \eqn{\le} \code{trees$n_trees}). Candidate ensemble size(s).
#' @param intermediate_regions Integer or integer vector (each \eqn{\ge 2}).
#'   Candidate number(s) of intermediate regions per cut-set. Default 50.
#' @param final_regions Integer or integer vector (each \eqn{\ge 2}). Candidate
#'   number(s) of final consensus clusters \emph{k}.
#' @param partitions_per_tree Integer or integer vector (each \eqn{\ge 1}).
#'   Candidate number(s) of independent cut-sets applied per tree. Default 1.
#' @param weight_exponent Numeric or numeric vector (each \eqn{\ge 0}).
#'   Candidate cut-probability exponent(s) \eqn{\gamma}. Default 1.
#' @param fuzziness Numeric \eqn{> 1}. Fuzzifier. Default 2.
#' @param metric Character. One of \code{"borda"} (default; rank aggregate),
#'   \code{"fuzzy_calinski_harabasz"}, \code{"Fukuyama_Sugeno"},
#'   \code{"dunn5"}, \code{"fuzzy_hyper_volume"}, \code{"stability"}.
#' @param crisp Logical. Produce hard memberships? Default \code{FALSE}.
#' @param strategy Character. \code{"grid"} (default), \code{"sequential"}, or
#'   \code{"iterative"}. See Details.
#' @param stability_B Integer. Bootstrap replicates for the stability
#'   sub-metric. Set to 0 to skip. Default 25.
#' @param prior_strength Non-negative numeric scalar \eqn{\lambda_E}. Passed to
#'   the membership computation. \code{Inf} short-circuits to an exact
#'   reproduction of \code{prior_typology}. Default 0. Not a tuning axis.
#' @param prior_typology Optional. Vertex attribute name or label vector; used
#'   only by the \code{prior_strength = Inf} short-circuit.
#' @param max_iter Integer. Maximum passes for \code{"iterative"}. Default 10.
#' @param max_ensembles Integer or \code{Inf}. Guard against an accidentally
#'   enormous search: the call fails before any work if the number of distinct
#'   ensemble triples exceeds this. Default 64.
#' @param large_n_threshold Integer. CLARA/PAM cutoff. Default 20000.
#' @param seed Optional integer seed; the caller's RNG state is restored on
#'   exit.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return An object of class \code{"regions_tuning"}: a list with
#'   \code{best_result} (a \code{"fuzzy_clusters"} object), \code{best_config}
#'   (one-row data frame of the winning values on all five axes), the individual
#'   slots \code{best_intermediate_regions}, \code{best_final_regions},
#'   \code{best_partitions_per_tree}, \code{best_weight_exponent},
#'   \code{best_n_st}, plus \code{best_score}, \code{tuning_log} (one row per
#'   configuration, with all five axis columns and all CVI columns),
#'   \code{candidates}, \code{n_ensembles_built}, \code{metric},
#'   \code{strategy}, and \code{prior_strength}.
#'
#' @export
tune_regions <- function(graph, 
                         trees,
                         n_st,
                         intermediate_regions = 50,
                         final_regions,
                         partitions_per_tree = 1L,
                         weight_exponent     = 1,
                         fuzziness  = 2,
                         metric = c(
                                 "borda",
                                 "fuzzy_calinski_harabasz",
                                 "Fukuyama_Sugeno",
                                 "dunn5",
                                 "fuzzy_hyper_volume",
                                 "stability"
                         ),
                         crisp = FALSE,
                         strategy = c("grid", "sequential", "iterative"),
                         stability_B       = 25L,
                         prior_strength    = 0,
                         prior_typology    = NULL,
                         max_iter          = 10,
                         size_influence    = 2,
                         max_ensembles     = 64L,
                         large_n_threshold = 20000,
                         seed              = NULL,
                         verbose           = TRUE) {
        
        # --- Checks -----------------------------------------------------*
        stopifnot(inherits(trees, "spanning_trees"))
        stopifnot(is.numeric(prior_strength), length(prior_strength) == 1L,
                  prior_strength >= 0)
        
        # --- Argument Matching  -----------------------------------------*
        strategy <- match.arg(strategy)
        metric   <- match.arg(metric)
        
        # --- Exact prior reproduction short-circuit --------------------------*
        if (is.infinite(prior_strength)) {
                if (is.null(prior_typology)) {
                        stop("`prior_strength = Inf` requires `prior_typology`.", call. = FALSE)
                }
                prior_labels <- .resolve_prior_typology(graph, prior_typology)
                return(.exact_prior_result(prior_labels, verbose = verbose))
        }
        
        # --- Normalize candidate sets ----------------------------------------*
        # weight_exponent stays double; the other four are counts.
        intermediate_regions <- sort(unique(as.integer(intermediate_regions)))
        final_regions        <- sort(unique(as.integer(final_regions)))
        partitions_per_tree  <- sort(unique(as.integer(partitions_per_tree)))
        n_st                 <- sort(unique(as.integer(n_st)))
        weight_exponent      <- sort(unique(as.numeric(weight_exponent)))
        
        # --- Checks -----------------------------------------------------*
        stopifnot(length(intermediate_regions) >= 1L, all(intermediate_regions >= 2),
                  all(intermediate_regions <= trees$n_vertices))
        stopifnot(length(final_regions) >= 1L, all(final_regions >= 2))
        stopifnot(length(partitions_per_tree) >= 1L, all(partitions_per_tree >= 1))
        stopifnot(length(n_st) >= 1L, all(n_st >= 2), all(n_st <= trees$n_trees))
        stopifnot(length(weight_exponent) >= 1L, all(is.finite(weight_exponent)),
                  all(weight_exponent >= 0))
        stopifnot(is.numeric(max_ensembles), length(max_ensembles) == 1L,
                  max_ensembles >= 1)
        
        cand <- list(
                intermediate_regions = intermediate_regions,
                weight_exponent      = weight_exponent,
                partitions_per_tree  = partitions_per_tree,
                n_st                 = n_st,
                final_regions        = final_regions
        )
        active_axes <- .TUNE_AXES[vapply(cand[.TUNE_AXES], length, integer(1)) > 1L]
        
        # Distinct ensemble triples: the actual cost of the search.
        n_ens_triples <- length(intermediate_regions) *
                length(weight_exponent) * length(partitions_per_tree)
        if (strategy == "grid" && n_ens_triples > max_ensembles) {
                stop(sprintf(paste0(
                        "Grid requires %d distinct ensemble computations (limit `max_ensembles` = %s). ",
                        "Use strategy = \"sequential\" / \"iterative\", shrink the ",
                        "intermediate_regions / weight_exponent / partitions_per_tree ",
                        "candidate sets, or raise `max_ensembles` deliberately."),
                        n_ens_triples, format(max_ensembles)), call. = FALSE)
        }
        
        if (verbose && "n_st" %in% active_axes) {
                message("Note: most validity indices drift monotonically with ensemble ",
                        "size; every ensemble will be built at n_st = ",
                        max(n_st), ".")
        }
        
        # --- RNG handling ----------------------------------------------------*
        if (!is.null(seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                set.seed(seed)
        }
        
        # --- Single-slot ensemble cache --------------------------------------*
        # Key is the ensemble triple only. n_st is served by slicing columns of a
        # matrix built once at max(n_st): with ppt cut-sets per tree the columns
        # are grouped tree-major, so the first n_st * ppt columns are exactly the
        # ensemble that n_st trees would have produced.
        n_st_build <- max(n_st)
        ens_cache  <- new.env(parent = emptyenv())
        ens_cache$key   <- NULL
        ens_cache$mat   <- NULL
        ens_cache$p0c   <- NULL
        ens_cache$built <- 0L

        get_memb <- function(sr, we, ppt, ns) {
                
                key <- paste(sr, sprintf("%.17g", we), ppt, sep = "|")
                
                if (!identical(key, ens_cache$key)) {
                        ens <- compute_ensemble_memberships(
                                trees                = trees,
                                n_st                 = n_st_build,
                                graph                = graph,
                                intermediate_regions = sr,
                                prior_strength       = prior_strength,
                                partitions_per_tree  = ppt,
                                weight_exponent      = we,
                                size_influence       = size_influence,
                                verbose              = verbose
                        )
                        # Drop the previous matrix before holding two at once.
                        ens_cache$mat   <- NULL
                        ens_cache$mat   <- ens$memberships
                        ens_cache$p0c   <- ens$chance_coassoc_by_col
                        ens_cache$key   <- key
                        ens_cache$built <- ens_cache$built + 1L
                }
                n_col <- ns * ppt
                mat <- if (n_col >= ncol(ens_cache$mat)) {
                        ens_cache$mat
                } else {
                        ens_cache$mat[, seq_len(n_col), drop = FALSE]
                }
                .new_ensemble_memberships(
                        memb                 = mat,
                        n_st                 = ns,
                        intermediate_regions = sr,
                        partitions_per_tree  = ppt,
                        weight_exponent      = we,
                        prior_strength       = prior_strength,
                        n_vertices           = nrow(mat),
                        p0_by_col            = ens_cache$p0c[seq_len(n_col)]
                )
        }
        
        # Canonical column sets (kept fixed so rbind() works across batches).
        axis_cols <- .TUNE_AXES
        diag_cols <- c("P0", "N_EFF")
        cvi_cols  <- c("CH_D", "SIL_F","STAB")
        col_order <- c(axis_cols, diag_cols, cvi_cols)
        
        # --- Score one configuration -> one row of CVI columns ---------------*
        score_config <- function(cfg) {
                
                memb_mat <- get_memb(cfg$intermediate_regions,
                                     cfg$weight_exponent,
                                     cfg$partitions_per_tree,
                                     cfg$n_st)
                
                res <- cluster_consensus(
                        memb_mat, 
                        graph = graph, 
                        k = cfg$final_regions,
                        fuzziness = fuzziness, 
                        crisp = crisp,
                        large_n_threshold = large_n_threshold,
                        verbose = FALSE
                )
                cvis <- .compute_fuzzy_cvis(
                        res                      = res,
                        memb_mat                 = memb_mat$memberships,
                        stability_B              = stability_B,
                        stability_subsample_frac = 0.8,
                        stability_seed           = seed,
                        large_n_threshold        = large_n_threshold,
                        verbose                  = verbose
                )

                list (
                        out = cbind(
                                cfg[, axis_cols, drop = FALSE],
                                data.frame(
                                        P0 = memb_mat$chance_coassoc,
                                        N_EFF = 1 / memb_mat$chance_coassoc
                                ),
                                cvis
                        ), 
                        typicality = res$typicality,
                        entropy = res$entropy
                )

        }
        # --- tuning_log accumulator ------------------------------------------*
        tuning_log <- NULL
        
        rebuild_score_column <- function() {
                # Borda ranks are computed globally over all rows evaluated so far, so
                # the score column must be rebuilt after every batch.
                if (is.null(tuning_log) || nrow(tuning_log) == 0L) return(invisible(NULL))
                tuning_log$score <<- .derive_score_column(tuning_log, metric)
                invisible(NULL)
        }
        
        evaluate_configs <- function(cfg_df, tick = NULL) {
                cfg_df <- cfg_df[, axis_cols, drop = FALSE]
                seen   <- if (is.null(tuning_log)) character(0) else .tune_key(tuning_log)
                keys   <- .tune_key(cfg_df)
                new_rows <- vector("list", nrow(cfg_df))
                for (i in seq_len(nrow(cfg_df))) {
                        if (keys[i] %in% seen) {          # already evaluated
                                if (!is.null(tick)) tick(i)
                                next
                        }
                        new_rows[[i]] <- score_config(cfg_df[i, , drop = FALSE])$out
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
                others <- setdiff(axis_cols, axis)
                keep   <- rep(TRUE, nrow(tuning_log))
                for (o in others) keep <- keep & (tuning_log[[o]] == cur[[o]])
                idx <- which(keep)
                if (length(idx) == 0L) return(cur[[axis]])
                sc <- tuning_log$score[idx]
                if (all(is.na(sc))) return(cur[[axis]])
                tuning_log[[axis]][idx][which.max(sc)]
        }
        
        sweep_axis <- function(axis, label) {
                cfgs <- axis_configs(axis)
                # Descending n_st means the max-size ensemble is built on the first
                # visit and every smaller candidate is a free column slice.
                if (axis == "n_st") cfgs <- cfgs[order(-cfgs$n_st), , drop = FALSE]
                if (verbose) message(sprintf("  Sweeping %s (%d candidates)",
                                             axis, nrow(cfgs)))
                tick <- if (verbose) .make_progress(nrow(cfgs), label) else NULL
                evaluate_configs(cfgs, tick)
                cur[[axis]] <<- best_along(axis)
                if (verbose) message(sprintf("  >> Best %s = %s", axis,
                                             format(cur[[axis]])))
                invisible(NULL)
        }
        
        # --- Run the chosen strategy -----------------------------------------*
        if (verbose) {
                message(sprintf("Tuning (%s strategy, metric = %s) ...", strategy, metric))
                message(sprintf("  Active axes: %s",
                                if (length(active_axes)) paste(active_axes, collapse = ", ")
                                else "none (single configuration)"))
        }
        
        if (length(active_axes) == 0L) {
                # Degenerate case: nothing to search, evaluate the one configuration.
                evaluate_configs(cur_df())
                
        } else if (strategy == "grid") {
                # expand.grid varies its first argument fastest; ordering it this way
                # keeps configurations that share an ensemble contiguous, and visits
                # n_st descending so the max-size build happens first in each block.
                pairs <- expand.grid(
                        final_regions        = final_regions,
                        n_st                 = rev(n_st),
                        intermediate_regions = intermediate_regions,
                        weight_exponent      = weight_exponent,
                        partitions_per_tree  = partitions_per_tree,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
                )
                if (verbose) message(sprintf(
                        "  Grid: %d configurations across %d ensemble computations",
                        nrow(pairs), n_ens_triples))
                tick <- if (verbose) .make_progress(nrow(pairs), "Grid search") else NULL
                evaluate_configs(pairs, tick)
                
        } else if (strategy == "sequential") {
                evaluate_configs(cur_df())             # seed the incumbent
                for (axis in active_axes) sweep_axis(axis, axis)
                
        } else if (strategy == "iterative") {
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
        tuning_log$score <- .derive_score_column(tuning_log, metric)
        
        if (all(is.na(tuning_log$score))) {
                stop("All hyper-parameter combinations produced an undefined score; ",
                     "check the metric and inputs.", call. = FALSE)
        }
        
        best_idx <- which.max(tuning_log$score)
        best_cfg <- tuning_log[best_idx, axis_cols, drop = FALSE]
        rownames(best_cfg) <- NULL
        best_sc  <- tuning_log$score[best_idx]
        
        if (verbose) {
                message(sprintf(paste0(
                        "\n>> Best: intermediate_regions = %d, weight_exponent = %g, ",
                        "partitions_per_tree = %d, n_st = %d, final_regions = %d ",
                        "(%s score = %.4f)"),
                        best_cfg$intermediate_regions, best_cfg$weight_exponent,
                        best_cfg$partitions_per_tree, best_cfg$n_st,
                        best_cfg$final_regions, metric, best_sc))
                message("  Computing final memberships ...")
        }
        
        memb_mat <- get_memb(best_cfg$intermediate_regions,
                             best_cfg$weight_exponent,
                             best_cfg$partitions_per_tree,
                             best_cfg$n_st)
        best_res <- cluster_consensus(
                memb_mat, graph = graph, k = best_cfg$final_regions,
                fuzziness = fuzziness, crisp = crisp,
                large_n_threshold = large_n_threshold, verbose = verbose
        )
        
        structure(
                list(
                        best_result               = best_res,
                        best_config               = best_cfg,
                        best_intermediate_regions = best_cfg$intermediate_regions,
                        best_final_regions        = best_cfg$final_regions,
                        best_partitions_per_tree  = best_cfg$partitions_per_tree,
                        best_weight_exponent      = best_cfg$weight_exponent,
                        best_n_st                 = best_cfg$n_st,
                        best_score                = best_sc,
                        tuning_log                = tuning_log,
                        candidates                = cand,
                        active_axes               = active_axes,
                        n_ensembles_built         = ens_cache$built,
                        metric                    = metric,
                        strategy                  = strategy,
                        prior_strength            = prior_strength,
                        ensemble                  = best_res$ensemble_meta$memb_mat
                ),
                class = "regions_tuning"
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
    "memberships", "typicality", "entropy", "hard_clusters", "k",
    "d_to_medoids", "medoid_idx", "fuzziness",
    "config", "tuned", "score", "metric", "strategy", "tuning_log",
    "candidates", "n_ensembles_built",
    "prior_strength", "coverage", "fit", "trees", "ensemble",
    "p0", "n_eff"
)

#' Fuzzy memberships from an N x k distance-to-medoid matrix.
#'
#' \eqn{u_{ij} \propto d_{ij}^{-2/(m-1)}}, rows normalised to 1. Factored out of
#' cluster_consensus() so that recompute_memberships() cannot drift from it:
#' one implementation, two callers.
#'
#' Two degenerate rows are handled explicitly. A node that IS a medoid has an
#' exact zero distance and gets a one-hot row. A node whose smallest distance
#' underflows the power transform to Inf would otherwise produce NaN after
#' normalisation, so it also collapses to one-hot on its nearest medoid.
#'
#' @param d Numeric N x k matrix of distances, non-negative.
#' @param fuzziness Numeric > 1.
#' @return Numeric N x k matrix, rows sum to 1.
#' @keywords internal
#' @noRd
.memberships_from_d <- function(d, fuzziness) {
    stopifnot(is.matrix(d), is.numeric(fuzziness),
              length(fuzziness) == 1L, fuzziness > 1)
    n <- nrow(d); k <- ncol(d)
    U <- matrix(0, n, k)
    
    dp   <- d ^ (-2 / (fuzziness - 1))
    bad  <- (rowSums(d == 0) > 0L) | !is.finite(rowSums(dp))
    good <- !bad
    
    if (any(good)) {
        g <- dp[good, , drop = FALSE]
        U[good, ] <- g / rowSums(g)
    }
    if (any(bad)) {
        idx <- which(bad)
        # max.col on the negated distances is the argmin.
        U[cbind(idx, max.col(-d[idx, , drop = FALSE], ties.method = "first"))] <- 1
    }
    U
}

#' Row-wise Shannon entropy, normalised by log(k).
#'
#' Rows are renormalised before the calculation, so this is defined for a
#' possibilistic matrix as well, but for such a matrix it measures only the row
#' SHAPE and discards the overall level. Entropy of a typicality matrix is
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
.as_config <- function(x) {
    if (is.list(x) && !is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
    stopifnot(is.data.frame(x), nrow(x) == 1L)
    miss <- setdiff(.TUNE_AXES, names(x))
    if (length(miss))
        stop("configuration is missing axis column(s): ",
             paste(miss, collapse = ", "), call. = FALSE)
    x <- x[, .TUNE_AXES, drop = FALSE]
    for (nm in setdiff(.TUNE_AXES, "weight_exponent")) x[[nm]] <- as.integer(x[[nm]])
    x$weight_exponent <- as.numeric(x$weight_exponent)
    rownames(x) <- NULL
    x
}

#' Construct a "pulse_regions" object
#'
#' The single exit point for get_regions(). Every branch (tuned, untuned, exact
#' prior) routes through here, so the returned shape cannot drift between them.
#' Derived quantities (typicality, entropy, k) are computed once, here, rather
#' than independently in each branch.
#'
#' @param fit A "fuzzy_clusters" object from cluster_consensus(), or NULL for
#'   the exact-prior degenerate case, in which case \code{memberships} and
#'   \code{hard_clusters} must be supplied directly.
#' @param config One-row data frame or named list giving the five tuning axes.
#' @param trees The "spanning_trees" object used, or NULL.
#' @param prior_strength Numeric lambda_E actually applied.
#' @param tuned Logical. TRUE if a search was run.
#' @param score,metric,strategy Tuning metadata. NA in the untuned branch: the
#'   hardened silhouette is deliberately NOT substituted, because it is not
#'   comparable to a tuned score and a populated-but-incomparable slot is worse
#'   than an empty one.
#' @param tuning_log,candidates,n_ensembles_built Search provenance, or NULL/NA.
#' @param coverage Length-N logical, exact-prior case only. NULL otherwise.
#' @param memberships,hard_clusters Supplied directly only when \code{fit} is
#'   NULL. Ignored otherwise.
#' @param fuzziness Fuzzifier. Taken from \code{fit} when available.
#' @param validate Logical. Run the shape assertions. TRUE except in tight
#'   internal loops.
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
                               metric            = NA_character_,
                               strategy          = NA_character_,
                               tuning_log        = NULL,
                               candidates        = NULL,
                               n_ensembles_built = NA_integer_,
                               coverage          = NULL,
                               memberships       = NULL,
                               hard_clusters     = NULL,
                               fuzziness         = NA_real_,
                               validate          = TRUE,
                               ensemble          = NULL) {
    
    degenerate <- is.null(fit)
    
    # --- Pull the primitives ---------------------------------------------*
    if (!degenerate) {
        need <- c("memberships", "hard_clusters", "d_to_medoids",
                  "medoid_idx", "fuzziness")
        miss <- setdiff(need, names(fit))
        if (length(miss))
            stop("cluster_consensus() did not return: ",
                 paste(miss, collapse = ", "),
                 ". The pulse_regions contract cannot be met.", call. = FALSE)
        
        U             <- fit$memberships
        hard_clusters <- fit$hard_clusters
        d_to_medoids  <- fit$d_to_medoids
        medoid_idx    <- fit$medoid_idx
        fuzziness     <- fit$fuzziness
    } else {
        if (is.null(memberships) || is.null(hard_clusters))
            stop("`fit = NULL` requires `memberships` and `hard_clusters`.",
                 call. = FALSE)
        U            <- memberships
        d_to_medoids <- NULL      # no medoids exist in the exact-prior case
        medoid_idx   <- NULL
    }
    .p0v <- if (degenerate) NA_real_ else {
            v <- tryCatch(fit$ensemble_meta$chance_coassoc,
                          error = function(e) NULL)
            if (is.numeric(v) && length(v) > 0L) as.numeric(v) else NA_real_
    }
    # --- Derived quantities, computed exactly once ------------------------*
    k <- ncol(U)
    
    # Typicality is possibilistic: rows do NOT sum to 1. It is a separate
    # quantity from `memberships`, not a rescaling of it, and downstream code
    # must branch on which it is holding.
    typicality <- if (degenerate) {
        U                                   # one-hot; typicality == membership
    } else {
        tryCatch(compute_typicality(fit),
                 error = function(e) {
                     warning("compute_typicality() failed: ",
                             conditionMessage(e), "; returning NULL.",
                             call. = FALSE)
                     NULL
                 })
    }
    
    # Entropy is of the MEMBERSHIPS, always. See .row_entropy() on why the
    # typicality version is a different quantity and is not stored here.
    entropy <- .row_entropy(U)
    
    config <- .as_config(config)
    
    # --- Shape assertions -------------------------------------------------*
    if (isTRUE(validate)) {
        N <- nrow(U)
        stopifnot(is.matrix(U), is.numeric(U), N >= 1L, k >= 1L)
        stopifnot(length(hard_clusters) == N)
        stopifnot(length(entropy) == N)
        
        # Rows sum to 1 up to floating point. Uncovered vertices in the
        # exact-prior case are all-zero rows and are exempt.
        rs <- rowSums(U)
        chk <- if (!is.null(coverage)) coverage else rep(TRUE, N)
        if (any(abs(rs[chk] - 1) > 5e-3))
            warning("membership rows do not sum to 1 (max deviation ",
                    signif(max(abs(rs[chk] - 1)), 3), ").", call. = FALSE)
        
        if (!is.null(d_to_medoids)) {
            # This is the check that protects recompute_memberships(): if the
            # distances and the memberships disagree in shape, every
            # recomputed matrix would be silently wrong.
            stopifnot(is.matrix(d_to_medoids),
                      nrow(d_to_medoids) == N,
                      ncol(d_to_medoids) == k)
            stopifnot(length(medoid_idx) == k)
            stopifnot(is.numeric(fuzziness), length(fuzziness) == 1L,
                      fuzziness > 1)
        }
        if (!is.null(typicality))
            stopifnot(is.matrix(typicality), nrow(typicality) == N,
                      ncol(typicality) == k)
        if (!is.null(coverage)) stopifnot(length(coverage) == N)
        if (!degenerate && config$final_regions != k)
            warning("config$final_regions (", config$final_regions,
                    ") disagrees with ncol(memberships) (", k, ").",
                    call. = FALSE)
    }
    
    p0        = mean(.p0v, na.rm = TRUE)
    m <- mean(.p0v, na.rm = TRUE)
    n_eff <- if (is.finite(m) && m > 0) 1 / m else NA_real_ 
    
    out <- list(
        memberships       = U,
        typicality        = typicality,
        entropy           = entropy,
        hard_clusters     = as.integer(hard_clusters),
        k                 = as.integer(k),
        d_to_medoids      = d_to_medoids,
        medoid_idx        = medoid_idx,
        fuzziness         = fuzziness,
        config            = config,
        tuned             = isTRUE(tuned),
        score             = as.numeric(score),
        metric            = as.character(metric),
        strategy          = as.character(strategy),
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


#' Recompute fuzzy memberships at a different fuzzifier
#'
#' The membership transform \eqn{u_{ij} \propto d_{ij}^{-2/(m-1)}} is a
#' deterministic function of the stored distances to medoids, so changing the
#' fuzzifier requires no re-clustering and no re-sampling of spanning trees.
#'
#' Useful for sensitivity analysis (how much does the reported uncertainty
#' depend on \emph{m}?) and for calibration work, where \emph{m} is the natural
#' recalibration parameter: it rescales confidence monotonically, so it changes
#' calibration completely while leaving every rank-based diagnostic untouched.
#'
#' @param x A "pulse_regions" object.
#' @param fuzziness Numeric > 1.
#' @param update Logical. If TRUE (default) return a modified copy of \code{x}
#'   with memberships, entropy and fuzziness updated; if FALSE return the bare
#'   N x k matrix.
#' @return A "pulse_regions" object or a numeric matrix.
#' @export
recompute_memberships <- function(x, fuzziness, update = TRUE) {
    stopifnot(inherits(x, "pulse_regions"))
    if (is.null(x$d_to_medoids))
        stop("this object carries no `d_to_medoids` (exact-prior result?); ",
             "memberships cannot be recomputed.", call. = FALSE)
    U <- .memberships_from_d(x$d_to_medoids, fuzziness)
    if (!isTRUE(update)) return(U)
    x$memberships <- U
    x$entropy     <- .row_entropy(U)
    x$fuzziness   <- fuzziness
    # hard_clusters come from PAM/CLARA, not from argmax(U), so they are
    # unchanged by construction. typicality is left as-is: it is derived from
    # the same distances but under the possibilistic transform, which has its
    # own scale parameter.
    x
}

#' @export
print.pulse_regions <- function(x, ...) {
    cat(sprintf("<pulse_regions>  %d units x %d regions | m = %s | %s\n",
                nrow(x$memberships), x$k, format(x$fuzziness),
                if (x$tuned) sprintf("tuned (%s, %s)", x$metric, x$strategy)
                else "single configuration"))
    cfg <- x$config
    cat(sprintf("  ir = %d | gamma = %g | ppt = %d | n_st = %d\n",
                cfg$intermediate_regions, cfg$weight_exponent,
                cfg$partitions_per_tree, cfg$n_st))
    cat(sprintf("  mean normalised entropy = %.3f | lambda_E = %s\n",
                mean(x$entropy), format(x$prior_strength)))
    if (!is.null(x$coverage))
        cat(sprintf("  prior coverage: %d/%d vertices\n",
                    sum(x$coverage), length(x$coverage)))
    invisible(x)
}