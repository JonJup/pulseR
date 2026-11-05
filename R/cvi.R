#' ## -- CVIs -- ## 



#' Pseudo-F (Calinski-Harabasz) and R^2 of the hard partition on `env`.
#' `env` should be the same transformed predictors whose Euclidean distance is
#' the edge dissimilarity (e.g. the GDM I-spline-transformed predictor matrix),
#' so that this is a PERMANOVA-style criterion in the space the trees were cut in.
#' @keywords internal
#' @noRd
.pseudo_F_env <- function(env, hard, k) {
        X <- as.matrix(env)
        stopifnot(nrow(X) == length(hard))
        ok <- stats::complete.cases(X)
        X <- X[ok, , drop = FALSE]; h <- hard[ok]
        N <- nrow(X)
        n_k <- tabulate(h, k)
        G   <- matrix(0, k, ncol(X))
        rs  <- rowsum(X, h)
        G[as.integer(rownames(rs)), ] <- rs / pmax(n_k[as.integer(rownames(rs))], 1)
        SSW <- sum((X - G[h, , drop = FALSE])^2)
        SST <- sum(scale(X, scale = FALSE)^2)
        SSB <- SST - SSW
        k_eff <- sum(n_k > 0)
        c(CH_ENV = (SSB / (k_eff - 1)) / (SSW / (N - k_eff)),
          R2_ENV = SSB / SST)
}

#' Proportion of ambiguous co-association (Senbabaoglu et al. 2014), estimated
#' on random unit pairs. Ensemble-level diagnostic: independent of k.
#' @keywords internal
#' @noRd
.pac <- function(memb, n_pairs = 20000L, lower = 0.1, upper = 0.9) {
        N <- nrow(memb)
        i <- sample.int(N, n_pairs, replace = TRUE)
        m <- sample.int(N, n_pairs, replace = TRUE)
        ok <- i != m
        ca <- rowMeans(memb[i[ok], , drop = FALSE] == memb[m[ok], , drop = FALSE])
        mean(ca > lower & ca < upper)
}

#' Unit-subsampling stability of the consensus (mean ARI vs. reference).
#' Refits with the SAME consensus procedure, so it measures how much the
#' regionalisation depends on which units are present.
#' @keywords internal
#' @noRd
.stability_consensus <- function(memb, hard, k, B = 25L, subsample_frac = 0.8,
                                 large_n_threshold = 20000, sample_size = 5000L,
                                 n_runs = 4L, n_cores = 1L, seed = NULL) {
        if (B <= 0L) return(NA_real_)
        N <- nrow(memb); n_sub <- floor(N * subsample_frac)
        if (n_sub < 2L * k) return(NA_real_)
        .local_seed(seed)
        aris <- vapply(seq_len(B), function(b) {
                idx   <- sort(sample.int(N, n_sub))
                lab_b <- tryCatch(
                        .consensus_partition(memb[idx, , drop = FALSE], k, large_n_threshold,
                                             sample_size, refine_iter = 3L, n_runs, n_cores,
                                             verbose = FALSE)$labels,
                        error = function(e) NULL)
                if (is.null(lab_b)) return(NA_real_)
                .adjusted_rand_index(lab_b, hard[idx])
        }, numeric(1))
        mean(aris, na.rm = TRUE)
}

#' One row of validity columns for a fitted consensus.
#'
#' Columns and their direction:
#'   CH_ENV, R2_ENV  external (environment)          higher better   NA if env = NULL
#'   STAB            unit-subsampling ARI            higher better
#'   SIL_CA          mean co-association silhouette  higher better   ensemble-internal
#'   CERT            mean own-region Pi membership   higher better   ensemble-internal
#'   EVI             expected VI(draw, consensus)    lower better    ensemble-internal
#'   PAC             ambiguous co-association        lower better    ensemble-level, k-free
#' @keywords internal
#' @noRd
.consensus_cvis <- function(res, memb, env = NULL,
                            stability_B = 25L, stability_subsample_frac = 0.8,
                            stability_seed = NULL, large_n_threshold = 20000,
                            sample_size = 5000L, n_cores = 1L, verbose = FALSE) {
        k <- res$k; hard <- res$hard_clusters
        fe <- if (is.null(env)) c(CH_ENV = NA_real_, R2_ENV = NA_real_) else .pseudo_F_env(env, hard, k)
        if (verbose && stability_B > 0L) message(sprintf("    stability (B = %d) ...", stability_B))
        stab <- .stability_consensus(memb, hard, k, B = stability_B,
                                     subsample_frac = stability_subsample_frac,
                                     large_n_threshold = large_n_threshold,
                                     sample_size = sample_size, n_cores = n_cores,
                                     seed = stability_seed)
        data.frame(CH_ENV = unname(fe["CH_ENV"]),
                   R2_ENV = unname(fe["R2_ENV"]),
                   STAB   = stab,
                   SIL_CA = res$avg_sil_width,
                   CERT   = mean(res$certainty, na.rm = TRUE),
                   EVI    = res$expected_vi,
                   PAC    = .pac(memb))
}

#' Direction of each validity column (for Borda).
#' @keywords internal
#' @noRd
.CVI_DIRECTION <- c(CH_ENV = "max", R2_ENV = "max", STAB = "max",
                    SIL_CA = "max", CERT = "max", EVI = "min", PAC = "min")

#' Columns entering the Borda aggregate. Only ensemble-EXTERNAL criteria by
#' default: the co-association indices are functions of the ensemble and
#' reward degenerate ensembles (gamma -> Inf) on the ensemble-shape axes.
#' @keywords internal
#' @noRd
.BORDA_COLS <- c("CH_ENV", "STAB")

# REPLACES the existing .derive_score_column()
.derive_score_column <- function(df, metric, voters = .BORDA_COLS) {
        switch(metric,
               borda = {
                       cols <- voters[vapply(voters, function(cc) !all(is.na(df[[cc]])), logical(1))]
                       if (!length(cols)) {
                               warning("No external criterion available (env = NULL, stability_B = 0); ",
                                       "Borda falls back to SIL_CA + EVI, which are only meaningful ",
                                       "when comparing k within one ensemble.", call. = FALSE)
                               cols <- c("SIL_CA", "EVI")
                       }
                       -.borda_score(df[, cols, drop = FALSE], unname(.CVI_DIRECTION[cols]))
               },
               ch_env        = df$CH_ENV,
               r2_env        = df$R2_ENV,
               stability     = df$STAB,
               silhouette_ca = df$SIL_CA,
               expected_vi   = -df$EVI,
               stop("Unknown metric: ", metric, call. = FALSE))
}










#' #' Adjusted Rand Index (Hubert & Arabie 1985)
#' #'
#' #' Pure-base implementation. Permutation-invariant; ARI = 1 for identical
#' #' partitions (up to label permutation), expected ARI = 0 under random
#' #' assignment.
#' #'
#' #' @param a,b Integer or character label vectors of equal length.
#' #' @return Numeric scalar in approximately \code{[-0.5, 1]}; \code{NA_real_}
#' #'   if input length < 2.
#' #'
#' #' @details Used by \code{.stability_fuzzy()} to score how well two hardened
#' #'   partitions agree on the overlap of resampled objects. We use the crisp ARI
#' #'   on argmax labels rather than a fuzzy-ARI generalisation because (i)
#' #'   stability of *structure* is what matters here, not stability of *softness*,
#' #'   and (ii) it has no dependencies. If a non-trivial fraction of objects have
#' #'   a top-two membership margin below ~0.1, consider substituting a
#' #'   Hüllermeier–Rifqi or Campello fuzzy Rand index (not implemented).
#' #'
#' #' @keywords internal
#' #' @noRd
#' .adjusted_rand_index <- function(a, b) {
#'         stopifnot(length(a) == length(b))
#'         if (length(a) < 2) return(NA_real_)
#'         tab <- table(a, b)
#'         n   <- sum(tab)
#'         if (n < 2) return(NA_real_)
#'         sum_nij <- sum(choose(tab, 2))
#'         ai <- rowSums(tab); bj <- colSums(tab)
#'         sum_ai <- sum(choose(ai, 2))
#'         sum_bj <- sum(choose(bj, 2))
#'         expected <- (sum_ai * sum_bj) / choose(n, 2)
#'         max_idx  <- 0.5 * (sum_ai + sum_bj)
#'         denom    <- max_idx - expected
#'         # Degenerate case: both partitions are singletons or a single block.
#'         if (denom == 0) return(1)
#'         (sum_nij - expected) / denom
#' }
#' 
#' 
#' #' Borda rank-aggregation across sub-metric columns
#' #'
#' #' For each column of \code{df}, ranks rows (1 = best per the direction); NAs
#' #' receive the worst rank within that column; fully-NA columns are dropped;
#' #' remaining ranks are summed row-wise. The optimal row MINIMISES the returned
#' #' Borda score.
#' #'
#' #' @param df Data frame of sub-metric values (one column per voter).
#' #' @param directions Character vector of length \code{ncol(df)}; each element
#' #'   \code{"min"} or \code{"max"}.
#' #' @return Numeric vector of length \code{nrow(df)}; all-NA if every voter is NA.
#' #'
#' #' @details Uses \code{ties.method = "average"} so genuine ties do not falsely
#' #'   separate candidates. Weighted Borda is intentionally not used: each voter
#' #'   is a different *aspect* of validity (compactness/separation, fuzziness-
#' #'   aware silhouette, stability), and over-weighting one would re-introduce the
#' #'   bias direction we are trying to balance.
#' #'
#' #' @keywords internal
#' #' @noRd
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
#' 
#' # -----------------------------------------------------------------------------*
#' # -----------------------------------------------------------------------------*
#' # -----------------------------------------------------------------------------*
#' 
#' #' Calinski-Harabasz index computed from a distance matrix
#' #'
#' #' Ratio of between-group to within-group dispersion, corrected for the number
#' #' of groups. Needs only a distance matrix and a membership matrix; no
#' #' coordinates and no centroids are required, so the same estimator works on
#' #' consensus Hamming distances (level 1) and Jensen-Shannon distances (level 3).
#' #'
#' #' @param D Distance matrix or \code{dist} object over N objects.
#' #' @param U N x k membership matrix.
#' #' @param fuzziness Fuzzifier m. Ignored when \code{crisp = TRUE}.
#' #' @param square Logical. If TRUE, \code{D} is squared before use. Set this so
#' #'   that the quantity entering the formula is the SQUARED distance in the
#' #'   sense of an ordinary geometric layout. See details in the notes below.
#' #' @param crisp Logical. Use one-hot argmax weights instead of u^m.
#' #' @return Numeric scalar; higher is better. NA_real_ if undefined.
#' #' @keywords internal
#' #' @noRd
#' .ch_distance <- function(D, U, fuzziness = 2, square = T, crisp = FALSE){
#'         
#'         if(inherits(D, "dist")) D <- as.matrix(D)
#'         N <- nrow(U) 
#'         k <- ncol(U)
#'         if (N != nrow(D)) stop("Dimensions of D and N do not agree", call. = FALSE)
#'         if (k < 2 || N < k) return(NA_real_)
#'         D2 <- if(square) D*D else D
#'         if (crisp){
#'                 Wt <- matrix(0, N, k)
#'                 Wt[cbind(seq_len(N), max.col(U, ties.method = "first"))] <- 1
#'         } else {
#'                 Wt <- U ^ fuzziness
#'         }
#'         
#'         # Group size-like notion. Excatlcy the group size for crisp. 
#'         mass <- colSums(Wt)
#'         if (any(!is.finite(mass)) || any(mass<=0)) return(NA_real_)
#'         
#'         # Within group distance. Wt provides fuzzy notion of within-groupies. 
#'         # The 0.5 corrects for counting every relationship twice (once from each side). 
#'         within_k <- 0.5 * colSums(Wt * (D2 %*% Wt)) / mass
#'         W <- sum(within_k)
#'         
#'         # Total spread 
#'         Tot <- 0.5 * sum(D2) / N
#'         # Between groups 
#'         s   <- rowSums(D2)
#'         h   <- rowSums(Wt)
#'         M   <- sum(mass)
#'         B   <- sum(h * s) / N - W - M * Tot / N
#'         if (!is.finite(B) || !is.finite(W) || B <= 0 || W <= 0) return(NA_real_)
#'         (B/(k - 1))/(W/(N-k))
#' }
#' 
#' 
#' 
#' #' Compute fuzzy-aware Cluster Validity Indices on a final fuzzy partition
#' #'
#' #' Returns a 1-row data frame of values, NOT a single scalar. The caller derives
#' #' a \code{score} column according to the chosen \code{metric} (see
#' #' \code{.derive_score_column}).
#' #'
#' #' @param res results of cluster_consensus()
#' #' @param env Environmental data passed through to \code{.compute_fuzzy_cvis()}.
#' #'  The same data used to weight the edges in \code{\link{add_edge_weight}}.
#' #' @param env_id Name of any ID columns in `env`, which should not be used in 
#' #'   the computation of cluster validity metrics. If `NULL` (default) no columns
#' #'   are removed from `env`.
#' #' @param d_to_medoids N x k matrix of distances from each object to each
#' #'   medoid (from \code{cluster_consensus}).
#' #' @param medoid_dist k x k matrix of inter-medoid distances.
#' #' @param memb_mat Integer matrix N x n_trees. Required if \code{SIL_F} or
#' #'   \code{STAB} are to be computed; may be \code{NULL} otherwise (those entries
#' #'   become \code{NA}).
#' #' @param fuzziness Fuzzifier \emph{m}; used in the XB* numerator weighting.
#' #' @param silhouette_hard Numeric scalar; hardened average silhouette width from
#' #'   CLARA/PAM (passed through as \code{SILH_HARD}).
#' #' @param stability_B Integer; bootstrap replicates for STAB. \code{0} skips it.
#' #' @param stability_subsample_frac Fraction of rows per stability replicate.
#' #' @param stability_seed Optional seed.
#' #' @param large_n_threshold CLARA/PAM cutoff for stability re-fits.
#' #' @param clara_samples,clara_sampsize CLARA controls for stability re-fits.
#' #' @param verbose Logical.
#' #' @return One-row \code{data.frame} with columns PC, PE, MPC, PEN, XB_star,
#' #'   SIL_F, STAB, SILH_HARD.
#' #' @keywords internal
#' #' @noRd
#' .compute_fuzzy_cvis <- function(res,
#'                                 # env,
#'                                 # env_id = NULL,
#'                                 memb_mat = NULL,
#'                                 
#'                                 stability_B              = 0L,
#'                                 stability_subsample_frac = 0.8,
#'                                 stability_seed           = NULL,
#'                                 silf_alpha               = 2,
#'                                 large_n_threshold        = 20000L,
#'                                 clara_samples            = 20L,
#'                                 clara_sampsize           = NULL,
#'                                 verbose                  = FALSE) {
#'         
#'         
#'         if (is.null(memb_mat)) stop("No membership matrix provided", call. = FALSE)
#'         fuzziness = res$fuzziness
#'         # U = Membership matrix 
#'         U <- res$memberships
#'         D <- parallelDist::parallelDist(memb_mat, method = "hamming")
#'         SIL_F <- .fuzzy_silhouette(D, U, alpha = silf_alpha)
#'         # Hamming already plays the role of the squared distance, so square = FALSE.
#'         CH_D <- .ch_distance(D, U, fuzziness = fuzziness, square = FALSE)
#'         # --- STAB (bootstrap ARI) ---- *
#'         STAB <- NA_real_
#'         if (stability_B > 0) {
#'                 hard_clusters <- res$hard_clusters
#'                 
#'                 STAB <- tryCatch(
#'                         .stability_fuzzy(memb_mat = memb_mat,
#'                                          hard_clusters = hard_clusters, 
#'                                          k = res$k, 
#'                                          fuzziness = res$fuzziness, 
#'                                          B = stability_B, 
#'                                          subsample_frac = stability_subsample_frac, 
#'                                          large_n_threshold = large_n_threshold,
#'                                          clara_samples = clara_samples, 
#'                                          clara_sampsize = clara_sampsize,
#'                                          seed = stability_seed, 
#'                                          verbose  = verbose),
#'                         error = function(e) {
#'                                 if (verbose) message("    [STAB] failed: ", conditionMessage(e))
#'                                 NA_real_
#'                         }
#'                 )
#'         }
#'         
#'         
#'         data.frame(SIL_F = SIL_F, 
#'                    CH_D = CH_D,
#'                    STAB = STAB)
#'         
#'         
#'         
#' }
#' 
#' #' Crisp Rousseeuw silhouette from a dissimilarity matrix
#' #'
#' #' Vectorised over the cluster axis; O(N k) after the N x N matrix is
#' #' materialised. Intended for a subsample of moderate size (a few thousand
#' #' rows), not the full graph.
#' #'
#' #' @param D Dissimilarity, either a \code{dist} object or a full N x N matrix.
#' #' @param labels Integer vector of cluster labels of length N.
#' #' @return Numeric vector of per-object silhouette widths.
#' #' @keywords internal
#' #' @noRd
#' .rousseeuw_silhouette <- function(D, labels) {
#'         if (inherits(D, "dist")) D <- as.matrix(D)
#'         N <- length(labels)
#'         clusters <- sort(unique(labels))
#'         K <- length(clusters)
#'         if (K < 2L) return(rep(0, N))
#'         cluster_sizes <- as.integer(table(factor(labels, levels = clusters)))
#'         
#'         # N x K matrix of summed distances from each point to each cluster.
#'         sums <- matrix(0, nrow = N, ncol = K)
#'         for (j in seq_len(K)) {
#'                 in_j <- labels == clusters[j]
#'                 sums[, j] <- rowSums(D[, in_j, drop = FALSE])
#'         }
#'         
#'         own_idx <- match(labels, clusters)
#'         s   <- numeric(N)
#'         eps <- .Machine$double.eps
#'         for (i in seq_len(N)) {
#'                 o <- own_idx[i]
#'                 if (cluster_sizes[o] <= 1L) { s[i] <- 0; next }
#'                 means <- sums[i, ] / cluster_sizes
#'                 # Own-cluster mean excludes self (self-distance is 0, so subtract nothing
#'                 # from the sum but divide by size - 1).
#'                 means[o] <- sums[i, o] / (cluster_sizes[o] - 1L)
#'                 a_i <- means[o]
#'                 b_i <- min(means[-o])
#'                 s[i] <- (b_i - a_i) / max(a_i, b_i, eps)
#'         }
#'         s
#' }
#'