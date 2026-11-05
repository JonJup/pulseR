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

#' Split-half ensemble agreement of the consensus (STAB)
#'
#' Splits the ensemble columns at random into two disjoint halves, runs the
#' same consensus procedure at \code{k} on each half, and records the ARI
#' between the two half-ensemble partitions; repeated \code{B} times. The
#' halves share no trees and no cuts, so STAB measures whether independent
#' sets of trees and cuts give the same regions. Each half has half the
#' columns, so STAB is a conservative estimate for the full ensemble.
#'
#' @param memb Integer label matrix, N x M.
#' @param k Target region count.
#' @param B Number of random splits. 0 returns \code{NA}.
#' @return Mean ARI over splits, or \code{NA} if \code{M < 4}.
#' @keywords internal
#' @noRd
.stability_consensus <- function(memb, k, B = 25L,
                                 large_n_threshold = 20000, sample_size = 5000L,
                                 n_runs = 4L, n_cores = 1L, seed = NULL) {
        if (B <= 0L) return(NA_real_)
        M <- ncol(memb)
        if (M < 4L || nrow(memb) < 2L * k) return(NA_real_)
        .local_seed(seed)
        half <- M %/% 2L
        fit_half <- function(cols) tryCatch(
                .consensus_partition(memb[, cols, drop = FALSE], k, large_n_threshold,
                                     sample_size, refine_iter = 3L, n_runs, n_cores,
                                     verbose = FALSE)$labels,
                error = function(e) NULL)
        aris <- vapply(seq_len(B), function(b) {
                perm <- sample.int(M)
                la <- fit_half(sort(perm[seq_len(half)]))
                lb <- fit_half(sort(perm[half + seq_len(half)]))
                if (is.null(la) || is.null(lb)) return(NA_real_)
                .adjusted_rand_index(la, lb)
        }, numeric(1))
        mean(aris, na.rm = TRUE)
}

#' One row of validity columns for a fitted consensus.
#'
#' Columns and their direction:
#'   CH_ENV, R2_ENV  external (environment)          higher better   NA if env = NULL
#'   STAB            split-half ensemble ARI         higher better
#'   SIL_CA          mean co-association silhouette  higher better   ensemble-internal
#'   CERT            mean own-region Pi membership   higher better   ensemble-internal
#'   EVI             expected VI(draw, consensus)    lower better    ensemble-internal
#'   PAC             ambiguous co-association        lower better    ensemble-level, k-free
#' @keywords internal
#' @noRd
.consensus_cvis <- function(res, memb, env = NULL,
                            stability_B = 25L,
                            stability_seed = NULL, large_n_threshold = 20000,
                            sample_size = 5000L, n_cores = 1L, verbose = FALSE) {
        k <- res$k; hard <- res$hard_clusters
        fe <- if (is.null(env)) c(CH_ENV = NA_real_, R2_ENV = NA_real_) else .pseudo_F_env(env, hard, k)
        if (verbose && stability_B > 0L) message(sprintf("    stability (B = %d) ...", stability_B))
        stab <- .stability_consensus(memb, k, B = stability_B,
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

#' Borda rank aggregation across voter columns
#'
#' Ranks each column of \code{df} (1 = best for its direction) and returns the
#' row-wise mean rank. NAs share the average of the unused worst ranks;
#' all-NA columns are dropped. Lower is better.
#'
#' @param df Data frame, one column per voter.
#' @param directions Character vector, \code{"min"} or \code{"max"} per column;
#'   if named, used to select and order the columns of \code{df}.
#' @return Numeric vector of length \code{nrow(df)}; all \code{NA} if no voter
#'   has data.
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
