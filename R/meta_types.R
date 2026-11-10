# =============================================================================*
# pulseR — LEVEL 2b: meta-types
#
# Why this exists
# ---------------*
# A local type belongs to exactly one region by construction, so knowing a
# unit's local type tells you its region with certainty. Any signature
# built on local types therefore encodes region identity, and the mosaic level
# recovers the region map (ARI 0.79 in the test case) while its distance matrix
# saturates, because two neighbourhoods in different regions share no classes at
# all.
#
# Merging local types ACROSS regions by how much their environmental densities
# overlap breaks that determinism. A meta-type occurs in several regions, so two
# distant neighbourhoods can share classes, signatures overlap, and what is left
# for the mosaic level to find is configuration -- which types sit next to which
# -- rather than location.
#
# The test of whether it worked is not cluster quality. It is whether meta-types
# cross regions: see the `crossing` element and the print method.
#
# Requires: stats, cluster. Optional: ggplot2 (for plot_meta_types).
# =============================================================================*


# -----------------------------------------------------------------------------*
# 1. Distances between local types
# -----------------------------------------------------------------------------*

#' @keywords internal
#' @noRd
.local_gaussians <- function(local) {
        if (is.null(local$models))
                stop("This `pulse_local_types` object carries no `models`; refit with ",
                     "keep_models = TRUE.", call. = FALSE)
        mu <- list(); Sg <- list(); nm <- character(0); reg <- integer(0)
        for (r in seq_along(local$models)) {
                fit <- local$models[[r]]
                if (is.null(fit$parameters$mean))
                        stop(sprintf("Model for region %d has no parameters.", r), call. = FALSE)
                ## d x G. With one predictor mclust stores the means as a
                ## length-G vector, which as.matrix() would turn into G x 1.
                M <- matrix(fit$parameters$mean, nrow = length(local$variables))
                V <- fit$parameters$variance
                G <- ncol(M)
                for (g in seq_len(G)) {
                        mu[[length(mu) + 1L]] <- M[, g]
                        Sg[[length(Sg) + 1L]] <- .mclust_sigma(V, g, nrow(M))
                }
                nm  <- c(nm, local$key$name[local$key$region == r])
                reg <- c(reg, rep(r, G))
        }
        if (length(mu) != nrow(local$key))
                stop("Model components do not match the local-type key.", call. = FALSE)
        list(mu = mu, sigma = Sg, names = nm, region = reg)
}

#' @keywords internal
#' @noRd
.mclust_sigma <- function(V, g, d) {
        ## [[ ]] instead of $: `V$sigma` partially matches `sigmasq`, the slot
        ## mclust uses for one-predictor fits.
        if (!is.null(V[["sigma"]])) return(matrix(V[["sigma"]][, , g], d, d))
        if (!is.null(V[["Sigma"]])) return(matrix(V[["Sigma"]], d, d))  # shared across G
        if (!is.null(V[["cholsigma"]])) {
                R <- V[["cholsigma"]][, , g]; return(crossprod(R))
        }
        if (!is.null(V[["sigmasq"]])) {                    # univariate "E" / "V"
                s2 <- V[["sigmasq"]]
                return(matrix(if (length(s2) == 1L) s2 else s2[g], 1L, 1L))
        }
        stop("Unsupported mclust variance structure; refit with model_name = \"VVV\".",
             call. = FALSE)
}

#' @keywords internal
#' @noRd
.ridge <- function(S, eps = 1e-8) S + diag(eps * mean(diag(S)) + 1e-12, nrow(S))

#' Bhattacharyya distance between two Gaussians
#' @keywords internal
#' @noRd
.bhattacharyya <- function(m1, S1, m2, S2) {
        S1 <- .ridge(S1)
        S2 <- .ridge(S2)
        Sm <- (S1 + S2) / 2
        d  <- m1 - m2
        ld <- function(A) 2 * sum(log(diag(chol(A))))
        term1 <- drop(crossprod(d, solve(Sm, d))) / 8
        term2 <- (ld(Sm) - 0.5 * (ld(S1) + ld(S2))) / 2
        max(term1 + term2, 0)
}

#' 2-Wasserstein (Bures) distance between two Gaussians
#' @keywords internal
#' @noRd
.wasserstein2 <- function(m1, S1, m2, S2) {
        S1 <- .ridge(S1); S2 <- .ridge(S2)
        e  <- eigen(S2, symmetric = TRUE)
        R  <- e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% t(e$vectors)
        ev <- eigen(R %*% S1 %*% R, symmetric = TRUE)$values
        cross <- sum(sqrt(pmax(ev, 0)))
        sqrt(max(sum((m1 - m2)^2) + sum(diag(S1)) + sum(diag(S2)) - 2 * cross, 0))
}

#' Pairwise distances between local types
#'
#' @param local A \code{"pulse_local_types"} object.
#' @param method \code{"bhattacharyya"} (default) turns the overlap of the two
#'   fitted densities into a Hellinger distance, \eqn{\sqrt{1 - e^{-B}}}, bounded
#'   in \eqn{\\[0,1\\]} and directly readable as "how much do these two types
#'   describe the same conditions". \code{"wasserstein"} is unbounded and keeps
#'   discriminating between types that barely overlap, which matters if you
#'   later reuse this matrix as a ground metric. \code{"empirical"} weights the
#'   comparison by where units actually are in predictor space, so two
#'   types whose Gaussians overlap in a region of the space that is empty do not
#'   count as similar; it needs \code{env}.
#' @param env Predictor matrix for \code{method = "empirical"}, one row per unit,
#'   columns \code{local$variables}.
#' @param n_sample Rows sampled for the empirical comparison. Default 20000.
#' @param seed Seed for that sample, restored on exit. \code{NULL} (default)
#'   uses the current RNG stream.
#' @return A \code{dist} over the local types, with type names as labels.
#' @export
compute_local_type_distance <- function(local,
                                method = c("bhattacharyya", "wasserstein", "empirical"),
                                env = NULL, n_sample = 20000L, seed = NULL) {
        method <- match.arg(method)
        gs <- .local_gaussians(local)
        C  <- length(gs$mu)
        
        if (method == "empirical") {
                if (is.null(env)) stop("`method = \"empirical\"` needs `env`.", call. = FALSE)
                if (!is.numeric(n_sample) || length(n_sample) != 1L || is.na(n_sample) ||
                    n_sample < 2 || n_sample != round(n_sample))
                        stop("`n_sample` must be a single integer >= 2.", call. = FALSE)
                N <- nrow(env)
                if (N > n_sample) {
                        .local_seed(seed)
                        env <- env[sort(sample.int(N, n_sample)), , drop = FALSE]
                }
                X <- as.matrix(env[, local$variables, drop = FALSE])
                if (!is.null(local$scaling))
                        X <- sweep(sweep(X, 2L, local$scaling$center, "-"), 2L, local$scaling$scale, "/")
                ## Column c = the density of type c over the sampled units,
                ## normalised to a distribution. Bhattacharyya coefficient between
                ## columns then measures overlap where the data actually live.
                L <- vapply(seq_len(C), function(c)
                        .log_dmvn(X, gs$mu[[c]], gs$sigma[[c]]), numeric(nrow(X)))
                L <- exp(sweep(L, 2, apply(L, 2, max), "-"))
                L <- sweep(L, 2, colSums(L), "/")
                BC <- crossprod(sqrt(L))
                D  <- sqrt(pmax(1 - BC, 0))
                diag(D) <- 0
        } else {
                D <- matrix(0, C, C)
                for (i in seq_len(C - 1L)) for (j in (i + 1L):C) {
                        v <- if (method == "bhattacharyya") {
                                sqrt(max(1 - exp(-.bhattacharyya(gs$mu[[i]], gs$sigma[[i]],
                                                                 gs$mu[[j]], gs$sigma[[j]])), 0))
                        } else {
                                .wasserstein2(gs$mu[[i]], gs$sigma[[i]],
                                              gs$mu[[j]], gs$sigma[[j]])
                        }
                        D[i, j] <- D[j, i] <- v
                }
        }
        dimnames(D) <- list(gs$names, gs$names)
        attr(D, "region") <- gs$region
        stats::as.dist(D)
}

#' @keywords internal
#' @noRd
.log_dmvn <- function(X, mu, S) {
        S <- .ridge(S); R <- chol(S); d <- ncol(X)
        Z <- backsolve(R, t(sweep(X, 2, mu, "-")), transpose = TRUE)
        -0.5 * colSums(Z^2) - sum(log(diag(R))) - 0.5 * d * log(2 * pi)
}


# -----------------------------------------------------------------------------*
# 2. Meta-types
# -----------------------------------------------------------------------------*

#' Merge local types into cross-cutting meta-types
#'
#' Groups the local types of all regions by how much their fitted densities
#' overlap, producing a smaller set of classes that recur across regions.
#'
#' @param local A \code{"pulse_local_types"} object.
#' @param n_meta_types Integer vector. One value fixes the number of
#'   meta-types. Several values (default \code{2:12}) are candidates: the
#'   largest one whose confined mass is at most \code{max_confined} is chosen.
#'   Confined mass is the share of unit membership held by meta-types that
#'   occur in one region only. Coarse cuts always cross regions, so the rule
#'   takes the finest cut that still does. Warns if no candidate passes or the
#'   largest one does.
#' @param max_confined Threshold for candidate selection, in \eqn{[0, 1]}.
#'   Default 0.2. Some confinement is real, so 0 is usually too strict.
#' @param sensitivity Thresholds at which the selection is re-read from
#'   \code{selection}, without refitting. Default \code{c(0.1, 0.2, 0.3)}.
#' @param method,env,n_sample,seed Passed to \code{\link{compute_local_type_distance}}.
#' @param linkage Passed to \code{\link[stats]{hclust}}. Default
#'   \code{"average"}, which merges on mean overlap and does not chain the way
#'   single linkage would.
#' @param distance Optional precomputed \code{dist} over local types, in place of
#'   recomputing one.
#' @return An object of class \code{"pulse_meta_types"}:
#'   \describe{
#'     \item{\code{memberships}}{N x n_meta_types matrix. Rows sum to one, exactly, by
#'       summing local-type columns.}
#'     \item{\code{hard_labels}, \code{certainty}}{Row argmax and its value.}
#'     \item{\code{entropy}}{Row entropy of \code{memberships}, normalized by
#'       \eqn{\log} \code{n_meta_types}.}
#'     \item{\code{key}}{Local type, its region, and its meta-type.}
#'     \item{\code{crossing}}{Regions represented in each meta-type, and the ARI
#'       between the meta-type and region labelings OF THE LOCAL TYPES. Low ARI
#'       is the goal: it means meta-types are not region proxies.}
#'     \item{\code{n_meta_types}}{Number of meta-types.}
#'     \item{\code{selection}}{Per candidate: confined mass, share of
#'       single-region meta-types, regions per meta-type, ARI with region,
#'       silhouette, and whether it passes. \code{NULL} for a single value
#'       outside 2 to (number of local types - 1).}
#'     \item{\code{sensitivity}}{Per threshold: chosen n_meta_types, confined mass,
#'       unweighted and mass-weighted ARI, silhouette. Attribute \code{stable}
#'       is TRUE if every threshold picks the same cut. \code{NULL} if
#'       \code{n_meta_types} was fixed.}
#'     \item{\code{rule}}{How \code{n_meta_types} was chosen.}
#'     \item{\code{hclust}, \code{distance}}{The tree and the distance used.}
#'     \item{\code{method}}{Distance method used.}
#'     \item{\code{parent}}{The local-type key of \code{local}.}
#'   }
#' @importFrom stats hclust cutree as.dist
#' @importFrom cluster silhouette
#' @export
get_meta_types <- function(local,
                           n_meta_types     = 2:12,
                           max_confined = 0.2,
                           sensitivity  = c(0.1, 0.2, 0.3),
                           method     = c("bhattacharyya", "wasserstein", "empirical"),
                           env        = NULL,
                           n_sample   = 20000L,
                           seed       = NULL,
                           linkage    = c("average","ward.D", "single", "complete",
                                          "mcquitty", "median", "centroid", "ward.D2"),
                           distance   = NULL) {
        
        stopifnot(inherits(local, "pulse_local_types"))
        method <- match.arg(method)
        linkage <- match.arg(linkage)
        
        d <- if (is.null(distance)) compute_local_type_distance(local, method, env, n_sample, seed) else distance
        C <- attr(d, "Size")
        reg <- local$key$region
        hc  <- stats::hclust(d, method = linkage)
        
        
        if (!is.numeric(n_meta_types) || !length(n_meta_types) || anyNA(n_meta_types) ||
            any(n_meta_types != round(n_meta_types)))
                stop("`n_meta_types` must be one or more integers.", call. = FALSE)
        n_meta_types <- sort(unique(as.integer(n_meta_types)))
        auto   <- length(n_meta_types) > 1L
        cand   <- n_meta_types[n_meta_types >= 2L & n_meta_types < C]
        if (auto && !length(cand))
                stop(sprintf("No candidate in `n_meta_types` lies in 2..%d (the number of local types minus one).",
                             C - 1L), call. = FALSE)
        ## Unit mass of each local type: share of the N rows' membership.
        w <- colSums(local$memberships) / nrow(local$memberships)
        sel <- if (!length(cand)) NULL else do.call(rbind, lapply(cand, function(M) {
                cl <- stats::cutree(hc, k = M)
                nr <- vapply(split(reg, cl), function(v) length(unique(v)), integer(1))
                confined <- as.integer(names(nr)[nr == 1L])
                data.frame(n_meta_types = M,
                           mass_confined = sum(w[cl %in% confined]),
                           frac_single_region = mean(nr == 1L),
                           mean_regions_per_meta = mean(nr),
                           min_regions_per_meta = min(nr),
                           ari_with_region = .adjusted_rand_index(cl, reg),
                           ari_with_region_w = .weighted_ari(cl, reg, w * nrow(local$memberships)),
                           silhouette = mean(cluster::silhouette(cl, d)[, 3]))
        }))
        if (!is.null(sel)) sel$passes <- sel$mass_confined <= max_confined + 1e-9
        
        .pick <- function(t) {
                p <- sel$mass_confined <= t + 1e-9
                if (any(p)) max(sel$n_meta_types[p])
                else max(sel$n_meta_types[sel$mass_confined <= min(sel$mass_confined) + 1e-12])
        }
        sens <- if (auto && !is.null(sel)) {
                s <- data.frame(max_confined = sort(unique(c(sensitivity, max_confined))))
                s$n_meta_types <- vapply(s$max_confined, .pick, numeric(1))
                s <- merge(s, sel[, c("n_meta_types", "mass_confined", "ari_with_region",
                                      "ari_with_region_w", "silhouette")], by = "n_meta_types")
                s <- s[order(s$max_confined), ]
                attr(s, "stable") <- length(unique(s$n_meta_types)) == 1L
                s
        } else NULL
        
        if (auto) {
                if (length(unique(reg)) < 2L)
                        stop("Crossing is undefined with one region; give one value in `n_meta_types`.",
                             call. = FALSE)
                if (any(sel$passes)) {
                        ## Finest cut that still crosses regions.
                        M <- max(sel$n_meta_types[sel$passes])
                        if (M == max(cand))
                                warning(sprintf(paste0(
                                        "n_meta_types = %d is the largest candidate and passes; a finer ",
                                        "cut may also pass. Extend `n_meta_types`."), M), call. = FALSE)
                } else {
                        best <- min(sel$mass_confined)
                        M <- max(sel$n_meta_types[sel$mass_confined <= best + 1e-12])
                        warning(sprintf(paste0(
                                "No candidate keeps confined mass <= %.2f (best %.2f at n_meta_types = %d). ",
                                "Meta-types will track regions; check `selection`."),
                                max_confined, best, M), call. = FALSE)
                }
                rule <- sprintf("crossing: largest n_meta_types with mass_confined <= %.2f",
                                max_confined)
        } else {
                M <- n_meta_types
                if (M < 1L || M > C)
                        stop("`n_meta_types` must lie in 1..", C, " (the number of local types).", call. = FALSE)
                rule <- "user"
        }
        cl <- stats::cutree(hc, k = M)
        
        
        ## Meta memberships: sum the local-type columns of each meta-type. Exact,
        ## because local-type rows already sum to one.
        U <- local$memberships
        idx <- split(seq_len(ncol(U)), cl)
        Um  <- vapply(idx, function(j) rowSums(U[, j, drop = FALSE]), numeric(nrow(U)))
        colnames(Um) <- paste0("meta", names(idx))
        
        key <- data.frame(local_type = local$key$name,
                          region     = reg,
                          meta       = as.integer(cl),
                          meta_name  = paste0("meta", cl),
                          stringsAsFactors = FALSE)
        nr <- vapply(split(reg, cl), function(v) length(unique(v)), integer(1))
        
        hard <- max.col(Um, ties.method = "first")
        out <- list(
                memberships  = Um,
                hard_labels = as.integer(hard),
                certainty    = Um[cbind(seq_len(nrow(Um)), hard)],
                entropy      = .row_entropy(Um),
                key          = key,
                n_meta_types       = M,
                crossing     = list(regions_per_meta = nr,
                                    ari_with_region  = .adjusted_rand_index(cl, reg),
                                    ari_with_region_w = .weighted_ari(cl, reg, w * nrow(U)),
                                    frac_single_region = mean(nr == 1)),
                selection    = sel,
                sensitivity  = sens,
                rule         = rule,
                hclust       = hc,
                distance     = d,
                method       = method,
                parent       = local$key)
        structure(out, class = "pulse_meta_types")
}

#' @rdname get_meta_types
#' @param x An object of class \code{"pulse_meta_types"}.
#' @param ... Unused.
#' @return \code{print.pulse_meta_types} returns \code{x} invisibly.
#' @export
print.pulse_meta_types <- function(x, ...) {
        cat(sprintf("<pulse_meta_types>  %d local types -> %d meta-types (%s, %s linkage)\n",
                    nrow(x$key), x$n_meta_types, x$method, x$hclust$method))
        cat(sprintf("  regions per meta-type: median %.0f, min %d | %.0f%% confined to one region\n",
                    stats::median(x$crossing$regions_per_meta),
                    min(x$crossing$regions_per_meta),
                    100 * x$crossing$frac_single_region))
        cat(sprintf("  ARI(meta-type, region) = %.3f  -- lower is better here\n",
                    x$crossing$ari_with_region))
        if (x$crossing$ari_with_region > 0.5)
                cat("  WARNING: meta-types still track regions closely. Merging has not\n",
                    "  decoupled class from location; the mosaic level will keep\n",
                    "  reproducing the region map.\n")
        cat(sprintf("  n_meta_types chosen by %s\n", x$rule))
        cat(sprintf("  mean certainty = %.3f\n", mean(x$certainty)))
        if (!is.null(x$sensitivity) && !isTRUE(attr(x$sensitivity, "stable")))
                cat(sprintf("  Choice depends on max_confined: n_meta_types = %s across %s\n",
                            paste(x$sensitivity$n_meta_types, collapse = "/"),
                            paste(x$sensitivity$max_confined, collapse = "/")))
        invisible(x)
}

# -----------------------------------------------------------------------------*
# 4. Did it work?
# -----------------------------------------------------------------------------*

#' Cross-tabulate meta-types against regions
#'
#' @param meta A \code{"pulse_meta_types"} object.
#' @return A table of local-type counts, meta-types by region.
#' @export
tabulate_meta_regions <- function(meta) table(meta = meta$key$meta, region = meta$key$region)