## Dirichlet-multinomial mixture (DMM) for Level 3 mosaic types.
## Base R only (stats for kmeans). Model: Holmes, Harris & Quince 2012, PLoS ONE 7:e30126.
## Alpha update: Minka (2000) fixed point, which never decreases the likelihood.

#' Convert signatures back to weighted pair counts
#'
#' @param sigs Output of \code{compute_all_signatures()}, rows summing to 1.
#' @param totals Per-row totals; defaults to \code{attr(sigs, "totals")}.
#' @param n_scale Multiplier on the totals (sensitivity analysis: try 0.5 and 2).
#' @return Numeric matrix of non-negative, possibly fractional counts, NA rows dropped.
#' @keywords internal
.sig_to_counts <- function(sigs, totals = attr(sigs, "totals"), n_scale = 1) {
        if (is.null(totals))
                stop("No totals: add `attr(sigs, \"totals\") <- rs` in compute_all_signatures().",
                     call. = FALSE)
        ok <- stats::complete.cases(sigs)
        X <- unclass(sigs)[ok, , drop = FALSE] * (totals[ok] * n_scale)
        attr(X, "valid") <- which(ok)
        X
}

## Grouped sums over contiguous blocks: one cumsum, about 5x faster than rowsum().
## `ends` = cumsum(tabulate(group)); empty groups return 0.
.gsum <- function(v, ends) {
        cs <- c(0, cumsum(v))
        cs[ends + 1L] - cs[c(0L, ends[-length(ends)]) + 1L]
}

## Non-zero triplet form of X: zero cells contribute nothing to either step.
## which() returns column-major order, so column blocks are contiguous; `orow` gives row blocks.
.dmm_triplets <- function(X) {
        nz <- which(X > 0, arr.ind = TRUE)
        tr <- list(i = nz[, 1], j = nz[, 2], x = X[nz], n = rowSums(X),
                   N = nrow(X), C = ncol(X))
        tr$orow  <- order(tr$i, method = "radix")
        tr$rends <- cumsum(tabulate(tr$i, tr$N))
        tr$cends <- cumsum(tabulate(tr$j, tr$C))
        tr$lcoef <- lgamma(tr$n + 1) - .gsum(lgamma(tr$x + 1)[tr$orow], tr$rends)
        tr
}

## N x K matrix of log DM(x_i | alpha_k), multinomial coefficient included.
.dmm_logdens <- function(tr, alpha) {
        K <- nrow(alpha)
        out <- matrix(0, tr$N, K)
        for (k in seq_len(K)) {
                a <- alpha[k, ]
                A <- sum(a)
                v <- lgamma(tr$x + a[tr$j]) - lgamma(a[tr$j])
                out[, k] <- lgamma(A) - lgamma(tr$n + A) + .gsum(v[tr$orow], tr$rends) + tr$lcoef
        }
        out
}

## Row-wise log-sum-exp via max.col (no apply()).
.lse_rows <- function(L) {
        m <- L[cbind(seq_len(nrow(L)), max.col(L, ties.method = "first"))]
        m + log(rowSums(exp(L - m)))
}

## Precision gradient (and Hessian) of the weighted DM log-likelihood in A = sum(alpha),
## mean m fixed. Zero cells contribute nothing, so only the non-zero triplets are summed.
.dmm_grad_A <- function(tr, r, m, A, hess = FALSE, r_min = 1e-8) {
        keep <- which(r[tr$i] > r_min)
        rk <- r[tr$i[keep]]; xk <- tr$x[keep]; mj <- m[tr$j[keep]]; amj <- A * mj
        rows <- r > r_min; nr <- tr$n[rows]; rr <- r[rows]
        g <- sum(rr * (digamma(A) - digamma(nr + A))) +
                sum(rk * mj * (digamma(xk + amj) - digamma(amj)))
        if (!hess) return(g)
        h <- sum(rr * (trigamma(A) - trigamma(nr + A))) +
                sum(rk * mj^2 * (trigamma(xk + amj) - trigamma(amj)))
        c(g, h)
}

## M-step for one component, weights r. Each inner step does
##  (i) a safeguarded Newton search for the precision A in log A, mean fixed. Newton
##      is used only where the likelihood is concave in log A; elsewhere it takes a
##      unit step uphill. Far above the optimum the likelihood is convex in log A,
##      and plain Newton then runs to infinity.
## (ii) one Minka (2000) fixed-point step on the full vector, which never decreases
##      the likelihood.
## Step (i) replaces the slow part of the old update: the fixed point alone moves A
## only a little per iteration, so a large or infinite optimum took hundreds of EM
## iterations to approach. A is capped at A_max (the multinomial limit).
.dmm_update_alpha <- function(tr, r, a, inner = 2L, newton = 5L, alpha_min = 1e-4,
                              A_max = 1e6, r_min = 1e-8) {
        keep <- which(r[tr$i] > r_min)
        rk   <- r[tr$i[keep]]; xk <- tr$x[keep]; jk <- tr$j[keep]
        rows <- r > r_min; nr <- tr$n[rows]; rr <- r[rows]
        w <- numeric(length(tr$x))
        for (s in seq_len(inner)) {
                m <- a / sum(a)
                for (t in seq_len(newton)) {
                        A  <- sum(a)
                        gh <- .dmm_grad_A(tr, r, m, A, hess = TRUE, r_min = r_min)
                        d1 <- A * gh[1]; d2 <- d1 + A^2 * gh[2]       # derivatives in log A
                        step <- if (is.finite(d2) && d2 < 0) -d1 / d2 else sign(d1)
                        step <- max(min(step, 2), -2)
                        A_new <- min(A * exp(step), A_max)
                        a <- pmax(A_new * m, alpha_min)
                        if (abs(step) < 1e-4 || A_new >= A_max) break
                }
                A   <- sum(a)
                w[keep] <- rk * (digamma(xk + a[jk]) - digamma(a[jk]))
                num <- .gsum(w, tr$cends)
                den <- sum(rr * (digamma(nr + A) - digamma(A)))
                a   <- pmax(a * num / den, alpha_min)
                if (sum(a) > A_max) a <- pmax(a * (A_max / sum(a)), alpha_min)
        }
        a
}

## Multinomial mixture EM: the DMM in the limit A -> Inf. One matrix product per
## iteration. Used to initialise the DMM and as its exact solution when the data
## show no overdispersion.
.fit_mm <- function(X, lab, K, max_iter = 200L, tol = 1e-7) {
        N <- nrow(X); R <- matrix(0, N, K); R[cbind(seq_len(N), lab)] <- 1
        ll_old <- -Inf
        for (it in seq_len(max_iter)) {
                Th   <- crossprod(R, X) + 1e-10; Th <- Th / rowSums(Th)
                pi_k <- pmax(colMeans(R), .Machine$double.eps)
                L    <- X %*% t(log(Th)) + rep(log(pi_k), each = N)
                lse  <- .lse_rows(L); R <- exp(L - lse); ll <- sum(lse)
                if (abs(ll - ll_old) < tol * abs(ll)) break
                ll_old <- ll
        }
        list(R = R, theta = Th)
}

#' Fit a Dirichlet-multinomial mixture
#'
#' Each start: k-means on the Hellinger embedding, then a multinomial mixture
#' (cheap), then DMM EM from there. If the precision gradient is still positive
#' at \code{A_max} for every component under the multinomial fit, the data show
#' no overdispersion: the multinomial mixture is the DMM solution and EM is
#' skipped (\code{limit = "multinomial"}).
#'
#' @param X Count matrix (N x C), from \code{.sig_to_counts()}. Fractional counts are fine.
#' @param K Number of mosaic types.
#' @param n_start Random restarts; best log-likelihood kept.
#' @param max_iter,tol EM limits; \code{tol} is on the relative log-likelihood change.
#' @param A_init Starting precision sum(alpha) for the DMM EM.
#' @param alpha_min Floor on alpha. A pair absent from a type would otherwise get
#'   alpha -> 0, which forbids any catchment containing that pair from joining it.
#' @param A_max Cap on the precision. At the cap the DMM equals a multinomial
#'   mixture to numerical precision for any realistic count total.
#' @param seed Optional integer seed.
#' @return Object of class \code{pulse_dmm}: memberships \code{U} (N x K), hard
#'   \code{cluster}, \code{alpha} (K x C), mixing weights \code{pi}, \code{loglik},
#'   \code{BIC}, \code{ICL} (both lower is better), \code{limit}, \code{trace},
#'   \code{converged}.
#' @keywords internal
.fit_dmm <- function(X, K, n_start = 5L, max_iter = 500L, tol = 1e-5,
                     A_init = 50, alpha_min = 1e-4, A_max = 1e6, seed = NULL) {
        X <- as.matrix(X)
        if (any(X < 0, na.rm = TRUE) || anyNA(X)) stop("`X` must be non-negative with no NA.", call. = FALSE)
        if (any(rowSums(X) <= 0)) stop("Every row of `X` needs a positive total.", call. = FALSE)
        if (K < 1L || K > nrow(X)) stop("`K` must be between 1 and nrow(X).", call. = FALSE)
        .local_seed(seed)
        
        tr  <- .dmm_triplets(X)
        H   <- sqrt(X / tr$n)                         # Hellinger embedding for initialisation
        eps <- .Machine$double.eps
        mstep <- function(R, alpha) {
                for (k in seq_len(K))
                        alpha[k, ] <- .dmm_update_alpha(tr, R[, k], alpha[k, ],
                                                        alpha_min = alpha_min, A_max = A_max)
                alpha
        }
        best <- NULL
        
        for (st in seq_len(n_start)) {
                lab <- if (K == 1L) rep(1L, tr$N) else
                        stats::kmeans(H, centers = K, nstart = 1L, iter.max = 50L,
                                      algorithm = "MacQueen")$cluster
                mm   <- .fit_mm(X, lab, K)
                R    <- mm$R
                pi_k <- pmax(colMeans(R), eps)
                
                ## No overdispersion anywhere: the multinomial mixture is the solution.
                limit <- all(vapply(seq_len(K), function(k)
                        .dmm_grad_A(tr, R[, k], mm$theta[k, ], A_max) > 0, logical(1)))
                if (limit) {
                        alpha <- pmax(mm$theta * A_max, alpha_min)
                        L <- .dmm_logdens(tr, alpha) + rep(log(pi_k), each = tr$N)
                        lse <- .lse_rows(L)
                        fit <- list(R = exp(L - lse), alpha = alpha, pi = pi_k, loglik = sum(lse),
                                    trace = sum(lse), converged = TRUE, start = st,
                                    limit = "multinomial")
                } else {
                        alpha <- mstep(R, pmax(mm$theta * A_init, alpha_min))
                        trace <- numeric(0); converged <- FALSE
                        for (it in seq_len(max_iter)) {
                                L   <- .dmm_logdens(tr, alpha) + rep(log(pi_k), each = tr$N)
                                lse <- .lse_rows(L)
                                R   <- exp(L - lse)
                                trace[it] <- sum(lse)
                                if (it > 1L && abs(trace[it] - trace[it - 1L]) < tol * abs(trace[it])) {
                                        converged <- TRUE; break
                                }
                                pi_k  <- pmax(colMeans(R), eps)
                                alpha <- mstep(R, alpha)
                        }
                        fit <- list(R = R, alpha = alpha, pi = pi_k, loglik = trace[length(trace)],
                                    trace = trace, converged = converged, start = st,
                                    limit = "dirichlet-multinomial")
                }
                if (is.null(best) || fit$loglik > best$loglik) best <- fit
        }
        
        ## At the cap the concentrations are not free parameters.
        df  <- if (best$limit == "multinomial") K * (tr$C - 1L) + (K - 1L) else K * tr$C + (K - 1L)
        BIC <- -2 * best$loglik + df * log(tr$N)
        R   <- best$R
        EN  <- -sum(R[R > 0] * log(R[R > 0]))
        colnames(best$alpha) <- colnames(X)
        structure(list(U = R, cluster = max.col(R, ties.method = "first"),
                       alpha = best$alpha, pi = best$pi, K = K,
                       loglik = best$loglik, BIC = BIC, ICL = BIC + 2 * EN, df = df,
                       concentration = rowSums(best$alpha), limit = best$limit,
                       trace = best$trace, converged = best$converged,
                       valid = attr(X, "valid")),
                  class = "pulse_dmm")
}

#' Memberships and log predictive density for new rows (held-out blocks, unsampled catchments)
#' @param fit A \code{pulse_dmm} object.
#' @param X_new Count matrix with the same columns as the fitted \code{X}.
#' @return List with \code{U}, \code{cluster}, per-row \code{logpred} and its sum \code{loglik}.
#' @keywords internal
.score_dmm <- function(fit, X_new) {
        X_new <- as.matrix(X_new)
        if (ncol(X_new) != ncol(fit$alpha)) stop("Column mismatch with the fitted model.", call. = FALSE)
        tr  <- .dmm_triplets(X_new)
        L   <- .dmm_logdens(tr, fit$alpha) + rep(log(fit$pi), each = tr$N)
        lse <- .lse_rows(L)
        U   <- exp(L - lse)
        list(U = U, cluster = max.col(U, ties.method = "first"),
             logpred = lse, loglik = sum(lse))
}

#' @export
print.pulse_dmm <- function(x, ...) {
        cat(sprintf("Dirichlet-multinomial mixture: K = %d, N = %d, C = %d\n",
                    x$K, nrow(x$U), ncol(x$alpha)))
        cat(sprintf("logLik %.1f | BIC %.1f | ICL %.1f (lower is better) | %s after %d iterations\n",
                    x$loglik, x$BIC, x$ICL,
                    if (x$converged) "converged" else "NOT converged", length(x$trace)))
        if (identical(x$limit, "multinomial"))
                cat("Multinomial limit: no overdispersion detected; concentrations sit at the cap.\n")
        cat("Sizes:        ", tabulate(x$cluster, x$K), "\n")
        cat("Concentration:", signif(x$concentration, 3), "\n")
        invisible(x)
}