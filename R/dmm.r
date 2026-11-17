## Dirichlet-multinomial mixture (DMM) for Level 3 mosaic types.
## Base R only (stats for kmeans). Model: Holmes, Harris & Quince 2012, PLoS ONE 7:e30126.
## Alpha update: Minka (2000) fixed point, which never decreases the likelihood.

#' Convert signatures back to weighted pair counts
#'
#' @param signatures Output of \code{compute_all_signatures()}, rows summing to 1.
#' @param totals Per-row totals; defaults to \code{attr(signatures, "totals")}.
#' @param n_scale Multiplier on the totals (sensitivity analysis: try 0.5 and 2).
#' @return Numeric matrix of non-negative, possibly fractional counts, NA rows dropped.
#' @keywords internal
#' @noRd
.sig_to_counts <- function(signatures, totals = attr(signatures, "totals"), n_scale = 1) {
        if (is.null(totals))
                stop("No totals: add `attr(signatures, \"totals\") <- rs` in compute_all_signatures().",
                     call. = FALSE)
        ok <- stats::complete.cases(signatures)
        X <- unclass(signatures)[ok, , drop = FALSE] * (totals[ok] * n_scale)
        attr(X, "valid") <- which(ok)
        X
}
## n_scale = "auto": 1 / mean neighbourhoods per edge, so the summed
## pseudo-counts equal the number of distinct contributing edges.
.resolve_n_scale <- function(n_scale, signatures) {
        if (identical(n_scale, "auto")) {
                m <- attr(signatures, "edge_multiplicity")
                if (is.null(m) || !is.finite(m) || m <= 0)
                        stop("n_scale = \"auto\" needs attr(signatures, \"edge_multiplicity\"). ",
                             "Recompute with compute_all_signatures(), or pass a number.", call. = FALSE)
                return(1 / m)
        }
        if (!is.numeric(n_scale) || length(n_scale) != 1L || !is.finite(n_scale) || n_scale <= 0)
                stop("`dmm_control$n_scale` must be \"auto\" or one positive number.", call. = FALSE)
        n_scale
}
## Temperature tau on the per-unit evidence, r_ik ∝ pi_k * DM(x_i | alpha_k)^tau.
## "auto" and "calibrate" both start at 1 / mean neighbourhoods per edge; "calibrate"
## later replaces it by .calibrate_dmm_tau() for the selected k.
.resolve_temperature <- function(temperature, signatures) {
        if (identical(temperature, "auto") || identical(temperature, "calibrate")) {
                m <- attr(signatures, "edge_multiplicity")
                if (is.null(m) || !is.finite(m) || m <= 0)
                        stop("temperature = \"", temperature, "\" needs attr(signatures, ",
                             "\"edge_multiplicity\"). Recompute with compute_all_signatures(), ",
                             "or pass a number.", call. = FALSE)
                return(min(1, 1 / m))
        }
        if (!is.numeric(temperature) || length(temperature) != 1L || !is.finite(temperature) ||
            temperature <= 0 || temperature > 1)
                stop("`dmm_control$temperature` must be \"auto\", \"calibrate\" or one number in (0, 1].",
                     call. = FALSE)
        temperature
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

## DMM EM from given parameters. Always untempered: tempering the E-step is
## deterministic annealing, and below a critical tau all components merge into one
## (planted 4-zone fixture: ARI 0.93 at tau = 0.5, 0 at tau = 0.13).
.dmm_em <- function(tr, alpha, pi_k, max_iter = 500L, tol = 1e-5,
                    alpha_min = 1e-4, A_max = 1e6) {
        eps <- .Machine$double.eps
        trace <- numeric(0); converged <- FALSE
        for (it in seq_len(max_iter)) {
                L   <- .dmm_logdens(tr, alpha) + rep(log(pi_k), each = tr$N)
                lse <- .lse_rows(L)
                R   <- exp(L - lse)
                trace[it] <- sum(lse)
                if (it > 1L && abs(trace[it] - trace[it - 1L]) < tol * abs(trace[it])) {
                        converged <- TRUE; break
                }
                pi_k <- pmax(colMeans(R), eps)
                for (k in seq_len(nrow(alpha)))
                        alpha[k, ] <- .dmm_update_alpha(tr, R[, k], alpha[k, ],
                                                        alpha_min = alpha_min, A_max = A_max)
        }
        list(R = R, alpha = alpha, pi = pi_k, loglik = trace[length(trace)],
             trace = trace, converged = converged)
}

## Assemble a pulse_dmm object from untempered EM output (alpha, pi, loglik, trace,
## converged, limit). Temperature tau acts after fitting, as a power on each unit's
## evidence: U_ik ∝ pi_k DM(x_i | alpha_k)^tau. The criteria get the composite-
## likelihood magnitude adjustment (Pauli, Racugno & Ventura 2011): BIC uses
## tau * loglik, and ICL adds the entropy of the tempered U. loglik stays untempered.
.new_pulse_dmm <- function(fit, tr, K, tau, colnames_X = NULL, valid = NULL) {
        df  <- if (fit$limit == "multinomial") K * (tr$C - 1L) + (K - 1L) else K * tr$C + (K - 1L)
        BIC <- -2 * tau * fit$loglik + df * log(tr$N)
        L   <- tau * .dmm_logdens(tr, fit$alpha) + rep(log(fit$pi), each = tr$N)
        R   <- exp(L - .lse_rows(L))
        EN  <- -sum(R[R > 0] * log(R[R > 0]))
        colnames(fit$alpha) <- colnames_X
        structure(list(U = R, cluster = max.col(R, ties.method = "first"),
                       alpha = fit$alpha, pi = fit$pi, K = K, tau = tau,
                       loglik = fit$loglik, BIC = BIC, ICL = BIC + 2 * EN, df = df,
                       concentration = rowSums(fit$alpha), limit = fit$limit,
                       trace = fit$trace, converged = fit$converged, valid = valid),
                  class = "pulse_dmm")
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
#'   alpha -> 0, which forbids any unit containing that pair from joining it.
#' @param A_max Cap on the precision. At the cap the DMM equals a multinomial
#'   mixture to numerical precision for any realistic count total.
#' @param tau Temperature on the per-unit evidence, in [0, 1], applied after
#'   fitting: alpha and pi are the untempered estimates, memberships are
#'   U_ik ∝ pi_k DM(x_i | alpha_k)^tau (tau -> 0 sends them to \code{pi}), and
#'   BIC uses tau * loglik. 1 is the plain DMM.
#' @param seed Optional integer seed.
#' @return Object of class \code{pulse_dmm}: memberships \code{U} (N x K), hard
#'   \code{cluster}, \code{alpha} (K x C), mixing weights \code{pi}, \code{loglik},
#'   \code{BIC}, \code{ICL} (both lower is better; magnitude-adjusted by tau),
#'   \code{tau}, \code{limit}, \code{trace}, \code{converged}.
#' @keywords internal
#' @noRd
.fit_dmm <- function(X, K, n_start = 5L, max_iter = 500L, tol = 1e-5,
                     A_init = 50, alpha_min = 1e-4, A_max = 1e6, tau = 1, seed = NULL) {
        X <- as.matrix(X)
        if (any(X < 0, na.rm = TRUE) || anyNA(X)) stop("`X` must be non-negative with no NA.", call. = FALSE)
        if (any(rowSums(X) <= 0)) stop("Every row of `X` needs a positive total.", call. = FALSE)
        if (K < 1L || K > nrow(X)) stop("`K` must be between 1 and nrow(X).", call. = FALSE)
        if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0 || tau > 1)
                stop("`tau` must be one number in [0, 1].", call. = FALSE)
        n_u <- nrow(unique(sqrt(X / rowSums(X))))
        if (K > n_u) stop("`K` (", K, ") exceeds the ", n_u,
                          " distinct row compositions of `X`.", call. = FALSE)
        .local_seed(seed)
        
        tr  <- .dmm_triplets(X)
        H   <- sqrt(X / tr$n)                         # Hellinger embedding for initialisation
        eps <- .Machine$double.eps
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
                                    trace = sum(lse), converged = TRUE)
                } else {
                        alpha <- pmax(mm$theta * A_init, alpha_min)
                        for (k in seq_len(K))
                                alpha[k, ] <- .dmm_update_alpha(tr, R[, k], alpha[k, ],
                                                                alpha_min = alpha_min, A_max = A_max)
                        fit <- .dmm_em(tr, alpha, pi_k, max_iter = max_iter, tol = tol,
                                       alpha_min = alpha_min, A_max = A_max)
                }
                fit$start <- st
                fit$limit <- if (limit) "multinomial" else "dirichlet-multinomial"
                if (is.null(best) || fit$loglik > best$loglik) best <- fit
        }
        .new_pulse_dmm(best, tr, K, tau, colnames(X), attr(X, "valid"))
}

## Same fit, new temperature: recomputes U, cluster, BIC and ICL; alpha and pi unchanged.
.temper_dmm <- function(fit, X, tau) {
        tr <- .dmm_triplets(as.matrix(X))
        .new_pulse_dmm(list(alpha = unname(fit$alpha), pi = fit$pi, loglik = fit$loglik,
                            trace = fit$trace, converged = fit$converged, limit = fit$limit),
                       tr, fit$K, tau, colnames(fit$alpha), fit$valid)
}

## Simulate counts from a fitted DMM with given row totals `n`. Fractional totals:
## draw round(n) (at least 1) multinomial counts, then rescale each row to n.
.dmm_simulate <- function(alpha, pi_k, n) {
        N <- length(n); C <- ncol(alpha)
        z <- sample.int(nrow(alpha), N, replace = TRUE, prob = pi_k)
        G <- matrix(stats::rgamma(N * C, shape = alpha[z, , drop = FALSE]), N, C)
        empty <- rowSums(G) <= 0                       # all shapes tiny: fall back to the mode
        if (any(empty)) G[cbind(which(empty), max.col(alpha[z[empty], , drop = FALSE]))] <- 1
        size <- pmax(1L, round(n))
        X <- t(vapply(seq_len(N), function(i) as.numeric(stats::rmultinom(1L, size[i], G[i, ])),
                      numeric(C)))
        list(X = X * (n / size), z = z)
}

#' Temperature by parametric-bootstrap ECE
#'
#' Simulates \code{B} count matrices from \code{fit} with the observed row
#' totals, refits each by (untempered) EM started at the generating
#' parameters, so labels stay aligned with the simulated truth without
#' matching. Each refit is rescored at every grid temperature and the
#' top-class reliability bins are pooled over the \code{B} maps.
#'
#' The simulator is the fitted DMM itself: rows are independent and the model
#' is correct by construction. The search corrects sharpness from estimating
#' alpha and pi. It cannot see misspecification or dependence between rows.
#' @param fit A \code{pulse_dmm} object.
#' @param n Row totals of the fitted count matrix.
#' @param tau_grid Candidate temperatures in (0, 1].
#' @param B Simulated maps.
#' @return List with the selected \code{tau} and a \code{table} (tau, ece,
#'   mean confidence, accuracy) over the grid.
#' @keywords internal
#' @noRd
.calibrate_dmm_tau <- function(fit, n, tau_grid = exp(seq(log(0.01), 0, length.out = 25L)),
                               B = 50L, n_bins = 10L, max_iter = 200L, tol = 1e-5,
                               alpha_min = 1e-4, A_max = 1e6, seed = NULL) {
        .local_seed(seed)
        G <- length(tau_grid)
        s_conf <- s_acc <- n_bin <- matrix(0, G, n_bins)
        alpha <- unname(fit$alpha)
        for (b in seq_len(B)) {
                sim <- .dmm_simulate(alpha, fit$pi, n)
                tr  <- .dmm_triplets(sim$X)
                em  <- .dmm_em(tr, alpha, fit$pi, max_iter = max_iter, tol = tol,
                               alpha_min = alpha_min, A_max = A_max)
                ld  <- .dmm_logdens(tr, em$alpha)
                lp  <- rep(log(em$pi), each = tr$N)
                for (g in seq_len(G)) {
                        L    <- tau_grid[g] * ld + lp
                        R    <- exp(L - .lse_rows(L))
                        top  <- max.col(R, ties.method = "first")
                        conf <- R[cbind(seq_len(tr$N), top)]
                        bin  <- pmin(n_bins, 1L + floor(conf * n_bins))
                        s_conf[g, ] <- s_conf[g, ] + .gsum(conf[order(bin)], cumsum(tabulate(bin, n_bins)))
                        s_acc[g, ]  <- s_acc[g, ] + tabulate(bin[top == sim$z], n_bins)
                        n_bin[g, ]  <- n_bin[g, ] + tabulate(bin, n_bins)
                }
        }
        tot <- rowSums(n_bin)
        tab <- data.frame(tau = tau_grid, ece = rowSums(abs(s_acc - s_conf)) / tot,
                          mean_conf = rowSums(s_conf) / tot, accuracy = rowSums(s_acc) / tot)
        list(tau = tab$tau[which.min(tab$ece)], table = tab, B = B)
}

#' Memberships and log predictive density for new rows (held-out blocks, unsampled units)
#' @param fit A \code{pulse_dmm} object.
#' @param X_new Count matrix with the same columns as the fitted \code{X}.
#' @param tau Temperature; defaults to the one the model was fitted at.
#' @return List with \code{U}, \code{cluster}, per-row \code{logpred} and its sum \code{loglik}.
#' @keywords internal
#' @noRd
.score_dmm <- function(fit, X_new, tau = if (is.null(fit$tau)) 1 else fit$tau) {
        X_new <- as.matrix(X_new)
        if (ncol(X_new) != ncol(fit$alpha)) stop("Column mismatch with the fitted model.", call. = FALSE)
        tr  <- .dmm_triplets(X_new)
        ld  <- .dmm_logdens(tr, fit$alpha)
        lp  <- rep(log(fit$pi), each = tr$N)
        lse <- .lse_rows(ld + lp)                     # predictive density: always untempered
        L   <- tau * ld + lp
        U   <- exp(L - .lse_rows(L))
        list(U = U, cluster = max.col(U, ties.method = "first"),
             logpred = lse, loglik = sum(lse))
}

#' @rdname get_mosaic_types
#' @return \code{print.pulse_dmm} returns \code{x} invisibly.
#' @export
print.pulse_dmm <- function(x, ...) {
        cat(sprintf("Dirichlet-multinomial mixture: K = %d, N = %d, C = %d\n",
                    x$K, nrow(x$U), ncol(x$alpha)))
        cat(sprintf("logLik %.1f | BIC %.1f | ICL %.1f (lower is better) | %s after %d iterations\n",
                    x$loglik, x$BIC, x$ICL,
                    if (x$converged) "converged" else "NOT converged", length(x$trace)))
        if (identical(x$limit, "multinomial"))
                cat("Multinomial limit: no overdispersion detected; concentrations sit at the cap.\n")
        if (!is.null(x$tau) && x$tau != 1)
                cat(sprintf("Temperature:   %.3g (memberships, BIC and ICL; logLik untempered)\n", x$tau))
        cat("Sizes:        ", tabulate(x$cluster, x$K), "\n")
        cat("Concentration:", signif(x$concentration, 3), "\n")
        invisible(x)
}


## ---------------------------------------------------------------------------*
## Two-channel mixture: DMM on signature counts (composition) x diagonal
## Gaussian on configuration features (arrangement). One label per unit, shared
## by both channels. Tempered composite evidence (memberships, BIC, ICL):
##   tau * log DM(x_i | alpha_k) + w * log N(y_i | mu_k, diag(s2_k)) + log pi_k
## tau is the DMM temperature (edge overcounting); w weights the arrangement
## channel against the tempered composition evidence. As in .dmm_em(), EM runs
## untempered: it maximises log DM + (w / tau) log N, the same channel balance at
## temperature 1. Tempering inside EM is deterministic annealing and merges all
## components once the configuration channel carries no zone signal.
## ---------------------------------------------------------------------------*

## N x K matrix of log N(y_i | mu_k, diag(s2_k)).
.gauss_logdens <- function(Y, mu, s2) {
        K <- nrow(mu); out <- matrix(0, nrow(Y), K)
        for (k in seq_len(K)) {
                Z <- sweep(Y, 2, mu[k, ])^2
                out[, k] <- -0.5 * (rowSums(sweep(Z, 2, s2[k, ], "/")) +
                                    sum(log(2 * pi * s2[k, ])))
        }
        out
}

## Weighted means and variances per component; s2_min floors the variances.
.gauss_mstep <- function(Y, R, s2_min) {
        nk <- pmax(colSums(R), .Machine$double.eps)
        mu <- crossprod(R, Y) / nk
        s2 <- matrix(0, ncol(R), ncol(Y))
        for (k in seq_len(ncol(R))) {
                Z <- sweep(Y, 2, mu[k, ])^2
                s2[k, ] <- colSums(R[, k] * Z) / nk[k]
        }
        list(mu = mu, s2 = sweep(s2, 2, s2_min, pmax))
}

#' Fit the two-channel (composition x configuration) mixture
#'
#' @param X Count matrix (N x C) from \code{.sig_to_counts()}.
#' @param Y Configuration features (N x L), no NA.
#' @param K Number of mosaic types.
#' @param tau DMM temperature, in (0, 1].
#' @param w Weight of the configuration channel, > 0.
#' @param n_start Starts. The first three start from k-means on both channels,
#'   on configuration only and on composition only; the rest on both.
#' @return A \code{pulse_dmm} object with extra \code{mu}, \code{s2}, \code{w},
#'   \code{channels = "composition+configuration"}. \code{loglik} is the
#'   untempered composite log-likelihood log DM + (w / tau) log N; BIC uses
#'   tau * loglik, as for the one-channel DMM.
#' @keywords internal
#' @noRd
.fit_dmm2 <- function(X, Y, K, tau = 1, w = 1, n_start = 5L, max_iter = 500L, tol = 1e-5,
                      A_init = 50, alpha_min = 1e-4, A_max = 1e6, seed = NULL) {
        X <- as.matrix(X); Y <- as.matrix(Y)
        if (nrow(X) != nrow(Y)) stop("`X` and `Y` need the same rows.", call. = FALSE)
        if (anyNA(Y) || anyNA(X)) stop("`X` and `Y` must have no NA.", call. = FALSE)
        if (any(rowSums(X) <= 0)) stop("Every row of `X` needs a positive total.", call. = FALSE)
        if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau > 1)
                stop("`tau` must be one number in (0, 1].", call. = FALSE)
        if (!is.numeric(w) || length(w) != 1L || !is.finite(w) || w <= 0)
                stop("The configuration weight must be one positive number.", call. = FALSE)
        .local_seed(seed)

        tr <- .dmm_triplets(X)
        N  <- tr$N; eps <- .Machine$double.eps
        s2_min <- pmax(1e-3 * apply(Y, 2, stats::var), 1e-8)

        ## Initialisation embeddings, each block scaled to total variance 1.
        H  <- sqrt(X / tr$n)
        Hs <- H / sqrt(max(sum(apply(H, 2, stats::var)), eps))
        Ys <- scale(Y); Ys[!is.finite(Ys)] <- 0; Ys <- Ys / sqrt(ncol(Y))
        emb <- list(both = cbind(Hs, Ys), conf = Ys, comp = Hs)
        n_dist <- function(M) nrow(unique(round(M, 10)))

        best <- NULL
        for (st in seq_len(n_start)) {
                E <- emb[[c("both", "conf", "comp")[min(st, 3L)]]]
                if (n_dist(E) < K) E <- emb$both
                lab <- if (K == 1L) rep(1L, N) else
                        stats::kmeans(E, centers = K, nstart = 1L, iter.max = 50L,
                                      algorithm = "MacQueen")$cluster
                R <- matrix(0, N, K); R[cbind(seq_len(N), lab)] <- 1
                Th <- crossprod(R, X) + 1e-10; Th <- Th / rowSums(Th)
                alpha <- pmax(Th * A_init, alpha_min)
                gm <- .gauss_mstep(Y, R, s2_min)
                pi_k <- pmax(colMeans(R), eps)
                trace <- numeric(0); converged <- FALSE
                for (it in seq_len(max_iter)) {
                        ld <- .dmm_logdens(tr, alpha); lg <- .gauss_logdens(Y, gm$mu, gm$s2)
                        L  <- ld + (w / tau) * lg + rep(log(pi_k), each = N)
                        lse <- .lse_rows(L); R <- exp(L - lse)
                        trace[it] <- sum(lse)
                        if (it > 1L && abs(trace[it] - trace[it - 1L]) < tol * abs(trace[it])) {
                                converged <- TRUE; break
                        }
                        pi_k <- pmax(colMeans(R), eps)
                        for (k in seq_len(K))
                                alpha[k, ] <- .dmm_update_alpha(tr, R[, k], alpha[k, ],
                                                                alpha_min = alpha_min, A_max = A_max)
                        gm <- .gauss_mstep(Y, R, s2_min)
                }
                fit <- list(R = R, alpha = alpha, pi = pi_k, mu = gm$mu, s2 = gm$s2,
                            loglik = trace[length(trace)], trace = trace,
                            converged = converged, start = st)
                if (is.null(best) || fit$loglik > best$loglik) best <- fit
        }
        .new_pulse_dmm2(best, tr, Y, K, tau, w, colnames(X), colnames(Y), attr(X, "valid"))
}

## Assemble the two-channel fit. df counts alpha, Gaussian means and variances,
## and mixing weights.
.new_pulse_dmm2 <- function(fit, tr, Y, K, tau, w, colnames_X = NULL, colnames_Y = NULL,
                            valid = NULL) {
        L   <- tau * .dmm_logdens(tr, fit$alpha) + w * .gauss_logdens(Y, fit$mu, fit$s2) +
                rep(log(fit$pi), each = tr$N)
        R   <- exp(L - .lse_rows(L))
        df  <- K * tr$C + 2L * K * ncol(Y) + (K - 1L)
        BIC <- -2 * tau * fit$loglik + df * log(tr$N)
        EN  <- -sum(R[R > 0] * log(R[R > 0]))
        colnames(fit$alpha) <- colnames_X
        colnames(fit$mu) <- colnames(fit$s2) <- colnames_Y
        structure(list(U = R, cluster = max.col(R, ties.method = "first"),
                       alpha = fit$alpha, pi = fit$pi, mu = fit$mu, s2 = fit$s2,
                       K = K, tau = tau, w = w, loglik = fit$loglik, BIC = BIC,
                       ICL = BIC + 2 * EN, df = df,
                       concentration = rowSums(fit$alpha), limit = "dirichlet-multinomial",
                       channels = "composition+configuration",
                       trace = fit$trace, converged = fit$converged, valid = valid),
                  class = "pulse_dmm")
}

## Memberships for new rows under a two-channel fit.
.score_dmm2 <- function(fit, X_new, Y_new) {
        X_new <- as.matrix(X_new); Y_new <- as.matrix(Y_new)
        if (ncol(X_new) != ncol(fit$alpha) || ncol(Y_new) != ncol(fit$mu))
                stop("Column mismatch with the fitted model.", call. = FALSE)
        tr <- .dmm_triplets(X_new)
        L  <- fit$tau * .dmm_logdens(tr, fit$alpha) + fit$w * .gauss_logdens(Y_new, fit$mu, fit$s2) +
                rep(log(fit$pi), each = tr$N)
        lse <- .lse_rows(L); U <- exp(L - lse)
        list(U = U, cluster = max.col(U, ties.method = "first"), logpred = lse, loglik = sum(lse))
}
