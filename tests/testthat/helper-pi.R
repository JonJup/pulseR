# Test helpers for the Pi membership tests. Loaded automatically by testthat.

# Brute-force Pi: loop over every unit and every partition, count piece-mates
# by region, average the shares. Slow and obvious; the reference for the
# sparse implementation in .ensemble_stats().
pi_bruteforce <- function(memb, lab, k) {
        N  <- nrow(memb)
        Pi <- matrix(0, N, k)
        Mi <- integer(N)
        for (t in seq_len(ncol(memb))) {
                for (i in seq_len(N)) {
                        mates <- setdiff(which(memb[, t] == memb[i, t]), i)
                        if (!length(mates)) next
                        Pi[i, ] <- Pi[i, ] + tabulate(lab[mates], k) / length(mates)
                        Mi[i]   <- Mi[i] + 1L
                }
        }
        Pi <- Pi / pmax(Mi, 1)
        Pi[Mi == 0, ] <- NA_real_
        list(Pi = Pi, pi_support = Mi)
}

# Labels for contiguous pieces on a chain 1..N, given cut positions
# (a cut at c separates unit c from unit c + 1).
chain_labels <- function(N, cuts) {
        lab <- rep(1L, N)
        for (c in sort(cuts)) lab[(c + 1L):N] <- lab[(c + 1L):N] + 1L
        lab
}

# Random contiguous chain ensemble: M columns, R pieces each.
chain_ensemble <- function(N, M, R, seed = 1) {
        set.seed(seed)
        sapply(seq_len(M), function(t) chain_labels(N, sample.int(N - 1L, R - 1L)))
}
