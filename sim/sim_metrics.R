# =========================================================================== *
# sim_metrics.R
#
# Label-agnostic agreement measures and membership-quality measures.
# Estimated and true labels never need one-to-one matching: memberships are
# FOLDED many-to-one (each estimated class goes to the true class it overlaps
# most), which handles over- and under-segmentation without discarding mass.
# =========================================================================== *

#' Adjusted Rand index. NA where undefined (one class on both sides).
ari <- function(a, b) {
    ok <- !is.na(a) & !is.na(b)
    tab <- table(a[ok], b[ok])
    n   <- sum(tab)
    if (n < 2L) return(NA_real_)
    c2  <- function(x) sum(x * (x - 1) / 2)
    sij <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab))
    e   <- sa * sb / (n * (n - 1) / 2)
    den <- (sa + sb) / 2 - e
    if (den == 0) return(NA_real_)
    (sij - e) / den
}

#' Homogeneity, completeness and V-measure (Rosenberg & Hirschberg 2007).
#' truth = C, est = K. Homogeneity: each estimated class holds one true class.
#' Completeness: each true class sits in one estimated class.
v_measure <- function(truth, est) {
    ok  <- !is.na(truth) & !is.na(est)
    tab <- table(truth[ok], est[ok]); n <- sum(tab)
    H   <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
    HC  <- H(rowSums(tab) / n); HK <- H(colSums(tab) / n)
    HCK <- H(tab / n) - HK                       # H(C | K)
    HKC <- H(tab / n) - HC                       # H(K | C)
    h <- if (HC == 0) 1 else 1 - HCK / HC
    c <- if (HK == 0) 1 else 1 - HKC / HK
    v <- if (h + c == 0) 0 else 2 * h * c / (h + c)
    c(homogeneity = h, completeness = c, v_measure = v)
}

#' Map each estimated label to the true label it overlaps most.
#' Returns a named vector: names = estimated labels, values = true labels.
fold_map <- function(est, truth) {
    tab <- table(est, truth)
    stats::setNames(colnames(tab)[max.col(tab, ties.method = "first")], rownames(tab))
}

#' Fold an N x k membership matrix onto the true classes.
#' @param P N x k matrix, columns in the order of `est_levels`.
#' @param map Output of fold_map().
#' @param true_levels Column order of the result.
fold_columns <- function(P, map, est_levels, true_levels) {
    F <- matrix(0, nrow(P), length(true_levels), dimnames = list(NULL, true_levels))
    for (j in seq_along(est_levels)) {
        tgt <- map[as.character(est_levels[j])]
        if (!is.na(tgt)) F[, tgt] <- F[, tgt] + P[, j]
    }
    F
}

#' AUC: probability that a random error unit scores higher than a random
#' correct unit (Mann-Whitney). score = uncertainty, y = TRUE for errors.
auc <- function(score, y) {
    ok <- is.finite(score) & !is.na(y)
    score <- score[ok]; y <- as.logical(y[ok])
    n1 <- sum(y); n0 <- sum(!y)
    if (n1 == 0L || n0 == 0L) return(NA_real_)
    r <- rank(score)
    (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Multiclass Brier score on folded memberships.
brier <- function(F, truth_idx) {
    ok <- stats::complete.cases(F)
    Y  <- matrix(0, nrow(F), ncol(F)); Y[cbind(seq_len(nrow(F)), truth_idx)] <- 1
    mean(rowSums((F[ok, , drop = FALSE] - Y[ok, , drop = FALSE])^2))
}

#' Expected calibration error of the top-class confidence, equal-width bins.
ece <- function(F, truth_idx, n_bins = 10L) {
    ok   <- stats::complete.cases(F)
    conf <- apply(F[ok, , drop = FALSE], 1, max)
    acc  <- max.col(F[ok, , drop = FALSE], "first") == truth_idx[ok]
    bin  <- cut(conf, seq(0, 1, length.out = n_bins + 1L), include.lowest = TRUE)
    gap  <- tapply(acc, bin, mean) - tapply(conf, bin, mean)
    w    <- tabulate(bin, n_bins) / length(conf)
    sum(w * abs(gap), na.rm = TRUE)
}

#' Reliability table (for pooled reliability diagrams across replicates).
reliability_bins <- function(conf, correct, n_bins = 10L) {
    bin <- cut(conf, seq(0, 1, length.out = n_bins + 1L), include.lowest = TRUE)
    data.frame(bin = levels(bin), n = tabulate(bin, n_bins),
               conf = as.numeric(tapply(conf, bin, mean)),
               acc  = as.numeric(tapply(correct, bin, mean)))
}

#' Mean ARI of `est` against `truth` inside each group, size-weighted.
ari_within <- function(est, truth, group) {
    gs <- split(seq_along(group), group)
    a  <- vapply(gs, function(i) ari(est[i], truth[i]), numeric(1))
    w  <- lengths(gs)
    stats::weighted.mean(a, w, na.rm = TRUE)
}

#' Share of pieces of size 1, averaged over ensemble columns.
singleton_share <- function(E) {
    if (is.null(E)) return(NA_real_)
    mean(apply(E, 2, function(v) { s <- tabulate(v); s <- s[s > 0]; mean(s == 1L) }))
}

onehot <- function(z, levels = sort(unique(z))) {
    M <- matrix(0, length(z), length(levels), dimnames = list(NULL, levels))
    M[cbind(seq_along(z), match(z, levels))] <- 1
    M
}

#' Extra connected pieces summed over classes. 0 = every class contiguous.
excess_pieces <- function(lab, g) {
    sp <- split(seq_along(lab), lab)
    sum(vapply(sp, function(i)
        as.numeric(igraph::components(igraph::induced_subgraph(g, i))$no) - 1, numeric(1)))
}
