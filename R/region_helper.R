#' Chance co-association per ensemble column
#'
#' For each column of the ensemble membership matrix, the probability that two units drawn uniformly at random
#' share an tree region. Averaged over columns this is the floor
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
