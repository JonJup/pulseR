# Minimal pulseR objects for crosswalk tests: only the fields
# crosswalk_types() reads.
cw_onehot <- function(h, K = max(h)) {
        U <- matrix(0, length(h), K)
        U[cbind(seq_along(h), h)] <- 1
        U
}
cw_regions <- function(U, hard = max.col(U, ties.method = "first"))
        structure(list(memberships = U, hard_clusters = hard), class = "pulse_regions")
cw_local <- function(U, region, hard = max.col(U, ties.method = "first"))
        structure(list(memberships = U, hard_classes = hard,
                       key = data.frame(name = paste0("local", region, "_",
                                                      stats::ave(region, region, FUN = seq_along)),
                                        region = region)),
                  class = "pulse_local_types")
cw_mosaic <- function(U, valid, hard = max.col(U, ties.method = "first"))
        structure(list(memberships = U, hard_clusters = hard, valid = valid),
                  class = "mosaic_types")
cw_ext <- function(id, lab, ...) data.frame(ID = id, type = lab, ..., stringsAsFactors = FALSE)
