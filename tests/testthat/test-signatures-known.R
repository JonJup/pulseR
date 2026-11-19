## Level 3 known answers (round-5 review, sections 4.2 and 4.3).

.path_P <- function() rbind(c(1, 0), c(0.5, 0.5), c(0, 1), c(0.25, 0.75))

test_that("signatures on a 4-node path match a hand count", {
        g <- igraph::make_lattice(4)                       # path 1-2-3-4
        P <- .path_P()
        # Edge contributions to (c1:c1, c1:c2, c2:c2):
        #   e12 = (.5, .5, 0), e23 = (0, .5, .5), e34 = (0, .25, .75).
        # An edge counts for focal f when both ends are in N(f).
        s <- compute_all_signatures(P, g, k_order = 1, verbose = FALSE)
        C1 <- rbind(c(0.5, 0.5,  0),                        # N(1) = {1,2}:   e12
                    c(0.5, 1,    0.5),                      # N(2) = {1,2,3}: e12 + e23
                    c(0,   0.75, 1.25),                     # N(3) = {2,3,4}: e23 + e34
                    c(0,   0.25, 0.75))                     # N(4) = {3,4}:   e34
        expect_identical(colnames(s), c("c1:c1", "c1:c2", "c2:c2"))
        expect_equal(unclass(s)[, ], C1 / rowSums(C1), ignore_attr = TRUE)
        expect_equal(attr(s, "totals"), rowSums(C1))
        expect_equal(attr(s, "edge_multiplicity"), 2)       # every edge in 2 neighbourhoods

        s2 <- compute_all_signatures(P, g, k_order = 2, verbose = FALSE)
        C2 <- rbind(c(0.5, 1,    0.5),                      # {1,2,3}: e12 + e23
                    c(0.5, 1.25, 1.25),                     # all edges
                    c(0.5, 1.25, 1.25),                     # all edges
                    c(0,   0.75, 1.25))                     # {2,3,4}: e23 + e34
        expect_equal(unclass(s2)[, ], C2 / rowSums(C2), ignore_attr = TRUE)
        expect_equal(attr(s2, "totals"), rowSums(C2))
        expect_equal(attr(s2, "edge_multiplicity"), 10 / 3)
})

test_that("swapping class columns swaps the matching signature bins", {
        g  <- igraph::make_lattice(4)
        s  <- compute_all_signatures(.path_P(), g, k_order = 1, verbose = FALSE)
        sw <- compute_all_signatures(.path_P()[, 2:1], g, k_order = 1, verbose = FALSE)
        # new c1:c1 = old c2:c2, c1:c2 unchanged
        expect_equal(unclass(sw)[, c(3, 2, 1)], unclass(s)[, ], ignore_attr = TRUE)
})

test_that("renumbering the units permutes the signature rows", {
        set.seed(1)
        g <- igraph::make_lattice(c(6, 6))
        P <- matrix(stats::runif(36 * 3), 36, 3); P <- P / rowSums(P)
        perm <- sample.int(36)                           # vertex i becomes perm[i]
        g2 <- igraph::permute(g, perm)
        P2 <- P; P2[perm, ] <- P
        s  <- compute_all_signatures(P,  g,  k_order = 2, verbose = FALSE)
        s2 <- compute_all_signatures(P2, g2, k_order = 2, verbose = FALSE)
        expect_equal(unclass(s2)[perm, ], unclass(s)[, ], ignore_attr = TRUE)
        expect_equal(attr(s2, "totals")[perm], attr(s, "totals"))
})

test_that("check_signature_saturation() counts pairs near the metric ceiling", {
        ceil <- sqrt(log(2))
        d <- stats::as.dist(matrix(c(0,    ceil,        0.1,
                                     ceil, 0,           0.96 * ceil,
                                     0.1,  0.96 * ceil, 0), 3))
        expect_equal(check_signature_saturation(d, "js"), 2 / 3)
        expect_equal(check_signature_saturation(d / ceil, "hellinger"), 2 / 3)
        expect_equal(check_signature_saturation(d, "js", thresh = 0.99), 1 / 3)
        expect_true(is.na(check_signature_saturation(d, "euclidean")))
})

test_that("compute_mosaic_contrast() matches a hand count and skips NA rows", {
        g <- igraph::make_lattice(4)
        S <- rbind(c(1, 0), c(0, 1), c(1, 0), c(1, 0))
        # mean |s_i - mean of neighbours' s|: 1, 1, 0.5, 0
        expect_equal(compute_mosaic_contrast(S, g), c(1, 1, 0.5, 0))
        S[4, ] <- NA                                     # unit 3 now sees only unit 2
        expect_equal(compute_mosaic_contrast(S, g), c(1, 1, 1, NA))
})
