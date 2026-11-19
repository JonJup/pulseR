## .refine_regions() on non-planar graphs (round-5 review, section 4.2).
## The connectivity check was only tested on rook lattices. Queen lattices and
## dense random graphs have crossing edges, where "neighbours linked around i"
## is not the same as "region stays connected".

.queen_lattice <- function(n) {
        g  <- igraph::make_lattice(c(n, n))
        id <- function(x, y) x + (y - 1L) * n            # make_lattice: x varies fastest
        xy <- expand.grid(x = seq_len(n - 1L), y = seq_len(n - 1L))
        d  <- rbind(id(xy$x, xy$y), id(xy$x + 1L, xy$y + 1L),
                    id(xy$x + 1L, xy$y), id(xy$x, xy$y + 1L))
        igraph::add_edges(g, as.vector(d))           # pairs (row 1, 2) and (row 3, 4)
}

.dense_random <- function(n, p, seed) {
        set.seed(seed)
        g <- igraph::sample_gnp(n, p)
        g <- igraph::add_edges(g, as.vector(rbind(seq_len(n), c(seq_len(n)[-1], 1L))))
        igraph::simplify(g)                              # ring: always connected
}

# Contiguous start: graph-Voronoi cells around k random seeds, repaired.
.start_partition <- function(g, k, seed) {
        set.seed(seed)
        s <- sample.int(igraph::vcount(g), k)
        z <- max.col(-t(igraph::distances(g, v = s)), ties.method = "random")
        .enforce_contiguity(g, z, diag(k)[z, , drop = FALSE])
}

.ssw <- function(X, z) {
        G <- rowsum(X, z) / tabulate(z)
        sum((X - G[z, , drop = FALSE])^2)
}

test_that("boundary refinement keeps regions connected on non-planar graphs", {
        graphs <- list(queen = .queen_lattice(12), random = .dense_random(120, 0.06, 1))
        moves  <- 0L
        for (nm in names(graphs)) {
                g <- graphs[[nm]]
                n <- igraph::vcount(g)
                expect_gt(igraph::ecount(g), 3 * n - 6)  # more edges than any planar graph
                for (s in 1:4) for (ms in c(1L, 8L)) {
                        set.seed(100 + s)
                        X  <- cbind(stats::rnorm(n), stats::rnorm(n))
                        z0 <- .start_partition(g, 6, s)
                        expect_true(.check_contiguity(g, z0)$contiguous)
                        r  <- .refine_regions(z0, g, X, min_size = ms)
                        info <- sprintf("%s, seed %d, min_size %d", nm, s, ms)
                        expect_true(.check_contiguity(g, r$z)$contiguous, info = info)
                        expect_setequal(unique(r$z), 1:6)
                        expect_true(all(tabulate(r$z, 6) >= pmin(tabulate(z0, 6), ms)), info = info)
                        expect_lte(.ssw(X, r$z), .ssw(X, z0) + 1e-9)
                        moves <- moves + r$moves
                }
        }
        expect_gt(moves, 0L)                             # the test exercised real moves
})
