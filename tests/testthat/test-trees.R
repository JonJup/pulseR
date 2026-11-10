test_that("sample_spanning_trees() works without eid and is reproducible", {
        set.seed(1)
        g <- igraph::make_lattice(c(8, 8))
        igraph::E(g)$weight <- stats::runif(igraph::ecount(g))   # no eid on purpose
        set.seed(99); before <- .Random.seed
        a <- sample_spanning_trees(g, n = 5, seed = 42, verbose = FALSE)
        expect_identical(.Random.seed, before)
        b <- sample_spanning_trees(g, n = 5, seed = 42, verbose = FALSE)
        ids <- function(x) lapply(seq_len(x$n_trees), function(j) sort(x$edge_ids[, j]))
        expect_identical(ids(a), ids(b))
        expect_true(all(lengths(ids(a)) == igraph::vcount(g) - 1L))
})

test_that("a non-positional eid is overwritten, not used to map edges", {
        set.seed(2)
        g <- igraph::make_lattice(c(6, 6))
        igraph::E(g)$weight <- stats::runif(igraph::ecount(g))
        ref <- sample_spanning_trees(g, n = 3, seed = 1, verbose = FALSE)   # no eid
        igraph::E(g)$eid <- rev(seq_len(igraph::ecount(g)))
        tr  <- sample_spanning_trees(g, n = 3, seed = 1, verbose = FALSE)   # reversed eid
        expect_identical(tr$edge_ids, ref$edge_ids)
})

test_that(".subtree_sizes() matches brute force on random spanning trees", {
        set.seed(11)
        for (rep in 1:5) {
                g0 <- igraph::sample_gnp(40, 0.15)
                cc <- igraph::components(g0)
                g  <- igraph::induced_subgraph(g0, which(cc$membership == which.max(cc$csize)))
                igraph::E(g)$eid <- seq_len(igraph::ecount(g))
                tr  <- igraph::mst(g, weights = stats::runif(igraph::ecount(g)))
                ep  <- igraph::ends(g, igraph::E(g), names = FALSE)
                eid <- igraph::E(tr)$eid
                got <- .subtree_sizes(eid, ep, igraph::vcount(g))
                root <- ep[eid[1], 1]
                ref <- vapply(seq_along(eid), function(j) {
                        cmp <- igraph::components(igraph::delete_edges(tr, j))$membership
                        sum(cmp != cmp[root])
                }, numeric(1))
                expect_equal(got, ref)
        }
})

test_that("a disconnected graph stops", {
        g <- igraph::disjoint_union(igraph::make_lattice(c(4, 4)), igraph::make_ring(3))
        igraph::E(g)$weight <- 1
        expect_error(sample_spanning_trees(g, n = 2, verbose = FALSE),
                     "2 connected components")
})
make_eid <- function(g) { igraph::E(g)$eid <- seq_len(igraph::ecount(g)); g }

test_that("eta calibration returns 0 and warns on a tree graph", {
        .local_seed(1)
        g <- make_eid(igraph::make_tree(20, children = 2, mode = "undirected"))
        w <- stats::runif(igraph::ecount(g))
        expect_warning(eta <- .calibrate_eta_empirical(g, w, s = 1, target = 0.1), "tree")
        expect_identical(eta, 0)
})

test_that("eta calibration warns when the target is unreachable", {
        .local_seed(1)
        g <- make_eid(igraph::make_ring(30))   # 1 cycle: max turnover 1/29
        w <- stats::runif(igraph::ecount(g))
        expect_warning(.calibrate_eta_empirical(g, w, s = 1, target = 0.5),
                       "cannot exceed")
})

test_that("eta calibration is silent when it converges", {
        .local_seed(1)
        g <- make_eid(igraph::make_lattice(c(8, 8)))
        w <- stats::runif(igraph::ecount(g))
        expect_no_warning(.calibrate_eta_empirical(g, w, s = stats::mad(w), target = 0.1))
})