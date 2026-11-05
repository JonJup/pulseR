test_that("sample_spanning_trees() works without eid and is reproducible", {
        set.seed(1)
        g <- igraph::make_lattice(c(8, 8))
        igraph::E(g)$weight <- stats::runif(igraph::ecount(g))   # no eid on purpose
        set.seed(99); before <- .Random.seed
        a <- sample_spanning_trees(g, n = 5, seed = 42, verbose = FALSE)
        expect_identical(.Random.seed, before)
        b <- sample_spanning_trees(g, n = 5, seed = 42, verbose = FALSE)
        ids <- function(x) lapply(x$trees, function(t) sort(igraph::E(t)$eid))
        expect_identical(ids(a), ids(b))
        expect_true(all(lengths(ids(a)) == igraph::vcount(g) - 1L))
})

test_that("a non-positional eid is overwritten, not used to map edges", {
        set.seed(2)
        g <- igraph::make_lattice(c(6, 6))
        igraph::E(g)$weight <- stats::runif(igraph::ecount(g))
        igraph::E(g)$eid <- rev(seq_len(igraph::ecount(g)))
        tr <- sample_spanning_trees(g, n = 3, seed = 1, verbose = FALSE)
        for (t in tr$trees) {
                ep_tree <- igraph::ends(t, igraph::E(t), names = FALSE)
                ep_par  <- tr$edge_endpoints[igraph::E(t)$eid, , drop = FALSE]
                expect_equal(t(apply(ep_tree, 1, sort)), t(apply(ep_par, 1, sort)))
        }
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
