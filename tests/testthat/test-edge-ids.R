## tests/testthat/test-edge-ids.R
## pulse_trees stores parent-edge ids, not igraph trees.

test_that("edge-id format reproduces pre-refactor memberships under a fixed seed", {
        ref <- readRDS(test_path("fixtures", "memb_golden.rds"))
        skip_if_not(identical(attr(ref, "igraph"),
                              as.character(utils::packageVersion("igraph"))),
                    "fixture built with another igraph version (mst edge order may differ)")
        expect_identical(.golden_run(), ref)
})

test_that("pulse_trees holds an integer (N - 1) x n edge-id matrix and no igraph objects", {
        x  <- .golden_fixture(side = 10, n_tree = 7)
        id <- x$tr$edge_ids
        expect_null(x$tr$trees)
        expect_identical(typeof(id), "integer")
        expect_identical(dim(id), c(99L, 7L))
        expect_false(any(vapply(x$tr, inherits, NA, "igraph")))
        expect_true(all(apply(id, 2, anyDuplicated) == 0L))
        for (j in seq_len(ncol(id)))                       # each column spans the graph
                expect_identical(max(.labels_from_edges(
                        x$tr$edge_endpoints[id[, j], , drop = FALSE], 100L)), 1L)
})

test_that("an old-format pulse_trees stops with a re-run hint", {
        x <- .golden_fixture(side = 6, n_tree = 2)
        x$tr$edge_ids <- NULL
        expect_error(.compute_ensemble_internal(x$g, x$tr, 2, 2), "Re-run sample_spanning_trees")
})
