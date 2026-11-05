test_that("STAB is split-half agreement: high for a stable ensemble, low for noise", {
        set.seed(4)
        N <- 120; base <- rep(1:3, each = 40)
        # Every column refines the same 3 blocks: halves agree.
        stable <- sapply(1:20, function(t) (base - 1L) * 2L + sample.int(2L, N, TRUE))
        # Columns are unrelated random partitions: halves disagree.
        noise  <- sapply(1:20, function(t) sample.int(6L, N, TRUE))
        s1 <- .stability_consensus(stable, 3L, B = 3L, seed = 1)
        s0 <- .stability_consensus(noise,  3L, B = 3L, seed = 1)
        expect_gt(s1, 0.9)
        expect_lt(s0, 0.3)
})

test_that("STAB uses disjoint halves and handles edge cases", {
        memb <- matrix(rep(1:2, each = 10), 20, 3)
        expect_true(is.na(.stability_consensus(memb, 2L, B = 2L)))       # M < 4
        expect_true(is.na(.stability_consensus(cbind(memb, memb), 2L, B = 0L)))
        seen <- list()
        local_mocked_bindings(.consensus_partition = function(memb, k, ...) {
                seen[[length(seen) + 1L]] <<- memb[1, ]
                list(labels = rep(1:2, length.out = nrow(memb)))
        })
        m <- matrix(rep(1:6, each = 10), 10, 6)              # column t holds value t
        .stability_consensus(m, 2L, B = 1L, seed = 3)
        expect_length(seen, 2L)
        expect_length(intersect(seen[[1]], seen[[2]]), 0L)
        expect_setequal(c(seen[[1]], seen[[2]]), 1:6)
})
