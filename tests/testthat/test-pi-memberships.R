# Pi memberships (Section 6 of PULSE_the_method_as_it_stands.md)
#
# Pi_ik = share of unit i's piece-mates that sit in consensus region k,
# averaged over the partitions in which i has at least one piece-mate.
# Cohesion is kept as the raw, non-normalised
# co-association and is NOT a probability.
#
# Run with: devtools::test(filter = "pi-memberships")
# The snapshot file tests/testthat/_snaps/pi-memberships.md shows the worked
# examples as printed tables.

# ---------------------------------------------------------------------------
# 1. The definition, checked by hand
# ---------------------------------------------------------------------------

test_that("Pi matches a hand calculation on a 5-unit, 2-partition ensemble", {
        # Partition 1: {1,2,3} {4,5}    Partition 2: {1,2} {3,4,5}
        # Consensus regions: A = {1,2,3}, B = {4,5}
        memb <- cbind(c(1, 1, 1, 2, 2), c(1, 1, 2, 2, 2))
        lab  <- c(1L, 1L, 1L, 2L, 2L)
        st   <- .ensemble_stats(memb, lab, k = 2L, pi = TRUE)
        
        # Unit 3: P1 mates {1,2} -> (1, 0); P2 mates {4,5} -> (0, 1); mean (.5, .5)
        # Unit 4: P1 mates {5}   -> (0, 1); P2 mates {3,5} -> (.5, .5); mean (.25, .75)
        expected <- rbind(c(1, 0), c(1, 0), c(.5, .5), c(.25, .75), c(.25, .75))
        expect_equal(st$Pi, expected)
        expect_equal(st$pi_support, rep(2L, 5))
})

test_that("Pi equals a brute-force loop on random ensembles, for any block size", {
        for (seed in 1:3) {
                set.seed(seed)
                N <- 40; M <- 15; R <- 6; k <- 4
                memb <- matrix(sample.int(R, N * M, replace = TRUE), N, M)
                lab  <- sample(rep_len(seq_len(k), N))
                ref  <- pi_bruteforce(memb, lab, k)
                
                one_block   <- .ensemble_stats(memb, lab, k, pi = TRUE)
                many_blocks <- .ensemble_stats(memb, lab, k, pi = TRUE, budget = 2 * N)
                
                expect_equal(one_block$Pi, ref$Pi)
                expect_equal(one_block$pi_support, ref$pi_support)
                expect_equal(many_blocks$Pi, one_block$Pi)
        }
})

test_that("Pi rows sum to 1; cohesion rows do not", {
        memb <- chain_ensemble(N = 60, M = 30, R = 8)
        lab  <- chain_labels(60, c(20, 45))                 # 3 regions: 20, 25, 15
        st   <- .ensemble_stats(memb, lab, k = 3L, pi = TRUE)
        coh  <- .cohesion_from_stats(st, seq_len(60), lab, 60)
        
        expect_equal(rowSums(st$Pi), rep(1, 60))
        expect_true(all(st$Pi >= 0 & st$Pi <= 1))
        expect_gt(max(abs(rowSums(coh) - 1)), 0.1)
})

# ---------------------------------------------------------------------------
# 2. Why Pi: the small-region bias of normalised cohesion
# ---------------------------------------------------------------------------

test_that("Pi removes the small-region bias of normalised cohesion", {
        # Chain of 45 units. Region L = units 1-40, region S = units 41-45.
        # In every partition, units 35-45 form one piece: 6 L units and 5 S units.
        # The rest of the chain is cut at random.
        N <- 45; M <- 4
        set.seed(42)
        memb <- sapply(seq_len(M), function(t) chain_labels(N, c(sort(sample.int(33, 3)), 34)))
        lab  <- c(rep(1L, 40), rep(2L, 5))
        st   <- .ensemble_stats(memb, lab, k = 2L, pi = TRUE)
        coh  <- .cohesion_from_stats(st, seq_len(N), lab, N)
        U_old <- coh / rowSums(coh)                         # pre-Pi memberships
        
        # Unit 40 (last L unit): piece-mates are 5 L + 5 S, so Pi = (.5, .5).
        expect_equal(st$Pi[40, ], c(.5, .5))
        # Cohesion: 5 of 39 L members always co-occur (5/39); all 5 S members do (1).
        expect_equal(coh[40, ], c(5 / 39, 1))
        # Normalising that gives about (0.11, 0.89): an 89% vote for the small region.
        expect_equal(round(U_old[40, ], 2), c(.11, .89))
        
        shown <- data.frame(
                unit            = c(40, 41),
                label           = c("L", "S"),
                mates_L         = c(5, 6),
                mates_S         = c(5, 4),
                cohesion_L      = round(coh[40:41, 1], 3),
                cohesion_S      = round(coh[40:41, 2], 3),
                old_U_L         = round(U_old[40:41, 1], 3),
                old_U_S         = round(U_old[40:41, 2], 3),
                Pi_L            = round(st$Pi[40:41, 1], 3),
                Pi_S            = round(st$Pi[40:41, 2], 3)
        )
        expect_snapshot(print(shown, row.names = FALSE))
})

# ---------------------------------------------------------------------------
# 3. Edge cases
# ---------------------------------------------------------------------------

test_that("a unit that is a one-unit piece in every partition gets an NA row", {
        # Unit 5 is alone in both partitions.
        memb <- cbind(c(1, 1, 2, 2, 3), c(1, 2, 2, 2, 3))
        lab  <- c(1L, 1L, 2L, 2L, 2L)
        st   <- .ensemble_stats(memb, lab, k = 2L, pi = TRUE)
        
        expect_equal(st$pi_support, c(1L, 2L, 2L, 2L, 0L))
        expect_true(all(is.na(st$Pi[5, ])))
        expect_false(anyNA(st$Pi[1:4, ]))
        expect_equal(st$Pi, pi_bruteforce(memb, lab, 2L)$Pi)
})

test_that("Pi refuses a partial reference labelling", {
        memb <- chain_ensemble(N = 20, M = 5, R = 4)
        expect_error(
                .ensemble_stats(memb, c(1L, 2L), k = 2L, ref_idx = c(1L, 20L), pi = TRUE),
                "full reference")
        # Without Pi the sampled-path call still works as before.
        expect_silent(.ensemble_stats(memb, c(1L, 2L), k = 2L, ref_idx = c(1L, 20L)))
})

test_that("certainty, margin and entropy are read off Pi against the final labels", {
        U <- rbind(c(.70, .20, .10),
                   c(.40, .45, .15),        # label 1, but region 2 is larger
                   c(1/3, 1/3, 1/3),
                   c(NA,  NA,  NA))
        hard <- c(1L, 1L, 2L, 3L)
        sm <- .membership_summaries(U, hard)
        
        expect_equal(sm$certainty, c(.70, .40, 1/3, NA))
        expect_equal(sm$margin,    c(.50, .05, 0, NA))
        expect_equal(sm$entropy[3], 1)
        expect_true(is.na(sm$entropy[4]))
})

# ---------------------------------------------------------------------------
# 4. compute_consensus(): Pi is recomputed AFTER contiguity repair
# ---------------------------------------------------------------------------

test_that("memberships and summaries come from the repaired labels", {
        # Chain of 12 units; true regions 1-8 and 9-12. The consensus step is
        # mocked to return an exclave: unit 3 labelled region 2.
        N <- 12
        g <- igraph::make_lattice(N)
        memb  <- chain_ensemble(N, M = 40, R = 4, seed = 7)
        truth <- c(rep(1L, 8), rep(2L, 4))
        bad   <- truth; bad[3] <- 2L
        local_mocked_bindings(.consensus_partition = mock_partition(bad))
        
        fit <- compute_consensus(memb, final_regions = 2L, graph = g, verbose = FALSE)
        
        # Repair moved unit 3 back.
        expect_equal(fit$hard_labels, truth)
        
        # Pi is the one for the repaired labels, not the mocked ones.
        pi_after  <- .ensemble_stats(memb, truth, 2L, pi = TRUE)$Pi
        pi_before <- .ensemble_stats(memb, bad,   2L, pi = TRUE)$Pi
        expect_equal(fit$memberships, pi_after)
        expect_false(isTRUE(all.equal(fit$memberships, pi_before)))
        
        # Summaries agree with the returned memberships and labels.
        sm <- .membership_summaries(fit$memberships, fit$hard_labels)
        expect_equal(fit$certainty, sm$certainty)
        expect_equal(fit$margin,    sm$margin)
        expect_equal(fit$entropy,   sm$entropy)
        
        expect_true(is.matrix(fit$cohesion))
        
        shown <- data.frame(
                unit             = 1:6,
                label_mocked     = bad[1:6],
                label_repaired   = fit$hard_labels[1:6],
                Pi1_if_unrepaired = round(pi_before[1:6, 1], 3),
                Pi1_returned     = round(fit$memberships[1:6, 1], 3),
                certainty        = round(fit$certainty[1:6], 3)
        )
        expect_snapshot(print(shown, row.names = FALSE))
})

test_that("a unit without piece-mates is filled, one-hot without a graph", {
        memb <- cbind(c(1, 1, 2, 2, 3), c(1, 2, 2, 2, 3))
        local_mocked_bindings(.consensus_partition = mock_partition(c(1L, 1L, 2L, 2L, 2L)))
        expect_warning(fit <- compute_consensus(memb, final_regions = 2L, verbose = FALSE), "pi_support")
        expect_equal(fit$pi_support[5], 0L)
        expect_equal(fit$memberships[5, ], c(0, 1))
        g <- igraph::make_ring(5, circular = FALSE)
        expect_warning(fit <- compute_consensus(memb, final_regions = 2L, graph = g, verbose = FALSE), "pi_support")
        expect_equal(fit$memberships[5, ], fit$memberships[4, ])
})

test_that(".fill_unsupported_pi() fills in synchronous passes", {
        Pi  <- rbind(c(1, 0), c(NA, NA), c(NA, NA), c(0, 1))
        out <- .fill_unsupported_pi(Pi, c(1L, 1L, 2L, 2L), igraph::make_ring(4, circular = FALSE))
        expect_equal(out[2, ], c(1, 0)); expect_equal(out[3, ], c(0, 1))
})

# ---------------------------------------------------------------------------
# 5. Downstream contract
# ---------------------------------------------------------------------------

test_that("pulse_regions and Level 2 accessors read the new fields", {
        N <- 12
        memb <- chain_ensemble(N, M = 40, R = 4, seed = 7)
        truth <- c(rep(1L, 8), rep(2L, 4))
        local_mocked_bindings(.consensus_partition = mock_partition(truth))
        fit <- compute_consensus(memb, final_regions = 2L, verbose = FALSE)
        
        cfg <- list(final_regions = 2L, regions_per_partition = 4L, N_tree = 40L,
                    partitions_per_tree = 1L, weight_exponent = 1)
        cfg <- cfg[.TUNE_AXES]
        reg <- expect_silent(.new_pulse_regions(fit, config = cfg))
        
        expect_equal(reg$memberships, fit$memberships)
        expect_equal(reg$entropy, fit$entropy)
        expect_equal(.region_weights(reg, "memberships"), fit$memberships)
        expect_equal(.region_weights(reg, "cohesion"),  fit$cohesion)
})

test_that("contiguity repair regrows freed units by Pi, not by cohesion", {
        # Chain of 45. Regions: 3 = units 1-5, 1 (L) = 6-39, 2 (S) = 41-45.
        # Unit 40 is mocked into region 3, so it is an exclave and gets freed.
        # Its neighbours are 39 (L) and 41 (S). In every partition units 34-45
        # form one piece: unit 40 has 6 L piece-mates and 5 S piece-mates.
        N <- 45
        g <- igraph::make_lattice(N)
        set.seed(3)
        memb <- sapply(1:6, function(t) chain_labels(N, c(sort(sample(6:32, 3)), 33)))
        bad  <- c(rep(3L, 5), rep(1L, 34), 3L, rep(2L, 5))
        
        st  <- .ensemble_stats(memb, bad, 3L, pi = TRUE)
        coh <- .cohesion_from_stats(st, seq_len(N), bad, N)
        
        # The two scores disagree for unit 40: cohesion votes S, Pi votes L.
        expect_equal(coh[40, 1:2], c(6 / 34, 1))
        expect_equal(st$Pi[40, ], c(6 / 11, 5 / 11, 0))
        expect_equal(.enforce_contiguity(g, bad, coh)[40],    2L)
        expect_equal(.enforce_contiguity(g, bad, st$Pi)[40],  1L)
        
        # compute_consensus() uses Pi.
        local_mocked_bindings(.consensus_partition = mock_partition(bad))
        fit <- compute_consensus(memb, final_regions = 3L, graph = g, verbose = FALSE)
        expect_equal(fit$hard_labels[40], 1L)
        
        shown <- data.frame(
                score          = c("cohesion", "Pi"),
                to_L           = round(c(coh[40, 1], st$Pi[40, 1]), 3),
                to_S           = round(c(coh[40, 2], st$Pi[40, 2]), 3),
                regrown_into   = c("S", "L")
        )
        expect_snapshot(print(shown, row.names = FALSE))
})