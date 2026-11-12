## tests/testthat/helper-golden.R
## Fixed-seed memberships used to guard the edge-id refactor.
.golden_fixture <- function(prior = FALSE, side = 20, n_tree = 15) {
        set.seed(1)
        g   <- igraph::make_lattice(c(side, side))
        xy  <- expand.grid(x = seq_len(side), y = seq_len(side))
        env <- cbind(cut(xy$x, 3, labels = FALSE) * 2 + stats::rnorm(side^2, sd = 0.5),
                    stats::rnorm(side^2))
        ep  <- igraph::ends(g, igraph::E(g), names = FALSE)
        igraph::E(g)$weight <- sqrt(rowSums((env[ep[, 1], ] - env[ep[, 2], ])^2))
        g   <- igraph::set_graph_attr(g, "env", env)
        pr  <- if (prior) ifelse(xy$y <= side / 2, 1L, 2L) else NULL
        list(g = g, tr = sample_spanning_trees(g, n = n_tree, seed = 1, verbose = FALSE,
                                                prior_typology = pr))
}

.golden_run <- function() {
        out <- list()
        for (p in c(FALSE, TRUE)) {
                x <- .golden_fixture(prior = p)
                for (rule in c("propensity", "ssd")) for (k in c(2, 6, 12)) {
                        set.seed(7)
                        out[[paste(p, rule, k)]] <- .compute_ensemble_internal(
                                x$g, x$tr, 15, k,
                                prior_strength = if (p) 1 else 0, size_influence = 1,
                                cut_rule = rule, min_size = if (rule == "ssd") 3L else 1L)
                }
                ccp <- check_cut_placement(x$tr, 6, weight_exponent = 1, seed = 3)$cuts
                # weight_exponent column added after the fixture was built;
                # values of all other columns must stay identical.
                out[[paste(p, "ccp")]] <- ccp[setdiff(names(ccp), "weight_exponent")]
        }
        structure(out, igraph = as.character(utils::packageVersion("igraph")))
}
