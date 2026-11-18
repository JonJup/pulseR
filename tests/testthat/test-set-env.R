## B9: graphs whose edge weights come from elsewhere (e.g. GDM).

.b9_graphs <- function() {
        x   <- make_grid_sf()
        net <- suppressMessages(build_network(x, verbose = FALSE))
        g   <- add_edge_weight(net, id_col = "ID", verbose = FALSE)
        ## "Foreign" weights: same graph, env matrix removed.
        g0  <- igraph::delete_graph_attr(g, "env")
        list(x = x, env = sf::st_drop_geometry(x), g = g, g0 = g0)
}

.collect_messages <- function(expr) {
        msgs <- character()
        res  <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                                    message = function(m) {
                                            msgs <<- c(msgs, conditionMessage(m))
                                            invokeRestart("muffleMessage")
                                    })
        list(result = res, messages = msgs)
}

test_that("set_env() reproduces the matrix add_edge_weight() stores", {
        b  <- .b9_graphs()
        g1 <- set_env(b$g0, b$env, scale = TRUE, verbose = FALSE)
        expect_equal(unname(igraph::graph_attr(g1, "env")),
                     unname(igraph::graph_attr(b$g, "env")))
        expect_identical(colnames(igraph::graph_attr(g1, "env")), c("temp", "prec"))
})

test_that("rows are matched by ID; order and extra rows do not matter", {
        b   <- .b9_graphs()
        ref <- igraph::graph_attr(set_env(b$g0, b$env, verbose = FALSE), "env")
        set.seed(2)
        shuf  <- b$env[sample(nrow(b$env)), ]
        extra <- rbind(shuf, data.frame(ID = 1000:1002, temp = 0, prec = 0))
        got   <- igraph::graph_attr(set_env(b$g0, extra, env_id = "ID",
                                            verbose = FALSE), "env")
        expect_identical(got, ref)
        mat <- as.matrix(shuf)
        expect_identical(igraph::graph_attr(set_env(b$g0, mat, verbose = FALSE), "env"), ref)
        expect_identical(igraph::graph_attr(set_env(b$g0, b$x, verbose = FALSE), "env"), ref)
})

test_that("set_env() stops on bad input with a message that names the problem", {
        b <- .b9_graphs()
        expect_error(set_env(b$g0, b$env[-1, ], verbose = FALSE), "no environmental row")
        dup <- b$env; dup$ID[2] <- dup$ID[1]
        expect_error(set_env(b$g0, dup, verbose = FALSE), "duplicated ID")
        na <- b$env; na$temp[3] <- NA
        expect_error(set_env(b$g0, na, verbose = FALSE), "1 missing value.*temp")
        expect_error(set_env(b$g0, b$env, variables = "nope", verbose = FALSE), "not found")
        chr <- b$env; chr$lab <- "a"
        expect_error(set_env(b$g0, chr, variables = c("temp", "lab"), verbose = FALSE),
                     "not numeric: lab")
        expect_error(set_env(b$g0, b$env, env_id = "nope", verbose = FALSE), "env_id")
        expect_error(set_env(b$g0, list(1), verbose = FALSE), "data frame")
})

test_that("without IDs, row order is assumed and a warning says so", {
        b  <- .b9_graphs()
        g2 <- igraph::delete_vertex_attr(b$g0, "polygon_id")
        expect_warning(set_env(g2, b$env, verbose = FALSE), "polygon_id")
        expect_warning(set_env(b$g0, b$env[, c("temp", "prec")], verbose = FALSE),
                       "ID column")
})

test_that("set_env() reports whether the null baseline will accept the weights", {
        b <- .b9_graphs()
        m_own <- .collect_messages(set_env(b$g0, b$env, scale = TRUE))$messages
        expect_true(any(grepl("Correlation.*: 1\\.000\\.\\s*$", m_own)))
        gf <- b$g0
        igraph::E(gf)$weight <- stats::runif(igraph::ecount(gf))
        m_gdm <- .collect_messages(set_env(gf, b$env, scale = TRUE))$messages
        expect_true(any(grepl("will stop", m_gdm)))
        expect_true(any(grepl("Replacing", .collect_messages(set_env(b$g, b$env))$messages)))
})

test_that("the defaults stop before any sampling when the graph has no env", {
        b   <- .b9_graphs()
        out <- .collect_messages(get_regions(b$g0, final_regions = 2, N_tree = 8,
                                             verbose = TRUE))
        expect_s3_class(out$result, "error")
        expect_match(conditionMessage(out$result), "set_env\\(\\)")
        expect_match(conditionMessage(out$result),
                     "cut_rule = \"ssd\" and refine_boundary = TRUE need")
        expect_false(any(grepl("spanning tree", out$messages, ignore.case = TRUE)))
        expect_error(get_regions(b$g0, final_regions = 2, N_tree = 8,
                                 cut_rule = "propensity", verbose = FALSE),
                     "refine_boundary = TRUE needs")
        tr <- sample_spanning_trees(b$g0, n = 8, seed = 1, verbose = FALSE)
        expect_error(tune_regions(b$g0, tr, N_tree = 8, final_regions = 2:3,
                                  verbose = FALSE), "set_env")
        expect_error(compute_ensemble(tr, N_tree = 8, graph = b$g0,
                                      regions_per_partition = 6, verbose = FALSE),
                     "cut_rule = \"ssd\" needs")
        ens <- compute_ensemble(tr, N_tree = 8, graph = b$g0, regions_per_partition = 6,
                                cut_rule = "propensity", verbose = FALSE)
        expect_error(compute_consensus(ens, 2, graph = b$g0, verbose = FALSE),
                     "refine_boundary = TRUE needs")
})

test_that("a malformed stored env stops with a specific message", {
        b <- .b9_graphs()
        X <- igraph::graph_attr(b$g, "env")
        bad_rows <- igraph::set_graph_attr(b$g0, "env", X[-1, ])
        expect_error(.require_env(bad_rows, "ssd"), "63 rows, `graph` has 64")
        bad_df <- igraph::set_graph_attr(b$g0, "env", as.data.frame(X))
        expect_error(.require_env(bad_df, "ssd"), "numeric matrix, not data.frame")
        X[2, 1] <- NA
        bad_na <- igraph::set_graph_attr(b$g0, "env", X)
        expect_error(.require_env(bad_na, refine_boundary = TRUE), "1 missing")
        expect_null(.require_env(b$g0, "propensity", FALSE))
})

test_that("weights from elsewhere run under the defaults after set_env()", {
        b  <- .b9_graphs()
        gf <- b$g0
        ep <- igraph::ends(gf, igraph::E(gf), names = FALSE)
        igraph::E(gf)$weight <- abs(b$env$temp[ep[, 1]] - b$env$temp[ep[, 2]]) + 0.01
        gf  <- set_env(gf, b$env, scale = TRUE, verbose = FALSE)
        res <- suppressWarnings(get_regions(gf, final_regions = 2, N_tree = 8,
                                            seed = 1, verbose = FALSE))
        expect_s3_class(res, "pulse_regions")
        expect_equal(unname(rowSums(res$memberships)), rep(1, 64), tolerance = 1e-8)
        ## Edge weights only: no env needed.
        res2 <- suppressWarnings(get_regions(b$g0, final_regions = 2, N_tree = 8,
                                             cut_rule = "propensity",
                                             refine_boundary = FALSE,
                                             seed = 1, verbose = FALSE))
        expect_s3_class(res2, "pulse_regions")
})
