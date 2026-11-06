test_that("identical crisp typologies give perfect agreement and 1:1 relations", {
        h  <- rep(1:4, each = 25)
        cw <- crosswalk_types(cw_regions(cw_onehot(h)), cw_ext(1:100, letters[h]),
                              label = "type", units = 1:100)
        a  <- stats::setNames(cw$agreement$crisp, cw$agreement$index)
        expect_equal(unname(a[c("ARI", "AMI", "NMI", "purity")]), c(1, 1, 1, 1))
        expect_equal(unname(a["VI"]), 0)
        expect_true(all(cw$types$pulse$relation == "1:1"))
        expect_true(all(cw$units$match))
        expect_equal(cw$agreement$fuzzy[3:9], cw$agreement$crisp[3:9])
})

test_that("units are matched by id, not row order", {
        .local_seed(1)
        h   <- sample(1:3, 60, TRUE); e <- sample(c("x", "y"), 60, TRUE)
        obj <- cw_regions(cw_onehot(h))
        ext <- cw_ext(1:60, e)
        a1 <- crosswalk_types(obj, ext, "type", units = 1:60)
        a2 <- crosswalk_types(obj, ext[sample(60), ], "type", units = 1:60)
        expect_equal(a1$overlap, a2$overlap)
        expect_equal(a1$agreement, a2$agreement)
        ## character vs integer IDs compare as strings
        ext$ID <- as.character(ext$ID)
        expect_equal(crosswalk_types(obj, ext, "type", units = 1:60)$overlap, a1$overlap)
})

test_that("refinement shows up in the directional VI parts and relations", {
        h_p <- rep(1:6, each = 20)                 # PULSE splits each external type in two
        h_e <- c("A", "B", "C")[ceiling(h_p / 2)]
        cw  <- crosswalk_types(cw_regions(cw_onehot(h_p)), cw_ext(1:120, h_e), "type",
                               units = 1:120)
        a   <- stats::setNames(cw$agreement$crisp, cw$agreement$index)
        expect_equal(unname(a["H_ext_given_pulse"]), 0)
        expect_equal(unname(a["H_pulse_given_ext"]), log(2))
        cr <- function(d) d$relation[d$basis == "crisp"]
        expect_true(all(cr(cw$types$pulse) == "nested"))
        expect_true(all(cr(cw$types$external) == "contains"))
})

test_that("independent labels give chance-level ARI and AMI", {
        .local_seed(2)
        n  <- 2000
        cw <- crosswalk_types(cw_regions(cw_onehot(sample(1:5, n, TRUE))),
                              cw_ext(seq_len(n), sample(letters[1:4], n, TRUE)),
                              "type", units = seq_len(n))
        a  <- stats::setNames(cw$agreement$crisp, cw$agreement$index)
        expect_lt(abs(a["ARI"]), 0.01)
        expect_lt(abs(a["AMI"]), 0.01)
})

test_that("exact EMI equals the mean MI over all permutations", {
        x <- c(1, 1, 1, 2, 2, 3, 3); y <- c(1, 1, 2, 2, 2, 2, 3)
        perms <- function(v) if (length(v) <= 1L) list(v) else
                do.call(c, lapply(seq_along(v), function(i)
                        lapply(perms(v[-i]), function(p) c(v[i], p))))
        mi <- function(a, b) { C <- table(a, b) / length(a)
                .cw_h(rowSums(C)) + .cw_h(colSums(C)) - .cw_h(C) }
        brute <- mean(vapply(perms(y), function(p) mi(x, p), numeric(1)))
        expect_equal(.cw_emi_exact(as.numeric(table(x)), as.numeric(table(y)), 7), brute,
                     tolerance = 1e-12)
})

test_that("AMI falls back to the approximation above ami_max_terms", {
        .local_seed(3)
        h  <- sample(1:3, 300, TRUE)
        e  <- ifelse(stats::runif(300) < 0.7, h, sample(1:3, 300, TRUE))
        f  <- function(m) crosswalk_types(cw_regions(cw_onehot(h)), cw_ext(1:300, e), "type",
                                          units = 1:300, ami_max_terms = m)
        ex <- f(5e7); ap <- f(1)
        expect_equal(ex$settings$ami_method, "exact")
        expect_equal(ap$settings$ami_method, "approx")
        get <- function(cw) cw$agreement$crisp[cw$agreement$index == "AMI"]
        expect_equal(get(ex), get(ap), tolerance = 0.01)
})

test_that("one-hot fuzzy external input reproduces the crisp-label path", {
        .local_seed(4)
        U   <- matrix(stats::runif(150), 50); U <- U / rowSums(U)
        e   <- sample(1:2, 50, TRUE)
        V   <- cw_onehot(e, 2)
        ext <- data.frame(ID = 1:50, type = c("p", "q")[e], p = V[, 1], q = V[, 2])
        a <- crosswalk_types(cw_regions(U), ext, "type", units = 1:50)
        b <- crosswalk_types(cw_regions(U), ext, c("p", "q"), units = 1:50)
        expect_equal(unname(a$overlap$fuzzy), unname(b$overlap$fuzzy))
        expect_equal(a$units$support, b$units$support)
        expect_true(b$settings$external_fuzzy)
})

test_that("fuzzy and crisp overlap carry the total weight; ARI ignores weights", {
        .local_seed(5)
        U   <- matrix(stats::runif(120), 40); U <- U / rowSums(U)
        ext <- cw_ext(1:40, sample(c("a", "b"), 40, TRUE), area = stats::runif(40, 1, 100))
        cw  <- crosswalk_types(cw_regions(U), ext, "type", units = 1:40, weights = "area")
        un  <- crosswalk_types(cw_regions(U), ext, "type", units = 1:40)
        expect_equal(sum(cw$overlap$fuzzy), sum(ext$area))
        expect_equal(sum(cw$overlap$crisp), sum(ext$area))
        expect_equal(sum(cw$overlap$counts), 40)
        ari <- function(z) z$agreement$crisp[z$agreement$index == "ARI"]
        expect_equal(ari(cw), ari(un))
})

test_that("mosaic units without a type are excluded and counted", {
        valid <- rep(c(TRUE, FALSE), c(30, 5))
        h  <- rep(1:3, each = 10)
        cw <- crosswalk_types(cw_mosaic(cw_onehot(h), valid), cw_ext(1:35, rep("a", 35)),
                              "type", units = 1:35)
        expect_equal(cw$coverage$n_matched, 30)
        expect_equal(cw$coverage$n_no_pulse_type, 5)
        expect_setequal(cw$units$id, as.character(1:30))
})

test_that("within_region equals separate crosswalks per region", {
        .local_seed(6)
        reg_of_type <- c(1, 1, 2, 2, 2)
        h <- c(sample(1:2, 40, TRUE), sample(3:5, 60, TRUE))
        e <- sample(c("x", "y", "z"), 100, TRUE)
        cw <- crosswalk_types(cw_local(cw_onehot(h, 5), reg_of_type), cw_ext(1:100, e),
                              "type", units = 1:100, local_scope = "within_region")
        br <- cw$by_region[cw$by_region$index == "ARI", ]
        expect_equal(br$crisp[br$region == 1], .adjusted_rand_index(h[1:40], e[1:40]))
        expect_equal(br$crisp[br$region == 2], .adjusted_rand_index(h[41:100], e[41:100]))
        expect_equal(cw$agreement$crisp[cw$agreement$index == "ARI"],
                     sum(br$crisp * c(40, 60)) / 100)
        expect_equal(cw$scope, "within_region")
})

test_that("bad inputs fail with clear messages", {
        obj <- cw_regions(cw_onehot(rep(1:2, 5)))
        expect_error(crosswalk_types(obj, cw_ext(c(1, 1:9), rep("a", 10)), "type", units = 1:10),
                     "duplicated ID")
        expect_error(crosswalk_types(obj, cw_ext(1:10, rep("a", 10)), "type"),
                     "`units` is required")
        expect_error(crosswalk_types(obj, cw_ext(1:10, rep("a", 10)), "type", units = 1:9),
                     "has 9 IDs")
        expect_error(crosswalk_types(obj, cw_ext(1:10, rep("a", 10)), "nope", units = 1:10),
                     "not found")
        expect_warning(crosswalk_types(obj, cw_ext(1:9, rep("a", 9)), "type", units = 1:10),
                       "Only 90.0%")
})

test_that("a named list crosswalks every level", {
        h1 <- rep(1:2, each = 30); h2 <- rep(1:6, each = 10)
        set <- crosswalk_types(list(regions = cw_regions(cw_onehot(h1)),
                                    local = cw_local(cw_onehot(h2), rep(1:2, each = 3))),
                               cw_ext(1:60, letters[h1]), "type", units = 1:60)
        expect_s3_class(set, "pulse_crosswalk_set")
        expect_equal(set$by_level$level, c("regions", "local"))
        expect_equal(set$by_level$AMI[1], 1)
        expect_output(print(set), "Compare levels on AMI")
        expect_output(print(set$results$local), "global scope includes the level 1 signal")
})

test_that("sf external and igraph units work together", {
        g   <- make_grid_sf(4L, 4L)
        net <- igraph::set_vertex_attr(igraph::make_ring(16), "polygon_id", value = rev(g$ID))
        g$type <- ifelse(g$temp > stats::median(g$temp), "warm", "cold")
        h  <- ifelse(rev(g$type) == "warm", 1L, 2L)        # vertex order is reversed
        cw <- crosswalk_types(cw_regions(cw_onehot(h)), g, "type", units = net)
        expect_equal(cw$agreement$crisp[cw$agreement$index == "ARI"], 1)
})
