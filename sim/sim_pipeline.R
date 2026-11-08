# =========================================================================== *
# sim_pipeline.R
#
# Adapter between simulate_landscape() and the pulseR API, plus scoring.
# This is the ONLY file that calls pulseR. When the API changes, fix it here.
#
# Stages (each timed; errors and warnings are recorded as data, not thrown):
#   graph -> weights -> regions -> [floor] -> local (est) -> local (oracle)
#   -> meta -> mosaic (est) -> mosaic (oracle)
#
# Oracles isolate each level from upstream error:
#   local oracle  : get_local_types() on one-hot TRUE regions
#   mosaic oracle : signatures from one-hot TRUE meta-types
# Handoff loss = oracle score - estimated score.
# =========================================================================== *

#' Pipeline settings. Anything not given takes these values.
default_cfg <- function() list(
    weights_on            = "all",     # "all" columns or "regional" only for L1
    final_regions         = NULL,      # NULL = R_true (fixed k); vector = tuned
    regions_per_partition = 50L,
    N_tree                = 200L,
    weight_exponent       = 1,
    size_influence        = 2,
    target_edge_turnover  = 0.2,
    tuning_strategy       = "grid",
    stability_B           = 10L,
    large_n_threshold     = 20000,
    n_local_types         = 1:6,
    local_cutoff          = 0.1,
    meta_method           = "bhattacharyya",
    k_order               = 2L,
    mosaic_engine         = "dmm",
    n_mosaic              = 2:8,
    mosaic_sample_size    = 5000L,
    floor                 = FALSE,     # contiguity-only floor (one extra L1 fit)
    null_B                = 0L,        # >0: compute_null_baseline() replicates
    oracle                = TRUE,
    baselines = TRUE,
    refine = FALSE,
    cg_alpha = "cross",
    levels                = 1:3)

#' Run the pipeline on one simulated landscape and score it.
#' @return list(metrics = 1-row data.frame, units = data.frame (one row per
#'   unit), stages = data.frame, reliability = data.frame)
run_pipeline <- function(sim, cfg = list(), seed = 1L) {
    cfg <- utils::modifyList(default_cfg(), cfg)
    tr  <- sim$truth
    N   <- length(tr$z)
    log <- new.env()
    m   <- list()                                 # scalar metrics
    u   <- data.frame(unit = seq_len(N), z = tr$z, d_boundary = tr$d_boundary,
                      ustar_max = apply(tr$ustar, 1, max), region_size = tr$region_size,
                      zone = tr$zone, arrangement = tr$arrangement,
                      shared = tr$shared)

    vars_all <- unlist(sim$env_cols, use.names = FALSE)
    vars_l1  <- if (cfg$weights_on == "regional") sim$env_cols$reg else vars_all
    env_df   <- sf::st_drop_geometry(sim$polygons)

    ## ---- Graph and weights -------------------------------------------------*
    g  <- .stage(log, "graph", pulseR::polygon_to_network(sim$polygons, verbose = FALSE))
    gw <- .stage(log, "weights", if (.ok(g))
        pulseR::add_edge_weight(g, variables = vars_l1, scale = TRUE, verbose = FALSE))
    if (!.ok(gw)) return(.finish(m, u, log, NULL))
    m$n_bridges <- igraph::ecount(g$graph) - nrow(g$edge_data)

    ## ---- L1 -----------------------------------------------------------------*
    k_req <- if (is.null(cfg$final_regions)) tr$R_true else cfg$final_regions
    fit_l1 <- function(graph) pulseR::get_regions(
        graph, final_regions = k_req,
        regions_per_partition = cfg$regions_per_partition, N_tree = cfg$N_tree,
        weight_exponent = cfg$weight_exponent, size_influence = cfg$size_influence,
        target_edge_turnover = cfg$target_edge_turnover,
        tuning_strategy = cfg$tuning_strategy, stability_B = cfg$stability_B,
        large_n_threshold = cfg$large_n_threshold, verbose = FALSE, seed = seed)
    reg <- .stage(log, "regions", fit_l1(gw))
    if (.ok(reg)) {
        sc <- score_l1(reg, tr)
        m  <- c(m, sc$metrics); u <- cbind(u, sc$units)
        rel_l1 <- sc$reliability
    } else rel_l1 <- NULL

    if (.ok(reg)) m$l1_excess_pieces <- excess_pieces(reg$hard_clusters, gw)
    
    if (.ok(reg) && isTRUE(cfg$refine)) {
        rf <- .stage(log, "refine", .refine_regions(reg$hard_clusters, gw,
                                                    igraph::graph_attr(gw, "env"),
                                                    min_size = getOption("pulseR.min_size", 1L)))
        if (.ok(rf)) {
            m$l1_ari_refined    <- ari(rf$z, tr$z)
            m$l1_refined_pieces <- excess_pieces(rf$z, gw)
            m$l1_refine_moves   <- rf$moves
        }
    }
    
    if (isTRUE(cfg$baselines)) m <- c(m, .run_baselines(log, gw, tr, cfg$cg_alpha))
    
    if (isTRUE(cfg$floor)) {
        gp <- .permute_graph(gw)
        rf <- .stage(log, "floor", fit_l1(gp))
        m$ari_l1_floor <- if (.ok(rf)) ari(rf$hard_clusters, tr$z) else NA_real_
    }
    if (cfg$null_B > 0L && .ok(reg) && !reg$tuned) {
        nl <- .stage(log, "null", .null_test(reg, gw, cfg, seed))
        if (.ok(nl)) m <- c(m, nl)
    }
    if (!(2L %in% cfg$levels)) return(.finish(m, u, log, rel_l1))

    ## ---- L2 -----------------------------------------------------------------*
    lt <- .stage(log, "local", if (.ok(reg)) pulseR::get_local_types(
        env_df, reg, variables = vars_all, n_local_types = cfg$n_local_types,
        cutoff = cfg$local_cutoff))
    if (.ok(lt)) {
        sc <- score_l2(lt, tr, "l2"); m <- c(m, sc$metrics); u <- cbind(u, sc$units)
    }
    if (isTRUE(cfg$oracle)) {
        lt_o <- .stage(log, "local_oracle", pulseR::get_local_types(
            env_df, onehot(tr$z), variables = vars_all,
            n_local_types = cfg$n_local_types, cutoff = cfg$local_cutoff))
        if (.ok(lt_o)) {
            sc <- score_l2(lt_o, tr, "l2o"); m <- c(m, sc$metrics)
            G  <- lt_o$n_local_types
            m$l2o_G_exact <- mean(G == tr$K_r)
            m$l2o_G_err   <- mean(G - tr$K_r)
            m$l2_handoff  <- m$l2o_ari - .na(m$l2_ari)
        }
    }
    m$l2_ari_bayes <- ari(paste(tr$z, tr$bayes_l2), tr$joint)

    ## ---- L2b meta-types -------------------------------------------------------*
    meta <- .stage(log, "meta", if (.ok(lt))
        pulseR::get_meta_types(lt, method = cfg$meta_method))
    if (.ok(meta)) {
        m$meta_ari       <- ari(meta$hard_classes, tr$type_id)
        m$meta_n         <- meta$n_meta
        m$meta_n_err     <- meta$n_meta - tr$n_meta_true
        m$meta_ari_reg   <- meta$crossing$ari_with_region
        m$meta_frac_single <- meta$crossing$frac_single_region
        mp <- fold_map(meta$hard_classes, tr$type_id)
        u$meta_correct <- mp[as.character(meta$hard_classes)] == tr$type_id
    }
    if (!(3L %in% cfg$levels)) return(.finish(m, u, log, rel_l1))

    ## ---- L3 -----------------------------------------------------------------*
    run_mosaic <- function(P) {
        sg <- pulseR::compute_all_signatures(pulseR::as_mosaic_P(P), gw,
                                             k = cfg$k_order, verbose = FALSE)
        list(sigs = sg, fit = pulseR::get_mosaic_types(
            sg, n_clusters = cfg$n_mosaic, method = cfg$mosaic_engine,
            sample_size = cfg$mosaic_sample_size, seed = seed, verbose = FALSE))
    }
    mo <- .stage(log, "mosaic", if (.ok(meta)) run_mosaic(meta))
    if (.ok(mo)) {
        sc <- score_l3(mo$fit, tr, "l3"); m <- c(m, sc$metrics); u$l3_hard <- sc$hard
    }
    if (isTRUE(cfg$oracle)) {
        mo_o <- .stage(log, "mosaic_oracle", run_mosaic(onehot(tr$type_id)))
        if (.ok(mo_o)) {
            sc <- score_l3(mo_o$fit, tr, "l3o"); m <- c(m, sc$metrics)
            m$l3_handoff <- m$l3o_ari - .na(m$l3_ari)
        }
    }
    .finish(m, u, log, rel_l1)
}

# --------------------------------------------------------------------------- *
# Scoring
# --------------------------------------------------------------------------- *

score_l1 <- function(reg, tr) {
    hard <- reg$hard_clusters
    Pi   <- reg$memberships
    R    <- tr$R_true
    mp   <- fold_map(hard, tr$z)
    F    <- fold_columns(Pi, mp, seq_len(ncol(Pi)), as.character(seq_len(R)))
    F[!stats::complete.cases(Pi), ] <- NA
    wrong <- as.integer(mp[as.character(hard)]) != tr$z
    cert  <- apply(Pi, 1, function(r) if (anyNA(r)) NA_real_ else max(r))
    coh   <- reg$fit$cohesion
    coh_c <- if (is.null(coh)) rep(NA_real_, length(hard)) else
        apply(coh / rowSums(coh), 1, max)
    interior <- tr$ustar[cbind(seq_along(hard), tr$z)] > 0.99
    vm   <- v_measure(tr$z, hard)
    Fc   <- F[stats::complete.cases(F), , drop = FALSE]
    zc   <- tr$z[stats::complete.cases(F)]

    metrics <- list(
        l1_k_hat        = reg$k,
        l1_k_err        = reg$k - R,
        l1_ari          = ari(hard, tr$z),
        l1_homogeneity  = vm[["homogeneity"]],
        l1_completeness = vm[["completeness"]],
        l1_vmeasure     = vm[["v_measure"]],
        l1_ari_bayes    = ari(max.col(tr$bayes_l1, "first"), tr$z),
        l1_err_rate     = mean(wrong),
        l1_auc          = auc(1 - cert, wrong),
        l1_auc_cohesion = auc(1 - coh_c, wrong),
        l1_brier        = brier(Fc, zc),
        l1_ece          = ece(Fc, zc),
        l1_mean_cert    = mean(cert, na.rm = TRUE),
        l1_mae_ustar    = mean(abs(F - tr$ustar), na.rm = TRUE),
        l1_size_cor_pi  = suppressWarnings(stats::cor(cert[interior], log(tr$region_size[interior]),
                                                      method = "spearman", use = "complete.obs")),
        l1_size_cor_coh = suppressWarnings(stats::cor(coh_c[interior], log(tr$region_size[interior]),
                                                      method = "spearman", use = "complete.obs")),
        l1_frac_na_pi   = mean(is.na(cert)),
        l1_singleton    = singleton_share(reg$ensemble),
        l1_n_eff        = reg$n_eff,
        l1_tuned        = reg$tuned)
    units <- data.frame(l1_hard = hard, l1_wrong = wrong, l1_cert = cert,
                        l1_cert_coh = coh_c)
    rel <- reliability_bins(apply(Fc, 1, max), max.col(Fc, "first") == zc)
    list(metrics = metrics, units = units, reliability = rel)
}

score_l2 <- function(lt, tr, prefix) {
    hard  <- lt$hard_classes
    mp    <- fold_map(hard, tr$joint)
    wrong <- mp[as.character(hard)] != tr$joint
    metrics <- list(
        ari        = ari(hard, tr$joint),
        ari_within = ari_within(hard, tr$local, tr$z),
        auc        = auc(1 - lt$certainty, wrong),
        err_rate   = mean(wrong),
        n_types    = ncol(lt$memberships),
        mean_cert  = mean(lt$certainty))
    names(metrics) <- paste0(prefix, "_", names(metrics))
    units <- stats::setNames(data.frame(hard, wrong, lt$certainty),
                             paste0(prefix, c("_hard", "_wrong", "_cert")))
    list(metrics = metrics, units = units)
}

score_l3 <- function(fit, tr, prefix) {
    hard <- rep(NA_integer_, length(tr$zone))
    hard[fit$valid] <- fit$hard_clusters
    v    <- fit$validity
    kcol <- intersect(c("k", "n_clusters", "K"), names(v))[1]
    row  <- if (is.na(kcol)) NULL else v[v[[kcol]] == fit$best_k, , drop = FALSE]
    get  <- function(col) if (!is.null(row) && col %in% names(row)) row[[col]][1] else NA_real_
    metrics <- list(
        ari         = ari(hard, tr$zone),
        ari_arrange = ari(hard, tr$arrangement),
        k_hat       = fit$best_k,
        k_err       = fit$best_k - tr$M_true,
        R2_JS       = get("R2_JS"),
        saturation  = if (is.numeric(fit$saturation)) fit$saturation[1] else NA_real_,
        frac_valid  = mean(fit$valid))
    names(metrics) <- paste0(prefix, "_", names(metrics))
    list(metrics = metrics, hard = hard)
}

# Baselines at k = R_true on the L1 graph, weights and scaled variables.
# SKATER: one MST, cut by spdep::skater (the ensemble's ablation).
# ClustGeo: Ward on D0 = env distance, D1 = centroid distance; alpha where the
#   normalised feature and spatial homogeneity curves cross. Crisp only.
.run_baselines <- function(log, gw, tr, cg_alpha = "cross") {
    k <- tr$R_true
    X <- igraph::graph_attr(gw, "env")
    m <- list()
    sk <- .stage(log, "skater", {
        tree <- igraph::mst(gw, weights = igraph::E(gw)$weight)
        spdep::skater(edges = igraph::ends(tree, igraph::E(tree), names = FALSE),
                      data = X, ncuts = k - 1L)$groups
    })
    if (.ok(sk)) m <- c(m, score_crisp(sk, tr, gw, "skater"))
    cg <- .stage(log, "clustgeo", {
        D0 <- stats::dist(X); D1 <- stats::dist(tr$xy)
        al <- if (is.numeric(cg_alpha)) cg_alpha else {
            ca <- ClustGeo::choicealpha(D0, D1, range.alpha = seq(0, 1, 0.1),
                                        K = k, graph = FALSE)
            ca$range.alpha[which.min(abs(ca$Qnorm[, 1] - ca$Qnorm[, 2]))]
        }
        list(hard = stats::cutree(ClustGeo::hclustgeo(D0, D1, alpha = al), k = k),
             alpha = al)
    })
    if (.ok(cg)) m <- c(m, score_crisp(cg$hard, tr, gw, "clustgeo"),
                        list(clustgeo_alpha = cg$alpha))
    m
}

score_crisp <- function(hard, tr, g, prefix) {
    vm <- v_measure(tr$z, hard)
    m  <- list(ari = ari(hard, tr$z), homogeneity = vm[["homogeneity"]],
               completeness = vm[["completeness"]],
               excess_pieces = excess_pieces(hard, g))
    stats::setNames(m, paste0(prefix, "_", names(m)))
}
# --------------------------------------------------------------------------- *
# Helpers
# --------------------------------------------------------------------------- *

# Evaluate expr, record time, warnings and errors under `name` in `log`.
.stage <- function(log, name, expr) {
    warns <- character()
    t0 <- proc.time()[["elapsed"]]
    val <- withCallingHandlers(
        tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)),
                                                     class = "stage_error")),
        warning = function(w) { warns <<- c(warns, conditionMessage(w))
                                invokeRestart("muffleWarning") },
        message = function(msg) invokeRestart("muffleMessage"))
    log[[name]] <- data.frame(
        stage = name, seconds = proc.time()[["elapsed"]] - t0,
        ok = !inherits(val, "stage_error") && !is.null(val),
        error = if (inherits(val, "stage_error")) val$msg else NA_character_,
        n_warn = length(warns),
        warnings = if (length(warns)) paste(unique(warns), collapse = " | ") else NA_character_)
    val
}

.ok <- function(x) !is.null(x) && !inherits(x, "stage_error")
.na <- function(x) if (is.null(x)) NA_real_ else x

.finish <- function(m, u, log, rel) {
    st <- do.call(rbind, mget(ls(log), envir = log))
    rownames(st) <- NULL
    m  <- lapply(m, function(v) if (is.null(v) || length(v) != 1L) NA else v)
    list(metrics = as.data.frame(m), units = u, stages = st, reliability = rel)
}

# Permutation null in the edge-weight space (B2-safe: env from the graph).
# Observed CVIs use the package's own scorer, so null and observed match.
.null_test <- function(reg, gw, cfg, seed) {
    env <- igraph::graph_attr(gw, "env")
    obs <- pulseR:::.consensus_cvis(reg$fit, reg$ensemble, env = env,
                                    stability_B = cfg$stability_B, stability_seed = seed,
                                    large_n_threshold = cfg$large_n_threshold)
    obs$final_regions <- reg$k
    nb <- pulseR::compute_null_baseline(
        gw, final_regions = reg$k, regions_per_partition = cfg$regions_per_partition,
        N_tree = cfg$N_tree, weight_exponent = cfg$weight_exponent,
        size_influence = cfg$size_influence, target_edge_turnover = cfg$target_edge_turnover,
        B = cfg$null_B, stability_B = cfg$stability_B, observed = obs,
        large_n_threshold = cfg$large_n_threshold, seed = seed, verbose = FALSE)
    s <- nb$summary
    list(null_p_ch   = s$p[s$criterion == "CH_ENV"],
         null_p_stab = s$p[s$criterion == "STAB"],
         null_excess_ch   = s$excess[s$criterion == "CH_ENV"],
         null_excess_stab = s$excess[s$criterion == "STAB"])
}

# Same graph, edge weights recomputed from row-permuted environment: keeps
# contiguity, destroys environmental signal.
.permute_graph <- function(gw) {
    X  <- igraph::graph_attr(gw, "env")
    Xp <- X[sample.int(nrow(X)), , drop = FALSE]
    ep <- igraph::ends(gw, igraph::E(gw), names = FALSE)
    igraph::E(gw)$weight <- sqrt(rowSums((Xp[ep[, 1], , drop = FALSE] -
                                              Xp[ep[, 2], , drop = FALSE])^2))
    rownames(Xp) <- rownames(X)
    igraph::set_graph_attr(gw, "env", Xp)
}

# Contiguity-preserving SSD descent on a hard partition. Moves a boundary
# unit to an adjacent region if total within-region SSD drops and the old
# region stays connected. Same objective as the "ssd" cut rule.
.refine_regions <- function(z, graph, X, max_pass = 20L) {
    z   <- match(z, sort(unique(z))); k <- max(z)
    nb  <- lapply(igraph::adjacent_vertices(graph, igraph::V(graph)), as.integer)
    n_r <- tabulate(z, k)
    S   <- rowsum(X, z, reorder = TRUE)
    moves <- 0L
    for (pass in seq_len(max_pass)) {
        moved <- 0L
        bnd <- which(vapply(seq_along(z), function(i) any(z[nb[[i]]] != z[i]), logical(1)))
        for (i in bnd) {
            a <- z[i]; cand <- setdiff(unique(z[nb[[i]]]), a)
            if (!length(cand) || n_r[a] <= 1L) next
            out  <- n_r[a] / (n_r[a] - 1) * sum((X[i, ] - S[a, ] / n_r[a])^2)
            gain <- vapply(cand, function(b)
                out - n_r[b] / (n_r[b] + 1) * sum((X[i, ] - S[b, ] / n_r[b])^2), numeric(1))
            if (max(gain) <= 0) next
            b    <- cand[which.max(gain)]
            same <- nb[[i]][z[nb[[i]]] == a]
            if (length(same) > 1L) {              # is i a cut vertex of region a?
                ia <- setdiff(which(z == a), i)
                cm <- igraph::components(igraph::induced_subgraph(graph, ia))$membership
                if (length(unique(cm[match(same, ia)])) > 1L) next
            }
            z[i] <- b
            S[a, ] <- S[a, ] - X[i, ]; S[b, ] <- S[b, ] + X[i, ]
            n_r[a] <- n_r[a] - 1L;     n_r[b] <- n_r[b] + 1L
            moved <- moved + 1L
        }
        moves <- moves + moved
        if (!moved) break
    }
    list(z = z, moves = moves, passes = pass)
}