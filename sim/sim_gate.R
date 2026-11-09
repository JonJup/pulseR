# =========================================================================== *
# sim_gate.R
#
# Fast go/no-go gate for iterative development. Local desktop, any OS
# (PSOCK workers, so Windows works). Target: under 20 min on 4 cores,
# hard cap `budget_min`.
#
#   source("sim_gate.R")
#   v <- run_gate(reps = 1, cores = 4, pkg_path = "~/pulseR")   # pilot
#   v <- run_gate(reps = 6, cores = 4, pkg_path = "~/pulseR")   # gate
#
# Verdicts per rule:
#   PASS    claim holds by more than 2 MCSE
#   STOP    claim fails by more than 2 MCSE, or no data (stage failed)
#   WARN    same, for rules that should not block development
#   UNSURE  within 2 MCSE of the threshold: add reps before deciding
# Overall: STOP if any STOP, CHECK if any UNSURE or WARN, else GO.
#
# Every call writes to a NEW time-stamped out_dir by default. run_one()
# skips existing files, so reusing an out_dir after a code change would
# silently report the OLD code. Pass the same out_dir only to resume.
#
# Fixed master_seed = same 30 landscapes every time, so verdicts are
# comparable across code versions. Rotate it before trusting a GO: tuning
# against fixed landscapes is teaching to the test.
# =========================================================================== *

.gate_files <- c("sim_dgp.R", "sim_metrics.R", "sim_pipeline.R", "sim_run.R")
for (f in .gate_files) source(f)

# --------------------------------------------------------------------------- *
# Design: 5 cells, L1 dose-response plus null, one L2 cell
# --------------------------------------------------------------------------- *

gate_design <- function(cfg = list()) {
    fast <- utils::modifyList(list(N_tree = 100L, regions_per_partition = 20L,
                                   stability_B = 2L, n_local_types = 1:4,
                                   baselines = TRUE), cfg)
    base <- list(N = 800L, R_true = 4L, K_types = 2L, share_types = 1)
    l1   <- c(fast, list(levels = 1L))
    d <- list(
        cell("G", "null", dgp = c(base, list(contrast = 0)),
             cfg = c(l1, list(floor = TRUE))),
        cell("G", "c1",   dgp = c(base, list(contrast = 1)), cfg = l1),
        cell("G", "c2",   dgp = c(base, list(contrast = 2)), cfg = l1),
        cell("G", "c4",   dgp = c(base, list(contrast = 4)), cfg = l1),
        cell("G", "l2",   dgp = c(base, list(contrast = 2.5, share_types = 0.5)),
             cfg = c(fast, list(levels = 1:2, oracle = TRUE))))
    stats::setNames(d, vapply(d, `[[`, "", "id"))
}

#' Thresholds. Fix these BEFORE looking at results.
gate_thresholds <- function() list(
    recover   = 0.80,   # c4: pulseR ARI / Bayes ARI >= 0.80
    vs_base   = -0.05,  # c2+c4: pulseR ARI - baseline ARI >= -0.05 (paired)
    dose      = 0.20,   # ARI(c4) - ARI(c1) >= 0.20
    null_ari  = 0.10,   # null: pulseR ARI - permuted-floor ARI <= 0.10
    conf_drop = 0.05,   # mean certainty(c4) - mean certainty(null) >= 0.05
    auc       = 0.65,   # c2: AUC of (1 - certainty) for errors >= 0.65
    handoff   = 0.15)   # l2: oracle L2 ARI - estimated L2 ARI <= 0.15

# --------------------------------------------------------------------------- *
# Runner
# --------------------------------------------------------------------------- *

run_gate <- function(reps = 6L, cores = 4L, budget_min = 45,
                     out_dir = file.path("gate_out", format(Sys.time(), "%Y%m%d-%H%M%S")),
                     pkg_path = NULL, master_seed = 20261105L,
                     thr = gate_thresholds(), cfg = list(), opts = list()) {
    gate_preflight(pkg_path, cores)
    if (!is.null(pkg_path)) pkgload::load_all(pkg_path, quiet = TRUE)
    d  <- gate_design(cfg)
    old_opts <- options(opts); on.exit(options(old_opts), add = TRUE)
    jb <- job_list(d, reps)                       # rep-major: a cut-off leaves
    saveRDS(d, file.path(.mkdir(out_dir), "design.rds"))   # balanced reps
    deadline <- Sys.time() + 60 * budget_min      # no NEW job starts after this
    t0 <- Sys.time()
    status <- if (cores > 1L) {
        cl <- parallel::makePSOCKcluster(cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterCall(cl, .gate_init, getwd(), .gate_files, pkg_path, opts)
        unlist(parallel::parLapplyLB(cl, seq_len(nrow(jb)), .gate_job, d = d,
                                     jb = jb, out_dir = out_dir,
                                     master_seed = master_seed, deadline = deadline))
    } else {
        vapply(seq_len(nrow(jb)), .gate_job, "", d = d, jb = jb, out_dir = out_dir,
               master_seed = master_seed, deadline = deadline)
    }
    res <- collect_results(out_dir)
    v   <- gate_verdict(res, thr)
    v$run <- list(out_dir = out_dir, reps = reps, cores = cores,
                  wall_min = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  status = table(sub(":.*", "", status)),
                  proj_min_6reps = sum(v$timing$per_run_s) * 6 / cores / 60)
    saveRDS(v, file.path(out_dir, "verdict.rds"))
    print_gate(v)
    invisible(v)
}

gate_preflight <- function(pkg_path, cores) {
    need <- c("sf", "spdep", "igraph", "RANN", "ClustGeo",
              if (is.null(pkg_path)) "pulseR" else "pkgload")
    miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss))
        stop("Missing packages. Run:\n  install.packages(c(",
             paste0("\"", miss, "\"", collapse = ", "), "))", call. = FALSE)
    if (is.null(pkg_path) && cores > 1L && "pulseR" %in% loadedNamespaces() &&
        requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("pulseR"))
        stop("pulseR is load_all()-ed in this session, but workers would load the ",
             "INSTALLED version. Pass pkg_path = \"<path to pulseR source>\".",
             call. = FALSE)
    invisible(TRUE)
}

.gate_init <- function(wd, files, pkg_path, opts = list()) {
    setwd(wd)
    if (!is.null(pkg_path)) pkgload::load_all(pkg_path, quiet = TRUE)
    options(opts)
    for (f in files) source(f)
    NULL
}

.gate_job <- function(i, d, jb, out_dir, master_seed, deadline) {
    if (Sys.time() > deadline) return("skipped")
    tryCatch({ run_one(d[[jb$cell[i]]], jb$rep[i], out_dir, master_seed); "done" },
             error = function(e) paste0("error: ", conditionMessage(e)))
}

# --------------------------------------------------------------------------- *
# Verdict
# --------------------------------------------------------------------------- *

gate_verdict <- function(res, thr = gate_thresholds()) {
    M <- res$metrics
    if (is.null(M)) stop("No finished runs in out_dir.", call. = FALSE)
    v <- function(cell, col) {
        r <- M[M$cell_id == paste0("G_", cell), , drop = FALSE]
        if (col %in% names(r)) as.numeric(r[[col]]) else numeric()
    }
    pr <- function(x) {                           # one sample (or paired diffs)
        x <- x[is.finite(x)]; n <- length(x)
        c(est = if (n) mean(x) else NA, se = if (n > 1L) stats::sd(x) / sqrt(n) else NA, n = n)
    }
    up <- function(a, b) {                        # two independent samples
        a <- a[is.finite(a)]; b <- b[is.finite(b)]
        if (!length(a) || !length(b)) return(c(est = NA, se = NA, n = 0))
        se <- sqrt(stats::var(a) / length(a) + stats::var(b) / length(b))
        c(est = mean(a) - mean(b), se = se, n = min(length(a), length(b)))
    }
    paired2 <- function(col_a, col_b)             # pooled over c2 and c4
        pr(c(v("c2", col_a) - v("c2", col_b), v("c4", col_a) - v("c4", col_b)))

    st  <- res$stages
    bad <- if (is.null(st)) NULL else st[!st$ok, c("cell_id", "rep", "stage", "error")]
    n_bad <- if (is.null(st)) NA else nrow(bad)
    pieces <- M$l1_excess_pieces[is.finite(M$l1_excess_pieces)]

    rules <- rbind(
        .judge("health",      "all",   c(est = n_bad, se = 0, n = nrow(M)), 0, "<=", "STOP",
               "every stage of every run succeeded"),
        .judge("contiguous",  "all",   c(est = if (length(pieces)) max(pieces) else NA,
                                         se = 0, n = length(pieces)), 0, "<=", "STOP",
               "consensus regions are contiguous"),
        .judge("recover",     "c4",    pr(v("c4", "l1_ari") / v("c4", "l1_ari_bayes")),
               thr$recover, ">=", "STOP", "pulseR ARI near the Bayes ceiling at strong signal"),
        .judge("vs_skater",   "c2+c4", paired2("l1_ari", "skater_ari"),
               thr$vs_base, ">=", "STOP", "ensemble not worse than its single-tree ablation"),
        .judge("vs_clustgeo", "c2+c4", paired2("l1_ari", "clustgeo_ari"),
               thr$vs_base, ">=", "WARN", "not worse than ClustGeo"),
        .judge("dose",        "c4-c1", up(v("c4", "l1_ari"), v("c1", "l1_ari")),
               thr$dose, ">=", "STOP", "ARI rises with contrast"),
        .judge("null_ari",    "null",  pr(v("null", "l1_ari") - v("null", "ari_l1_floor")),
               thr$null_ari, "<=", "STOP", "no structure beyond the permuted floor"),
        .judge("conf_drop",   "c4-null", up(v("c4", "l1_mean_cert"), v("null", "l1_mean_cert")),
               thr$conf_drop, ">=", "STOP", "certainty falls when signal is absent"),
        .judge("auc",         "c2",    pr(v("c2", "l1_auc")),
               thr$auc, ">=", "STOP", "low certainty flags errors"),
        .judge("handoff",     "l2",    pr(v("l2", "l2_handoff")),
               thr$handoff, "<=", "WARN", "L1 errors cost L2 little"))

    overall <- if (any(rules$verdict == "STOP")) "STOP" else
        if (any(rules$verdict %in% c("UNSURE", "WARN"))) "CHECK" else "GO"

    per_run  <- tapply(as.numeric(M$seconds), M$cell_id, stats::median)
    per_stage <- if (is.null(st)) NULL else tapply(st$seconds, st$stage, stats::median)
    list(overall = overall, rules = rules, failed_stages = bad,
         baselines = .baseline_table(M),
         timing = list(per_run_s = per_run, per_stage_s = per_stage))
}

.judge <- function(id, cell, s, thr, dir, sev, claim) {
    est <- unname(s[["est"]]); se <- unname(s[["se"]])
    gap <- if (dir == ">=") est - thr else thr - est
    verdict <- if (!is.finite(gap)) "STOP" else if (!is.finite(se)) "UNSURE" else
        if (gap >= 2 * se) "PASS" else if (gap < -2 * se) sev else "UNSURE"
    data.frame(verdict = verdict, rule = id, cell = cell, est = round(est, 3),
               se = round(se, 3), n = unname(s[["n"]]), need = paste(dir, thr),
               claim = claim)
}

# Mean ARI per cell for pulseR, both baselines and the Bayes reference.
.baseline_table <- function(M) {
    cols <- intersect(c("l1_ari", "skater_ari", "clustgeo_ari", "l1_ari_bayes",
                        "clustgeo_alpha"), names(M))
    if (!length(cols)) return(NULL)
    stats::aggregate(M[cols], by = list(cell = M$cell_id),
                     FUN = function(x) round(mean(as.numeric(x), na.rm = TRUE), 3))
}

print_gate <- function(v) {
    cat(sprintf("\nGATE: %s   (%s, %.1f min wall, %d cores)\n\n", v$overall,
                paste(names(v$run$status), v$run$status, collapse = " / "),
                v$run$wall_min, v$run$cores))
    print(v$rules[, c("verdict", "rule", "cell", "est", "se", "n", "need")], row.names = FALSE)
    cat("\nMean ARI by cell:\n"); print(v$baselines, row.names = FALSE)
    if (!is.null(v$failed_stages) && nrow(v$failed_stages)) {
        cat("\nFailed stages:\n"); print(utils::head(v$failed_stages, 10), row.names = FALSE)
    }
    cat("\nMedian seconds per stage:\n"); print(round(v$timing$per_stage_s, 1))
    cat(sprintf("\nProjected wall time for reps = 6 on %d cores: %.0f min\n",
                v$run$cores, v$run$proj_min_6reps))
    cat("Output:", v$run$out_dir, "\n")
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
