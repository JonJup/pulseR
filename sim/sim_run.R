# =========================================================================== *
# sim_run.R
#
# Design, runner and collector.
#
#   design <- make_design()                       # M1-M5 cells, Tier I
#   run_design(design, reps = 20, out_dir = "out", cores = 8)
#   res <- collect_results("out")                 # metrics, stages, units
#
# On Slurm (one array task per chunk):
#   run_chunk(design, reps, chunk = task_id, n_chunks = n_tasks, out_dir)
#
# Seeds. The landscape seed depends on (dgp settings, rep) only, so cells that
# differ only in pipeline settings see the SAME landscapes (common random
# numbers: knob effects are paired, not confounded with DGP noise). The
# pipeline seed depends on (cell, rep).
#
# Resume. One .rds per (cell, rep). Existing files are skipped, so a killed
# job is restarted by running the same call again.
# =========================================================================== *

# --------------------------------------------------------------------------- *
# Design
# --------------------------------------------------------------------------- *

core_dgp <- function() list(N = 2000L, R_true = 6L, contrast = 2.5, tau = 1,
                            K_types = 3L, share_types = 0.5, delta = 2,
                            zones = c("block", "fine", "banded", "random"))

#' One design cell.
#' @param module "M1".."M10".
#' @param label Short readable id; must be unique within the design.
cell <- function(module, label, dgp = list(), cfg = list()) {
    list(id = paste0(module, "_", label), module = module, label = label,
         dgp = utils::modifyList(core_dgp(), dgp), cfg = cfg)
}

#' Tier I design for the core-claim modules M1-M5.
#' Core scenario plus one-factor-at-a-time departures; crossings only where
#' the interaction is the claim (contrast x tau, contrast x handoff).
make_design <- function() {
    d <- list()
    add <- function(x) d[[x$id]] <<- x

    ## M1: L1 recovery and null. contrast x N x R_true, L1 only.
    ## share_types = 1: confined local types carry region signal, so with
    ## share_types < 1 contrast = 0 is NOT a null (smoke test: ARI 0.73).
    for (ct in c(0, 0.5, 1, 1.5, 2, 3, 4)) for (n in c(500L, 2000L)) for (R in c(4L, 8L))
        add(cell("M1", sprintf("c%.1f_N%d_R%d", ct, n, R),
                 dgp = list(contrast = ct, N = n, R_true = R, share_types = 1),
                 cfg = list(levels = 1L, floor = TRUE, null_B = 19L)))

    ## M2: Pi as uncertainty. tau x contrast, L1 only.
    for (ta in c(0, 0.5, 1, 2, 4)) for (ct in c(1.5, 2.5, 4))
        add(cell("M2", sprintf("tau%.1f_c%.1f", ta, ct),
                 dgp = list(tau = ta, contrast = ct, share_types = 1),
                 cfg = list(levels = 1L)))

    ## M3: L2 types and handoff. K and delta one at a time; contrast for handoff.
    for (K in c(2L, 3L, 5L))
        add(cell("M3", sprintf("K%d", K), dgp = list(K_types = K), cfg = list(levels = 1:2)))
    for (de in c(0.5, 1, 3, 4))
        add(cell("M3", sprintf("delta%.1f", de), dgp = list(delta = de), cfg = list(levels = 1:2)))
    for (ct in c(1.5, 4))
        add(cell("M3", sprintf("c%.1f", ct), dgp = list(contrast = ct), cfg = list(levels = 1:2)))

    ## M4: meta-types. share_types x distance method.
    for (sh in c(0, 1 / 3, 2 / 3, 1)) for (me in c("bhattacharyya", "wasserstein"))
        add(cell("M4", sprintf("share%.2f_%s", sh, substr(me, 1, 4)),
                 dgp = list(share_types = sh), cfg = list(levels = 1:2, meta_method = me)))

    ## M5: L3 configuration. Types fully shared so zones cut across regions
    ## cleanly. "null" = all zones random, so zones are unrecoverable by design.
    zsets <- list(mixed = c("block", "fine", "banded", "random"),
                  null  = rep("random", 4))
    for (zs in names(zsets)) for (ko in 1:3) for (en in c("pam", "dmm"))
        add(cell("M5", sprintf("%s_k%d_%s", zs, ko, en),
                 dgp = list(zones = zsets[[zs]], share_types = 1),
                 cfg = list(levels = 1:3, k_order = ko, mosaic_engine = en)))
    d
}

#' Flatten a design to one row per cell (vectors collapsed with "+").
design_table <- function(design) {
    flat <- function(x) vapply(x, function(v) paste(v, collapse = "+"), character(1))
    rows <- lapply(design, function(c) {
        data.frame(cell_id = c$id, module = c$module, label = c$label,
                   t(flat(c$dgp)), t(flat(c$cfg)), check.names = FALSE)
    })
    cols <- unique(unlist(lapply(rows, names)))
    do.call(rbind, lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
}

# --------------------------------------------------------------------------- *
# Runner
# --------------------------------------------------------------------------- *

#' Run one (cell, rep) and save it. Returns the file path.
run_one <- function(cell, rep, out_dir, master_seed = 20261105L, overwrite = FALSE) {
    dir.create(file.path(out_dir, cell$id), recursive = TRUE, showWarnings = FALSE)
    f <- file.path(out_dir, cell$id, sprintf("rep_%03d.rds", rep))
    if (file.exists(f) && !overwrite) return(invisible(f))

    dgp_seed  <- .hash_int(list(master_seed, rep, cell$dgp))
    pipe_seed <- .hash_int(list(master_seed, rep, cell$id))
    t0 <- proc.time()[["elapsed"]]
    sim <- tryCatch(do.call(simulate_landscape, c(cell$dgp, list(seed = dgp_seed))),
                    error = function(e) e)
    res <- if (inherits(sim, "error")) {
        list(metrics = data.frame(dgp_error = conditionMessage(sim)),
             units = NULL, stages = NULL, reliability = NULL)
    } else {
        set.seed(pipe_seed)                       # stream for floor permutation
        run_pipeline(sim, cell$cfg, seed = pipe_seed)
    }
    res$id <- list(cell_id = cell$id, module = cell$module, label = cell$label,
                   rep = rep, dgp_seed = dgp_seed, pipe_seed = pipe_seed,
                   seconds = proc.time()[["elapsed"]] - t0,
                   pulseR_version = as.character(utils::packageVersion("pulseR")),
                   R_version = R.version.string, time = format(Sys.time()))
    saveRDS(res, f)
    invisible(f)
}

#' All (cell, rep) jobs, in an order that spreads cells across chunks.
job_list <- function(design, reps) {
    jb <- expand.grid(rep = seq_len(reps), cell = names(design), stringsAsFactors = FALSE)
    jb[order(jb$rep, jb$cell), ]
}

#' Run every job. cores > 1 uses forked workers (Linux, macOS).
run_design <- function(design, reps, out_dir, cores = 1L, master_seed = 20261105L) {
    jb  <- job_list(design, reps)
    fun <- function(i) tryCatch(
        run_one(design[[jb$cell[i]]], jb$rep[i], out_dir, master_seed),
        error = function(e) { message("job ", i, " failed: ", conditionMessage(e)); NA })
    saveRDS(design, file.path(.mkdir(out_dir), "design.rds"))
    if (cores > 1L && .Platform$OS.type == "unix")
        invisible(parallel::mclapply(seq_len(nrow(jb)), fun, mc.cores = cores,
                                     mc.preschedule = FALSE))
    else invisible(lapply(seq_len(nrow(jb)), fun))
}

#' Run one chunk of the job list (for Slurm arrays).
run_chunk <- function(design, reps, chunk, n_chunks, out_dir, master_seed = 20261105L) {
    jb  <- job_list(design, reps)
    idx <- which((seq_len(nrow(jb)) - 1L) %% n_chunks + 1L == chunk)
    saveRDS(design, file.path(.mkdir(out_dir), "design.rds"))
    for (i in idx) run_one(design[[jb$cell[i]]], jb$rep[i], out_dir, master_seed)
    invisible(idx)
}

# --------------------------------------------------------------------------- *
# Collector
# --------------------------------------------------------------------------- *

#' Read all results. Returns metrics (one row per run, joined to the design
#' table), stages (one row per run x stage) and, if units = TRUE, per-unit rows.
collect_results <- function(out_dir, units = FALSE) {
    files <- list.files(out_dir, pattern = "^rep_.*\\.rds$", recursive = TRUE,
                        full.names = TRUE)
    rs <- lapply(files, readRDS)
    idf <- function(r) as.data.frame(r$id[c("cell_id", "rep", "seconds")])
    met <- .rbind_fill(lapply(rs, function(r) cbind(idf(r), r$metrics)))
    stg <- .rbind_fill(lapply(rs, function(r) if (!is.null(r$stages)) cbind(idf(r)[1:2], r$stages)))
    dsg <- file.path(out_dir, "design.rds")
    if (file.exists(dsg)) met <- merge(design_table(readRDS(dsg)), met, by = "cell_id")
    out <- list(metrics = met, stages = stg)
    if (units) out$units <- .rbind_fill(lapply(rs, function(r)
        if (!is.null(r$units)) cbind(idf(r)[1:2], r$units)))
    out
}

#' Monte Carlo SE and replicates needed for a target MCSE, per cell.
mcse_table <- function(metrics, metric, target = 0.03) {
    sp <- split(metrics[[metric]], metrics$cell_id)
    out <- data.frame(cell_id = names(sp),
                      n    = vapply(sp, function(v) sum(is.finite(v)), integer(1)),
                      mean = vapply(sp, function(v) mean(v, na.rm = TRUE), numeric(1)),
                      sd   = vapply(sp, function(v) stats::sd(v, na.rm = TRUE), numeric(1)),
                      row.names = NULL)
    out$mcse        <- out$sd / sqrt(out$n)
    out$reps_needed <- ceiling((out$sd / target)^2)
    out
}

# --------------------------------------------------------------------------- *
# Helpers
# --------------------------------------------------------------------------- *

# djb2 string hash of any R object, as a positive integer seed.
.hash_int <- function(x) {
    s <- utf8ToInt(paste(deparse(x), collapse = ""))
    h <- 5381
    for (ch in s) h <- (h * 33 + ch) %% 2147483629
    as.integer(h)
}

.mkdir <- function(p) { dir.create(p, recursive = TRUE, showWarnings = FALSE); p }

.rbind_fill <- function(lst) {
    lst <- Filter(Negate(is.null), lst)
    if (!length(lst)) return(NULL)
    cols <- unique(unlist(lapply(lst, names)))
    do.call(rbind, lapply(lst, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}
