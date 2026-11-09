# =========================================================================== *
# sim_smoke.R
#
# End-to-end check of the harness on a tiny design. Run before every Tier I
# launch and after every pulseR update. About 4 minutes on one core.
#
#   source("sim_smoke.R"); smoke_test()
#
# Passes when: every stage of every run is ok, the key metric columns exist,
# a rerun skips finished jobs (resume), and collect_results() joins the design.
# =========================================================================== *

for (f in c("sim_dgp.R", "sim_metrics.R", "sim_pipeline.R", "sim_run.R"))
    source(f)

smoke_design <- function() {
    fast <- list(N_tree = 60L, regions_per_partition = 20L, stability_B = 2L,
                 n_local_types = 1:4, n_mosaic = 2:5)
    small <- list(N = 600L, R_true = 4L, K_types = 2L)
    d <- list(
        cell("S", "l1_null",   dgp = c(small, contrast = 0, share_types = 1),
             cfg = c(fast, levels = 1L, floor = TRUE, null_B = 3L)),
        cell("S", "l1_signal", dgp = c(small, contrast = 3),
             cfg = c(fast, levels = 1L, floor = TRUE)),
        cell("S", "l1_tuned",  dgp = small,
             cfg = c(fast, levels = 1L, list(final_regions = 3:5))),
        cell("S", "l2",        dgp = small, cfg = c(fast, levels = list(1:2))),
        cell("S", "l3",        dgp = c(small, share_types = 1,
                                       list(zones = c("block", "random"))),
             cfg = c(fast, levels = list(1:3), mosaic_engine = "pam")))
    stats::setNames(d, vapply(d, `[[`, "", "id"))
}

smoke_test <- function(out_dir = file.path(tempdir(), "pulse_smoke")) {
    unlink(out_dir, recursive = TRUE)
    d  <- smoke_design()
    t1 <- system.time(run_design(d, reps = 1L, out_dir = out_dir))[["elapsed"]]
    t2 <- system.time(run_design(d, reps = 1L, out_dir = out_dir))[["elapsed"]]
    r  <- collect_results(out_dir)

    bad <- r$stages[!r$stages$ok, c("cell_id", "stage", "error")]
    need <- c("l1_ari", "l1_auc", "ari_l1_floor", "null_p_ch", "l2_ari",
              "l2o_ari", "l2_handoff", "meta_ari", "l3_ari", "l3o_ari")
    miss <- setdiff(need, names(r$metrics))
    checks <- c(all_stages_ok  = nrow(bad) == 0L,
                metrics_present = length(miss) == 0L,
                one_row_per_run = nrow(r$metrics) == length(d),
                resume_skips   = t2 < 5,
                design_joined  = "module" %in% names(r$metrics))
    cat(sprintf("first run %.0f s, rerun %.1f s\n", t1, t2))
    print(checks)
    if (nrow(bad)) print(bad)
    if (length(miss)) cat("missing metrics:", miss, "\n")
    invisible(list(checks = checks, results = r))
}
