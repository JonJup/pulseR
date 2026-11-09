
library(devtools)
load_all()
setwd("../simulation_tests/")
source("sim_gate.R")

v <- run_gate(reps = 6, cores = 4, pkg_path = "../pulseR/", master_seed = 2L,
              cfg  = list(weight_exponent = 4, size_influence = 0,
                          regions_per_partition = 4L),
              opts = list(pulseR.cut_rule = "ssd", pulseR.refine = TRUE))

# v <- run_gate(reps = 6, cores = 4, pkg_path = "../pulseR/", master_seed = 1L,
#               cfg  = list(weight_exponent = 4, size_influence = 0,
#                           regions_per_partition = 4L),
#               opts = list(pulseR.cut_rule = "ssd"))
# v <- run_gate(reps = 6, cores = 4, pkg_path = "../pulseR/", master_seed = 1L,
#               cfg  = list(weight_exponent = 4, size_influence = 0,
#                           regions_per_partition = 4L, refine = TRUE),
#               opts = list(pulseR.cut_rule = "ssd"))
# r <- collect_results(v$run$out_dir)
# aggregate(cbind(l1_ari, l1_ari_refined, l1_refined_pieces, l1_refine_moves) ~ cell_id,
#           r$metrics, mean)


r <- collect_results("gate_out/20261106-143057", units = TRUE)
u <- r$units[r$units$cell_id == "G_c4", ]
u$band <- cut(u$d_boundary, c(-Inf, 1, 2, 4, 8, Inf))
rbind(n      = table(u$band),
      err    = round(tapply(u$l1_wrong, u$band, mean), 3),
      share  = round(tapply(u$l1_wrong, u$band, sum) / sum(u$l1_wrong), 3))
aggregate(l1_ari ~ rep, r$metrics[r$metrics$cell_id == "G_c4", ], mean)

for (f in c("sim_pipeline.R", "sim_metrics.R", "sim_gate.R"))
        cat(f, grep(".refine_regions <- function", readLines(f), fixed = TRUE), "\n")

# test experimental setup  ------------------------------------------------
options(pulseR.cut_rule = "ssd")
gd <- gate_design(); sw <- list()
for (cl in c("c4", "null")) for (g in c(1, 2, 4)) {
        x <- gd[[paste0("G_", cl)]]
        x$cfg <- modifyList(x$cfg, list(weight_exponent = g, size_influence = 0,
                                        regions_per_partition = 4L,
                                        baselines = FALSE, floor = FALSE))
        x$label <- sprintf("%s_ssd_g%g", cl, g); x$id <- paste0("D_", x$label)
        sw[[x$id]] <- x
}
out <- file.path("gate_out", "diag_ssd")
for (r in 1:5) for (x in sw) run_one(x, r, out)
options(pulseR.cut_rule = NULL)
res <- collect_results(out)
aggregate(cbind(l1_ari, l1_mean_cert, l1_n_eff) ~ cell_id, res$metrics, mean)



# original ----------------------------------------------------------------



v <- run_gate(reps = 25, cores = 4, pkg_path = "../pulseR/")




# check up 1 --------------------------------------------------------------

d0   <- gate_design()$G_c4
grid <- expand.grid(si = c(0, 1, 2), rpp = c(4L, 20L))
sw   <- list()
for (j in seq_len(nrow(grid))) {
        x <- d0
        x$cfg   <- modifyList(x$cfg, list(size_influence = grid$si[j],
                                          regions_per_partition = grid$rpp[j],
                                          baselines = FALSE))
        x$label <- sprintf("c4_si%g_rpp%d", grid$si[j], grid$rpp[j])
        x$id    <- paste0("D_", x$label)
        sw[[x$id]] <- x
}
out <- file.path("gate_out", "diag_size")
for (r in 1:5) for (x in sw) run_one(x, r, out)
res <- collect_results(out)
aggregate(cbind(l1_ari, l1_mean_cert) ~ cell_id, res$metrics, mean)



# check up two ------------------------------------------------------------

d0 <- gate_design()$G_c4
sw <- list()
for (g in c(1, 4, 8, 16)) {
        x <- d0
        x$cfg   <- modifyList(x$cfg, list(weight_exponent = g, size_influence = 0,
                                          regions_per_partition = 4L, baselines = FALSE))
        x$label <- sprintf("c4_g%g", g)
        x$id    <- paste0("D_", x$label)
        sw[[x$id]] <- x
}
out <- file.path("gate_out", "diag_gamma")
for (r in 1:5) for (x in sw) run_one(x, r, out)
res <- collect_results(out)
aggregate(cbind(l1_ari, l1_mean_cert, l1_n_eff) ~ cell_id, res$metrics, mean)
