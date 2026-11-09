# Entry point for one Slurm array task. Usage:
#   Rscript sim_slurm.R <out_dir> <reps>
# Chunk = SLURM_ARRAY_TASK_ID, n_chunks = SLURM_ARRAY_TASK_COUNT.
suppressMessages({ library(sf); library(pulseR) })
for (f in c("sim_dgp.R", "sim_metrics.R", "sim_pipeline.R", "sim_run.R")) source(f)
a <- commandArgs(trailingOnly = TRUE)
run_chunk(make_design(), reps = as.integer(a[2]),
          chunk    = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")),
          n_chunks = as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT")),
          out_dir  = a[1])
