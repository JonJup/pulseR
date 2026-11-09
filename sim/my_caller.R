source("sim_dgp.R")
source("sim_metrics.R")
source("sim_pipeline.R")
source("sim_run.R")
source("sim_smoke.R")




design <- make_design()                       # M1-M5 cells, Tier I
run_design(design, reps = 20, out_dir = "out", cores = 8)
res <- collect_results("out")                 # metrics, stages, units