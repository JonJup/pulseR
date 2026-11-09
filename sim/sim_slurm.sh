#!/bin/bash
# Tier I on Turso. Edit partition, time and image path, then:
#   sbatch --array=1-200 sim_slurm.sh
# One task = 1/200 of all (cell, rep) jobs; finished jobs are skipped on rerun.
#SBATCH --job-name=pulse_tier1
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%A_%a.out

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # no hidden BLAS threading
OUT=$WRKDIR/pulse_sim/tier1
mkdir -p "$OUT" logs
singularity exec "$WRKDIR/images/pulse.sif" Rscript sim_slurm.R "$OUT" 20
