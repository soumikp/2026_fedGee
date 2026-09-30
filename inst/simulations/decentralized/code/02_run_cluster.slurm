#!/bin/bash
#
# FedGEE Simulation — Pitt CRC SLURM submission
#
# Usage:
#   Step 1: Time ONE task locally first (the 3 decentralized methods each
#           cost O(dec_n_iter * K) local GEE evaluations per call -- much
#           slower per call than the closed-form centralized methods):
#             Rscript 02_run_cluster.R local 1 test
#
#   Step 2: Check grid size / method count for THIS pipeline
#           Rscript 02_run_cluster.R count
#
#   Step 3: Submit
#           sbatch run_sim.slurm
#
#   Step 4: Monitor
#           squeue -u $USER
#           sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS
#
###############################################################################
#SBATCH --job-name=fedgee_decentralized
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --array=1-300
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-06:00:00
#SBATCH --output=/ihome/spurkayastha/fec54/Fed-GEE/One-shot/Decentralized/simulations/log/sim_%A_%a.out
#SBATCH --error=/ihome/spurkayastha/fec54/Fed-GEE/One-shot/Decentralized/simulations/log/sim_%A_%a.err
#SBATCH --mail-user=FEC54@pitt.edu
#SBATCH --mail-type=END,FAIL

###############################################################################
# Notes on parameter choices (updated for the CURRENT pipeline: run_one_rep()
# in 00_simulation_config.R now fits 6 centralized benchmark methods PLUS 3
# decentralized Dec-Fed-GEE methods -- one per network topology (dec_hub,
# dec_ring, dec_visn) via decentralized_fedgee() in 01_decentralized_fedgee.R
# -- 9 result rows per replication, not 6, not 8. There is no `dec_configs`
# list and no "complete" topology in this pipeline -- those belonged to an
# earlier iteration of the decentralized code and no longer apply):
#
# --cluster=smp        Single-node R jobs, no MPI needed.
#
# --array=1-300         PLACEHOLDER. With n_reps = 200 (set in
#                      02_run_cluster.R) and 9 result rows/rep, re-run
#                      `Rscript 02_run_cluster.R count` for the CURRENT grid
#                      and set --array to whatever it recommends, then
#                      adjust the %concurrency-cap suffix (e.g.
#                      --array=1-300%50) if you want to be polite to a
#                      shared partition.
#
# --mem=8G             R + geepack + glmer/lme4 can spike; the decentralized
#                      loop itself is mostly small-matrix (p x p) operations
#                      per site so it should NOT meaningfully raise memory
#                      use -- 8G should still be fine, but confirm via
#                      `sacct --format=...,MaxRSS` on your first completed
#                      array task before trusting this for a full run.
#
# --time=0-06:00:00     Each of the 3 decentralized methods costs
#                      O(dec_n_iter * K) local GEE evaluations per call
#                      (default dec_n_iter = 100, set in run_one_rep()'s
#                      signature in 00_simulation_config.R), run once per
#                      topology (hub, ring, visn) on top of the 6 closed-
#                      form centralized methods -- meaningfully more CPU
#                      time per task than the closed-form-only pipeline.
#                      6 hours is a SAFER STARTING GUESS, not a measured
#                      value. Time one real task first
#                      (`Rscript 02_run_cluster.R local 1 test`, or submit a
#                      single array index and check
#                      `sacct ... --format=Elapsed`) and shrink/grow this
#                      accordingly. If a chunk would exceed the partition's
#                      time limit, either reduce the chunk size (larger
#                      --array) or add:
#                        #SBATCH --qos=long
#
# %100                 (If added to --array) polite concurrency cap; raise
#                      if the cluster is quiet, lower if you are asked to.
###############################################################################

# -----------------------------------------------------------------------------
# PRE-FLIGHT CHECK (do this once before your first sbatch on a new cluster):
#   02_run_cluster.R sources THREE files by hard-coded ABSOLUTE path (there
#   is no single BASE_DIR variable in this script -- unlike an earlier
#   iteration of this pipeline). Open 02_run_cluster.R and confirm all three
#   paths actually point to real files:
#     00_simulation_config.R
#     01_fedgee.R                 (centralized/pooled FedGEE -- defines
#                                   fedgee(), used by the fed_site /
#                                   fed_site_md / fed_patient methods)
#     01_decentralized_fedgee.R   (defines decentralized_fedgee(), used by
#                                   the dec_hub / dec_ring / dec_visn methods)
#
#   NOTE: the SLURM-array-mode source() block and the local-parallel-mode
#   clusterEvalQ() source() block in 02_run_cluster.R currently point to TWO
#   DIFFERENT directories (this predates these changes -- not introduced by
#   them). Since THIS submission script only ever runs SLURM-array mode
#   (see the Rscript call at the bottom), only that block's paths matter for
#   sbatch -- but if you also plan to run `Rscript 02_run_cluster.R local N`
#   directly (not via sbatch), check the other block's paths too.
#
#   01_fedgee.R calls `library(Rcpp)` at the top (even though no compiled
#   code is actually invoked at runtime) -- confirm the Rcpp package is
#   installed in the R module you `module load` below, or sourcing will
#   fail immediately with "there is no package called 'Rcpp'".
# -----------------------------------------------------------------------------

# Create log directory if it doesn't exist
mkdir -p log

# Load R module (check available versions with: module spider R)
module load gcc
module load r

# Print job info for debugging
echo "============================================"
echo "SLURM Job ID:    ${SLURM_JOB_ID}"
echo "Array Task ID:   ${SLURM_ARRAY_TASK_ID}"
echo "Node:            ${SLURMD_NODENAME}"
echo "Working Dir:     $(pwd)"
echo "Start Time:      $(date)"
echo "============================================"

# Set working directory to where scripts live.
cd ${SLURM_SUBMIT_DIR}

# Run the simulation
Rscript 02_run_cluster.R ${SLURM_ARRAY_TASK_ID}

echo "End Time: $(date)"
echo "Exit Code: $?"
