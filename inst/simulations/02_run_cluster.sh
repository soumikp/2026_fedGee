#!/bin/bash
#
# FedGEE Simulation — Pitt CRC SLURM submission
#
# Usage:
#   Step 1: Run the setup script to count tasks
#           Rscript 02_run_cluster.R count
#
#   Step 2: Submit
#           sbatch run_sim.slurm
#
#   Step 3: Monitor
#           squeue -u $USER
#           sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS
#
###############################################################################

#SBATCH --job-name=fedgee_sim
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --array=1-300
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-02:00:00
#SBATCH --output=/ihome/spurkayastha/soumik/2026_fedGee/simulations/log/sim_%A_%a.out
#SBATCH --error=/ihome/spurkayastha/soumik/2026_fedGee/simulations/log/sim_%A_%a.err
#SBATCH --mail-user=soumik@pitt.edu
#SBATCH --mail-type=END,FAIL

###############################################################################
# Notes on parameter choices:
#
# --cluster=smp        Single-node R jobs, no MPI needed
# --array=1-1000%100   1000 array tasks, max 100 running concurrently
#                      Each task handles ceil(total_tasks / 1000) reps
#                      Adjust 1000 to match your grid (see 02_run_cluster.R count)
# --mem=8G             R + geeglm + glmer can spike; 8G is safe
#                      Monitor with sacct and reduce if wasteful
# --time=0-04:00:00    4 hours per array task; adjust after timing one job
#                      If > 3 days, add: #SBATCH --qos=long
# %100                 Polite concurrency cap; increase if cluster is quiet
###############################################################################

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

# Set working directory to where scripts live
cd ${SLURM_SUBMIT_DIR}

# Run the simulation
Rscript 02_run_cluster.R ${SLURM_ARRAY_TASK_ID}

echo "End Time: $(date)"
echo "Exit Code: $?"
