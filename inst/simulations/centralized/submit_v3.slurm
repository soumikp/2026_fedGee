#!/bin/bash
###############################################################################
# Fed-GEE simulation v3 (experiment EV3) -- Pitt CRC SLURM submission
#
# Adds, relative to the calibration run:
#   * TRUE site-level Fay-Graubard (the v2 "fg" was observation-level HC2)
#   * Bell-McCaffrey df in score space, and FG's own df (saws, dtildeH)
#   * t(K-1) reference alongside z / t(K-p) / t(nu)
#   * UNEQUAL site sizes in the small-K regime (CV ~0.78 VHA-matched, ~1.0)
#
# Grid: 117 scenarios x 1000 reps = 117,000 tasks.
#   K        in {10, 18, 30, 50}
#   p total  in {2, 5, 10, 20}   (feasible cells K - p >= 5)
#   beta_0   in {-3, -1, 0}      (~6% / ~27% / ~50% prevalence)
#   sizes    in {uniform, cv078, cv100}, mean ~20 patients/site
#   rho_H = 0.05, rho_P = 0.30
# Measured locally: 0.3-0.5 s/task at K <= 18, ~2 s/task at K = 50, p = 20.
# 117000 / 499 = 235 tasks/job -> ~2-8 min/job.
#
# BEFORE FIRST LAUNCH (once):  module load gcc r;  Rscript -e 'install.packages("saws")'
#
# Usage from this directory:
#   cd /ihome/spurkayastha/soumik/2026_fedGee/sim_v3
#   sbatch submit_v3.slurm
#   squeue -M smp -u $USER
#   bash check_progress.sh
#   FEDGEE_PROJECT_ROOT=$(pwd) Rscript code/collect_results.R   # when done
###############################################################################

#SBATCH --job-name=fedgee_v3
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --account=spurkayastha
#SBATCH --array=1-499
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/ihome/spurkayastha/soumik/2026_fedGee/sim_v3/logs/v3_%A_%a.out
#SBATCH --error=/ihome/spurkayastha/soumik/2026_fedGee/sim_v3/logs/v3_%A_%a.err
#SBATCH --mail-user=soumik@pitt.edu
#SBATCH --mail-type=END,FAIL

PROJECT_ROOT=/ihome/spurkayastha/soumik/2026_fedGee/sim_v3

export FEDGEE_SIM_DIR=${PROJECT_ROOT}/code
export FEDGEE_OUT_DIR=${PROJECT_ROOT}/output
export FEDGEE_EXPS=EV3
export FEDGEE_N_REPS=${FEDGEE_N_REPS:-1000}

mkdir -p "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/output"

module load gcc
module load r

echo "============================================"
echo "SLURM Job ID:  ${SLURM_JOB_ID}"
echo "Array Task:    ${SLURM_ARRAY_TASK_ID} / ${SLURM_ARRAY_TASK_MAX}"
echo "Node:          ${SLURMD_NODENAME}"
echo "Experiments:   ${FEDGEE_EXPS}"
echo "Reps:          ${FEDGEE_N_REPS}"
echo "Start:         $(date)"
echo "============================================"

cd "${PROJECT_ROOT}"
Rscript "${FEDGEE_SIM_DIR}/02_run_cluster_v3.R" "${SLURM_ARRAY_TASK_ID}"
STATUS=$?

echo "End:           $(date)"
echo "Exit Code:     ${STATUS}"
exit ${STATUS}
