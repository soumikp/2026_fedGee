#!/bin/bash
# Monitor the v3 array job.
#   bash check_progress.sh
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
n_done=$(ls "${PROJECT_ROOT}/output"/results_[0-9]*.rds 2>/dev/null | wc -l)
echo "Result chunks written : ${n_done} / 499"
echo "Queue:"
squeue -M smp -u "$USER" -o "%.10i %.12j %.8T %.10M %.6D %R" 2>/dev/null | head -20
echo
echo "Recent errors (non-empty .err files):"
find "${PROJECT_ROOT}/logs" -name 'v3_*.err' -size +0c 2>/dev/null | head -5
