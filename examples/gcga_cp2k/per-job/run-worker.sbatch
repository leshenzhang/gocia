#!/bin/bash
# One GCGA worker per scheduler job (16 cores on a shared node; no whole-node request). Submit N copies for N workers:
#   for i in $(seq N); do sbatch per-job/run-worker.sbatch; done
# input.py: cp2k_cmd = 'mpirun -np 16 <CP2K_ROOT>/bin/cp2k.psmp'. Core count per worker: docs/SPEED_TUNING.md.
#SBATCH --job-name=gocia_cp2k_w
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
##SBATCH --partition=REPLACE
##SBATCH --constraint=REPLACE_NODE_FEATURE     # pin the node type chosen by the speed tests
cd "$SLURM_SUBMIT_DIR" || exit 1
source per-job/env.sh
${GOCIA_PY:-python3} -u per-job/ga-loop.py ${TOTCONF:-2000} ${MINCONF:-200}
