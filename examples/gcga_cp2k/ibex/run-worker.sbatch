#!/bin/bash
# Ibex: one GCGA worker per job (16 cores, shared Turin node). Submit N copies: for i in $(seq N); do sbatch ibex/run-worker.sbatch; done
# input.py: cp2k_cmd = 'mpirun -np 16 /ibex/user/reny0b/zls/soft/cp2k-2026.2-lpb/install_z4/bin/cp2k.psmp'
#SBATCH --job-name=gocia_cp2k_w
#SBATCH --partition=batch
#SBATCH --constraint=cpu_amd_epyc_9655
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
cd "$SLURM_SUBMIT_DIR" || exit 1
source ibex/env.sh
${GOCIA_PY:-python3} -u ibex/ga-loop.py ${TOTCONF:-2000} ${MINCONF:-200}
