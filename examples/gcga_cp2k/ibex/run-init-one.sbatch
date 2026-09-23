#!/bin/bash
# Ibex: relax ONE initial structure (3 stages) and write it into gcga.db. From the search root:
#   for d in s0*/; do sbatch --chdir=$d ibex/run-init-one.sbatch; done      (one 16-core job per structure; path from the root)
#SBATCH --job-name=gocia_cp2k_init
#SBATCH --partition=batch
#SBATCH --constraint=cpu_amd_epyc_9655
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=init-%j.out
source ../ibex/env.sh
${GOCIA_PY:-python3} -u ../init-worker.py > worker.log 2>&1
