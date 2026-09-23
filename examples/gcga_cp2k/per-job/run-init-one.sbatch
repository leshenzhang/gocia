#!/bin/bash
# Relax ONE initial structure (3 stages) and write it into gcga.db. From the search root, one job per structure:
#   for d in s0*/; do sbatch --chdir=$d per-job/run-init-one.sbatch; done
#SBATCH --job-name=gocia_cp2k_init
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=init-%j.out
##SBATCH --partition=REPLACE
##SBATCH --constraint=REPLACE_NODE_FEATURE
source ../per-job/env.sh
${GOCIA_PY:-python3} -u ../init-worker.py > worker.log 2>&1
