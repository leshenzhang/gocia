#!/bin/bash
# Relax the initial population with CP2K, 4 init-workers per node (Shaheen wangc0i pack-4 rule).
#   python boxSample.py substrate.vasp 40 && python db2dirs.py init.db && sbatch run-init.sbatch
# Each s0000N/ runs ../init-worker.py which writes its relaxed structure into ../gcga.db (done=1).
#SBATCH --job-name=01__gocia_init
#SBATCH --account=k10175
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
module switch PrgEnv-cray PrgEnv-gnu
module load cray-fftw/3.3.10.10
module load cray-python
R=/scratch/wangc0i/zls/soft/cp2k-2026.2
export CP2K_DATA_DIR=$R/src/cp2k-2026.2/data
export OMP_NUM_THREADS=2
export PYTHONPATH=/scratch/wangc0i/zls/soft/gocia-cp2k:/scratch/wangc0i/zls/soft/pylib:$PYTHONPATH
NPAR=$(( 4 * SLURM_JOB_NUM_NODES ))
cd "$SLURM_SUBMIT_DIR" || exit 1
one() { ( cd $1 && python -u ../init-worker.py > worker.log 2>&1; echo "### $1 rc=$?" ); }
for d in s0*/; do
  d=${d%/}
  [ -f $d/FAIL ] && rm -f $d/FAIL
  grep -q "written:" $d/worker.log 2>/dev/null && { echo "### skip $d (done)"; continue; }
  one $d &
  while (( $(jobs -rp | wc -l) >= NPAR )); do wait -n; done
done
wait
echo "### ALL DONE"
