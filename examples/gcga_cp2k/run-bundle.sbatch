#!/bin/bash
# Shaheen (wangc0i, k10175) bundled GCGA job: 1 node = 192 cores = 4 CP2K workers x 48 cores.
# Scale: --nodes=N and nworker=4*N. Job name keeps the shared-budget prefix 01__.
#SBATCH --job-name=01__gocia_cp2k
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
R=/scratch/wangc0i/zls/soft/cp2k-2026.2
export CP2K_DATA_DIR=$R/src/cp2k-2026.2/data
export OMP_NUM_THREADS=2
# python with ase+numpy+natsort visible from the compute nodes (/scratch, not /home: Cray DVS)
export PYTHONPATH=/scratch/wangc0i/zls/soft/gocia-cp2k:$PYTHONPATH
PY=/scratch/reny0b/iops/sw/envs/gs/bin/python
NW=$(( 4 * SLURM_JOB_NUM_NODES ))
cd "$SLURM_SUBMIT_DIR" || exit 1
$PY -u ga-bundle.py $NW 400 100
echo "### ALL DONE"
