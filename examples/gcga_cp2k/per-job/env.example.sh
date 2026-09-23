# Runtime for GOCIA CP2K workers, sourced by per-job/*.sbatch. Copy to per-job/env.sh and fill in for your cluster.
# How the settings below were chosen and how to re-check them on new hardware: docs/SPEED_TUNING.md.
CP2K_ROOT=/REPLACE/cp2k                     # install with bin/cp2k.psmp and lib/
# MPI / compiler / MKL environment of the CP2K build, e.g. `module load ...` or `source .../mpi/env/vars.sh`
# REPLACE: module load <compiler> <mpi> <mkl>
export CP2K_DATA_DIR=$CP2K_ROOT/data
export LD_LIBRARY_PATH=$CP2K_ROOT/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# Optional, AMD CPUs with AVX-512 only, after `tools/speedtest/speedtest.py ab` showed a gain with identical energies:
# export LD_PRELOAD=/REPLACE/libmkl_amd_shim.so MKL_ENABLE_INSTRUCTIONS=AVX512     # tools/speedtest/mkl_amd_shim.c
# python with ase + natsort and this repository on the compute nodes
export PYTHONPATH=${GOCIA_REPO:-/REPLACE/gocia}:${GOCIA_PYLIB:-}:$PYTHONPATH
ulimit -s unlimited
