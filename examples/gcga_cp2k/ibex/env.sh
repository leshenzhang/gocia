# Ibex CP2K runtime for GOCIA workers (sourced by the sbatch files). Measured choices (ai-reconstr-pre D-010/D-011):
# install_z4 = gcc 13.2 -march=znver4 + ELPA 2025.06.001 + the redistribution fix; MKL AVX-512 kernels on AMD via the
# libfakeintel.so shim (mkl_serv_intel_cpu_true -> 1): 0.80x wall of the x86-64-v3 ScaLAPACK build, energies identical.
R=/ibex/user/reny0b/zls/soft/cp2k-2026.2-lpb
source /sw/rl9g/intel/2022/mpi/2021.7.1/env/vars.sh
export CP2K_DATA_DIR=$R/src/cp2k-2026.2/data
export LD_LIBRARY_PATH=/sw/rl9g/intel/2022/mkl/2022.2.1/lib/intel64:/sw/rl9g/fftw/3.3.10/rl9_gnu11.3_ompi4.1.4_dp/lib:$R/install_z4/lib:/sw/rl9g/gcc/13.2.0/rl9_binary/lib64:$LD_LIBRARY_PATH
export I_MPI_HYDRA_BOOTSTRAP=fork OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MKL_ENABLE_INSTRUCTIONS=AVX512 LD_PRELOAD=$R/tools/libfakeintel.so
# python with ase + natsort + this repo on the compute nodes: set GOCIA_PY (e.g. a --target pylib); not provided here.
export PYTHONPATH=${GOCIA_REPO:-/ibex/user/wangc0i/zls/soft/gocia-cp2k}:${GOCIA_PYLIB:-}:$PYTHONPATH
ulimit -s unlimited
