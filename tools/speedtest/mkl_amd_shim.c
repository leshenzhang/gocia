/* MKL chooses its code path from the CPU vendor. On AMD CPUs with AVX-512 (e.g. Zen 4/5) it may stay on AVX2
 * kernels. This shim makes MKL's Intel check return true so the AVX-512 kernels are used:
 *     gcc -O2 -shared -fPIC -o libmkl_amd_shim.so mkl_amd_shim.c
 *     export LD_PRELOAD=/path/libmkl_amd_shim.so MKL_ENABLE_INSTRUCTIONS=AVX512
 * It changes only MKL's dispatch. Before adopting it, confirm with `speedtest.py ab` that energies are unchanged
 * (0.0 meV) and that the run is faster on your nodes; on Intel CPUs it does nothing. */
int mkl_serv_intel_cpu_true(void) { return 1; }
