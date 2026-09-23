# Tuned VASP stage INCARs for GOCIA relaxations

Drop-in `INCAR-1/2/3` for the GOCIA VASP workflow (`examples/gcga`), from the speed campaign in
`docs/SPEED_TUNING.md` (section 7). Physics (functional, smearing, cutoff, `EDIFFG` 0.5 / 0.2 / 0.1 eV/A, `PREC`)
is that of the project protocol; change it for your system. The performance settings are:

| setting | where | why |
|---|---|---|
| `NCORE = 8`, `NSIM = 16` | all stages | measured best at 16 MPI ranks (NSIM 16: 0.87x wall, identical results). NCORE depends on the rank count: re-test with `tools/speedtest` when the ranks change |
| `EDIFF = 1e-4` | all stages | 0.90x wall at 0.41 meV / 7e-4 A after 5 steps (stages 2-3 were 1e-5) |
| `IBRION = 3`, `IOPT = 7` (FIRE), `POTIM = 0`, `MAXMOVE = 0.2`, `TIMESTEP = 0.1` | all stages | about half the SCF work of CG in stage 2; `IBRION = 1` diverged far from a minimum and crawled after a FIRE stage. Needs a VASP build with the VTST optimisers; without it use `IBRION = 2`, `POTIM = 0.5` (CG) |
| `LWAVE`/`LCHARG = .TRUE.` in stage 2, `ISTART = 1`, `ICHARG = 1` in stage 3 | stages 2-3 | stage 3 starts from the stage-2 wavefunction and charge |

Run VASP pure MPI (OpenMP was 1.7-2.3x slower). On AMD CPUs with AVX-512, test the MKL code-path shim
(`tools/speedtest/mkl_amd_shim.c`; 0.90x in our runs, identical energies).
