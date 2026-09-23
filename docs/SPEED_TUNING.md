# Speed tuning of VASP / CP2K relaxations for GOCIA on a new cluster

Chinese version: `SPEED_TUNING.zh.md`. Tools: `tools/speedtest/`. Worked reference result: section 7.

GOCIA spends almost all of its time in local relaxations of offspring structures. The quantity to maximise is
**relaxed structures per hour for a given core budget, at fixed accuracy**. A setting is adopted only if it is
faster *and* passes the accuracy gates of section 4. Everything below is scheduler- and cluster-independent; the
only cluster-specific file is the JSON config of `tools/speedtest/`.

## 1. What you need

- **Stage templates** of the production protocol: `INCAR-1..3` (VASP) or `cp2k-1..3.inp` (CP2K, GOCIA template
  contract with `@INCLUDE` / `@SET`).
- **One pre-relaxed structure** of the production system (per-step tests, section 3 phases B-G).
- **At least four unrelaxed GA offspring** of the same system: two mutation children (e.g. rattle) and two crossover
  children, i.e. realistic far-from-equilibrium starts (per-structure tests, phases H-I). Structures that need a long
  rearrangement exist in every GA run; keep one if you have it, it decides the optimiser comparison.
- The scheduler facts: node types (CPU model, cores, memory per node), per-user core cap, whether small jobs share
  nodes, whether whole-node jobs queue for long.

## 2. Measurement rules

1. **Metric = wall seconds per force evaluation** (job wall time / number of force evaluations) for per-step tests,
   and **wall minutes per structure** (sum of all stages) for relaxations. Program-internal timers are for diagnosis
   only; if you use one, reconcile it with the wall time on one run first (a gap > 5 % means the timer misses part
   of the work). Count force evaluations the same way for both codes: VASP `NSW=5` gives 5, CP2K `MAX_ITER 5` gives 6.
2. **Compare variants inside one allocation.** On shared nodes the neighbour load changes run times by 10-30 %.
   `speedtest.py ab` puts every variant into the same job and runs them in rotated order; use >= 2 jobs x 2 rounds
   (>= 4 samples per variant) and read the **within-job ratio** to the reference variant.
3. **Use deterministic metrics to separate cause from noise**: SCF iterations per force evaluation, the energy of the
   first step (same geometry), and the final energy and geometry after N optimiser steps do not depend on node load.
4. **Time single jobs, not packed runs.** If several runs share one allocation in parallel, verify the core binding
   of every run (`ps -o psr` per process); a mistake here overlaps all runs on the same cores and looks like load.
5. **Record the host of every run** and keep the raw outputs; a result without its host and input is not reusable.
6. **Plan with co-located numbers.** GOCIA runs many workers side by side; measure 10-12 identical jobs submitted at
   once and use their mean, not the best lone run.

## 3. Procedure

Each phase: what to vary, what to measure, when to accept. Phases B-G use the pre-relaxed structure with a 5-step
geometry optimisation; H-I use the offspring.

**A. Inventory and builds.** List node types and limits; build or locate each code. CP2K: build for the target CPU
(compiler flags for the micro-architecture), with ELPA if possible; confirm the build with a regression test. Check
the known build/runtime traps of section 5 before timing anything.

**B. Node type x cores.** Single jobs of 8, 16, 24, 48 (96) cores on every node type, both codes, pure MPI. Pick the
fastest node type, then the largest core count whose parallel efficiency relative to 8 cores stays >= 90 %; beyond
that, extra cores per structure cost throughput. Very large rank counts for Gamma-point slabs can also lose SCF
convergence; do not go beyond what you measured.

**C. Co-location.** 10-12 identical jobs of the chosen size at once; ratio of their mean to the lone run. Use the
co-located value in all planning.

**D. Parallel layout.**
- VASP: `NCORE` over the divisors of the rank count (sqrt(ranks) is only a starting guess), `NSIM` 4/8/16,
  `KPAR` only with several k-points, MPI x OpenMP (e.g. 8x2, 4x4).
- CP2K: MPI x OpenMP (1/2/4 threads).

**E. Libraries and builds.** Same inputs, different binaries or libraries: diagonalisation library (ScaLAPACK vs
ELPA, `&GLOBAL PREFERRED_DIAG_LIBRARY`), generic vs CPU-specific build, and on AMD CPUs with AVX-512 the MKL code path
(`tools/speedtest/mkl_amd_shim.c`). Accept only with identical energies (<= 0.01 meV).

**F. SCF settings, in two regimes.**
- CP2K: `EPS_SCF` (1e-5 / 1e-4 / 1e-3), Kerker damping `BETA` (1.0 / 1.5 / 2.5) in Broyden mixing, `ALPHA`,
  `NBROYDEN`, Pulay mixing, wavefunction extrapolation, `ADDED_MOS`, FFTW plan type. Diagonalisation + Fermi-Dirac
  smearing for metals (OT converges to minimiser-dependent integer-occupation states on metal slabs).
- VASP: `EDIFF`, `ALGO` (Fast / VeryFast / Normal), `NELMIN` only if your protocol allows it.
- Regime 1: the pre-relaxed structure, 5 steps (small steps near a minimum).
- Regime 2: one offspring, the first 10 steps at stage-1 settings (large steps). Looser SCF saves much less when the
  geometry changes a lot per step, and some settings behave differently there. A setting may be right for stage 1
  and wrong for stages 2-3.

**G. Cutoff and grid.** Single points of three structures at the production cutoff and at a high reference cutoff;
compare **relative** energies (absolute energies are not the criterion). Any change to the FFT grid options needs
the same check (section 5, item 7).

**H. Optimiser per stage.** Full multi-stage relaxations of all offspring (`speedtest.py pilot`), varying only the
optimiser of one stage at a time: CP2K BFGS (trust radius), LBFGS, CG; VASP CG (`IBRION 2`), quasi-Newton
(`IBRION 1`), VTST FIRE / LBFGS (`IBRION 3`, `IOPT 7` / `1`). Record per stage the steps and wall time, and check
that the stage-1 energy trace goes down (an optimiser that accepts uphill steps turns the pre-relaxation into a random
kick). The optimiser usually decides the cost per structure more than any per-step setting.

**I. Final pilot and throughput.** All chosen settings together on all offspring; mean minutes per structure;
throughput = floor(core cap / cores per job) x 60 / mean minutes (`speedtest.py throughput`). Compare codes on the
same structures, and check that both reached the same kind of minimum (section 4, last gate) before crediting a code
with a speed advantage.

## 4. Accuracy gates

| change | gate |
|---|---|
| build / library / MKL path | first and final energy identical (<= 0.01 meV) |
| SCF tolerance, mixing (stages 2-3) | same-geometry first-step dE < 1 meV; after 5 steps from a minimum dE_final < 1 meV and max displacement < 1e-3 A |
| SCF tolerance (stage 1, pre-relaxation) | same-geometry first-step dE < 10 meV |
| cutoff / FFT grid | relative energies within 5-10 meV of the high-cutoff reference |
| smearing (CP2K) | no "Add more MOs for proper smearing" warning |
| optimiser, full relaxations | converged in every stage; final structures compared on one PES (same code, grid and tight SCF single points); the same offspring can land in basins ~1 eV apart after small numerical changes, so do not rank settings by final energy on a handful of structures |
| code vs code | same site assignment of the adsorbates in the final structures (energies of different functionals are not comparable) |

## 5. Pitfalls met in practice

1. **Timer spans.** A parser that sums SCF-iteration times for one code and whole ionic steps for the other
   understated the first code by 8-18 %. Use wall time per force evaluation for both.
2. **Packed runs overlap.** Core slices computed in a subshell never advanced; every run landed on cores 0..N-1.
   Check `ps -o psr`.
3. **Intel MPI 2021 wrappers prepend the environment `FCFLAGS`/`CFLAGS`** to their own arguments. configure exports
   `FCFLAGS`, so autoconf's module test saw `-J` twice ("gfortran: Only one -J option allowed"). Fix: small wrapper
   scripts that `unset FCFLAGS FFLAGS CFLAGS CXXFLAGS` and exec the real `mpif90` / `mpicc`.
4. **Intel MPI 2021 ships `mpi.mod` only up to gfortran 11**; with gfortran >= 12 the wrapper falls back to the
   Intel-compiled module and `use mpi` fails. Fix: add `-I<impi>/include/gfortran/11.1.0` (those modules load in
   gfortran 13).
5. **Prebuilt CP2K (conda) may need a newer glibc** than the cluster provides; build from source.
6. **CP2K 2026.2 diagonalisation redistribution.** `cp_fm_redistribute_init` and `cp_fm_redistribute_work_finalize`
   reset the redistribution settings with the type constructor after assigning them, so every diagonalisation runs on
   a subset of ranks and `&FM_DIAG_SETTINGS` has no effect. Check with `&GLOBAL &FM_DIAG_SETTINGS
   PRINT_FM_REDISTRIBUTE` (all ranks should take part); fix by applying the reset before the assignments (measured:
   diagonalisation 235 -> 114 s at 16 ranks).
7. **`EXTENDED_FFT_LENGTHS` changes the grid.** One `&GLOBAL` keyword present in one input family and absent in the
   other shifted absolute energies by 0.75 eV (250 Ry) and 1.25 eV (350 Ry) and relative energies by 30-60 meV. When
   two runs of the "same" settings disagree at the first step, diff the full inputs before explaining anything.
8. **CP2K BFGS as a pre-relaxation optimiser.** With the default 0.25 A trust radius it accepted uphill steps on GA
   offspring (energy swinging by ~10 eV over 25 steps). LBFGS converged in 17-24 evaluations. For stages near a
   minimum, BFGS with its model Hessian was the fastest; LBFGS there spent extra line-search evaluations.
9. **VASP `IBRION 1` far from a minimum** diverged (energies jumping by hundreds of eV) and, after a FIRE stage,
   crawled for 86-342 steps. FIRE in every stage was the fastest robust choice.
10. **Loose SCF helps less at large steps**: SCF iterations per force evaluation dropped from ~16 to ~9 at stage-1
    settings, but only when the tolerance was relaxed to 1e-3; mixing changes alone gave ~10 %.
11. **Whole-node jobs on shared clusters** may wait for days while 16-core jobs start at once; the per-job worker
    launcher (`examples/gcga_cp2k/per-job/`) exists for that case.

## 6. Tools (`tools/speedtest/`)

```bash
cp tools/speedtest/cluster.example.json cluster.json          # scheduler header, modules, MPI launcher, binaries
# per-step A/B (phases B-G): a 5-step benchmark input from a stage template and the pre-relaxed structure
python tools/speedtest/speedtest.py bench-input --code cp2k --template examples/gcga_cp2k/cp2k-2.inp \
       --structure relaxed.vasp --out inputs/
python tools/speedtest/speedtest.py ab cluster.json --code cp2k --inputs inputs/ \
       --variants tools/speedtest/examples/variants_cp2k.json --out runs/scf --jobs 2 --rounds 2
for j in runs/scf/j*/job.sh; do sbatch $j; done
python tools/speedtest/speedtest.py parse-ab runs/scf --ref base --csv scf.csv
# relaxations (phases H-I): templates as shipped vs stage 1 switched back to BFGS
python tools/speedtest/speedtest.py pilot cluster.json --code cp2k --stages examples/gcga_cp2k \
       --structures kid1.vasp kid2.vasp kid3.vasp kid4.vasp --out runs/pilot_ref
python tools/speedtest/speedtest.py pilot cluster.json --code cp2k --stages examples/gcga_cp2k \
       --structures kid1.vasp kid2.vasp kid3.vasp kid4.vasp --out runs/pilot_bfgs1 \
       --edits tools/speedtest/examples/edits_stage1_bfgs.json
for j in runs/pilot_*/*/job.sh; do sbatch $j; done
python tools/speedtest/speedtest.py parse-pilot runs/pilot_* --csv pilots.csv
python tools/speedtest/speedtest.py throughput --minutes 24.3 --cap 1300 --cores-per-job 16
```

Variants are exact text edits of the main input (`edits`: each target must occur exactly once), extra lines
(`append`), and optional `ranks`, `omp`, `exe`, `env` per variant (another binary, `LD_PRELOAD`, ...).

## 7. Reference result

Hardware: 192-core AMD EPYC 9655 nodes (Zen 5, AVX-512), shared between jobs, 16 cores per structure, pure MPI.
System: Cu(100) 6x6x4 slab + 4 CO + 10 H (162 atoms, bottom two layers fixed, Gamma point, vacuum); VASP RPBE, CP2K
PBE (GTH / MOLOPT-SR). Three stages with force criteria 0.5 / 0.2 / 0.1 eV/A.

**CP2K** (templates `examples/gcga_cp2k/cp2k-1..3.inp`)

| setting | value |
|---|---|
| build | CPU-specific build with ELPA; `&GLOBAL PREFERRED_DIAG_LIBRARY ELPA`; MKL AVX-512 path on AMD |
| SCF | diagonalisation + Fermi-Dirac 300 K, Broyden ALPHA 0.4 / NBROYDEN 8 + Kerker BETA 1.5, ADDED_MOS = max(30, ceil(N_atoms/2)) |
| stage 1 | LBFGS, CUTOFF 250 / REL_CUTOFF 40, EPS_SCF 1e-3, MAX_ITER 25, writes the wavefunction |
| stages 2-3 | BFGS, 350 / 50 Ry, EPS_SCF 1e-4, restart from the previous stage's wavefunction |
| off | `EXTENDED_FFT_LENGTHS`, OT, OpenMP |

**VASP** (`examples/gcga_vasp_tuned/INCAR-1..3`)

| setting | value |
|---|---|
| parallel | 16 MPI ranks, `NCORE 8`, `NSIM 16`, MKL AVX-512 path on AMD; no OpenMP |
| SCF | `ALGO Fast`, `EDIFF 1e-4` in all stages |
| optimiser | FIRE in all stages (`IBRION 3`, `IOPT 7`, `POTIM 0`, `MAXMOVE 0.2`, `TIMESTEP 0.1`; VTST build) |
| restart | stage 2 writes WAVECAR/CHGCAR, stage 3 reads them (`ISTART 1`, `ICHARG 1`) |

**Measured** (per structure, four offspring; means)

| | before tuning | after tuning | structures/h at 1300 cores |
|---|---|---|---|
| CP2K | 95.8 min | 24.3 min | 200 |
| VASP | 155.6 min | 77.8 min | 62 |

Per-step levers (5-step benchmark, within-job ratios): CP2K ELPA 0.90x, ELPA + MKL AVX-512 path 0.80x, EPS_SCF 1e-4
0.64x (dE 0.06 meV), + Kerker 1.5 0.57x; CPU-specific compilation alone, fewer ADDED_MOS, FFTW MEASURE, REL_CUTOFF
40, the DGEMM grid backend, Pulay mixing and other extrapolations gave nothing or were slower. VASP NSIM 16 0.87x,
MKL AVX-512 path 0.90x, EDIFF 1e-4 0.90x (dE 0.41 meV); OpenMP 1.7-2.3x slower, ALGO VeryFast no gain.
Per-structure levers: CP2K stage-1 LBFGS 0.67x; VASP FIRE halved the SCF work of stage 2 relative to CG. The
hardest offspring (a ~1 eV downhill rearrangement) took CP2K 36 min and VASP 178 min, and both reached the same
adsorbate configuration.
