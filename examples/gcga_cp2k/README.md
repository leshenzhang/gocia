# GOCIA x CP2K (branch `cp2k-interface`)

GOCIA core untouched. New module `gocia/utils/cp2k.py` mirrors `gocia/utils/vasp.py`; results
enter the population through the existing `pop.add_aseResult(atoms, workdir)`.

## Files

| file | role | VASP counterpart |
|---|---|---|
| `cp2k-1.inp` `cp2k-2.inp` `cp2k-3.inp` | 3-stage GEO_OPT templates (loose -> production) | INCAR-1/2/3 |
| `cp2k-gce.inp` | route B: GEO_OPT at constant potential (`&SCCS DEBYE_LENGTH` + `&SCF&GCE`) | — (new) |
| `cp2k-sc.inp` | route C: fixed-charge single point with electrolyte | INCAR-sc |
| `input.py` | popSize / subsPot / chemPotDict / zLim / route / cp2k_cmd | input.py |
| `init-worker.py` | relax one initial structure, write `done=1` row | init-worker.py |
| `ga-worker.py` | offspring -> Hookean preopt -> CP2K -> db | ga-worker.py |
| `ga-bundle.py` + `run-bundle.sbatch` | ONE sbatch job, N nodes, 4 workers/node | ga-slurm.py + slurm-vasp.sh |
| `sc-worker.py` | charge scan -> `sc.dat` -> `parabola.dat` | do_surfChrg_batch + get_parabola |

Template contract: `@SET PROJ`, `@SET CHG` rewritten per run; `@INCLUDE cell.inc / coord.inc /
constraint.inc` written from the ASE Atoms (FixAtoms -> `&FIXED_ATOMS LIST`). Add `&KIND` blocks
for your elements.

## Energies

| route | per structure | Omega fed to GOCIA | needs |
|---|---|---|---|
| A canonical/CHE | 1 GEO_OPT, CHARGE 0 | E − subsPot − Σµ'N (µ'_H carries pH, U by CHE) | stock CP2K |
| B constant potential | stages 1-2 neutral + 1 GEO_OPT with `&GCE TARGET_POTENTIAL` | Ω_el(U)/nsides − subsPot − Σµ'N, Ω_el = E − E_F·N_e from the last `GCE\| trace` line | ai-cp2k-lpb build |
| C charge scan (a posteriori) | 5 single points, CHARGE −2..2 | fit Ω(U)=aU²+bU+c, `add_result_SC(pop, atoms, u_she, nsides=2)` | ai-cp2k-lpb build (fractional: `NELECTRON_EXCESS`) |

`W_f = −E_F` because the bulk electrolyte potential is 0 by construction (DEBYE_LENGTH > 0);
`U = W_f − PHI_SHE`. For B/C the slab must be z-symmetrised and `subsPot` must be the bare
slab's Ω_el(U)/nsides at the same U.

## Run (Shaheen wangc0i, k10175)

```bash
# 0. gocia + this branch visible from compute nodes (/scratch, not /home)
git clone -b cp2k-interface <fork> /scratch/wangc0i/zls/soft/gocia-cp2k
# 1. sample + relax the initial population (boxSample.py from examples/init, then per s0*:)
mkdir s0001 && cp s0001.vasp s0001/POSCAR && (cd s0001 && python ../init-worker.py)
# 2. GCGA, bundled (nworker = 4 x nodes, job-name 01__gocia_cp2k, <=10 concurrent 01__ jobs)
sbatch run-bundle.sbatch          # touch STOP to end
# 3. a posteriori constant potential on the ensemble (route C)
(cd 000042 && python ../sc-worker.py -2 -1 0 1 2)
```

Offline tests: `python tests/test_cp2k_interface.py` (6 tests, synthetic CP2K output).
