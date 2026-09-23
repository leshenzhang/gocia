# GOCIA x CP2K (branch `cp2k-interface`)

GOCIA core untouched. `gocia/utils/cp2k.py` mirrors `gocia/utils/vasp.py`; every energy enters the
population through the existing `pop.add_aseResult(atoms, workdir)`. Route A needs only stock CP2K.
The constant-potential layer (`electrolyte/`) is optional and needs the separate implicit-electrolyte
build https://github.com/leshenzhang/cp2k-implicit-electrolyte — the model is not part of this repo.

## Files

| file | role | VASP counterpart |
|---|---|---|
| `input.py` | popSize / subsPot / chemPotDict / sampleSpecies / zLim / xyzLims / route / cp2k_cmd | input.py |
| `boxSample.py` | initial population by box sampling -> `init.db` | init/boxSample.py |
| `db2dirs.py` | `init.db` -> `s00001/POSCAR` ... | init/db2vasp.py |
| `init-worker.py` + `run-init.sbatch` | relax each initial structure (3 stages), write `done=1` rows into `gcga.db`; 4 per node | init-worker.py + slurm-vasp-init.sh |
| `dedup_db.py` | mark duplicate initial structures `alive=0` | init/collectVASP.py |
| `ga-worker.py` | offspring -> Hookean preopt -> CP2K (route A or B) -> db | ga-worker.py |
| `ga-bundle.py` + `run-bundle.sbatch` | ONE sbatch job, N nodes, 4 workers/node, `touch STOP` to end | ga-slurm.py + slurm-vasp.sh |
| `ibex/` | Ibex: `run-worker.sbatch` (one 16-core job = one worker, submit N) + `ga-loop.py` (race-free kid numbering) + `ga-prep.py` (run once) + `run-init-one.sbatch` + `env.sh` (install_z4, ELPA, MKL AVX-512 path) | ga-slurm.py |
| `electrolyte/sc-worker.py` | route C charge scan: `--prepare` / packed run / `--collect` -> `sc.dat`, `parabola.dat` (separate electrolyte build) | do_surfChrg_batch + get_parabola |
| `cp2k-1/2/3.inp` | 3-stage GEO_OPT, lean: diagonalisation + Fermi-Dirac 300 K (Multiwfn) + Kerker BETA 1.5; stage 1 LBFGS, 250/40 Ry, EPS_SCF 1e-3; stages 2-3 BFGS, 350/50 Ry, EPS_SCF 1e-4; force criteria 0.5/0.2/0.1 eV/A (= EDIFFG); wavefunction chained 1 -> 2 -> 3 (atomic guess if atoms were removed); `PREFERRED_DIAG_LIBRARY ELPA` (ScaLAPACK where ELPA is absent); only the final structure (FINAL xyz) is written | INCAR-1/2/3 |
| `electrolyte/cp2k-gce.inp` | route B: GEO_OPT at constant potential (`&SCCS DEBYE_LENGTH` + `&SCF&GCE`), separate electrolyte build | — |
| `electrolyte/cp2k-sc.inp` | route C: fixed-charge single point with electrolyte (5-point parabola), separate electrolyte build | INCAR-sc |

Template contract (all written from the ASE Atoms by `write_includes`): `@INCLUDE cell.inc`,
`coord.inc`, `constraint.inc` (FixAtoms -> `&FIXED_ATOMS LIST 1..8`), `kind.inc` (one `&KIND` per
element, `DZVP-MOLOPT-SR-GTH-qN` / `GTH-PBE-qN` from `KIND_Q`, override via `kind_overrides`);
`@SET PROJ`, `@SET CHG`, `@SET MULT` (UKS parity from the valence count and the charge; GCGA adds
and removes atoms so it flips per structure), `@SET NVAL`, `@SET ADDMOS` rewritten per run. Smearing
settings follow the group's Multiwfn 3.8 reference: Fermi-Dirac 300 K (raise if SCF stalls), Broyden
ALPHA 0.4 / NBROYDEN 8, ADDED_MOS = max(30, ceil(N_atoms/2)) ("n n" and UKS via `@IF ${MULT} == 2` for odd electron counts).

## Energies

| route | per structure | Omega fed to GOCIA | needs |
|---|---|---|---|
| A canonical/CHE | 3-stage GEO_OPT, CHARGE 0 | E − subsPot − Σµ'N (µ'_H carries pH, U by CHE) | stock CP2K |
| B constant potential | stages 1-2 neutral + GEO_OPT with `&GCE TARGET_POTENTIAL` | Ω_el(U)/nsides − subsPot − Σµ'N, Ω_el = E − E_F·N_e from the last `GCE\| trace` line | electrolyte build |
| C charge scan (a posteriori) | 5 single points, CHARGE −2..2 (or fractional `NELECTRON_EXCESS`) | fit Ω(U)=aU²+bU+c, `add_result_SC(pop, atoms, u_she, nsides=2)` | electrolyte build |

`W_f = −E_F` (bulk electrolyte potential 0 by construction, DEBYE_LENGTH > 0), `U = W_f − PHI_SHE`.
Route C validated on the project's stored Pt(111) 4x4x4 five-point outputs: Ω−Ω0, U0, R² identical.
The structure set is fixed across U; only Ω_i(U) is re-evaluated, so the 5-point scan is done once per
structure. For B/C the slab must be z-symmetrised and `subsPot` = bare-slab Ω_el(U)/nsides at the same U.
`&GCE` needs diagonalization + Fermi-Dirac smearing (no OT), hence the two SCF flavours in the templates.

## Run (Shaheen wangc0i, k10175)

```bash
# 0. code + python visible from compute nodes (/scratch, /home is not mounted there)
git clone -b cp2k-interface https://github.com/leshenzhang/gocia /scratch/wangc0i/zls/soft/gocia-cp2k
module load cray-python && python -m pip install --target /scratch/wangc0i/zls/soft/pylib ase natsort
# 1. initial population
python boxSample.py substrate.vasp 40 && python db2dirs.py init.db && sbatch run-init.sbatch
python dedup_db.py gcga.db
# 2. GCGA (job-name 01__gocia_cp2k, nworker = 4 x nodes, <=10 concurrent 01__ jobs)
sbatch run-bundle.sbatch
# 3. a posteriori constant potential on the ensemble (route C), inside a kid directory
python ../electrolyte/sc-worker.py --prepare -2 -1 0 1 2   # then run q_*/sc.inp packed, then
python ../electrolyte/sc-worker.py --collect -2 -1 0 1 2
```

Offline tests: `python tests/test_cp2k_interface.py` (7 tests, synthetic CP2K output).

## Run (Ibex, Turin nodes, never a whole node)

Measured on four GOCIA offspring of Cu(100) 6x6x4 + 4 CO + 10 H (162 atoms, vacuum), 16 cores per structure
(ai-reconstr-pre D-004..D-015): 24 min per 3-stage relaxation with these templates on `install_z4`, against 96 min
with the previous settings; about 200 relaxed structures/h per 1300-core account (81 workers).

```bash
# prerequisites: python with ase + natsort on the compute nodes (GOCIA_PY / GOCIA_PYLIB in ibex/env.sh), this repo
# at GOCIA_REPO, input.py cp2k_cmd = the Ibex line
python boxSample.py substrate.vasp 40 && python db2dirs.py init.db
for d in s0*/; do sbatch --chdir=$d ibex/run-init-one.sbatch; done
python dedup_db.py gcga.db && python ibex/ga-prep.py
for i in $(seq 25); do TOTCONF=2000 MINCONF=200 sbatch ibex/run-worker.sbatch; done     # 25 workers; touch STOP to end
```

Keep `EXTENDED_FFT_LENGTHS` off: at 350 Ry it moves relative energies by 30-60 meV (3-5 meV without, D-012).

