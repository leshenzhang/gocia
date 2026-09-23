# speedtest: speed tests of VASP / CP2K relaxations on any cluster

Method, measurement rules, accuracy gates and a worked result: `docs/SPEED_TUNING.md` (Chinese: `SPEED_TUNING.zh.md`).

| file | role |
|---|---|
| `speedtest.py` | `bench-input`, `ab`, `pilot`, `parse-ab`, `parse-pilot`, `throughput` (standard library; CP2K input writing uses ase + `gocia.utils.cp2k`) |
| `cluster.example.json` | the only cluster-specific file: scheduler header, setup lines, MPI launcher, binaries, POTCAR directory |
| `examples/variants_*.json` | A/B variants written against the templates `examples/gcga_cp2k/cp2k-2.inp` and `examples/gcga_vasp_tuned/INCAR-2` |
| `examples/edits_*.json` | stage edits for `pilot` (stage-1 BFGS instead of LBFGS; VASP stage-2 CG instead of FIRE) |
| `mkl_amd_shim.c` | MKL AVX-512 code path on AMD CPUs (`LD_PRELOAD`); adopt only with identical energies |

```bash
cp tools/speedtest/cluster.example.json cluster.json    # fill in
S=tools/speedtest/speedtest.py
# 1. per-step A/B on a pre-relaxed structure (5 optimiser steps, all variants in one allocation)
python $S bench-input --code cp2k --template examples/gcga_cp2k/cp2k-2.inp --structure relaxed.vasp --out inputs_cp2k
python $S ab cluster.json --code cp2k --inputs inputs_cp2k --variants tools/speedtest/examples/variants_cp2k.json --out runs/ab_cp2k
for j in runs/ab_cp2k/j*/job.sh; do sbatch $j; done
python $S parse-ab runs/ab_cp2k --ref base --csv ab_cp2k.csv
# 2. full relaxations of GA offspring (templates as shipped, and with one stage changed)
python $S pilot cluster.json --code cp2k --stages examples/gcga_cp2k --structures kid*.vasp --out runs/pilot_ref
python $S pilot cluster.json --code cp2k --stages examples/gcga_cp2k --structures kid*.vasp --out runs/pilot_bfgs1 \
       --edits tools/speedtest/examples/edits_stage1_bfgs.json
for j in runs/pilot_*/*/job.sh; do sbatch $j; done
python $S parse-pilot runs/pilot_ref runs/pilot_bfgs1 --csv pilots.csv
# 3. throughput for your core cap
python $S throughput --minutes 24.3 --cap 1300 --cores-per-job 16
```

A/B layout: `RUNDIR/j<k>/<variant>/` with outputs per round (`cp2k_r<R>.out` + `final_r<R>.xyz`, or `OUTCAR_r<R>` +
`CONTCAR_r<R>`) and `RUNDIR/j<k>/timings.txt`. Pilot layout: `RUNDIR/<structure>/stages.txt` with one
`STAGE i RC WALL STEPS CONVERGED` line per stage.
