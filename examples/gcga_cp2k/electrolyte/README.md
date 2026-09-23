# Constant-potential layer (optional)

Everything in this directory needs the separate implicit-electrolyte CP2K build
https://github.com/leshenzhang/cp2k-implicit-electrolyte (keywords `&SCCS DEBYE_LENGTH`, `METHOD ISODENSITY`,
`&SCF&GCE`). Stock CP2K runs route A (`../cp2k-1/2/3.inp`) without any of this.

| file | role |
|---|---|
| `cp2k-gce.inp` | route B: GEO_OPT at constant potential (`&GCE TARGET_POTENTIAL`), one run per structure |
| `cp2k-sc.inp` | route C: fixed-charge single point, 5 points per structure -> parabola |
| `sc-worker.py` | `--prepare` / packed run / `--collect` -> `sc.dat`, `parabola.dat` |

The parsing side (`gocia.utils.cp2k.gce_result`, `collect_surfChrg`, `get_parabola`, `add_result_SC`) lives in
the generic module and is a no-op without these outputs.
