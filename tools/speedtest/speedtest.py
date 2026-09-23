#!/usr/bin/env python3
"""Cluster-agnostic speed tests for VASP / CP2K relaxations (method: docs/SPEED_TUNING.md).

  speedtest.py bench-input --code cp2k|vasp --template STAGE_FILE --structure POSCAR --out DIR [--steps 5] [--cfg CFG]
  speedtest.py ab          CFG --code cp2k|vasp --inputs DIR --variants VARIANTS.json --out RUNDIR
                               [--ranks 16] [--jobs 2] [--rounds 2] [--name TAG]
  speedtest.py pilot       CFG --code cp2k|vasp --stages DIR --structures POSCAR [POSCAR ...] --out RUNDIR
                               [--ranks 16] [--edits EDITS.json] [--name TAG] [--charge 0]
  speedtest.py parse-ab    RUNDIR [--ref base] [--csv OUT.csv]
  speedtest.py parse-pilot RUNDIR [RUNDIR ...] [--csv OUT.csv]
  speedtest.py throughput  --minutes M [M ...] --cap CORES --cores-per-job N [--structures 8000]

ab     one scheduler job runs every variant (in rotated order, `--rounds` times) on the same allocation, so the
       node's neighbour load is shared by all variants; ratios are then taken inside a job.
pilot  full multi-stage relaxations (stage templates cp2k-1..N.inp or INCAR-1..N) of real starting structures.
Everything is standard library except `pilot --code cp2k`, which uses ase + gocia.utils.cp2k to fill the
@INCLUDE files and @SET values of the GOCIA CP2K templates. Cluster specifics live only in the JSON config
(see cluster.example.json).
"""
import argparse
import csv
import glob
import json
import os
import re
import shutil
import statistics as st
import sys

HA = 27.211386245988
CP2K_E = re.compile(r'ENERGY\| Total FORCE_EVAL \( QS \) energy \[hartree\]\s+(-?\d+\.\d+)')
VASP_E = re.compile(r'free  energy   TOTEN  =\s+(-?\d+\.\d+)')


# ----------------------------------------------------------------------------------------------- job scripts
def header(cfg, name, ntasks, output):
    return '\n'.join(l.format(name=name, ntasks=ntasks, walltime=cfg.get('walltime', '04:00:00'), output=output)
                     for l in cfg['header']) + '\n' + '\n'.join(cfg.get('setup', [])) + '\n'


def launch(cfg, code, ranks, omp, exe=None, env=None):
    c = cfg[code]
    env = c.get('env', '') if env is None else env
    return f"env OMP_NUM_THREADS={omp} {env} {cfg['mpi'].format(np=ranks)} {exe or c['exe']}".replace('  ', ' ')


def apply_edits(text, edits, where):
    for old, new in edits:
        if text.count(old) != 1:
            sys.exit(f'{where}: edit target found {text.count(old)} times (need exactly 1): {old!r}')
        text = text.replace(old, new)
    return text


def write_job(path, text):
    with open(path, 'w') as f:
        f.write(text)
    os.chmod(path, 0o755)


# ------------------------------------------------------------------------------------------------------- ab
def cmd_ab(a):
    cfg = json.load(open(a.cfg))
    var = json.load(open(a.variants))
    main = 'cp2k.inp' if a.code == 'cp2k' else 'INCAR'
    if not os.path.isfile(os.path.join(a.inputs, main)):
        sys.exit(f'{a.inputs} has no {main}')
    names = list(var)
    tag = a.name or os.path.basename(os.path.normpath(a.out))
    for j in range(1, a.jobs + 1):
        jd = os.path.join(a.out, f'j{j}')
        if os.path.exists(jd):
            sys.exit(f'{jd} exists: refusing to overwrite')
        k = (j - 1) // 2 % len(names)
        order = names[k:] + names[:k]
        order = order if j % 2 else order[::-1]
        body = []
        for v in names:
            spec = var[v]
            vd = os.path.join(jd, v)
            shutil.copytree(a.inputs, vd)
            text = apply_edits(open(os.path.join(vd, main)).read(), spec.get('edits', []), f'{v}/{main}')
            if spec.get('append'):
                text = text.rstrip('\n') + '\n' + '\n'.join(spec['append']) + '\n'
            open(os.path.join(vd, main), 'w').write(text)
        for v in order:
            spec = var[v]
            run = launch(cfg, a.code, spec.get('ranks', a.ranks), spec.get('omp', 1), spec.get('exe'), spec.get('env'))
            if a.code == 'cp2k':
                body += [f'cd {v}; rm -f *.wfn *.wfn.bak-* *.Hessian; T0=$(date +%s)',
                         f'{run} -i cp2k.inp -o cp2k_r$R.out > run_r$R.log 2>&1',
                         f'echo "{v} R=$R RC=$? WALL=$(( $(date +%s) - T0 ))" >> ../timings.txt',
                         'f=$(ls -t *FINAL*.xyz *-pos-1.xyz 2>/dev/null | head -1); [ -n "$f" ] && cp "$f" final_r$R.xyz',
                         'cd ..']
            else:
                body += [f'cd {v}; rm -f WAVECAR CHGCAR CHG; T0=$(date +%s)',
                         f'{run} > run_r$R.log 2>&1',
                         f'echo "{v} R=$R RC=$? WALL=$(( $(date +%s) - T0 ))" >> ../timings.txt',
                         'cp OUTCAR OUTCAR_r$R; cp CONTCAR CONTCAR_r$R; rm -f WAVECAR CHGCAR CHG vasprun.xml', 'cd ..']
        ntasks = max(var[v].get('ranks', a.ranks) * var[v].get('omp', 1) for v in names)
        script = header(cfg, f'{tag}_j{j}', ntasks, os.path.abspath(os.path.join(jd, 'job.%j.out')))
        script += f'cd "{os.path.abspath(jd)}" || exit 1\nfor R in $(seq 1 {a.rounds}); do\n'
        script += '\n'.join('  ' + l for l in body) + '\ndone\n'
        write_job(os.path.join(jd, 'job.sh'), script)
        print(os.path.join(jd, 'job.sh'), ' '.join(order))


# ---------------------------------------------------------------------------------------------------- pilot
XYZ2COORD = '''import sys
L = open(sys.argv[1]).read().splitlines()
n = int(L[0].split()[0]); frame = L[-n:]
assert all(len(l.split()) == 4 for l in frame), "last xyz frame incomplete"
open(sys.argv[2], "w").write("".join("  %s %s %s %s\\n" % tuple(l.split()) for l in frame))
'''


def stage_files(code, d):
    pat = r'cp2k-(\d+)\.inp$' if code == 'cp2k' else r'INCAR-(\d+)$'
    got = sorted((int(re.search(pat, f).group(1)), f) for f in os.listdir(d) if re.search(pat, f))
    if not got or [i for i, _ in got] != list(range(1, len(got) + 1)):
        sys.exit(f'{d}: need stage templates numbered 1..N ({pat})')
    return got


def potcar_for(poscar, potcar_dir, out):
    L = open(poscar).read().splitlines()
    elems = L[5].split()
    with open(out, 'w') as f:
        for e in elems:
            f.write(open(os.path.join(potcar_dir, e, 'POTCAR')).read())


def cmd_bench_input(a):
    """5-step geometry-optimisation benchmark input from a stage template + structure (all N steps always run)."""
    if os.path.exists(a.out):
        sys.exit(f'{a.out} exists: refusing to overwrite')
    os.makedirs(a.out)
    n = a.steps
    if a.code == 'cp2k':
        from ase.io import read
        from gocia.utils import cp2k
        atoms = read(a.structure)
        cp2k.write_includes(atoms, a.out)
        out = os.path.join(a.out, 'cp2k.inp')
        cp2k.make_input(a.template, out, 'bench', charge=a.charge, sets=cp2k.default_sets(atoms, a.charge))
        t = open(out).read()
        t = re.sub(r'(?m)^(\s*MAX_ITER)\s+\d+', rf'\1 {n}', t)
        t = re.sub(r'(?m)^(\s*(?:MAX|RMS)_FORCE)\s+\S+', r'\1 1.0E-5', t)
        t = re.sub(r'(?m)^\s*(?:SCF_GUESS\s+RESTART|WFN_RESTART_FILE_NAME\s+\S+)\s*\n', '', t)
        open(out, 'w').write(t)
    else:
        cfg = json.load(open(a.cfg)) if a.cfg else {}
        t = open(a.template).read()
        for key, val in (('NSW', str(n)), ('EDIFFG', '-0.001'), ('ISTART', '0'), ('ICHARG', '2'),
                         ('LWAVE', '.FALSE.'), ('LCHARG', '.FALSE.')):
            t, k = re.subn(rf'(?m)^\s*{key}\s*=.*$', f'{key} = {val}', t)
            if not k:
                t = t.rstrip('\n') + f'\n{key} = {val}\n'
        open(os.path.join(a.out, 'INCAR'), 'w').write(t)
        shutil.copy(a.structure, os.path.join(a.out, 'POSCAR'))
        kp = os.path.join(os.path.dirname(os.path.abspath(a.template)), 'KPOINTS')
        if os.path.isfile(kp):
            shutil.copy(kp, a.out)
        if cfg.get('vasp', {}).get('potcar_dir'):
            potcar_for(a.structure, cfg['vasp']['potcar_dir'], os.path.join(a.out, 'POTCAR'))
    print('benchmark input in', a.out)


def cmd_pilot(a):
    cfg = json.load(open(a.cfg))
    stages = stage_files(a.code, a.stages)
    edits = json.load(open(a.edits)) if a.edits else {}
    tag = a.name or os.path.basename(os.path.normpath(a.out))
    run = launch(cfg, a.code, a.ranks, 1)
    for poscar in a.structures:
        d = os.path.join(a.out, os.path.splitext(os.path.basename(poscar))[0])
        if os.path.exists(d):
            sys.exit(f'{d} exists: refusing to overwrite')
        os.makedirs(d)
        if a.code == 'cp2k':
            from ase.io import read
            from gocia.utils import cp2k
            atoms = read(poscar)
            cp2k.write_includes(atoms, d)
            sets = cp2k.default_sets(atoms, a.charge)
            for i, f in stages:
                out = os.path.join(d, f'stage{i}.inp')
                cp2k.make_input(os.path.join(a.stages, f), out, f'stage{i}', charge=a.charge, sets=sets)
                ed = edits.get('all', []) + edits.get(str(i), [])
                text = apply_edits(open(out).read(), ed, out)      # read before reopening for writing (truncates)
                open(out, 'w').write(text)
            open(os.path.join(d, 'xyz2coord.py'), 'w').write(XYZ2COORD)
            body = [f'for i in $(seq 1 {len(stages)}); do', '  T0=$(date +%s)',
                    f'  {run} -i stage$i.inp -o stage$i.out > run-$i.log 2>&1; rc=$?',
                    '  n=$(grep -c "ENERGY| Total FORCE_EVAL" stage$i.out); c=$(grep -c "GEOMETRY OPTIMIZATION COMPLETED" stage$i.out)',
                    '  echo "STAGE $i RC=$rc WALL=$(( $(date +%s) - T0 )) STEPS=$n CONVERGED=$c" >> stages.txt',
                    '  f=$(ls -t stage$i-FINAL*.xyz stage$i-pos-1.xyz 2>/dev/null | head -1); [ -s "$f" ] && python3 xyz2coord.py "$f" coord.inc',
                    'done', 'rm -f *.wfn *.wfn.bak-*']
        else:
            shutil.copy(poscar, os.path.join(d, 'POSCAR'))
            shutil.copy(poscar, os.path.join(d, 'POSCAR.init'))
            for f in os.listdir(a.stages):
                if not re.match(r'INCAR-\d+$', f):
                    shutil.copy(os.path.join(a.stages, f), d)
            if not os.path.isfile(os.path.join(d, 'POTCAR')):
                potcar_for(poscar, cfg['vasp']['potcar_dir'], os.path.join(d, 'POTCAR'))
            for i, f in stages:
                ed = edits.get('all', []) + edits.get(str(i), [])
                text = apply_edits(open(os.path.join(a.stages, f)).read(), ed, f)
                open(os.path.join(d, f'INCAR-{i}'), 'w').write(text)
            body = [f'for i in $(seq 1 {len(stages)}); do', '  cp INCAR-$i INCAR; T0=$(date +%s)',
                    f'  {run} > run-$i.log 2>&1; rc=$?', '  cp OUTCAR OUTCAR-$i',
                    '  n=$(grep -c "LOOP+" OUTCAR); c=$(grep -c "reached required accuracy" OUTCAR)',
                    '  echo "STAGE $i RC=$rc WALL=$(( $(date +%s) - T0 )) STEPS=$n CONVERGED=$c" >> stages.txt',
                    '  [ -s CONTCAR ] && cp CONTCAR CONTCAR-$i && cp CONTCAR POSCAR', 'done',
                    'rm -f WAVECAR CHGCAR CHG vasprun.xml']
        script = header(cfg, f'{tag}_{os.path.basename(d)}', a.ranks, os.path.abspath(os.path.join(d, 'job.%j.out')))
        script += f'cd "{os.path.abspath(d)}" || exit 1\n' + '\n'.join(body) + '\n'
        write_job(os.path.join(d, 'job.sh'), script)
        print(os.path.join(d, 'job.sh'))


# ---------------------------------------------------------------------------------------------------- parse
def read_cp2k(out):
    t = open(out, errors='replace').read()
    scf = [int(x) for x in re.findall(r'SCF run converged in\s+(\d+)\s+steps', t)]
    scf += [-1] * len(re.findall(r'SCF run NOT converged', t))
    return scf, [float(x) * HA for x in CP2K_E.findall(t)]


def read_vasp(outcar):
    scf, n = [], 0
    with open(outcar, errors='replace') as f:
        text = f.read()
    for l in text.splitlines():
        if 'LOOP+:' in l:
            scf.append(n)
            n = 0
        elif 'LOOP:' in l:
            n += 1
    return scf, [float(x) for x in VASP_E.findall(text)]


def read_xyz(f):
    L = open(f).read().splitlines()
    n = int(L[0].split()[0])
    return None, [[float(x) for x in l.split()[1:4]] for l in L[-n:]]


def read_contcar(f):
    L = open(f).read().splitlines()
    s = float(L[1])
    A = [[s * float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]
    n = sum(int(x) for x in L[6].split())
    k = 8 if L[7].strip()[:1] in 'Ss' else 7
    frac = [[float(x) for x in L[k + 1 + i].split()[:3]] for i in range(n)]
    return A, frac


def max_disp(g1, g2):
    """Largest atomic displacement (A). Fractional input (CONTCAR) is wrapped with the minimum image."""
    A, p1 = g1
    _, p2 = g2
    m = 0.0
    for a, b in zip(p1, p2):
        d = [x - y for x, y in zip(a, b)]
        if A is not None:
            d = [x - round(x) for x in d]
            d = [sum(d[i] * A[i][j] for i in range(3)) for j in range(3)]
        m = max(m, sum(x * x for x in d) ** 0.5)
    return m


def cmd_parse_ab(a):
    W, runs = {}, {}
    for jd in sorted(glob.glob(os.path.join(a.rundir, 'j*'))):
        tf = os.path.join(jd, 'timings.txt')
        for l in open(tf) if os.path.isfile(tf) else []:
            m = re.match(r'(\S+)(?: R=(\d+))? RC=(\d+) WALL=(\d+)', l)
            if m and m.group(3) == '0':
                W[(jd, m.group(1), m.group(2) or '1')] = int(m.group(4))
        for key in [k for k in W if k[0] == jd]:
            _, v, r = key
            vd = os.path.join(jd, v)
            if os.path.isfile(os.path.join(vd, f'cp2k_r{r}.out')):
                scf, e = read_cp2k(os.path.join(vd, f'cp2k_r{r}.out'))
                g = os.path.join(vd, f'final_r{r}.xyz')
                geo = read_xyz(g) if os.path.isfile(g) else None
            elif os.path.isfile(os.path.join(vd, f'OUTCAR_r{r}')):
                scf, e = read_vasp(os.path.join(vd, f'OUTCAR_r{r}'))
                g = os.path.join(vd, f'CONTCAR_r{r}')
                geo = read_contcar(g) if os.path.isfile(g) else None
            else:
                continue
            runs[key] = dict(scf=scf, e=e, geo=geo)
    if not W:
        sys.exit(f'no finished runs under {a.rundir}')
    variants = sorted({k[1] for k in W}, key=lambda v: (v != a.ref, v))
    rows = []
    for v in variants:
        walls = [W[k] for k in sorted(W) if k[1] == v]
        ratios, de1, def_, dr, scf = [], [], [], [], None
        for k in sorted(W):
            if k[1] != v:
                continue
            ref = (k[0], a.ref, k[2])
            if ref in W:
                ratios.append(W[k] / W[ref])
            r, rr = runs.get(k), runs.get(ref)
            if r and scf is None:
                scf = r['scf']
            if r and rr and r['e'] and rr['e']:
                de1.append((r['e'][0] - rr['e'][0]) * 1000)
                def_.append((r['e'][-1] - rr['e'][-1]) * 1000)
            if r and rr and r['geo'] and rr['geo']:
                dr.append(max_disp(r['geo'], rr['geo']))
        f = lambda xs, fmt: fmt % max(xs, key=abs) if xs else '-'
        rows.append(dict(variant=v, samples=len(walls), wall_mean=round(st.mean(walls), 1),
                         ratio_mean=round(st.mean(ratios), 3) if ratios else '', ratios=' '.join('%.2f' % x for x in ratios),
                         scf_per_force_eval=' '.join(map(str, scf or [])), dE_first_meV=f(de1, '%.2f'),
                         dE_final_meV=f(def_, '%.2f'), max_dr_A=f(dr, '%.1e')))
    print('%-14s %3s %8s %7s %-24s %-28s %9s %9s %9s' % ('variant', 'n', 'wall_s', 'ratio', 'ratios (same job+round)',
                                                         'SCF per force evaluation', 'dE1 meV', 'dEf meV', 'max dr A'))
    for r in rows:
        print('%-14s %3d %8.1f %7s %-24s %-28s %9s %9s %9s' % (r['variant'], r['samples'], r['wall_mean'], r['ratio_mean'],
              r['ratios'][:24], r['scf_per_force_eval'][:28], r['dE_first_meV'], r['dE_final_meV'], r['max_dr_A']))
    if a.csv:
        with open(a.csv, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return rows


def cmd_parse_pilot(a):
    rows = []
    for root in a.rundirs:
        for st_file in sorted(glob.glob(os.path.join(root, '*', 'stages.txt'))):
            d = os.path.dirname(st_file)
            L = re.findall(r'STAGE (\d+) RC=(\d+) WALL=(\d+) (?:STEPS|FORCE_EVALS|IONIC)=(\d+) CONVERGED=(\d+)', open(st_file).read())
            nst = len(glob.glob(os.path.join(d, 'stage*.inp'))) or len(glob.glob(os.path.join(d, 'INCAR-*')))
            e = None
            if L and len(L) == nst:
                last = L[-1][0]
                if os.path.isfile(os.path.join(d, f'stage{last}.out')):
                    x = read_cp2k(os.path.join(d, f'stage{last}.out'))[1]
                elif os.path.isfile(os.path.join(d, f'OUTCAR-{last}')):
                    x = read_vasp(os.path.join(d, f'OUTCAR-{last}'))[1]
                else:
                    x = []
                e = x[-1] if x else None
            rows.append(dict(run=os.path.basename(os.path.normpath(root)), structure=os.path.basename(d),
                             steps='+'.join(s[3] for s in L), wall_s='+'.join(s[2] for s in L),
                             total_min=round(sum(int(s[2]) for s in L) / 60, 1), converged='+'.join(s[4] for s in L),
                             failed=sum(s[1] != '0' for s in L), final_E_eV='' if e is None else '%.4f' % e,
                             complete=int(bool(L) and len(L) == nst)))
    if not rows:
        sys.exit('no stages.txt found')
    for r in rows:
        print('%-16s %-12s %-14s %8.1f min  converged %-8s failed %d  E %s%s' % (
            r['run'], r['structure'], r['steps'], r['total_min'], r['converged'], r['failed'], r['final_E_eV'],
            '' if r['complete'] else '  (incomplete)'))
    for run in sorted({r['run'] for r in rows}):
        done = [r['total_min'] for r in rows if r['run'] == run and r['complete']]
        n = sum(r['run'] == run for r in rows)
        if done:
            print(f'{run}: mean {st.mean(done):.1f} min over {len(done)}/{n} complete structures')
    if a.csv:
        with open(a.csv, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return rows


def cmd_throughput(a):
    workers = a.cap // a.cores_per_job
    for m in a.minutes:
        per_h = workers * 60 / m
        print(f'{m:6.1f} min/structure x {workers} workers ({a.cores_per_job} cores each, cap {a.cap}): '
              f'{per_h:.0f} structures/h, {a.structures} structures in {a.structures / per_h / 24:.1f} days')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd', required=True)
    x = sp.add_parser('ab')
    x.add_argument('cfg'); x.add_argument('--code', choices=['cp2k', 'vasp'], required=True)
    x.add_argument('--inputs', required=True); x.add_argument('--variants', required=True); x.add_argument('--out', required=True)
    x.add_argument('--ranks', type=int, default=16); x.add_argument('--jobs', type=int, default=2)
    x.add_argument('--rounds', type=int, default=2); x.add_argument('--name')
    x = sp.add_parser('bench-input')
    x.add_argument('--code', choices=['cp2k', 'vasp'], required=True); x.add_argument('--template', required=True)
    x.add_argument('--structure', required=True); x.add_argument('--out', required=True)
    x.add_argument('--steps', type=int, default=5); x.add_argument('--cfg'); x.add_argument('--charge', type=int, default=0)
    x = sp.add_parser('pilot')
    x.add_argument('cfg'); x.add_argument('--code', choices=['cp2k', 'vasp'], required=True)
    x.add_argument('--stages', required=True); x.add_argument('--structures', nargs='+', required=True)
    x.add_argument('--out', required=True); x.add_argument('--ranks', type=int, default=16)
    x.add_argument('--edits'); x.add_argument('--name'); x.add_argument('--charge', type=int, default=0)
    x = sp.add_parser('parse-ab')
    x.add_argument('rundir'); x.add_argument('--ref', default='base'); x.add_argument('--csv')
    x = sp.add_parser('parse-pilot')
    x.add_argument('rundirs', nargs='+'); x.add_argument('--csv')
    x = sp.add_parser('throughput')
    x.add_argument('--minutes', type=float, nargs='+', required=True); x.add_argument('--cap', type=int, required=True)
    x.add_argument('--cores-per-job', type=int, required=True); x.add_argument('--structures', type=int, default=8000)
    a = p.parse_args(argv)
    return {'ab': cmd_ab, 'bench-input': cmd_bench_input, 'pilot': cmd_pilot, 'parse-ab': cmd_parse_ab,
            'parse-pilot': cmd_parse_pilot, 'throughput': cmd_throughput}[a.cmd](a)


if __name__ == '__main__':
    main()
