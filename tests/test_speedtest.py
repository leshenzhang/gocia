"""Tests for tools/speedtest/speedtest.py (run: python tests/test_speedtest.py)."""
import json
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'tools', 'speedtest'))
sys.path.insert(0, ROOT)
import speedtest  # noqa: E402

EX = os.path.join(ROOT, 'tools', 'speedtest', 'examples')
CFG = {'header': ['#!/bin/bash', '#SBATCH --job-name={name}', '#SBATCH --ntasks={ntasks}', '#SBATCH --time={walltime}',
                  '#SBATCH --output={output}'], 'setup': ['echo setup'], 'mpi': 'mpirun -np {np}',
       'cp2k': {'exe': '/opt/cp2k.psmp', 'env': ''}, 'vasp': {'exe': '/opt/vasp_gam', 'env': '', 'potcar_dir': 'POT'}}
POSCAR = """Cu H
1.0
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 12.0
Cu H
2 1
Selective dynamics
Direct
  0.0 0.0 0.0 F F F
  0.5 0.5 0.0 F F F
  0.25 0.25 0.12 T T T
"""


def _setup(tmp):
    os.chdir(tmp)
    json.dump(CFG, open('cfg.json', 'w'))
    open('s1.vasp', 'w').write(POSCAR)
    for e in ('Cu', 'H'):
        os.makedirs(f'POT/{e}')
        open(f'POT/{e}/POTCAR', 'w').write(f'PAW {e}\n')


def test_apply_edits():
    assert speedtest.apply_edits('a b a', [['b', 'c']], 'x') == 'a c a'
    for bad in (['z', 'y'], ['a', 'y']):
        try:
            speedtest.apply_edits('a b a', [bad], 'x')
        except SystemExit:
            continue
        raise AssertionError('edit target that is missing or repeated must stop')


def test_bench_input_vasp_and_variants(tmp):
    _setup(tmp)
    tpl = os.path.join(ROOT, 'examples', 'gcga_vasp_tuned', 'INCAR-2')
    speedtest.main(['bench-input', '--code', 'vasp', '--template', tpl, '--structure', 's1.vasp', '--out', 'bv', '--cfg', 'cfg.json'])
    inc = open('bv/INCAR').read()
    for line in ('NSW = 5', 'EDIFFG = -0.001', 'ISTART = 0', 'ICHARG = 2', 'LWAVE = .FALSE.', 'LCHARG = .FALSE.'):
        assert line in inc, line
    assert open('bv/POTCAR').read() == 'PAW Cu\nPAW H\n'
    for v, spec in json.load(open(os.path.join(EX, 'variants_vasp.json'))).items():
        speedtest.apply_edits(inc, spec.get('edits', []), v)          # every example edit matches exactly once
    speedtest.main(['ab', 'cfg.json', '--code', 'vasp', '--inputs', 'bv', '--variants', os.path.join(EX, 'variants_vasp.json'),
                    '--out', 'runs/ab', '--jobs', '2', '--rounds', '2'])
    j1, j2 = open('runs/ab/j1/job.sh').read(), open('runs/ab/j2/job.sh').read()
    assert 'mpirun -np 8' in j1 and 'OMP_NUM_THREADS=2' in j1 and 'for R in $(seq 1 2)' in j1
    first = lambda s: s.split('echo "')[1].split(' ')[0]
    assert first(j1) != first(j2), 'second job must run the variants in a different order'
    assert 'NSIM = 4' in open('runs/ab/j1/nsim4/INCAR').read() and 'NCORE' not in open('runs/ab/j1/omp2/INCAR').read()


def test_bench_input_cp2k_and_pilot(tmp):
    _setup(tmp)
    stages = os.path.join(ROOT, 'examples', 'gcga_cp2k')
    speedtest.main(['bench-input', '--code', 'cp2k', '--template', os.path.join(stages, 'cp2k-2.inp'), '--structure', 's1.vasp', '--out', 'bc'])
    inp = open('bc/cp2k.inp').read()
    assert 'MAX_ITER 5' in inp and 'MAX_FORCE 1.0E-5' in inp and 'SCF_GUESS RESTART' not in inp and 'WFN_RESTART' not in inp
    for f in ('cell.inc', 'coord.inc', 'constraint.inc', 'kind.inc'):
        assert os.path.isfile(os.path.join('bc', f)), f
    for v, spec in json.load(open(os.path.join(EX, 'variants_cp2k.json'))).items():
        speedtest.apply_edits(inp, spec.get('edits', []), v)
    speedtest.main(['pilot', 'cfg.json', '--code', 'cp2k', '--stages', stages, '--structures', 's1.vasp', '--out', 'pl',
                    '--edits', os.path.join(EX, 'edits_stage1_bfgs.json')])
    assert 'OPTIMIZER BFGS' in open('pl/s1/stage1.inp').read() and 'OPTIMIZER BFGS' in open('pl/s1/stage2.inp').read()
    job = open('pl/s1/job.sh').read()
    assert 'seq 1 3' in job and 'xyz2coord.py' in job


def test_pilot_vasp_with_edits(tmp):
    _setup(tmp)
    stages = os.path.join(ROOT, 'examples', 'gcga_vasp_tuned')
    open(os.path.join(tmp, 'KPOINTS'), 'w').write('Gamma\n0\nG\n1 1 1\n')
    sd = os.path.join(tmp, 'stages')
    os.makedirs(sd)
    for i in (1, 2, 3):
        open(os.path.join(sd, f'INCAR-{i}'), 'w').write(open(os.path.join(stages, f'INCAR-{i}')).read())
    open(os.path.join(sd, 'KPOINTS'), 'w').write('Gamma\n0\nG\n1 1 1\n')
    speedtest.main(['pilot', 'cfg.json', '--code', 'vasp', '--stages', sd, '--structures', 's1.vasp', '--out', 'pv',
                    '--edits', os.path.join(EX, 'edits_vasp_cg_stage2.json')])
    assert 'IBRION = 2' in open('pv/s1/INCAR-2').read() and 'IOPT = 7' in open('pv/s1/INCAR-1').read()
    assert os.path.isfile('pv/s1/KPOINTS') and open('pv/s1/POTCAR').read() == 'PAW Cu\nPAW H\n'


def _cp2k_out(e_ha, scf):
    lines = []
    for e, n in zip(e_ha, scf):
        lines += [f'  *** SCF run converged in    {n} steps ***',
                  f' ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]              {e:.12f}']
    return '\n'.join(lines) + '\n'


def test_parse_ab_and_pilot(tmp):
    os.chdir(tmp)
    for j, (wb, wv) in ((1, (100, 80)), (2, (120, 90))):
        for v, off in (('base', 0.0), ('fast', 1e-6)):
            d = f'ab/j{j}/{v}'
            os.makedirs(d)
            open(f'{d}/cp2k_r1.out', 'w').write(_cp2k_out([-10.0 + off, -10.1 + off], [20, 8]))
            open(f'{d}/final_r1.xyz', 'w').write('2\n\nH 0 0 0\nH 0 0 %.4f\n' % (0.74 + (0.0005 if v == 'fast' else 0)))
        open(f'ab/j{j}/timings.txt', 'w').write(f'base R=1 RC=0 WALL={wb}\nfast R=1 RC=0 WALL={wv}\n')
    rows = {r['variant']: r for r in speedtest.main(['parse-ab', 'ab', '--ref', 'base'])}
    assert rows['fast']['ratios'] == '0.80 0.75' and abs(rows['fast']['ratio_mean'] - 0.775) < 1e-9
    assert rows['fast']['dE_final_meV'] == '0.03' and rows['fast']['max_dr_A'] == '5.0e-04'
    assert rows['fast']['scf_per_force_eval'] == '20 8'
    os.makedirs('pil/k1')
    open('pil/k1/stage1.inp', 'w').write('x'); open('pil/k1/stage2.inp', 'w').write('x')
    open('pil/k1/stage2.out', 'w').write(_cp2k_out([-10.0], [9]))
    open('pil/k1/stages.txt', 'w').write('STAGE 1 RC=0 WALL=600 STEPS=22 CONVERGED=1\nSTAGE 2 RC=0 WALL=300 STEPS=10 CONVERGED=1\n')
    r = speedtest.main(['parse-pilot', 'pil'])[0]
    assert r['total_min'] == 15.0 and r['steps'] == '22+10' and r['complete'] == 1 and r['final_E_eV'] == '%.4f' % (-10.0 * speedtest.HA)


def test_vasp_readers(tmp):
    os.chdir(tmp)
    open('OUTCAR', 'w').write('LOOP: a\nLOOP: b\nLOOP+: x\n  free  energy   TOTEN  =      -5.00000000 eV\nLOOP: c\nLOOP+: y\n'
                              '  free  energy   TOTEN  =      -5.10000000 eV\n')
    scf, e = speedtest.read_vasp('OUTCAR')
    assert scf == [2, 1] and e == [-5.0, -5.1]
    open('C1', 'w').write(POSCAR)
    open('C2', 'w').write(POSCAR.replace('  0.25 0.25 0.12 T T T', '  0.99 0.25 0.12 T T T'))
    d = speedtest.max_disp(speedtest.read_contcar('C1'), speedtest.read_contcar('C2'))
    assert abs(d - 5.0 * 0.26) < 1e-9, d      # minimum image: 0.25 -> 0.99 is 0.26 of a 5 A cell, not 0.74


if __name__ == '__main__':
    import tempfile
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for t in tests:
        with tempfile.TemporaryDirectory() as tmp:
            here = os.getcwd()
            try:
                t(tmp) if t.__code__.co_argcount else t()
            finally:
                os.chdir(here)
        print('PASS', t.__name__)
    print(f'{len(tests)} tests passed')
