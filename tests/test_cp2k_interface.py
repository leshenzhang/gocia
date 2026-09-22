"""Offline tests for gocia.utils.cp2k (synthetic CP2K output; no CP2K binary needed).
Run:  python tests/test_cp2k_interface.py
"""
import os
import sys
import tempfile
import numpy as np
from ase.build import fcc100, add_adsorbate
from ase.constraints import FixAtoms
from ase.db import connect
from ase.io import write, read

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gocia.utils import cp2k
from gocia.ga.popGrandCanon import PopulationGrandCanonical

HA = cp2k.HA2EV

FAKE_OUT = """ SCF run converged in    18 steps
 Fermi energy:                                        -0.18000000
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:              -100.5000000000
 SCF run converged in    12 steps
 Fermi energy:                                        -0.20000000
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:              -101.0000000000
 GCE| trace  step 2  N_e 0.2500 e  W_f 5.4400 eV  U 1.0100 V vs SHE  U_target 1.0000 V  dU 0.0100 V
 ***                    GEOMETRY OPTIMIZATION COMPLETED                      ***
  PROGRAM ENDED AT                               2026-09-22 15:00:00.000
"""


def slab():
    a = fcc100('Pt', size=(2, 2, 3), vacuum=6.0)
    a.set_constraint(FixAtoms(indices=[i for i in range(len(a)) if a[i].tag >= 2]))
    return a


def test_ranges():
    assert cp2k._ranges_1based([0, 1, 2, 5, 6, 9]) == '1..3 6..7 10'
    assert cp2k._ranges_1based([]) == ''


def test_includes_and_input(tmp):
    a = slab()
    cp2k.write_includes(a, tmp)
    cons = open(os.path.join(tmp, 'constraint.inc')).read()
    assert '&FIXED_ATOMS' in cons and 'LIST 1..8' in cons, cons
    assert len(open(os.path.join(tmp, 'coord.inc')).read().splitlines()) == 12
    kind = open(os.path.join(tmp, 'kind.inc')).read()
    assert '&KIND Pt' in kind and 'DZVP-MOLOPT-SR-GTH-q18' in kind and 'GTH-PBE-q18' in kind
    assert cp2k.valence_electrons(a) == 12 * 18 and cp2k.multiplicity(a) == 1 and cp2k.multiplicity(a, 1) == 2
    from ase import Atoms
    ref = Atoms('Mo24S46O92H185')                       # the group's reference input: UKS MULTIPLICITY 2
    assert cp2k.valence_electrons(ref) == 1349 and cp2k.multiplicity(ref) == 2
    assert cp2k.kind_block(['Mo', 'S']).count('&KIND') == 2
    tpl = os.path.join(tmp, 'tpl.inp')
    open(tpl, 'w').write('@SET PROJ x\n@SET CHG 0\n&GLOBAL\n PROJECT ${PROJ}\n&END GLOBAL\n')
    out = os.path.join(tmp, 'run.inp')
    cp2k.make_input(tpl, out, 'kid1', charge=-1, sets={'USHE': -0.3})
    t = open(out).read()
    assert '@SET PROJ kid1' in t and '@SET CHG -1' in t and t.startswith('@SET USHE -0.3')


def test_parse(tmp):
    fn = os.path.join(tmp, 'x.out')
    open(fn, 'w').write(FAKE_OUT)
    d = cp2k.parse_output(fn)
    assert abs(d['energy_eV'] + 101.0 * HA) < 1e-6
    assert abs(d['fermi_eV'] + 0.2 * HA) < 1e-6
    assert d['n_force_eval'] == 2 and d['ended'] and d['geo_converged']
    assert d['gce'] == dict(step=2, N_e=0.25, W_f=5.44, U_SHE=1.01)
    assert cp2k.is_success(fn, require_geo=True)
    assert cp2k.parse_output(os.path.join(tmp, 'missing.out')) is None


def test_final_geometry(tmp):
    a = slab()
    proj = 'stage1'
    b = a.copy()
    b.positions[-1] += [0.0, 0.0, 0.3]
    with open(os.path.join(tmp, f'{proj}-pos-1.xyz'), 'w') as f:
        for frame in (a, b):
            f.write(f'{len(frame)}\n i =        1, E =      -100.5\n')
            for s, p in zip(frame.get_chemical_symbols(), frame.positions):
                f.write(f' {s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n')
    new = cp2k.read_final_geometry(proj, a, dirname=tmp)
    assert abs(new.positions[-1, 2] - a.positions[-1, 2] - 0.3) < 1e-5
    assert len(new.constraints) == 1 and list(new.constraints[0].get_indices()) == list(a.constraints[0].get_indices())
    assert np.allclose(new.cell, a.cell)


def test_grand_potential_and_parabola(tmp):
    # Omega = E - E_F*N_e ; identical to gocia.utils.vasp.pb_calc with FERMI_SHIFT=0
    assert abs(cp2k.grand_potential_el(-10.0, -4.5, 1.0) - (-10.0 + 4.5)) < 1e-12
    U = np.array([-0.6, 0.1, 0.8, 1.5, 2.1])
    a0, b0, c0 = -0.74, 1.16, -0.45
    with open(os.path.join(tmp, 'sc.dat'), 'w') as f:
        for u in U:
            f.write(f'{0}\t{u}\t{a0*u*u+b0*u+c0}\n')
    a, b, c, r2 = cp2k.get_parabola(tmp)
    assert max(abs(a - a0), abs(b - b0), abs(c - c0)) < 1e-9 and r2 > 0.999999
    assert abs(cp2k.parabola_energy(0.5, os.path.join(tmp, 'parabola.dat')) - (a0*0.25+b0*0.5+c0)) < 1e-9


def test_add_result_SC(tmp):
    os.chdir(tmp)
    subs = slab()
    write('substrate.vasp', subs, format='vasp')
    init = subs.copy(); add_adsorbate(init, 'H', 1.0, 'hollow')
    init.set_constraint(subs.constraints[0])
    with connect('gcga.db') as db:
        db.write(init, name='s0001', mag=0, eV=-50.0, grandPot=-50.0 - (-40.0) - (-3.4),
                 mated=0, done=1, alive=1, label='0 0 init')
    pop = PopulationGrandCanonical(gadb='gcga.db', substrate='substrate.vasp', popSize=5,
                                   subsPot=-40.0, chemPotDict={'H': -3.4}, zLim=[8, 14])
    os.makedirs('000001', exist_ok=True)
    kid = subs.copy(); add_adsorbate(kid, 'H', 1.0, 'ontop'); kid.set_constraint(subs.constraints[0])
    open('000001/label', 'w').write('1 1 grow')
    open('000001/parabola.dat', 'w').write('-0.6\t1.0\t-110.0')  # Omega_el(U) for the 2-sided slab
    from ase.calculators.singlepoint import SinglePointCalculator
    kid.calc = SinglePointCalculator(kid, energy=-49.0)
    g = cp2k.add_result_SC(pop, kid, u_she=-0.5, workdir='000001', nsides=2)
    ene_sc = -0.6*0.25 + 1.0*(-0.5) - 110.0
    expect = ene_sc/2 - (-40.0) - (-3.4)
    assert abs(g - expect) < 1e-9, (g, expect)
    row = pop.gadb.get(id=2)
    assert row.done == 1 and row.alive == 1 and abs(row.sc_U + 0.5) < 1e-12 and abs(row.a + 0.6) < 1e-12
    assert open('gmid').read().strip() == '2'     # lower Omega than the init row -> new GM
    # gce_result path
    open('000001/gce.out', 'w').write(FAKE_OUT)
    res, info = cp2k.gce_result('000001/gce.out', kid, nsides=2)
    omega = (-101.0*HA) - (-5.44)*0.25
    assert abs(res.get_potential_energy() - omega/2) < 1e-6 and res.info['sc_U'] == 1.01


def test_prepare_collect(tmp):
    os.chdir(tmp)
    write('POSCAR', slab(), format='vasp')
    open('tpl.inp', 'w').write('@SET PROJ x\n@SET CHG 0\n&GLOBAL\n PROJECT ${PROJ}\n&END GLOBAL\n')
    dirs = cp2k.prepare_surfChrg([-1, 0, 1], template='tpl.inp')
    assert dirs == ['q_-1', 'q_+0', 'q_+1'] and '@SET CHG -1' in open('q_-1/sc.inp').read()
    assert '@SET MULT 1' in open('q_+0/sc.inp').read() and '@SET MULT 2' in open('q_+1/sc.inp').read()
    assert os.path.isfile('q_+1/coord.inc')
    for q, ef, e in ((-1, -0.17, -100.2), (0, -0.19, -100.0), (1, -0.21, -100.3)):
        open(f'{cp2k._sc_dirname(q, False)}/sc.out', 'w').write(
            f' Fermi energy:   {ef}\n ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]   {e}\n PROGRAM ENDED AT x\n')
    rows = cp2k.collect_surfChrg([-1, 0, 1], phi_she=4.43)
    assert len(rows) == 3 and abs(rows[0][0] - 1) < 1e-12          # dN = -q
    assert abs(rows[1][1] - (0.19 * HA - 4.43)) < 1e-9              # U = -E_F - phi_SHE
    assert abs(rows[0][2] - ((-100.2 * HA) - (-0.17 * HA) * 1.0)) < 1e-6
    assert len(open('sc.dat').read().splitlines()) == 3
    fr = cp2k.prepare_surfChrg([-0.5], template='tpl.inp', fractional=True)
    assert fr == ['q_-0.50'] and '@SET NEX 0.5' in open('q_-0.50/sc.inp').read()


if __name__ == '__main__':
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
