"""
CP2K interface for GOCIA — mirrors gocia.utils.vasp.

Design rule: GOCIA core (Interface / Population* / GA) is untouched. This
module only (i) writes CP2K inputs from an ASE Atoms + a template, (ii) runs
cp2k, (iii) parses the output back into an Atoms carrying an energy, so the
existing population API `pop.add_aseResult(atoms, workdir)` can be used.

Three ways to get an energy per structure:
  A. canonical / CHE          : GEO_OPT at CHARGE 0          -> E
  B. constant potential (GCE) : GEO_OPT with &SCF&GCE        -> Omega_el(U) = E - E_F*N_e
  C. fixed-charge scan        : N single points, CHARGE q    -> sc.dat -> parabola.dat -> Omega_el(U)
B and C need the implicit electrolyte build (DEBYE_LENGTH > 0): with it the
bulk-electrolyte potential is 0 by construction, so W_f = -E_F and U = W_f - PHI_SHE.

Template convention (see examples/gcga_cp2k/cp2k-*.inp):
  @SET PROJ <name>      -> replaced per run
  @SET CHG  <int>       -> replaced per run (fixed-charge)
  @INCLUDE cell.inc     inside &SUBSYS&CELL
  @INCLUDE coord.inc    inside &SUBSYS&COORD
  @INCLUDE constraint.inc inside &MOTION&CONSTRAINT   (FixAtoms -> &FIXED_ATOMS)
Any other '@SET KEY value' line can be overridden through make_input(sets={...}).
"""
import os
import re
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.calculators.singlepoint import SinglePointCalculator
from gocia.interface import Interface
from gocia.geom import get_fragments, del_freeMol, is_bonded, detect_bond_between_adsFrag
from gocia.geom.frag import read_frag, update_frag_del

HA2EV = 27.211386245988


# ---------------------------------------------------------------- input side
def _fixed_indices(atoms):
    idx = []
    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            idx += [int(i) for i in c.get_indices()]
    return sorted(set(idx))


def _ranges_1based(idx):
    """[0,1,2,5,6] -> '1..3 6..7' (CP2K LIST syntax, 1-based)."""
    if not idx:
        return ''
    out, start, prev = [], idx[0], idx[0]
    for i in idx[1:]:
        if i == prev + 1:
            prev = i
            continue
        out.append(f'{start+1}..{prev+1}' if prev > start else f'{start+1}')
        start = prev = i
    out.append(f'{start+1}..{prev+1}' if prev > start else f'{start+1}')
    return ' '.join(out)


def write_includes(atoms, dirname='.'):
    """Write cell.inc / coord.inc / constraint.inc for the template @INCLUDEs."""
    cell = np.array(atoms.get_cell())
    with open(os.path.join(dirname, 'cell.inc'), 'w') as f:
        for lab, v in zip('ABC', cell):
            f.write(f'  {lab} {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n')
    with open(os.path.join(dirname, 'coord.inc'), 'w') as f:
        for s, p in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            f.write(f'  {s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n')
    fixed = _fixed_indices(atoms)
    with open(os.path.join(dirname, 'constraint.inc'), 'w') as f:
        if fixed:
            f.write(f'  &FIXED_ATOMS\n    LIST {_ranges_1based(fixed)}\n  &END FIXED_ATOMS\n')


def make_input(template, out, project, charge=None, sets=None):
    """Copy template -> out, rewriting '@SET PROJ', '@SET CHG' and any keys in `sets`."""
    sets = dict(sets or {})
    sets['PROJ'] = project
    if charge is not None:
        sets['CHG'] = charge
    lines = open(template).read().splitlines()
    seen = set()
    for i, l in enumerate(lines):
        m = re.match(r'^\s*@SET\s+(\S+)\s+(.*)$', l)
        if m and m.group(1) in sets:
            lines[i] = f'@SET {m.group(1)} {sets[m.group(1)]}'
            seen.add(m.group(1))
    missing = [k for k in sets if k not in seen]
    header = [f'@SET {k} {sets[k]}' for k in missing]
    with open(out, 'w') as f:
        f.write('\n'.join(header + lines) + '\n')
    return out


# --------------------------------------------------------------- output side
def parse_output(out):
    """Parse a CP2K output. Energies in eV. Returns None if the file lacks a FORCE_EVAL energy."""
    if not os.path.isfile(out):
        return None
    t = open(out, errors='replace').read()
    e = re.findall(r'ENERGY\| Total FORCE_EVAL[^\n]*?\s(-?\d+\.\d+)', t)
    if not e:
        return None
    ef = re.findall(r'Fermi energy:\s+(-?\d+\.\d+)', t)
    d = dict(
        energy_eV=float(e[-1]) * HA2EV,
        fermi_eV=float(ef[-1]) * HA2EV if ef else None,
        n_force_eval=len(e),
        ended='PROGRAM ENDED AT' in t,
        geo_converged='GEOMETRY OPTIMIZATION COMPLETED' in t,
        geo_max_steps='MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED' in t,
        scf_not_converged=('SCF run NOT converged' in t),
        gce=None,
    )
    # constant-potential trace (ai-cp2k-lpb GCE build): one line per ionic step
    g = re.findall(r'GCE\| trace\s+step\s+(\d+)\s+N_e\s+(-?\d+\.\d+)\s+e\s+W_f\s+(-?\d+\.\d+)\s+eV\s+U\s+(-?\d+\.\d+)\s+V', t)
    if g:
        step, ne, wf, u = g[-1]
        d['gce'] = dict(step=int(step), N_e=float(ne), W_f=float(wf), U_SHE=float(u))
    return d


def is_success(out, require_geo=False):
    d = parse_output(out)
    if d is None or not d['ended']:
        return False
    if require_geo and not (d['geo_converged'] or d['geo_max_steps']):
        return False
    return True


def grand_potential_el(energy_eV, fermi_eV, n_excess):
    """Omega_el = E - E_F * N_e (Goodpaster form; phi_bulk = 0 with the implicit electrolyte).
    Identical to gocia.utils.vasp.pb_calc's  G = E + dN * W_f  since W_f = -E_F here."""
    return energy_eV - fermi_eV * n_excess


def read_final_geometry(project, ref_atoms, dirname='.'):
    """Last frame of <project>-pos-1.xyz with the cell/pbc/constraints of ref_atoms.
    Falls back to ref_atoms (single point runs write no trajectory)."""
    fn = os.path.join(dirname, f'{project}-pos-1.xyz')
    new = ref_atoms.copy()
    if os.path.isfile(fn):
        fr = read(fn, index=-1)
        if len(fr) != len(ref_atoms):
            raise RuntimeError(f'{fn}: {len(fr)} atoms, expected {len(ref_atoms)}')
        new.set_positions(fr.get_positions())
    return new


def result_atoms(atoms, energy_eV):
    a = atoms.copy()
    a.calc = SinglePointCalculator(a, energy=energy_eV)
    return a


def run_cp2k(cp2k_cmd, inp, out):
    """cp2k_cmd already contains the launcher, e.g. 'srun -n 24 -c 2 /path/cp2k.psmp'."""
    rc = os.system(f'{cp2k_cmd} -i {inp} -o {out}')
    return rc


# ------------------------------------------------- multi-stage local optimisation
def do_multiStep_opt(step=3, cp2k_cmd='', poscar='POSCAR', template='../cp2k-%i.inp',
                     charge=0, zLim=None, substrate='../substrate.vasp', chkMol=False,
                     fn_frag='fragments', list_keep=[0], has_fragList=False,
                     has_fragSurfBond=False, check_rxn_frags=False, rmAtomsNotInBond=[],
                     sets=None, require_geo=False, clean=True):
    """Mirror of gocia.utils.vasp.do_multiStep_opt for CP2K.

    template: python %-format with the stage number -> '../cp2k-1.inp', ...
    Returns an Atoms with a SinglePointCalculator energy (eV) of the last stage,
    or None on failure (a FAIL file is also touched, like the VASP path).
    """
    if read_frag(fn=fn_frag) is not None:
        has_fragList = True
    counter = 1
    energy = None
    while counter <= step:
        print(f'Optimization step: {counter}')
        atoms = read(poscar)
        write_includes(atoms)
        proj = f'stage{counter}'
        make_input(template % counter, f'{proj}.inp', proj, charge=charge, sets=sets)
        run_cp2k(cp2k_cmd, f'{proj}.inp', f'{proj}.out')
        if not is_success(f'{proj}.out', require_geo=require_geo):
            os.system('touch FAIL')
            return None
        res = parse_output(f'{proj}.out')
        energy = res['energy_eV']
        final = read_final_geometry(proj, atoms)
        write(f'out-{counter}.vasp', final, format='vasp')
        write(poscar, final, format='vasp')
        atom_tmp = final.copy()
        counter += 1

        # --- geometry / connectivity checks, same order as the VASP path ---
        if zLim is not None:
            surf = Interface(read(poscar), substrate, zLim=zLim)
            if surf.has_outsideBox():
                if has_fragList:
                    surf.del_outsideBox_frag(fn_frag)
                else:
                    surf.del_outsideBox()
                surf.write(poscar)

        if has_fragList:
            _remove_broken_frags(poscar, fn_frag)

        if has_fragList and has_fragSurfBond:
            my_fragList = read_frag(fn=fn_frag)
            struct = read(poscar)
            list1 = list(range(len(read(substrate))))
            list_del = []
            for i in range(len(my_fragList)):
                if not is_bonded(struct, list1, my_fragList[i]):
                    list_del += [j for j in my_fragList[i] if j not in list_del]
            if list_del:
                list_del.sort()
                print('Remove associated fragments containing:', list_del)
                update_frag_del(list_del, fn=fn_frag)
                del struct[list_del]
                write(poscar, struct, format='vasp')

        if has_fragList and check_rxn_frags:
            my_fragList = read_frag(fn=fn_frag)
            struct = read(poscar)
            ads_bonds = detect_bond_between_adsFrag(struct, my_fragList)
            if ads_bonds:
                import random
                list_del = []
                for b in ads_bonds:
                    i = random.choice(b[0:2])
                    list_del += [j for j in my_fragList[i] if j not in list_del]
                if list_del:
                    list_del.sort()
                    update_frag_del(list_del, fn=fn_frag)
                    del struct[list_del]
                    write(poscar, struct, format='vasp')

        if chkMol:
            geom_tmp, list_del = del_freeMol(read(poscar), list_keep=list_keep)
            write(poscar, geom_tmp, format='vasp')
            if has_fragList and list_del:
                list_del.sort()
                update_frag_del(list_del, fn=fn_frag)

        if rmAtomsNotInBond:
            struct = read(poscar)
            list_del = []
            for p in rmAtomsNotInBond:
                list_a0 = [a.index for a in struct if a.symbol == p[0]]
                list_a1 = [a.index for a in struct if a.symbol == p[1]]
                for a0 in list_a0:
                    if not is_bonded(struct, [a0], list_a1):
                        list_del.append(a0)
            if list_del:
                list_del.sort()
                print(f'Atoms {list_del} are removed because not in required bonds.')
                del struct[list_del]
                write(poscar, struct, format='vasp')

        if has_fragList:
            _remove_broken_frags(poscar, fn_frag)

        natoms_removed = len(atom_tmp) - len(read(poscar))
        if natoms_removed > 0 and counter > step:
            print(f'Redo the last opt step due to removal of {natoms_removed} atoms')
            counter -= 1
            continue

    if clean:
        os.system('rm -f *.wfn *.wfn.bak-* *.cube *-RESTART.kp')
    return result_atoms(read(poscar), energy)


def _remove_broken_frags(poscar, fn_frag):
    my_fragList = read_frag(fn=fn_frag)
    struct = read(poscar)
    my_fragAtoms = [struct[f] for f in my_fragList]
    list_del = []
    for i in range(len(my_fragList)):
        if len(get_fragments(my_fragAtoms[i])) != 1:
            list_del += [j for j in my_fragList[i] if j not in list_del]
    if list_del:
        list_del.sort()
        print('Remove broken fragments containing:', list_del)
        update_frag_del(list_del, fn=fn_frag)
        del struct[list_del]
        write(poscar, struct, format='vasp')


# ------------------------------------------ constant potential, route B (GCE)
def gce_result(out, atoms, nsides=1):
    """Omega_el(U) from a &GCE run: E - E_F*N_e with E_F = -W_f.
    nsides=2 for a z-symmetrised slab (both faces charged), as in GOCIA's add_vaspResult_SC.
    Returns (atoms_with_energy, info) or (None, None)."""
    d = parse_output(out)
    if d is None or d['gce'] is None:
        return None, None
    g = d['gce']
    omega = grand_potential_el(d['energy_eV'], -g['W_f'], g['N_e'])
    a = result_atoms(atoms, omega / nsides)
    a.info.update(sc_U=g['U_SHE'], sc_eV=omega, N_e=g['N_e'], W_f=g['W_f'], E_eV=d['energy_eV'])
    return a, d


# ------------------------------------- constant potential, route C (charge scan)
def _sc_dirname(q, fractional):
    return f'q_{q:+.2f}' if fractional else f'q_{int(q):+d}'


def prepare_surfChrg(list_charge, template='../cp2k-sc.inp', poscar='POSCAR', sets=None, fractional=False):
    """Write one sub-directory per charge with sc.inp + includes; nothing is run.
    fractional=True sets NEX (NELECTRON_EXCESS, ai-cp2k-lpb build) instead of the integer CHARGE."""
    atoms = read(poscar)
    dirs = []
    for q in list_charge:
        d = _sc_dirname(q, fractional)
        os.makedirs(d, exist_ok=True)
        write_includes(atoms, d)
        s = dict(sets or {})
        if fractional:
            s['NEX'] = -q
            make_input(template, os.path.join(d, 'sc.inp'), 'sc', charge=0, sets=s)
        else:
            make_input(template, os.path.join(d, 'sc.inp'), 'sc', charge=int(q), sets=s)
        dirs.append(d)
    return dirs


def collect_surfChrg(list_charge, phi_she=4.43, fractional=False, out='sc.out'):
    """Parse the charge sub-directories -> sc.dat rows (dN, U_SHE, Omega_el).
    dN = excess electrons = -q (same sign convention as gocia.utils.vasp sc.dat).
    Omega_el = E - E_F*dN with W_f = -E_F (bulk electrolyte potential = 0)."""
    rows = []
    for q in list_charge:
        d = _sc_dirname(q, fractional)
        r = parse_output(os.path.join(d, out))
        if r is None or r['fermi_eV'] is None or not r['ended']:
            print(f'{d}: no energy/Fermi level or not ended, skipped')
            continue
        dN = -q
        wf = -r['fermi_eV']
        u = wf - phi_she
        rows.append((dN, u, grand_potential_el(r['energy_eV'], r['fermi_eV'], dN)))
    with open('sc.dat', 'w') as f:
        for dN, u, om in rows:
            f.write(f'{dN}\t{u}\t{om}\n')
    return rows


def do_surfChrg_batch(list_charge, cp2k_cmd, template='../cp2k-sc.inp', poscar='POSCAR',
                      phi_she=4.43, sets=None, fractional=False):
    """Serial prepare -> run -> collect. For packed parallel runs call prepare_surfChrg,
    launch cp2k in each q_* directory yourself, then collect_surfChrg."""
    home = os.getcwd()
    for d in prepare_surfChrg(list_charge, template, poscar, sets, fractional):
        os.chdir(d)
        run_cp2k(cp2k_cmd, 'sc.inp', 'sc.out')
        os.chdir(home)
    return collect_surfChrg(list_charge, phi_she, fractional)


def get_parabola(dir_sc='.'):
    """Quadratic fit Omega(U) = a U^2 + b U + c on sc.dat -> parabola.dat (same 3-number
    contract as gocia.utils.vasp.get_parabola, without the matplotlib/sklearn dependency)."""
    data = np.loadtxt(os.path.join(dir_sc, 'sc.dat'), ndmin=2)
    x, y = data[:, 1], data[:, 2]
    a, b, c = np.polyfit(x, y, 2)
    yfit = a * x**2 + b * x + c
    ss_res = float(((y - yfit) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    with open(os.path.join(dir_sc, 'parabola.dat'), 'w') as f:
        f.write(f'{a}\t{b}\t{c}')
    with open(os.path.join(dir_sc, 'parabola_fit.txt'), 'w') as f:
        f.write(f'a={a:.6f} b={b:.6f} c={c:.6f} R2={r2:.6f} C_eff={-2*a:.6f} eV/V^2 U0={-b/(2*a):.4f} V\n')
    return a, b, c, r2


def parabola_energy(u_she, parabola='parabola.dat'):
    a, b, c = np.loadtxt(parabola)
    return a * u_she**2 + b * u_she + c


def add_result_SC(pop, atoms, u_she, workdir='.', parabola='parabola.dat', nsides=2, isAlive=1):
    """Mirror of PopulationGrandCanonical.add_vaspResult_SC for a CP2K charge scan.
    Uses only the population's public methods; writes the same db fields."""
    if not os.path.isabs(parabola):
        parabola = os.path.join(workdir, parabola)
    a, b, c = np.loadtxt(parabola)
    ene_sc = a * u_she**2 + b * u_she + c
    # read the neutral energy BEFORE calc_grandPot: Interface() resets cell/pbc on the
    # Atoms in place, which invalidates a SinglePointCalculator's stored state
    try:
        ene_eV = atoms.get_potential_energy()
    except Exception:
        ene_eV = None
    grndPot = pop.calc_grandPot(atoms, ene_sc / nsides)
    myLabel = open(f'{workdir}/label').read()
    name = os.path.basename(os.path.abspath(workdir))
    if name.isdigit():          # ase.db rejects int-looking strings as values
        name = 'kid' + name
    print('\n%s IS BORN with G = %.3f eV\t[%s]' % (name, grndPot, myLabel))
    unique = pop.is_uniqueInAll(atoms, grndPot)
    if unique and grndPot < pop.get_GMrow()['grandPot'] and isAlive == 1:
        print(f' |- {name} is the new GM!')
        with open(os.path.join(workdir, '..', 'gmid'), 'w') as f:
            f.write(str(len(pop) + 1))   # id the new row gets below (upstream writes len(pop): off by one)
    pop.gadb.write(atoms, name=name, mag=0, eV=ene_eV if ene_eV is not None else ene_sc,
                   sc_U=u_she, sc_eV=ene_sc, grandPot=grndPot, a=a, b=b, c=c,
                   mated=0, done=1, alive=isAlive if unique else 0, label=myLabel)
    if not unique:
        print(f' |- {name} is a duplicate!')
    return grndPot


# TODO: PopulationGrandCanonicalPoly (fragments) - add_aseResult(fn_frag=...) already exists upstream
