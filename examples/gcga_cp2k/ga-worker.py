# GCGA worker (CP2K). Run inside its own sub-directory of the search root:
#   root/gcga.db  root/substrate.vasp  root/input.py  root/cp2k-*.inp  root/000123/  <- cwd
import os
import sys
import numpy as np
sys.path.insert(0, '..')
import input
from ase.io import read
from gocia.ga.popGrandCanon import PopulationGrandCanonical
from gocia.utils import cp2k

pop = PopulationGrandCanonical(
    gadb='../gcga.db',
    substrate='../substrate.vasp',
    popSize=input.popSize,
    subsPot=input.subsPot,
    chemPotDict=input.chemPotDict,
    zLim=input.zLim,
)

# --- child generation (identical to the VASP example) ---
kid = None
while kid is None:
    kid = pop.gen_offspring_box(
        mutRate=0.4,
        xyzLims=np.array(input.xyzLims),
        bondRejList=getattr(input, 'bondRejList', None),
        constrainTop=True,
        transVec=[[-3, 3], [-3, 3]],
    )
kid.preopt_hooke(cutoff=1.2, toler=0.1)
kid.write('POSCAR')

# --- local optimisation with CP2K ---
if input.route == 'A':
    atoms = cp2k.do_multiStep_opt(step=3, cp2k_cmd=input.cp2k_cmd, template='../cp2k-%i.inp',
                                  zLim=input.zLim, substrate='../substrate.vasp')
    if atoms is None:
        sys.exit('CP2K failed (FAIL touched)')
    pop.add_aseResult(atoms, workdir='.')

elif input.route == 'B':
    # stages 1-2 neutral (cheap), final stage at constant potential
    atoms = cp2k.do_multiStep_opt(step=2, cp2k_cmd=input.cp2k_cmd, template='../cp2k-%i.inp',
                                  zLim=input.zLim, substrate='../substrate.vasp', clean=False)
    if atoms is None:
        sys.exit('CP2K failed (FAIL touched)')
    cp2k.write_includes(atoms)
    cp2k.make_input('../cp2k-gce.inp', 'gce.inp', 'gce',
                    sets=dict(cp2k.default_sets(atoms, 0), USHE=input.u_she, PHISHE=input.phi_she))
    cp2k.run_cp2k(input.cp2k_cmd, 'gce.inp', 'gce.out')
    final = cp2k.read_final_geometry('gce', atoms)
    res, info = cp2k.gce_result('gce.out', final, nsides=input.nsides)
    if res is None:
        os.system('touch FAIL')
        sys.exit('no GCE trace in gce.out')
    from ase.io import write
    write('out-gce.vasp', res, format='vasp')
    pop.add_aseResult(res, workdir='.')     # grandPot = Omega_el/nsides - subsPot - sum(mu'N)

pop.natural_selection()
os.system('rm -f *.wfn *.wfn.bak-* *.cube')
