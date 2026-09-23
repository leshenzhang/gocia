# Relax one initial-population structure with CP2K and write it straight into gcga.db
# with done=1 (GOCIA's get_GMrow needs >=1 relaxed row before the GA starts; D-13 lesson).
# Run inside root/s0000N/ containing POSCAR (from db2vasp.py of the sampled *.db).
import os
import sys
sys.path.insert(0, '..')
import input
from ase.db import connect
from ase.io import read
from gocia.utils import cp2k
from gocia.interface import Interface

atoms = cp2k.do_multiStep_opt(step=3, cp2k_cmd=input.cp2k_cmd, template='../cp2k-%i.inp',
                              zLim=input.zLim, substrate='../substrate.vasp')
if atoms is None:
    sys.exit('CP2K failed (FAIL touched)')

subs = read('../substrate.vasp')
ene = atoms.get_potential_energy()
grandPot = ene - input.subsPot
for s in Interface(atoms, subs).get_adsAtoms().get_chemical_symbols():
    if s in input.chemPotDict:
        grandPot -= input.chemPotDict[s]
name = os.path.basename(os.getcwd())
with connect('../gcga.db') as db:
    db.write(atoms, name=name, mag=0, eV=ene, grandPot=grandPot,
             mated=0, done=1, alive=1, label='0 0 init')
print(f'{name} written: E = {ene:.4f} eV, Omega = {grandPot:.4f} eV')
os.system('rm -f *.wfn *.wfn.bak-* *.cube')
