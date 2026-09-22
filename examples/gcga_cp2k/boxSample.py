# Initial population by box sampling (engine-agnostic; identical to examples/init/boxSample.py
# but driven by input.py). Usage, from the search root:
#   python boxSample.py substrate.vasp 40        -> init.db (done=0), N random structures
# input.py must define:
#   xyzLims                       sampling box [[x0,x1],[y0,y1],[z0,z1]]
#   zLim                          z window of the Interface
#   sampleSpecies = {'H': (4, 12)} element -> (min, max) count per structure
#   bondRejList = [['H','H']]     optional
import random
import sys
import numpy as np
import ase.io as ai
from ase.db import connect
from gocia.interface import Interface
from gocia.geom import build
import input

surfName = sys.argv[1]
nSample = int(sys.argv[2])
surf = Interface(tags=surfName.split('.')[0], allAtoms=ai.read(surfName), subAtoms=ai.read(surfName), zLim=input.zLim)
surf.print()

db = connect('init.db', append=False)
n_ok = 0
while n_ok < nSample:
    elems = []
    for el, (lo, hi) in input.sampleSpecies.items():
        elems += [el] * random.randint(lo, hi)
    new = build.boxSample_adatom(
        surf, elems,
        xyzLims=np.array(input.xyzLims),
        toler_BLmin=-0.2, toler_CNmax=3,
        bondRejList=getattr(input, 'bondRejList', None),
        doShuffle=True, rattle=True, rattleStdev=0.05,
    )
    if new is None:
        continue
    new.preopt_hooke(cutoff=1.2, toler=0.1)
    db.write(new.get_allAtoms(), done=0)
    n_ok += 1
    print(f'sample {n_ok}: {new.get_allAtoms().get_chemical_formula()}')
print('init.db written:', n_ok)
