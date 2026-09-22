# init.db -> s00001/POSCAR ... (one directory per sampled structure, VASP format keeps FixAtoms)
#   python db2dirs.py init.db
import os
import sys
from ase.db import connect
from ase.io import write

db = connect(sys.argv[1] if len(sys.argv) > 1 else 'init.db')
n = 0
for row in db.select():
    n += 1
    d = 's%05d' % n
    os.makedirs(d, exist_ok=True)
    write(os.path.join(d, 'POSCAR'), row.toatoms(), format='vasp')
print(f'{n} directories written (s00001..s{n:05d})')
