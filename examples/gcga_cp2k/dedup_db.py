# After run-init: drop duplicates among the relaxed initial structures in gcga.db
# (same test as the GA uses: identical formula + sorted-distance similarity + |dOmega| < eneCut).
#   python dedup_db.py gcga.db        -> duplicates get alive=0 (rows are kept, never deleted)
import sys
from ase.db import connect
from gocia.ensemble.comparator import srtDist_similar_zz

fn = sys.argv[1] if len(sys.argv) > 1 else 'gcga.db'
eneCut = 0.05
db = connect(fn)
rows = [r for r in db.select('done=1')]
kept = []
n_dup = 0
for r in sorted(rows, key=lambda r: r.grandPot):
    a = r.toatoms()
    dup = any(k.toatoms().get_chemical_formula() == a.get_chemical_formula()
              and abs(k.grandPot - r.grandPot) < eneCut
              and srtDist_similar_zz(a, k.toatoms()) for k in kept)
    if dup:
        db.update(r.id, alive=0)
        n_dup += 1
    else:
        kept.append(r)
print(f'{len(rows)} relaxed rows, {n_dup} duplicates set alive=0, {len(kept)} unique')
