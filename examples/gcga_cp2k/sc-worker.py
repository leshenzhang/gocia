# Route C: a posteriori constant-potential energetics for ONE relaxed structure (thesis eq 5.5).
# Run inside a directory holding POSCAR (relaxed, z-symmetrised slab) and ../cp2k-sc.inp.
#   python sc-worker.py -2 -1 0 1 2        (integer charges; add --frac for NELECTRON_EXCESS)
# Produces sc.dat, parabola.dat, parabola_fit.txt. Feed into GOCIA with
#   gocia.utils.cp2k.add_result_SC(pop, atoms, u_she, nsides=2)   (mirror of add_vaspResult_SC)
import sys
sys.path.insert(0, '..')
import input
from gocia.utils import cp2k

args = [a for a in sys.argv[1:] if a != '--frac']
frac = '--frac' in sys.argv
charges = [float(a) if frac else int(a) for a in args] or [-2, -1, 0, 1, 2]
rows = cp2k.do_surfChrg_batch(charges, input.cp2k_cmd, template='../cp2k-sc.inp',
                              phi_she=input.phi_she, fractional=frac)
if len(rows) >= 3:
    a, b, c, r2 = cp2k.get_parabola('.')
    print(f'parabola a={a:.6f} b={b:.6f} c={c:.6f} R2={r2:.5f}  C_eff={-2*a:.4f} eV/V^2  U0={-b/(2*a):.4f} V')
else:
    sys.exit('fewer than 3 charge points succeeded; no fit')
