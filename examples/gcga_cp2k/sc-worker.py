# Route C: a posteriori constant-potential energetics for ONE relaxed structure (thesis eq 5.5).
# Run inside a directory holding POSCAR (relaxed, z-symmetrised slab) and ../cp2k-sc.inp.
#   python sc-worker.py -2 -1 0 1 2               serial: prepare, run, collect, fit
#   python sc-worker.py --prepare -2 -1 0 1 2     only write q_*/sc.inp (then run them packed)
#   python sc-worker.py --collect -2 -1 0 1 2     parse q_*/sc.out -> sc.dat -> parabola.dat
# add --frac for fractional charges via NELECTRON_EXCESS.
# Feed into GOCIA with gocia.utils.cp2k.add_result_SC(pop, atoms, u_she, nsides=2).
import sys
sys.path.insert(0, '..')
import input
from gocia.utils import cp2k

flags = {a for a in sys.argv[1:] if a.startswith('--')}
args = [a for a in sys.argv[1:] if not a.startswith('--')]
frac = '--frac' in flags
charges = [float(a) if frac else int(a) for a in args] or [-2, -1, 0, 1, 2]
template = getattr(input, 'sc_template', '../cp2k-sc.inp')

if '--prepare' in flags:
    print('prepared:', cp2k.prepare_surfChrg(charges, template=template, fractional=frac))
    sys.exit(0)
if '--collect' in flags:
    rows = cp2k.collect_surfChrg(charges, phi_she=input.phi_she, fractional=frac)
else:
    rows = cp2k.do_surfChrg_batch(charges, input.cp2k_cmd, template=template,
                                  phi_she=input.phi_she, fractional=frac)
for dN, u, om in rows:
    print(f'dN={dN:+.2f}  U={u:+.4f} V vs SHE  Omega={om:.6f} eV')
if len(rows) >= 3:
    a, b, c, r2 = cp2k.get_parabola('.')
    print(f'parabola a={a:.6f} b={b:.6f} c={c:.6f} R2={r2:.6f}  C_eff={-2*a:.4f} eV/V^2  U0={-b/(2*a):.4f} V')
else:
    sys.exit('fewer than 3 charge points succeeded; no fit')
