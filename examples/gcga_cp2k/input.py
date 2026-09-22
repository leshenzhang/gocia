# GOCIA x CP2K example input (python data file, imported by the workers)

popSize = 10

# Reference for the fixed part (group B): energy of the bare substrate.vasp relaxed with
# cp2k-3.inp (route A, eV). For route B (constant potential) it must instead be
# Omega_el,B(U)/nsides of the bare symmetrised slab at the SAME target potential.
subsPot = 0.0            # <-- fill from the bare-slab run

# chemical potentials mu' (eV) of the exchangeable species (group A). H via CHE:
# mu'_H = 0.5*E_H2 + dE_gas - ln(10)*kT*pH - |e|*U_SHE - dE_ads
chemPotDict = {
    'H': -3.40,          # <-- placeholder, recompute for your level of theory / conditions
}

# z window of the sampling region (A); below = buffer, above = vacuum/electrolyte
zLim = [12.0, 16.0]
# sampling box for growMut_box / boxSample_adatom: [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
xyzLims = [[0.0, 11.1], [0.0, 9.6], [12.0, 16.0]]

# route: 'A' canonical/CHE (cp2k-1..3.inp), 'B' constant potential (cp2k-gce.inp)
route = 'A'
u_she = -0.3             # route B target potential (V vs SHE); also enters mu'_H via CHE
phi_she = 4.43           # absolute SHE potential used by the CP2K build (D-131 convention)
nsides = 2               # 2 = z-symmetrised slab (both faces charged), route B/C only

# launcher + binary. Shaheen wangc0i pack-4 rule (192 cores/node, 4 x 48): 24 MPI x 2 OMP each
cp2k_cmd = 'srun --exact --mem=90000 --hint=nomultithread -n 24 -c 2 /scratch/wangc0i/zls/soft/cp2k-2026.2/install/bin/cp2k.psmp'
