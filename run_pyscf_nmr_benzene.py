import json
import numpy as np
from pyscf import gto, scf
from pyscf.prop.nmr.rhf import NMR

# 0) Create conda environment and install conda dependancies, pip install pyscf and properties package (for nmr)
# conda create -n pyscf_nmr_env python=3.9 numpy matplotlib -c conda-forge
# conda activate pyscf_nmr_env
# pip install --upgrade pyscf pyscf-properties

# 1) Build benzene molecule
mol = gto.M(
    atom = '''
      C   0.000000   1.402720   0.000000
      C   1.214790   0.701360   0.000000
      C   1.214790  -0.701360   0.000000
      C   0.000000  -1.402720   0.000000
      C  -1.214790  -0.701360   0.000000
      C  -1.214790   0.701360   0.000000
      H   0.000000   2.490290   0.000000
      H   2.156660   1.245150   0.000000
      H   2.156660  -1.245150   0.000000
      H   0.000000  -2.490290   0.000000
      H  -2.156660  -1.245150   0.000000
      H  -2.156660   1.245150   0.000000
    ''',
    basis   = '6-31G(d)',
    verbose = 1,
).build()

model_name = 'benzene'

# 2) Run RHF
mf = scf.RHF(mol)
mf.kernel()

# 3) Compute GIAO NMR shielding tensors
nmr = NMR(mf)
shield_list = nmr.kernel()  
# shield_list is a list of length mol.natm, each a 3×3 array

# 4) Compute isotropic shieldings
iso_shieldings = []
for idx, tensor in enumerate(shield_list):
    arr = np.array(tensor, dtype=float)
    if arr.shape != (3,3):
        raise ValueError(f"Atom {idx}: expected 3×3, got {arr.shape}")
    iso = float(np.trace(arr) / 3.0)
    iso_shieldings.append(iso)

# 5) Write JSON using mol.atom_symbol(i)
output = []
for i in range(mol.natm):
    sym = mol.atom_symbol(i)
    output.append({
        'index':         i,
        'symbol':        sym,
        'iso_shielding': iso_shieldings[i]
    })

with open(f'nmr_shieldings_{model_name}.json','w') as f:
    json.dump(output, f, indent=2)

print("✔ Written isotropic shieldings to nmr_shieldings.json")