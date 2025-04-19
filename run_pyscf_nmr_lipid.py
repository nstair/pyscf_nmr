import json
import numpy as np
from pyscf import gto, scf
from pyscf.prop.nmr.rhf import NMR

# 1) Build lipid molecule
mol = gto.M(
    atom = '''
# 8‑carbon chain (approximate zig‑zag) + COOH
 C    0.0000    0.0000    0.0000
 C    1.5400    0.0000    0.0000
 C    3.0800    1.2600    0.0000
 C    4.6200    1.2600    0.0000
 C    6.1600    0.0000    0.0000
 C    7.7000    0.0000    0.0000
 C    9.2400    1.2600    0.0000
 C   10.7800    1.2600    0.0000
# Carboxyl group
 O   12.1300    0.0000    0.0000
 O   10.7800    2.5500    0.0000

# Terminal methyl protons
 H   -0.5200   -0.5200    0.9000
 H   -0.5200   -0.5200   -0.9000
 H    1.5400    0.9000    0.9000
 H    1.5400   -0.9000   -0.9000

# CH2 protons along the chain
 H    3.0800    2.1600    0.0000
 H    4.6200    0.3600    0.9000
 H    4.6200    0.3600   -0.9000
 H    6.1600   -0.9000    0.9000
 H    6.1600    0.9000   -0.9000
 H    7.7000    0.9000    0.0000
 H    9.2400    2.1600    0.0000

# Carboxylic‑OH proton
 H   10.7800    2.5500    0.9000
''',
    basis = '6-31G(d)',
    verbose=4,
).build()

model_name = 'lipid'

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