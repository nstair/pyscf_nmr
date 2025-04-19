import json
import numpy as np
import matplotlib.pyplot as plt

model_name = 'benzene'

def lorentzian(x, x0, fwhm, intensity=1):
    gamma = fwhm/2
    return intensity*(1/np.pi)*(gamma/((x-x0)**2 + gamma**2))

# 1) Load the JSON data
with open(f'nmr_shieldings_{model_name}.json','r') as f:
    atoms = json.load(f)

# 2) Select only protons and convert to chemical shifts δ = σ_ref – σ_iso
sigma_ref = 31.88  # TMS reference in ppm
shifts = [sigma_ref - atom['iso_shielding']
          for atom in atoms if atom['symbol'] == 'H']

# 3) For benzene, the 6 protons are equivalent:
mean_shift = float(np.mean(shifts))
total_int  = len(shifts)

# 4) Simulate a Lorentzian spectrum
linewidth = 0.1  # ppm
x = np.linspace(mean_shift + 2, mean_shift - 2, 2000)
spectrum = lorentzian(x, mean_shift, linewidth, total_int)

# 5) Plot
plt.figure(figsize=(6,3))
plt.plot(x, spectrum, lw=2)
plt.xlabel('Chemical Shift (ppm)')
plt.ylabel('Intensity (a.u.)')
plt.title(f'Simulated $^1$H NMR Spectrum of {model_name}')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()