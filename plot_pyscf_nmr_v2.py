import json
import numpy as np
import matplotlib.pyplot as plt

model_name = 'lipid'

def lorentzian(x, x0, fwhm, intensity=1):
    """Lorentzian line shape: center x0, FWHM fwhm, scaled by intensity."""
    gamma = fwhm/2.0
    return intensity * (1/np.pi) * (gamma / ((x - x0)**2 + gamma**2))

def simulate_spectrum(chemical_shifts, intensities, linewidth, x_min, x_max, num_points=2000):
    x = np.linspace(x_min, x_max, num_points)
    spectrum = np.zeros_like(x)
    for cs, inten in zip(chemical_shifts, intensities):
        spectrum += lorentzian(x, cs, linewidth, inten)
    return x, spectrum

# 1) Load the shielding data
with open(f"nmr_shieldings_{model_name}.json","r") as f:
    atoms = json.load(f)

# 2) Extract *all* proton shifts
sigma_ref = 31.88  # TMS reference in ppm
shifts = []
for entry in atoms:
    if entry["symbol"].upper() == "H":
        sigma_iso = entry["iso_shielding"]
        shifts.append(sigma_ref - sigma_iso)

# 3) Set each proton to unit intensity
intensities = [1]*len(shifts)

# 4) Determine plotting range ±1 ppm around the extrema
margin = 1.0
x_min = min(shifts) - margin
x_max = max(shifts) + margin

# 5) Simulate
linewidth = 0.1  # ppm
x, spectrum = simulate_spectrum(shifts, intensities, linewidth, x_min, x_max)

# 6) Plot
plt.figure(figsize=(7,3))
plt.plot(x, spectrum, lw=2)
plt.xlabel('Chemical Shift (ppm)')
plt.ylabel('Intensity (a.u.)')
plt.title('Simulated $^1$H NMR Spectrum')
plt.gca().invert_xaxis()

# Optional: mark the individual peak centers
for cs in shifts:
    plt.axvline(cs, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()