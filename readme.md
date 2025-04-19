# NMR Shielding Simulation with PySCF and Matplotlib

This repository demonstrates how to compute NMR shielding constants for a molecule using PySCF, export the results to a JSON file, and then simulate and plot the corresponding NMR spectrum with NumPy and Matplotlib.

## Repository Structure
- `run_pyscf_nmr.py`: PySCF script to compute shieldings and write `nmr_shieldings.json`  
- `simulate_nmr.py`: Python script to read `nmr_shieldings.json` and plot the NMR spectrum  
- `nmr_shieldings.json`: Generated output (ignored by Git)  
- `README.md`: This file  

## Prerequisites
- Anaconda (if not already installed): download from https://www.anaconda.com/products/distribution  
- Git (optional, for cloning the repo)  

## Installation

1. **Install Anaconda (if needed)**  
    Download and install Anaconda from the official site (link above).

2. **Clone this repository**  
    git clone https://github.com/YourUser/YourRepo.git  
    cd YourRepo

3. **Create and activate the Conda environment**  
    conda create -n pyscf_nmr_env python=3.9 pyscf numpy matplotlib -c conda-forge  
    conda activate pyscf_nmr_env

4. **Install PySCF with the NMR properties add‑on**  
    pip install --upgrade pyscf pyscf-properties  
    (Or: pip install --upgrade "pyscf[properties]")

## Usage

1. **Compute NMR Shieldings**  
    python run_pyscf_nmr.py  
    Runs an RHF calculation and writes isotropic shieldings to `nmr_shieldings.json`.

2. **Simulate and Plot the NMR Spectrum**  
    python simulate_nmr.py  
    Reads the JSON file, converts shieldings to chemical shifts, and displays the simulated 1H NMR spectrum.

## Notes
- Replace the benzene geometry in `run_pyscf_nmr.py` with your molecule of interest.  
- Adjust `sigma_ref` and `linewidth` in `simulate_nmr.py` to fine‑tune the spectrum.  
- Ensure the Conda environment is active before running the scripts.

## License
This project is licensed under the MIT License. Feel free to use and modify as needed.