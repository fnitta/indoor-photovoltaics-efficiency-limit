# Efficiency Limit of Solar Cells Under Various Indoor Lighting Conditions
This repository contains MATLAB code to calculate the fundamental performance limits of single-junction solar cells with a realistic analysis based on the Tiedje-Yablonovitch model and including defect-assisted Shockley-Read-Hall (SRH) recombination. It calculates current density-voltage (J-V) characteristics and key photovoltaic parameters: short-circuit current (J<sub>sc</sub>), open-circuit voltage (V<sub>oc</sub>), fill factor (FF), and power conversion efficiency (eff). The MATLAB code also allows for visual representation of these parameters through plots to better understand the efficiency trends with varying thicknesses. Our work focuses on the efficiency limits of multilayer MoS<sub>2</sub>, MoSe<sub>2</sub>, WS<sub>2</sub>, and WSe<sub>2</sub> solar cells under compact fluorescent lamp (CFL), light-emitting diode (LED) lamp, halogen lamp, and low-light AM 1.5 G illumination. However, this code can be applied to calculate the efficiency limits of any material under any spectrum.

This README serves two purposes:
1. It provides an overview of the MATLAB code we employed to produce the results presented in our paper.
2. It guides users on the necessary modifications and data files required to calculate the efficiency limits of (a) their materials of interest and/or (b) under other spectra.

## Overview of the MATLAB Code
All of the data needed to reproduce the results in our paper is in this repository. The MATLAB script is divided into four main parts:

0. Parameter Initialization and Directory Setup
  - Defines constants and directory path of the main folder.
  - Initializes electrical parameters (band gap, effective electron mass, and effective hole mass) in Lines 21 - 23, and optical constants (n and k) in Lines 35 - 39.
  - Defines spectra in Lines 52 - 53.
  - Initializes recombination parameters in Lines 68 - 70.
  - Range of thicknesses defined in Line 82. Thickness(es) to generate J-V curves for are defined in Line 83 (note: make sure that the thicknesses are in ascending order, and are in the range of thicknesses).
  - Defines directory paths to store output data and figures.

1. Efficiency Limits via Shockley-Queisser (SQ) Model
  - Generates J-V characteristics and calculates J<sub>sc</sub>, V<sub>oc</sub>, FF, and eff using the SQ model.

2. Efficiency Limits via Extended Tiedje-Yablonovitch (TY) Model
  - Uses Equation 11 from Supplementary Note 1 (extended TY model) to generate J-V characteristics and calculate V<sub>oc</sub>, J<sub>sc</sub>, FF, and eff for a range of thicknesses.
  - Analysis of recombination magnitude is conducted. Data includes recombination units and recombination mechanism lifetimes in units of mA cm<sup>-2</sup>, and all recombination mechanism lifetimes in units of seconds.

3. Data Visualization
  - Plots for each parameter versus material thickness.
  - Saves plots in both .fig and .jpg formats. Each curve, corresponding to a different illumination intensity (or 𝜏<sub>SRH</sub> value), has a different color (or line style).

## Calculating the Efficiency Limits Under Other Conditions
We hope this code will be useful for calculating the efficiency limits under other conditions. Here is how you can modify the code to calculate the efficiency limits for other materials of interest and for other spectra.

### For Other Materials of Interest
- Initialize the electrical parameters for the material: band gap, effective electron mass, and effective hole mass in Lines 21 - 23. These parameters are used to calculate the effective conduction band density of states, valence band density of states, and intrinsic carrier concentration in Lines 27 - 32.
  - If the material's band gap is high, you may need to extend the range of the voltage in the SQ model in Line 133. 
- Define the optical constants (n and k) in Lines 35 - 39.
  - To minimize the amount of code to edit, define your material name with the variable Material_Name in Line 18 (this will be used for the plots); then, put your n and k data files in the folder "Data/Material/[Material_Name]" (without the square brackets) and name them "[Material_Name]-n.txt" and "[Material_Name]-k.txt", respectively. See how it is done for the example data (MoS<sub>2</sub>, MoSe<sub>2</sub>, WS<sub>2</sub>, and WSe<sub>2</sub>), which is in the repository.
  - Note: Make sure that your k data is defined for all **wavelengths** of your spectrum. If your k data is defined for energy (such as eV or J), use Planck's equation (λ = hc / E) to convert energy into wavelength. For AM 1.5 G, this means having k data for all wavelengths between 280 nm and 4000 nm, but for other spectra, this can be redefined in Lines 57 - 58. This may mean that you will have to extrapolate some of your k data. Make sure that your n data is defined at the band gap energy of your material, converted to wavelength. Also, make sure your n data is defined in terms of wavelengths, not energy.
- Define the recombination parameters in Lines 68 - 70.

### For Other Spectra
- Define the illumination spectra and their intensities in Lines 52 - 53.
    - If you use more than four illumination values, you will need to add more colors to the `Colors_All` variable in Lines 77 - 78.
    - See the note above for the k data. The spectrum should also be defined for all **wavelengths** rather than energy, and if it's defined for energy, Planck's equation should be used to convert energy into wavelength.
 
## Citing Our Work
If you find this model useful for your research, please cite our paper: Nitta, F. U., Nazif, K. N., & Pop, E. (2024). Transition Metal Dichalcogenide Solar Cells for Indoor Energy Harvesting. 
_arXiv preprint arXiv:2411.02642_.
