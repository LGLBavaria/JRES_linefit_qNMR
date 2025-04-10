# JRES_linefit_qNMR
This R code provides functions to evaluate J-Resolution 2D NMR-spectra and use the received signal information to create a Line-Fit algorithmn for 1D NMR-spectra. Via PULCON it is possible to quantify analytes within complex matrices by this code combination.

Please cite the following manuscripts if you use the package:
publication will be inserted as soon as it is published

# Usage
## JRES_import()
function explained here
## JRES_peak_picking()
function explained here
## JRES_signal_identification
### JRES_find_singulet()
function explained here
### JRES_find_even_signal()
function and subfunction explained here
### JRES_find_odd_signal()
function and subfunction explained here
## JRES_plot_spectrum()
function explained here

## read_1d_rnmrdata(path)
To read bruker NMR files the package rnmrfit (https://github.com/ssokolen/rnmrfit) is used. 
`read_1d_rnmrdata()` is a modified function to create NMRData1D objects, supplemented by the slot ‘DATA’ which contains the real numbers of the spectrum intensity, adjusted by the scaling factor nc.proc and additional acquisitiondata. The input value is the `path` to the spectra defined as `.../exp_number/`.

## extract_ppm(Spec)
`extract_ppm()` calculate the step size for the ppm axis of the NMRData1D object by the following formula: $spectral width (=Spec.sw)/(number of data points (=Spec.size) -1)$.
Using the calculated step size, the ppm axis is constructed from the maximum ppm value (Spec.maxppm) to the minimum ppm value (Spec.minppm).

## fit_signal()
The function `fit_signal()` is used to precisely analyse NMR spectra by fitting a (or multiple) pseudo-Voigt curve(s) into a baseline corrected spectral area.

$\eta \cdot \left( \frac{a}{1 + \left( \frac{x - x_0}{\gamma} \right)^2 } \right) + (1 - \eta) \cdot \left( a \cdot \exp \left( \frac{-(x - x_0)^2}{2 \sigma^2} \right) \right)$

<sub>η: Ratio of Lorentzian contribution to the pseudo-Voigt function.</sub>  
<sub>1 - η: Ratio of Gaussian contribution to the pseudo-Voigt function.</sub>  
<sub>a: Amplitude of the fit.</sub>  
<sub>x: signal center.</sub>  
<sub>x₀: Peak shift relative to the signal center; the combination of x and x₀ determines the actual peak position.</sub>  
<sub>γ: Half-width of the Lorentzian component.</sub>  
<sub>σ: Half-width of the Gaussian component.</sub>

To optimize the fitting parameters, a genetic algorithm (GA) (https://github.com/luca-scr/GA) was employed. 

The function is expanded to fit multiple peaks by defining their estimated appearance beforehand.

Therefore `fit_signal()` is based on the following required input values:  
`sampleID`: sample ID for visualisation  
`x`: chemical shift of the relevant spectral region  
`y`: signal intensity of this spectral region  
`peaknumber`: defining the amount of peaks fittet; choose from "one", "two", "three", "four", "five", "six", "seven", "eight"  
`peak_center`: center of the signal fit  
`tolerance`: tolerance accepted of the center  
`SF`: spectral frequency of the spectrometer to calculate Hz to ppm  
`peak_shifts_Hz`: Peak position to peak_center in Hz (from lowfield (high ppm) to highfield (low ppm)) without any sign (+/-)  
`area_ratios`: estimated integral ratio (from lowfield (high ppm) to highfield (low ppm))  

optional input values:  
`coupling_tolerance = 0.05`: accepted tolerance of peak_shifts_Hz, default: 0,05 = 5%  
`error_threshold = 0.35 * max(y)`: relative threshold to emphazise smaller errors, default: 35% of maximal signal intensity  
`error_factor = 5`: factor to minimize (<1) or emphazise (>1) errorcalculation of small errors, default: 5  
`minor_signal = F`: possible option to emphazise errorcalculation for minor signals, default: FALSE  
`minor_signal_factor = 0.1`: factor to minimize (<1) or emphazise (>1) errorcalculation for minor signals, default: 0.1  
`minor_signal_threshold = 2`: thershold for usage of minor_signal_factor, default: 2 (double max(y))  
`overfit_cor = F`: possible option to use an extra penalty for signals overfitting the spectral line, default: FALSE  
`overfit_factor = 10`: factor for extra penalty (overfitting the spectral line), default: low = 1, medium = 10, high = 1000  

The output of `fit_signal()` function is a ggplot of the fit in the defined sectral region and a list `fit_result` which contains:  
`$best_params`: optimised fit-parameters    
`fit-values` and `residuals` for the plot  
`integral_value`: total integral  
`mse`: mean squared error of the fit  
`plot`: a ggplot for visualising the fit in the defines spectral region with the spectral line and it's residuals



# Example Spectra and Usage
```
xxx
```
