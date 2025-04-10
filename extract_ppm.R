#extract_ppm Function

# Description: resulting ppm variable is an array representing the chemical shift axis in ppm units corresponding to the spectrum data.
# Procedure: 
# Calculate the step size for the ppm axis: spectral width (=Spec.sw) / (number of data points (=Spec.size) -1)
# Using the calculated step size, the ppm axis is constructed from the maximum ppm value (Spec.maxppm) to the minimum ppm value (Spec.minppm)


extract_ppm <- function(Spec) {
  # results ppm axis for spectrum
  step <- Spec@parameters$sw / (Spec@parameters$size - 1)
  ppm <- seq(Spec@parameters$maxppm, Spec@parameters$minppm, by = -step)
  return(ppm)
}
