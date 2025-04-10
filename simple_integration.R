# Load necessary R-Packages
library(ggplot2)

# Function for simple integration of 1D NMR signals

# necessary input:
# x: vector with chemical shift values of selected spectral region
# y: vector with corresponding intensity values in selected spectral region
# peak_center: position of signal center on chem. shift axis (x) in ppm
# left_range: distance to left integral border in ppm
# right_range: distance to right integral border in ppm

# optional input:
# sampleID: ID-number of the examined sample, used for the plot heading (Default: NA)
# analyte: name of the examined compound, used for the plot heading (Default: NA)


integrate_signal <- function(x, y, peak_center, left_range, right_range, sampleID = NA, analyte  = NA){
  

  # Definition of the integral range
  integral_x <- which((
    x > (peak_center - right_range) &
      x < (peak_center + left_range)),
    arr.ind = T)
  
  # Correction of the baseline by miniumum in the integral range
  correction_range_x <- which((
    x > (peak_center - right_range) &
      x < (peak_center + left_range)),
    arr.ind = T)
  
  min_y <- min(y[correction_range_x])
  y_corrected <- y - min_y
  
  # Calculation of the integral in the selected area
  integral <- sum(y_corrected[integral_x])
  
  # Creating data frame for plotting
  b <- data.frame(x = x, y_corrected = y_corrected)
  integral_lower_border <- min(x[integral_x])
  integral_upper_border <- max(x[integral_x])
  
  p <- ggplot(b, aes(x = x)) + 
    geom_line(aes(y = y_corrected, color = "spectral data"), linewidth = 0.75) +
    geom_ribbon(data = subset(b, x > (peak_center - right_range) & x < (peak_center + left_range)), aes(ymin = 0, ymax = y_corrected, fill = "integral area"), alpha = 0.3) +
    geom_vline(xintercept = c(integral_lower_border, integral_upper_border), linetype = "dashed") +
    scale_x_reverse() +
    scale_color_manual("", values = c("spectral data" = "black")) +
    scale_fill_manual("", values = c("integral area" = "blue")) +
    labs(title = paste(sampleID, "-", analyte, "Integral"),
         x = "chemical shift (ppm)",
         y = "signal intensity") +
    theme_minimal() +
    theme(legend.position = "right", legend.text = element_text(size = 10))
  
  results_integration <- list(
    integral_value = integral,
    integral_x = integral_x,
    integral_y = y_corrected,
    plot = p
  )
  
  return(results_integration)
}
