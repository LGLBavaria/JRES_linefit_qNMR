# Load necessary R-Packages
library(dplyr)
library(ggplot2)
library(mrbin)
library(purrr)
library(tidyr)

# 1. Function for importing a 2D spectrum in bruker format (mrbin-dependent) -----------------
# automatic referencing on internal standard (TSP-d4 in our case) for complete comparability
# therefore Intenselimit can be adjusted for the peakpicking of the reference signal

# filepath to sample spectrum & sample name necessary
JRES_import <- function(path, IntenseLimit = 1000) {
  
  # Import of a 2D spectrum with mrbin function
  JRES2D_Import <- readBruker(
    folder = paste(path, "pdata", "1", sep = "/"),
    dimension = "2D"
  )
  
  JRES2D_Data <- as.data.frame(JRES2D_Import$currentSpectrum)
  
  #  Identification of the reference signal in JRES spectrum
  JRES_Peaks <- JRES_peak_picking(JRES2D_Data, 0.30, - 0.30, IntenseLimit = IntenseLimit)
  
  # Extraction Peak identification J-resolved
  sub_spectrum <- JRES_Peaks$sub_spectrum
  JRESPeaklist <- JRES_Peaks$JRESPeaklist
  
  # Identification of the analyte signal
  cat("\n")
  cat("Adjustment on reference signal\n")
  result_list <- JRES_find_singulet(JRESPeaklist, 0.000)
  
  # Extracting results of the J-resolved signal identification
  mark_signal_data <- result_list$mark_signal_data
  
  F1_TSP <- mark_signal_data$F1ppm
  F2_TSP <- mark_signal_data$F2ppm
  
  # Reading out the existing axes
  couple_vector <- as.numeric(rownames(JRES2D_Data)) # Kopplung: Zeilennnamen als PPM-Werte
  shift_vector <- as.numeric(colnames(JRES2D_Data))  # Verschiebung: Spaltennamen als PPM-Werte
  
  # Creating the adjusted axes
  couple_vector_new <- couple_vector - F1_TSP
  shift_vector_new <- shift_vector - F2_TSP
  
  # Assigning the adjusted axes
  rownames(JRES2D_Data) <- couple_vector_new
  colnames(JRES2D_Data) <- shift_vector_new
  
  # Assigning the adjusted spectrum to the imported spectrum
  JRES2D_Import$currentSpectrum <- JRES2D_Data
  
  # Return import
  return(JRES2D_Import)
}

# 2. Function for identifying all peaks in the defined section ----------------------------------------

# Default: Inspection area of 21x21 entries (10 left, 1 center, 10 right/10 top, 1 center, 10 bottom) in DF, intensity limit 1000
# Info: Naming for JRES: F1 = coupling axis in ppm; F2 = chem. shift axis in ppm

JRES_peak_picking <- function(spectrum, ppm_max, ppm_min,
                              horiz_range = 20, vert_range = 20,
                              IntenseLimit = 1000) {
  
  # Extract chemical shift axis (F2) from column names
  shift_vector <- as.numeric(colnames(spectrum))
  # Extract coupling axis (F1) from row names
  couple_vector <- as.numeric(rownames(spectrum))
  
  # Select relevant F2 region by finding indices within the specified ppm range
  ppm_indices <- which(shift_vector >= ppm_min & shift_vector <= ppm_max)
  # Extract the submatrix of interest for processing
  sub_spectrum <- spectrum[, ppm_indices, drop = FALSE]
  
  # Get dimensions of the sub-spectrum
  nrow_s <- nrow(sub_spectrum)
  ncol_s <- ncol(sub_spectrum)
  
  # Initialize a list to store found peaks (efficient for growing data)
  result_list <- list()
  count <- 1
  
  # Calculate half-widths for horizontal (F2) and vertical (F1) peak detection windows
  h_half <- floor(horiz_range / 2)
  v_half <- floor(vert_range / 2)
  
  # Loop over each point in the sub-spectrum, including border areas
  for (i in 1:nrow_s) {
    for (j in 1:ncol_s) {
      current_value <- sub_spectrum[i, j]
      # Skip if current intensity is below the defined threshold
      if (current_value <= IntenseLimit) next
      
      # Dynamically define the boundaries of the local window, adjusting for spectrum borders
      i_min <- max(1, i - v_half) # Ensure window does not go below 1st row
      i_max <- min(nrow_s, i + v_half) # Ensure window does not exceed last row
      j_min <- max(1, j - h_half) # Ensure window does not go below 1st column
      j_max <- min(ncol_s, j + h_half) # Ensure window does not exceed last column
      
      # Extract the local window around the current point
      window <- sub_spectrum[i_min:i_max, j_min:j_max, drop = FALSE]
      
      # Determine the center's position within the extracted local window
      center_in_window_i <- i - i_min + 1
      center_in_window_j <- j - j_min + 1
      
      # Create a temporary window with the center value set to NA for comparison
      window_center_removed <- window
      window_center_removed[center_in_window_i, center_in_window_j] <- NA
      
      # Check if the current point is a local maximum (i.e., greater than all its neighbors)
      if (current_value > max(window_center_removed, na.rm = TRUE)) {
        # Store the F1 ppm, F2 ppm, and Intensity of the detected peak
        result_list[[count]] <- c(
          F1ppm = couple_vector[i],
          F2ppm = shift_vector[ppm_indices[j]],
          Intensity = current_value
        )
        count <- count + 1
      }
    }
  }
  
  # Handle the case where no peaks are found
  if (length(result_list) == 0) {
    JRESPeaklist <- data.frame(F1ppm = numeric(0), F2ppm = numeric(0), Intensity = numeric(0))
  } else {
    # Combine all found peaks from the list into a data frame
    JRESPeaklist <- do.call(rbind, result_list)
    JRESPeaklist <- as.data.frame(JRESPeaklist)
  }
  # Ensure correct column names for the peak list
  colnames(JRESPeaklist) <- c("F1ppm", "F2ppm", "Intensity")
  
  # Return both the processed sub-spectrum and the detected peak list
  return(list(
    sub_spectrum = sub_spectrum,
    JRESPeaklist = JRESPeaklist
  ))
}

