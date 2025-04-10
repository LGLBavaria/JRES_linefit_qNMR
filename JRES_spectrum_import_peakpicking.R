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

JRES_peak_picking <- function(spectrum, ppm_max, 
                            ppm_min,
                            horiz_range = 20,
                            vert_range = 20,
                            IntenseLimit = 1000) {
  
  # Saving the axis information as vectors
  shift_vector <- as.numeric(colnames(spectrum))  # Shift: Column names as ppm values
  couple_vector <- as.numeric(rownames(spectrum)) # Coupling: Line names as ppm values
  
  # Find the indices that correspond to the ppm range
  ppm_indices <- which(shift_vector >= ppm_min & shift_vector <= ppm_max)
  
  # Extract submatrix containing the ppm range
  sub_spectrum <- spectrum[, ppm_indices]
  
  
  # Create an empty dataframe for the results
  JRESPeaklist <- data.frame(
    "F1ppm" = numeric(),
    "F2ppm" = numeric(),
    Intensity = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(sub_spectrum)) {
    for (j in 1:ncol(sub_spectrum)) {
      # Current point (intensity)
      current_value <- sub_spectrum[i, j]
      
      # Define the boundaries of the neighboring area taking into account the dataframe border
      i_min <- max(1, i - vert_range/2)
      i_max <- min(nrow(sub_spectrum), i + vert_range/2)
      j_min <- max(1, j - horiz_range/2)
      j_max <- min(ncol(sub_spectrum), j + horiz_range/2)
      
      # Extract the area of the neighboring cells
      neighbors <- sub_spectrum[i_min:i_max, j_min:j_max]
      
      # Remove the center (the current point) from the neighboring cell view
      neighbors_without_center <- neighbors
      center_i <- i - i_min + 1
      center_j <- j - j_min + 1
      neighbors_without_center[center_i, center_j] <- NA
      
      # Check whether the current point is a local maximum
      if (current_value > IntenseLimit &&
          current_value > max(neighbors_without_center, na.rm = TRUE)) {
        
        # Find original position
        i2 <- rownames(sub_spectrum)[i]
        j2 <- colnames(sub_spectrum)[j]
        
        # Add the local maximum to JRESPeaklist
        JRESPeaklist <- rbind(JRESPeaklist,
                              data.frame(
                                "F1ppm" = as.numeric(i2),
                                "F2ppm" = as.numeric(j2),
                                Intensity = current_value
                              ))
      }
    }
  }
  
  JRES_Peaks <- list(
    sub_spectrum = sub_spectrum,
    JRESPeaklist = JRESPeaklist
  )
  
  return(JRES_Peaks)
}
