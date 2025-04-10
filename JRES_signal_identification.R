# Load necessary R-Packages
library(dplyr)
library(ggplot2)
library(mrbin)
library(purrr)
library(tidyr)

# 1. singulet identification----------------------------------------------
# singulets separated from odd numbered multiplicity, because of its special properties

# necessary input:
# JRESPeaklist: Dataframe with the peaks picked by JRES_peak_picking()
# F2shifttarget: setpoint for assumed signal position on F2-axis
# F1shifttol: tolerance granted for middle signals to shift up/down in ppm (Default: 0.0005 #[ppm])
# F2shifttol: tolerance granted for the peaks of one signal to differ in their F2 value (Default: 0.002 #[ppm])

JRES_find_singulet <- function(JRESPeaklist, F2shifttarget, F1shifttol = 0.0005, F2shifttol = 0.001){
  
  cat("\n")
  
  # Identification of the center peak ---------------------------------------
  # Search for signals on/near the zero line (F1)
  middle_peak <- which((
    JRESPeaklist$F1ppm > ( - F1shifttol) &
      JRESPeaklist$F1ppm < F1shifttol
  ),
  arr.ind = T
  )
  
  # Check for coupling peaks -------------------------------------------
  # Creating the deletion vector
  delete_vector <- c(rep(FALSE, length(middle_peak)))
  
  # Search for potentially coupling peaks
  for (n in 1:length(middle_peak)) {
    row <- middle_peak[n]
    
    # Extract the corresponding values from JRESPeakList
    F2value <- JRESPeaklist[row, "F2ppm"]
    
    couplesignals <- which((
      JRESPeaklist$F2ppm > ( F2value - F2shifttol) &
        JRESPeaklist$F2ppm < F2value + F2shifttol
    ),
    arr.ind = T
    )
    
    # Check whether there is more than one peak at this chemical shift
    if (length(couplesignals) > 1) {
      delete_vector[n] <- TRUE
    }
  }
  # Delete the selected items
  middle_peak <- middle_peak[!delete_vector]
  
  # Predefining error messages
  Error_JRES_singulet0 <- FALSE
  Error_JRES_singulet1 <- FALSE
  
  # Creating the data frame for marking the signal in the plot function
  mark_signal_data <- data.frame(
    "F1ppm" = numeric(),
    "F2ppm" = numeric(),
    "Intensity" = as.numeric(),
    stringsAsFactors = FALSE
  )
  
  # Results section: Checking the number of selected signals and outputting the results ------------
  
  # Check whether a singlet has been selected
  if (length(middle_peak) == 0) {
    F1_singuletppm <- 0
    F2_singuletppm <- F2shifttarget
    
    Error_JRES_singulet0 <- TRUE
    
    cat("JRES: No singulet was identified in the evaluated region!")
  } else {
    if (length(middle_peak) == 1) {
      
      # Output of the selected signal
      F1_singuletppm <- JRESPeaklist$F1ppm[middle_peak] 
      F2_singuletppm <- JRESPeaklist$F2ppm[middle_peak]
      singulet_Intensity <- JRESPeaklist$Intensity[middle_peak]
      
      # Appending the peak data to mark_signal_data for marking the selected signal in Plot2D
      mark_signal_data <- rbind(mark_signal_data, data.frame(
        "F1ppm" = as.numeric(F1_singuletppm),
        "F2ppm" = as.numeric(F2_singuletppm),
        "Intensity" = as.numeric(singulet_Intensity)
      ))
      
      cat("JRES: Exactly one singulet was identified in the evaluated region!")
      cat(paste("\nchem. shift:", round(
        F2_singuletppm, digits = 4
      ), "ppm", sep = " "))
    }
    
    if (length(middle_peak) >= 2) {
      vector_shiftdistance <- c()
      
      for (l in 1:length(middle_peak)) {
        shift_distance <- abs(JRESPeaklist$F2ppm[middle_peak[l]] - F2shifttarget)
        vector_shiftdistance <- c(vector_shiftdistance, shift_distance)
      }
      
      nearest_signal <- as.numeric(which.min(vector_shiftdistance))
      
      # Output of the selected signal
      F1_singuletppm <- JRESPeaklist$F1ppm[middle_peak[nearest_signal]] 
      F2_singuletppm <- JRESPeaklist$F2ppm[middle_peak[nearest_signal]]
      singulet_Intensity <- JRESPeaklist$Intensity[middle_peak[nearest_signal]]
      
      # Appending the peak data to mark_signal_data for marking the selected signal in Plot2D
      mark_signal_data <- rbind(mark_signal_data, data.frame(
        "F1ppm" = as.numeric(F1_singuletppm),
        "F2ppm" = as.numeric(F2_singuletppm),
        "Intensity" = as.numeric(singulet_Intensity)
      ))
      
      cat("JRES: More than one singulet was identified in the evaluated region! Signal with the least shift to the F2 target was selected."
      )
      cat(paste("\nchem. shift:", round(
        F2_singuletppm, digits = 4
      ), "ppm", sep = " "))
      
      Error_JRES_singulet1 <- TRUE
    }
  }
  
  # Creating result list of data to be output
  result_list <- list(
    mark_signal_data = mark_signal_data,
    peak_center = F2_singuletppm,
    Error_JRES0 = Error_JRES_singulet0,
    Error_JRES1 = Error_JRES_singulet1
  )
  
  cat("\n")
  
  return(result_list)
  
}

# 2. Identification of even-numbered coupled signals (2, 4, 6, 8, etc.) ----------------------------------------------

# necessary input:
# JRESPeaklist: Dataframe with the peaks picked by JRES_peak_picking()
# F2shifttarget: setpoint for assumed signal position on F2-axis
# coupleconst: vector with assumed coupling distance to signal center in ppm (for example c(0.1, 0.1))
# multiplicity: multiplicity of the signal (2 ,4 ,5 ,6 etc.)
# coupletol: tolerance granted when evaluating the coupleconstant (Default: 0.05 = 5%)
# F1shifttol: tolerance granted for middle signals to shift up/down in ppm (Default: 0.0005 #[ppm])
# F2shifttol: tolerance granted for the peaks of one signal to differ in their F2 value (Default: 0.002 #[ppm])
# check_height: turns height check of signals on and off (TRUE = On, False = Off; Default: FALSE)

# optional input:
# heightratio: vector with realtive height information, length has to equal multiplicity (Default: NA when not used)
# heighttol: tolerance granted when evaluating height information (Default: 0.1 = 10%)

JRES_find_even_signal <- function(JRESPeaklist,
                               F2shifttarget,
                               coupleconst,
                               multiplicity,
                               coupletol = 0.05,
                               F1shifttol = 0.0005,
                               F2shifttol = 0.002,
                               check_height = FALSE,
                               heightratio = NA,
                               heighttol = 0.1) {
  
  cat("\n")
  
  # 1. Excluding odd-numbered multiplicities -------------------------
  
  even_couplings <- cut_odd_couplings(JRESPeaklist, F1shifttol = F1shifttol, F2shifttol = F2shifttol)
  
  # 2. Structuring and pre-selection of possible couplings -------------------
  
  diff_even_coupling <- filter_even_couplings(even_couplings, F2shifttol = F2shifttol)
  
  # 3. Selection according to desired even-numbered multiplicity -----------------
  
  multiplicity_checked <- find_even_coupled_signals(diff_even_coupling,
                                                    multiplicity = multiplicity,
                                                    F2shifttol = F2shifttol)
  
  # 4. checking the coupling information -----------------------------------
  
  couple_checked <- check_even_coupling(multiplicity_checked, coupleconst, multiplicity = multiplicity, coupletol = coupletol)
 
  # 5.  Optional: Checking the peak height ratio ------------------------------
  
  height_checked <- check_height_even_signal(couple_checked,
                                             check_height = check_height,
                                             heightratio,
                                             heighttol = heighttol)
 
  # 6. results section ---------------------------------------------------------
  
  result_list <- results_even_signal(height_checked, multiplicity, F2shifttarget, coupleconst)
  
  cat("\n")
  
  return(result_list)
}


# 3. Identification of odd-numbered couplings (3, 5, 7, 9, etc.) -----------

# necessary input:
# JRESPeaklist: Dataframe with the peaks picked by JRES_peak_picking()
# F2shifttarget: setpoint for assumed signal position on F2-axis
# coupleconst: vector with assumed coupling distance to signal center in ppm (for example c(0.1, 0.1))
# multiplicity: multiplicity of the signal (3, 5, 7, 9, etc.)
# coupletol: tolerance granted when evaluating the coupleconstant (Default: 0.05 = 5%)
# F1shifttol: tolerance granted for middle signals to shift up/down in ppm (Default: 0.0005 #[ppm])
# F2shifttol: tolerance granted for the peaks of one signal to differ in their F2 value (Default: 0.003 #[ppm])
# check_height: turns height check of signals on and off (TRUE = On, False = Off; Default: FALSE)

# optional input:
# heightratio: vector with realtive height information, length has to equal multiplicity (Default: NA when not used)
# heighttol: tolerance granted when evaluating height information (Default: 0.1 = 10%)

JRES_find_odd_signal <- function(JRESPeaklist,
                              F2shifttarget,
                              coupleconst,
                              multiplicity,
                              coupletol = 0.05,
                              F1shifttol = 0.0005,
                              F2shifttol = 0.003,
                              check_height = FALSE,
                              heightratio = NA,
                              heighttol = 0.1) {
  
  cat("\n")
  
  # 1. Exclude even multiplicities -------------------------
  
  middle_peak <- cut_even_couplings(JRESPeaklist, F1shifttol = F1shifttol, F2shifttol = F2shifttol)
  
  # 2. Structuring and preselection of possible couplings -------------------
  
  diff_odd_coupling <- filter_odd_couplings(middle_peak, F2shifttol = F2shifttol)
  
  # 3. Selection according to desired odd multiplicity -----------------
  
  multiplicity_checked <- find_odd_coupled_signals(diff_odd_coupling,
                                                      multiplicity = multiplicity,
                                                      F2shifttol = F2shifttol)
  
  # 4. Checking the coupling information -----------------------------------
  
  couple_checked <- check_odd_coupling(multiplicity_checked, coupleconst, multiplicity = multiplicity, coupletol = coupletol)
  
  # 5. Optional: Checking the peak height ratio ------------------------------
  
  height_checked <- check_height_odd_signal(couple_checked,
                                               check_height = check_height,
                                               heightratio,
                                               heighttol = heighttol)
  
  # 6. Results section ---------------------------------------------------------
  
  result_list <- results_odd_signal(height_checked, multiplicity, F2shifttarget, coupleconst)
  
  cat("\n")
  
  return(result_list)
  
}

