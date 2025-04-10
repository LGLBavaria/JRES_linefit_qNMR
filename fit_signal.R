#LineFit function:
#fit_signal()
#Written and Copyright by Sabine Milbert, Bavarian Health and Food Safety authority, 2024


#Description of parameters, which have to be defined at function usage:
#sampleID                              sample ID for visualisation
#x                                     chemical shift of the relevant spectral region
#y                                     signal intensity of this spectral region
#peaknumber                            defining the amount of peaks fittet; choose from "one", "two", "three", "four", "five", "six", "seven", "eight"
#peak_center                           center of the signal fit
#tolerance                             tolerance accepted of the center
#SF                                    spectral frequency of the spectrometer to calculate Hz to ppm 
#peak_shifts_Hz                        Peak position to peak_center in Hz (from lowfield (high ppm) to highfield (low ppm)) without any sign (+/-)   
#area_ratios                           estimated integral ratio (from lowfield (high ppm) to highfield (low ppm))                         

#optional parameters, which may be defined at function usage:
#coupling_tolerance = 0.05             accepted tolerance of peak_shifts_Hz, default: 0,05 = 5%
#error_threshold = 0.35 * max(y)       relative threshold to emphazise smaller errors, default: 35% of maximal signal intensity
#error_factor = 5                      factor to minimize (<1) or emphazise (>1) errorcalculation of small errors, default: 5
#minor_signal = F                      possible option to emphazise errorcalculation for minor signals, default: FALSE
#minor_signal_factor = 0.1             factor to minimize (<1) or emphazise (>1) errorcalculation for minor signals, default: 0.1
#minor_signal_threshold = 2            thershold for usage of minor_signal_factor, default: 2 (double max(y))
#overfit_cor = F                       possible option to use an extra penalty for signals overfitiing the spectral line, default: FALSE
#overfit_factor = 10                   factor for extra penalty (overfitting the spectral line), default: low = 1, medium = 10, high = 1000

#----------------------------------------- definition of peak-functions ------------------------------------------------------

#Definition of the Pseudo-Voigt-function for single peaks
#x: signal center
#a: aplitude of the fit
#x0: peak shift to signal center --> x and x0 resultion in the actual peak position
#gamma: half width of the Lorentz-function
#sigma: half width of the Gauss-function 
#eta = ratio of Lorentz in Pseudo-Voigt
#1 - eta = ratio of Gauss in Pseudo-Voigt
pseudo_voigt_single <- function(x, a, x0, sigma, gamma, eta) { 
  eta * (a / (1 + ((x - x0) / gamma)^2)) + 
    (1 - eta) * (a * exp(-((x - x0)^2) / (2 * sigma^2)))
}

# Definition of the Pseudo-Voigt-function for two peaks
pseudo_voigt_twoPeaks <- function(x, a1, x0, sigma, gamma, eta, shifts, ratios) {
  peak1 <- pseudo_voigt_single(x, a1 * ratios[1], x0 + shifts[1], sigma, gamma, eta)
  peak2 <- pseudo_voigt_single(x, a1 * ratios[2], x0 - shifts[2], sigma, gamma, eta)
  peak1 + peak2
}

# Definition of the Pseudo-Voigt-function for three peaks
pseudo_voigt_threePeaks <- function(x, a1, x0, sigma, gamma, eta, shifts, ratios) {
  peak1 <- pseudo_voigt_single(x, a1 * ratios[1], x0 + shifts[1], sigma, gamma, eta)
  peak2 <- pseudo_voigt_single(x, a1 * ratios[2], x0 + shifts[2], sigma, gamma, eta)
  peak3 <- pseudo_voigt_single(x, a1 * ratios[3], x0 - shifts[3], sigma, gamma, eta)
  peak1 + peak2 + peak3
}

# Definition of the Pseudo-Voigt-function for four peaks
pseudo_voigt_fourPeaks <- function(x, a1, x0, sigma, gamma, eta, shifts, ratios) {
  peak1 <- pseudo_voigt_single(x, a1 * ratios[1], x0 + shifts[1], sigma, gamma, eta)
  peak2 <- pseudo_voigt_single(x, a1 * ratios[2], x0 + shifts[2], sigma, gamma, eta)
  peak3 <- pseudo_voigt_single(x, a1 * ratios[3], x0 - shifts[3], sigma, gamma, eta)
  peak4 <- pseudo_voigt_single(x, a1 * ratios[4], x0 - shifts[4], sigma, gamma, eta)
  peak1 + peak2 + peak3 + peak4
}

# Definition of the Pseudo-Voigt-function for five peaks
pseudo_voigt_fivePeaks <- function(x, a1, x0, sigma, gamma, eta, shifts, ratios) {
  peak1 <- pseudo_voigt_single(x, a1 * ratios[1], x0 + shifts[1], sigma, gamma, eta)
  peak2 <- pseudo_voigt_single(x, a1 * ratios[2], x0 + shifts[2], sigma, gamma, eta)
  peak3 <- pseudo_voigt_single(x, a1 * ratios[3], x0 + shifts[3], sigma, gamma, eta)
  peak4 <- pseudo_voigt_single(x, a1 * ratios[4], x0 - shifts[4], sigma, gamma, eta)
  peak5 <- pseudo_voigt_single(x, a1 * ratios[5], x0 - shifts[5], sigma, gamma, eta)
  peak1 + peak2 + peak3 + peak4 + peak5
}

# Definition of the Pseudo-Voigt-function for six peaks
pseudo_voigt_sixPeaks <- function(x, a1, x0, sigma, gamma, eta, shifts, ratios) {
  peak1 <- pseudo_voigt_single(x, a1 * ratios[1], x0 + shifts[1], sigma, gamma, eta)
  peak2 <- pseudo_voigt_single(x, a1 * ratios[2], x0 + shifts[2], sigma, gamma, eta)
  peak3 <- pseudo_voigt_single(x, a1 * ratios[3], x0 + shifts[3], sigma, gamma, eta)
  peak4 <- pseudo_voigt_single(x, a1 * ratios[4], x0 - shifts[4], sigma, gamma, eta)
  peak5 <- pseudo_voigt_single(x, a1 * ratios[5], x0 - shifts[5], sigma, gamma, eta)
  peak6 <- pseudo_voigt_single(x, a1 * ratios[6], x0 - shifts[6], sigma, gamma, eta)
  peak1 + peak2 + peak3 + peak4 + peak5 + peak6
}

#  Definition of the Pseudo-Voigt-function for seven peaks
pseudo_voigt_sevenPeaks <- function(x, a1, x0, sigma, gamma, eta, shifts, ratios) {
  peak1 <- pseudo_voigt_single(x, a1 * ratios[1], x0 + shifts[1], sigma, gamma, eta)
  peak2 <- pseudo_voigt_single(x, a1 * ratios[2], x0 + shifts[2], sigma, gamma, eta)
  peak3 <- pseudo_voigt_single(x, a1 * ratios[3], x0 + shifts[3], sigma, gamma, eta)
  peak4 <- pseudo_voigt_single(x, a1 * ratios[4], x0 + shifts[4], sigma, gamma, eta)
  peak5 <- pseudo_voigt_single(x, a1 * ratios[5], x0 - shifts[5], sigma, gamma, eta)
  peak6 <- pseudo_voigt_single(x, a1 * ratios[6], x0 - shifts[6], sigma, gamma, eta)
  peak7 <- pseudo_voigt_single(x, a1 * ratios[7], x0 - shifts[7], sigma, gamma, eta)
  peak1 + peak2 + peak3 + peak4 + peak5 + peak6 + peak7
}

# Definition of the Pseudo-Voigt-function for eight peaks
pseudo_voigt_eightPeaks <- function(x, a1, x0, sigma, gamma, eta, shifts, ratios) {
  peak1 <- pseudo_voigt_single(x, a1 * ratios[1], x0 + shifts[1], sigma, gamma, eta)
  peak2 <- pseudo_voigt_single(x, a1 * ratios[2], x0 + shifts[2], sigma, gamma, eta)
  peak3 <- pseudo_voigt_single(x, a1 * ratios[3], x0 + shifts[3], sigma, gamma, eta)
  peak4 <- pseudo_voigt_single(x, a1 * ratios[4], x0 + shifts[4], sigma, gamma, eta)
  peak5 <- pseudo_voigt_single(x, a1 * ratios[5], x0 - shifts[5], sigma, gamma, eta)
  peak6 <- pseudo_voigt_single(x, a1 * ratios[6], x0 - shifts[6], sigma, gamma, eta)
  peak7 <- pseudo_voigt_single(x, a1 * ratios[7], x0 - shifts[7], sigma, gamma, eta)
  peak8 <- pseudo_voigt_single(x, a1 * ratios[8], x0 - shifts[8], sigma, gamma, eta)
  peak1 + peak2 + peak3 + peak4 + peak5 + peak6 + peak7 + peak8
}

#----------------------------------------- baseline_correction ------------------------------------------------------

# Preprocessing: baseline_correction; setting the min y-Signal to 0
baseline_correction <- function(y) {
  min_y <- min(y)
  y_corrected <- y - min_y
  return(y_corrected)
}

#----------------------------------------- fit_signal function ------------------------------------------------------

# function to model the signal fit, based on function parameters
fit_signal <- function(sampleID, x, y, peaknumber, peak_center, tolerance, SF, peak_shifts_Hz = NULL, area_ratios = NULL, 
                       coupling_tolerance = 0.05, error_threshold = 0.35 * max(y), error_factor = 0.1, minor_signal = F, minor_signal_factor = 2, 
                       minor_signal_threshold = 5, overfit_factor = 10, overfit_cor = F) {
  
  set.seed(1997) #seed for reproducability
  
  # data preprocessing
  y_corrected <- baseline_correction(y)
  
  # calculate peaks_shifts from Hz to ppm
  if (!is.null(peak_shifts_Hz)) {
    peak_shifts <- peak_shifts_Hz / SF
  } else {
    peak_shifts <- NULL
  }
  
  # fitness-function of the GA optimizer
  if (peaknumber == "one") {
    fitness_function <- function(params) {
      a <- params[1]
      x0 <- params[2]
      sigma <- params[3]
      gamma <- params[4]
      eta <- params[5]
      fit_values <- pseudo_voigt_single(x, a, x0, sigma, gamma, eta)
      
      # Test for parameter validity
      if (any(is.na(params)) || any(is.infinite(params)) || sigma <= 0 || gamma <= 0 || eta < 0 || eta > 1) {
        cat("Invalid parameter: ", params, "\n")
        return(NA)
      }
      
      # Error calculation
      ERROR <- y_corrected - fit_values
      
      if (overfit_cor == T){
        # If fit at chemical shift value is exceeding spectral line (fit_values > y_corrected)
        overfit <- fit_values > y_corrected
        ERROR[overfit] <- ERROR[overfit] * overfit_factor  # Emphasizing error of overfitted chemical shifts
      }
      
      if (minor_signal) {
        #emphazise errorcalculation for minor signals
        igr <- fit_values > max(y_corrected) * minor_signal_threshold
        ERROR[igr] <- ERROR[igr] * minor_signal_factor
      }
      
      #Emphazising the smaller errors (error_threshold)
      im1 <- (ERROR < error_threshold)
      im2 <- (ERROR >= error_threshold)
      
      #summarize sqaured errors
      D <- sum(ERROR[im1]^2) * error_factor + sum(ERROR[im2]^2)
      #calculate mean error
      D <- D / length(ERROR)
      #return negative mean error
      -D  # fitness-function return has to be negative D, as GA optimizer is minimizing
    }
    
    # fitness-function return has to be negative D, as GA optimizer is minimizing
    lower_bounds <- c(0, peak_center - tolerance, 0, 0, 0)
    upper_bounds <- c(max(y_corrected), peak_center + tolerance, 2 / SF, 2 / SF, 1)
    
  } 
  else if (peaknumber == "two") {
    fitness_function <- function(params) {
      a1 <- params[1]
      x0 <- params[2]
      sigma <- params[3]
      gamma <- params[4]
      eta <- params[5]
      shift1 <- params[6]  # changed shift for peak 1
      shift2 <- params[7]  # changed shift for peak 2
      
      # Test for parameter validity
      if (any(is.na(params)) || any(is.infinite(params)) || sigma <= 0 || gamma <= 0 || eta < 0 || eta > 1) {
        cat("Invalid parameter: ", params, "\n")
        return(NA)
      }
      
      shifts_adjusted <- c(shift1, shift2)
      
      fit_values <- pseudo_voigt_twoPeaks(x, a1, x0, sigma, gamma, eta, shifts_adjusted, area_ratios)
      
      # Error calculation
      ERROR <- y_corrected - fit_values
      
      if (overfit_cor == T){
        # If fit at chemical shift value is exceeding spectral line (fit_values > y_corrected)
        overfit <- fit_values > y_corrected
        ERROR[overfit] <- ERROR[overfit] * overfit_factor  # Emphasizing error of overfitted chemical shifts
      }
      
      if (minor_signal) {
        #emphazise errorcalculation for minor signals
        igr <- fit_values > max(y_corrected) * minor_signal_threshold
        ERROR[igr] <- ERROR[igr] * minor_signal_factor
      }
      
      im1 <- (ERROR < error_threshold)
      im2 <- (ERROR >= error_threshold)
      
      D <- sum(ERROR[im1]^2) * error_factor + sum(ERROR[im2]^2)
      D <- D / length(ERROR)
      
      -D  # fitness-function return has to be negative D, as GA optimizer is minimizing
    }
    
    # fitness-function return has to be negative D, as GA optimizer is minimizing
    lower_bounds <- c(0, peak_center - tolerance, 0, 0, 0, peak_shifts * (1 - coupling_tolerance))
    upper_bounds <- c(max(y_corrected), peak_center + tolerance, 2 / SF, 2 / SF, 1, peak_shifts * (1 + coupling_tolerance))
    
  } 
  else if (peaknumber == "three") {
    fitness_function <- function(params) {
      a1 <- params[1]
      x0 <- params[2]
      sigma <- params[3]
      gamma <- params[4]
      eta <- params[5]
      shift1 <- params[6]  # changed shift for peak 1
      shift2 <- params[7]  # changed shift for peak 2
      shift3 <- params[8]  # changed shift for peak 3
      
      # Test for parameter validity
      if (any(is.na(params)) || any(is.infinite(params)) || sigma <= 0 || gamma <= 0 || eta < 0 || eta > 1) {
        cat("Invalid parameter: ", params, "\n")
        return(NA)
      }
      
      shifts_adjusted <- c(shift1, shift2, shift3)
      
      fit_values <- pseudo_voigt_threePeaks(x, a1, x0, sigma, gamma, eta, shifts_adjusted, area_ratios)
      
      # Error calculation
      ERROR <- y_corrected - fit_values
      
      if (overfit_cor == T){
        # If fit at chemical shift value is exceeding spectral line (fit_values > y_corrected)
        overfit <- fit_values > y_corrected
        ERROR[overfit] <- ERROR[overfit] * overfit_factor  # Emphasizing error of overfitted chemical shifts
      }
      
      if (minor_signal) {
        #emphazise errorcalculation for minor signals
        igr <- fit_values > max(y_corrected) * minor_signal_threshold
        ERROR[igr] <- ERROR[igr] * minor_signal_factor
      }
      
      im1 <- (ERROR < error_threshold)
      im2 <- (ERROR >= error_threshold)
      
      D <- sum(ERROR[im1]^2) * error_factor + sum(ERROR[im2]^2)
      D <- D / length(ERROR)
      
      -D  # fitness-function return has to be negative D, as GA optimizer is minimizing
    }
    
    # fitness-function return has to be negative D, as GA optimizer is minimizing
    lower_bounds <- c(0, peak_center - tolerance, 0, 0, 0, peak_shifts * (1 - coupling_tolerance))
    upper_bounds <- c(max(y_corrected), peak_center + tolerance, 2 / SF, 2 / SF, 1, peak_shifts * (1 + coupling_tolerance))
    
  } 
  else if (peaknumber == "four") {
    fitness_function <- function(params) {
      a1 <- params[1]
      x0 <- params[2]
      sigma <- params[3]
      gamma <- params[4]
      eta <- params[5]
      shift1 <- params[6]  # changed shift for peak 1
      shift2 <- params[7]  # changed shift for peak 2
      shift3 <- params[8]  # changed shift for peak 3
      shift4 <- params[9]  # changed shift for peak 4
      
      # Test for parameter validity
      if (any(is.na(params)) || any(is.infinite(params)) || sigma <= 0 || gamma <= 0 || eta < 0 || eta > 1) {
        cat("Invalid parameter: ", params, "\n")
        return(NA)
      }
      
      shifts_adjusted <- c(shift1, shift2, shift3, shift4)
      
      fit_values <- pseudo_voigt_fourPeaks(x, a1, x0, sigma, gamma, eta, shifts_adjusted, area_ratios)
      
      # Error calculation
      ERROR <- y_corrected - fit_values
      
      if (overfit_cor == T){
        # If fit at chemical shift value is exceeding spectral line (fit_values > y_corrected)
        overfit <- fit_values > y_corrected
        ERROR[overfit] <- ERROR[overfit] * overfit_factor  # Emphasizing error of overfitted chemical shifts
      }
      
      if (minor_signal) {
        #emphazise errorcalculation for minor signals
        igr <- fit_values > max(y_corrected) * minor_signal_threshold
        ERROR[igr] <- ERROR[igr] * minor_signal_factor
      }
      
      im1 <- (ERROR < error_threshold)
      im2 <- (ERROR >= error_threshold)
      
      D <- sum(ERROR[im1]^2) * error_factor + sum(ERROR[im2]^2)
      D <- D / length(ERROR)
      
      -D  # fitness-function return has to be negative D, as GA optimizer is minimizing
    }
    
    # fitness-function return has to be negative D, as GA optimizer is minimizing
    lower_bounds <- c(0, peak_center - tolerance, 0, 0, 0, peak_shifts * (1 - coupling_tolerance))
    upper_bounds <- c(max(y_corrected), peak_center + tolerance, 2 / SF, 2 / SF, 1, peak_shifts * (1 + coupling_tolerance))
    
  } 
  else if (peaknumber == "five") {
    fitness_function <- function(params) {
      a1 <- params[1]
      x0 <- params[2]
      sigma <- params[3]
      gamma <- params[4]
      eta <- params[5]
      shift1 <- params[6]  # changed shift for peak 1
      shift2 <- params[7]  # changed shift for peak 2
      shift3 <- params[8]  # changed shift for peak 3
      shift4 <- params[9]  # changed shift for peak 4
      shift5 <- params[10] # changed shift for peak 5
      
      # Test for parameter validity
      if (any(is.na(params)) || any(is.infinite(params)) || sigma <= 0 || gamma <= 0 || eta < 0 || eta > 1) {
        cat("Invalid parameter: ", params, "\n")
        return(NA)
      }
      
      shifts_adjusted <- c(shift1, shift2, shift3, shift4, shift5)
      
      fit_values <- pseudo_voigt_fivePeaks(x, a1, x0, sigma, gamma, eta, shifts_adjusted, area_ratios)
      
      # Error calculation
      ERROR <- y_corrected - fit_values
      
      if (overfit_cor == T){
        # If fit at chemical shift value is exceeding spectral line (fit_values > y_corrected)
        overfit <- fit_values > y_corrected
        ERROR[overfit] <- ERROR[overfit] * overfit_factor  # Emphasizing error of overfitted chemical shifts
      }
      
      if (minor_signal) {
        #emphazise errorcalculation for minor signals
        igr <- fit_values > max(y_corrected) * minor_signal_threshold
        ERROR[igr] <- ERROR[igr] * minor_signal_factor
      }
      
      im1 <- (ERROR < error_threshold)
      im2 <- (ERROR >= error_threshold)
      
      D <- sum(ERROR[im1]^2) * error_factor + sum(ERROR[im2]^2)
      D <- D / length(ERROR)
      
      -D  # fitness-function return has to be negative D, as GA optimizer is minimizing
    }
    
    # fitness-function return has to be negative D, as GA optimizer is minimizing
    lower_bounds <- c(0, peak_center - tolerance, 0, 0, 0, peak_shifts * (1 - coupling_tolerance))
    upper_bounds <- c(max(y_corrected), peak_center + tolerance, 2 / SF, 2 / SF, 1, peak_shifts * (1 + coupling_tolerance))
    
  }
  else if (peaknumber == "six") {
    fitness_function <- function(params) {
      a1 <- params[1]
      x0 <- params[2]
      sigma <- params[3]
      gamma <- params[4]
      eta <- params[5]
      shift1 <- params[6]  # changed shift for peak 1
      shift2 <- params[7]  # changed shift for peak 2
      shift3 <- params[8]  # changed shift for peak 3
      shift4 <- params[9]  # changed shift for peak 4
      shift5 <- params[10] # changed shift for peak 5
      shift6 <- params[11] # changed shift for peak 6
      
      # Test for parameter validity
      if (any(is.na(params)) || any(is.infinite(params)) || sigma <= 0 || gamma <= 0 || eta < 0 || eta > 1) {
        cat("Invalid parameter: ", params, "\n")
        return(NA)
      }
      
      shifts_adjusted <- c(shift1, shift2, shift3, shift4, shift5, shift6)
      
      fit_values <- pseudo_voigt_fourPeaks(x, a1, x0, sigma, gamma, eta, shifts_adjusted, area_ratios)
      
      # Error calculation
      ERROR <- y_corrected - fit_values
      
      if (overfit_cor == T){
        # If fit at chemical shift value is exceeding spectral line (fit_values > y_corrected)
        overfit <- fit_values > y_corrected
        ERROR[overfit] <- ERROR[overfit] * overfit_factor  # Emphasizing error of overfitted chemical shifts
      }
      
      if (minor_signal) {
        #emphazise errorcalculation for minor signals
        igr <- fit_values > max(y_corrected) * minor_signal_threshold
        ERROR[igr] <- ERROR[igr] * minor_signal_factor
      }
      
      im1 <- (ERROR < error_threshold)
      im2 <- (ERROR >= error_threshold)
      
      D <- sum(ERROR[im1]^2) * error_factor + sum(ERROR[im2]^2)
      D <- D / length(ERROR)
      
      -D  # fitness-function return has to be negative D, as GA optimizer is minimizing
    }
    
    # fitness-function return has to be negative D, as GA optimizer is minimizing
    lower_bounds <- c(0, peak_center - tolerance, 0, 0, 0, peak_shifts * (1 - coupling_tolerance))
    upper_bounds <- c(max(y_corrected), peak_center + tolerance, 2 / SF, 2 / SF, 1, peak_shifts * (1 + coupling_tolerance))
    
  }
  else if (peaknumber == "seven") {
    fitness_function <- function(params) {
      a1 <- params[1]
      x0 <- params[2]
      sigma <- params[3]
      gamma <- params[4]
      eta <- params[5]
      shift1 <- params[6]  # changed shift for peak 1
      shift2 <- params[7]  # changed shift for peak 2
      shift3 <- params[8]  # changed shift for peak 3
      shift4 <- params[9]  # changed shift for peak 4
      shift5 <- params[10] # changed shift for peak 5
      shift6 <- params[11] # changed shift for peak 6
      shift7 <- params[12] # changed shift for peak 7
      
      # Test for parameter validity
      if (any(is.na(params)) || any(is.infinite(params)) || sigma <= 0 || gamma <= 0 || eta < 0 || eta > 1) {
        cat("Invalid parameter: ", params, "\n")
        return(NA)
      }
      
      shifts_adjusted <- c(shift1, shift2, shift3, shift4, shift5, shift6, shift7)
      
      fit_values <- pseudo_voigt_fourPeaks(x, a1, x0, sigma, gamma, eta, shifts_adjusted, area_ratios)
      
      # Error calculation
      ERROR <- y_corrected - fit_values
      
      if (overfit_cor == T){
        # If fit at chemical shift value is exceeding spectral line (fit_values > y_corrected)
        overfit <- fit_values > y_corrected
        ERROR[overfit] <- ERROR[overfit] * overfit_factor  # Emphasizing error of overfitted chemical shifts
      }
      
      if (minor_signal) {
        #emphazise errorcalculation for minor signals
        igr <- fit_values > max(y_corrected) * minor_signal_threshold
        ERROR[igr] <- ERROR[igr] * minor_signal_factor
      }
      
      im1 <- (ERROR < error_threshold)
      im2 <- (ERROR >= error_threshold)
      
      D <- sum(ERROR[im1]^2) * error_factor + sum(ERROR[im2]^2)
      D <- D / length(ERROR)
      
      -D  # fitness-function return has to be negative D, as GA optimizer is minimizing
    }
    
    # fitness-function return has to be negative D, as GA optimizer is minimizing
    lower_bounds <- c(0, peak_center - tolerance, 0, 0, 0, peak_shifts * (1 - coupling_tolerance))
    upper_bounds <- c(max(y_corrected), peak_center + tolerance, 2 / SF, 2 / SF, 1, peak_shifts * (1 + coupling_tolerance))
    
  }
  else if (peaknumber == "eight") {
    fitness_function <- function(params) {
      a1 <- params[1]
      x0 <- params[2]
      sigma <- params[3]
      gamma <- params[4]
      eta <- params[5]
      shift1 <- params[6]  # changed shift for peak 1
      shift2 <- params[7]  # changed shift for peak 2
      shift3 <- params[8]  # changed shift for peak 3
      shift4 <- params[9]  # changed shift for peak 4
      shift5 <- params[10] # changed shift for peak 5
      shift6 <- params[11] # changed shift for peak 6
      shift7 <- params[12] # changed shift for peak 7
      shift8 <- params[13] # changed shift for peak 8
      
      # Test for parameter validity
      if (any(is.na(params)) || any(is.infinite(params)) || sigma <= 0 || gamma <= 0 || eta < 0 || eta > 1) {
        cat("Invalid parameter: ", params, "\n")
        return(NA)
      }
      
      shifts_adjusted <- c(shift1, shift2, shift3, shift4, shift5, shift6, shift7, shift8)
      
      fit_values <- pseudo_voigt_fourPeaks(x, a1, x0, sigma, gamma, eta, shifts_adjusted, area_ratios)
      
      # Error calculation
      ERROR <- y_corrected - fit_values
      
      if (overfit_cor == T){
        # If fit at chemical shift value is exceeding spectral line (fit_values > y_corrected)
        overfit <- fit_values > y_corrected
        ERROR[overfit] <- ERROR[overfit] * overfit_factor  # Emphasizing error of overfitted chemical shifts
      }
      
      if (minor_signal) {
        #emphazise errorcalculation for minor signals
        igr <- fit_values > max(y_corrected) * minor_signal_threshold
        ERROR[igr] <- ERROR[igr] * minor_signal_factor
      }
      
      im1 <- (ERROR < error_threshold)
      im2 <- (ERROR >= error_threshold)
      
      D <- sum(ERROR[im1]^2) * error_factor + sum(ERROR[im2]^2)
      D <- D / length(ERROR)
      
      -D  # fitness-function return has to be negative D, as GA optimizer is minimizing
    }
    
    # fitness-function return has to be negative D, as GA optimizer is minimizing
    lower_bounds <- c(0, peak_center - tolerance, 0, 0, 0, peak_shifts * (1 - coupling_tolerance))
    upper_bounds <- c(max(y_corrected), peak_center + tolerance, 2 / SF, 2 / SF, 1, peak_shifts * (1 + coupling_tolerance))
    
  }
  else {
    stop("Invalid entry for peaknumber. Has to be 'one', 'two', 'three', 'four', 'five', 'six', 'seven' oder 'eight'.")
  }
  
  # GA optimizing
  ga_result <- ga(type = "real-valued", 
                  fitness = fitness_function, 
                  lower = lower_bounds, 
                  upper = upper_bounds, 
                  popSize = 50, 
                  maxiter = 1000, 
                  run = 150,
                  monitor = FALSE # no print in console
  ) 
  
  best_params <- ga_result@solution
  integral_value <- NA  # Initializing
  mse <- NA  # Initializing mse
  
  if (peaknumber == "one") {
    a <- best_params[1]
    x0 <- best_params[2]
    sigma <- best_params[3]
    gamma <- best_params[4]
    eta <- best_params[5]
    fit_values <- pseudo_voigt_single(x, a, x0, sigma, gamma, eta)
    
    # Calculation integral of Pseudo-Voigt-function for one peak
    pseudo_voigt_integral_singlet <- function(a, sigma, gamma, eta) {
      lorentz_integral <- pi * a * gamma
      gauss_integral <- a * sigma * sqrt(2 * pi)
      eta * lorentz_integral + (1 - eta) * gauss_integral
    }
    integral_value <- pseudo_voigt_integral_singlet(a, sigma, gamma, eta)
    
  } else if (peaknumber == "two") {
    a1 <- best_params[1]
    x0 <- best_params[2]
    sigma <- best_params[3]
    gamma <- best_params[4]
    eta <- best_params[5]
    shift1 <- best_params[6]
    shift2 <- best_params[7]
    fit_values <- pseudo_voigt_twoPeaks(x, a1, x0, sigma, gamma, eta, c(shift1, shift2), area_ratios)
    
    # Calculation integral of Pseudo-Voigt-function for for two peaks
    pseudo_voigt_integral_twoPeaks <- function(a1, sigma, gamma, eta, ratios) {
      lorentz_integral <- pi * a1 * gamma
      gauss_integral <- a1 * sigma * sqrt(2 * pi)
      single_integral <- eta * lorentz_integral + (1 - eta) * gauss_integral
      sum(single_integral * ratios)
    }
    integral_value <- pseudo_voigt_integral_twoPeaks(a1, sigma, gamma, eta, area_ratios)
    
  } else if (peaknumber == "three") {
    a1 <- best_params[1]
    x0 <- best_params[2]
    sigma <- best_params[3]
    gamma <- best_params[4]
    eta <- best_params[5]
    shift1 <- best_params[6]
    shift2 <- best_params[7]
    shift3 <- best_params[8]
    fit_values <- pseudo_voigt_threePeaks(x, a1, x0, sigma, gamma, eta, c(shift1, shift2, shift3), area_ratios)
    
    # Calculation integral of Pseudo-Voigt-function for three peaks
    pseudo_voigt_integral_threePeaks <- function(a1, sigma, gamma, eta, ratios) {
      lorentz_integral <- pi * a1 * gamma
      gauss_integral <- a1 * sigma * sqrt(2 * pi)
      single_integral <- eta * lorentz_integral + (1 - eta) * gauss_integral
      sum(single_integral * ratios)
    }
    integral_value <- pseudo_voigt_integral_threePeaks(a1, sigma, gamma, eta, area_ratios)
    
  } else if (peaknumber == "four"){
    a1 <- best_params[1]
    x0 <- best_params[2]
    sigma <- best_params[3]
    gamma <- best_params[4]
    eta <- best_params[5]
    shift1 <- best_params[6]
    shift2 <- best_params[7]
    shift3 <- best_params[8]
    shift4 <- best_params[9]
    fit_values <- pseudo_voigt_triplet_fourPeaks(x, a1, x0, sigma, gamma, eta, c(shift1, shift2, shift3, shift4), area_ratios)
    
    # Calculation integral of Pseudo-Voigt-function for four peaks
    pseudo_voigt_integral_fourPeaks <- function(a1, sigma, gamma, eta, ratios) {
      lorentz_integral <- pi * a1 * gamma
      gauss_integral <- a1 * sigma * sqrt(2 * pi)
      single_integral <- eta * lorentz_integral + (1 - eta) * gauss_integral
      sum(single_integral * ratios)
    }
    integral_value <- pseudo_voigt_integral_fourPeaks(a1, sigma, gamma, eta, area_ratios)
  } else if (peaknumber == "five"){
    a1 <- best_params[1]
    x0 <- best_params[2]
    sigma <- best_params[3]
    gamma <- best_params[4]
    eta <- best_params[5]
    shift1 <- best_params[6]
    shift2 <- best_params[7]
    shift3 <- best_params[8]
    shift4 <- best_params[9]
    shift5 <- best_params[10]
    fit_values <- pseudo_voigt_fivePeaks(x, a1, x0, sigma, gamma, eta, c(shift1, shift2, shift3, shift4, shift5), area_ratios)
    
    # Calculation integral of Pseudo-Voigt-function for five peaks
    pseudo_voigt_integral_fivePeaks <- function(a1, sigma, gamma, eta, ratios) {
      lorentz_integral <- pi * a1 * gamma
      gauss_integral <- a1 * sigma * sqrt(2 * pi)
      single_integral <- eta * lorentz_integral + (1 - eta) * gauss_integral
      sum(single_integral * ratios)
    }
    integral_value <- pseudo_voigt_integral_fivePeaks(a1, sigma, gamma, eta, area_ratios)
  } else if (peaknumber == "six"){
    a1 <- best_params[1]
    x0 <- best_params[2]
    sigma <- best_params[3]
    gamma <- best_params[4]
    eta <- best_params[5]
    shift1 <- best_params[6]
    shift2 <- best_params[7]
    shift3 <- best_params[8]
    shift4 <- best_params[9]
    shift5 <- best_params[10]
    shift6 <- best_params[11]
    fit_values <- pseudo_voigt_sixPeaks(x, a1, x0, sigma, gamma, eta, c(shift1, shift2, shift3, shift4, shift5, shift6), area_ratios)
    
    # Calculation integral of Pseudo-Voigt-function for six peaks
    pseudo_voigt_integral_sixPeaks <- function(a1, sigma, gamma, eta, ratios) {
      lorentz_integral <- pi * a1 * gamma
      gauss_integral <- a1 * sigma * sqrt(2 * pi)
      single_integral <- eta * lorentz_integral + (1 - eta) * gauss_integral
      sum(single_integral * ratios)
    }
    integral_value <- pseudo_voigt_integral_sixPeaks(a1, sigma, gamma, eta, area_ratios)
  }
  else if (peaknumber == "seven"){
    a1 <- best_params[1]
    x0 <- best_params[2]
    sigma <- best_params[3]
    gamma <- best_params[4]
    eta <- best_params[5]
    shift1 <- best_params[6]
    shift2 <- best_params[7]
    shift3 <- best_params[8]
    shift4 <- best_params[9]
    shift5 <- best_params[10]
    shift6 <- best_params[11]
    shift7 <- best_params[12]
    fit_values <- pseudo_voigt_sevenPeaks(x, a1, x0, sigma, gamma, eta, c(shift1, shift2, shift3, shift4, shift5, shift6, shift7), area_ratios)
    
    # Calculation integral of Pseudo-Voigt-function for seven peaks
    pseudo_voigt_integral_sevenPeaks <- function(a1, sigma, gamma, eta, ratios) {
      lorentz_integral <- pi * a1 * gamma
      gauss_integral <- a1 * sigma * sqrt(2 * pi)
      single_integral <- eta * lorentz_integral + (1 - eta) * gauss_integral
      sum(single_integral * ratios)
    }
    integral_value <- pseudo_voigt_integral_sevenPeaks(a1, sigma, gamma, eta, area_ratios)
  } else {
    a1 <- best_params[1]
    x0 <- best_params[2]
    sigma <- best_params[3]
    gamma <- best_params[4]
    eta <- best_params[5]
    shift1 <- best_params[6]
    shift2 <- best_params[7]
    shift3 <- best_params[8]
    shift4 <- best_params[9]
    shift5 <- best_params[10]
    shift6 <- best_params[11]
    shift7 <- best_params[12]
    shift8 <- best_params[13]
    fit_values <- pseudo_voigt_eightPeaks(x, a1, x0, sigma, gamma, eta, c(shift1, shift2, shift3, shift4, shift5, shift6, shift7, shift8), area_ratios)
    
    # Calculation integral of Pseudo-Voigt-function for eight peaks
    pseudo_voigt_integral_eightPeaks <- function(a1, sigma, gamma, eta, ratios) {
      lorentz_integral <- pi * a1 * gamma
      gauss_integral <- a1 * sigma * sqrt(2 * pi)
      single_integral <- eta * lorentz_integral + (1 - eta) * gauss_integral
      sum(single_integral * ratios)
    }
    integral_value <- pseudo_voigt_integral_eightPeaks(a1, sigma, gamma, eta, area_ratios)
  }
  
  # residuals calculation
  residuals <- y_corrected - fit_values
  
  # calculation of mean squared error (MSE)
  mse <- mean(residuals^2)
  
  # data frame for ggplot2 visualisation
  data <- data.frame(x = x, y_corrected = y_corrected, fit_values = fit_values, residuals = residuals)
  
  # print the spectral data, the fit and the r miesiduals using ggplot2
  p <- ggplot(data, aes(x = x)) +
    geom_line(aes(y = y_corrected, color = "spectral data"), linewidth = 0.75) +
    geom_line(aes(y = residuals, color = "residuals"), linewidth = 0.75) +
    geom_ribbon(aes(ymin = 0, ymax = fit_values, fill = "fit"), alpha = 0.3) +
    scale_x_reverse() +
    scale_color_manual("", values = c("spectral data" = "black", "residuals" = "red")) +
    scale_fill_manual("", values = c("fit" = "blue")) +
    labs(title = paste(sampleID, "-",peaknumber, "Peak", "- Pseudo-Voigt-Fit (GA Optimizer)"),
         x = "chemical shift (ppm)",
         y = "signal intensity") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.text = element_text(size = 10))
  
  print(p)
  
  # Print the patrameters of the best fit and it's integral
  if (peaknumber == "one") {
    cat("best fit for one peak:\n")
    cat("a:", a, "\n")
    cat("x0:", x0, "\n")
    cat("sigma:", sigma, "\n")
    cat("gamma:", gamma, "\n")
    cat("eta:", eta, "\n")
  } else if (peaknumber == "two") {
    cat("best fit for two peaks:\n")
    cat("a1:", a1, "\n")
    cat("x0:", x0, "\n")
    cat("sigma:", sigma, "\n")
    cat("gamma:", gamma, "\n")
    cat("eta:", eta, "\n")
    cat("peak positions (ppm):", c(shift1, shift2), "\n")
    cat("area ratios:", area_ratios, "\n")
  } else if (peaknumber == "three") {
    cat("best fit for three peaks:\n")
    cat("a1:", a1, "\n")
    cat("x0:", x0, "\n")
    cat("sigma:", sigma, "\n")
    cat("gamma:", gamma, "\n")
    cat("eta:", eta, "\n")
    cat("peak positions (ppm):", c(shift1, shift2, shift3), "\n")
    cat("area ratios:", area_ratios, "\n")
  } else if (peaknumber == "four"){
    cat("best fit for four peaks:\n")
    cat("a1:", a1, "\n")
    cat("x0:", x0, "\n")
    cat("sigma:", sigma, "\n")
    cat("gamma:", gamma, "\n")
    cat("eta:", eta, "\n")
    cat("peak positions (ppm):", c(shift1, shift2, shift3, shift4), "\n")
    cat("area ratios:", area_ratios, "\n")
  } else if (peaknumber == "five"){
    cat("best fit for five peaks:\n")
    cat("a1:", a1, "\n")
    cat("x0:", x0, "\n")
    cat("sigma:", sigma, "\n")
    cat("gamma:", gamma, "\n")
    cat("eta:", eta, "\n")
    cat("peak positions (ppm):", c(shift1, shift2, shift3, shift4, shift5), "\n")
    cat("area ratios:", area_ratios, "\n")
  }else if (peaknumber == "six"){
    cat("best fit for six peaks:\n")
    cat("a1:", a1, "\n")
    cat("x0:", x0, "\n")
    cat("sigma:", sigma, "\n")
    cat("gamma:", gamma, "\n")
    cat("eta:", eta, "\n")
    cat("peak positions (ppm):", c(shift1, shift2, shift3, shift4, shift5, shift6), "\n")
    cat("area ratios:", area_ratios, "\n")
  }else if (peaknumber == "seven"){
    cat("best fit for seven peaks:\n")
    cat("a1:", a1, "\n")
    cat("x0:", x0, "\n")
    cat("sigma:", sigma, "\n")
    cat("gamma:", gamma, "\n")
    cat("eta:", eta, "\n")
    cat("peak positions (ppm):", c(shift1, shift2, shift3, shift4, shift5, shift6, shift7), "\n")
    cat("area ratios:", area_ratios, "\n")
  }else {
    cat("best fit for eight peaks:\n")
    cat("a1:", a1, "\n")
    cat("x0:", x0, "\n")
    cat("sigma:", sigma, "\n")
    cat("gamma:", gamma, "\n")
    cat("eta:", eta, "\n")
    cat("peak positions (ppm):", c(shift1, shift2, shift3, shift4, shift5, shift6, shift7, shift8), "\n")
    cat("area ratios:", area_ratios, "\n")
  }
  cat("Integral of the Pseudo-Voigt-Fits:", integral_value, "\n")
  cat("Mean Squared Error (MSE) of the Fits:", mse, "\n")
  
  # returning the fit results as object
  fit_result <- list(
    best_params = best_params,
    fit_values = fit_values,
    residuals = residuals,
    integral_value = integral_value,
    mse = mse,
    plot = p
  )
  
  return(fit_result)
}
