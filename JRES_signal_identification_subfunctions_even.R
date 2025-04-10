# subscripts for selecting even numbered signals

# Load necessary R-Packages
library(dplyr)
library(purrr)
library(tidyr)

# Identification of the center peaks and removing uneven numbered multiplicities --------------

cut_odd_couplings <- function(JRESPeaklist,
                                 F1shifttol = 0.0005,
                                 F2shifttol = 0.001) {
  middle_signal <- which((
    JRESPeaklist$F1ppm > (-F1shifttol) &
      JRESPeaklist$F1ppm < F1shifttol
  ),
  arr.ind = T)
  if (length(middle_signal) > 0) {
    # Search for potentially coupling peaks of middle peaks
    row_couplesignals <- c()
    
    for (n in 1:length(middle_signal)) {
      # Extract the corresponding values from JRESPeakList
      F2value <- JRESPeaklist[middle_signal[n], "F2ppm"]
      
      # Search for coupling partners to the middle signals
      couplesignals <- c(which((
        JRESPeaklist$F2ppm > (F2value - F2shifttol) &
          JRESPeaklist$F2ppm < (F2value + F2shifttol)
      ), arr.ind = T))
      
      row_couplesignals <- c(row_couplesignals, couplesignals)
    }
    
    even_couplings <- JRESPeaklist[-row_couplesignals, ]
  }
  if (length(middle_signal) == 0) {
    even_couplings <- JRESPeaklist
  }
  return(even_couplings)
}


# Structuring and preselection of possible couplings --------

filter_even_couplings <- function(even_couplings, F2shifttol = 0.001) {
  
  # Create an empty data frame for the difference information
  diff_even_coupling <- data.frame(
    "row_Peak1" = numeric(),
    "row_Peak2" = numeric(),
    Differenceppm = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (nrow(even_couplings) >= 2) {
    
    # Breakdown of positive and negative values
    position_positive <- even_couplings$F1ppm > 0
    position_negative <- even_couplings$F1ppm < 0
    
    # Use of the original rows
    rownames_positive <- rownames(even_couplings)[position_positive]
    rownames_negative <- rownames(even_couplings)[position_negative]
    
    # Positive and negative coupling indices
    pos_indices <- as.numeric(rownames_positive)
    neg_indices <- as.numeric(rownames_negative)
    
    # Create the possible combinations of peaks
    combinations <- expand.grid(row_Peak1 = pos_indices, row_Peak2 = neg_indices)
    
    # Calculate the differences and add to the empty dataframe
    for (n in 1:nrow(combinations)) {
      rowPeak1 <- as.numeric(combinations[n, 1])
      rowPeak2 <- as.numeric(combinations[n, 2])
      difference <- abs(as.numeric(JRESPeaklist$F1ppm[rowPeak1]) - as.numeric(JRESPeaklist$F1ppm[rowPeak2]))
      
      diff_even_coupling <- rbind(
        diff_even_coupling,
        data.frame(
          "row_Peak1" = rowPeak1,
          "row_Peak2" = rowPeak2,
          Differenceppm = difference
        )
      )
    }
    
    # Exclusion of couplings that are not in the same position
    diff_delete_rows <- rep(FALSE, nrow(diff_even_coupling))
    
    # Iteration over the rows of diff_even_coupling
    for (i in 1:nrow(diff_even_coupling)) {
      row1 <- diff_even_coupling[i, 1]
      row2 <- diff_even_coupling[i, 2]
      
      # Extract the corresponding values from JRESPeakList
      value1 <- JRESPeaklist[row1, "F2ppm"]
      value2 <- JRESPeaklist[row2, "F2ppm"]
      
      # Matching and marking couplings to be removed
      if (abs(value1 - value2) > F2shifttol) {
        diff_delete_rows[i] <- TRUE
      }
    }
    
    # Remove the couplings that are not in the same position
    diff_even_coupling <- diff_even_coupling[!diff_delete_rows, ]
  }
  return(diff_even_coupling)
}


# Selection of the desired even numbered multiplicity ---------------------

find_even_coupled_signals <- function(diff_even_coupling, multiplicity, F2shifttol = 0.001) {
 
  if (nrow(diff_even_coupling) >= 1){
   # Sort the determined couplings into F2 groups of the same position
  # Add the chemical shift to the data
  diff_even_coupling$F2ppm_Peak1 <- JRESPeaklist$F2ppm[diff_even_coupling$row_Peak1]
  diff_even_coupling$F2ppm_Peak2 <- JRESPeaklist$F2ppm[diff_even_coupling$row_Peak2]
  
  # Empty lists for saving groups and IDs
  groups <- list()
  group_ids <- integer(nrow(diff_even_coupling))  # Initialize a column for group IDs
  group_id <- 1  # Start of the group ID
  
  # Loop through each row of the data set
  for (i in 1:nrow(diff_even_coupling)) {
    current_peak1 <- diff_even_coupling$F2ppm_Peak1[i]
    current_peak2 <- diff_even_coupling$F2ppm_Peak2[i]
    
    # Check whether Peak1 or Peak2 are already in an existing group
    in_group <- FALSE
    for (j in seq_along(groups)) {
      if (any(abs(groups[[j]] - current_peak1) <= F2shifttol) ||
          any(abs(groups[[j]] - current_peak2) <= F2shifttol)) {
        
        # Add the peaks to the existing group
        groups[[j]] <- c(groups[[j]], current_peak1, current_peak2)
        group_ids[i] <- j  # Assigns the group ID to the line
        in_group <- TRUE
        break
      }
    }
    
    # If the peak does not yet belong to a group, create a new group
    if (!in_group) {
      groups[[group_id]] <- c(current_peak1, current_peak2)
      group_ids[i] <- group_id  # Assign new group ID
      group_id <- group_id + 1
    }
  }
  
  # Assign the group IDs to the source data set
  diff_even_coupling$GroupID <- group_ids
  
  # Calculate the number of each group
  group_counts <- table(diff_even_coupling$GroupID)
  
  # Filter the groups that occur exactly X times (number of possible couplings in even multiplicities: x = (multiplicity/2)^2)
  groups_with_x <- as.numeric(names(group_counts[group_counts == (multiplicity/2)^2]))
  
  # Create a new dataframe that only contains the rows whose GroupID occurs X times
  multiplicity_checked <- diff_even_coupling[diff_even_coupling$GroupID %in% groups_with_x, ]
 
  }
  if (nrow(diff_even_coupling) < 1) {
    multiplicity_checked <- diff_even_coupling %>% mutate(GroupID = NA)
  }
  
  return(multiplicity_checked)
} 

# Checking the coupling information --------------------------------------

check_even_coupling <- function(multiplicity_checked,
                                coupleconst,
                                multiplicity,
                                coupletol = 0.01) {
  diff_checked <- c()
  for (m in 1:(length(coupleconst) / 2)) {
    # Check whether couplings of the left half of the vector are fulfilled
    diff_checked_intern <- which((
      multiplicity_checked$Differenceppm > coupleconst[m] * 2 - coupletol * coupleconst[m] *
        2 &
        multiplicity_checked$Differenceppm < coupleconst[m] * 2 + coupletol * coupleconst[m] *
        2
    ),
    arr.ind = T
    )
    
    # Add to result vector
    diff_checked <- c(diff_checked, diff_checked_intern)
  }
 
   # Preparation of the coupling peaks
  couple_checked <- multiplicity_checked[diff_checked, ]  
  
  # Calculate the number of each group
  group_counts <- table(couple_checked$GroupID)
  
  # Filter the groups that occur exactly X times (number of possible linkages in even multiplicities after linkage check: x = (multiplicity/2))
  groups_with_x <- as.numeric(names(group_counts[group_counts == (multiplicity/2)]))
  
  # Create a new dataframe that only contains the rows whose GroupID occurs X times
  couple_checked <- couple_checked[couple_checked$GroupID %in% groups_with_x, ]
  
  return(couple_checked)
}


# Checking height ratios -----------------------------------------------
# Caveat: Signal multiplicity is derived from vector heightratio
# Height ratio check implemented as optional; to activate check_height = TRUE
check_height_even_signal <- function(couple_checked,
                                     check_height = FALSE,
                                     heightratio,
                                     heighttol = 0.1) {

  
  if (check_height == TRUE &&
      nrow(couple_checked) >= 1)  {
    # Defining the boundaries
    upper_bound <- heightratio + heightratio * heighttol
    lower_bound <- heightratio - heightratio * heighttol
    
    # Create a list with the position vectors of all peaks per signal
    list_of_vectors <- couple_checked %>%
      group_by(GroupID) %>%
      group_split() %>%
      map( ~ sort(c(.x$row_Peak1, .x$row_Peak2)))
    
    # Using the GroupIDs as list headers
    names(list_of_vectors) <- unique(couple_checked$GroupID)
    
    # Reference vector for later adjustment and deletion of unsuitable signals
    check_ref_vector <- unique(couple_checked$GroupID)
    check_delete_vector <- rep(FALSE, length(check_ref_vector))
    
    
    # Checking all identified signals
    for (i in seq_along(list_of_vectors)) {
      # Getting the peak positions of the i. signal
      postion_vector <- list_of_vectors[[i]]
      
      # Reading out the corresponding intensities
      height_absolute <- JRESPeaklist$Intensity[postion_vector]
      
      # Standardize the height information to the highest value
      max_value <- max(height_absolute)
      height_complete_intern <- height_absolute/max_value
      
      # Testing the boundaries
      height_intern <- all(height_complete_intern >= lower_bound &
                             height_complete_intern <= upper_bound)
      
      # If all peaks have the correct height ratio, the signal is marked
      if (height_intern == FALSE) {
        check_delete_vector[i] <- TRUE
      }
      
      # Delete the unsuitable signals in the check_vector
      check_ref_vector <- check_ref_vector[!check_delete_vector]
      
      # Reduce the input data frame
      height_checked <- couple_checked %>%
        filter(GroupID %in% check_ref_vector)
      
    }
    
  } else {
    height_checked <- couple_checked
  }
  return(height_checked)
}


# result section ------------------------------------------------------------
# Check the number of signals in the selected range and select by chem. shift if necessary
# Output of relevant information and error messages

results_even_signal <- function(height_checked, multiplicity, F2shifttarget, coupleconst) {
  # Predefining error messages
  Error_JRES0 <- FALSE
  Error_JRES1 <- FALSE
  
  # Definition of signal type-specific parameters
  if (multiplicity == 2) {
    name_signal <- "doublet"
    nrow_target <- 1
  }
  if
  (multiplicity == 4) {
    name_signal <- "quadruplet"
    nrow_target <- 2
  }
  if
  (multiplicity == 6) {
    name_signal <- "sextet"
    nrow_target <- 3
  }
  if
  (multiplicity == 8) {
    name_signal <- "octet"
    nrow_target <- 4
  }
  if
  (multiplicity > 8) {
    name_signal <- "multiplet"
    nrow_target <- multiplicity/2
  }
  
  # Creating the data frame for marking the signal in the plot function
  mark_signal_data <- data.frame(
    "F1ppm" = numeric(),
    "F2ppm" = numeric(),
    "Intensity" = as.numeric(),
    stringsAsFactors = FALSE
  )
  
  if (nrow(height_checked) < nrow_target) {
    detected_coupling_ppm <- coupleconst
    peak_center <- F2shifttarget
    Error_JRES0 <- TRUE
    
    cat(paste("JRES: No suitable", name_signal, "was identified in the evaluated region!\n", sep = " "))
    
  }
  if (nrow(height_checked) == nrow_target) {
    # Create a list with the position vectors of all peaks per signal
    list_of_vectors <- height_checked %>%
      group_by(GroupID) %>%
      group_split() %>%
      map(~ unique(sort(c(
        .x$row_Peak1, .x$row_Peak2
      ))))
    
    # Creation of coupling vector
    detected_coupling_ppm <- c()
    
    for (n in 1:nrow(height_checked)) {
      coupling_intern <- height_checked$Differenceppm[n] / 2
      detected_coupling_ppm <- c(detected_coupling_ppm, coupling_intern)
    }
    
    detected_coupling_ppm <- c(detected_coupling_ppm, rev(detected_coupling_ppm))
    peak_center <- mean(JRESPeaklist$F2ppm[list_of_vectors[[1]]])
    
    # Appending the peak data to the mark_signal_data dataframe for marking the selected signal in Plot2D
    mark_signal_data <- rbind(
      mark_signal_data,
      data.frame(
        "F1ppm" = as.numeric(c(JRESPeaklist$F1ppm[list_of_vectors[[1]]])),
        "F2ppm" = as.numeric(c(JRESPeaklist$F2ppm[list_of_vectors[[1]]])),
        "Intensity" = as.numeric(c(JRESPeaklist$Intensity[list_of_vectors[[1]]]))
      )
    )
    
    # Output of the results in the console
    cat(
      paste(
        "JRES: Exactly one suitable",
        name_signal,
        "was identified in the evaluated region!\n",
        sep = " "
      )
    )
    cat(paste("chem. shift:\n"))
    
    cat(paste(round(mean(
      mark_signal_data$F2ppm
    ), digits = 4), "[ppm]", "\n"))
    
    cat("Coupling to the signal center:\n")
    
    cat(paste(round(detected_coupling_ppm, digits = 5), "[ppm]"))
    
    cat("\n")
    
    }
  
  if (nrow(height_checked) > nrow_target) {
    # Create a list with the position vectors of all peaks per signal (defined as a group)
    list_of_vectors <- height_checked %>%
      group_by(GroupID) %>%
      group_split() %>%
      map(~ unique(sort(c(
        .x$row_Peak1, .x$row_Peak2
      ))))
    
    # Checking the chem. shift and determination of the next signal
    vector_shiftdistance <- c()
    
    for (i in seq_along(list_of_vectors)) {
      shift_distance <- abs(mean(JRESPeaklist$F2ppm[list_of_vectors[[i]]]) - F2shifttarget)
      vector_shiftdistance <- c(vector_shiftdistance, shift_distance)
    }
    
    nearest_signal <- which.min(vector_shiftdistance)
    
    # Read out and store the difference information as a list by group
    list_of_vectors_diff <- height_checked %>%
      group_by(GroupID) %>%
      group_split() %>%
      map( ~ c(.x$Differenceppm))
    
    detected_coupling_ppm <- c(list_of_vectors_diff[[nearest_signal]]) / 2
    detected_coupling_ppm <- c(detected_coupling_ppm, rev(detected_coupling_ppm))
    peak_center <- mean(JRESPeaklist$F2ppm[list_of_vectors[[nearest_signal]]])
    
    # Appending the peak data to the mark_signal_data dataframe for marking the selected signal in Plot2D
    mark_signal_data <- rbind(
      mark_signal_data,
      data.frame(
        "F1ppm" = as.numeric(c(JRESPeaklist$F1ppm[list_of_vectors[[nearest_signal]]])),
        "F2ppm" = as.numeric(c(JRESPeaklist$F2ppm[list_of_vectors[[nearest_signal]]])),
        "Intensity" = as.numeric(c(JRESPeaklist$Intensity[list_of_vectors[[nearest_signal]]]))
      )
    )
    
    # Output of the results in the console
    cat(
      "JRES: More then one suitable",
      name_signal,
      "was identified in the evaluated region!",
      "\nSignal with the least shift to the F2 target was selected.\n"
    )
    cat(paste("chem. shift:\n"))
    
    cat(paste(round(mean(
      mark_signal_data$F2ppm
    ), digits = 4), "[ppm]", "\n"))
    
    cat("Coupling to the signal center:\n")
    
    cat(paste(round(detected_coupling_ppm, digits = 5), "[ppm]"))
    
    cat("\n")
    
    Error_JRES1 <- TRUE
    
  }
  
  # Assignment of the data to be output
  result_list <- list(
    mark_signal_data = mark_signal_data,
    detected_coupling_ppm = detected_coupling_ppm,
    peak_center = peak_center,
    Error_JRES0 = Error_JRES0,
    Error_JRES1 = Error_JRES1
    )

  return(result_list)
  
}



