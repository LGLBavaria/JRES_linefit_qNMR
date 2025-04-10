# subscripts for selecting odd numbered signals

# Load necessary R-Packages
library(dplyr)
library(purrr)
library(tidyr)

# Identification of the center peaks and removing singulets ---------------------------------------

cut_even_couplings <- function(JRESPeaklist,
                             F1shifttol = 0.0005,
                             F2shifttol = 0.001) {
  
  # Search for signals on/near the zero line (F1 = 0)
  middle_peak <- which((
    JRESPeaklist$F1ppm > (-F1shifttol) &
      JRESPeaklist$F1ppm < F1shifttol
  ),
  arr.ind = T)
  # Creating the deletion vector
  delete_vector <- c(rep(FALSE, length(middle_peak)))
  
  # Searching for and removing singulets
  for (n in 1:length(middle_peak)) {
    # Extract the corresponding values from JRESPeakList
    F2value <- JRESPeaklist[middle_peak[n], "F2ppm"]
    
    couplesignals <- which((
      JRESPeaklist$F2ppm > (F2value - F2shifttol) &
        JRESPeaklist$F2ppm < (F2value + F2shifttol)
    ), arr.ind = T)
    
    # Check whether there is only one signal at this displacement; if yes, mark for removal of the singulet
    if (length(couplesignals) == 1) {
      delete_vector[n] <- TRUE
    }
  }
  
  # Delete the selected items
  middle_peak <- middle_peak[!delete_vector]
  
  return(middle_peak)
}


# Structuring and preselection of possible couplings ---------------------------


filter_odd_couplings <- function(middle_peak, F2shifttol = 0.001) {
  # Create an empty data frame for the difference information
  diff_odd_coupling <- data.frame(
    "row_outer_Peak" = numeric(),
    "row_middle_peak" = numeric(),
    Differenceppm = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (length(middle_peak) > 0) {
    row_couplesignals <- c()
    
    for (n in 1:length(middle_peak)) {
      # Extract the corresponding values from JRESPeakList
      F2value <- JRESPeaklist[middle_peak[n], "F2ppm"]
      
      couplesignals <- c(which((
        JRESPeaklist$F2ppm > (F2value - F2shifttol) &
          JRESPeaklist$F2ppm < (F2value + F2shifttol)
      ), arr.ind = T))
      
      row_couplesignals <- c(row_couplesignals, couplesignals)
    }
    
    
    # Removing the center positions of row_couplesignals
    row_couplesignals <- row_couplesignals[!row_couplesignals %in% middle_peak]
    
    odd_couplings <- JRESPeaklist[row_couplesignals, ]
    
    # Breakdown of positive and negative values
    position_positive <- odd_couplings$F1ppm > 0
    position_negative <- odd_couplings$F1ppm < 0
    
    # Use of the original rows
    rownames_positive <- rownames(odd_couplings)[position_positive]
    rownames_negative <- rownames(odd_couplings)[position_negative]
    
    # Positive and negative coupling indices
    pos_indices <- as.numeric(rownames_positive)
    neg_indices <- as.numeric(rownames_negative)
    
    # reate the possible combinations of peaks
    pos_combinations <- expand.grid(row_Peak1 = pos_indices, middle_peak = middle_peak)
    neg_combinations <- expand.grid(row_Peak1 = neg_indices, middle_peak = middle_peak)
    
    # Calculate the differences and add to the empty dataframe (positive part)
    for (n in 1:nrow(pos_combinations)) {
      rowPeak1 <- as.numeric(pos_combinations[n, 1])
      rowPeak2 <- as.numeric(pos_combinations[n, 2])
      difference <- abs(as.numeric(JRESPeaklist$F1ppm[rowPeak1]) - as.numeric(JRESPeaklist$F1ppm[rowPeak2]))
      
      diff_odd_coupling <- rbind(
        diff_odd_coupling,
        data.frame(
          "row_outer_Peak" = rowPeak1,
          "row_middle_peak" = rowPeak2,
          Differenceppm = difference
        )
      )
    }
    # Calculate the differences and add to the empty dataframe (negative part)
    for (n in 1:nrow(neg_combinations)) {
      rowPeak1 <- as.numeric(neg_combinations[n, 1])
      rowPeak2 <- as.numeric(neg_combinations[n, 2])
      difference <- abs(as.numeric(JRESPeaklist$F1ppm[rowPeak1]) - as.numeric(JRESPeaklist$F1ppm[rowPeak2]))
      
      diff_odd_coupling <- rbind(
        diff_odd_coupling,
        data.frame(
          "row_outer_Peak" = rowPeak1,
          "row_middle_peak" = rowPeak2,
          Differenceppm = difference
        )
      )
    }
    
    # Exclusion of lines with NAs; can occur at low intensities if peaks are completely missing on one side
    diff_odd_coupling <- diff_odd_coupling[complete.cases(diff_odd_coupling), ]
    
    # Exclusion of couplings that are not in the same position
    diff_delete_rows <- rep(FALSE, nrow(diff_odd_coupling))
    
    # Iteration over the rows of diff_odd_coupling
    for (i in 1:nrow(diff_odd_coupling)) {
      row1 <- diff_odd_coupling[i, 1]
      row2 <- diff_odd_coupling[i, 2]
      
      # Extract the corresponding values from JRESPeakList
      value1 <- JRESPeaklist[row1, "F2ppm"]
      value2 <- JRESPeaklist[row2, "F2ppm"]
      
      # Comparison and marking
      if (abs(value1 - value2) > F2shifttol) {
        diff_delete_rows[i] <- TRUE
      }
    }
    
    # Remove the couplings that are not in the same position
    diff_odd_coupling <- diff_odd_coupling[!diff_delete_rows, ]
    
    # Remove dublicate rows in diff_odd_coupling
    diff_odd_coupling$row_outer_Peak <- round(diff_odd_coupling$row_outer_Peak)
    
    diff_odd_coupling <- diff_odd_coupling %>% distinct(row_outer_Peak,
                                                row_middle_peak,
                                                Differenceppm,
                                                .keep_all = TRUE)
    
  
  }
  return(diff_odd_coupling)
}


# Selection of the desired even multiplicity ---------------------

find_odd_coupled_signals <- function(diff_odd_coupling, multiplicity, F2shifttol = 0.001) {
  
  if (nrow(diff_odd_coupling) >= 2){
  # Sort the determined couplings into F2 groups of the same position
  # Add the chemical shift to the data
  diff_odd_coupling$F2ppm_outer_Peak <- JRESPeaklist$F2ppm[diff_odd_coupling$row_outer_Peak]
  diff_odd_coupling$F2ppm_middle_peak <- JRESPeaklist$F2ppm[diff_odd_coupling$row_middle_peak]
  
  # Empty lists for saving groups and IDs
  groups <- list()
  group_ids <- integer(nrow(diff_odd_coupling))  # Initialize a column for group IDs
  group_id <- 1  # Start of the group ID
  
  # Loop through each row of the data set
  for (i in 1:nrow(diff_odd_coupling)) {
    current_peak2 <- diff_odd_coupling$F2ppm_middle_peak[i]
    
    # Check whether outer_peak or middle_peak are already in an existing group
    in_group <- FALSE
    for (j in seq_along(groups)) {
      if (any(abs(groups[[j]] - current_peak2) <= F2shifttol)) {
        # Add the peaks to the existing group
        groups[[j]] <- c(groups[[j]], current_peak2)
        group_ids[i] <- j  # Assigns the group ID to the line
        in_group <- TRUE
        break
      }
    }
    # If the peak does not yet belong to a group, create a new group
    if (!in_group) {
      groups[[group_id]] <- c(current_peak2)
      group_ids[i] <- group_id  # Assign new group ID
      group_id <- group_id + 1
    }
  }
  
  
  # for (i in 1:nrow(diff_odd_coupling)) {
  #   current_peak1 <- diff_odd_coupling$F2ppm_outer_Peak[i]
  #   current_peak2 <- diff_odd_coupling$F2ppm_middle_peak[i]
  #   
  #   # Check whether outer_peak or middle_peak are already in an existing group
  #   in_group <- FALSE
  #   for (j in seq_along(groups)) {
  #     if (any(abs(groups[[j]] - current_peak1) <= F2shifttol) ||
  #         any(abs(groups[[j]] - current_peak2) <= F2shifttol)) {
  #       
  #       # Add the peaks to the existing group
  #       groups[[j]] <- c(groups[[j]], current_peak1, current_peak2)
  #       group_ids[i] <- j  # Assigns the group ID to the line
  #       in_group <- TRUE
  #       break
  #     }
  #   }
  #   
  #   # If the peak does not yet belong to a group, create a new group
  #   if (!in_group) {
  #     groups[[group_id]] <- c(current_peak1, current_peak2)
  #     group_ids[i] <- group_id  # Assign new group ID
  #     group_id <- group_id + 1
  #   }
  # }
  
  # Assign the group IDs to the source data set
  diff_odd_coupling$GroupID <- group_ids
  
  # Calculate the number of each group
  group_counts <- table(diff_odd_coupling$GroupID)
  
  # Filter the groups that occur exactly X times (number of possible couplings in odd numbered multiplicities: x = (multiplicity-1)
  groups_with_x <- as.numeric(names(group_counts[group_counts == (multiplicity - 1)]))
  
  # Create a new dataframe that only contains the rows whose GroupID occurs X times
  multiplicity_checked <- diff_odd_coupling[diff_odd_coupling$GroupID %in% groups_with_x, ]
  }
  if (nrow(diff_odd_coupling) < 2) {
    multiplicity_checked <- diff_odd_coupling %>% mutate(GroupID = NA)
  }
  return(multiplicity_checked)
}

# Checking the coupling information --------------------------------------

check_odd_coupling <- function(multiplicity_checked,
                                coupleconst,
                                multiplicity,
                                coupletol = 0.01) {
  diff_checked <- c()
  for (m in 1:(length(coupleconst) / 2)) {
    # Check whether external coupling is fulfilled
    diff_checked_intern <- which((
      multiplicity_checked$Differenceppm > coupleconst[m] - coupletol * coupleconst[m] &
        multiplicity_checked$Differenceppm < coupleconst[m] + coupletol * coupleconst[m]
    ),
    arr.ind = T
    )
    
    # Add to result vector
    diff_checked <- c(diff_checked, diff_checked_intern)
  }
  
  # Output of the coupling peaks
  couple_checked <- multiplicity_checked[diff_checked, ]
  
  # Calculate the number of each group
  group_counts <- table(couple_checked$GroupID)
  
  # Filter the groups that occur exactly X times (number of possible couplings in odd numbered multiplicities: x = (multiplicity-1)
  groups_with_x <- as.numeric(names(group_counts[group_counts == (multiplicity - 1)]))
  
  # Create a new dataframe that only contains the rows whose GroupID occurs X times
  couple_checked <- couple_checked[couple_checked$GroupID %in% groups_with_x, ]
  
  return(couple_checked)
}

# Checking height ratios -----------------------------------------------
# Caveat: Signal multiplicity is derived from vector heightratio
# Height ratio check implemented as optional; to activate check_height = TRUE
check_height_odd_signal <- function(couple_checked,
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
      map( ~ unique(sort(c(.x$row_outer_Peak, .x$row_middle_peak))))
    
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

results_odd_signal <- function(height_checked, multiplicity, F2shifttarget, coupleconst) {
  # Predefining error messages
  Error_JRES0 <- FALSE
  Error_JRES1 <- FALSE
  
  # Definition of signal type-specific parameters
  if (multiplicity == 3) {
    name_signal <- "triplet"
    nrow_target <- 2
  }
  if
  (multiplicity == 5) {
    name_signal <- "quintet"
    nrow_target <- 4
  }
  if
  (multiplicity == 7) {
    name_signal <- "septet"
    nrow_target <- 6
  }
  if
  (multiplicity == 9) {
    name_signal <- "nonet"
    nrow_target <- 8
  }
  if
  (multiplicity > 9) {
    name_signal <- "multiplet"
    nrow_target <- multiplicity - 1
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
    
  } else {
    if (nrow(height_checked) == nrow_target) {
      
      # Create a list with the position vectors of all peaks per signal
      list_of_vectors <- height_checked %>%
        group_by(GroupID) %>%
        group_split() %>%
        map( ~ unique(sort(c(.x$row_outer_Peak, .x$row_middle_peak))))
      
      # Creation of coupling vector
      left_detected_coupling_ppm <- c()
      
      for (n in seq(from = 1, to = nrow(height_checked), by = 2)) {
        coupling_intern <- height_checked$Differenceppm[n]
        left_detected_coupling_ppm <- c(left_detected_coupling_ppm, coupling_intern)
      }
      
      right_detected_coupling_ppm <- c()
      for (n in seq(from = 2, to = nrow(height_checked), by = 2)) {
        coupling_intern <- height_checked$Differenceppm[n]
        right_detected_coupling_ppm <- c(right_detected_coupling_ppm, coupling_intern)
      }
      
      detected_coupling_ppm <- c(left_detected_coupling_ppm, 0, rev(right_detected_coupling_ppm))
      peak_center <- mean(JRESPeaklist$F2ppm[list_of_vectors[[1]]])
      
      # Appending the peak data to the mark_signal_data dataframe for marking the selected signal in Plot2D
      mark_signal_data <- rbind(mark_signal_data, data.frame(
        "F1ppm" = as.numeric(c(JRESPeaklist$F1ppm[list_of_vectors[[1]]])),
        "F2ppm" = as.numeric(c(JRESPeaklist$F2ppm[list_of_vectors[[1]]])),
        "Intensity" = as.numeric(c(JRESPeaklist$Intensity[list_of_vectors[[1]]]))
      ))
      
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
        map( ~ unique(sort(c(.x$row_outer_Peak, .x$row_middle_peak))))
      
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
        map(~ c(.x$Differenceppm))
      
      list_coupling <- c(list_of_vectors_diff[[nearest_signal]])
      
      # Creation of coupling vector
      left_detected_coupling_ppm <- c()
      
      for (n in seq(from = 1, to = length(list_coupling), by = 2)) {
        coupling_intern <- list_coupling[n]
        left_detected_coupling_ppm <- c(left_detected_coupling_ppm, coupling_intern)
      }
      
      right_detected_coupling_ppm <- c()
      
      for (n in seq(from = 2, to = length(list_coupling), by = 2)) {
        coupling_intern <- list_coupling[n]
        right_detected_coupling_ppm <- c(right_detected_coupling_ppm, coupling_intern)
      }
      
      detected_coupling_ppm <- c(left_detected_coupling_ppm, 0, rev(right_detected_coupling_ppm))
      peak_center <- mean(JRESPeaklist$F2ppm[list_of_vectors[[nearest_signal]]])
      
      # Appending the peak data to the mark_signal_data dataframe for marking the selected signal in Plot2D
      mark_signal_data <- rbind(mark_signal_data, data.frame(
        "F1ppm" = as.numeric(c(JRESPeaklist$F1ppm[list_of_vectors[[nearest_signal]]])),
        "F2ppm" = as.numeric(c(JRESPeaklist$F2ppm[list_of_vectors[[nearest_signal]]])),
        "Intensity" = as.numeric(c(JRESPeaklist$Intensity[list_of_vectors[[nearest_signal]]]))
      ))
      
      # Output of the results in the console
      cat(
        "JRES: More then one suitable", name_signal, "was identified in the evaluated region!",
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
  }
  
  # Assignment of the data to be output
  result_list <- list(
    mark_signal_data = mark_signal_data,
    peak_center = peak_center,
    detected_coupling_ppm = detected_coupling_ppm,
    Error_JRES0 = Error_JRES0,
    Error_JRES1 = Error_JRES1
  )
  
  return(result_list)
  
}

