# Load necessary R-Packages
library(dplyr)
library(tidyr)
library(stats) # For dist function

#' Calculate Centers of Mass for Coherent Peaks in a Spectral Region
#'
#' This function identifies coherent signal regions (peaks) within a specified
#' rectangular region of a 2D NMR spectrum using a Connected-Component-Algorithm (CCA).
#' It calculates the center of mass (CoM) for each identified coherent region
#' and compiles a list of all these CoMs and their associated base area points.
#' This function does not perform target selection.
#'
#' @param spectrum A 2D numeric matrix representing the JRES spectrum.
#'   Row names should be F1 (coupling) ppm values and column names F2 (chemical shift) ppm values.
#' @param ppm_f2_min Numeric, minimum F2 ppm value for the region of interest.
#' @param ppm_f2_max Numeric, maximum F2 ppm value for the region of interest.
#' @param ppm_f1_min Numeric, minimum F1 ppm value for the region of interest.
#'   (Default: -0.1, often the minimal possible F1 value in JRES).
#' @param ppm_f1_max Numeric, maximum F1 ppm value for the region of interest.
#'   (Default: 0.1, often the maximal possible F1 value in JRES).
#' @param intensity_threshold Numeric, only intensities above this threshold will be considered for CCA and CoM calculation.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{sub_spectrum}: A dataframe representation of the extracted spectral sub-region (matrix).
#'   \item \code{JRES_CoMList}: A dataframe where each row represents a distinct coherent component's
#'         Center of Mass, including F1ppm, F2ppm, total 'Intensity' (sum of intensities in the component),
#'         and a unique 'ComponentID'. Returns an empty dataframe if no components are found.
#'   \item \code{all_base_areas_df}: A dataframe of all individual points (F1ppm, F2ppm, Intensity, ComponentID)
#'         that belong to any coherent component identified.
#'         This can be used for plotting the base areas of all found components.
#' }
#' @importFrom dplyr filter mutate rename select bind_rows
#' @importFrom tidyr pivot_longer
#' @export

JRES_calculate_CoM <- function(spectrum,
                               ppm_f2_min,
                               ppm_f2_max,
                               ppm_f1_min = -0.1, # Corrected default as per documentation
                               ppm_f1_max = 0.1,  # Corrected default as per documentation
                               intensity_threshold = 1000) {
  
  if (!is.matrix(spectrum)) {
    spectrum <- as.matrix(spectrum)
  }
  
  couple_vector <- as.numeric(rownames(spectrum))
  shift_vector <- as.numeric(colnames(spectrum))
  
  f1_indices <- which(couple_vector >= ppm_f1_min & couple_vector <= ppm_f1_max)
  f2_indices <- which(shift_vector >= ppm_f2_min & shift_vector <= ppm_f2_max)
  
  # Initialize empty dataframes for results
  JRES_CoMList <- data.frame(F1ppm = numeric(), F2ppm = numeric(), Intensity = numeric(), ComponentID = integer())
  all_base_areas_df <- data.frame(F1ppm = numeric(), F2ppm = numeric(), Intensity = numeric(), ComponentID = integer())
  
  if (length(f1_indices) == 0 || length(f2_indices) == 0) {
    # warning("CoM calculation: No data points found within the specified PPM range. Returning empty results.")
    return(list(sub_spectrum = data.frame(), JRES_CoMList = JRES_CoMList, all_base_areas_df = all_base_areas_df))
  }
  
  peak_region_data <- spectrum[f1_indices, f2_indices, drop = FALSE]
  
  # Get actual ppm values for the extracted sub-region
  sub_couple_vector <- as.numeric(rownames(peak_region_data))
  sub_shift_vector <- as.numeric(colnames(peak_region_data))
  
  # Convert the peak region to a long format and add grid indices (1-based within peak_region_data)
  all_points <- as.data.frame(as.table(peak_region_data)) %>%
    dplyr::rename(F1_label = Var1, F2_label = Var2, Intensity = Freq) %>%
    dplyr::mutate(
      F1ppm = as.numeric(as.character(F1_label)),
      F2ppm = as.numeric(as.character(F2_label)),
      row_idx = match(F1ppm, sub_couple_vector),
      col_idx = match(F2ppm, sub_shift_vector)
    ) %>%
    dplyr::select(F1ppm, F2ppm, Intensity, row_idx, col_idx)
  
  points_for_cca <- all_points %>%
    dplyr::filter(Intensity > intensity_threshold)
  
  if (nrow(points_for_cca) == 0) {
    # warning("CoM calculation: No positive intensities above threshold found in the specified region. Returning empty results.")
    return(list(sub_spectrum = as.data.frame(peak_region_data), JRES_CoMList = JRES_CoMList, all_base_areas_df = all_base_areas_df))
  }
  
  # Prepare for Connected Components Algorithm (CCA)
  # point_map: key is "row_idx_col_idx", value is intensity
  points_for_cca$key <- paste(points_for_cca$row_idx, points_for_cca$col_idx, sep = "_")
  point_map <- points_for_cca$Intensity
  names(point_map) <- points_for_cca$key
  
  # visited_keys: To track visited points during BFS
  visited_keys <- new.env(hash = TRUE, parent = emptyenv())
  
  all_components_coms_list <- list()
  component_id_counter <- 0
  all_coherent_points_list <- list()
  
  # Iterate through all potential points to find components
  for (k in 1:nrow(points_for_cca)) {
    current_row <- points_for_cca$row_idx[k]
    current_col <- points_for_cca$col_idx[k]
    current_key <- points_for_cca$key[k]
    
    if (is.null(visited_keys[[current_key]])) { # If point not yet visited
      component_id_counter <- component_id_counter + 1
      current_component_points_df <- data.frame(F1ppm = numeric(), F2ppm = numeric(), Intensity = numeric())
      
      # BFS queue stores keys (row_idx_col_idx)
      q <- list(current_key)
      visited_keys[[current_key]] <- TRUE
      
      head_idx <- 1
      while (head_idx <= length(q)) {
        dequeue_key <- q[[head_idx]]
        head_idx <- head_idx + 1
        
        parts <- as.numeric(strsplit(dequeue_key, "_")[[1]])
        drow <- parts[1]
        dcol <- parts[2]
        
        # Check if the key exists in point_map before accessing.
        # (Defensive programming; dequeue_key should always be valid as it originates from point_map.)
        if (! (dequeue_key %in% names(point_map))) {
          # warning(paste("Dequeued key '", dequeue_key, "' not found in point_map. Skipping.", sep = ""))
          next
        }
        dintensity <- point_map[[dequeue_key]]
        
        # Add point to current component (with original PPMs)
        current_component_points_df <- dplyr::bind_rows(current_component_points_df,
                                                        data.frame(F1ppm = sub_couple_vector[drow],
                                                                   F2ppm = sub_shift_vector[dcol],
                                                                   Intensity = dintensity))
        
        # Check 8 neighbors (offsets -1, 0, 1 for row_idx and col_idx)
        for (row_offset in c(-1, 0, 1)) {
          for (col_offset in c(-1, 0, 1)) {
            if (row_offset == 0 && col_offset == 0) next # Skip self
            
            neighbor_row <- drow + row_offset
            neighbor_col <- dcol + col_offset
            
            # Check if neighbor is within peak_region_data bounds
            if (neighbor_row < 1 || neighbor_row > nrow(peak_region_data) ||
                neighbor_col < 1 || neighbor_col > ncol(peak_region_data)) {
              next # Skip neighbor if out of bounds
            }
            
            neighbor_key <- paste(neighbor_row, neighbor_col, sep = "_")
            
            # First, check if neighbor key exists in point_map names
            if (neighbor_key %in% names(point_map)) {
              if (is.null(visited_keys[[neighbor_key]])) { # Then, check if visited
                visited_keys[[neighbor_key]] <- TRUE
                q[[length(q) + 1]] <- neighbor_key
              }
            }
          }
        }
      } # End of BFS for current component
      
      # Calculate CoM for this coherent component
      if (nrow(current_component_points_df) > 0 && sum(current_component_points_df$Intensity) > 0) {
        comp_f1_center <- sum(current_component_points_df$F1ppm * current_component_points_df$Intensity) / sum(current_component_points_df$Intensity)
        comp_f2_center <- sum(current_component_points_df$F2ppm * current_component_points_df$Intensity) / sum(current_component_points_df$Intensity)
        comp_total_intensity <- sum(current_component_points_df$Intensity) # CoM intensity is sum of component intensities
        
        all_components_coms_list[[length(all_components_coms_list) + 1]] <- data.frame(
          F1ppm = comp_f1_center,
          F2ppm = comp_f2_center,
          Intensity = comp_total_intensity,
          ComponentID = component_id_counter
        )
        
        # Store points for plotting all base areas
        current_component_points_df$ComponentID <- component_id_counter
        all_coherent_points_list[[length(all_coherent_points_list) + 1]] <- current_component_points_df
      }
    }
  }
  
  if (length(all_components_coms_list) > 0) {
    JRES_CoMList <- dplyr::bind_rows(all_components_coms_list)
  } else {
    # warning("CoM calculation: No coherent peaks found above threshold. Returning empty results.")
    return(list(sub_spectrum = as.data.frame(peak_region_data), JRES_CoMList = JRES_CoMList, all_base_areas_df = all_base_areas_df))
  }
  
  # Combine all coherent points for base_areas_df plotting
  if (length(all_coherent_points_list) > 0) {
    all_base_areas_df <- dplyr::bind_rows(all_coherent_points_list)
  }
  
  return(list(
    sub_spectrum = as.data.frame(peak_region_data), # Return the extracted sub-spectrum as a dataframe
    JRES_CoMList = JRES_CoMList,
    all_base_areas_df = all_base_areas_df
  ))
}


#' Select the Closest Center of Mass to a Target Position
#'
#' This function takes a dataframe of multiple Centers of Mass (CoMs) and
#' a dataframe of all associated base area points, then selects the CoM
#' that is closest to a specified target F1 and F2 PPM position. It also
#' returns the base area points specifically for the selected cluster.
#'
#' @param JRES_CoMList A dataframe containing the F1ppm, F2ppm, Intensity, and ComponentID
#'   for all identified coherent components' CoMs (e.g., from `JRES_calculate_CoM`).
#' @param all_base_areas_df A dataframe containing all individual points (F1ppm, F2ppm, Intensity, ComponentID)
#'   that belong to any coherent component (e.g., from `JRES_calculate_CoM`).
#' @param target_f1_val Numeric, the F1 ppm value of the target position for selection.
#' @param target_f2_val Numeric, the F2 ppm value of the target position for selection.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{selected_com_point_df}: A dataframe with one row, containing the F1ppm, F2ppm
#'         coordinates of the selected center of mass and its total 'Intensity'. This is
#'         conceptually similar to 'mark_signal_data' for the chosen cluster. Returns an
#'         empty dataframe if no CoMs are provided.
#'   \item \code{selected_base_area_df}: A dataframe of all points (F1ppm, F2ppm, Intensity, ComponentID)
#'         that belong *only* to the selected coherent component. This can be used for plotting.
#'         Returns an empty dataframe if no CoMs are provided.
#' }
#' @importFrom dplyr filter select
#' @importFrom stats dist
#' @export

JRES_select_target_CoM <- function(JRES_CoMList, all_base_areas_df, target_f1_val, target_f2_val) {
  
  # Initialize selected dataframes as truly empty
  selected_com_point_df <- data.frame(F1ppm = numeric(), F2ppm = numeric(), Intensity = numeric())
  selected_base_area_df <- data.frame(F1ppm = numeric(), F2ppm = numeric(), Intensity = numeric(), ComponentID = integer())

  Error_JRES0 <- FALSE
  Error_JRES1 <- FALSE
  
  if (nrow(JRES_CoMList) == 0) {
    Error_JRES0 <- TRUE
    return(list(selected_com_point_df = selected_com_point_df, selected_base_area_df = selected_base_area_df, Error_JRES0 = Error_JRES0, Error_JRES1 = Error_JRES1))
  }
  
    if (nrow(JRES_CoMList) > 1) {
    Error_JRES1 <- TRUE
  }
  
  # Calculate distances to target position and select the closest CoM
  target_coords <- c(target_f1_val, target_f2_val)
  
  # Calculate Euclidean distance for each component's CoM to the target_coords
  distances <- apply(JRES_CoMList[, c("F1ppm", "F2ppm")], 1, function(row_coords) {
    stats::dist(rbind(row_coords, target_coords))[1]
  })
  
  closest_com_idx <- which.min(distances)
  selected_com_point_df <- JRES_CoMList[closest_com_idx, ] %>%
    dplyr::select(F1ppm, F2ppm, Intensity) # Remove ComponentID for final output
  
  # Filter the base areas to include only points from the selected component
  selected_component_id <- JRES_CoMList[closest_com_idx, "ComponentID"]
  selected_base_area_df <- all_base_areas_df %>%
    dplyr::filter(ComponentID == selected_component_id)
  
  return(list(
    selected_com_point_df = selected_com_point_df,
    selected_base_area_df = selected_base_area_df,
    Error_JRES0 = Error_JRES0,
    Error_JRES1 = Error_JRES1
  ))
}
