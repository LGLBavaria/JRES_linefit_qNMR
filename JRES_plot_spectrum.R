# Load necessary R-Packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(mrbin)
library(purrr)
library(tidyr)


# 1. Plot a spectral section with optional marking of picked Signal  ------

# necessary input:
# sub_spectrum: dataframe with subregion of J-resolved spectrum selected generated in JRES_Ident_Peaks

# optional input:
# mark_signal_data: dataframe generated in the signal identification procedure, part of the output result_list 
# sampleID: ID-number of the examined sample, used for the plot heading (Default: NA)
# analyte: name of the examined compound, used for the plot heading (Default: NA)
# mark_signal: turns the marking of signals on and off (TRUE = On, False = Off; Default: FALSE)
# mark_color: color of the signal markings and information (Default: red)
# mark_rel_x_shift: defines the relative offset between signal marking and information box (Default: 0.05 = 5%)


JRES_plot_spectrum <- function(sub_spectrum,
                              mark_signal_data,
                              sampleID = "NA",
                              analyte = "NA",
                              mark_signal = FALSE,
                              mark_color = "red",
                              mark_rel_x_shift = 0.05) {
  
  # Converting the spectrum to a long format (necessary for ggplot2)
  sub_spectrum_long <- sub_spectrum %>%
    pivot_longer(cols = everything(),
                 names_to = "F2ppm",
                 values_to = "Intensity") %>%
    mutate(F1ppm = rep(rownames(sub_spectrum), each = ncol(sub_spectrum)))
  
  sub_spectrum_long$F1ppm <- as.numeric(as.character(sub_spectrum_long$F1ppm))
  sub_spectrum_long$F2ppm <- as.numeric(as.character(sub_spectrum_long$F2ppm))
  
  # Creation of the base plot of the spectral region
  p <- ggplot(sub_spectrum_long, aes(x = F2ppm, y = F1ppm, z = Intensity)) +
    geom_contour_filled(aes(fill = after_stat(level)), show.legend = FALSE) +  # Filled contours with automatic breaks
    geom_contour(color = "black", linewidth = 0.2) +
    scale_fill_grey(start = 1, end = 0) +  # Grayscale for discrete contours (start = 1 for black, end = 0 for white)
    labs(
      title = paste(sampleID, analyte, "J-resolved", sep = "_"),
      x = "F2 [ppm]",
      y = "F1 [ppm]"
    ) +
    theme_minimal() +
    scale_x_reverse(
      expand = c(0, 0),
      breaks = seq(round(min(
        sub_spectrum_long$F2ppm
      ), digits = 2), round(max(
        sub_spectrum_long$F2ppm
      ), digits = 2), by = 0.025),
      minor_breaks = seq(round(min(
        sub_spectrum_long$F2ppm
      ), digits = 2), round(max(
        sub_spectrum_long$F2ppm
      ), digits = 2), by = 0.005),
      guide = guide_axis(minor.ticks = TRUE)
    ) +  # Prevents additional spacing on the x-axis
    scale_y_reverse(
      expand = c(0, 0),
      breaks = seq(round(min(
        sub_spectrum_long$F1ppm
      ), digits = 2), round(max(
        sub_spectrum_long$F1ppm
      ), digits = 2), by = 0.025),
      minor_breaks = seq(round(min(
        sub_spectrum_long$F1ppm
      ), digits = 2), round(max(
        sub_spectrum_long$F1ppm
      ), digits = 2), by = 0.005),
      guide = guide_axis(minor.ticks = TRUE)
    ) + # Prevents additional spacing on the y-axis
    theme(
      panel.grid = element_blank(),
      # Removes the grid
      legend.background = element_rect(fill = "white"),
      # Set background to white
      axis.line = element_line(color = "black"),
      # Intensifies the axis lines
      axis.ticks = element_line(color = "black"),
      # Intensifies the axis ticks
      axis.text = element_text(size = 10),
      # Larger font size for axis labels
      axis.title = element_text(size = 12) # Larger and bolder font for axis titles
    )
  
  
  # Optional: Insert markers for the selected signals ----------
  
  # Creating the data point label shift
  x_shift <- diff(range(sub_spectrum_long$F2ppm)) * mark_rel_x_shift
  
  # If selected peaks are to be marked, function input mark_signal is queried. If TRUE, marking is carried out.
  if (isTRUE(mark_signal)) {
    # Markierungspunkte annehmen, dass sie bereits korrekt in einer DataFrame-Struktur vorliegen
    p <- p +
      geom_point(
        data = mark_signal_data,
        aes(x = F2ppm, y = F1ppm),
        color = mark_color,
        size = 3,
        shape = 4
      ) +
      # geom_text(data = mark_signal_data, aes(label = paste("F2[ppm] =", round(F2ppm, 3))), position = position_dodge(width = 0.1)) +
      # geom_text(data = mark_signal_data, aes(label = paste("F1[ppm] =", round(F1ppm, 3))), position = position_dodge(width = 0.1), vjust = -0.05)
    
      geom_label_repel(
        data = mark_signal_data,
        aes(
          x = F2ppm,
          y = F1ppm,
          label = paste(
            "F2[ppm] =", round(F2ppm, 4), "\nF1[ppm] =", round(F1ppm, 4)
          )
        ),
        color = mark_color,  # text color
        fill = scales::alpha("white", 0.5),     # background color of the frame
        direction = "y",     # Move text only vertically
        nudge_x = x_shift,   # Horizontal offset to the right
        hjust = 0,           # Align to the right of the points
        box.padding = 0.3,   # Distance between labels and points
        point.padding = 0.2, # spacing between dots and text
        max.overlaps = Inf,  # Show all texts
        fontface = "bold"    # Increase font thickness
      )
  }
  
  # return plot
  print(p)
  
  return(p)
}
