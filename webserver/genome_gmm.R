# Keep every requested biological component numerically present in automatic
# fits.  This is intentionally a very small floor: it prevents a harmonic from
# collapsing to effectively zero (as in virginalis) without forcing weak peaks
# to carry the same mass as clear peaks.  The optional n+1 nuisance component
# is not floored.
stabilize_biological_alpha <- function(alpha, required_components,
                                       minimum_share = 0.005) {
  alpha <- as.numeric(alpha)
  total <- sum(alpha)
  required_components <- min(as.integer(required_components), length(alpha))
  if (!is.finite(total) || total <= 0 || required_components < 1L) return(alpha)
  proportions <- pmax(alpha / total, 0)
  proportions[seq_len(required_components)] <- pmax(
    proportions[seq_len(required_components)], minimum_share
  )
  proportions <- proportions / sum(proportions)
  proportions * total
}

# Build the same continuous reference curve that is drawn in grey in the PDF.
# Keeping this in one helper prevents the optimiser from fitting integer-bin
# heights while the user judges a different, spline-interpolated curve.
make_gmm_reference_curve <- function(x, y, dense_points = NULL) {
  valid <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  ordered <- order(x)
  x <- x[ordered]
  y <- y[ordered]
  if (anyDuplicated(x)) {
    grouped <- split(y, x)
    x <- as.numeric(names(grouped))
    y <- vapply(grouped, mean, numeric(1))
    ordered <- order(x)
    x <- x[ordered]
    y <- y[ordered]
  }
  if (is.null(dense_points)) {
    # The PDF uses 1000 points; using the identical resolution also keeps the
    # automatic candidate comparison practical for large test collections.
    dense_points <- 1000L
  }
  dense_points <- max(50L, as.integer(dense_points))
  if (length(x) > 3L) {
    reference <- stats::spline(x, y, n = dense_points)
  } else if (length(x) > 1L) {
    dense_x <- seq(min(x), max(x), length.out = dense_points)
    reference <- list(
      x = dense_x,
      y = stats::approx(x, y, xout = dense_x, rule = 2)$y
    )
  } else {
    reference <- list(x = x, y = y)
  }
  reference$y[!is.finite(reference$y)] <- 0
  reference
}

# Refine an integer-bin local maximum on the continuous grey reference.  The
# correction is deliberately limited to two bins: the five-bin detection
# smoother can shift a skewed peak by more than one bin, while two bins remain
# far below the biological peak spacing and cannot jump to a neighbouring mode.
estimate_continuous_peak <- function(x, y, peak_x, max_shift_bins = 2.0) {
  if (!length(peak_x) || !is.finite(peak_x)) return(peak_x)
  x_steps <- diff(sort(unique(as.numeric(x))))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  # Histograms can have a very long, sparse high-frequency tail.  Constructing
  # one uniform spline over that whole range would undersample the genomic peak.
  neighbourhood <- which(
    x >= peak_x - (max_shift_bins + 3) * bin_width &
      x <= peak_x + (max_shift_bins + 3) * bin_width
  )
  if (length(neighbourhood) < 4L) return(as.numeric(peak_x))
  reference <- make_gmm_reference_curve(
    x[neighbourhood], y[neighbourhood], dense_points = 801L
  )
  local <- which(
    reference$x >= peak_x - max_shift_bins * bin_width &
      reference$x <= peak_x + max_shift_bins * bin_width
  )
  if (!length(local)) return(as.numeric(peak_x))
  reference$x[local[which.max(reference$y[local])]]
}

estimate_reference_height <- function(x, y, centre) {
  if (!length(centre) || !is.finite(centre)) return(NA_real_)
  x_steps <- diff(sort(unique(as.numeric(x))))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  neighbourhood <- which(
    x >= centre - 5 * bin_width & x <= centre + 5 * bin_width
  )
  if (length(neighbourhood) < 4L) {
    return(y[which.min(abs(x - centre))])
  }
  reference <- make_gmm_reference_curve(
    x[neighbourhood], y[neighbourhood], dense_points = 801L
  )
  stats::approx(reference$x, reference$y, xout = centre, rule = 2)$y
}

# Internal Gaussian-mixture fitting engine. Its numerical arguments are
# candidate values supplied only by the automatic search; they are not part of
# the user-facing interface.
.GenomeGMM_fit <- function(data_file, n, main_peak_ratio = NULL, symmetry_tolerance = 0.2, 
                      main_peak_weight = 1, use_calibration = FALSE, species_name = "", 
                      max_sigma_ratio = 1.35,
                      sigma_mode = c("legacy", "adaptive", "common", "trend"),
                      position_mode = c("legacy", "adaptive", "relaxed"),
                      joint_refinement = c("off", "auto", "always"),
                      save_plot = TRUE,
                      em_max_iter = NULL, em_tolerance = 1e-8,
                      extra_component_mode = c("legacy", "always", "never"),
                      .preloaded_data = NULL, .analysis_max_x = NULL, ...) {
  extra_component_mode <- match.arg(extra_component_mode)
  sigma_mode <- match.arg(sigma_mode)
  position_mode <- match.arg(position_mode)
  joint_refinement <- match.arg(joint_refinement)
  
  if (missing(n)) {
    stop("Parameter 'n' (number of Gaussian components) must be specified.")
  }
  
  if (!is.numeric(n) || n <= 0 || n != round(n)) {
    stop("Parameter 'n' must be a positive integer.")
  }
  
  # Check if data file exists
  if (!file.exists(data_file)) {
    stop(paste("Data file '", data_file, "' does not exist, please check the file path."))
  }
  
  # Core logic: If main_peak_ratio is specified by the user, automatically enable calibration
  if (!is.null(main_peak_ratio)) {
    use_calibration <- TRUE
    cat("main_peak_ratio specified (", main_peak_ratio, "), automatically enabling calibration\n")
  } else {
    # If user did not specify main_peak_ratio, use default value
    main_peak_ratio <- 0.3
  }
  
  # Extract file name without path and extension
  file_basename <- tools::file_path_sans_ext(basename(data_file))
  
  # Read full data
  full_data <- if (is.null(.preloaded_data)) {
    read.table(data_file, header = FALSE)
  } else {
    .preloaded_data
  }
  
  # Find start point using gradient change
  find_start_by_gradient <- function(y_values, min_gradient = 0.005) {
    n <- length(y_values)
    gradients <- diff(y_values) / y_values[-n]
    
    for (i in seq_len(n - 1)) {
      if (i + 2 <= n) {
        avg_gradient <- mean(gradients[i:(i + 2)], na.rm = TRUE)
        if (avg_gradient > min_gradient) {
          # The historical gradient trigger can fire one bin before the
          # low-frequency error tail reaches its valley. Keep that trigger,
          # then finish at the lowest count in a short forward neighbourhood.
          valley_end <- min(n, i + 3L)
          valley_window <- i:valley_end
          return(valley_window[which.min(y_values[valley_window])])
        }
      }
    }
    return(1)
  }
  
  # Determine start point
  start_index <- find_start_by_gradient(full_data[, 2])
  cat("Dynamically determined start index:", start_index, "\n")
  cat("Corresponding k-mer frequency:", full_data[start_index, 1], "\n")
  
  # Filter data
  filtered_data <- full_data[start_index:nrow(full_data), ]
  if (!is.null(.analysis_max_x) && is.finite(.analysis_max_x)) {
    filtered_data <- filtered_data[filtered_data[, 1] <= .analysis_max_x, ]
  }
  x_filtered <- filtered_data[, 1]
  y_filtered <- filtered_data[, 2]
  
  
  find_peaks_properly <- function(x_values, y_values, min_peak_prominence = 0.1) {
    # Smooth data for peak detection
    y_smooth <- stats::filter(y_values, rep(1, 5)/5, sides = 2)
    y_smooth[is.na(y_smooth)] <- y_values[is.na(y_smooth)]
    
    # Find local maxima
    is_peak <- diff(sign(diff(y_smooth))) < 0
    peak_indices <- which(is_peak) + 1
    
    if (length(peak_indices) == 0) {
      # If no peaks detected, use global maximum
      peak_indices <- which.max(y_values)
    }
    
    peak_x <- x_values[peak_indices]
    peak_y <- y_values[peak_indices]
    
    # Calculate peak prominence to filter noise
    if (length(peak_x) > 1) {
      prominences <- numeric(length(peak_x))
      for (i in 1:length(peak_x)) {
        # Simple prominence calculation
        left_min <- if (i > 1) min(peak_y[1:(i-1)]) else min(y_values)
        right_min <- if (i < length(peak_x)) min(peak_y[(i+1):length(peak_x)]) else min(y_values)
        prominences[i] <- peak_y[i] - max(left_min, right_min, na.rm = TRUE)
      }
      
      # Filter by prominence
      significant <- prominences > (max(peak_y) * min_peak_prominence)
      peak_x <- peak_x[significant]
      peak_y <- peak_y[significant]
    }
    
    # Remove duplicates and sort by x
    unique_peaks <- !duplicated(peak_x)
    peak_x <- peak_x[unique_peaks]
    peak_y <- peak_y[unique_peaks]
    
    if (length(peak_x) == 0) {
      return(list(
        first_peak_x = NULL,
        first_peak_y = NULL,
        highest_peak_x = NULL,
        highest_peak_y = NULL,
        all_peaks_x = NULL,
        all_peaks_y = NULL
      ))
    }
    
    # Sort by x to ensure first peak is truly the first
    sorted_indices <- order(peak_x)
    peak_x <- peak_x[sorted_indices]
    peak_y <- peak_y[sorted_indices]
    
    # Highest peak
    highest_idx <- which.max(peak_y)
    highest_peak_x <- peak_x[highest_idx]
    highest_peak_y <- peak_y[highest_idx]
    
    # First peak (leftmost significant peak)
    first_peak_x <- peak_x[1]
    first_peak_y <- peak_y[1]
    
    return(list(
      first_peak_x = first_peak_x,
      first_peak_y = first_peak_y,
      highest_peak_x = highest_peak_x,
      highest_peak_y = highest_peak_y,
      all_peaks_x = peak_x,
      all_peaks_y = peak_y
    ))
  }
  
  # Detect peaks properly
  peak_info <- find_peaks_properly(x_filtered, y_filtered)
  first_peak_x_bin <- peak_info$first_peak_x
  highest_peak_x_bin <- peak_info$highest_peak_x
  first_peak_x <- estimate_continuous_peak(
    x_filtered, y_filtered, first_peak_x_bin
  )
  highest_peak_x <- estimate_continuous_peak(
    x_filtered, y_filtered, highest_peak_x_bin
  )
  first_peak_y <- estimate_reference_height(
    x_filtered, y_filtered, first_peak_x
  )
  highest_peak_y <- estimate_reference_height(
    x_filtered, y_filtered, highest_peak_x
  )
  
  cat("First peak position:", first_peak_x, ", height:", first_peak_y, "\n")
  cat("Highest peak position:", highest_peak_x, ", height:", highest_peak_y, "\n")
  
  # Use highest peak as main peak for symmetry processing
  first_main_peak_x <- highest_peak_x
  first_main_peak_y <- highest_peak_y
  
  # Define main peak focus region
  main_peak_range_ratio <- main_peak_ratio
  main_peak_left <- round(first_main_peak_x * (1 - main_peak_range_ratio))
  main_peak_right <- round(first_main_peak_x * (1 + main_peak_range_ratio))
  
  # Ensure range is within data bounds
  main_peak_left <- max(main_peak_left, min(x_filtered))
  main_peak_right <- min(main_peak_right, max(x_filtered))
  
  cat("Main peak focus region:", main_peak_left, "-", main_peak_right, "\n")
  
  # Extract main peak region data
  main_peak_indices <- which(x_filtered >= main_peak_left & x_filtered <= main_peak_right)
  main_peak_x <- x_filtered[main_peak_indices]
  main_peak_y <- y_filtered[main_peak_indices]
  
  # Enhanced symmetry check and data completion
  symmetric_main_peak_data <- complete_main_peak_symmetry_enhanced(
    main_peak_x, main_peak_y, first_main_peak_x, first_main_peak_y, 
    symmetry_tolerance
  )
  
  main_peak_x_symmetric <- symmetric_main_peak_data$x
  main_peak_y_symmetric <- symmetric_main_peak_data$y
  
  cat("Processed main peak data points:", length(main_peak_x_symmetric), "\n")
  
  # Determine fit right boundary
  focus_range_end <- first_main_peak_x * (n + 2)
  
  # Smooth data for peak detection
  y_smooth_peaks <- stats::filter(y_filtered, rep(1, 5)/5, sides = 2)
  y_smooth_peaks[is.na(y_smooth_peaks)] <- y_filtered[is.na(y_smooth_peaks)]
  
  # Find local maxima with improved filtering
  is_peak <- diff(sign(diff(y_smooth_peaks))) < 0
  peak_indices <- which(is_peak) + 1
  
  # Filter peaks: remove noise and insignificant peaks
  peak_x <- x_filtered[peak_indices]
  peak_y <- y_filtered[peak_indices]
  
  # First limit to focus range
  in_focus_range <- peak_x <= focus_range_end
  peak_x_focused <- peak_x[in_focus_range]
  peak_y_focused <- peak_y[in_focus_range]
  
  # 1. Remove peaks that are too close to each other (keep the higher one)
  min_peak_distance <- first_peak_x * 0.5  # Minimum distance between peaks
  if (length(peak_x_focused) > 1) {
    keep_peak <- rep(TRUE, length(peak_x_focused))
    for (i in 2:length(peak_x_focused)) {
      if (peak_x_focused[i] - peak_x_focused[i-1] < min_peak_distance) {
        # Keep the higher peak
        if (peak_y_focused[i] > peak_y_focused[i-1]) {
          keep_peak[i-1] <- FALSE
        } else {
          keep_peak[i] <- FALSE
        }
      }
    }
    significant_peaks <- keep_peak
  } else {
    significant_peaks <- rep(TRUE, length(peak_x_focused))  
  }
  
  # 2. Improved significance filtering: require peak height to be at least 1.5x of adjacent valleys
  if (length(peak_x_focused) > 1) {
    keep_peak_significant <- rep(TRUE, length(peak_x_focused))
    for (i in seq_along(peak_x_focused)) {
      # Skip main peak from significance filtering
      if (abs(peak_x_focused[i] - highest_peak_x_bin) <= 0.5) next
      
      # Find left valley
      left_valley <- if (i > 1) {
        # Look for minimum between current peak and previous peak
        valley_left_indices <- which(x_filtered < peak_x_focused[i] & x_filtered > peak_x_focused[i-1])
        if (length(valley_left_indices) > 0) {
          min(y_filtered[valley_left_indices])
        } else {
          # If no data between peaks, use the lower of the two adjacent points
          min(y_filtered[which(x_filtered == peak_x_focused[i-1])], 
              y_filtered[which(x_filtered == peak_x_focused[i])])
        }
      } else {
        # For first peak, look for minimum to the left
        valley_left_indices <- which(x_filtered < peak_x_focused[i] & x_filtered >= min(x_filtered))
        if (length(valley_left_indices) > 0) {
          min(y_filtered[valley_left_indices])
        } else {
          y_filtered[which(x_filtered == peak_x_focused[i])] * 0.5  # Conservative estimate
        }
      }
      
      # Find right valley
      right_valley <- if (i < length(peak_x_focused)) {
        # Look for minimum between current peak and next peak
        valley_right_indices <- which(x_filtered > peak_x_focused[i] & x_filtered < peak_x_focused[i+1])
        if (length(valley_right_indices) > 0) {
          min(y_filtered[valley_right_indices])
        } else {
          # If no data between peaks, use the lower of the two adjacent points
          min(y_filtered[which(x_filtered == peak_x_focused[i])], 
              y_filtered[which(x_filtered == peak_x_focused[i+1])])
        }
      } else {
        # For last peak, look for minimum to the right
        valley_right_indices <- which(x_filtered > peak_x_focused[i] & x_filtered <= max(x_filtered))
        if (length(valley_right_indices) > 0) {
          min(y_filtered[valley_right_indices])
        } else {
          y_filtered[which(x_filtered == peak_x_focused[i])] * 0.5  # Conservative estimate
        }
      }
      
      # Use the higher valley as reference (more conservative filtering)
      max_valley <- max(left_valley, right_valley, na.rm = TRUE)
      
      # If peak height is less than 1.05x of the higher valley, filter it out
      if (peak_y_focused[i] < max_valley * 1.05) {
        keep_peak_significant[i] <- FALSE
      }
    }
    significant_peaks <- significant_peaks & keep_peak_significant
  }
  
  # Apply filtering
  peak_x_final <- peak_x_focused[significant_peaks]
  peak_y_final <- peak_y_focused[significant_peaks]
  continuous_peak_x <- vapply(peak_x_final, function(value) {
    estimate_continuous_peak(x_filtered, y_filtered, value)
  }, numeric(1))
  peak_x_final <- continuous_peak_x

  # A local maximum counts as an explicit biological peak only when it belongs
  # to one of the requested 1..n harmonic slots.  A weak fifth-harmonic bump
  # must not occupy the unresolved fourth component (the floridensis failure).
  # Keep all maxima in peak_x_final so an optional n+1 component can still use
  # them as nuisance-position evidence.
  harmonic_ratio <- continuous_peak_x / first_peak_x
  harmonic_slot <- round(harmonic_ratio)
  biological_mask <- is.finite(harmonic_ratio) &
    harmonic_slot >= 1L & harmonic_slot <= n &
    abs(harmonic_ratio - harmonic_slot) <= 0.35
  biological_indices <- which(biological_mask)
  if (length(biological_indices)) {
    # If noise creates two maxima in one slot, retain the one closest to its
    # biological harmonic rather than inflating the explicit-peak count.
    biological_indices <- unlist(lapply(
      sort(unique(harmonic_slot[biological_indices])), function(slot) {
        candidates <- biological_indices[
          harmonic_slot[biological_indices] == slot
        ]
        candidates[which.min(abs(harmonic_ratio[candidates] - slot))]
      }
    ), use.names = FALSE)
  }
  observed_peak_x <- continuous_peak_x[biological_indices]
  observed_peak_count <- length(observed_peak_x)
  
  cat("Final detected peak positions in focus range:\n")
  print(peak_x_final)
  cat("Number of peaks detected in focus range:", length(peak_x_final), "\n")
  cat("Biologically matched harmonic peak positions (1..n):\n")
  print(observed_peak_x)
  cat("Number of biologically matched peaks:", observed_peak_count, "\n")
  
  # Shoulder detection function
  detect_shoulder_peaks <- function(x_values, y_values, first_peak_x, first_peak_y, n_shoulders, start_x) {
    cat("Detecting shoulder peaks before first peak...\n")
    search_start <- start_x * 1.1
    search_end <- first_peak_x * 0.9
    search_indices <- which(x_values >= search_start & x_values <= search_end)
    
    if (length(search_indices) < 5) return(NULL)
    
    search_x <- x_values[search_indices]
    search_y <- y_values[search_indices]
    
    y_smooth <- stats::filter(search_y, rep(1, 3)/3, sides = 2)
    y_smooth[is.na(y_smooth)] <- search_y[is.na(y_smooth)]
    
    dy_dx <- diff(y_smooth) / diff(search_x)
    x_deriv <- search_x[-length(search_x)] + diff(search_x)/2
    
    if (length(dy_dx) > 2) {
      d2y_dx2 <- diff(dy_dx) / diff(x_deriv)
      x_2deriv <- x_deriv[-length(x_deriv)] + diff(x_deriv)/2
    } else {
      return(NULL)
    }
    
    shoulder_candidates <- which(d2y_dx2 > 0)
    if (length(shoulder_candidates) == 0) return(NULL)
    
    significant_shoulders <- numeric(0)
    for (candidate_idx in shoulder_candidates) {
      if (candidate_idx < 2 || candidate_idx >= length(d2y_dx2)) next
      
      prev_2deriv <- if (candidate_idx > 1) d2y_dx2[candidate_idx - 1] else NA
      
      is_shoulder <- TRUE
      if (!is.na(prev_2deriv) && prev_2deriv > 0) {
        found_negative <- FALSE
        for (i in (candidate_idx-1):max(1, candidate_idx-5)) {
          if (d2y_dx2[i] < 0) {
            found_negative <- TRUE
            break
          }
        }
        is_shoulder <- found_negative
      }
      
      if (is_shoulder) {
        shoulder_x <- x_2deriv[candidate_idx]
        if (shoulder_x >= search_start && shoulder_x <= search_end) {
          left_min <- min(y_smooth[1:candidate_idx], na.rm = TRUE)
          right_min <- min(y_smooth[candidate_idx:length(y_smooth)], na.rm = TRUE)
          prominence <- y_smooth[which(search_x == search_x[which.min(abs(search_x - shoulder_x))])] - max(left_min, right_min)
          
          if (prominence > max(y_smooth, na.rm = TRUE) * 0.005) {
            significant_shoulders <- c(significant_shoulders, shoulder_x)
          }
        }
      }
    }
    
    if (length(significant_shoulders) > 1) {
      significant_shoulders <- sort(significant_shoulders)
      keep <- rep(TRUE, length(significant_shoulders))
      for (i in 2:length(significant_shoulders)) {
        if (significant_shoulders[i] - significant_shoulders[i-1] < first_peak_x * 0.1) {
          keep[i] <- FALSE
        }
      }
      significant_shoulders <- significant_shoulders[keep]
    }
    
    shoulder_count <- length(significant_shoulders)
    if (shoulder_count > n_shoulders) shoulder_count <- n_shoulders
    
    if (shoulder_count > 0) {
      shoulder_positions <- numeric(shoulder_count)
      for (i in 1:shoulder_count) {
        shoulder_positions[i] <- first_peak_x / (shoulder_count + 1) * i
      }
      return(shoulder_positions)
    } else {
      return(NULL)
    }
  }
  
  # Helper function: determine data decay rate using expected peak positions
  is_slow_decaying_data <- function(y_values, x_values, first_peak_x, n_components, peak_x_final) {
    expected_positions <- first_peak_x * 1:n_components
    expected_y_values <- numeric(n_components)
    for (i in 1:n_components) {
      pos <- expected_positions[i]
      closest_index <- which.min(abs(x_values - pos))
      expected_y_values[i] <- y_values[closest_index]
    }
    
    fast_decay_detected <- FALSE
    for (i in 1:(n_components - 1)) {
      current_y <- expected_y_values[i]
      next_y <- expected_y_values[i + 1]
      
      if (current_y > 0) {
        decay_ratio <- next_y / current_y
        if (decay_ratio < 0.2) {
          fast_decay_detected <- TRUE
          break
        }
        if (next_y < max(y_values) * 0.01) {
          fast_decay_detected <- TRUE
          break
        }
      }
    }
    
    if (fast_decay_detected) return(FALSE) else return(TRUE)
  }
  
  # Determine truncation point
  if (length(peak_x_final) < n) {
    missing_peaks <- n - length(peak_x_final)
    shoulder_peaks <- detect_shoulder_peaks(x_filtered, y_filtered, first_peak_x, first_peak_y, missing_peaks, full_data[start_index, 1])
    
    if (!is.null(shoulder_peaks) && length(shoulder_peaks) > 0) {
      peak_x_final <- sort(c(peak_x_final, shoulder_peaks))
      peak_y_final <- c(peak_y_final, rep(NA, length(shoulder_peaks)))
      for (i in 1:length(peak_x_final)) {
        peak_idx <- which.min(abs(x_filtered - peak_x_final[i]))
        peak_y_final[i] <- y_filtered[peak_idx]
      }
      
      first_peak_x <- peak_x_final[1]
      first_peak_y <- peak_y_final[1]
      highest_idx <- which.max(peak_y_final)
      first_main_peak_x <- peak_x_final[highest_idx]
      first_main_peak_y <- peak_y_final[highest_idx]
    }
    
    if (length(peak_x_final) < n) {
      if (is_slow_decaying_data(y_filtered, x_filtered, first_peak_x, n, peak_x_final)) {
        peak_multiplier <- n + 2
      } else {
        peak_multiplier <- n + 1
      }
      right_value <- peak_x_final[1] * peak_multiplier  
    } else {
      nth_peak_x <- peak_x_final[n]
      search_start <- nth_peak_x
      search_end <- nth_peak_x + first_main_peak_x
      search_indices <- which(x_filtered >= search_start & x_filtered <= search_end)
      
      if (length(search_indices) > 0) {
        search_x <- x_filtered[search_indices]
        search_y <- y_filtered[search_indices]
        nth_peak_y <- y_filtered[which.min(abs(x_filtered - nth_peak_x))]
        threshold_y <- nth_peak_y * 0.1
        significant_drop_indices <- which(search_y < threshold_y)
        
        if (length(significant_drop_indices) > 0) {
          right_value <- search_x[significant_drop_indices[1]]
        } else {
          right_value <- nth_peak_x + first_main_peak_x
        }
      } else {
        right_value <- nth_peak_x * 1.5
      }
    }
  } else {
    nth_peak_x <- peak_x_final[n]
    if (length(peak_x_final) == n) {
      search_start <- nth_peak_x
      search_end <- nth_peak_x + first_main_peak_x
      search_indices <- which(x_filtered >= search_start & x_filtered <= search_end)
      
      if (length(search_indices) > 0) {
        search_x <- x_filtered[search_indices]
        search_y <- y_filtered[search_indices]
        nth_peak_y <- y_filtered[which.min(abs(x_filtered - nth_peak_x))]
        threshold_y <- nth_peak_y * 0.1
        significant_drop_indices <- which(search_y < threshold_y)
        
        if (length(significant_drop_indices) > 0) {
          right_value <- search_x[significant_drop_indices[1]]
        } else {
          right_value <- nth_peak_x + first_main_peak_x
        }
      } else {
        right_value <- nth_peak_x * 1.5
      }
    } else {
      search_end <- nth_peak_x + first_main_peak_x
      search_indices <- which(x_filtered >= nth_peak_x & x_filtered <= search_end)
      if (length(search_indices) > 0) {
        min_index_in_range <- search_indices[which.min(y_filtered[search_indices])]
        right_value <- x_filtered[min_index_in_range]
      } else {
        right_value <- nth_peak_x * 1.5
      }
    }
  }
  
  right_value <- min(right_value, max(x_filtered))
  right_index <- which(full_data[, 1] >= right_value)[1]
  
  if (!is.na(right_index)) {
    fit_range <- start_index:right_index
  } else {
    fit_range <- start_index:nrow(full_data)
  }
  
  cat("Fit range: frequency", full_data[start_index, 1], "to", right_value, "\n")
  
  # Extract data for fitting
  fitted_data <- full_data[fit_range, ]
  x_fit <- fitted_data[, 1]
  y_fit <- fitted_data[, 2]
  
  # Create processed fit data with weights
  processed_fit_data <- create_processed_fit_data_with_weight(
    x_fit, y_fit, main_peak_x_symmetric, main_peak_y_symmetric,
    main_peak_left, main_peak_right, main_peak_weight
  )
  
  x_fit_processed <- processed_fit_data$x
  y_fit_processed <- processed_fit_data$y
  weight_fit_processed <- processed_fit_data$weights
  
  cat("Main peak region weight:", main_peak_weight, "\n")
  cat("Main peak data points:", length(main_peak_x_symmetric), "\n")
  cat("Total data points for fitting:", length(x_fit_processed), "\n")
  
  # Select whether to use calibration based on parameters
  if (use_calibration) {
    cat("Using EM algorithm with post-fitting calibration\n")
    # Execute EM algorithm with calibration
    EM_result <- EM_algorithm_with_calibration(
      x_fit_processed, y_fit_processed, weight_fit_processed, n, 
      first_main_peak_x, first_peak_x, peak_x_final,
      main_peak_x_symmetric, main_peak_y_symmetric,
      max_sigma_ratio = max_sigma_ratio,
      sigma_mode = sigma_mode,
      position_mode = position_mode,
      observed_peak_count = observed_peak_count,
      max_iter = em_max_iter,
      tolerance = em_tolerance,
      extra_component_mode = extra_component_mode
    )
  } else {
    cat("Using original EM algorithm without calibration\n")
    # Execute original EM algorithm without calibration
    EM_result <- EM_algorithm(
      x_fit_processed, y_fit_processed, weight_fit_processed, n, 
      first_main_peak_x, first_peak_x, peak_x_final,
      max_sigma_ratio = max_sigma_ratio,
      sigma_mode = sigma_mode,
      position_mode = position_mode,
      observed_peak_count = observed_peak_count,
      max_iter = em_max_iter,
      tolerance = em_tolerance,
      extra_component_mode = extra_component_mode
    )
  }
  
  alpha <- EM_result$alpha
  miu <- EM_result$miu
  sigma <- EM_result$sigma
  curve_refinement <- refine_gmm_curve_widths(
    x = x_fit, y = y_fit, alpha = alpha, miu = miu, sigma = sigma,
    total_samples_fit = sum(y_fit_processed), first_peak_x = first_peak_x,
    sigma_reliability = EM_result$sigma_reliability,
    sigma_mode = sigma_mode
  )
  alpha <- curve_refinement$alpha
  sigma <- curve_refinement$sigma
  if (!identical(sigma_mode, "legacy")) {
    alpha <- stabilize_biological_alpha(alpha, n)
  }
  EM_result$alpha <- alpha
  EM_result$sigma <- sigma
  if (!identical(sigma_mode, "legacy")) {
    cat("Direct curve width refinement:",
        if (isTRUE(curve_refinement$refined)) "applied" else "kept EM result",
        "\n")
    cat("Original-data curve scale:",
        round(curve_refinement$curve_scale, 6), "\n")
  }
  # The guarded comparison is attempted for every final fit.  A sample no
  # longer has to cross one brittle trigger threshold: baseline, standard
  # joint, width-balanced and flexible multi-start candidates compete under
  # the same original-curve and biological-sanity checks.
  joint_triggered <- !identical(joint_refinement, "off")
  joint_accepted <- FALSE
  joint_profile <- "baseline"
  joint_candidates <- NULL
  joint_score_before <- NA_real_
  joint_score_after <- NA_real_
  if (joint_triggered) {
    baseline_joint_score <- calculate_gmm_fit_score(
      x_observed = x_fit, y_observed = y_fit,
      alpha = alpha, miu = miu, sigma = sigma,
      total_samples_fit = sum(y_fit_processed),
      main_peak_left = main_peak_left, main_peak_right = main_peak_right,
      sigma_mode = sigma_mode,
      sigma_reliability = EM_result$sigma_reliability,
      detected_peak_x = observed_peak_x,
      first_peak_x = first_peak_x
    )
    joint_score_before <- baseline_joint_score$score
    candidate_values <- list(baseline = list(
      alpha = alpha, miu = miu, sigma = sigma,
      score_info = baseline_joint_score, refined = TRUE,
      eligible = TRUE
    ))
    overlap_case <- n >= 5L && observed_peak_count < n &&
      mean(EM_result$sigma_reliability, na.rm = TRUE) < 0.35
    low_overlap_case <- n <= 3L && observed_peak_count < n &&
      mean(EM_result$sigma_reliability, na.rm = TRUE) < 0.35
    total_peak_case <- n >= 2L && observed_peak_count >= n &&
      baseline_joint_score$primary_total_peak_position_error > 0.02
    joint_profiles <- c("standard", "width_balanced", "flexible")
    # The two extra optimisers are residual-triggered, not species-triggered.
    # Clean fits avoid unnecessary work; any sample with the corresponding
    # visible defect can enter the same guarded comparison.
    if (baseline_joint_score$left_overshoot_error > 0.002 ||
        baseline_joint_score$left_relative_error > 0.06) {
      joint_profiles <- c(joint_profiles, "flank_balanced")
    }
    if (baseline_joint_score$region_max_nrmse > 0.01) {
      joint_profiles <- c(joint_profiles, "region_balanced")
    }
    if (overlap_case && baseline_joint_score$global_overshoot_nrmse > 0.006) {
      joint_profiles <- c(joint_profiles, "overlap_balanced")
    }
    if (low_overlap_case &&
        baseline_joint_score$global_overshoot_nrmse > 0.006) {
      joint_profiles <- c(joint_profiles, "low_ploidy_balanced")
    }
    if (total_peak_case) {
      joint_profiles <- c(joint_profiles, "total_peak_aligned")
    }
    for (profile in joint_profiles) {
      candidate_start <- if (identical(profile, "low_ploidy_balanced") &&
          !is.null(candidate_values$flank_balanced) &&
          isTRUE(candidate_values$flank_balanced$refined)) {
        candidate_values$flank_balanced
      } else list(alpha = alpha, miu = miu, sigma = sigma)
      candidate <- refine_gmm_curve_joint(
        x = x_fit, y = y_fit,
        alpha = candidate_start$alpha,
        miu = candidate_start$miu,
        sigma = candidate_start$sigma,
        total_samples_fit = sum(y_fit_processed), first_peak_x = first_peak_x,
        sigma_reliability = EM_result$sigma_reliability,
        # Synthetic shoulder positions may guide the optimiser, but they are
        # deliberately excluded from n+1 decisions and fit diagnostics.
        detected_peak_x = peak_x_final, profile = profile,
        max_iterations = if (profile %in%
          c("flank_balanced", "region_balanced", "overlap_balanced")) {
          350L
        } else if (identical(profile, "low_ploidy_balanced")) {
          350L
        } else if (identical(profile, "total_peak_aligned")) {
          450L
        } else 600L
      )
      candidate_score <- if (isTRUE(candidate$refined)) {
        candidate$alpha <- stabilize_biological_alpha(candidate$alpha, n)
        calculate_gmm_fit_score(
          x_observed = x_fit, y_observed = y_fit,
          alpha = candidate$alpha, miu = candidate$miu,
          sigma = candidate$sigma,
          total_samples_fit = sum(y_fit_processed),
          main_peak_left = main_peak_left, main_peak_right = main_peak_right,
          sigma_mode = sigma_mode,
          sigma_reliability = EM_result$sigma_reliability,
          detected_peak_x = observed_peak_x,
          first_peak_x = first_peak_x
        )
      } else NULL
      score_improvement <- if (!is.null(candidate_score)) {
        (joint_score_before - candidate_score$score) /
          max(joint_score_before, 1e-8)
      } else -Inf
      main_peak_guard <- !is.null(candidate_score) &&
        candidate_score$main_peak_nrmse <=
          baseline_joint_score$main_peak_nrmse * 1.05 + 0.001
      height_guard <- !is.null(candidate_score) &&
        candidate_score$clear_peak_height_error <=
          baseline_joint_score$clear_peak_height_error + 0.006
      peak_ratio_floor <- if (is.finite(
        baseline_joint_score$peak_height_ratio
      )) min(0.985, baseline_joint_score$peak_height_ratio - 0.003) else 0.985
      peak_ratio_guard <- !is.null(candidate_score) &&
        is.finite(candidate_score$peak_height_ratio) &&
        candidate_score$peak_height_ratio >= peak_ratio_floor
      clear_underfit_guard <- !is.null(candidate_score) &&
        candidate_score$clear_peak_underfit <=
          max(0.02, baseline_joint_score$clear_peak_underfit + 0.003)
      new_balanced_profile <- profile %in%
        c("flank_balanced", "region_balanced", "overlap_balanced",
          "low_ploidy_balanced", "total_peak_aligned")
      curve_guard <- !is.null(candidate_score) &&
        candidate_score$nrmse <= baseline_joint_score$nrmse *
          if (new_balanced_profile) 1.02 else 1.01
      position_guard <- !is.null(candidate_score) &&
        (!new_balanced_profile ||
          candidate_score$primary_peak_position_error <= max(
            0.03,
            baseline_joint_score$primary_peak_position_error + 0.015
          ))
      local_guard <- !is.null(candidate_score)
      if (identical(profile, "flank_balanced")) {
        local_guard <- candidate_score$left_relative_error <=
          baseline_joint_score$left_relative_error * 0.92 + 0.002
      } else if (identical(profile, "region_balanced")) {
        local_guard <- candidate_score$region_max_nrmse <=
          baseline_joint_score$region_max_nrmse * 0.94 + 0.0005
      } else if (identical(profile, "overlap_balanced")) {
        local_guard <- candidate_score$global_overshoot_nrmse <=
          baseline_joint_score$global_overshoot_nrmse * 0.90
      } else if (identical(profile, "low_ploidy_balanced")) {
        local_guard <- candidate_score$global_overshoot_nrmse <=
          baseline_joint_score$global_overshoot_nrmse * 0.90
      } else if (identical(profile, "total_peak_aligned")) {
        local_guard <- candidate_score$primary_total_peak_position_error <=
          min(0.02,
              baseline_joint_score$primary_total_peak_position_error * 0.70)
      }
      eligible <- is.finite(score_improvement) &&
        score_improvement >= if (new_balanced_profile) 0.005 else 0.015 &&
        main_peak_guard && height_guard && peak_ratio_guard &&
        clear_underfit_guard && curve_guard && position_guard && local_guard

      # In high-ploidy histograms, several unresolved harmonics often merge
      # into one broad late peak.  The ordinary clear-peak guard can then reject
      # a candidate that visibly improves both wings merely because the height
      # of one weak, overlapping harmonic changes by about one per cent.  Allow
      # only a tightly guarded Pareto exception; explicit-peak and low-ploidy
      # samples never enter this branch.
      candidate_width_ratio <- if (!is.null(candidate_score) &&
          length(candidate$sigma) > 1L) {
        ordered_sigma <- candidate$sigma[order(candidate$miu)]
        max(pmax(
          ordered_sigma[-1L] / ordered_sigma[-length(ordered_sigma)],
          ordered_sigma[-length(ordered_sigma)] / ordered_sigma[-1L]
        ))
      } else Inf
      overlap_pareto <- profile %in%
        c("region_balanced", "overlap_balanced") &&
        overlap_case && !is.null(candidate_score) &&
        candidate_score$nrmse <= baseline_joint_score$nrmse * 1.01 &&
        candidate_score$global_overshoot_nrmse <=
          baseline_joint_score$global_overshoot_nrmse * 0.85 &&
        candidate_score$region_max_nrmse <=
          baseline_joint_score$region_max_nrmse * 0.85 &&
        candidate_score$main_peak_nrmse <=
          baseline_joint_score$main_peak_nrmse * 1.05 + 0.001 &&
        candidate_score$primary_peak_position_error <= max(
          0.03,
          baseline_joint_score$primary_peak_position_error + 0.015
        ) &&
        candidate_score$clear_peak_height_error <= 0.015 &&
        candidate_score$clear_peak_underfit <= 0.015 &&
        is.finite(candidate_score$peak_height_ratio) &&
        candidate_score$peak_height_ratio >= 0.98 &&
        candidate_score$peak_height_ratio <= 1.015 &&
        candidate_width_ratio <= 1.35
      eligible <- eligible || overlap_pareto
      low_ploidy_pareto <- identical(profile, "low_ploidy_balanced") &&
        low_overlap_case && !is.null(candidate_score) &&
        candidate_score$nrmse <= baseline_joint_score$nrmse * 1.01 &&
        candidate_score$global_overshoot_nrmse <=
          baseline_joint_score$global_overshoot_nrmse * 0.85 &&
        candidate_score$global_underfit_nrmse <=
          baseline_joint_score$global_underfit_nrmse * 1.10 + 0.001 &&
        candidate_score$region_max_nrmse <=
          baseline_joint_score$region_max_nrmse * 0.95 + 0.0005 &&
        candidate_score$main_peak_nrmse <=
          baseline_joint_score$main_peak_nrmse * 1.05 + 0.001 &&
        candidate_score$primary_peak_position_error <= max(
          0.03,
          baseline_joint_score$primary_peak_position_error + 0.015
        ) &&
        candidate_score$clear_peak_height_error <= max(
          0.015,
          baseline_joint_score$clear_peak_height_error + 0.003
        ) &&
        candidate_score$clear_peak_underfit <= max(
          0.015,
          baseline_joint_score$clear_peak_underfit + 0.003
        ) &&
        is.finite(candidate_score$peak_height_ratio) &&
        candidate_score$peak_height_ratio >= 0.985 &&
        candidate_score$peak_height_ratio <= 1.015 &&
        candidate_score$width_identifiability_penalty <=
          baseline_joint_score$width_identifiability_penalty + 0.002
      eligible <- eligible || low_ploidy_pareto
      total_peak_pareto <- identical(profile, "total_peak_aligned") &&
        total_peak_case && !is.null(candidate_score) &&
        candidate_score$primary_total_peak_position_error <=
          min(0.02,
              baseline_joint_score$primary_total_peak_position_error * 0.70) &&
        candidate_score$total_peak_position_error <=
          baseline_joint_score$total_peak_position_error * 0.80 &&
        candidate_score$nrmse <= baseline_joint_score$nrmse * 1.03 &&
        candidate_score$main_peak_nrmse <=
          baseline_joint_score$main_peak_nrmse * 1.05 + 0.001 &&
        candidate_score$primary_peak_position_error <= max(
          0.04,
          baseline_joint_score$primary_peak_position_error + 0.015
        ) &&
        candidate_score$clear_peak_height_error <=
          baseline_joint_score$clear_peak_height_error + 0.0075 &&
        candidate_score$clear_peak_underfit <= max(
          0.02,
          baseline_joint_score$clear_peak_underfit + 0.004
        ) &&
        is.finite(candidate_score$peak_height_ratio) &&
        candidate_score$peak_height_ratio >= 0.985 &&
        candidate_score$peak_height_ratio <= 1.02 &&
        candidate_score$region_max_nrmse <=
          baseline_joint_score$region_max_nrmse * 1.20 + 0.001 &&
        candidate_score$global_overshoot_nrmse <=
          baseline_joint_score$global_overshoot_nrmse * 1.15 + 0.002 &&
        candidate_score$global_underfit_nrmse <=
          baseline_joint_score$global_underfit_nrmse * 1.15 + 0.002 &&
        candidate_score$width_identifiability_penalty <=
          baseline_joint_score$width_identifiability_penalty + 0.004
      eligible <- eligible || total_peak_pareto
      if (identical(profile, "total_peak_aligned")) {
        eligible <- total_peak_pareto
      }
      candidate$score_info <- candidate_score
      candidate$eligible <- eligible
      candidate$overlap_pareto <- overlap_pareto
      candidate$low_ploidy_pareto <- low_ploidy_pareto
      candidate$total_peak_pareto <- total_peak_pareto
      candidate$score_improvement <- score_improvement
      candidate_values[[profile]] <- candidate
    }
    for (width_name in c("width_soft", "width_medium")) {
      balance_strength <- if (identical(width_name, "width_soft")) 0.008 else 0.025
      candidate <- refine_gmm_width_compromise(
        x = x_fit, y = y_fit, alpha = alpha, miu = miu, sigma = sigma,
        total_samples_fit = sum(y_fit_processed), first_peak_x = first_peak_x,
        detected_peak_x = peak_x_final,
        balance_strength = balance_strength
      )
      candidate_score <- if (isTRUE(candidate$refined)) {
        candidate$alpha <- stabilize_biological_alpha(candidate$alpha, n)
        calculate_gmm_fit_score(
          x_observed = x_fit, y_observed = y_fit,
          alpha = candidate$alpha, miu = candidate$miu,
          sigma = candidate$sigma,
          total_samples_fit = sum(y_fit_processed),
          main_peak_left = main_peak_left, main_peak_right = main_peak_right,
          sigma_mode = sigma_mode,
          sigma_reliability = EM_result$sigma_reliability,
          detected_peak_x = observed_peak_x,
          first_peak_x = first_peak_x
        )
      } else NULL
      score_improvement <- if (!is.null(candidate_score)) {
        (joint_score_before - candidate_score$score) /
          max(joint_score_before, 1e-8)
      } else -Inf
      main_peak_guard <- !is.null(candidate_score) &&
        candidate_score$main_peak_nrmse <=
          baseline_joint_score$main_peak_nrmse * 1.05 + 0.001
      height_guard <- !is.null(candidate_score) &&
        candidate_score$clear_peak_height_error <=
          baseline_joint_score$clear_peak_height_error + 0.006
      peak_ratio_floor <- if (is.finite(
        baseline_joint_score$peak_height_ratio
      )) min(0.985, baseline_joint_score$peak_height_ratio - 0.003) else 0.985
      peak_ratio_guard <- !is.null(candidate_score) &&
        is.finite(candidate_score$peak_height_ratio) &&
        candidate_score$peak_height_ratio >= peak_ratio_floor
      clear_underfit_guard <- !is.null(candidate_score) &&
        candidate_score$clear_peak_underfit <=
          max(0.02, baseline_joint_score$clear_peak_underfit + 0.003)
      curve_guard <- !is.null(candidate_score) &&
        candidate_score$nrmse <= baseline_joint_score$nrmse * 1.03
      eligible <- is.finite(score_improvement) &&
        score_improvement >= 0.01 && main_peak_guard && height_guard &&
        peak_ratio_guard && clear_underfit_guard && curve_guard
      candidate$score_info <- candidate_score
      candidate$eligible <- eligible
      candidate$score_improvement <- score_improvement
      candidate_values[[width_name]] <- candidate
    }
    joint_candidates <- do.call(rbind, lapply(names(candidate_values),
      function(profile) {
        value <- candidate_values[[profile]]
        info <- value$score_info
        data.frame(
          profile = profile,
          refined = isTRUE(value$refined),
          eligible = isTRUE(value$eligible),
          overlap_pareto = isTRUE(value$overlap_pareto),
          low_ploidy_pareto = isTRUE(value$low_ploidy_pareto),
          total_peak_pareto = isTRUE(value$total_peak_pareto),
          score = if (is.null(info)) Inf else info$score,
          nrmse = if (is.null(info)) Inf else info$nrmse,
          main_peak_nrmse = if (is.null(info)) Inf else info$main_peak_nrmse,
          peak_height_error = if (is.null(info)) Inf else
            info$peak_height_error,
          peak_height_ratio = if (is.null(info)) NA_real_ else
            info$peak_height_ratio,
          clear_peak_height_error = if (is.null(info)) Inf else
            info$clear_peak_height_error,
          clear_peak_underfit = if (is.null(info)) Inf else
            info$clear_peak_underfit,
          peak_position_error = if (is.null(info)) Inf else
            info$peak_position_error,
          primary_peak_position_error = if (is.null(info)) Inf else
            info$primary_peak_position_error,
          total_peak_position_error = if (is.null(info)) Inf else
            info$total_peak_position_error,
          primary_total_peak_position_error = if (is.null(info)) Inf else
            info$primary_total_peak_position_error,
          primary_total_peak_signed_error = if (is.null(info)) Inf else
            info$primary_total_peak_signed_error,
          flank_overshoot_error = if (is.null(info)) Inf else
            info$flank_overshoot_error,
          left_overshoot_error = if (is.null(info)) Inf else
            info$left_overshoot_error,
          right_overshoot_error = if (is.null(info)) Inf else
            info$right_overshoot_error,
          left_relative_error = if (is.null(info)) Inf else
            info$left_relative_error,
          global_overshoot_nrmse = if (is.null(info)) Inf else
            info$global_overshoot_nrmse,
          global_underfit_nrmse = if (is.null(info)) Inf else
            info$global_underfit_nrmse,
          region_mean_nrmse = if (is.null(info)) Inf else
            info$region_mean_nrmse,
          region_max_nrmse = if (is.null(info)) Inf else
            info$region_max_nrmse,
          width_penalty = if (is.null(info)) Inf else
            info$width_identifiability_penalty,
          max_adjacent_sigma_ratio = if (length(value$sigma) > 1L) {
            max(pmax(value$sigma[-1L] / value$sigma[-length(value$sigma)],
                     value$sigma[-length(value$sigma)] / value$sigma[-1L]))
          } else 1,
          means = paste(round(value$miu, 5), collapse = ";"),
          sigmas = paste(round(value$sigma, 5), collapse = ";"),
          stringsAsFactors = FALSE
        )
      }))
    eligible_profiles <- names(candidate_values)[vapply(
      candidate_values, function(value) isTRUE(value$eligible), logical(1)
    )]
    eligible_scores <- vapply(eligible_profiles, function(profile) {
      candidate_values[[profile]]$score_info$score
    }, numeric(1))
    joint_profile <- eligible_profiles[which.min(eligible_scores)]

    # A severe rising-flank defect is a separate visible failure mode.  If the
    # flank candidate Pareto-dominates the composite-score winner on both the
    # local flank and the global curve, do not let an unrelated region term
    # overturn that improvement (the diploid shoulder failure).
    if ("flank_balanced" %in% eligible_profiles &&
        baseline_joint_score$left_relative_error > 0.04) {
      flank_info <- candidate_values$flank_balanced$score_info
      winner_info <- candidate_values[[joint_profile]]$score_info
      if (flank_info$left_relative_error <=
            0.92 * winner_info$left_relative_error &&
          flank_info$nrmse <= 1.02 * winner_info$nrmse) {
        joint_profile <- "flank_balanced"
      }
    }
    if ("region_balanced" %in% eligible_profiles &&
        isTRUE(candidate_values$region_balanced$overlap_pareto)) {
      joint_profile <- "region_balanced"
    }
    if ("overlap_balanced" %in% eligible_profiles &&
        isTRUE(candidate_values$overlap_balanced$overlap_pareto)) {
      overlap_info <- candidate_values$overlap_balanced$score_info
      reference_info <- candidate_values[[joint_profile]]$score_info
      if (overlap_info$global_overshoot_nrmse <=
            0.92 * reference_info$global_overshoot_nrmse &&
          overlap_info$nrmse <= 1.03 * reference_info$nrmse &&
          overlap_info$region_max_nrmse <=
            1.05 * reference_info$region_max_nrmse &&
          overlap_info$global_underfit_nrmse <=
            1.08 * reference_info$global_underfit_nrmse + 0.001) {
        joint_profile <- "overlap_balanced"
      }
    }
    if ("low_ploidy_balanced" %in% eligible_profiles &&
        isTRUE(candidate_values$low_ploidy_balanced$low_ploidy_pareto)) {
      low_info <- candidate_values$low_ploidy_balanced$score_info
      reference_info <- candidate_values[[joint_profile]]$score_info
      if (low_info$global_overshoot_nrmse <=
            0.92 * reference_info$global_overshoot_nrmse &&
          low_info$nrmse <= 1.03 * reference_info$nrmse &&
          low_info$region_max_nrmse <=
            1.05 * reference_info$region_max_nrmse &&
          low_info$global_underfit_nrmse <=
            1.10 * reference_info$global_underfit_nrmse + 0.001 &&
          low_info$primary_peak_position_error <= max(
            0.03,
            reference_info$primary_peak_position_error + 0.01
          )) {
        joint_profile <- "low_ploidy_balanced"
      }
    }
    if ("total_peak_aligned" %in% eligible_profiles &&
        isTRUE(candidate_values$total_peak_aligned$total_peak_pareto)) {
      total_info <- candidate_values$total_peak_aligned$score_info
      reference_info <- candidate_values[[joint_profile]]$score_info
      if (total_info$primary_total_peak_position_error <=
            0.75 * reference_info$primary_total_peak_position_error &&
          total_info$total_peak_position_error <=
            0.85 * reference_info$total_peak_position_error &&
          total_info$nrmse <= 1.03 * reference_info$nrmse &&
          total_info$main_peak_nrmse <=
            1.05 * reference_info$main_peak_nrmse + 0.001 &&
          total_info$region_max_nrmse <=
            1.20 * reference_info$region_max_nrmse + 0.001) {
        joint_profile <- "total_peak_aligned"
      }
    }
    selected_joint <- candidate_values[[joint_profile]]
    joint_accepted <- !identical(joint_profile, "baseline")
    joint_score_after <- selected_joint$score_info$score
    if (joint_accepted) {
      alpha <- selected_joint$alpha
      miu <- selected_joint$miu
      sigma <- selected_joint$sigma
      EM_result$alpha <- alpha
      EM_result$miu <- miu
      EM_result$sigma <- sigma
    }
    cat("Guarded joint candidate comparison: completed\n")
    cat("Selected refinement profile:", joint_profile, "\n")
    cat("Joint refinement accepted:", joint_accepted, "\n")
    cat("Joint score before/after:", round(joint_score_before, 6), "/",
        round(joint_score_after, 6), "\n")
  }
  converged <- EM_result$converged
  main_peak_component <- EM_result$main_peak_component
  
  # ================== New Code Block ==================
  # Extract extra fitting information
  fit_extra_component <- if (!is.null(EM_result$fit_extra_component)) {
    EM_result$fit_extra_component
  } else {
    FALSE
  }
  
  internal_n <- if (!is.null(EM_result$internal_n)) {
    EM_result$internal_n
  } else {
    n
  }
  
  # Print relevant information
  cat("[Extra Component Info]\n")
  cat("Fitted extra components:", fit_extra_component, "\n")
  cat("Internal fitted components count:", internal_n, "\n")
  cat("Original requested components count:", n, "\n")
  
  if (fit_extra_component) {
    cat(sprintf("Note: Extra components fitted (%d total), components from %d onwards will be shown in light gray\n", 
                internal_n, n + 1))
  }
  # ====================================================
  
  # Output parameters
  cat("\n[Parameters after EM iteration]\n")
  cat("Mixing coefficients alpha:", round(alpha, 4), "\n")
  cat("Means miu:", round(miu, 4), "\n")
  cat("Standard deviations sigma:", round(sigma, 4), "\n")
  
  if (!converged) {
    cat("Warning: EM algorithm did not converge within maximum iteration", EM_result$max_iter, "!\n")
  }
  
  # Calculate genome size
  # Preserve the original genome-size convention: use the biologically
  # assigned main component as the divisor (for example the 2x homozygous peak
  # in a heterozygous diploid), while the plot can still label the 1x harmonic.
  current_main_component <- EM_result$main_peak_component
  genome_size_result <- calculate_genome_size(
    x_fit_processed, y_fit_processed, alpha, miu, sigma, 
    full_data, start_index, right_value, 
    first_peak_component = current_main_component
  )
  
  total_kmers <- genome_size_result$total_kmers
  kmer_depth <- genome_size_result$kmer_depth
  genome_size_bp <- genome_size_result$genome_size_bp
  component_sizes_bp <- genome_size_result$component_sizes_bp # Get component sizes

  # Score the fitted curve against the ORIGINAL histogram.  This deliberately
  # does not score against the symmetry-completed data, otherwise an artificial
  # main peak could look better than the observation it is meant to explain.
  fit_score_info <- calculate_gmm_fit_score(
    x_observed = x_fit,
    y_observed = y_fit,
    alpha = alpha,
    miu = miu,
    sigma = sigma,
    total_samples_fit = sum(y_fit_processed),
    main_peak_left = main_peak_left,
    main_peak_right = main_peak_right,
    sigma_mode = sigma_mode,
    sigma_reliability = EM_result$sigma_reliability,
    detected_peak_x = observed_peak_x,
    first_peak_x = first_peak_x
  )
  
  # === Modification: Add residual component calculation ===
  n_components_fitted <- length(component_sizes_bp)
  total_fitted_size <- sum(component_sizes_bp)
  residual_size_bp <- genome_size_bp - total_fitted_size
  
  # Add residual component to component list
  if (residual_size_bp > 0) {
    component_sizes_bp <- c(component_sizes_bp, residual_size_bp)
  }
  
  cat("\n[Final calculation results]\n")
  cat("Estimated genome size:", genome_size_bp, "bp\n")
  cat("Kmer depth:", kmer_depth, "\n")
  
  # Print individual component sizes
  cat("\n[Component Sizes]\n")
  for(i in 1:length(component_sizes_bp)) {
    if (i <= n_components_fitted) {
      cat(sprintf("Component %d: %.2f Mb\n", i, component_sizes_bp[i] / 1e6))
    } else {
      cat(sprintf("Component %d (residual): %.2f Mb\n", i, component_sizes_bp[i] / 1e6))
    }
  }
  
  # Save only the final selected fit; automatic candidate fits skip this block.
  filename <- NULL
  if (isTRUE(save_plot)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- paste0("em_fitted_gmm_", file_basename, "_n", n, "_", timestamp, ".pdf")
    
    pdf_width <- 7.08 
    pdf_height <- 5.3 
    
    # Set PDF with consistent font
    pdf(filename, width = pdf_width, height = pdf_height, family = "Helvetica")
    draw_plot(
      x_fit, y_fit,  # Use original data for plotting, but processed data for fitting
      alpha, miu, sigma, genome_size_bp, component_sizes_bp, # Pass component sizes
      main_peak_x_symmetric, main_peak_y_symmetric,
      x_fit_processed, y_fit_processed,
      full_data, right_value,
      fit_extra_component = fit_extra_component,  # Pass extra parameter
      internal_n = internal_n,                     # Pass internal fitted count
      first_peak_x = first_peak_x,
      n = n,
      species_name = species_name
    )
    dev.off()
    cat("Plot saved as:", filename, "\n")
  }
  cat("Automatic fit score (lower is better):", round(fit_score_info$score, 6), "\n")
  cat("Algorithm execution completed!\n")
  
  return(invisible(list(
    genome_size_bp = genome_size_bp,
    component_sizes_bp = component_sizes_bp,
    alpha = alpha,
    miu = miu,
    sigma = sigma,
    converged = converged,
    iterations = if (!is.null(EM_result$iterations)) EM_result$iterations else NA_integer_,
    use_calibration = use_calibration,
    main_peak_weight = main_peak_weight,
    fit_extra_component = fit_extra_component,
    internal_n = internal_n,
    max_sigma_ratio = if (sigma_mode == "legacy") max_sigma_ratio else NA_real_,
    sigma_mode = sigma_mode,
    position_mode = position_mode,
    joint_refinement = joint_refinement,
    joint_triggered = joint_triggered,
    joint_accepted = joint_accepted,
    joint_profile = joint_profile,
    joint_candidates = joint_candidates,
    joint_score_before = joint_score_before,
    joint_score_after = joint_score_after,
    sigma_reliability = EM_result$sigma_reliability,
    curve_refined = isTRUE(curve_refinement$refined),
    curve_scale = curve_refinement$curve_scale,
    symmetry_tolerance = symmetry_tolerance,
    extra_component_mode = extra_component_mode,
    n = n,
    main_peak_ratio = main_peak_ratio,
    fit_score = fit_score_info$score,
    fit_metrics = fit_score_info,
    right_value = right_value,
    pdf_file = filename
  )))
}

#' Create processed fit data with optional main peak weighting
create_processed_fit_data_with_weight <- function(x_fit, y_fit, main_peak_x_symmetric, main_peak_y_symmetric,
                                                  main_peak_left, main_peak_right, main_peak_weight = 1.0) {
  
  # Remove main peak region from original data
  non_main_peak_indices <- which(x_fit < main_peak_left | x_fit > main_peak_right)
  x_non_main <- x_fit[non_main_peak_indices]
  y_non_main <- y_fit[non_main_peak_indices]
  
  # Merge processed main peak data with other data
  x_processed <- c(main_peak_x_symmetric, x_non_main)
  y_processed <- c(main_peak_y_symmetric, y_non_main)
  
  # Apply weights: main peak region gets higher weight, others get 1.0
  weights_processed <- c(
    rep(main_peak_weight, length(main_peak_x_symmetric)),  # Main peak region
    rep(1.0, length(x_non_main))                           # Other regions
  )
  
  # Sort by x values
  sorted_indices <- order(x_processed)
  
  return(list(
    x = x_processed[sorted_indices],
    y = y_processed[sorted_indices],
    weights = weights_processed[sorted_indices]
  ))
}

#' Enhanced symmetry completion using both sides information
complete_main_peak_symmetry_enhanced <- function(main_peak_x, main_peak_y, peak_center, peak_height, tolerance = 0.2) {
  
  # Find peak position
  peak_index <- which.min(abs(main_peak_x - peak_center))
  
  # Split left and right data
  left_indices <- which(main_peak_x < peak_center)
  right_indices <- which(main_peak_x > peak_center)
  
  left_x <- main_peak_x[left_indices]
  left_y <- main_peak_y[left_indices]
  right_x <- main_peak_x[right_indices]
  right_y <- main_peak_y[right_indices]
  
  cat("Main peak symmetry analysis - left points:", length(left_x), ", right points:", length(right_x), "\n")
  
  # If data volume differs significantly on both sides, perform symmetry completion
  if (abs(length(left_x) - length(right_x)) / max(length(left_x), length(right_x)) > tolerance) {
    
    if (length(left_x) > length(right_x)) {
      # Insufficient right data, complete with left symmetry
      completed_data <- complete_peak_side(left_x, left_y, right_x, right_y, peak_center, "right")
    } else {
      # Insufficient left data, complete with right symmetry  
      completed_data <- complete_peak_side(right_x, right_y, left_x, left_y, peak_center, "left")
    }
    
    processed_x <- completed_data$x
    processed_y <- completed_data$y
    cat("Symmetry completion done, added", length(processed_x) - length(main_peak_x), "data points\n")
    
  } else {
    # Data is basically symmetric, no need for completion
    processed_x <- main_peak_x
    processed_y <- main_peak_y
    cat("Main peak data is basically symmetric, no completion needed\n")
  }
  
  # Sort by x values
  sorted_indices <- order(processed_x)
  return(list(x = processed_x[sorted_indices], y = processed_y[sorted_indices]))
}

#' Complete one side of peak using symmetry from the other side
complete_peak_side <- function(source_x, source_y, target_x, target_y, peak_center, side) {
  
  # Generate symmetric x coordinates
  if (side == "right") {
    symmetric_x <- 2 * peak_center - source_x
    # Keep only points on the right side of the peak and not in existing data
    missing_x <- symmetric_x[symmetric_x > peak_center & !(symmetric_x %in% target_x)]
  } else {
    symmetric_x <- 2 * peak_center - source_x  
    # Keep only points on the left side of the peak and not in existing data
    missing_x <- symmetric_x[symmetric_x < peak_center & !(symmetric_x %in% target_x)]
  }
  
  # Use y values from source data (symmetric)
  source_y_mapping <- approx(source_x, source_y, xout = 2 * peak_center - missing_x, rule = 2)$y
  missing_y <- source_y_mapping
  
  # Merge data
  all_x <- c(source_x, target_x, missing_x)
  all_y <- c(source_y, target_y, missing_y)
  
  return(list(x = all_x, y = all_y))
}

# Measure how much the histogram itself identifies the width of each expected
# harmonic peak. Deep valleys, a local maximum close to the expected position,
# and adequate signal all increase reliability. Weak or overlapping peaks are
# deliberately given little freedom so they cannot absorb neighbouring peaks.
estimate_sigma_reliability <- function(x, y, first_peak_x, component_count) {
  reliability <- rep(0, component_count)
  valid <- is.finite(x) & is.finite(y) & y >= 0
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 5 || !is.finite(first_peak_x) || first_peak_x <= 0) {
    return(reliability)
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  y_smooth <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  y_smooth[!is.finite(y_smooth)] <- y[!is.finite(y_smooth)]
  global_height <- max(y_smooth, na.rm = TRUE)
  if (!is.finite(global_height) || global_height <= 0) return(reliability)

  for (j in seq_len(component_count)) {
    centre <- first_peak_x * j
    cell <- which(x >= centre - 0.5 * first_peak_x &
                  x <= centre + 0.5 * first_peak_x)
    peak_window <- which(x >= centre - 0.18 * first_peak_x &
                         x <= centre + 0.18 * first_peak_x)
    left_window <- which(x >= centre - 0.5 * first_peak_x &
                         x <= centre - 0.08 * first_peak_x)
    right_window <- which(x >= centre + 0.08 * first_peak_x &
                          x <= centre + 0.5 * first_peak_x)
    if (length(cell) < 3 || length(peak_window) == 0 ||
        length(left_window) == 0 || length(right_window) == 0) next

    local_index <- peak_window[which.max(y_smooth[peak_window])]
    peak_height <- y_smooth[local_index]
    valley_height <- max(min(y_smooth[left_window], na.rm = TRUE),
                         min(y_smooth[right_window], na.rm = TRUE))
    if (!is.finite(peak_height) || peak_height <= 0 ||
        !is.finite(valley_height)) next

    prominence <- max(0, min(1, (peak_height - valley_height) / peak_height))
    displacement <- abs(x[local_index] - centre) / first_peak_x
    alignment <- exp(-0.5 * (displacement / 0.18)^2)
    signal <- sqrt(max(0, min(1, peak_height / global_height)))
    reliability[j] <- max(0, min(1, sqrt(prominence) * alignment * signal))
  }
  reliability
}

# Apply either the historical hard adjacent-ratio rule or a data-adaptive soft
# width model. The adaptive model learns a power trend from identifiable peaks
# and shrinks only poorly identified peaks toward it. Bounds are relative to
# harmonic spacing and are numerical anti-collapse guards, not biological
# adjacent-width assumptions.
regularize_sigma_widths <- function(sigma, miu, first_peak_x, reliability,
                                    sigma_mode = "adaptive",
                                    max_sigma_ratio = 1.35) {
  if (length(sigma) < 1) return(sigma)
  sorted_indices <- order(miu)
  sigma_sorted <- as.numeric(sigma[sorted_indices])

  if (identical(sigma_mode, "legacy")) {
    if (length(sigma_sorted) > 1) {
      for (i in 2:length(sigma_sorted)) {
        sigma_sorted[i] <- min(
          max(sigma_sorted[i], sigma_sorted[i - 1] / max_sigma_ratio),
          sigma_sorted[i - 1] * max_sigma_ratio
        )
      }
    }
    sigma[sorted_indices] <- sigma_sorted
    return(sigma)
  }

  spacing <- max(as.numeric(first_peak_x), .Machine$double.eps)
  minimum_sigma <- max(spacing * 0.015, .Machine$double.eps)
  maximum_sigma <- spacing * 0.80
  sigma_sorted <- pmin(maximum_sigma, pmax(minimum_sigma, sigma_sorted))

  rel <- as.numeric(reliability)
  if (length(rel) != length(sigma_sorted)) rel <- rep(0, length(sigma_sorted))
  rel <- pmin(1, pmax(0, rel[seq_along(sigma_sorted)]))
  log_sigma <- log(sigma_sorted)
  weights <- 0.15 + 0.85 * rel

  if (identical(sigma_mode, "common") || length(sigma_sorted) == 1) {
    target <- rep(stats::weighted.mean(log_sigma, weights), length(log_sigma))
    regularized <- target
  } else {
    log_rank <- log(seq_along(sigma_sorted))
    rank_centre <- stats::weighted.mean(log_rank, weights)
    sigma_centre <- stats::weighted.mean(log_sigma, weights)
    denominator <- sum(weights * (log_rank - rank_centre)^2)
    slope <- if (denominator > .Machine$double.eps) {
      sum(weights * (log_rank - rank_centre) *
            (log_sigma - sigma_centre)) / denominator
    } else 0

    # Two genuinely visible peaks are enough to support a width trend. With
    # overlapping diploid peaks the evidence is small and the slope naturally
    # shrinks toward a common width; with sturgeon, clear early peaks anchor the
    # trend used for weak late peaks.
    slope_support <- min(1, sum(rel) / 2)
    slope <- slope * slope_support
    intercept <- stats::weighted.mean(log_sigma - slope * log_rank, weights)
    target <- intercept + slope * log_rank

    if (identical(sigma_mode, "trend")) {
      regularized <- target
    } else {
      # Clear peaks remain almost data-driven; uncertain peaks receive stronger
      # shrinkage. Retaining a moderate data contribution is essential for
      # sturgeon-like broad shoulders, where several weak harmonic components
      # jointly describe one visibly wide region.
      data_fraction <- 0.70 + 0.30 * rel
      regularized <- data_fraction * log_sigma +
        (1 - data_fraction) * target
    }
  }

  sigma_sorted <- pmin(maximum_sigma,
                       pmax(minimum_sigma, exp(regularized)))
  sigma[sorted_indices] <- sigma_sorted
  sigma
}

# Final curve refinement for automatic width modes. EM provides biologically
# constrained means and a stable starting point; this step directly minimises
# the displayed red-vs-grey curve error by adjusting only component weights and
# widths. Harmonic means are never changed.
refine_gmm_curve_widths <- function(x, y, alpha, miu, sigma,
                                    total_samples_fit, first_peak_x,
                                    sigma_reliability,
                                    sigma_mode = "adaptive",
                                    max_iterations = 300L) {
  if (identical(sigma_mode, "legacy")) {
    return(list(alpha = alpha, sigma = sigma, curve_scale = 1,
                refined = FALSE))
  }
  valid <- is.finite(x) & is.finite(y) & y >= 0
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  component_count <- length(alpha)
  if (length(x) < 5 || component_count < 1 ||
      !is.finite(total_samples_fit) || total_samples_fit <= 0) {
    return(list(alpha = alpha, sigma = sigma, curve_scale = 1,
                refined = FALSE))
  }

  x_steps <- diff(sort(unique(x)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps) > 0) stats::median(x_steps) else 1
  spacing <- max(first_peak_x, .Machine$double.eps)
  minimum_sigma <- max(spacing * 0.015, .Machine$double.eps)
  maximum_sigma <- spacing * 0.80
  sigma <- pmin(maximum_sigma, pmax(minimum_sigma, sigma))
  rel <- pmin(1, pmax(0, as.numeric(sigma_reliability)))
  if (length(rel) != component_count) rel <- rep(0, component_count)

  # With only two or three strongly overlapping components, widths and weights
  # are not separately identifiable. Freeze their regularised widths (so a
  # bird-like yellow component cannot expand), but still refine weights and the
  # global curve scale to remove the systematic peak-height deficit.
  freeze_widths <- component_count <= 3 && mean(rel) < 0.35

  alpha <- pmax(as.numeric(alpha), 1e-10)
  alpha <- alpha / sum(alpha)
  alpha_parameters <- if (component_count > 1) {
    log(alpha[seq_len(component_count - 1L)] / alpha[component_count])
  } else numeric(0)

  log_rank <- log(seq_len(component_count))
  width_start <- if (freeze_widths) numeric(0) else switch(
      sigma_mode,
      common = stats::weighted.mean(log(sigma), 0.15 + 0.85 * rel),
      trend = {
        if (component_count == 1) {
          c(log(sigma[1]), 0)
        } else {
          fit <- stats::lm.wfit(cbind(1, log_rank), log(sigma),
                                w = 0.15 + 0.85 * rel)
          as.numeric(fit$coefficients)
        }
      },
      adaptive = log(sigma),
      log(sigma)
    )

  y_scale <- max(y, na.rm = TRUE)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1
  point_weights <- 0.25 + 0.75 * sqrt(pmax(y, 0) / y_scale)
  log_y_scale <- max(log1p(y_scale), 1)

  # Use the original smoothed histogram to define every data-supported peak
  # height. Reliability continuously downweights shoulders/noise instead of
  # applying a brittle yes/no peak rule.
  y_smooth <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  y_smooth[!is.finite(y_smooth)] <- y[!is.finite(y_smooth)]
  peak_indices <- vapply(seq_len(component_count), function(j) {
    local <- which(x >= miu[j] - 0.18 * spacing &
                   x <= miu[j] + 0.18 * spacing)
    if (length(local) == 0) which.min(abs(x - miu[j])) else
      local[which.max(y[local])]
  }, integer(1))
  # The PDF's grey spline interpolates the original histogram points. Using a
  # five-bin average here systematically lowered narrow peaks in simulations.
  peak_observed <- y[peak_indices]
  peak_weights <- rel^2

  left_flank_penalty <- function(prediction, centres = miu) {
    rising_mask <- rep(FALSE, length(x))
    for (centre in centres) {
      rising_mask <- rising_mask |
        (x >= centre - 0.45 * spacing & x <= centre - 0.08 * spacing)
    }
    if (sum(rising_mask) < 3L) return(0)
    mean(pmax((prediction[rising_mask] - y_smooth[rising_mask]) /
                y_scale, 0)^2)
  }

  initial_design <- vapply(seq_len(component_count), function(j) {
    stats::dnorm(x, miu[j], sigma[j]) * total_samples_fit * bin_width
  }, numeric(length(x)))
  initial_prediction <- as.numeric(initial_design %*% alpha)
  scale_denominator <- sum(point_weights * initial_prediction^2)
  initial_scale <- if (scale_denominator > .Machine$double.eps) {
    sum(point_weights * y * initial_prediction) / scale_denominator
  } else 1
  initial_scale <- min(4, max(0.25, initial_scale))

  decode_parameters <- function(parameters) {
    alpha_parameter_count <- length(alpha_parameters)
    if (component_count == 1) {
      fitted_alpha <- 1
      next_parameter <- 1L
    } else {
      alpha_logits <- c(parameters[seq_len(alpha_parameter_count)], 0)
      alpha_logits <- alpha_logits - max(alpha_logits)
      fitted_alpha <- exp(alpha_logits) / sum(exp(alpha_logits))
      next_parameter <- alpha_parameter_count + 1L
    }
    fitted_scale <- exp(parameters[next_parameter])
    width_parameters <- if (next_parameter < length(parameters)) {
      parameters[(next_parameter + 1L):length(parameters)]
    } else numeric(0)
    fitted_sigma <- if (freeze_widths) sigma else switch(
      sigma_mode,
      common = rep(exp(width_parameters[1]), component_count),
      trend = exp(width_parameters[1] + width_parameters[2] * log_rank),
      adaptive = exp(width_parameters),
      exp(width_parameters)
    )
    list(alpha = fitted_alpha, sigma = fitted_sigma, scale = fitted_scale)
  }

  objective <- function(parameters) {
    decoded <- decode_parameters(parameters)
    fitted_sigma <- decoded$sigma
    if (any(!is.finite(fitted_sigma)) ||
        any(fitted_sigma < minimum_sigma) ||
        any(fitted_sigma > maximum_sigma)) return(1e6)
    design <- vapply(seq_len(component_count), function(j) {
      stats::dnorm(x, miu[j], fitted_sigma[j]) *
        total_samples_fit * bin_width
    }, numeric(length(x)))
    prediction <- as.numeric(design %*% decoded$alpha) * decoded$scale
    scaled_error <- (prediction - y) / y_scale
    curve_error <- sum(point_weights * scaled_error^2) / sum(point_weights)
    log_error <- mean(((log1p(pmax(prediction, 0)) - log1p(y)) /
                       log_y_scale)^2)
    reliable <- is.finite(peak_weights) & peak_weights > 0 &
      is.finite(peak_observed) & peak_observed > 0
    relative_peak_error <- if (any(reliable)) {
      (prediction[peak_indices[reliable]] - peak_observed[reliable]) /
        pmax(peak_observed[reliable], 0.05 * y_scale)
    } else numeric(0)
    peak_error <- if (length(relative_peak_error)) {
      stats::weighted.mean(relative_peak_error^2, peak_weights[reliable])
    } else 0
    peak_underfit_error <- if (length(relative_peak_error)) {
      stats::weighted.mean(
        pmax(-relative_peak_error, 0)^2, peak_weights[reliable]
      )
    } else 0

    width_penalty <- 0
    if (identical(sigma_mode, "adaptive") && component_count > 1) {
      log_width <- log(fitted_sigma)
      trend_fit <- stats::lm.wfit(cbind(1, log_rank), log_width,
                                  w = 0.15 + 0.85 * rel)
      trend_target <- as.numeric(cbind(1, log_rank) %*%
                                   trend_fit$coefficients)
      width_penalty <- mean((1 - rel) * (log_width - trend_target)^2)
      width_penalty <- width_penalty +
        mean((1 - rel) * pmax(fitted_sigma / spacing - 0.45, 0)^2)
      pair_reliability <- pmin(rel[-1], rel[-component_count])
      width_penalty <- width_penalty + mean(
        (1 - pair_reliability) * diff(log_width)^2
      )
    }
    0.55 * curve_error + 0.15 * log_error + 0.30 * peak_error +
      0.30 * peak_underfit_error +
      0.12 * left_flank_penalty(prediction) + 0.020 * width_penalty
  }

  initial <- c(alpha_parameters, log(initial_scale), width_start)
  alpha_lower <- rep(-12, length(alpha_parameters))
  alpha_upper <- rep(12, length(alpha_parameters))
  if (freeze_widths) {
    width_lower <- numeric(0)
    width_upper <- numeric(0)
  } else if (identical(sigma_mode, "common")) {
    width_lower <- log(minimum_sigma)
    width_upper <- log(maximum_sigma)
  } else if (identical(sigma_mode, "trend")) {
    maximum_slope <- if (component_count > 1) {
      log(maximum_sigma / minimum_sigma) / log(component_count)
    } else 0
    width_lower <- c(log(minimum_sigma), -maximum_slope)
    width_upper <- c(log(maximum_sigma), maximum_slope)
  } else {
    width_lower <- rep(log(minimum_sigma), component_count)
    width_upper <- rep(log(maximum_sigma), component_count)
  }

  refined <- tryCatch(
    stats::optim(
      initial, objective, method = "L-BFGS-B",
      lower = c(alpha_lower, log(0.25), width_lower),
      upper = c(alpha_upper, log(4), width_upper),
      control = list(maxit = as.integer(max_iterations), factr = 1e8)
    ),
    error = function(e) NULL
  )
  if (is.null(refined) || !is.finite(refined$value) ||
      refined$value > objective(initial) + 1e-10) {
    return(list(alpha = alpha, sigma = sigma, curve_scale = 1,
                refined = FALSE))
  }
  decoded <- decode_parameters(refined$par)
  # Absorb the fitted global scale into alpha so all existing genome-size,
  # scoring and plotting code uses the same original-data normalization.
  list(alpha = decoded$alpha * decoded$scale, sigma = decoded$sigma,
       curve_scale = decoded$scale, refined = TRUE,
       objective = refined$value, convergence = refined$convergence)
}

# Construct data-defined peak cells used by both scoring and optimisation.
# Flanks begin/end at local smoothed valleys, so low-frequency sequencing noise
# outside the peak is not mistaken for a Gaussian fitting error.
build_peak_shape_regions <- function(x, y, detected_peak_x, spacing) {
  detected_peak_x <- sort(unique(as.numeric(detected_peak_x)))
  detected_peak_x <- detected_peak_x[
    is.finite(detected_peak_x) & detected_peak_x >= min(x) &
      detected_peak_x <= max(x)
  ]
  if (!length(detected_peak_x) || !is.finite(spacing) || spacing <= 0) {
    return(NULL)
  }
  y_smooth <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  y_smooth[!is.finite(y_smooth)] <- y[!is.finite(y_smooth)]
  peak_indices <- vapply(detected_peak_x, function(value) {
    which.min(abs(x - value))
  }, integer(1))
  left_regions <- vector("list", length(detected_peak_x))
  right_regions <- vector("list", length(detected_peak_x))
  prominence <- numeric(length(detected_peak_x))
  for (i in seq_along(detected_peak_x)) {
    peak_index <- peak_indices[i]
    peak_value <- detected_peak_x[i]
    left_window <- which(x >= peak_value - 0.75 * spacing &
                         x <= peak_value)
    right_window <- which(x >= peak_value &
                          x <= peak_value + 0.75 * spacing)
    if (length(left_window)) {
      left_valley <- left_window[which.min(y_smooth[left_window])]
      left_regions[[i]] <- seq.int(min(left_valley, peak_index), peak_index)
    } else left_regions[[i]] <- peak_index
    if (length(right_window)) {
      right_valley <- right_window[which.min(y_smooth[right_window])]
      right_regions[[i]] <- seq.int(peak_index, max(right_valley, peak_index))
    } else right_regions[[i]] <- peak_index
    peak_height <- max(y_smooth[peak_index], 1)
    valley_height <- max(
      min(y_smooth[left_regions[[i]]], na.rm = TRUE),
      min(y_smooth[right_regions[[i]]], na.rm = TRUE)
    )
    prominence[i] <- pmin(1, pmax(0, (peak_height - valley_height) /
                                      peak_height))
  }
  list(
    x = x,
    peaks = detected_peak_x,
    peak_indices = peak_indices,
    left_regions = left_regions,
    right_regions = right_regions,
    weights = 0.20 + 0.80 * sqrt(prominence),
    y_smooth = y_smooth
  )
}

evaluate_peak_shape <- function(prediction, miu, regions, spacing) {
  if (is.null(regions) || !length(regions$peaks)) {
    return(list(position_error = 0, flank_overshoot_error = 0,
                primary_peak_position_error = 0,
                total_peak_position_error = 0,
                primary_total_peak_position_error = 0,
                primary_total_peak_signed_error = 0,
                left_overshoot_error = 0, right_overshoot_error = 0,
                left_relative_error = 0))
  }
  component_assignment <- vapply(regions$peaks, function(value) {
    which.min(abs(miu - value))
  }, integer(1))
  position_residual <- (miu[component_assignment] - regions$peaks) / spacing
  position_error <- sqrt(stats::weighted.mean(
    position_residual^2, regions$weights
  ))
  primary_peak_position_error <- abs(
    position_residual[which.max(regions$weights)]
  )
  total_peak_positions <- vapply(seq_along(regions$peaks), function(i) {
    local <- which(
      regions$x >= regions$peaks[i] - 0.30 * spacing &
        regions$x <= regions$peaks[i] + 0.30 * spacing
    )
    if (!length(local)) return(regions$peaks[i])
    regions$x[local[which.max(prediction[local])]]
  }, numeric(1))
  total_peak_residual <-
    (total_peak_positions - regions$peaks) / spacing
  total_peak_position_error <- sqrt(stats::weighted.mean(
    total_peak_residual^2, regions$weights
  ))
  primary_peak_index <- which.max(regions$weights)
  primary_total_peak_signed_error <- total_peak_residual[primary_peak_index]
  primary_total_peak_position_error <- abs(
    primary_total_peak_signed_error
  )
  peak_scale <- pmax(regions$y_smooth[regions$peak_indices],
                     0.05 * max(regions$y_smooth, na.rm = TRUE))
  side_error <- function(side_regions) {
    errors <- vapply(seq_along(side_regions), function(i) {
      indices <- side_regions[[i]]
      if (!length(indices)) return(0)
      mean(pmax(prediction[indices] - regions$y_smooth[indices], 0)) /
        peak_scale[i]
    }, numeric(1))
    stats::weighted.mean(errors, regions$weights)
  }
  left_error <- side_error(regions$left_regions)
  right_error <- side_error(regions$right_regions)

  # Global-peak normalisation hides a 30-50% error at a low shoulder because
  # that error is small compared with the tallest homozygous peak.  Measure the
  # rising flank once more relative to its local grey height, with a modest
  # floor to prevent near-zero noise bins from dominating the result.
  global_scale <- max(regions$y_smooth, na.rm = TRUE)
  local_floor <- max(1, 0.03 * global_scale)
  relative_left_errors <- vapply(
    seq_along(regions$left_regions), function(i) {
      indices <- regions$left_regions[[i]]
      if (!length(indices)) return(0)
      denominator <- pmax(regions$y_smooth[indices], local_floor)
      sqrt(mean(((prediction[indices] - regions$y_smooth[indices]) /
                   denominator)^2))
    }, numeric(1)
  )
  left_relative_error <- stats::weighted.mean(
    relative_left_errors, regions$weights
  )
  list(position_error = position_error,
       primary_peak_position_error = primary_peak_position_error,
       total_peak_position_error = total_peak_position_error,
       primary_total_peak_position_error =
         primary_total_peak_position_error,
       primary_total_peak_signed_error = primary_total_peak_signed_error,
       flank_overshoot_error = max(left_error, right_error),
       left_overshoot_error = left_error,
       right_overshoot_error = right_error,
       left_relative_error = left_relative_error)
}

# Balance the residual over harmonic cells.  A global RMSE can otherwise hide
# a visibly poor first or late peak behind the many well-fitted bins elsewhere.
evaluate_harmonic_regions <- function(x, reference_y, prediction, spacing) {
  if (!is.finite(spacing) || spacing <= 0 || length(x) < 3L) {
    return(list(mean_nrmse = 0, max_nrmse = 0))
  }
  y_scale <- max(reference_y, na.rm = TRUE)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1
  region_count <- max(1L, ceiling(max(x, na.rm = TRUE) / spacing))
  region_errors <- numeric(0)
  for (j in seq_len(region_count)) {
    indices <- which(
      x >= (j - 0.5) * spacing & x < (j + 0.5) * spacing
    )
    if (length(indices) < 3L ||
        max(reference_y[indices], na.rm = TRUE) < 0.01 * y_scale) next
    weights <- 0.25 + 0.75 * sqrt(pmax(reference_y[indices], 0) / y_scale)
    residual <- (prediction[indices] - reference_y[indices]) / y_scale
    region_errors <- c(region_errors, sqrt(
      sum(weights * residual^2) / sum(weights)
    ))
  }
  if (!length(region_errors)) return(list(mean_nrmse = 0, max_nrmse = 0))
  list(mean_nrmse = mean(region_errors), max_nrmse = max(region_errors))
}

# Joint curve refinement for low-identifiability histograms.  Unlike the
# preceding width-only pass, this candidate optimises ordered peak centres,
# widths and component areas together.  Spacings are expressed as smooth
# multiplicative deviations from the first-peak distance, so an overlapping
# sturgeon-like curve can move coherently without allowing individual
# components to chase noise.
refine_gmm_curve_joint <- function(x, y, alpha, miu, sigma,
                                   total_samples_fit, first_peak_x,
                                   sigma_reliability,
                                   detected_peak_x = numeric(0),
                                   profile = c("standard", "width_balanced",
                                               "flexible", "flank_balanced",
                                               "region_balanced",
                                               "overlap_balanced",
                                               "low_ploidy_balanced",
                                               "total_peak_aligned"),
                                   max_iterations = 600L) {
  profile <- match.arg(profile)
  valid <- is.finite(x) & is.finite(y) & y >= 0
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  component_count <- length(alpha)
  spacing <- as.numeric(first_peak_x)
  if (length(x) < 8L || component_count < 2L ||
      !is.finite(spacing) || spacing <= 0 ||
      !is.finite(total_samples_fit) || total_samples_fit <= 0) {
    return(list(alpha = alpha, miu = miu, sigma = sigma,
                refined = FALSE, objective_before = Inf,
                objective_after = Inf))
  }

  ordered <- order(miu)
  miu <- as.numeric(miu[ordered])
  sigma <- as.numeric(sigma[ordered])
  alpha <- as.numeric(alpha[ordered])
  rel <- as.numeric(sigma_reliability)
  if (length(rel) == component_count) rel <- rel[ordered] else
    rel <- rep(0, component_count)
  rel <- pmin(1, pmax(0, rel))

  x_steps <- diff(sort(unique(x)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  minimum_sigma_fraction <- switch(profile,
    standard = 0.06, width_balanced = 0.08, flexible = 0.05,
    flank_balanced = 0.06, region_balanced = 0.05,
    overlap_balanced = 0.05, low_ploidy_balanced = 0.05,
    total_peak_aligned = 0.05)
  maximum_sigma_fraction <- switch(profile,
    standard = 0.72, width_balanced = 0.60, flexible = 0.80,
    flank_balanced = 0.72, region_balanced = 0.80,
    overlap_balanced = 0.80, low_ploidy_balanced = 0.75,
    total_peak_aligned = 0.80)
  gap_lower <- switch(profile,
    standard = 0.82, width_balanced = 0.85, flexible = 0.76,
    flank_balanced = 0.80, region_balanced = 0.76,
    overlap_balanced = 0.76, low_ploidy_balanced = 0.78,
    total_peak_aligned = 0.78)
  gap_upper <- switch(profile,
    standard = 1.18, width_balanced = 1.15, flexible = 1.24,
    flank_balanced = 1.20, region_balanced = 1.24,
    overlap_balanced = 1.24, low_ploidy_balanced = 1.22,
    total_peak_aligned = 1.22)
  minimum_sigma <- max(minimum_sigma_fraction * spacing,
                       .Machine$double.eps)
  maximum_sigma <- maximum_sigma_fraction * spacing
  sigma <- pmin(maximum_sigma, pmax(minimum_sigma, sigma))

  alpha <- pmax(alpha, 1e-10)
  initial_scale <- sum(alpha)
  alpha_proportion <- alpha / initial_scale
  alpha_parameters <- log(
    alpha_proportion[seq_len(component_count - 1L)] /
      alpha_proportion[component_count]
  )
  first_lower <- if (identical(profile, "overlap_balanced")) {
    0.975
  } else if (profile %in% c("flexible", "region_balanced")) 0.84 else 0.88
  first_upper <- if (identical(profile, "overlap_balanced")) {
    1.025
  } else if (profile %in% c("flexible", "region_balanced")) 1.16 else 1.12
  initial_first <- min(max(miu[1], first_lower * spacing),
                       first_upper * spacing)
  initial_gap_ratio <- pmin(gap_upper,
                            pmax(gap_lower, diff(miu) / spacing))
  initial_log_sigma <- log(sigma)

  reference <- make_gmm_reference_curve(x, y)
  diagnostic_x <- reference$x
  diagnostic_y <- pmax(reference$y, 0)
  y_scale <- max(diagnostic_y, na.rm = TRUE)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1
  log_y_scale <- max(log1p(y_scale), 1)
  point_weights <- 0.20 + 0.80 * sqrt(pmax(diagnostic_y, 0) / y_scale)

  detected_peak_x <- sort(unique(as.numeric(detected_peak_x)))
  detected_peak_x <- detected_peak_x[
    is.finite(detected_peak_x) & detected_peak_x >= min(x) &
      detected_peak_x <= max(x)
  ]
  detected_indices <- if (length(detected_peak_x)) {
    vapply(detected_peak_x, function(value) {
      local <- which(diagnostic_x >= value - 0.18 * spacing &
                     diagnostic_x <= value + 0.18 * spacing)
      if (length(local)) local[which.max(diagnostic_y[local])] else
        which.min(abs(diagnostic_x - value))
    }, integer(1))
  } else integer(0)
  detected_height <- diagnostic_y[detected_indices]
  left_flank_penalty <- function(decoded) {
    prediction <- diagnostic_prediction_from(decoded)
    rising_mask <- rep(FALSE, length(diagnostic_x))
    centres <- if (length(detected_peak_x)) detected_peak_x else decoded$miu
    for (centre in centres) {
      rising_mask <- rising_mask |
        (diagnostic_x >= centre - 0.55 * spacing &
           diagnostic_x <= centre - 0.04 * spacing)
    }
    if (sum(rising_mask) < 10L) return(0)
    mean(pmax((prediction[rising_mask] - diagnostic_y[rising_mask]) /
                y_scale, 0)^2)
  }
  relative_left_flank_penalty <- function(decoded) {
    prediction <- diagnostic_prediction_from(decoded)
    rising_mask <- rep(FALSE, length(diagnostic_x))
    centres <- if (length(detected_peak_x)) detected_peak_x else decoded$miu
    for (centre in centres) {
      rising_mask <- rising_mask |
        (diagnostic_x >= centre - 0.55 * spacing &
           diagnostic_x <= centre - 0.04 * spacing)
    }
    if (sum(rising_mask) < 10L) return(0)
    denominator <- pmax(diagnostic_y[rising_mask], 0.03 * y_scale)
    mean(((prediction[rising_mask] - diagnostic_y[rising_mask]) /
            denominator)^2)
  }
  decode <- function(parameters) {
    cursor <- 1L
    logits <- c(parameters[cursor:(cursor + component_count - 2L)], 0)
    cursor <- cursor + component_count - 1L
    logits <- logits - max(logits)
    proportions <- exp(logits) / sum(exp(logits))
    scale_value <- exp(parameters[cursor])
    cursor <- cursor + 1L
    first_value <- parameters[cursor]
    cursor <- cursor + 1L
    gap_ratios <- exp(parameters[cursor:(cursor + component_count - 2L)])
    cursor <- cursor + component_count - 1L
    fitted_miu <- c(first_value, first_value + cumsum(spacing * gap_ratios))
    fitted_sigma <- exp(parameters[cursor:(cursor + component_count - 1L)])
    list(alpha = proportions * scale_value, miu = fitted_miu,
         sigma = fitted_sigma, gap_ratios = gap_ratios)
  }

  prediction_from <- function(decoded) {
    design <- vapply(seq_len(component_count), function(j) {
      stats::dnorm(x, decoded$miu[j], decoded$sigma[j]) *
        total_samples_fit * bin_width
    }, numeric(length(x)))
    as.numeric(design %*% decoded$alpha)
  }

  diagnostic_prediction_from <- function(decoded) {
    design <- vapply(seq_len(component_count), function(j) {
      stats::dnorm(diagnostic_x, decoded$miu[j], decoded$sigma[j]) *
        total_samples_fit * bin_width
    }, numeric(length(diagnostic_x)))
    as.numeric(design %*% decoded$alpha)
  }

  objective <- function(parameters) {
    decoded <- decode(parameters)
    prediction <- diagnostic_prediction_from(decoded)
    if (any(!is.finite(prediction))) return(1e6)
    scaled_error <- (prediction - diagnostic_y) / y_scale
    curve_error <- sum(point_weights * scaled_error^2) / sum(point_weights)
    global_overshoot_error <- sum(
      point_weights * pmax(scaled_error, 0)^2
    ) / sum(point_weights)
    log_error <- mean(((log1p(pmax(prediction, 0)) - log1p(diagnostic_y)) /
                       log_y_scale)^2)
    slope_error <- if (length(prediction) > 1L) {
      mean((diff(prediction - diagnostic_y) / y_scale)^2)
    } else 0
    rising_overshoot <- left_flank_penalty(decoded)
    relative_rising_error <- relative_left_flank_penalty(decoded)
    region_imbalance <- evaluate_harmonic_regions(
      diagnostic_x, diagnostic_y, prediction, spacing
    )$max_nrmse^2

    peak_error <- 0
    peak_underfit_error <- 0
    position_anchor <- 0
    total_peak_slope_error <- 0
    if (length(detected_indices)) {
      height_scale <- pmax(detected_height, 0.05 * y_scale)
      relative_peak_error <- (prediction[detected_indices] -
                                detected_height) / height_scale
      peak_error <- mean(relative_peak_error^2)
      peak_underfit_error <- mean(pmax(-relative_peak_error, 0)^2)
      component_assignment <- pmin(
        component_count,
        pmax(1L, round(detected_peak_x / spacing))
      )
      position_anchor <- mean(
        (decoded$miu[component_assignment] - detected_peak_x)^2 / spacing^2
      )
      slope_terms <- vapply(seq_along(detected_indices), function(i) {
        index <- detected_indices[i]
        if (index <= 1L || index >= length(prediction)) return(0)
        local_slope <- (prediction[index + 1L] - prediction[index - 1L]) /
          (diagnostic_x[index + 1L] - diagnostic_x[index - 1L])
        local_slope * spacing / max(detected_height[i], 0.05 * y_scale)
      }, numeric(1))
      total_peak_slope_error <- mean(slope_terms^2)
    }

    log_gap <- log(decoded$gap_ratios)
    spacing_penalty <- mean(log_gap^2)
    spacing_smoothness <- if (length(log_gap) > 1L) mean(diff(log_gap)^2) else 0
    log_width <- log(decoded$sigma)
    width_smoothness <- if (length(log_width) > 1L) {
      pair_reliability <- pmin(rel[-1L], rel[-component_count])
      mean((0.25 + 0.75 * (1 - pair_reliability)) * diff(log_width)^2)
    } else 0
    broad_threshold <- if (identical(profile, "width_balanced")) 0.40 else 0.52
    broad_weak_penalty <- mean(
      (1 - rel) * pmax(decoded$sigma / spacing - broad_threshold, 0)^2
    )
    width_level_penalty <- 0
    if (identical(profile, "width_balanced")) {
      # When only one width is identifiable (bird-like overlap), the remaining
      # widths borrow strength from that anchor instead of expanding to absorb
      # a neighbouring component.  With several reliable peaks the target is a
      # smooth log-width trend rather than a fixed adjacent ratio.
      reliable_weights <- 0.05 + rel^2
      rank <- log(seq_len(component_count))
      if (component_count > 1L && sum(rel >= 0.25) >= 2L) {
        trend_fit <- stats::lm.wfit(cbind(1, rank), log_width,
                                    w = reliable_weights)
        width_target <- as.numeric(cbind(1, rank) %*%
                                     trend_fit$coefficients)
      } else {
        anchor <- which.max(rel)
        width_target <- rep(log_width[anchor], component_count)
      }
      width_level_penalty <- mean((1 - rel) *
                                    (log_width - width_target)^2)
    }

    weights <- switch(profile,
      standard = c(curve = 0.50, log = 0.12, peak = 0.13, slope = 0.05,
                   anchor = 0.07, rising = 0.12, underfit = 0.16,
                   spacing = 0.055,
                   spacing_smooth = 0.035, width_smooth = 0.035,
                   broad = 0.015, width_level = 0, region = 0,
                   relative_rising = 0, overshoot = 0, total_slope = 0),
      width_balanced = c(curve = 0.45, log = 0.10, peak = 0.12, slope = 0.04,
                         anchor = 0.07, rising = 0.12, underfit = 0.16,
                         spacing = 0.05,
                         spacing_smooth = 0.03, width_smooth = 0.10,
                         broad = 0.05, width_level = 0.12, region = 0,
                         relative_rising = 0, overshoot = 0,
                         total_slope = 0),
      flexible = c(curve = 0.56, log = 0.13, peak = 0.13, slope = 0.06,
                   anchor = 0.045, rising = 0.12, underfit = 0.18,
                   spacing = 0.025,
                   spacing_smooth = 0.025, width_smooth = 0.025,
                   broad = 0.01, width_level = 0, region = 0,
                   relative_rising = 0, overshoot = 0, total_slope = 0),
      flank_balanced = c(curve = 0.52, log = 0.12, peak = 0.14, slope = 0.05,
                         anchor = 0.08, rising = 0.12, underfit = 0.18,
                         spacing = 0.05, spacing_smooth = 0.03,
                         width_smooth = 0.012, broad = 0.02,
                         width_level = 0, region = 0,
                         relative_rising = 0.02, overshoot = 0,
                         total_slope = 0),
      region_balanced = c(curve = 0.52, log = 0.12, peak = 0.14, slope = 0.05,
                          anchor = 0.06, rising = 0.16, underfit = 0.18,
                          spacing = 0.03, spacing_smooth = 0.025,
                          width_smooth = 0.03, broad = 0.01,
                          width_level = 0, region = 0.30,
                          relative_rising = 0, overshoot = 0,
                          total_slope = 0),
      overlap_balanced = c(
        curve = 0.50, log = 0.12, peak = 0.16, slope = 0.05,
        anchor = 0.15, rising = 0.16, underfit = 0.22,
        spacing = 0.03, spacing_smooth = 0.025,
        width_smooth = 0.03, broad = 0.01,
        width_level = 0, region = 0.28,
        relative_rising = 0, overshoot = 0.45, total_slope = 0
      ),
      low_ploidy_balanced = c(
        curve = 0.50, log = 0.12, peak = 0.15, slope = 0.05,
        anchor = 0.18, rising = 0.14, underfit = 0.18,
        spacing = 0.05, spacing_smooth = 0.03,
        width_smooth = 0.025, broad = 0.02,
        width_level = 0, region = 0.15,
        relative_rising = 0.02, overshoot = 0.70, total_slope = 0
      ),
      total_peak_aligned = c(
        curve = 0.70, log = 0.12, peak = 0.30, slope = 0.05,
        anchor = 0.04, rising = 0.12, underfit = 0.30,
        spacing = 0.04, spacing_smooth = 0.03,
        width_smooth = 0.025, broad = 0.015,
        width_level = 0, region = 0.10,
        relative_rising = 0.02, overshoot = 0.10,
        total_slope = 0.030
      )
    )
    weights["curve"] * curve_error + weights["log"] * log_error +
      weights["peak"] * peak_error + weights["slope"] * slope_error +
      weights["underfit"] * peak_underfit_error +
      weights["anchor"] * position_anchor +
      weights["rising"] * rising_overshoot +
      weights["spacing"] * spacing_penalty +
      weights["spacing_smooth"] * spacing_smoothness +
      weights["width_smooth"] * width_smoothness +
      weights["broad"] * broad_weak_penalty +
      weights["width_level"] * width_level_penalty +
      weights["region"] * region_imbalance +
      weights["relative_rising"] * relative_rising_error +
      weights["overshoot"] * global_overshoot_error +
      weights["total_slope"] * total_peak_slope_error
  }

  initial <- c(
    alpha_parameters, log(initial_scale), initial_first,
    log(initial_gap_ratio), initial_log_sigma
  )
  alpha_lower <- rep(-12, component_count - 1L)
  alpha_upper <- rep(12, component_count - 1L)
  scale_lower <- log(max(initial_scale * 0.35, 1e-8))
  scale_upper <- log(max(initial_scale * 3.0, 1e-7))
  sigma_lower <- if (identical(profile, "total_peak_aligned")) {
    pmax(minimum_sigma, 0.995 * sigma)
  } else rep(minimum_sigma, component_count)
  sigma_upper <- if (identical(profile, "total_peak_aligned")) {
    pmin(maximum_sigma, 1.005 * sigma)
  } else rep(maximum_sigma, component_count)
  lower <- c(alpha_lower, scale_lower, first_lower * spacing,
             rep(log(gap_lower), component_count - 1L),
             log(sigma_lower))
  upper <- c(alpha_upper, scale_upper, first_upper * spacing,
             rep(log(gap_upper), component_count - 1L),
             log(sigma_upper))
  objective_before <- objective(initial)
  starts <- list(initial)
  # Deterministic alternative starts reduce local-optimum dependence without
  # making repeated automatic runs non-reproducible.
  if (!identical(profile, "total_peak_aligned")) {
    common_width <- stats::weighted.mean(initial_log_sigma, 0.05 + rel^2)
    common_start <- initial
    common_start[(length(common_start) - component_count + 1L):
                   length(common_start)] <- common_width
    starts[[length(starts) + 1L]] <- common_start
  }
  if (profile %in%
      c("flexible", "region_balanced", "overlap_balanced",
        "low_ploidy_balanced", "total_peak_aligned") &&
      component_count > 2L) {
    wave_start <- initial
    gap_start <- (component_count + 2L):(2L * component_count)
    wave_start[gap_start] <- pmin(log(gap_upper), pmax(log(gap_lower),
      wave_start[gap_start] + 0.04 * sin(seq_along(gap_start))))
    starts[[length(starts) + 1L]] <- wave_start
  }
  fits <- lapply(starts, function(start_value) tryCatch(
    stats::optim(
      start_value, objective, method = "L-BFGS-B", lower = lower,
      upper = upper,
      control = list(maxit = as.integer(max_iterations), factr = 1e7,
                     pgtol = 1e-7)
    ),
    error = function(e) NULL
  ))
  valid_fits <- vapply(fits, function(item) {
    !is.null(item) && is.finite(item$value)
  }, logical(1))
  refined <- if (any(valid_fits)) {
    fits[[which(valid_fits)[which.min(vapply(
      fits[valid_fits], function(item) item$value, numeric(1)
    ))]]]
  } else NULL
  if (is.null(refined) || !is.finite(refined$value) ||
      refined$value >= objective_before * 0.999) {
    return(list(alpha = alpha, miu = miu, sigma = sigma,
                refined = FALSE, objective_before = objective_before,
                objective_after = if (is.null(refined)) Inf else refined$value))
  }
  decoded <- decode(refined$par)
  list(alpha = decoded$alpha, miu = decoded$miu,
       sigma = decoded$sigma, refined = TRUE,
       profile = profile,
       objective_before = objective_before, objective_after = refined$value,
       convergence = refined$convergence)
}

# Fixed-centre width compromise.  This candidate is aimed at two/three highly
# overlapping components where a joint optimiser can move too many parameters
# at once.  It keeps the accepted peak positions, re-estimates component areas,
# and applies only a tunable soft penalty to adjacent log-width differences.
refine_gmm_width_compromise <- function(x, y, alpha, miu, sigma,
                                        total_samples_fit, first_peak_x,
                                        detected_peak_x = numeric(0),
                                        balance_strength = 0.02,
                                        max_iterations = 500L) {
  valid <- is.finite(x) & is.finite(y) & y >= 0
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  component_count <- length(alpha)
  spacing <- as.numeric(first_peak_x)
  if (length(x) < 8L || component_count < 2L ||
      !is.finite(spacing) || spacing <= 0) {
    return(list(alpha = alpha, miu = miu, sigma = sigma, refined = FALSE))
  }
  ordered <- order(miu)
  miu <- as.numeric(miu[ordered])
  sigma <- as.numeric(sigma[ordered])
  alpha <- pmax(as.numeric(alpha[ordered]), 1e-10)
  scale_start <- sum(alpha)
  proportions <- alpha / scale_start
  logits_start <- log(proportions[-component_count] /
                        proportions[component_count])
  sigma_minimum <- 0.06 * spacing
  sigma_maximum <- 0.70 * spacing
  sigma <- pmin(sigma_maximum, pmax(sigma_minimum, sigma))

  x_steps <- diff(sort(unique(x)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  y_scale <- max(y, na.rm = TRUE)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1
  log_y_scale <- max(log1p(y_scale), 1)
  point_weights <- 0.20 + 0.80 * sqrt(pmax(y, 0) / y_scale)
  y_smooth <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  y_smooth[!is.finite(y_smooth)] <- y[!is.finite(y_smooth)]
  detected_peak_x <- as.numeric(detected_peak_x)
  detected_peak_x <- detected_peak_x[
    is.finite(detected_peak_x) & detected_peak_x >= min(x) &
      detected_peak_x <= max(x)
  ]
  detected_indices <- vapply(detected_peak_x, function(value) {
    local <- which(x >= value - 0.18 * spacing &
                   x <= value + 0.18 * spacing)
    if (length(local)) local[which.max(y[local])] else
      which.min(abs(x - value))
  }, integer(1))
  left_flank_penalty <- function(prediction) {
    rising_mask <- rep(FALSE, length(x))
    for (centre in miu) {
      rising_mask <- rising_mask |
        (x >= centre - 0.45 * spacing & x <= centre - 0.08 * spacing)
    }
    if (sum(rising_mask) < 3L) return(0)
    mean(pmax((prediction[rising_mask] - y_smooth[rising_mask]) /
                y_scale, 0)^2)
  }
  decode <- function(parameters) {
    logits <- c(parameters[seq_len(component_count - 1L)], 0)
    logits <- logits - max(logits)
    fitted_proportions <- exp(logits) / sum(exp(logits))
    scale_value <- exp(parameters[component_count])
    fitted_sigma <- exp(parameters[(component_count + 1L):
                                     (2L * component_count)])
    list(alpha = fitted_proportions * scale_value, sigma = fitted_sigma)
  }
  objective <- function(parameters) {
    decoded <- decode(parameters)
    design <- vapply(seq_len(component_count), function(j) {
      stats::dnorm(x, miu[j], decoded$sigma[j]) *
        total_samples_fit * bin_width
    }, numeric(length(x)))
    prediction <- as.numeric(design %*% decoded$alpha)
    scaled_error <- (prediction - y) / y_scale
    curve_error <- sum(point_weights * scaled_error^2) / sum(point_weights)
    log_error <- mean(((log1p(pmax(prediction, 0)) - log1p(y)) /
                       log_y_scale)^2)
    relative_peak_error <- if (length(detected_indices)) {
      peak_scale <- pmax(y[detected_indices], 0.05 * y_scale)
      (prediction[detected_indices] - y[detected_indices]) /
        peak_scale
    } else numeric(0)
    peak_error <- if (length(relative_peak_error)) {
      mean(relative_peak_error^2)
    } else 0
    peak_underfit_error <- if (length(relative_peak_error)) {
      mean(pmax(-relative_peak_error, 0)^2)
    } else 0
    log_width <- log(decoded$sigma)
    width_difference <- mean(diff(log_width)^2)
    broad_penalty <- mean(pmax(decoded$sigma / spacing - 0.52, 0)^2)
    0.62 * curve_error + 0.13 * log_error + 0.20 * peak_error +
      0.16 * peak_underfit_error +
      0.12 * left_flank_penalty(prediction) +
      balance_strength * width_difference + 0.025 * broad_penalty
  }
  initial <- c(logits_start, log(scale_start), log(sigma))
  lower <- c(rep(-12, component_count - 1L),
             log(max(scale_start * 0.35, 1e-8)),
             rep(log(sigma_minimum), component_count))
  upper <- c(rep(12, component_count - 1L),
             log(max(scale_start * 3, 1e-7)),
             rep(log(sigma_maximum), component_count))
  objective_before <- objective(initial)
  refined <- tryCatch(stats::optim(
    initial, objective, method = "L-BFGS-B", lower = lower, upper = upper,
    control = list(maxit = as.integer(max_iterations), factr = 1e7,
                   pgtol = 1e-7)
  ), error = function(e) NULL)
  if (is.null(refined) || !is.finite(refined$value) ||
      refined$value >= objective_before * 0.9995) {
    return(list(alpha = alpha, miu = miu, sigma = sigma, refined = FALSE,
                objective_before = objective_before,
                objective_after = if (is.null(refined)) Inf else refined$value))
  }
  decoded <- decode(refined$par)
  list(alpha = decoded$alpha, miu = miu, sigma = decoded$sigma,
       refined = TRUE, balance_strength = balance_strength,
       objective_before = objective_before, objective_after = refined$value,
       convergence = refined$convergence)
}

#' EM algorithm with post-fitting main peak calibration
EM_algorithm_with_calibration <- function(x, y, weights, n, first_main_peak_x, first_peak_x, 
                                          peak_x_focused, main_peak_x_symmetric, main_peak_y_symmetric,
                                          max_sigma_ratio = 1.35,
                                          sigma_mode = "legacy",
                                          position_mode = "legacy",
                                          observed_peak_count = NULL,
                                          max_iter = NULL, tolerance = 1e-8,
                                          extra_component_mode = "legacy") {
  
  # Step 1: Perform normal EM fitting
  EM_result <- EM_algorithm(x, y, weights, n, first_main_peak_x, first_peak_x, peak_x_focused,
                            max_sigma_ratio = max_sigma_ratio,
                            sigma_mode = sigma_mode,
                            position_mode = position_mode,
                            observed_peak_count = observed_peak_count,
                            max_iter = max_iter, tolerance = tolerance,
                            extra_component_mode = extra_component_mode)
  
  alpha <- EM_result$alpha
  miu <- EM_result$miu
  sigma <- EM_result$sigma
  main_peak_component <- EM_result$main_peak_component
  
  # Step 2: Calibrate main peak component with original main peak data
  calibrated_result <- calibrate_main_peak_component(
    alpha, miu, sigma, main_peak_component,
    main_peak_x_symmetric, main_peak_y_symmetric,
    max_sigma_ratio = max_sigma_ratio,
    sigma_mode = sigma_mode,
    first_peak_x = first_peak_x,
    sigma_reliability = EM_result$sigma_reliability
  )
  
  # Return calibrated results
  return(list(
    alpha = calibrated_result$alpha,
    miu = calibrated_result$miu,
    sigma = calibrated_result$sigma,
    converged = EM_result$converged,
    iterations = EM_result$iterations,
    max_iter = EM_result$max_iter,
    main_peak_component = main_peak_component,
    fit_extra_component = EM_result$fit_extra_component,
    internal_n = EM_result$internal_n,
    original_n = EM_result$original_n,
    sigma_reliability = EM_result$sigma_reliability
  ))
}

#' Calibrate main peak component using original peak data
calibrate_main_peak_component <- function(alpha, miu, sigma, main_peak_component,
                                          main_peak_x_symmetric, main_peak_y_symmetric,
                                          max_sigma_ratio = 1.35,
                                          sigma_mode = "legacy",
                                          first_peak_x = NULL,
                                          sigma_reliability = NULL) {
  
  cat("=== Main peak calibration started ===\n")
  cat("Calibrating main peak component:", main_peak_component, "\n")
  cat("Original main peak position:", miu[main_peak_component], "\n")
  cat("Original main peak standard deviation:", sigma[main_peak_component], "\n")
  
  # Extract main peak data
  peak_x <- main_peak_x_symmetric
  peak_y <- main_peak_y_symmetric
  
  if (length(peak_x) < 3) {
    cat("Insufficient main peak data points, skipping calibration\n")
    return(list(alpha = alpha, miu = miu, sigma = sigma))
  }
  
  # Fit normal distribution to main peak data
  fit_main_peak <- fit_normal_to_peak(peak_x, peak_y)
  
  cat("Calibrated main peak position:", fit_main_peak$mean, "\n")
  cat("Calibrated main peak standard deviation:", fit_main_peak$sd, "\n")
  
  # Update parameters of main peak component
  miu_calibrated <- miu
  sigma_calibrated <- sigma
  alpha_calibrated <- alpha
  
  # Update main peak component with fitting results
  miu_calibrated[main_peak_component] <- fit_main_peak$mean
  sigma_calibrated[main_peak_component] <- fit_main_peak$sd

  # Calibration must obey the same width model used by EM. The manual interface
  # keeps its historical ratio projection; automatic fitting uses the evidence-
  # weighted common/trend/adaptive model instead.
  if (identical(sigma_mode, "legacy")) {
    sorted_components <- order(miu_calibrated)
    for (projection_step in seq_len(5L)) {
      for (i in 2:length(sorted_components)) {
        previous <- sorted_components[i - 1L]
        current <- sorted_components[i]
        sigma_calibrated[current] <- min(
          max(sigma_calibrated[current],
              sigma_calibrated[previous] / max_sigma_ratio),
          sigma_calibrated[previous] * max_sigma_ratio
        )
      }
      if (length(sorted_components) > 1L) {
        for (i in (length(sorted_components) - 1L):1L) {
          current <- sorted_components[i]
          following <- sorted_components[i + 1L]
          sigma_calibrated[current] <- min(
            max(sigma_calibrated[current],
                sigma_calibrated[following] / max_sigma_ratio),
            sigma_calibrated[following] * max_sigma_ratio
          )
        }
      }
    }
  } else {
    sigma_calibrated <- regularize_sigma_widths(
      sigma_calibrated, miu_calibrated, first_peak_x,
      sigma_reliability, sigma_mode, max_sigma_ratio
    )
  }
  
  # Adjust mixing coefficients to maintain total probability of 1
  total_alpha <- sum(alpha_calibrated)
  alpha_calibrated <- alpha_calibrated / total_alpha
  
  cat("=== Main peak calibration completed ===\n")
  
  return(list(
    alpha = alpha_calibrated,
    miu = miu_calibrated,
    sigma = sigma_calibrated
  ))
}

#' Fit normal distribution to peak data using maximum likelihood
fit_normal_to_peak <- function(x, y) {
  
  # Use weighted maximum likelihood estimation
  total_count <- sum(y)
  if (total_count == 0) {
    return(list(mean = mean(x), sd = sd(x)))
  }
  
  # Weighted mean
  weighted_mean <- sum(x * y) / total_count
  
  # Weighted standard deviation
  weighted_variance <- sum((x - weighted_mean)^2 * y) / total_count
  weighted_sd <- sqrt(weighted_variance)
  
  # Ensure reasonable standard deviation
  if (weighted_sd <= 0 || is.na(weighted_sd)) {
    weighted_sd <- sd(x) * 0.5
  }
  
  # Limit standard deviation range to avoid overfitting
  min_sd <- (max(x) - min(x)) * 0.05
  max_sd <- (max(x) - min(x)) * 0.5
  weighted_sd <- max(min_sd, min(weighted_sd, max_sd))
  
  return(list(mean = weighted_mean, sd = weighted_sd))
}

# Apply one of the peak-centre priors used by the EM M-step.  The adaptive
# modes deliberately distinguish a visible peak from a biologically expected
# but unresolved peak: the former is anchored to the observed maximum, while
# the latter keeps a harmonic target with a non-zero, evidence-dependent
# allowance for displacement.
constrain_peak_center <- function(new_miu, component, main_peak_component,
                                  first_main_peak_x, first_peak_x,
                                  observed_target = NA_real_,
                                  explicit_peak_set = FALSE,
                                  position_mode = "legacy") {
  harmonic_target <- first_peak_x * component

  if (identical(position_mode, "legacy")) {
    if (component == main_peak_component) {
      target <- first_main_peak_x
      limit <- first_peak_x * 0.10
      return(min(max(new_miu, target - limit), target + limit))
    }
    limit <- first_peak_x * 0.05
    if (abs(new_miu - harmonic_target) <= limit) return(new_miu)
    adjusted <- 0.20 * new_miu + 0.80 * harmonic_target
    return(min(max(adjusted, harmonic_target - limit),
               harmonic_target + limit))
  }

  is_relaxed <- identical(position_mode, "relaxed")
  has_observed_target <- is.finite(observed_target)

  if (component == main_peak_component) {
    # The main peak is itself directly observed.  Keep that strong anchor, but
    # avoid the old hard clipping at exactly ten per cent of one peak spacing.
    data_strength <- if (is_relaxed) 0.65 else 0.50
    target <- first_main_peak_x
    candidate <- (1 - data_strength) * new_miu + data_strength * target
    limit <- first_peak_x * if (is_relaxed) 0.30 else 0.20
    return(min(max(candidate, target - limit), target + limit))
  }

  if (explicit_peak_set && has_observed_target) {
    # n explicit maxima are available: the histogram determines the centre and
    # the harmonic relationship is only a weak biological prior.  In R, for
    # example, this lets the fourth centre approach 411 instead of being pulled
    # back to the exact harmonic value 392.
    data_strength <- if (is_relaxed) 0.75 else 0.58
    harmonic_strength <- if (is_relaxed) 0 else 0.03
    candidate <- (1 - data_strength) * new_miu +
      data_strength * observed_target
    candidate <- (1 - harmonic_strength) * candidate +
      harmonic_strength * harmonic_target
    limit <- first_peak_x * if (is_relaxed) 0.35 else 0.25
    return(min(max(candidate, observed_target - limit),
               observed_target + limit))
  }

  # An unresolved component remains harmonic-led, but never exactly locked.
  # A detected shoulder/local maximum expands the allowance from 10% to 15%
  # (18%/22% in the deliberately loose comparison mode).
  if (has_observed_target) {
    evidence_strength <- if (is_relaxed) 0.38 else 0.24
    harmonic_strength <- if (is_relaxed) 0.14 else 0.28
    limit_fraction <- if (is_relaxed) 0.22 else 0.15
    candidate <- (1 - evidence_strength) * new_miu +
      evidence_strength * observed_target
  } else {
    harmonic_strength <- if (is_relaxed) 0.20 else 0.38
    limit_fraction <- if (is_relaxed) 0.18 else 0.10
    candidate <- new_miu
  }
  candidate <- (1 - harmonic_strength) * candidate +
    harmonic_strength * harmonic_target
  limit <- first_peak_x * limit_fraction
  min(max(candidate, harmonic_target - limit), harmonic_target + limit)
}

#' EM algorithm with dynamic main peak assignment and peak distance constraints
EM_algorithm <- function(x, y, weights, n, first_main_peak_x, first_peak_x, peak_x_focused,
                         max_sigma_ratio = 1.35, max_iter = NULL,
                         tolerance = 1e-8,
                         extra_component_mode = "legacy",
                         sigma_mode = "legacy",
                         position_mode = "legacy",
                         observed_peak_count = NULL) {
  if (is.null(observed_peak_count) || !is.finite(observed_peak_count)) {
    observed_peak_count <- length(peak_x_focused)
  }
  observed_peak_count <- as.integer(observed_peak_count)
  # ================== New Logic ==================
  # Preserve the original rule: add n+1 only when fewer than n reliable peaks
  # were detected; data with n or more explicit peaks must fit exactly n.
  fit_one_extra <- switch(
    extra_component_mode,
    always = TRUE,
    never = FALSE,
    legacy = observed_peak_count < n,
    stop("Unknown extra_component_mode: ", extra_component_mode)
  )
  if (fit_one_extra) {
    internal_n <- n + 1  # Fit more peaks internally
    cat(sprintf("Fitting optional extra component (%d components internally)\n",
                internal_n))
    fit_extra_component <- TRUE  # Mark if extra components are fitted
  } else {
    internal_n <- n
    fit_extra_component <- FALSE  # No extra components fitted
  }
  # =============================================
  
  # Improved main peak assignment logic
  alpha <- rep(1/internal_n, internal_n)
  
  # Smarter main peak component assignment
  main_peak_component <- 1  # Default value
  
  # Calculate which component should represent the main peak (based on multiple relationship)
  expected_component <- round(first_main_peak_x / first_peak_x)
  expected_component <- min(max(expected_component, 1), internal_n)  # Limit to internal_n range
  
  cat("Expected main peak component based on ratio:", expected_component, "\n")
  
  if (length(peak_x_focused) > 0) {
    # Find the detected peak closest to the main peak
    peak_distances <- abs(peak_x_focused - first_main_peak_x)
    closest_peak_idx <- which.min(peak_distances)
    closest_peak <- peak_x_focused[closest_peak_idx]
    
    cat("Closest detected peak to main peak:", closest_peak, "at index:", closest_peak_idx, "\n")
    
    # Check if detected peak position matches expected position
    if (abs(closest_peak - first_main_peak_x) < first_peak_x * 0.3) {
      # If detected peak is close to main peak, verify if its position matches expected multiple relationship
      detected_ratio <- round(closest_peak / first_peak_x)
      if (detected_ratio == expected_component) {
        # Detected peak position matches expected multiple relationship, use expected component
        main_peak_component <- expected_component
        cat("Detected peak matches expected ratio, using component:", main_peak_component, "\n")
      } else {
        # If inconsistent, use assignment based on detected peaks
        main_peak_component <- expected_component  # Still prioritize expected component
        cat("Detected peak doesn't match expected ratio, but using expected component:", main_peak_component, "\n")
      }
    } else {
      # Large difference between detected peak and main peak position, use assignment based on multiple relationship
      main_peak_component <- expected_component
      cat("Using expected component due to large position difference:", main_peak_component, "\n")
    }
  } else {
    # When no peaks detected, use assignment based on multiple relationship
    main_peak_component <- expected_component
    cat("No peaks detected, using expected component:", main_peak_component, "\n")
  }
  
  cat("Final main peak component assignment:", main_peak_component, "\n")
  
  # Initialize means based on detected peaks count
  miu <- numeric(internal_n)
  selected_peaks <- numeric(0)
  position_targets <- rep(NA_real_, internal_n)
  explicit_peak_set <- observed_peak_count >= n
  
  # Determine initialization strategy based on detected peak count
  if (observed_peak_count < n) {
    # STRICT INITIALIZATION: Each peak at exact multiples of first_peak_x
    cat(sprintf("Using strict distance-based initialization (detected peaks < n, fitting %d peaks)\n", internal_n))
    for (j in 1:internal_n) {
      if (j == main_peak_component) {
        miu[j] <- first_main_peak_x
      } else {
        miu[j] <- first_peak_x * j
      }
    }
  } else {
    # FLEXIBLE INITIALIZATION: Use actual detected peak positions
    cat("Using flexible initialization with actual peak positions (detected peaks >= n)\n")
    
    # Sort detected peaks by position
    sorted_peaks <- sort(peak_x_focused)
    
    # If detected peaks > internal_n, select top internal_n most significant peaks
    if (length(sorted_peaks) > internal_n) {
      # Find height of each peak
      peak_heights <- numeric(length(sorted_peaks))
      for (i in 1:length(sorted_peaks)) {
        peak_idx <- which.min(abs(x - sorted_peaks[i]))
        peak_heights[i] <- y[peak_idx]
      }
      
      # Sort by height, select top internal_n peaks
      height_order <- order(peak_heights, decreasing = TRUE)
      selected_peaks <- sorted_peaks[height_order[1:internal_n]]
      selected_peaks <- sort(selected_peaks)  # Sort again by position
    } else {
      selected_peaks <- sorted_peaks
    }
    
    # Assign selected peak positions to miu
    for (j in 1:internal_n) {
      if (j <= length(selected_peaks)) {
        miu[j] <- selected_peaks[j]
      } else {
        # If selected peaks insufficient, supplement with multiples
        miu[j] <- first_peak_x * j
      }
    }
    
    # Ensure main peak component uses correct position
    if (main_peak_component <= length(selected_peaks)) {
      miu[main_peak_component] <- selected_peaks[main_peak_component]
    }
    position_targets[seq_along(selected_peaks)] <- selected_peaks
  }

  # When fewer than n peaks are visible, associate any usable local maxima or
  # shoulders with their nearest harmonic component.  Missing components keep
  # NA and are therefore governed mainly by the harmonic prior.
  if (!explicit_peak_set && length(peak_x_focused) > 0L) {
    for (observed_peak in peak_x_focused) {
      component <- round(observed_peak / first_peak_x)
      component <- min(max(component, 1L), internal_n)
      old_target <- position_targets[component]
      if (!is.finite(old_target) ||
          abs(observed_peak - first_peak_x * component) <
            abs(old_target - first_peak_x * component)) {
        position_targets[component] <- observed_peak
      }
    }
  }
  
  cat("Initialized means:", round(miu, 2), "\n")

  sigma_reliability <- estimate_sigma_reliability(
    x, y, first_peak_x, internal_n
  )
  cat("Peak-width data reliability:",
      paste(round(sigma_reliability, 3), collapse = " "), "\n")
  
  # Initialize sigma with reasonable values
  sigma <- rep(sd(x), internal_n)
  
  # EM algorithm parameters
  if (is.null(max_iter)) {
    max_iter <- ifelse(internal_n == 2, 50000, 100000)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 ||
      !is.finite(max_iter) || max_iter < 1) {
    stop("'max_iter' must be NULL or a positive number.")
  }
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tolerance) || length(tolerance) != 1 ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("'tolerance' must be a positive number.")
  }
  threshold <- tolerance
  converged <- FALSE
  
  # Safe handling of y values
  y_safe <- y
  y_safe[y_safe < 0] <- 0
  
  # Weighted y values
  y_weighted <- y_safe * weights
  
  # EM algorithm main loop
  prob <- matrix(0, nrow = length(x), ncol = internal_n)
  weight_matrix <- matrix(0, nrow = length(x), ncol = internal_n)
  
  for (step in 1:max_iter) {
    # E-step
    for (j in 1:internal_n) {
      prob_j <- dnorm(x, miu[j], sigma[j])
      prob_j[is.na(prob_j)] <- 1e-9
      weight_matrix[, j] <- alpha[j] * prob_j
    }
    
    row_weight <- rowSums(weight_matrix)
    row_weight[row_weight <= 0] <- 1e-9
    prob <- weight_matrix / row_weight
    
    old_alpha <- alpha
    old_miu <- miu
    old_sigma <- sigma
    
    # M-step with conditional peak distance constraints
    for (j in 1:internal_n) {
      resp_weights <- y_weighted * prob[, j]
      total_weight <- sum(resp_weights)
      
      if (total_weight < 1e-9) {
        # Reinitialize failed component
        if (observed_peak_count < n) {
          miu[j] <- first_peak_x * j  # Strict position
        } else {
          # Use detected peak positions to reinitialize
          if (j <= length(peak_x_focused)) {
            miu[j] <- peak_x_focused[j]
          } else {
            miu[j] <- first_peak_x * j
          }
        }
        sigma[j] <- if (identical(sigma_mode, "legacy")) {
          sd(x) * runif(1, 0.5, 2)
        } else {
          first_peak_x * runif(1, 0.15, 0.35)
        }
        next
      }
      
      # Update mean
      new_miu <- sum(x * resp_weights) / total_weight
      
      # Explicit maxima and unresolved harmonic components use different soft
      # priors in automatic fitting.  "legacy" reproduces the old 5% rule for
      # the comparison PDFs and for the unchanged manual interface.
      miu[j] <- constrain_peak_center(
        new_miu = new_miu,
        component = j,
        main_peak_component = main_peak_component,
        first_main_peak_x = first_main_peak_x,
        first_peak_x = first_peak_x,
        observed_target = position_targets[j],
        explicit_peak_set = explicit_peak_set,
        position_mode = position_mode
      )
      
      # Update standard deviation with width constraints
      new_sigma <- sqrt(sum((x - miu[j])^2 * resp_weights) / total_weight)
      
      # The manual interface retains its historical global-data bounds. Auto
      # fitting uses broad spacing-relative guards so histogram length does not
      # silently determine every peak width.
      if (identical(sigma_mode, "legacy")) {
        new_sigma <- max(new_sigma, sd(x) * 0.1)
        new_sigma <- min(new_sigma, sd(x) * 3)
      } else {
        new_sigma <- max(new_sigma, first_peak_x * 0.015)
        new_sigma <- min(new_sigma, first_peak_x * 0.80)
      }
      
      # Apply progressive width constraint: each peak's width can be at most max_sigma_ratio times the previous one
      if (identical(sigma_mode, "legacy") && j > 1) {
        # Find the previous component in sorted order
        sorted_indices <- order(miu)
        current_rank <- which(sorted_indices == j)
        
        if (current_rank > 1) {
          prev_component <- sorted_indices[current_rank - 1]
          previous_sigma <- sigma[prev_component]
          max_allowed_sigma <- previous_sigma * max_sigma_ratio
          min_allowed_sigma <- previous_sigma / max_sigma_ratio
          
          if (new_sigma > max_allowed_sigma) {
            #cat(sprintf("Component %d sigma constrained: %.4f -> %.4f (prev component %d sigma: %.4f, ratio: %.2f)\n", 
            #            j, new_sigma, max_allowed_sigma, prev_component, sigma[prev_component], max_sigma_ratio))
            new_sigma <- max_allowed_sigma
          }
          if (new_sigma < min_allowed_sigma) {
            new_sigma <- min_allowed_sigma
          }
        }
      }
      
      sigma[j] <- new_sigma
    }
    
    # Preserve the original one-sided height cap only for the manual legacy
    # interface. In automatic mode it caused a systematic downward bias because
    # over-height peaks were reduced while under-height peaks were never raised.
    if (identical(sigma_mode, "legacy")) for (j in 1:internal_n) {
      # Calculate theoretical height of current component at peak position
      gaussian_peak_height <- alpha[j] * dnorm(miu[j], miu[j], sigma[j])
      
      # Find actual height near peak position in original data
      peak_region_indices <- which(x >= miu[j] - sigma[j] & x <= miu[j] + sigma[j])
      if (length(peak_region_indices) > 0) {
        # Use max value in peak region as reference height
        actual_peak_height <- max(y[peak_region_indices] / sum(y * weights), na.rm = TRUE)
        
        # If theoretical height exceeds actual height, adjust alpha
        if (gaussian_peak_height > actual_peak_height * 1.002) {  # Allow 2% tolerance
          # Scale down alpha proportionally to ensure theoretical height doesn't exceed actual
          scale_factor <- actual_peak_height / gaussian_peak_height * 0.998  # Extra 2% safety margin
          old_alpha_j <- alpha[j]
          alpha[j] <- alpha[j] * scale_factor
        }
      }
    }
    
    # Ensure alpha sums to 1
    alpha <- alpha / sum(alpha)
    
    # Preserve the historical second harmonic pull only in the legacy baseline.
    if (identical(position_mode, "legacy")) {
      # Additional constraint: ensure peaks maintain expected spacing
      # Sort components by their means
      sorted_indices <- order(miu)
      miu_sorted <- miu[sorted_indices]
      sigma_sorted <- sigma[sorted_indices]
      
      # ENFORCE EXACT SPACING: Each peak should be at first_peak_x intervals
      for (i in 1:length(miu_sorted)) {
        expected_pos <- first_peak_x * i
        current_pos <- miu_sorted[i]
        
        # If position deviates more than 8%, adjust toward expected position
        if (abs(current_pos - expected_pos) > first_peak_x * 0.05) {
          adjustment_strength <- 0.5  # How strongly to adjust toward expected position
          adjusted_pos <- current_pos * (1 - adjustment_strength) + expected_pos * adjustment_strength
          miu_sorted[i] <- adjusted_pos
          #cat(sprintf("Strict spacing adjustment: component %d from %.2f to %.2f (expected: %.2f)\n",
          #            sorted_indices[i], current_pos, adjusted_pos, expected_pos))
        }
      }
      
      # Restore original order
      miu[sorted_indices] <- miu_sorted
      sigma[sorted_indices] <- sigma_sorted
    } else {
      # Adaptive modes retain ordering and a modest minimum separation without
      # pulling every component back to an exact integer multiple.
      sorted_indices <- order(miu)
      miu_sorted <- miu[sorted_indices]
      minimum_gap <- first_peak_x * 0.20
      if (length(miu_sorted) > 1L) {
        for (i in 2:length(miu_sorted)) {
          miu_sorted[i] <- max(miu_sorted[i], miu_sorted[i - 1L] + minimum_gap)
        }
      }
      miu[sorted_indices] <- miu_sorted
    }
    
    # Width regularisation is continuous rather than a fixed adjacent ratio in
    # automatic mode. Clear A/B/C-like peaks keep their own widths; overlapping
    # bird peaks and weak late sturgeon peaks shrink toward a learned trend.
    sigma <- regularize_sigma_widths(
      sigma, miu, first_peak_x, sigma_reliability,
      sigma_mode = sigma_mode, max_sigma_ratio = max_sigma_ratio
    )
    
    # Update mixing coefficients (recalculate as they might have been modified by peak height constraint)
    for (j in 1:internal_n) {
      resp_weights <- y_weighted * prob[, j]
      total_weight <- sum(resp_weights)
      alpha[j] <- total_weight / sum(y_weighted)
    }
    alpha <- alpha / sum(alpha)
    
    # The second historical cap is likewise manual-only. Automatic fitting uses
    # an unbiased all-clear-peak objective in refine_gmm_curve_widths().
    if (identical(sigma_mode, "legacy")) for (j in 1:internal_n) {
      gaussian_peak_height <- alpha[j] * dnorm(miu[j], miu[j], sigma[j])
      peak_region_indices <- which(x >= miu[j] - sigma[j] & x <= miu[j] + sigma[j])
      
      if (length(peak_region_indices) > 0) {
        actual_peak_height <- max(y[peak_region_indices] / sum(y * weights), na.rm = TRUE)
        
        if (gaussian_peak_height > actual_peak_height * 1.001) {  # Stricter 1% tolerance
          scale_factor <- actual_peak_height / gaussian_peak_height * 0.999
          old_alpha_j <- alpha[j]
          alpha[j] <- alpha[j] * scale_factor
        }
      }
    }
    
    # Finally ensure alpha sums to 1
    alpha <- alpha / sum(alpha)
    
    # Convergence check
    delta <- max(c(
      max(abs(alpha - old_alpha)),
      max(abs(miu - old_miu)),
      max(abs(sigma - old_sigma))
    ))
    
    if (delta < threshold) {
      converged <- TRUE
      break
    }
  }
  
  # ================== Logic Modification ==================
  # No longer filter components, return all components and add marker info
  cat(sprintf("\n=== EM Algorithm Completed, Fitted %d Gaussian Components ===\n", internal_n))
  cat("Fitted extra components:", fit_extra_component, "\n")
  cat("Original n value:", n, "\n")
  cat("Internal fitted components count:", internal_n, "\n")
  # ========================================================
  
  # Final assignment information
  cat("\nFinal component assignment:\n")
  sorted_final <- order(miu)
  for (idx in 1:length(sorted_final)) {
    j <- sorted_final[idx]
    role <- if (j == main_peak_component) " [MAIN PEAK]" else 
      if (j == 1) " [FIRST PEAK]" else ""
    
    # Display info based on detected peak count
    if (observed_peak_count < n) {
      expected_pos <- first_peak_x * idx
      distance_diff <- miu[j] - expected_pos
      spacing_info <- sprintf(" (expected: %.2f, diff: %.2f)", expected_pos, distance_diff)
    } else {
      spacing_info <- " (free positioning)"
    }
    
    spacing_accuracy <- if (idx > 1) {
      actual_spacing <- miu[j] - miu[sorted_final[idx-1]]
      expected_spacing <- first_peak_x
      sprintf(" (spacing: %.2f vs expected: %.2f, diff: %.2f%%)", 
              actual_spacing, expected_spacing, 
              abs(actual_spacing - expected_spacing)/expected_spacing * 100)
    } else ""
    
    sigma_ratio <- if (idx > 1) sprintf(" (sigma ratio: %.3f)", sigma[j] / sigma[sorted_final[idx-1]]) else ""
    
    # Check if extra component
    extra_marker <- if (fit_extra_component && idx > n) " [EXTRA - LIGHT GRAY]" else ""
    
    # Final peak height verification
    gaussian_peak_height <- alpha[j] * dnorm(miu[j], miu[j], sigma[j])
    peak_region_indices <- which(x >= miu[j] - sigma[j] & x <= miu[j] + sigma[j])
    actual_peak_height <- if (length(peak_region_indices) > 0) {
      max(y[peak_region_indices] / sum(y * weights), na.rm = TRUE)
    } else NA
    
    height_status <- if (!is.na(actual_peak_height)) {
      sprintf(" [height ratio: %.3f]", gaussian_peak_height / actual_peak_height)
    } else " [height: N/A]"
    
    cat(sprintf("Component %d: mean = %.2f, sigma = %.4f%s%s%s%s%s%s\n", 
                j, miu[j], sigma[j], sigma_ratio, spacing_info, spacing_accuracy, role, extra_marker, height_status))
  }
  
  return(list(
    alpha = alpha,
    miu = miu,
    sigma = sigma,
    converged = converged,
    iterations = step,
    max_iter = max_iter,
    main_peak_component = main_peak_component,
    fit_extra_component = fit_extra_component,  # Add marker
    internal_n = internal_n,                    # Record internal fitted count
    original_n = n,                             # Record original n value
    sigma_reliability = sigma_reliability
  ))
}

#' Calculate genome size from fitted parameters using FIRST peak
calculate_genome_size <- function(x_fit_processed, y_fit_processed, alpha, miu, sigma,
                                  full_data, start_index, right_value, first_peak_component = 1) {
  # Calculate probability density
  total_samples_filtered <- sum(y_fit_processed)
  bin_width <- ifelse(length(unique(diff(x_fit_processed))) == 1, 
                      diff(x_fit_processed)[1], 1)
  
  # Calculate fitted curve
  k <- length(alpha)
  x_vals <- seq(min(x_fit_processed) - bin_width, max(x_fit_processed) + bin_width, 
                length.out = 1000)
  
  # Calculate probability density function
  components_density <- lapply(1:k, function(j) alpha[j] * dnorm(x_vals, miu[j], sigma[j]))
  total_density <- Reduce(`+`, components_density)
  
  # Find FIRST peak position (changed from main peak)
  first_component <- components_density[[first_peak_component]]
  peak_index <- which.max(first_component)
  kmer_depth <- x_vals[peak_index]
  
  cat("FIRST peak position (kmer_depth):", kmer_depth, "\n")
  cat("First peak component:", first_peak_component, "\n")
  cat("All component means:", round(miu, 2), "\n")
  
  # For data within fitting range, use GMM estimated distribution
  fit_range_upper <- max(x_fit_processed)
  natural_x <- seq(min(x_fit_processed), fit_range_upper, by = 1)
  density_natural_x <- approx(x_vals, total_density, xout = natural_x)$y
  
  # Convert probability density to counts
  frequency_natural_x <- density_natural_x * total_samples_filtered * bin_width
  fit_weighted_frequencies <- natural_x * frequency_natural_x
  fit_total_kmers <- sum(fit_weighted_frequencies)
  
  # === Calculate Individual Component Sizes ===
  component_sizes_bp <- numeric(k)
  for (j in 1:k) {
    # 1. Calculate density for specific component j at natural x steps
    comp_density <- alpha[j] * dnorm(natural_x, miu[j], sigma[j])
    
    # 2. Convert to frequency counts (how many kmers have this frequency)
    comp_freq_counts <- comp_density * total_samples_filtered * bin_width
    # 3. Calculate "mass" (Total Kmers = Sum(Frequency * Count))
    comp_mass <- sum(natural_x * comp_freq_counts)
    
    # 4. Calculate Size in BP (Mass / Depth)
    component_sizes_bp[j] <- comp_mass / kmer_depth
  }
  
  # Calculate kmer count for remaining data
  remaining_data <- full_data[full_data[, 1] > fit_range_upper, ]
  if (nrow(remaining_data) > 0) {
    remaining_weighted_frequencies <- as.numeric(remaining_data[, 1]) * as.numeric(remaining_data[, 2])
    remaining_total_kmers <- sum(remaining_weighted_frequencies)
  } else {
    remaining_total_kmers <- 0
  }
  
  total_kmers <- fit_total_kmers + remaining_total_kmers
  genome_size_bp <- total_kmers / kmer_depth
  
  return(list(
    total_kmers = total_kmers,
    kmer_depth = kmer_depth,
    genome_size_bp = genome_size_bp,
    component_sizes_bp = component_sizes_bp # Return individual sizes
  ))
}

#' Draw plot with smooth curve overlay on original histogram
draw_plot <- function(x_original, y_original, alpha, miu, sigma, genome_size_bp, component_sizes_bp,
                      main_peak_x_symmetric, main_peak_y_symmetric,
                      x_fit_processed, y_fit_processed, full_data = NULL, right_value = NULL,
                      fit_extra_component = FALSE, internal_n = NULL, n = NULL,
                      species_name = "", first_peak_x = NULL) {
  
  k <- length(alpha)
  n_components_total <- length(component_sizes_bp)
  has_residual <- n_components_total > k
  
  bin_width <- ifelse(length(unique(diff(x_original))) == 1, diff(x_original)[1], 0.5)
  
  # ================== Use logarithmic decay for smoother size changes ==================
  if (!is.null(n)) {
    # Base size
    base_size <- 0.8
    
    # Logarithmic decay formula
    if (n <= 5) {
      legend_size <- base_size
    } else {
      # n increases from 5 to 20, size decreases from 0.8 to 0.4
      min_size <- 0.45
      max_size <- base_size
      
      # Use logarithmic function for smooth decay
      log_factor <- log10(n - 4) / log10(16)
      legend_size <- max_size - (max_size - min_size) * log_factor
      
      # Ensure size is within reasonable range
      legend_size <- max(min_size, min(max_size, legend_size))
    }
  } else {
    legend_size <- 0.8
  }
  # ================================================================
  
  # Build title
  if (species_name != "") {
    plot_title <- paste0("Genome size estimation of ", species_name)
  } else {
    plot_title <- "Genome size estimation of given species"
  }
  
  # ================== Critical fix: Determine x-axis range ==================
  # For diploids, keep the displayed grey curve inside the range that was
  # actually fitted. Otherwise the unfitted low-count tail looks like a GMM
  # failure even though it was deliberately excluded as background. Higher
  # ploidies retain the established harmonic display range.
  if (!is.null(n) && n <= 2L && !is.null(right_value) &&
      is.finite(right_value)) {
    x_limit <- right_value
  } else if (!is.null(first_peak_x) && !is.null(n)) {
    x_limit <- first_peak_x * (n + 2)
  } else if (!is.null(right_value)) {
    x_limit <- right_value
  } else {
    # If neither, use max value of original data
    x_limit <- max(x_original)
  }
  
  # Use complete data up to truncation point for histogram
  if (!is.null(full_data) && !is.null(right_value)) {
    # Get complete data from 0 to truncation point
    plot_indices <- which(full_data[, 1] <= x_limit)
    plot_x <- full_data[plot_indices, 1]
    plot_y <- full_data[plot_indices, 2]
  } else {
    # Fall back to original logic
    plot_x <- x_original
    plot_y <- y_original
  }
  
  # Calculate counts for GMM fitting
  total_samples_fit <- sum(y_fit_processed)
  
  # X range starts from 0, extends to full plot range
  x_vals_fit <- seq(0, x_limit + bin_width, length.out = 2000)
  
  # Calculate complete GMM distribution (starting from 0)
  components <- lapply(1:k, function(j) {
    # Use complete GMM formula, calculate from 0
    if (sigma[j] > 0) {
      alpha[j] * dnorm(x_vals_fit, miu[j], sigma[j]) * total_samples_fit * bin_width
    } else {
      rep(0, length(x_vals_fit))
    }
  })
  total_counts <- Reduce(`+`, components)
  
  # Calculate main peak height
  main_peak_height <- max(main_peak_y_symmetric, na.rm = TRUE)
  
  # Calculate y-axis range
  y_max <- main_peak_height * 1.6
  
  # Set graphical parameters for consistent font and size
  par(family = "Helvetica", cex = 1.0)
  
  # ================== Modify plot function call ==================
  plot(NULL, 
       xlim = c(0, x_limit),  # Use calculated x-axis length
       ylim = c(0, y_max), 
       main = plot_title,
       xlab = "kmer frequency", 
       ylab = "count"
  )
  # ===================================================
  
  # Histogram fill color - light gray
  histogram_fill <- "#E8E8E8"
  
  # Draw complete original histogram (Draw data only within x_limit range)
  if (length(plot_x) > 0) {
    rect(
      xleft = plot_x - bin_width/2, 
      ybottom = 0, 
      xright = plot_x + bin_width/2, 
      ytop = plot_y, 
      col = histogram_fill, 
      border = NA
    )
  }
  
  # Draw smooth curve of original data - dark gray
  if (length(plot_x) > 1) {
    # Limit within x_limit range
    smooth_x_max <- min(max(plot_x), x_limit)
    smooth_x <- seq(min(plot_x[-1]), smooth_x_max, length.out = 1000)
    
    # Method 1: Use the same spline reference used by fitting diagnostics.
    if (length(plot_x[-1]) > 3) {
      grey_reference <- make_gmm_reference_curve(
        plot_x[-1], plot_y[-1], dense_points = 1000L
      )
      smooth_x <- grey_reference$x
      smooth_y <- grey_reference$y
      # Draw data only within x_limit range
      lines(smooth_x[smooth_x <= x_limit], 
            smooth_y[smooth_x <= x_limit], 
            col = "#AAAAAA", lwd = 2, lty = 1)  # Dark gray
    } else {
      # Method 2: Use kernel smoothing
      smooth_y <- ksmooth(plot_x[-1], plot_y[-1], bandwidth = bin_width * 2, x.points = smooth_x)$y
      # Draw data only within x_limit range
      lines(smooth_x[smooth_x <= x_limit], 
            smooth_y[smooth_x <= x_limit], 
            col = "#AAAAAA", lwd = 2, lty = 1)  # Dark gray
    }
  }
  
  # Draw complete fitted curve (from 0 to x_limit)
  lines(x_vals_fit[x_vals_fit <= x_limit], 
        total_counts[x_vals_fit <= x_limit], 
        col = "#FF0000CC", lwd = 2)
  
  # Define component color scheme (line and fill colors)
  line_colors <- c(
    "#EA0000", "#FFDC6E", "#5ADBFE", "#861FAB", "#00A087", 
    "#3C5488", "#F39B7F", "#00468B", "#00468B", "#7E6148",
    "#42B540", "#922B21", "#B8860B", "#4527A0", "#556B2F",
    "#C2185B"
  )
  
  fill_colors <- c(
    "#FFBC98CC", "#FFF2CBCC", "#D0EBF3CC", "#C591D6CC", "#B2DFDBCC",
    "#C5CAE9CC", "#FFCCBCCC", "#BBDEFBCC", "#BBDEFBCC", "#D7CCC8CC",
    "#C8E6C9CC", "#E6B0AACC", "#F9E79FCC", "#D1C4E9CC", "#DCEDC8CC",
    "#F8BBD0CC"
  )
  
  # Check if light gray is needed for extra component color
  extra_color <- "#8491B4"  # Line color for extra components
  extra_fill <- "#CFD8DCCC"   # Fill color for extra components
  
  # Define extra component style function
  draw_component_with_fill <- function(x, y, line_color, fill_color, lty = 2, lwd = 1.5) {
    # Draw data only within x_limit range
    x_limited <- x[x <= x_limit]
    y_limited <- y[x <= x_limit]
    
    if (length(x_limited) > 0) {
      # Draw filled polygon
      polygon(c(x_limited, rev(x_limited)), c(y_limited, rep(0, length(y_limited))), 
              col = fill_color, border = NA)
      # Draw lines
      lines(x_limited, y_limited, col = line_color, lty = lty, lwd = lwd)
    }
  }
  
  # Draw components, from back to front (ensure correct overlap)
  for (j in k:1) {
    # Check if extra component (peak n+1 onwards)
    if (fit_extra_component && !is.null(n) && !is.null(internal_n) && 
        internal_n > n && j > n) {
      # Extra component: use light gray
      draw_component_with_fill(x_vals_fit, components[[j]], 
                               extra_color, extra_fill, lty = 2, lwd = 1.5)
    } else {
      # Normal component: use specified color scheme
      line_color_idx <- ((j-1) %% length(line_colors)) + 1
      fill_color_idx <- ((j-1) %% length(fill_colors)) + 1
      
      draw_component_with_fill(x_vals_fit, components[[j]], 
                               line_colors[line_color_idx], 
                               fill_colors[fill_color_idx], 
                               lty = 2, lwd = 1.5)
    }
  }
  
  # === Draw light gray vertical dashed lines for each peak position ===
  # Sort by miu (ascending)
  sorted_miu <- sort(miu)
  
  # Light gray dashed line parameters
  peak_line_color <- "#CCCCCC"  # Light gray
  peak_line_lty <- 2  # Dashed line
  peak_line_lwd <- 1.2
  
  # Draw vertical dashed line for each peak (only within x_limit)
  for (i in 1:length(sorted_miu)) {
    x_pos <- sorted_miu[i]
    
    # Check if x_pos is within plot range
    if (x_pos >= 0 && x_pos <= x_limit) {
      # Draw vertical dashed line
      abline(v = x_pos, col = peak_line_color, lty = peak_line_lty, lwd = peak_line_lwd)
    }
  }
  
  if (length(sorted_miu) > 0) {
    # Get position of first peak
    first_peak_x_draw <- sorted_miu[1]
    
    # Find y-value of first peak on fitted curve
    # Find component index for first peak
    first_peak_idx <- which.min(abs(miu - first_peak_x_draw))
    if (first_peak_idx <= k) {
      # Get curve values for that component
      component_curve <- components[[first_peak_idx]]
      # Find x-index closest to first_peak_x_draw
      x_idx <- which.min(abs(x_vals_fit - first_peak_x_draw))
      peak_y_value <- component_curve[x_idx]
    } else {
      # If not found, use default value
      peak_y_value <- y_max * 0.3
    }
    
    # Calculate label y-position (middle between peak y and ymax)
    label_y <- peak_y_value + (y_max - peak_y_value) / 3
    
    # Set label parameters
    label_cex <- 0.75  # Font size
    label_col <- "black"
    label_font <- 1.5
    
    # Calculate label x-position (right of dashed line, slight offset)
    offset <- x_limit * 0.01  # 1% offset
    label_x <- first_peak_x_draw + offset
    
    # Ensure label is within x_limit
    if (label_x <= x_limit) {
      # Draw label (right of dashed line, 1 decimal place, append 'x')
      text(x = label_x, y = label_y,
           labels = sprintf("%.1f×", first_peak_x_draw),
           col = label_col, cex = label_cex, font = label_font,
           adj = c(0, 0.5),  # Left align, vertically centered
           family = "Helvetica")
    }
  }
  
  # Format component labels
  comp_labels <- character(n_components_total)
  for (j in 1:n_components_total) {
    if (j <= k) {
      # Check if extra component
      if (fit_extra_component && !is.null(n) && j > n) {
        comp_labels[j] <- paste0("Component ", j, ": ", 
                                 sprintf("%.1f", component_sizes_bp[j] / 1e6), " Mb")
      } else {
        comp_labels[j] <- paste0("Component ", j, ": ", 
                                 sprintf("%.1f", component_sizes_bp[j] / 1e6), " Mb")
      }
    } else {
      comp_labels[j] <- paste0("Component > ", k, ": ", 
                               sprintf("%.1f", component_sizes_bp[j] / 1e6), " Mb")
    }
  }
  
  # Set legend colors - use line colors
  legend_line_colors <- character(k)
  for (j in 1:k) {
    if (fit_extra_component && !is.null(n) && j > n) {
      legend_line_colors[j] <- extra_color
    } else {
      line_color_idx <- ((j-1) %% length(line_colors)) + 1
      legend_line_colors[j] <- line_colors[line_color_idx]
    }
  }
  
  if (has_residual) {
    legend_line_colors <- c(legend_line_colors, NA)
  }
  
  # Legend items
  legend_items <- c("Original Data", "Original Smooth", "Total Fit", comp_labels)
  legend_cols <- c("#E0E0E0", "#AAAAAA", "#FF0000CC", legend_line_colors)
  legend_ltys <- c(NA, 1, 1, rep(2, k), if(has_residual) NA else NULL)
  legend_lwds <- c(NA, 2, 2, rep(1.5, k), if(has_residual) NA else NULL)
  legend_pchs <- c(15, NA, NA, rep(NA, n_components_total))
  
  # ================== New: Check if legend obstructs curve ==================
  max_attempts <- 20  # Max attempts
  attempt <- 0
  
  while (attempt < max_attempts) {
    attempt <- attempt + 1
    
    # First draw legend and save its position information
    legend_info <- legend("topright",
                          legend = legend_items,
                          col = legend_cols,
                          lty = legend_ltys,
                          lwd = legend_lwds,
                          pch = legend_pchs,
                          bty = "n",
                          cex = legend_size,
                          plot = FALSE
    )
    
    # === Calculate white background box position ===
    padding_x <- strwidth("M", cex = legend_size) * 1
    padding_y <- strheight("M", cex = legend_size) * 0.8
    
    # Bottom-left coordinates of white background
    bg_x1 <- legend_info$rect$left - padding_x
    bg_x2 <- legend_info$rect$left + legend_info$rect$w + padding_x
    
    # Calculate genome size info height
    # Use user input n as ploidy
    if (!is.null(n)) {
      # Calculate haploid genome size
      haploid_size_bp <- genome_size_bp / n
      haploid_size_mb <- haploid_size_bp / 1e6
      
      genome_text <- paste(
        "Full GSE: ", sprintf("%.1f", genome_size_bp / 1e6), " Mb\n",
        "Haploid GSE: ", sprintf("%.1f", haploid_size_mb), " Mb",
        sep = ""
      )
    } else {
      genome_text <- paste(
        "Full GS: ", sprintf("%.1f", genome_size_bp / 1e6), " Mb",
        sep = ""
      )
    }
    
    # Estimate position of genome size info
    genome_info <- legend("topright",
                          legend = genome_text,
                          col = "black",
                          lty = 0,
                          lwd = 0,
                          pch = NA,
                          bty = "n",
                          cex = legend_size,
                          plot = FALSE
    )
    
    # Count lines
    line_count <- length(strsplit(genome_text, "\n")[[1]])
    
    # Height per line
    line_height <- strheight("M", cex = legend_size) * 1.2
    
    # Total text height
    genome_total_height <- line_count * line_height
    
    # Text top position
    genome_top_y <- legend_info$rect$top - legend_info$rect$h - (strheight("M", cex = legend_size) * 0.5)
    
    # Text bottom position
    genome_bottom_y <- genome_top_y - genome_total_height
    
    # White bg bottom = text bottom - padding_y
    bg_y1 <- genome_bottom_y - padding_y
    bg_y2 <- legend_info$rect$top + padding_y
    
    # ================== Critical check: Does white bg box obstruct curve? ==================
    # 1. Calculate white bg box projection on x-axis
    bg_x_range <- c(bg_x1, bg_x2)
    
    # 2. Find original data points within this x range
    x_in_bg_range <- which(plot_x >= bg_x_range[1] & plot_x <= bg_x_range[2])
    
    if (length(x_in_bg_range) > 0) {
      # 3. Get y-values of these points
      y_in_bg_range <- plot_y[x_in_bg_range]
      
      # 4. Calculate max y-value (curve max) in this region
      max_y_in_bg_range <- max(y_in_bg_range, na.rm = TRUE)
      
      # 5. Add 1-character safety margin
      one_char_height <- strheight("M", cex = 1)
      safe_bottom <- max_y_in_bg_range + one_char_height
      
      # 6. Check if white bg bottom is lower than safe position
      if (bg_y1 < safe_bottom) {
        cat(sprintf("Attempt %d: Legend obstructs curve (bg bottom: %.2f, safe pos: %.2f), shrinking legend %.0f%%\n", 
                    attempt, bg_y1, safe_bottom, (legend_size * 0.95 / legend_size - 1) * 100))
        
        # Shrink legend size by 5%
        legend_size <- legend_size * 0.95
        
        # Ensure not shrinking infinitely
        if (legend_size < 0.3) {
          cat("Warning: Legend size shrunk to min 0.3\n")
          break
        }
        
        # Continue to next attempt
        next
      } else {
        cat(sprintf("Attempt %d: Legend position safe (bg bottom: %.2f, safe pos: %.2f)\n", 
                    attempt, bg_y1, safe_bottom))
        break
      }
    } else {
      # If no data points in white bg region, consider safe
      cat(sprintf("Attempt %d: No data points in white bg region, considered safe\n", attempt))
      break
    }
  }
  
  # If obstruction persists after max attempts, use min size
  if (attempt >= max_attempts) {
    cat("Warning: Could not avoid obstruction after", max_attempts, "attempts, using min size\n")
    legend_size <- 0.3
    # Recalculate legend_info
    legend_info <- legend("topright",
                          legend = legend_items,
                          col = legend_cols,
                          lty = legend_ltys,
                          lwd = legend_lwds,
                          pch = legend_pchs,
                          bty = "n",
                          cex = legend_size,
                          plot = FALSE
    )
  }
  
  # ================== Draw white background and legend ==================
  # Recalculate position (using final legend_size)
  padding_x <- strwidth("M", cex = legend_size) * 1
  padding_y <- strheight("M", cex = legend_size) * 0.8
  
  bg_x1 <- legend_info$rect$left - padding_x
  bg_y1 <- genome_bottom_y - padding_y
  bg_x2 <- legend_info$rect$left + legend_info$rect$w + padding_x
  bg_y2 <- legend_info$rect$top + padding_y
  
  # Ensure white background isn't too low
  min_bg_y1 <- par("usr")[3] + (par("usr")[4] - par("usr")[3]) * 0.1
  bg_y1 <- max(bg_y1, min_bg_y1)
  
  # Draw white background rectangle
  rect(bg_x1, bg_y1, bg_x2, bg_y2, 
       col = "white", border = NA)
  
  # Draw legend
  legend("topright",
         legend = legend_items,
         col = legend_cols,
         lty = legend_ltys,
         lwd = legend_lwds,
         pch = legend_pchs,
         bty = "n",
         cex = legend_size
  )
  
  # === Draw genome size text ===
  text(x = max(plot_x, x_limit), y = genome_top_y,
       labels = genome_text,
       col = "black", cex = legend_size, font = 2, adj = c(1, 1),
       family = "Helvetica")
  
  # Redraw black axes to ensure not obstructed by bg or lines
  axis_col <- "black"
  axis_lwd <- 1.5
  
  # Get current plot region boundaries
  plot_region <- par("usr")  # c(x1, x2, y1, y2)
  
  # Draw top border
  abline(h = plot_region[4], col = axis_col, lwd = axis_lwd)
  
  # Draw right border
  abline(v = plot_region[2], col = axis_col, lwd = axis_lwd)
}

#' Score a GMM curve against the observed histogram
#'
#' The score combines global curve error, main-peak error, peak-height error,
#' log-scale error and area error.  A small complexity penalty is included so
#' that automatic n selection does not always prefer more components.
calculate_gmm_fit_score <- function(x_observed, y_observed, alpha, miu, sigma,
                                    total_samples_fit, main_peak_left,
                                    main_peak_right,
                                    sigma_mode = "legacy",
                                    sigma_reliability = NULL,
                                    detected_peak_x = numeric(0),
                                    first_peak_x = NULL) {
  valid <- is.finite(x_observed) & is.finite(y_observed) & y_observed >= 0
  x_observed <- as.numeric(x_observed[valid])
  y_observed <- as.numeric(y_observed[valid])

  if (length(x_observed) < 3 || !is.finite(total_samples_fit) ||
      total_samples_fit <= 0 || any(!is.finite(sigma)) || any(sigma <= 0)) {
    return(list(
      score = Inf, nrmse = Inf, main_peak_nrmse = Inf,
      peak_height_error = Inf, log_rmse = Inf, area_error = Inf,
      peak_height_bias = Inf, peak_height_ratio = NA_real_,
      main_peak_underfit = Inf,
      clear_peak_height_error = Inf, clear_peak_height_bias = Inf,
      clear_peak_underfit = Inf,
      peak_position_error = Inf, flank_overshoot_error = Inf,
      primary_peak_position_error = Inf,
      total_peak_position_error = Inf,
      primary_total_peak_position_error = Inf,
      primary_total_peak_signed_error = Inf,
      left_overshoot_error = Inf, right_overshoot_error = Inf,
      left_relative_error = Inf,
      global_overshoot_nrmse = Inf, global_underfit_nrmse = Inf,
      region_mean_nrmse = Inf, region_max_nrmse = Inf,
      complexity_penalty = Inf, width_identifiability_penalty = Inf
    ))
  }

  x_steps <- diff(sort(unique(x_observed)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps) > 0) stats::median(x_steps) else 1

  component_counts <- vapply(seq_along(alpha), function(j) {
    alpha[j] * stats::dnorm(x_observed, miu[j], sigma[j]) *
      total_samples_fit * bin_width
  }, numeric(length(x_observed)))
  if (is.null(dim(component_counts))) {
    predicted <- component_counts
  } else {
    predicted <- rowSums(component_counts)
  }
  predicted[!is.finite(predicted)] <- 0

  # All visual diagnostics use the continuous grey curve shown in the PDF.
  # The integer-bin prediction is retained only for the histogram area check.
  reference <- make_gmm_reference_curve(x_observed, y_observed)
  diagnostic_x <- reference$x
  diagnostic_y <- pmax(reference$y, 0)
  dense_component_counts <- vapply(seq_along(alpha), function(j) {
    alpha[j] * stats::dnorm(diagnostic_x, miu[j], sigma[j]) *
      total_samples_fit * bin_width
  }, numeric(length(diagnostic_x)))
  diagnostic_prediction <- if (is.null(dim(dense_component_counts))) {
    dense_component_counts
  } else rowSums(dense_component_counts)
  diagnostic_prediction[!is.finite(diagnostic_prediction)] <- 0

  diagnostic_spacing <- if (!is.null(first_peak_x) &&
      is.finite(first_peak_x) && first_peak_x > 0) {
    first_peak_x
  } else if (length(miu) > 1L) {
    stats::median(diff(sort(miu)))
  } else max(diff(range(x_observed)), 1)
  shape_regions <- build_peak_shape_regions(
    diagnostic_x, diagnostic_y, detected_peak_x, diagnostic_spacing
  )
  shape_diagnostics <- evaluate_peak_shape(
    diagnostic_prediction, miu, shape_regions, diagnostic_spacing
  )
  region_diagnostics <- evaluate_harmonic_regions(
    diagnostic_x, diagnostic_y, diagnostic_prediction, diagnostic_spacing
  )

  y_scale <- max(diagnostic_y, na.rm = TRUE)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1

  # Preserve the visually important high-count bins without completely
  # ignoring valleys and later peaks.
  point_weights <- 0.25 + 0.75 * sqrt(pmax(diagnostic_y, 0) / y_scale)
  scaled_error <- (diagnostic_prediction - diagnostic_y) / y_scale
  nrmse <- sqrt(sum(point_weights * scaled_error^2) / sum(point_weights))
  global_overshoot_nrmse <- sqrt(
    sum(point_weights * pmax(scaled_error, 0)^2) / sum(point_weights)
  )
  global_underfit_nrmse <- sqrt(
    sum(point_weights * pmax(-scaled_error, 0)^2) / sum(point_weights)
  )

  peak_indices <- which(diagnostic_x >= main_peak_left &
                        diagnostic_x <= main_peak_right)
  if (length(peak_indices) >= 3) {
    main_peak_nrmse <- sqrt(mean(scaled_error[peak_indices]^2))
    observed_peak_height <- max(diagnostic_y[peak_indices], na.rm = TRUE)
    predicted_peak_height <- max(diagnostic_prediction[peak_indices], na.rm = TRUE)
    peak_height_bias <- (predicted_peak_height - observed_peak_height) /
      max(observed_peak_height, 1)
    peak_height_error <- abs(peak_height_bias)
    peak_height_ratio <- predicted_peak_height / max(observed_peak_height, 1)
  } else {
    main_peak_nrmse <- nrmse
    peak_height_error <- nrmse
    peak_height_bias <- NA_real_
    peak_height_ratio <- NA_real_
  }
  main_peak_underfit <- if (is.finite(peak_height_bias)) {
    pmax(-peak_height_bias, 0)
  } else nrmse

  # Evaluate signed height errors at every data-supported harmonic peak against
  # the same smoothed grey curve shown in the PDF. This exposes systematic
  # under-fitting that an absolute main-peak error alone cannot detect.
  rel <- as.numeric(sigma_reliability)
  clear_peak_height_error <- peak_height_error
  clear_peak_height_bias <- if (is.finite(peak_height_bias)) peak_height_bias else 0
  clear_peak_underfit <- main_peak_underfit
  if (length(rel) == length(alpha) && length(alpha) > 0) {
    harmonic_spacing <- if (length(miu) > 1) {
      stats::median(diff(sort(miu)))
    } else max(diff(range(diagnostic_x)), 1)
    local_peak_indices <- vapply(seq_along(miu), function(j) {
      local <- which(diagnostic_x >= miu[j] - 0.18 * harmonic_spacing &
                     diagnostic_x <= miu[j] + 0.18 * harmonic_spacing)
      if (length(local) == 0) which.min(abs(diagnostic_x - miu[j])) else
        local[which.max(diagnostic_y[local])]
    }, integer(1))
    observed_heights <- diagnostic_y[local_peak_indices]
    relative_biases <- (
      diagnostic_prediction[local_peak_indices] - observed_heights
    ) / pmax(observed_heights, 0.05 * y_scale)
    clear_weights <- pmin(1, pmax(0, rel))^2
    usable <- is.finite(relative_biases) & clear_weights > 0
    if (any(usable)) {
      clear_peak_height_error <- stats::weighted.mean(
        abs(relative_biases[usable]), clear_weights[usable]
      )
      clear_peak_height_bias <- stats::weighted.mean(
        relative_biases[usable], clear_weights[usable]
      )
      clear_peak_underfit <- stats::weighted.mean(
        pmax(-relative_biases[usable], 0), clear_weights[usable]
      )
    }
  }

  log_scale <- max(log1p(y_scale), 1)
  log_rmse <- sqrt(mean(((log1p(diagnostic_prediction) -
                           log1p(diagnostic_y)) /
                         log_scale)^2))
  area_error <- abs(sum(predicted) - sum(y_observed)) /
    max(sum(y_observed), 1)

  width_parameter_count <- switch(
    sigma_mode,
    common = 1,
    trend = min(2, length(alpha)),
    adaptive = length(alpha),
    legacy = length(alpha),
    length(alpha)
  )
  parameter_count <- 2 * length(alpha) - 1 + width_parameter_count
  complexity_penalty <- 0.01 * parameter_count *
    log(max(length(x_observed), 2)) / max(length(x_observed), 1)

  # In a weakly identifiable two/three-component fit, residual error alone can
  # be reduced by making one component swallow a neighbour (the bird failure).
  # The biologically relevant warning is a sudden adjacent-width change, not
  # merely that every overlapping peak is moderately broad.  Penalise smooth
  # log-width differences continuously and reserve the absolute-width term for
  # only genuinely near-swallowing components.  This avoids another arbitrary
  # fixed adjacent ratio such as 1.35.
  width_identifiability_penalty <- 0
  rel <- as.numeric(sigma_reliability)
  if (!identical(sigma_mode, "legacy") && length(alpha) <= 3 &&
      length(rel) == length(alpha) && mean(rel, na.rm = TRUE) < 0.35 &&
      length(miu) > 1) {
    harmonic_spacing <- stats::median(diff(sort(miu)))
    if (is.finite(harmonic_spacing) && harmonic_spacing > 0) {
      rel <- pmin(1, pmax(0, rel))
      sorted <- order(miu)
      sigma_sorted <- sigma[sorted]
      rel_sorted <- rel[sorted]
      pair_reliability <- pmin(rel_sorted[-1L],
                               rel_sorted[-length(rel_sorted)])
      adjacent_penalty <- mean(
        (1 - pair_reliability) * diff(log(sigma_sorted))^2
      )
      broad_excess <- mean(
        (1 - rel_sorted) *
          pmax(sigma_sorted / harmonic_spacing - 0.60, 0)^2
      )
      width_identifiability_penalty <-
        0.12 * adjacent_penalty + 0.15 * broad_excess
    }
  }

  score <- 0.40 * nrmse + 0.15 * main_peak_nrmse +
    0.10 * peak_height_error + 0.15 * clear_peak_height_error +
    0.15 * main_peak_underfit + 0.15 * clear_peak_underfit +
    0.15 * log_rmse +
    0.08 * shape_diagnostics$left_overshoot_error +
    0.04 * shape_diagnostics$left_relative_error +
    0.05 * region_diagnostics$max_nrmse +
    0.05 * area_error + complexity_penalty +
    width_identifiability_penalty

  list(
    score = score,
    nrmse = nrmse,
    main_peak_nrmse = main_peak_nrmse,
    peak_height_error = peak_height_error,
    peak_height_bias = peak_height_bias,
    peak_height_ratio = peak_height_ratio,
    main_peak_underfit = main_peak_underfit,
    clear_peak_height_error = clear_peak_height_error,
    clear_peak_height_bias = clear_peak_height_bias,
    clear_peak_underfit = clear_peak_underfit,
    peak_position_error = shape_diagnostics$position_error,
    primary_peak_position_error =
      shape_diagnostics$primary_peak_position_error,
    total_peak_position_error =
      shape_diagnostics$total_peak_position_error,
    primary_total_peak_position_error =
      shape_diagnostics$primary_total_peak_position_error,
    primary_total_peak_signed_error =
      shape_diagnostics$primary_total_peak_signed_error,
    flank_overshoot_error = shape_diagnostics$flank_overshoot_error,
    left_overshoot_error = shape_diagnostics$left_overshoot_error,
    right_overshoot_error = shape_diagnostics$right_overshoot_error,
    left_relative_error = shape_diagnostics$left_relative_error,
    global_overshoot_nrmse = global_overshoot_nrmse,
    global_underfit_nrmse = global_underfit_nrmse,
    region_mean_nrmse = region_diagnostics$mean_nrmse,
    region_max_nrmse = region_diagnostics$max_nrmse,
    log_rmse = log_rmse,
    area_error = area_error,
    complexity_penalty = complexity_penalty,
    width_identifiability_penalty = width_identifiability_penalty
  )
}

# Internal automatic search. Candidate fits are evaluated without plots and
# only the selected configuration is fitted again to create the final PDF.
.GenomeGMM_auto <- function(data_file, n = NULL, n_candidates = 2:8,
                           search = c("fast", "thorough"),
                           symmetry_tolerance = NULL,
                           position_mode = c("adaptive", "legacy", "relaxed"),
                           joint_refinement = c("auto", "off", "always"),
                           species_name = "",
                           seed = 1, save_search = TRUE,
                           candidate_max_iter = 3000,
                           candidate_tolerance = 1e-5,
                           verification_top = 2,
                           verification_max_iter = 100000,
                           final_max_iter = NULL,
                           save_plot = TRUE) {
  search <- match.arg(search)
  position_mode <- match.arg(position_mode)
  joint_refinement <- match.arg(joint_refinement)
  if (identical(n, "auto")) n <- NULL

  if (!file.exists(data_file)) {
    stop(paste("Data file '", data_file, "' does not exist, please check the file path."))
  }
  if (!is.null(n) && (!is.numeric(n) || length(n) != 1 ||
                      !is.finite(n) || n <= 0 || n != round(n))) {
    stop("Parameter 'n' must be a positive integer, NULL, or 'auto'.")
  }
  n_candidates <- sort(unique(as.integer(n_candidates)))
  n_candidates <- n_candidates[is.finite(n_candidates) & n_candidates > 0]
  if (is.null(n) && length(n_candidates) == 0) {
    stop("At least one positive integer must be supplied in n_candidates.")
  }

  # Histograms can contain hundreds of thousands of rows. Read once and reuse
  # it for every candidate; this changes no fitting logic or plotted content.
  cached_full_data <- read.table(data_file, header = FALSE)
  auto_analysis_limit <- local({
    y_values <- cached_full_data[, 2]
    row_count <- length(y_values)
    gradients <- diff(y_values) / pmax(y_values[-row_count], 1)
    start_index <- 1L
    if (length(gradients) >= 3L) {
      for (i in seq_len(length(gradients) - 2L)) {
        if (mean(gradients[i:(i + 2L)], na.rm = TRUE) > 0.005) {
          start_index <- i
          break
        }
      }
    }
    valley_end <- min(row_count, start_index + 3L)
    valley_window <- start_index:valley_end
    start_index <- valley_window[which.min(y_values[valley_window])]
    probe_end <- min(nrow(cached_full_data), start_index + 4095L)
    probe <- cached_full_data[start_index:probe_end, , drop = FALSE]
    smooth <- as.numeric(stats::filter(probe[, 2], rep(1 / 5, 5), sides = 2))
    smooth[!is.finite(smooth)] <- probe[!is.finite(smooth), 2]
    peak_indices <- which(diff(sign(diff(smooth))) < 0) + 1L
    if (length(peak_indices) == 0L) peak_indices <- which.max(smooth)
    strongest_peak <- probe[peak_indices[which.max(smooth[peak_indices])], 1]
    effective_n <- if (is.null(n)) max(n_candidates) else n
    max(strongest_peak * (effective_n + 3), probe[1, 1] + 100)
  })

  cat("\n=== GenomeGMM automatic parameter search ===\n")
  if (is.null(n)) {
    cat("Stage 1: comparing n candidates:", n_candidates, "\n")
  } else {
    cat("Using supplied n =", n, "and tuning curve parameters.\n")
  }

  cache <- new.env(parent = emptyenv())
  trial_counter <- 0L

  candidate_key <- function(n_value, ratio, weight, symmetry_value,
                            extra_mode, sigma_mode_value, stage) {
    ratio_text <- if (is.na(ratio)) "none" else sprintf("%.5f", ratio)
    paste(n_value, ratio_text, sprintf("%.5f", weight),
          sprintf("%.5f", symmetry_value), extra_mode, sigma_mode_value,
          stage, sep = "|")
  }

  run_candidate <- function(n_value, ratio = NA_real_, weight = 1,
                            symmetry_value = 0.2, extra_mode = "legacy",
                            sigma_mode_value = "adaptive",
                            stage = c("screening", "verification")) {
    stage <- match.arg(stage)
    key <- candidate_key(n_value, ratio, weight, symmetry_value,
                         extra_mode, sigma_mode_value, stage)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    trial_counter <<- trial_counter + 1L
    call_args <- list(
      data_file = data_file,
      n = as.integer(n_value),
      symmetry_tolerance = symmetry_value,
      main_peak_weight = weight,
      use_calibration = FALSE,
      species_name = species_name,
      sigma_mode = sigma_mode_value,
      position_mode = position_mode,
      joint_refinement = "off",
      extra_component_mode = extra_mode,
      .preloaded_data = cached_full_data,
      .analysis_max_x = auto_analysis_limit,
      save_plot = FALSE
    )
    if (stage == "screening") {
      call_args$em_max_iter <- candidate_max_iter
      call_args$em_tolerance <- candidate_tolerance
    } else if (!is.null(verification_max_iter)) {
      call_args$em_max_iter <- verification_max_iter
    }
    if (!is.na(ratio)) call_args$main_peak_ratio <- ratio

    fit <- NULL
    error_message <- NA_character_
    # Use the same random stream for every candidate.  Different seeds would
    # confound parameter quality with random component reinitialisation.
    set.seed(seed)
    tryCatch(
      invisible(utils::capture.output(
        fit <- do.call(.GenomeGMM_fit, call_args), type = "output"
      )),
      error = function(e) error_message <<- conditionMessage(e)
    )

    succeeded <- !is.null(fit) && is.finite(fit$fit_score)
    raw_fit_score <- if (succeeded) fit$fit_score else Inf
    score <- raw_fit_score
    if (succeeded && !isTRUE(fit$converged)) score <- score + 0.01

    row <- data.frame(
      trial = trial_counter,
      stage = stage,
      n = as.integer(n_value),
      calibration = !is.na(ratio),
      main_peak_ratio = if (is.na(ratio)) NA_real_ else ratio,
      main_peak_weight = weight,
      symmetry_tolerance = symmetry_value,
      extra_component = if (succeeded) isTRUE(fit$fit_extra_component) else FALSE,
      sigma_mode = sigma_mode_value,
      position_mode = position_mode,
      max_sigma_ratio = if (sigma_mode_value == "legacy") 1.35 else NA_real_,
      mean_sigma_reliability = if (succeeded) {
        mean(fit$sigma_reliability, na.rm = TRUE)
      } else NA_real_,
      fit_score = raw_fit_score,
      score = score,
      converged = if (succeeded) isTRUE(fit$converged) else FALSE,
      iterations = if (succeeded) fit$iterations else NA_integer_,
      status = if (succeeded) "ok" else error_message,
      stringsAsFactors = FALSE
    )
    value <- list(row = row, fit = fit)
    assign(key, value, envir = cache)
    value
  }

  all_rows <- function() {
    keys <- ls(cache, all.names = TRUE)
    rows <- lapply(keys, function(key) get(key, envir = cache)$row)
    table <- do.call(rbind, rows)
    table[order(table$trial), , drop = FALSE]
  }

  best_row <- function(n_value = NULL, stage = "screening") {
    table <- all_rows()
    if (!is.null(n_value)) table <- table[table$n == n_value, , drop = FALSE]
    if (!is.null(stage)) table <- table[table$stage == stage, , drop = FALSE]
    valid <- is.finite(table$score)
    if (!any(valid)) stop("All automatic candidate fits failed. Inspect the histogram and parameter ranges.")
    table[which.min(table$score), , drop = FALSE]
  }

  # When ploidy is unknown, first compare modest, identical configurations.
  # The score's complexity penalty prevents a systematic preference for large n.
  if (is.null(n)) {
    for (n_value in n_candidates) run_candidate(n_value)
    n <- best_row()$n[1]
    cat("Selected n =", n, "for detailed tuning.\n")
  }

  # Coarse exploration of the two modes used by the original function:
  # uncalibrated peak weighting, and calibrated peak-window selection.
  if (search == "fast") {
    coarse_weights <- c(0.85, 0.90, 0.95, 1.00, 1.05, 1.10)
    coarse_ratios <- c(0.30, 0.40, 0.50)
  } else {
    coarse_weights <- seq(0.80, 1.20, by = 0.05)
    coarse_ratios <- seq(0.25, 0.55, by = 0.05)
  }
  initial_symmetry <- if (is.null(symmetry_tolerance)) 0.2 else symmetry_tolerance
  # Preserve the original biological gate: visible peaks >= n means exactly n
  # components; only too few detected peaks trigger the historical n+1 model.
  for (extra_mode in "legacy") {
    for (weight in coarse_weights) {
      run_candidate(n, weight = weight, symmetry_value = initial_symmetry,
                    extra_mode = extra_mode)
    }
    for (ratio in coarse_ratios) {
      run_candidate(n, ratio = ratio, weight = 1,
                    symmetry_value = initial_symmetry,
                    extra_mode = extra_mode)
    }
  }

  current_best <- best_row(n)
  best_ratio <- current_best$main_peak_ratio[1]
  best_weight <- current_best$main_peak_weight[1]
  best_symmetry <- current_best$symmetry_tolerance[1]
  best_extra <- "legacy"

  # Local refinement around the best coarse peak weight.
  refined_weights <- sort(unique(pmin(1.40, pmax(
    0.60, best_weight + c(-0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04)
  ))))
  for (weight in refined_weights) {
    run_candidate(n, ratio = best_ratio, weight = weight,
                  symmetry_value = best_symmetry, extra_mode = best_extra)
  }

  current_best <- best_row(n)
  best_ratio <- current_best$main_peak_ratio[1]
  best_weight <- current_best$main_peak_weight[1]
  best_symmetry <- current_best$symmetry_tolerance[1]
  best_extra <- "legacy"

  # A calibrated fit also gets a local peak-window refinement.
  if (!is.na(best_ratio)) {
    ratio_offsets <- if (search == "fast") {
      c(-0.04, -0.02, 0, 0.02, 0.04)
    } else {
      seq(-0.05, 0.05, by = 0.01)
    }
    refined_ratios <- sort(unique(pmin(0.65, pmax(0.15, best_ratio + ratio_offsets))))
    for (ratio in refined_ratios) {
      run_candidate(n, ratio = ratio, weight = best_weight,
                    symmetry_value = best_symmetry, extra_mode = best_extra)
    }
  }

  current_best <- best_row(n)
  best_ratio <- current_best$main_peak_ratio[1]
  best_weight <- current_best$main_peak_weight[1]
  best_extra <- "legacy"

  symmetry_candidates <- if (!is.null(symmetry_tolerance)) {
    symmetry_tolerance
  } else if (search == "fast") {
    c(0.10, 0.20, 0.35, 0.50)
  } else {
    seq(0.05, 0.50, by = 0.05)
  }
  for (symmetry_value in symmetry_candidates) {
    run_candidate(n, ratio = best_ratio, weight = best_weight,
                  symmetry_value = symmetry_value, extra_mode = best_extra)
  }

  # Re-evaluate the original peak-count gate at the best nuisance parameters.
  current_best <- best_row(n)
  best_ratio <- current_best$main_peak_ratio[1]
  best_weight <- current_best$main_peak_weight[1]
  best_symmetry <- current_best$symmetry_tolerance[1]
  run_candidate(n, ratio = best_ratio, weight = best_weight,
                symmetry_value = best_symmetry, extra_mode = "legacy")

  # Compare parsimonious common-width and power-trend models with the adaptive
  # soft model at the best nuisance parameters. This is especially important
  # for Chinese sturgeon: the early visible peaks can support a trend while the
  # weak late peaks should not each spend an unconstrained width parameter.
  current_best <- best_row(n)
  best_ratio <- current_best$main_peak_ratio[1]
  best_weight <- current_best$main_peak_weight[1]
  best_symmetry <- current_best$symmetry_tolerance[1]
  best_extra <- "legacy"
  for (sigma_mode_value in c("common", "trend", "adaptive")) {
    run_candidate(n, ratio = best_ratio, weight = best_weight,
                  symmetry_value = best_symmetry, extra_mode = best_extra,
                  sigma_mode_value = sigma_mode_value)
  }

  # The constrained EM used by this project can continue moving even after a
  # few thousand iterations. Recheck a small finalist set with the original
  # full iteration budget so screening speed does not decide the winner.
  verification_top <- max(1L, as.integer(verification_top))
  screening_table <- all_rows()
  screening_table <- screening_table[
    screening_table$n == n & screening_table$stage == "screening" &
      is.finite(screening_table$score), , drop = FALSE
  ]
  screening_table <- screening_table[order(screening_table$score), , drop = FALSE]
  finalists <- head(screening_table, verification_top)

  # Also verify the immediate weight neighbours of the screening winner. This
  # catches narrow optima such as 0.93 versus 0.94 in the supplied example.
  winner <- screening_table[1, , drop = FALSE]
  neighbour_weights <- pmin(1.40, pmax(
    0.60, winner$main_peak_weight[1] + c(-0.01, 0.01)
  ))
  neighbours <- winner[rep(1, length(neighbour_weights)), , drop = FALSE]
  neighbours$main_peak_weight <- neighbour_weights
  finalists <- rbind(finalists, neighbours)
  finalist_key <- paste(finalists$n, finalists$calibration,
                        finalists$main_peak_ratio, finalists$main_peak_weight,
                        finalists$symmetry_tolerance,
                        finalists$extra_component, finalists$sigma_mode,
                        sep = "|")
  finalists <- finalists[!duplicated(finalist_key), , drop = FALSE]

  cat("Full-precision verification of", nrow(finalists), "finalists...\n")
  for (i in seq_len(nrow(finalists))) {
    run_candidate(
      finalists$n[i], ratio = finalists$main_peak_ratio[i],
      weight = finalists$main_peak_weight[i],
      symmetry_value = finalists$symmetry_tolerance[i],
      extra_mode = "legacy",
      sigma_mode_value = finalists$sigma_mode[i],
      stage = "verification"
    )
  }

  search_table <- all_rows()
  selected <- best_row(n, stage = "verification")
  search_table$selected <- search_table$trial == selected$trial[1]

  cat("\nBest automatic parameters:\n")
  cat("  n =", selected$n[1], "\n")
  cat("  calibration =", selected$calibration[1], "\n")
  if (selected$calibration[1]) {
    cat("  main_peak_ratio =", selected$main_peak_ratio[1], "\n")
  }
  cat("  main_peak_weight =", selected$main_peak_weight[1], "\n")
  cat("  symmetry_tolerance =", selected$symmetry_tolerance[1], "\n")
  cat("  optional n+1 component =", selected$extra_component[1], "\n")
  cat("  sigma width model =", selected$sigma_mode[1], "\n")
  cat("  peak position model =", position_mode, "\n")
  cat("  mean peak-width reliability =",
      round(selected$mean_sigma_reliability[1], 3), "\n")
  cat("  fit score =", round(selected$fit_score[1], 6), "\n")
  if (!selected$converged[1]) {
    cat("  selection score (includes non-convergence penalty) =",
        round(selected$score[1], 6), "\n")
  }
  cat("\n")

  # Refit once with normal output and create only the winning plot.
  final_args <- list(
    data_file = data_file,
    n = selected$n[1],
    symmetry_tolerance = selected$symmetry_tolerance[1],
    main_peak_weight = selected$main_peak_weight[1],
    use_calibration = FALSE,
    species_name = species_name,
    sigma_mode = selected$sigma_mode[1],
    position_mode = position_mode,
    joint_refinement = joint_refinement,
    extra_component_mode = "legacy",
    .preloaded_data = cached_full_data,
    .analysis_max_x = auto_analysis_limit,
    save_plot = save_plot,
    em_max_iter = final_max_iter
  )
  if (selected$calibration[1]) {
    final_args$main_peak_ratio <- selected$main_peak_ratio[1]
  }
  set.seed(seed)
  final_fit <- do.call(.GenomeGMM_fit, final_args)

  search_file <- NULL
  if (isTRUE(save_search)) {
    file_basename <- tools::file_path_sans_ext(basename(data_file))
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    search_file <- paste0("GenomeGMM_search_", file_basename,
                          "_n", selected$n[1], "_", timestamp, ".csv")
    utils::write.csv(search_table, search_file, row.names = FALSE, na = "")
    cat("Parameter search table saved as:", search_file, "\n")
  }

  final_fit$selected_parameters <- list(
    n = selected$n[1],
    use_calibration = selected$calibration[1],
    main_peak_ratio = selected$main_peak_ratio[1],
    main_peak_weight = selected$main_peak_weight[1],
    symmetry_tolerance = selected$symmetry_tolerance[1],
    fit_extra_component = selected$extra_component[1],
    sigma_mode = selected$sigma_mode[1],
    position_mode = position_mode,
    joint_refinement = joint_refinement,
    joint_triggered = isTRUE(final_fit$joint_triggered),
    joint_accepted = isTRUE(final_fit$joint_accepted),
    joint_profile = final_fit$joint_profile,
    mean_sigma_reliability = selected$mean_sigma_reliability[1],
    max_sigma_ratio = if (selected$sigma_mode[1] == "legacy") 1.35 else NA_real_,
    score = selected$fit_score[1],
    selection_score = selected$score[1]
  )
  final_fit$search_table <- search_table
  final_fit$search_file <- search_file
  invisible(final_fit)
}

#' Automatically fit a k-mer histogram with GenomeGMM
#'
#' Only the histogram path and known biological ploidy are required. All
#' numerical fitting parameters are selected automatically. `species_name` is
#' optional and only changes the plot title.
#'
#' @param data_file Path to a two-column k-mer histogram file.
#' @param n Known positive integer biological ploidy/component count.
#' @param species_name Optional species name shown in the plot title.
#'
#' @return An invisible list containing the fitted model, automatic parameter
#'         selection table, genome-size estimates, and output file paths.
#'
#' @examples
#' # GenomeGMM("sample.histo", n = 4)
#' # GenomeGMM("sample.histo", n = 4, species_name = "Example")
#'
#' @export
GenomeGMM <- function(data_file, n, species_name = "") {
  if (missing(n)) {
    stop("Parameter 'n' (known biological ploidy) must be specified.")
  }
  if (!is.character(species_name) || length(species_name) != 1L ||
      is.na(species_name)) {
    stop("Parameter 'species_name' must be one character string.")
  }

  .GenomeGMM_auto(
    data_file = data_file,
    n = n,
    species_name = species_name,
    search = "fast",
    candidate_max_iter = 1500,
    candidate_tolerance = 1e-5,
    verification_top = 4,
    verification_max_iter = 10000,
    final_max_iter = 10000,
    save_search = TRUE,
    save_plot = TRUE
  )
}

# ============================================================================
# Disabled experimental fitter (kept only as an internal rollback record).
# It is deliberately not evaluated and cannot replace the production fitter.
# ============================================================================
if (FALSE) {

# Non-negative weighted least squares by cyclic coordinate descent.  Keeping
# this implementation local avoids adding an external R-package dependency.
.ggmm_nnls <- function(design, response, weights = NULL, max_iter = 120,
                       tolerance = 1e-8) {
  design <- as.matrix(design)
  response <- as.numeric(response)
  if (is.null(weights)) weights <- rep(1, length(response))
  sqrt_weights <- sqrt(pmax(as.numeric(weights), 0))
  xw <- design * sqrt_weights
  yw <- response * sqrt_weights
  coefficients <- rep(0, ncol(xw))
  fitted <- rep(0, length(yw))
  denominators <- colSums(xw^2) + 1e-12

  for (iteration in seq_len(max_iter)) {
    previous <- coefficients
    for (j in seq_len(ncol(xw))) {
      partial_residual <- yw - fitted + xw[, j] * coefficients[j]
      updated <- max(0, sum(xw[, j] * partial_residual) / denominators[j])
      fitted <- fitted + xw[, j] * (updated - coefficients[j])
      coefficients[j] <- updated
    }
    if (max(abs(coefficients - previous)) < tolerance *
        max(1, max(abs(previous)))) break
  }
  coefficients
}

.ggmm_smooth <- function(values, window = 5L) {
  window <- max(3L, as.integer(window))
  if (window %% 2L == 0L) window <- window + 1L
  smoothed <- as.numeric(stats::filter(values, rep(1 / window, window),
                                       sides = 2))
  smoothed[!is.finite(smoothed)] <- values[!is.finite(smoothed)]
  smoothed
}

.ggmm_read_histogram <- function(data_file) {
  raw <- utils::read.table(data_file, header = FALSE,
                           colClasses = c("numeric", "numeric"))
  if (ncol(raw) < 2L) stop("The histogram must contain at least two columns.")
  raw <- raw[, 1:2]
  names(raw) <- c("x", "y")
  raw <- raw[is.finite(raw$x) & is.finite(raw$y) & raw$x >= 0 & raw$y >= 0, ]
  if (nrow(raw) < 10L) stop("The histogram contains fewer than 10 valid rows.")
  raw <- stats::aggregate(y ~ x, data = raw, FUN = sum)
  raw <- raw[order(raw$x), ]
  rownames(raw) <- NULL
  raw
}

# Locate the valley after the sequencing-error peak, then collect real genomic
# peak candidates.  No absolute coverage or count threshold is used.
.ggmm_prepare_data <- function(histogram, n) {
  x <- histogram$x
  y <- histogram$y
  number_of_rows <- length(y)

  relative_gradient <- diff(y) / pmax(y[-number_of_rows], 1)
  start_index <- 2L
  if (length(relative_gradient) >= 3L) {
    for (i in seq_len(length(relative_gradient) - 2L)) {
      local_gradient <- mean(relative_gradient[i:(i + 2L)], na.rm = TRUE)
      if (is.finite(local_gradient) && local_gradient > 0.005) {
        # i is the last falling bin; the following bin is the valley where the
        # genomic spectrum starts. Including bin i would keep part of the
        # sequencing-error wall and can dominate the fit (notably at low depth).
        start_index <- max(2L, i + 1L)
        break
      }
    }
  }

  segment_indices <- start_index:number_of_rows
  segment_y <- y[segment_indices]
  segment_x <- x[segment_indices]
  smooth_y <- .ggmm_smooth(segment_y, 5L)

  local_peak_indices <- which(diff(sign(diff(smooth_y))) < 0) + 1L
  # Moving-average edge artefacts cannot be biological peaks.
  local_peak_indices <- local_peak_indices[
    local_peak_indices > 2L & local_peak_indices < length(smooth_y) - 2L
  ]
  if (length(local_peak_indices) == 0L) local_peak_indices <- which.max(smooth_y)

  local_heights <- smooth_y[local_peak_indices]
  height_floor <- max(local_heights, na.rm = TRUE) * 0.002
  local_peak_indices <- local_peak_indices[local_heights >= height_floor]
  local_heights <- smooth_y[local_peak_indices]

  # Retain the strongest separated candidates. Closely spaced numerical maxima
  # are redundant for estimating a harmonic base depth.
  order_by_height <- order(local_heights, decreasing = TRUE)
  local_peak_indices <- local_peak_indices[order_by_height]
  kept <- integer(0)
  minimum_separation <- max(stats::median(diff(x)), 1) * 3
  for (index in local_peak_indices) {
    if (length(kept) == 0L ||
        all(abs(segment_x[index] - segment_x[kept]) >= minimum_separation)) {
      kept <- c(kept, index)
    }
    if (length(kept) >= max(12L, 2L * n + 4L)) break
  }
  local_peak_indices <- kept
  peak_x <- segment_x[local_peak_indices]
  peak_y <- smooth_y[local_peak_indices]

  strongest_peak_x <- peak_x[which.max(peak_y)]
  genomic_height <- max(peak_y)
  above_floor <- which(smooth_y >= genomic_height * 0.01)
  decay_end_x <- if (length(above_floor) > 0L) {
    segment_x[max(above_floor)]
  } else {
    strongest_peak_x * (n + 2)
  }
  # A common comparison range is essential. If every candidate ended at a
  # multiple of its own proposed depth, subharmonics could win merely by
  # ignoring the right-hand part of the observed curve.
  fit_end_x <- min(decay_end_x, strongest_peak_x * (n + 2))
  fit_end_x <- max(fit_end_x, strongest_peak_x * 1.25)

  list(
    start_index = start_index,
    start_x = x[start_index],
    segment_x = segment_x,
    segment_y = segment_y,
    smooth_y = smooth_y,
    peak_x = peak_x,
    peak_y = peak_y,
    fit_end_x = fit_end_x
  )
}

.ggmm_basis <- function(x, depth, base_sigma, width_exponent, components) {
  component_index <- seq_len(components)
  means <- depth * component_index
  sigmas <- base_sigma * component_index^width_exponent
  basis <- vapply(component_index, function(j) {
    exp(-0.5 * ((x - means[j]) / sigmas[j])^2)
  }, numeric(length(x)))
  if (components == 1L) basis <- matrix(basis, ncol = 1L)
  list(basis = basis, means = means, sigmas = sigmas)
}

.ggmm_curve_metrics <- function(observed, predicted) {
  observed <- pmax(as.numeric(observed), 0)
  predicted <- pmax(as.numeric(predicted), 0)
  scale <- max(observed, 1)
  observed_scaled <- observed / scale
  predicted_scaled <- predicted / scale
  importance <- 0.20 + 0.80 * sqrt(observed_scaled)
  nrmse <- sqrt(sum(importance * (predicted_scaled - observed_scaled)^2) /
                  sum(importance))
  log_rmse <- sqrt(mean((log1p(1000 * predicted_scaled) -
                         log1p(1000 * observed_scaled))^2)) / log(1001)
  area_error <- abs(sum(predicted) - sum(observed)) / max(sum(observed), 1)
  rss <- sum((predicted - observed)^2)
  tss <- sum((observed - mean(observed))^2)
  r_squared <- if (tss > 0) 1 - rss / tss else NA_real_
  score <- 0.65 * nrmse + 0.25 * log_rmse + 0.10 * area_error
  list(score = score, nrmse = nrmse, log_rmse = log_rmse,
       area_error = area_error, r_squared = r_squared)
}

.ggmm_fit_at_parameters <- function(parameters, x, y, components,
                                    weight_power) {
  depth <- parameters[1]
  base_sigma <- parameters[2]
  width_exponent <- parameters[3]
  model <- .ggmm_basis(x, depth, base_sigma, width_exponent, components)
  y_scale <- max(y, 1)
  y_scaled <- y / y_scale
  regression_weights <- (0.03 + y_scaled)^(-weight_power)
  amplitudes_scaled <- .ggmm_nnls(model$basis, y_scaled,
                                  regression_weights)
  predicted_scaled <- as.numeric(model$basis %*% amplitudes_scaled)
  predicted <- predicted_scaled * y_scale
  metrics <- .ggmm_curve_metrics(y, predicted)
  list(
    objective = metrics$score,
    predicted = predicted,
    amplitudes = amplitudes_scaled * y_scale,
    means = model$means,
    sigmas = model$sigmas,
    metrics = metrics
  )
}

.ggmm_optimize_candidate <- function(histogram, prepared, n, depth_initial,
                                     components, weight_power) {
  x_end <- min(max(histogram$x), prepared$fit_end_x)
  fit_indices <- which(histogram$x >= prepared$start_x & histogram$x <= x_end)
  if (length(fit_indices) < 10L) return(NULL)
  x_fit <- histogram$x[fit_indices]
  y_fit <- histogram$y[fit_indices]

  lower <- c(depth_initial * 0.82, max(0.35, depth_initial * 0.018), 0)
  upper <- c(depth_initial * 1.18, max(0.75, depth_initial * 0.65), 0.75)
  initial <- c(depth_initial, max(0.75, depth_initial * 0.12), 0.25)

  objective <- function(parameters) {
    result <- .ggmm_fit_at_parameters(parameters, x_fit, y_fit,
                                      components, weight_power)
    # Very broad first peaks erase biologically meaningful harmonic structure.
    # Keep each optimization in the peak-derived depth basin; multi-start
    # selection, rather than boundary drift, decides which harmonic is 1x.
    width_ratio <- parameters[2] / parameters[1]
    broadness_penalty <- 2.0 * max(0, width_ratio - 0.50)^2
    basin_penalty <- 0.02 * ((parameters[1] - depth_initial) /
                             (0.10 * depth_initial))^2
    result$objective + broadness_penalty + basin_penalty
  }

  optimized <- tryCatch(
    stats::optim(
      initial, objective, method = "L-BFGS-B", lower = lower, upper = upper,
      control = list(maxit = 120, factr = 1e8,
                     parscale = c(depth_initial, depth_initial * 0.1, 0.25))
    ),
    error = function(e) NULL
  )
  if (is.null(optimized) || any(!is.finite(optimized$par))) return(NULL)

  final <- .ggmm_fit_at_parameters(optimized$par, x_fit, y_fit,
                                   components, weight_power)
  # One nuisance harmonic is allowed, but only retained when it earns its cost.
  complexity_penalty <- 0.0025 * max(0, components - n)
  final$selection_score <- final$metrics$score + complexity_penalty
  final$depth <- optimized$par[1]
  final$base_sigma <- optimized$par[2]
  final$width_exponent <- optimized$par[3]
  final$components <- components
  final$weight_power <- weight_power
  final$x_fit <- x_fit
  final$y_fit <- y_fit
  final$x_end <- max(x_fit)
  final$convergence <- optimized$convergence
  final$optim_counts <- optimized$counts
  final
}

.ggmm_depth_candidates <- function(histogram, prepared, n) {
  peak_x <- prepared$peak_x
  peak_y <- prepared$peak_y
  if (length(peak_x) == 0L) stop("No genomic peak candidate could be detected.")

  candidates <- as.vector(outer(peak_x, seq_len(n), "/"))
  bin_width <- stats::median(diff(histogram$x))
  candidates <- candidates[
    is.finite(candidates) & candidates >= max(2 * bin_width, 1) &
      candidates <= max(peak_x) * 1.1
  ]
  candidates <- sort(candidates)
  # Merge candidates that differ by less than 3%; their optimization basins
  # overlap and repeating them adds time without information.
  merged <- numeric(0)
  for (candidate in candidates) {
    if (length(merged) == 0L ||
        min(abs(log(candidate / merged))) > log(1.03)) {
      merged <- c(merged, candidate)
    }
  }

  # Cheap preliminary fits retain the most plausible basins for full
  # optimization. Always retain divisions of the strongest observed peak.
  preliminary <- lapply(merged, function(depth) {
    x_end <- min(max(histogram$x), prepared$fit_end_x)
    indices <- which(histogram$x >= prepared$start_x & histogram$x <= x_end)
    if (length(indices) < 10L) return(Inf)
    x_fit <- histogram$x[indices]
    y_fit <- histogram$y[indices]
    trial <- .ggmm_fit_at_parameters(
      c(depth, max(0.75, 0.12 * depth), 0.25), x_fit, y_fit,
      n + 1L, 0.5
    )
    trial$metrics$score
  })
  preliminary <- unlist(preliminary)
  retained <- merged[order(preliminary)][seq_len(min(8L, length(merged)))]
  strongest_peak <- peak_x[which.max(peak_y)]
  anchor_candidates <- strongest_peak / seq_len(n)
  anchor_candidates <- anchor_candidates[anchor_candidates >= max(2 * bin_width, 1)]
  retained <- sort(unique(c(retained, anchor_candidates)))
  retained
}

.ggmm_fit_universal <- function(data_file, n, verbose = TRUE) {
  histogram <- .ggmm_read_histogram(data_file)
  prepared <- .ggmm_prepare_data(histogram, n)
  depth_candidates <- .ggmm_depth_candidates(histogram, prepared, n)

  if (verbose) {
    cat("Detected genomic fit start:", prepared$start_x, "\n")
    cat("Detected peak candidates:",
        paste(round(sort(prepared$peak_x), 2), collapse = ", "), "\n")
    cat("Optimizing", length(depth_candidates), "base-depth basins...\n")
  }

  fits <- list()
  fit_number <- 0L
  for (depth in depth_candidates) {
    for (components in unique(c(n, n + 1L))) {
      for (weight_power in c(0.25, 0.75)) {
        candidate <- .ggmm_optimize_candidate(
          histogram, prepared, n, depth, components, weight_power
        )
        if (!is.null(candidate) && is.finite(candidate$selection_score)) {
          fit_number <- fit_number + 1L
          candidate$initial_depth <- depth
          fits[[fit_number]] <- candidate
        }
      }
    }
  }
  if (length(fits) == 0L) stop("All constrained optimization attempts failed.")
  scores <- vapply(fits, function(fit) fit$selection_score, numeric(1))
  best <- fits[[which.min(scores)]]
  best$histogram <- histogram
  best$prepared <- prepared
  best$depth_candidates <- depth_candidates
  best$candidate_scores <- data.frame(
    initial_depth = vapply(fits, function(fit) fit$initial_depth, numeric(1)),
    depth = vapply(fits, function(fit) fit$depth, numeric(1)),
    base_sigma = vapply(fits, function(fit) fit$base_sigma, numeric(1)),
    width_exponent = vapply(fits, function(fit) fit$width_exponent, numeric(1)),
    components = vapply(fits, function(fit) fit$components, integer(1)),
    weight_power = vapply(fits, function(fit) fit$weight_power, numeric(1)),
    fit_score = vapply(fits, function(fit) fit$metrics$score, numeric(1)),
    selection_score = scores
  )
  best
}

.ggmm_genome_size <- function(fit) {
  bin_width <- stats::median(diff(fit$histogram$x))
  component_distinct_counts <- fit$amplitudes * fit$sigmas *
    sqrt(2 * pi) / bin_width
  component_kmer_mass <- component_distinct_counts * fit$means
  component_sizes_bp <- component_kmer_mass / fit$depth

  tail <- fit$histogram[fit$histogram$x > fit$x_end, ]
  tail_kmer_mass <- if (nrow(tail) > 0L) sum(tail$x * tail$y) else 0
  fitted_kmer_mass <- sum(component_kmer_mass)
  genome_size_bp <- (fitted_kmer_mass + tail_kmer_mass) / fit$depth
  residual_size_bp <- tail_kmer_mass / fit$depth

  list(
    genome_size_bp = genome_size_bp,
    component_sizes_bp = component_sizes_bp,
    residual_size_bp = residual_size_bp,
    fitted_kmer_mass = fitted_kmer_mass,
    tail_kmer_mass = tail_kmer_mass
  )
}

.ggmm_draw_universal_plot <- function(fit, genome, n, species_name = "") {
  histogram <- fit$histogram
  plot_limit <- min(max(histogram$x), fit$depth * (n + 2))
  plot_indices <- which(histogram$x >= fit$prepared$start_x &
                        histogram$x <= plot_limit)
  plot_x <- histogram$x[plot_indices]
  plot_y <- histogram$y[plot_indices]
  bin_width <- stats::median(diff(histogram$x))

  dense_x <- seq(0, plot_limit, length.out = 2500)
  dense_model <- .ggmm_basis(dense_x, fit$depth, fit$base_sigma,
                             fit$width_exponent, fit$components)
  components <- lapply(seq_len(fit$components), function(j) {
    dense_model$basis[, j] * fit$amplitudes[j]
  })
  total_curve <- Reduce(`+`, components)

  title <- if (nzchar(species_name)) {
    paste0("Genome size estimation of ", species_name)
  } else {
    "Automatic constrained genome GMM"
  }
  y_limit <- max(plot_y, total_curve, na.rm = TRUE) * 1.18
  graphics::plot(NULL, xlim = c(0, plot_limit), ylim = c(0, y_limit),
                 main = title, xlab = "kmer frequency", ylab = "count")
  graphics::rect(plot_x - bin_width / 2, 0, plot_x + bin_width / 2, plot_y,
                 col = "#E8E8E8", border = NA)
  if (length(plot_x) > 3L) {
    graphics::lines(plot_x, .ggmm_smooth(plot_y, 5L), col = "#999999", lwd = 1.5)
  }

  line_colors <- c("#EA0000", "#E69F00", "#00A6D6", "#861FAB",
                   "#00A087", "#3C5488", "#F39B7F", "#7E6148",
                   "#42B540", "#C2185B")
  fill_colors <- grDevices::adjustcolor(line_colors, alpha.f = 0.22)
  for (j in rev(seq_len(fit$components))) {
    color_index <- ((j - 1L) %% length(line_colors)) + 1L
    graphics::polygon(c(dense_x, rev(dense_x)),
                      c(components[[j]], rep(0, length(dense_x))),
                      col = fill_colors[color_index], border = NA)
    graphics::lines(dense_x, components[[j]],
                    col = line_colors[color_index], lty = 2, lwd = 1.2)
  }
  graphics::lines(dense_x, total_curve, col = "#E41A1C", lwd = 2.2)
  for (mean_value in fit$means) {
    if (mean_value <= plot_limit) {
      graphics::abline(v = mean_value, col = "#CCCCCC", lty = 3)
    }
  }

  component_labels <- vapply(seq_len(fit$components), function(j) {
    suffix <- if (j > n) " (nuisance)" else ""
    sprintf("Component %d%s: %.1f Mb", j, suffix,
            genome$component_sizes_bp[j] / 1e6)
  }, character(1))
  legend_labels <- c("Original data", "Smoothed data", "Total fit",
                     component_labels)
  legend_colors <- c("#E8E8E8", "#999999", "#E41A1C",
                     line_colors[seq_len(fit$components)])
  legend_lty <- c(NA, 1, 1, rep(2, fit$components))
  legend_pch <- c(15, NA, NA, rep(NA, fit$components))
  graphics::legend("topright", legend = legend_labels, col = legend_colors,
                   lty = legend_lty, pch = legend_pch, bty = "n",
                   cex = max(0.45, 0.82 - 0.035 * max(0, fit$components - 5)))

  info <- sprintf(
    "Depth: %.2fx\nFull GSE: %.1f Mb\nHaploid GSE: %.1f Mb\nR2: %.4f",
    fit$depth, genome$genome_size_bp / 1e6,
    genome$genome_size_bp / n / 1e6, fit$metrics$r_squared
  )
  graphics::legend("topleft", legend = info, bty = "n", cex = 0.78,
                   text.font = 2)
  graphics::box(lwd = 1.2)
}

#' Fully automatic, biologically constrained genome histogram fitting
#'
#' Only the histogram and biological ploidy/component count n are required.
#' Peak spacing is locked to integer multiples of an automatically inferred
#' base depth. Peak heights, base width, width-growth law, fit start, nuisance
#' component, loss weighting and all numerical initializations are selected
#' from the data by multi-start constrained optimization.
#'
#' @param data_file Path to a two-column k-mer histogram.
#' @param n Biological ploidy/component count. This is deliberately required:
#'        a histogram alone cannot identify biological ploidy reliably.
#' @param species_name Optional plot title species name.
#' @param save_plot Save the final PDF plot.
#' @param save_diagnostics Save all candidate fits as a CSV table.
#' @param verbose Print automatic detection and fit information.
#'
#' @return An invisible list containing the selected model, genome-size
#'         estimates, automatic parameters, diagnostics and output paths.
#'
#' @export
GenomeGMM_auto_experimental <- function(data_file, n, species_name = "", save_plot = TRUE,
                                        save_diagnostics = TRUE, verbose = TRUE) {
  if (missing(n) || !is.numeric(n) || length(n) != 1L ||
      !is.finite(n) || n < 1 || n != round(n)) {
    stop("'n' must be the known positive integer biological ploidy/component count.")
  }
  n <- as.integer(n)
  if (!file.exists(data_file)) stop("Histogram file does not exist: ", data_file)

  if (verbose) cat("\n=== Universal automatic constrained GenomeGMM ===\n")
  fit <- .ggmm_fit_universal(data_file, n, verbose = verbose)
  genome <- .ggmm_genome_size(fit)

  diagnostics <- fit$candidate_scores[order(fit$candidate_scores$selection_score), ]
  diagnostics$selected <- seq_len(nrow(diagnostics)) == 1L
  file_stem <- tools::file_path_sans_ext(basename(data_file))
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  pdf_file <- NULL
  diagnostics_file <- NULL

  if (isTRUE(save_plot)) {
    pdf_file <- paste0("auto_gmm_", file_stem, "_n", n, "_", timestamp, ".pdf")
    grDevices::pdf(pdf_file, width = 7.08, height = 5.3, family = "Helvetica")
    tryCatch(
      .ggmm_draw_universal_plot(fit, genome, n, species_name),
      finally = grDevices::dev.off()
    )
  }
  if (isTRUE(save_diagnostics)) {
    diagnostics_file <- paste0("auto_gmm_diagnostics_", file_stem,
                               "_n", n, "_", timestamp, ".csv")
    utils::write.csv(diagnostics, diagnostics_file, row.names = FALSE)
  }

  quality <- if (is.finite(fit$metrics$r_squared) && fit$metrics$r_squared >= 0.97) {
    "excellent"
  } else if (is.finite(fit$metrics$r_squared) && fit$metrics$r_squared >= 0.90) {
    "good"
  } else {
    "review"
  }
  warnings <- character(0)
  if (quality == "review") {
    warnings <- c(warnings,
                  "The constrained harmonic model explains less than 90% of curve variance; inspect this sample.")
  }
  if (fit$base_sigma / fit$depth > 0.45) {
    warnings <- c(warnings,
                  "The first peak is very broad relative to its spacing; peaks may be weakly identifiable.")
  }

  automatic_parameters <- list(
    fit_start = fit$prepared$start_x,
    base_depth = fit$depth,
    base_sigma = fit$base_sigma,
    width_exponent = fit$width_exponent,
    sigma_by_component = fit$sigmas,
    means = fit$means,
    amplitudes = fit$amplitudes,
    fitted_components = fit$components,
    nuisance_components = max(0L, fit$components - n),
    loss_weight_power = fit$weight_power
  )

  if (verbose) {
    cat("\nSelected automatic model:\n")
    cat("  base depth:", round(fit$depth, 4), "\n")
    cat("  harmonic means:", paste(round(fit$means, 3), collapse = ", "), "\n")
    cat("  component sigmas:", paste(round(fit$sigmas, 3), collapse = ", "), "\n")
    cat("  fitted components:", fit$components, "\n")
    cat("  R-squared:", round(fit$metrics$r_squared, 5), "\n")
    cat("  fit quality:", quality, "\n")
    cat("  full genome size:", round(genome$genome_size_bp / 1e6, 2), "Mb\n")
    cat("  haploid genome size:", round(genome$genome_size_bp / n / 1e6, 2), "Mb\n")
    if (!is.null(pdf_file)) cat("  plot:", pdf_file, "\n")
    if (!is.null(diagnostics_file)) cat("  diagnostics:", diagnostics_file, "\n")
    if (length(warnings) > 0L) cat("  warning:", paste(warnings, collapse = " "), "\n")
  }

  invisible(list(
    genome_size_bp = genome$genome_size_bp,
    haploid_genome_size_bp = genome$genome_size_bp / n,
    component_sizes_bp = genome$component_sizes_bp,
    residual_size_bp = genome$residual_size_bp,
    depth = fit$depth,
    n = n,
    automatic_parameters = automatic_parameters,
    fit_metrics = fit$metrics,
    fit_quality = quality,
    warnings = warnings,
    diagnostics = diagnostics,
    pdf_file = pdf_file,
    diagnostics_file = diagnostics_file
  ))
}
}
