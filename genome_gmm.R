# Integrate total k-mer mass with the fitted red curve only while it is needed
# to correct the unreliable leading part of the spectrum.  Once the fitted
# curve and the observed grey reference agree continuously, switch back to the
# observed histogram for every remaining coverage bin.
.GenomeGMM_integrate_red_then_observed <- function(
    natural_x, fitted_counts, alpha, miu, sigma, total_count_scale,
    full_data, fit_range_upper, biological_components, bin_width,
    agreement_tolerance = 0.05) {
  observed <- data.frame(
    x = suppressWarnings(as.numeric(full_data[, 1])),
    y = suppressWarnings(as.numeric(full_data[, 2]))
  )
  observed <- observed[
    is.finite(observed$x) & is.finite(observed$y) &
      observed$x >= min(natural_x) & observed$y >= 0,
    , drop = FALSE
  ]
  if (nrow(observed)) {
    observed <- stats::aggregate(y ~ x, data = observed, FUN = sum)
    observed <- observed[order(observed$x), , drop = FALSE]
  }

  switch_x <- NA_real_
  required_agreement_bins <- 5L
  grey_counts <- rep(NA_real_, length(natural_x))
  observed_in_fit <- observed[
    observed$x <= fit_range_upper, , drop = FALSE
  ]
  if (nrow(observed_in_fit) >= 4L) {
    grey_reference <- make_gmm_reference_curve(
      observed_in_fit$x, observed_in_fit$y,
      dense_points = min(20000L, max(1000L, 4L * length(natural_x)))
    )
    grey_counts <- stats::approx(
      grey_reference$x, pmax(grey_reference$y, 0),
      xout = natural_x, rule = 2
    )$y
    peak_scale <- max(grey_counts, na.rm = TRUE)
    valid_miu <- sort(miu[biological_components][
      is.finite(miu[biological_components]) & miu[biological_components] > 0
    ])
    harmonic_spacing <- if (length(valid_miu) >= 2L) {
      stats::median(diff(valid_miu), na.rm = TRUE)
    } else if (length(valid_miu) == 1L) {
      max(2 * sigma[biological_components[1]], 5 * bin_width)
    } else {
      5 * bin_width
    }
    required_agreement_bins <- max(
      5L, min(20L, as.integer(ceiling(
        0.25 * harmonic_spacing / max(bin_width, .Machine$double.eps)
      )))
    )
    # Five per cent local agreement is required.  The small peak-scale floor
    # keeps relative error stable in a low-count tail without allowing two
    # almost-zero curves to create an accidental switch.
    local_scale <- pmax(grey_counts, 0.02 * peak_scale, 1)
    relative_difference <- abs(fitted_counts - grey_counts) / local_scale
    search_start <- if (length(valid_miu)) min(valid_miu) else min(natural_x)
    agrees <- is.finite(relative_difference) &
      natural_x >= search_start &
      grey_counts >= 0.01 * peak_scale &
      relative_difference <= agreement_tolerance
    agreement_runs <- rle(agrees)
    run_ends <- cumsum(agreement_runs$lengths)
    run_starts <- run_ends - agreement_runs$lengths + 1L
    accepted_runs <- which(
      agreement_runs$values &
        agreement_runs$lengths >= required_agreement_bins
    )
    if (length(accepted_runs)) {
      switch_x <- natural_x[run_starts[accepted_runs[1]]]
    }
  }

  # If no stable agreement exists, retain the historical boundary: use the
  # fitted curve through the fitting range and observed counts only afterwards.
  observed_start_x <- if (is.finite(switch_x)) {
    switch_x
  } else {
    fit_range_upper + 0.5 * bin_width
  }
  red_indices <- which(natural_x < observed_start_x)
  red_total_kmers <- if (length(red_indices)) {
    sum(natural_x[red_indices] * fitted_counts[red_indices])
  } else 0
  observed_after_switch <- observed[
    observed$x >= observed_start_x, , drop = FALSE
  ]
  grey_total_kmers <- if (nrow(observed_after_switch)) {
    sum(observed_after_switch$x * observed_after_switch$y)
  } else 0

  # Keep component-size reporting consistent with the same integration rule.
  # Before the switch use fitted component counts.  Between the switch and the
  # fitting boundary allocate observed grey counts by GMM posterior responsibility;
  # data beyond the fitting boundary remains residual mass.
  component_masses <- numeric(length(alpha))
  grey_fit <- observed[
    observed$x >= observed_start_x &
      observed$x <= fit_range_upper,
    , drop = FALSE
  ]
  if (nrow(grey_fit)) {
    grey_component_density <- vapply(seq_along(alpha), function(j) {
      alpha[j] * stats::dnorm(grey_fit$x, miu[j], sigma[j])
    }, numeric(nrow(grey_fit)))
    if (is.null(dim(grey_component_density))) {
      grey_component_density <- matrix(
        grey_component_density, ncol = length(alpha)
      )
    }
    grey_density_total <- rowSums(grey_component_density)
  }
  for (j in seq_along(alpha)) {
    fitted_component_counts <- alpha[j] * stats::dnorm(
      natural_x, miu[j], sigma[j]
    ) * total_count_scale * bin_width
    red_component_mass <- if (length(red_indices)) {
      sum(natural_x[red_indices] * fitted_component_counts[red_indices])
    } else 0
    grey_component_mass <- 0
    if (nrow(grey_fit)) {
      responsibility <- ifelse(
        grey_density_total > 0,
        grey_component_density[, j] / grey_density_total,
        0
      )
      responsibility[!is.finite(responsibility)] <- 0
      grey_component_mass <- sum(
        grey_fit$x * grey_fit$y * responsibility
      )
    }
    component_masses[j] <- red_component_mass + grey_component_mass
  }

  list(
    total_kmers = red_total_kmers + grey_total_kmers,
    red_total_kmers = red_total_kmers,
    grey_total_kmers = grey_total_kmers,
    switch_x = switch_x,
    agreement_tolerance = agreement_tolerance,
    required_agreement_bins = required_agreement_bins,
    component_masses = component_masses
  )
}

# Require a data-supported genomic mode before either the known-n or unknown-n
# workflow is allowed to fit.  A median-smoothed log spectrum makes a single
# count spike harmless; the accepted structure must rise and fall over multiple
# consecutive bins.  Very broad, deeply overlapping modes remain eligible when
# they form a sustained high-curvature shoulder rather than a sharp maximum.
.GenomeGMM_no_peak_message <- function() {
  paste(
    "you need to increase sequencing coverage to estimate genome size",
    "using the k-mer frequencies."
  )
}

.GenomeGMM_write_no_peak_report <- function(data_file, n = NULL) {
  file_stem <- tools::file_path_sans_ext(basename(data_file))
  mode_stem <- if (is.null(n)) {
    ""
  } else {
    paste0("_n", as.integer(n))
  }
  result_file <- paste0(
    "em_fitted_gmm_", file_stem, mode_stem, "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"), "_result.txt"
  )
  message <- .GenomeGMM_no_peak_message()
  writeLines(
    enc2utf8(c(
      "GenomeGMM ploidy and genome-size inference",
      paste("Sample:", file_stem),
      "",
      "Status: no reliable genomic peak detected",
      message
    )),
    con = result_file, useBytes = TRUE
  )
  list(
    status = "no_reliable_peak",
    message = message,
    data_file = data_file,
    n = n,
    pdf_file = NULL,
    result_file = result_file,
    fitted = FALSE
  )
}

.GenomeGMM_assert_reliable_peak <- function(data_file) {
  failure_message <- .GenomeGMM_no_peak_message()
  if (!file.exists(data_file)) return(invisible(TRUE))
  histogram <- tryCatch(
    utils::read.table(data_file, header = FALSE),
    error = function(error) NULL
  )
  if (is.null(histogram) || ncol(histogram) < 2L || nrow(histogram) < 9L) {
    stop(failure_message, call. = FALSE)
  }
  spectrum <- data.frame(
    x = suppressWarnings(as.numeric(histogram[, 1])),
    y = suppressWarnings(as.numeric(histogram[, 2]))
  )
  # Zero/negative multiplicity cannot be a genomic coverage peak.
  spectrum <- spectrum[
    is.finite(spectrum$x) & is.finite(spectrum$y) &
      spectrum$x > 0 & spectrum$y >= 0,
    , drop = FALSE
  ]
  if (nrow(spectrum) < 9L) stop(failure_message, call. = FALSE)
  spectrum <- stats::aggregate(y ~ x, data = spectrum, FUN = sum)
  spectrum <- spectrum[order(spectrum$x), , drop = FALSE]
  x_steps <- diff(spectrum$x)
  positive_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  if (!length(positive_steps)) stop(failure_message, call. = FALSE)
  bin_width <- stats::median(positive_steps)

  # Median smoothing on log counts rejects isolated vertical spikes without
  # imposing an absolute count or species-specific coverage threshold.
  smooth_log_count <- stats::runmed(
    log1p(spectrum$y), k = 3L, endrule = "median"
  )
  smooth_count <- expm1(smooth_log_count)
  slopes <- diff(smooth_log_count) / pmax(x_steps, bin_width)
  peak_indices <- which(
    slopes[-length(slopes)] > 0 & slopes[-1L] <= 0
  ) + 1L
  if (!length(peak_indices)) stop(failure_message, call. = FALSE)

  # Exclude the first low-frequency error bin when defining the signal scale.
  reference_indices <- which(
    spectrum$x >= min(spectrum$x) + bin_width
  )
  reference_height <- max(smooth_count[reference_indices], na.rm = TRUE)
  if (!is.finite(reference_height) || reference_height <= 0) {
    stop(failure_message, call. = FALSE)
  }

  supported_peak <- vapply(peak_indices, function(index) {
    left <- index - 1L
    while (left > 1L && slopes[left - 1L] > 0) left <- left - 1L
    right <- index
    while (right <= length(slopes) && slopes[right] <= 0) {
      right <- right + 1L
    }
    rising_bins <- index - left
    falling_bins <- right - index
    prominence <- (smooth_count[index] - smooth_count[left]) /
      max(smooth_count[index], 1)
    height_fraction <- smooth_count[index] / reference_height

    # Ordinary explicit/shoulder peaks retain at least 0.3% of the post-error
    # signal and rise by at least 2%. A far-right peak below that level must be
    # exceptionally broad and prominent; this keeps genuine high-coverage
    # modes while rejecting one- or two-bin fluctuations in a count-poor tail.
    ordinary_support <- prominence >= 0.02 && height_fraction >= 0.003
    broad_deep_support <- prominence >= 0.50 && rising_bins >= 10L
    rising_bins >= 2L && falling_bins >= 2L &&
      (ordinary_support || broad_deep_support)
  }, logical(1))

  if (!any(supported_peak)) stop(failure_message, call. = FALSE)
  invisible(TRUE)
}

# Refine a recovered prefix with a fixed first-component height and width.
# Later widths are nondecreasing, including the optional tail component.
# Smoothness is selected from data; there is no fixed adjacent-width ratio.
.GenomeGMM_refine_prefix_residual <- function(
    x, y, alpha, miu, sigma, count_scale, spacing, max_iterations = 250L) {
  baseline <- list(alpha = alpha, miu = miu, sigma = sigma, accepted = FALSE)
  order_miu <- order(miu)
  first <- order_miu[1L]
  later <- order_miu[-1L]
  valid <- is.finite(x) & is.finite(y) & y >= 0 & x >= miu[first]
  x <- x[valid]
  y <- y[valid]
  if (length(later) < 1L || length(x) < 8L || !is.finite(spacing) ||
      spacing <= 0 || !is.finite(count_scale) || count_scale <= 0) return(baseline)
  # Every candidate is scored on the same observed domain. Moving mu1 must
  # never remove difficult bins from that domain or move the other means.
  local_scale <- pmax(y, 0.05 * max(y), 1)
  first_region <- x <= miu[later[1L]]
  second_region <- x >= (miu[first] + miu[later[1L]]) / 2 &
    x <= if (length(later) > 1L) mean(miu[later[1:2]]) else max(x)
  tail_region <- x > max(x[second_region], miu[later[1L]])
  regions <- list(first_region, second_region, tail_region)
  regions <- regions[vapply(regions, any, logical(1))]
  prediction <- function(a, centres, s) {
    as.numeric(vapply(order_miu, function(j) {
      a[j] * stats::dnorm(x, centres[j], s[j]) * count_scale
    }, numeric(length(x))) %*% rep(1, length(order_miu)))
  }
  metrics <- function(p) {
    error <- (p - y) / local_scale
    regional <- vapply(regions, function(mask) mean(error[mask]^2), numeric(1))
    list(loss = mean(regional), regional = regional,
         global = mean(((p - y) / max(y, 1))^2),
         under = mean(pmax(-error[second_region], 0)^2),
         flank_excess = mean(pmax(error[first_region] - 0.07, 0)^2),
         over = sum(pmax(p[first_region] - y[first_region], 0)) /
           max(sum(y[first_region]), 1),
         peak = (sum(alpha * stats::dnorm(miu[first], miu, sigma)) * count_scale))
  }
  before <- metrics(prediction(alpha, miu, sigma))
  count <- length(later)
  upper_sigma <- 0.80 * spacing
  width_span <- upper_sigma - sigma[first]
  if (!is.finite(width_span) || width_span <= 0) return(baseline)
  projected <- cummax(pmin(upper_sigma, pmax(sigma[first], sigma[order_miu])))
  increments <- diff(projected) / width_span
  shift_limit <- min(0.10 * spacing, 0.45 * (miu[later[1]] - miu[first]),
                     0.45 * miu[first])
  initial <- c(log(pmax(alpha[later], 1e-10)), increments, 0)
  decode <- function(parameters) {
    a <- alpha
    s <- sigma
    centres <- miu
    a[later] <- exp(parameters[seq_len(count)])
    growth <- parameters[count + seq_len(count)]
    # Nonnegative cumulative increments enforce monotonicity. The shared
    # normalization only enforces the pre-existing absolute anti-collapse
    # bound, not a per-neighbour growth cap.
    s[later] <- sigma[first] + width_span * cumsum(growth) / max(1, sum(growth))
    centres[first] <- miu[first] + shift_limit * parameters[2L * count + 1L]
    list(alpha = a, miu = centres, sigma = s)
  }
  width_roughness <- function(s) {
    growth <- diff(log(s[order_miu]))
    mean(growth^2) + if (length(growth) > 1L) mean(diff(growth)^2) else 0
  }
  objective <- function(parameters, strength) {
    candidate <- decode(parameters)
    m <- metrics(prediction(candidate$alpha, candidate$miu, candidate$sigma))
    peak <- sum(candidate$alpha * stats::dnorm(
      candidate$miu[first], candidate$miu, candidate$sigma)) * count_scale
    peak_excess <- pmax(abs(peak / max(before$peak, 1) - 1) - 0.03, 0)
    # Regional relative errors prevent the tall recovered first peak from
    # concealing a broad deficit around the second peak or in the tail.
    m$loss + m$under + m$over^2 + 100 * peak_excess^2 +
      1000 * m$flank_excess +
      1000 * pmax(m$over - 0.90 * min(0.05, max(before$over, 0.005)), 0)^2 +
      1000 * sum(pmax(m$regional - before$regional, 0)^2) +
      strength * width_roughness(candidate$sigma)
  }
  candidates <- lapply(c(0.003, 0.03, 0.3), function(strength) {
    solution <- tryCatch(stats::optim(initial, objective, strength = strength,
      method = "L-BFGS-B",
      lower = c(rep(log(1e-10), count), rep(0, count), -1),
      upper = c(rep(log(max(4 * sum(alpha), 1)), count), rep(1, count), 1),
      control = list(maxit = as.integer(max_iterations))), error = function(e) NULL)
    if (is.null(solution) || !is.finite(solution$value)) return(NULL)
    candidate <- decode(solution$par)
    after <- metrics(prediction(candidate$alpha, candidate$miu, candidate$sigma))
    new_peak <- sum(candidate$alpha * stats::dnorm(
      candidate$miu[first], candidate$miu, candidate$sigma)) * count_scale
    candidate$accepted <- after$loss < before$loss * 0.995 &&
      after$global <= before$global + 1e-8 &&
      all(after$regional <= before$regional + 1e-6) &&
      after$under <= max(before$under, 0.05^2) + 1e-6 &&
      after$over <= max(before$over, 0.005) + 1e-6 &&
      abs(new_peak / max(before$peak, 1) - 1) <= 0.04
    candidate$before <- before
    candidate$after <- after
    candidate$peak_ratio <- new_peak / max(before$peak, 1)
    candidate$width_strength <- strength
    candidate$width_roughness <- width_roughness(candidate$sigma)
    candidate$first_mean_shift <- candidate$miu[first] - miu[first]
    candidate
  })
  eligible <- Filter(function(value) !is.null(value) && isTRUE(value$accepted), candidates)
  if (!length(eligible)) {
    baseline$before <- before
    baseline$after <- before
    baseline$first_mean_shift <- 0
    return(baseline)
  }
  losses <- vapply(eligible, function(value) value$after$loss, numeric(1))
  # Prefer the smoother width sequence when its curve error is within 2% of
  # the best eligible curve. Do not compare different penalty strengths by
  # their penalized objectives, which are on different scales.
  close <- which(losses <= min(losses) * 1.02 + 1e-8)
  roughness <- vapply(eligible[close], function(value) value$width_roughness, numeric(1))
  eligible[[close[which.min(roughness)]]]
}

# Adapt the mass of an inferred first component to its observed right flank.
# The contaminated low-coverage prefix is deliberately excluded: the first
# centre and width describe the recovered geometry, whereas alpha remains an
# estimated quantity and must be allowed to move up or down with the data.
.GenomeGMM_refine_estimated_first_height <- function(
    x, y, alpha, miu, sigma, count_scale, spacing,
    biological_count, max_iterations = 120L) {
  baseline <- list(
    alpha = alpha, miu = miu, sigma = sigma,
    accepted = FALSE, triggered = FALSE, reason = "not_evaluated"
  )
  if (!is.finite(count_scale) || count_scale <= 0 ||
      !is.finite(spacing) || spacing <= 0 || biological_count < 2L) {
    baseline$reason <- "insufficient_input"
    return(baseline)
  }
  component_order <- order(miu)
  biological_count <- min(as.integer(biological_count), length(component_order))
  biological <- component_order[seq_len(biological_count)]
  first <- biological[1L]
  second <- biological[2L]
  right_limit <- 0.5 * (miu[first] + miu[second])
  valid <- is.finite(x) & is.finite(y) & y >= 0 &
    x >= miu[first] - 0.05 * spacing & x <= right_limit
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 6L) {
    baseline$reason <- "too_few_first_flank_bins"
    return(baseline)
  }
  ordered <- order(x)
  x <- x[ordered]
  y <- y[ordered]
  observed <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  observed[!is.finite(observed)] <- y[!is.finite(observed)]
  observed <- pmax(observed, 0)
  observed_max <- max(observed, na.rm = TRUE)
  if (!is.finite(observed_max) || observed_max <= 0) {
    baseline$reason <- "empty_first_flank"
    return(baseline)
  }

  component_basis <- vapply(seq_along(alpha), function(j) {
    stats::dnorm(x, miu[j], sigma[j]) * count_scale
  }, numeric(length(x)))
  if (is.null(dim(component_basis))) {
    component_basis <- matrix(component_basis, ncol = 1L)
  }
  other_prediction <- as.numeric(
    component_basis[, -first, drop = FALSE] %*% alpha[-first]
  )
  first_basis <- component_basis[, first]
  alpha_reference <- max(alpha[first], 1e-10)
  local_scale <- pmax(observed, 0.08 * observed_max, 1)
  centre_weight <- exp(-0.5 * ((x - miu[first]) /
    max(0.28 * spacing, .Machine$double.eps))^2)
  point_weight <- 0.35 + 0.65 * centre_weight

  evaluate <- function(alpha_first) {
    prediction <- other_prediction + alpha_first * first_basis
    relative_error <- (prediction - observed) / local_scale
    over <- sum(pmax(prediction - observed, 0)) / max(sum(observed), 1)
    under <- sum(pmax(observed - prediction, 0)) / max(sum(observed), 1)
    over_runs <- rle(relative_error > 0.035)
    under_runs <- rle(relative_error < -0.035)
    longest <- function(runs) {
      if (any(runs$values)) max(runs$lengths[runs$values]) else 0L
    }
    first_height <- stats::approx(
      x, observed, xout = miu[first], rule = 2
    )$y
    first_prediction <- stats::approx(
      x, prediction, xout = miu[first], rule = 2
    )$y
    list(
      loss = sum(point_weight * relative_error^2) / sum(point_weight) +
        0.30 * over^2 + 0.20 * under^2,
      over_mse = sum(point_weight * pmax(relative_error, 0)^2) /
        sum(point_weight),
      under_mse = sum(point_weight * pmax(-relative_error, 0)^2) /
        sum(point_weight),
      nrmse = sqrt(sum(point_weight *
        ((prediction - observed) / observed_max)^2) / sum(point_weight)),
      over = over,
      under = under,
      longest_over_run = as.integer(longest(over_runs)),
      longest_under_run = as.integer(longest(under_runs)),
      centre_ratio = first_prediction / max(first_height, 1)
    )
  }

  before <- evaluate(alpha_reference)
  trigger <- before$longest_over_run >= 2L ||
    before$longest_under_run >= 2L ||
    before$over >= 0.015 || before$under >= 0.020 ||
    before$centre_ratio < 0.98 || before$centre_ratio > 1.02
  baseline$before <- before
  baseline$after <- before
  baseline$triggered <- isTRUE(trigger)
  if (!isTRUE(trigger)) {
    baseline$reason <- "first_height_already_supported"
    return(baseline)
  }

  # When the fitted flank is predominantly above the observations, excess is
  # more informative than a small compensating deficit at one or two bins.
  # Derive the penalty from the measured over/under balance rather than using
  # a sample-specific lowering factor.
  overshoot_priority <- if (before$over > before$under) {
    min(4, max(1, before$over / max(before$under, 0.005)))
  } else 0
  selection_score <- function(metric, log_scale = 0) {
    metric$loss + overshoot_priority * metric$over_mse +
      0.20 * metric$under_mse + 0.004 * log_scale^2
  }
  objective <- function(log_scale) {
    metric <- evaluate(alpha_reference * exp(log_scale))
    selection_score(metric, log_scale)
  }
  solution <- tryCatch(stats::optimize(
    objective, interval = log(c(0.50, 1.60)),
    tol = 1e-5, maximum = FALSE
  ), error = function(error) NULL)
  if (is.null(solution) || !is.finite(solution$objective)) {
    baseline$reason <- "optimization_failed"
    return(baseline)
  }
  alpha_candidate <- alpha
  alpha_candidate[first] <- alpha_reference * exp(solution$minimum)
  after <- evaluate(alpha_candidate[first])
  changed <- abs(solution$minimum) >= 0.005
  before_score <- selection_score(before)
  after_score <- selection_score(after, solution$minimum)
  accepted <- isTRUE(changed) &&
    after_score <= before_score * 0.995 &&
    after$nrmse <= before$nrmse * 1.05 + 1e-6 &&
    after$centre_ratio >= 0.86 && after$centre_ratio <= 1.08 &&
    after$over <= max(before$over + 0.003, 0.01) &&
    after$under <= max(before$under + 0.006, 0.025)
  list(
    alpha = if (accepted) alpha_candidate else alpha,
    miu = miu, sigma = sigma,
    accepted = accepted, triggered = TRUE,
    reason = if (accepted) "accepted" else "guard_rejected",
    before = before, after = after,
    overshoot_priority = overshoot_priority,
    score_before = before_score, score_after = after_score,
    first_component = first,
    alpha_ratio = alpha_candidate[first] / alpha_reference
  )
}

# Final joint refinement for a recovered first Gaussian. Alpha, centre, and
# width are all inferred parameters; optimise them together against the first
# peak and its transition into component 2. Biological constraints remain as
# bounds, but no one parameter is treated as the fitting target itself.
.GenomeGMM_refine_estimated_first_shape <- function(
    x, y, alpha, miu, sigma, count_scale, spacing,
    biological_count, max_iterations = 220L,
    use_pdf_reference = FALSE, final_selection = FALSE) {
  baseline <- list(
    alpha = alpha, miu = miu, sigma = sigma,
    accepted = FALSE, triggered = FALSE, reason = "not_evaluated"
  )
  if (!is.finite(count_scale) || count_scale <= 0 ||
      !is.finite(spacing) || spacing <= 0 || biological_count < 2L) {
    baseline$reason <- "insufficient_input"
    return(baseline)
  }
  component_order <- order(miu)
  biological_count <- min(as.integer(biological_count), length(component_order))
  biological <- component_order[seq_len(biological_count)]
  first <- biological[1L]
  second <- biological[2L]
  harmonic_depth <- stats::median(
    miu[biological] / seq_along(biological), na.rm = TRUE
  )
  if (!is.finite(harmonic_depth) || harmonic_depth <= 0) {
    harmonic_depth <- spacing
  }
  centre_reference <- if (isTRUE(final_selection)) {
    harmonic_depth
  } else {
    miu[first]
  }
  input_valid <- is.finite(x) & is.finite(y) & y >= 0
  input_x <- as.numeric(x[input_valid])
  input_y <- as.numeric(y[input_valid])
  input_order <- order(input_x)
  input_x <- input_x[input_order]
  input_y <- input_y[input_order]
  left_limit <- centre_reference - 0.08 * spacing
  right_limit <- 0.5 * (centre_reference + miu[second]) +
    0.10 * spacing
  if (isTRUE(use_pdf_reference) && length(input_x) > 3L) {
    reference_rows <- if (length(input_x) > 1L) {
      2:length(input_x)
    } else seq_along(input_x)
    reference_x <- input_x[reference_rows]
    reference_y <- input_y[reference_rows]
    reference_left <- max(left_limit, min(reference_x))
    reference_right <- min(right_limit, max(reference_x))
    reference_grid <- seq(
      reference_left, reference_right, length.out = 240L
    )
    # Evaluate the same cubic spline as the PDF directly on a dense local
    # grid.  Creating 1000 points across the entire histogram first can leave
    # only one or two points in an early recovered peak when the tail is long.
    grey_reference <- if (length(reference_x) > 3L) {
      stats::spline(reference_x, reference_y, xout = reference_grid)
    } else {
      list(
        x = reference_grid,
        y = stats::approx(
          reference_x, reference_y, xout = reference_grid, rule = 2
        )$y
      )
    }
    x <- grey_reference$x
    y <- pmax(grey_reference$y, 0)
  } else {
    x <- input_x
    y <- input_y
  }
  valid <- is.finite(x) & is.finite(y) & y >= 0 &
    x >= left_limit & x <= right_limit
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 8L) {
    baseline$reason <- "too_few_first_peak_bins"
    return(baseline)
  }
  ordered <- order(x)
  x <- x[ordered]
  y <- y[ordered]
  observed <- if (isTRUE(use_pdf_reference)) {
    y
  } else {
    filtered <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
    filtered[!is.finite(filtered)] <- y[!is.finite(filtered)]
    filtered
  }
  observed <- pmax(observed, 0)
  observed_max <- max(observed, na.rm = TRUE)
  if (!is.finite(observed_max) || observed_max <= 0) {
    baseline$reason <- "empty_first_peak"
    return(baseline)
  }
  local_scale <- pmax(observed, 0.08 * observed_max, 1)
  point_weight <- 0.35 + 0.65 * exp(-0.5 * (
    (x - centre_reference) /
      max(0.38 * spacing, .Machine$double.eps)
  )^2)
  if (isTRUE(final_selection)) {
    # Bins preceding a reconstructed peak may still contain the descending
    # sequencing-error prefix.  They remain in the objective, but the peak
    # top and reliable right shoulder carry more evidence for the inferred
    # Gaussian shape.
    point_weight <- point_weight * ifelse(
      x < centre_reference, 0.70, 1.15
    )
  }

  prediction <- function(alpha_value, miu_value, sigma_value) {
    curves <- vapply(seq_along(alpha_value), function(j) {
      alpha_value[j] * stats::dnorm(
        x, miu_value[j], sigma_value[j]
      ) * count_scale
    }, numeric(length(x)))
    if (is.null(dim(curves))) as.numeric(curves) else rowSums(curves)
  }
  metrics <- function(predicted) {
    relative_error <- (predicted - observed) / local_scale
    over <- sum(pmax(predicted - observed, 0)) / max(sum(observed), 1)
    under <- sum(pmax(observed - predicted, 0)) / max(sum(observed), 1)
    observed_slope <- diff(observed) / observed_max
    predicted_slope <- diff(predicted) / observed_max
    list(
      mse = sum(point_weight * relative_error^2) / sum(point_weight),
      over_mse = sum(point_weight * pmax(relative_error, 0)^2) /
        sum(point_weight),
      under_mse = sum(point_weight * pmax(-relative_error, 0)^2) /
        sum(point_weight),
      nrmse = sqrt(sum(point_weight *
        ((predicted - observed) / observed_max)^2) / sum(point_weight)),
      slope_mse = mean((predicted_slope - observed_slope)^2),
      over = over, under = under,
      max_over = max(relative_error),
      max_under = max(-relative_error)
    )
  }

  before <- metrics(prediction(alpha, miu, sigma))
  baseline$before <- before
  trigger <- before$over >= 0.008 || before$under >= 0.015 ||
    before$max_over >= 0.035 || before$max_under >= 0.05
  baseline$triggered <- isTRUE(trigger)
  if (!isTRUE(trigger)) {
    baseline$after <- before
    baseline$reason <- "first_shape_already_supported"
    return(baseline)
  }
  overshoot_priority <- if (before$over > before$under) {
    min(4, max(1, before$over / max(before$under, 0.005)))
  } else 0
  alpha_reference <- max(alpha[first], 1e-10)
  sigma_reference <- max(sigma[first], 1e-8)
  sigma_lower <- max(
    if (isTRUE(final_selection)) 0.32 * spacing else 0.35 * spacing,
    0.75 * sigma_reference
  )
  sigma_upper <- min(
    if (isTRUE(final_selection)) 0.50 * spacing else 0.40 * spacing,
    1.25 * sigma_reference,
    if (isTRUE(final_selection)) sigma[second] else Inf
  )
  if (!is.finite(sigma_lower) || !is.finite(sigma_upper) ||
      sigma_upper <= sigma_lower) {
    baseline$after <- before
    baseline$reason <- "no_admissible_first_width"
    return(baseline)
  }
  decode <- function(parameters) {
    alpha_value <- alpha
    miu_value <- miu
    sigma_value <- sigma
    alpha_value[first] <- alpha_reference * exp(parameters[1L])
    miu_value[first] <- if (isTRUE(final_selection)) {
      centre_reference + spacing * parameters[2L]
    } else {
      miu[first] + spacing * parameters[2L]
    }
    sigma_value[first] <- sigma_reference * exp(parameters[3L])
    list(alpha = alpha_value, miu = miu_value, sigma = sigma_value)
  }
  score <- function(metric, parameters = c(0, 0, 0)) {
    if (isTRUE(final_selection)) {
      metric$mse +
        1.25 * metric$over_mse +
        0.20 * metric$under_mse +
        0.20 * (metric$over - metric$under)^2 +
        8 * metric$slope_mse +
        0.001 * parameters[1L]^2 +
        0.06 * parameters[2L]^2 +
        0.003 * parameters[3L]^2
    } else {
      metric$mse + overshoot_priority * metric$over_mse +
        0.20 * metric$under_mse + 12 * metric$slope_mse +
        0.006 * parameters[1L]^2 +
        0.20 * parameters[2L]^2 +
        0.015 * parameters[3L]^2
    }
  }
  objective <- function(parameters) {
    candidate <- decode(parameters)
    score(metrics(prediction(
      candidate$alpha, candidate$miu, candidate$sigma
    )), parameters)
  }
  lower_bound <- c(
    log(0.50), -0.08, log(sigma_lower / sigma_reference)
  )
  upper_bound <- c(
    log(1.60), 0.08, log(sigma_upper / sigma_reference)
  )
  baseline_parameters <- c(
    0,
    if (isTRUE(final_selection)) {
      (miu[first] - centre_reference) / spacing
    } else 0,
    0
  )
  baseline_parameters <- pmin(
    upper_bound, pmax(lower_bound, baseline_parameters)
  )
  starting_points <- list(baseline_parameters)
  if (isTRUE(final_selection)) {
    # A recovered peak can begin on either side of its initial estimate.
    # Multiple deterministic starts prevent a one-direction local optimum
    # from deciding the final PDF.
    start_grid <- expand.grid(
      alpha_shift = log(c(0.88, 1.00, 1.12)),
      centre_shift = c(-0.05, 0, 0.05),
      sigma_shift = log(c(0.90, 1.00)),
      KEEP.OUT.ATTRS = FALSE
    )
    starting_points <- c(starting_points, lapply(
      seq_len(nrow(start_grid)), function(index) {
        pmin(upper_bound, pmax(lower_bound, as.numeric(unlist(
          start_grid[index, ], use.names = FALSE
        ))))
      }
    ))
  }
  solutions <- lapply(starting_points, function(start_value) {
    tryCatch(stats::optim(
      start_value, objective, method = "L-BFGS-B",
      lower = lower_bound,
      upper = upper_bound,
      control = list(maxit = as.integer(max_iterations), factr = 1e8)
    ), error = function(error) NULL)
  })
  valid_solutions <- solutions[vapply(
    solutions,
    function(value) !is.null(value) && is.finite(value$value),
    logical(1)
  )]
  solution <- if (length(valid_solutions)) {
    valid_solutions[[which.min(vapply(
      valid_solutions, function(value) value$value, numeric(1)
    ))]]
  } else NULL
  if (is.null(solution) || !is.finite(solution$value)) {
    baseline$after <- before
    baseline$reason <- "optimization_failed"
    return(baseline)
  }
  candidate <- decode(solution$par)
  after <- metrics(prediction(
    candidate$alpha, candidate$miu, candidate$sigma
  ))
  before_score <- score(before, baseline_parameters)
  after_score <- score(after, solution$par)
  changed <- any(abs(solution$par) >= 0.004)
  accepted <- if (isTRUE(final_selection)) {
    isTRUE(changed) &&
      after_score <= before_score * 0.995 &&
      after$mse <= before$mse * 1.05 + 1e-8 &&
      after$nrmse <= before$nrmse * 1.05 + 1e-8 &&
      after$over + after$under <=
        before$over + before$under + 0.001 &&
      after$over <= before$over + 1e-8 &&
      after$max_over <= before$max_over + 1e-8 &&
      after$under <= max(before$under + 0.010, 0.040)
  } else {
    isTRUE(changed) &&
      after_score <= before_score * 0.99 &&
      after$mse <= before$mse * 1.02 + 1e-7 &&
      after$nrmse <= before$nrmse * 1.03 + 1e-6 &&
      after$over <= max(before$over + 0.003, 0.010) &&
      after$under <= max(before$under + 0.012, 0.030) &&
      after$max_over <= max(before$max_over + 0.01, 0.05)
  }
  list(
    alpha = if (accepted) candidate$alpha else alpha,
    miu = if (accepted) candidate$miu else miu,
    sigma = if (accepted) candidate$sigma else sigma,
    accepted = accepted, triggered = TRUE,
    reason = if (accepted) "accepted" else "guard_rejected",
    before = before, after = after,
    overshoot_priority = overshoot_priority,
    score_before = before_score, score_after = after_score,
    alpha_ratio = candidate$alpha[first] / alpha_reference,
    first_mean_shift = candidate$miu[first] - miu[first],
    first_sigma_ratio = candidate$sigma[first] / sigma_reference
  )
}

# Rebalance only component weights after hidden-prefix recovery. Peak centres,
# widths, component count, and the selected ploidy remain exactly fixed. This
# targets cases where overlapping Gaussian tails make one or more total-curve
# peaks too high even though the individual peak geometry is reasonable.
.GenomeGMM_refine_prefix_weights <- function(
    x, y, alpha, miu, sigma, count_scale, spacing,
    biological_count, max_iterations = 220L) {
  baseline <- list(
    alpha = alpha, miu = miu, sigma = sigma,
    accepted = FALSE, triggered = FALSE, reason = "not_evaluated"
  )
  if (!is.finite(count_scale) || count_scale <= 0 ||
      !is.finite(spacing) || spacing <= 0 || biological_count < 2L) {
    baseline$reason <- "insufficient_input"
    return(baseline)
  }
  component_order <- order(miu)
  biological_count <- min(
    as.integer(biological_count), length(component_order)
  )
  biological <- component_order[seq_len(biological_count)]
  first <- biological[1L]
  valid <- is.finite(x) & is.finite(y) & y >= 0 & x >= miu[first]
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 12L) {
    baseline$reason <- "too_few_observed_bins"
    return(baseline)
  }
  ordered <- order(x)
  x <- x[ordered]
  y <- y[ordered]
  observed <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  observed[!is.finite(observed)] <- y[!is.finite(observed)]
  observed <- pmax(observed, 0)
  observed_max <- max(observed, na.rm = TRUE)
  if (!is.finite(observed_max) || observed_max <= 0) {
    baseline$reason <- "empty_observed_curve"
    return(baseline)
  }

  unit_curves <- vapply(seq_along(alpha), function(j) {
    stats::dnorm(x, miu[j], sigma[j]) * count_scale
  }, numeric(length(x)))
  if (is.null(dim(unit_curves))) {
    unit_curves <- matrix(unit_curves, ncol = 1L)
  }
  predict_curve <- function(alpha_value) {
    prediction <- as.numeric(unit_curves %*% alpha_value)
    prediction[!is.finite(prediction)] <- 0
    prediction
  }

  sorted_centres <- miu[component_order]
  boundaries <- if (length(sorted_centres) > 1L) {
    0.5 * (sorted_centres[-length(sorted_centres)] +
             sorted_centres[-1L])
  } else numeric(0)
  region_id <- findInterval(x, boundaries) + 1L
  region_ids <- sort(unique(region_id))
  height_weight <- 0.25 + 0.75 * sqrt(observed / observed_max)
  region_weight <- rep(0, length(x))
  for (region in region_ids) {
    rows <- region_id == region
    region_weight[rows] <- height_weight[rows] /
      max(sum(height_weight[rows]), .Machine$double.eps)
  }
  # Every harmonic interval carries equal influence, preventing the recovered
  # first peak from concealing a high third or later component.
  region_weight <- region_weight / max(length(region_ids), 1L)
  local_scale <- pmax(observed, 0.05 * observed_max, 1)

  peak_windows <- lapply(biological, function(component) {
    abs(x - miu[component]) <= 0.16 * spacing
  })
  valley_centres <- 0.5 * (
    miu[biological[-length(biological)]] + miu[biological[-1L]]
  )
  valley <- rep(FALSE, length(x))
  for (centre in valley_centres) {
    valley <- valley | abs(x - centre) <= 0.10 * spacing
  }
  interpeak_windows <- lapply(
    seq_len(length(biological) - 1L), function(position) {
      x >= miu[biological[position]] &
        x <= miu[biological[position + 1L]]
    }
  )

  metrics <- function(prediction) {
    scaled_error <- (prediction - observed) / local_scale
    regional_mse <- vapply(region_ids, function(region) {
      rows <- region_id == region
      mean(scaled_error[rows]^2)
    }, numeric(1))
    regional_over <- vapply(region_ids, function(region) {
      rows <- region_id == region
      sum(pmax(prediction[rows] - observed[rows], 0)) /
        max(sum(observed[rows]), 1)
    }, numeric(1))
    regional_under <- vapply(region_ids, function(region) {
      rows <- region_id == region
      sum(pmax(observed[rows] - prediction[rows], 0)) /
        max(sum(observed[rows]), 1)
    }, numeric(1))
    peak_over <- vapply(peak_windows, function(rows) {
      if (!any(rows)) return(0)
      sum(pmax(prediction[rows] - observed[rows], 0)) /
        max(sum(observed[rows]), 1)
    }, numeric(1))
    peak_error <- vapply(peak_windows, function(rows) {
      if (!any(rows)) return(0)
      mean(abs(prediction[rows] - observed[rows])) /
        max(max(observed[rows]), 1)
    }, numeric(1))
    first_height <- stats::approx(
      x, observed, xout = miu[first], rule = 2
    )$y
    first_prediction <- stats::approx(
      x, prediction, xout = miu[first], rule = 2
    )$y
    overshoot_flag <- scaled_error > 0.035
    runs <- rle(overshoot_flag)
    longest_over_run <- if (any(runs$values)) {
      max(runs$lengths[runs$values])
    } else 0L
    slope_diagnostics <- lapply(interpeak_windows, function(rows) {
      if (sum(rows) < 3L) {
        return(c(mse = 0, artificial_rise = 0))
      }
      observed_slope <- diff(observed[rows]) / observed_max
      predicted_slope <- diff(prediction[rows]) / observed_max
      observed_rise <- max(c(observed_slope, 0), na.rm = TRUE)
      predicted_rise <- max(c(predicted_slope, 0), na.rm = TRUE)
      c(
        mse = mean((predicted_slope - observed_slope)^2),
        # A small numerical allowance prevents raw-bin jitter from being
        # interpreted as a biologically meaningful extra shoulder.
        artificial_rise = pmax(
          predicted_rise - observed_rise - 0.001, 0
        )
      )
    })
    slope_diagnostics <- do.call(rbind, slope_diagnostics)
    list(
      balanced_mse = sum(region_weight * scaled_error^2),
      global_nrmse = sqrt(mean(((prediction - observed) /
        observed_max)^2)),
      region_mse = regional_mse,
      max_region_mse = max(regional_mse),
      region_over = regional_over,
      max_region_over = max(regional_over),
      max_region_under = max(regional_under),
      peak_over = peak_over,
      max_peak_over = max(peak_over),
      peak_error = peak_error,
      max_peak_error = max(peak_error),
      first_peak_error = peak_error[1L],
      first_centre_ratio = first_prediction / max(first_height, 1),
      valley_over = if (any(valley)) {
        sum(pmax(prediction[valley] - observed[valley], 0)) /
          max(sum(observed[valley]), 1)
      } else 0,
      under_area = sum(pmax(observed - prediction, 0)) /
        max(sum(observed), 1),
      over_area = sum(pmax(prediction - observed, 0)) /
        max(sum(observed), 1),
      area_error = abs(sum(prediction) - sum(observed)) /
        max(sum(observed), 1),
      longest_over_run = as.integer(longest_over_run),
      slope_mse = mean(slope_diagnostics[, "mse"]),
      max_artificial_rise = max(
        slope_diagnostics[, "artificial_rise"]
      )
    )
  }

  before <- metrics(predict_curve(alpha))
  # A candidate is needed only for a sustained or region-level excess. Small
  # point noise is ignored and the ordinary fit is returned unchanged.
  trigger <- before$max_region_over >= 0.02 &&
    before$longest_over_run >= 3L
  baseline$before <- before
  baseline$after <- before
  baseline$triggered <- isTRUE(trigger)
  if (!isTRUE(trigger)) {
    baseline$reason <- "no_sustained_regional_overshoot"
    return(baseline)
  }

  alpha_reference <- pmax(alpha, 1e-10)
  decode <- function(parameters) alpha_reference * exp(parameters)
  objective <- function(parameters) {
    alpha_value <- decode(parameters)
    metric <- metrics(predict_curve(alpha_value))
    ridge <- mean(parameters^2)
    metric$balanced_mse + 2.00 * metric$global_nrmse^2 +
      0.50 * metric$max_region_mse +
      0.50 * metric$max_peak_over^2 +
      1.00 * metric$max_region_under^2 +
      0.50 * metric$valley_over^2 +
      1.00 * metric$under_area^2 +
      1.00 * metric$area_error^2 +
      10.0 * metric$slope_mse +
      500 * metric$max_artificial_rise^2 + 0.006 * ridge
  }
  initial <- rep(0, length(alpha))
  solution <- tryCatch(stats::optim(
    initial, objective, method = "L-BFGS-B",
    lower = rep(log(0.35), length(alpha)),
    upper = rep(log(2.50), length(alpha)),
    control = list(maxit = as.integer(max_iterations), factr = 1e8)
  ), error = function(error) NULL)
  if (is.null(solution) || !is.finite(solution$value)) {
    baseline$reason <- "optimization_failed"
    return(baseline)
  }
  candidate_alpha <- decode(solution$par)
  after <- metrics(predict_curve(candidate_alpha))
  shares <- candidate_alpha / max(sum(candidate_alpha), 1e-12)
  changed <- any(abs(solution$par) >= 0.005)
  overshoot_improved <- after$max_region_over <=
    before$max_region_over * 0.90 ||
    after$max_peak_over <= before$max_peak_over * 0.90
  accepted <- isTRUE(changed) && isTRUE(overshoot_improved) &&
    after$balanced_mse <= before$balanced_mse * 0.985 &&
    after$global_nrmse <= before$global_nrmse * 1.002 + 1e-6 &&
    after$max_region_mse <= before$max_region_mse + 1e-6 &&
    after$first_peak_error <= before$first_peak_error + 0.003 &&
    after$first_centre_ratio >= 0.86 &&
    after$first_centre_ratio <= 1.08 &&
    after$max_region_over <= max(
      0.07, before$max_region_over + 0.025
    ) &&
    after$max_artificial_rise <= max(
      0.002, before$max_artificial_rise * 0.55
    ) &&
    after$slope_mse <= before$slope_mse * 1.05 + 1e-8 &&
    after$valley_over <= max(
      before$valley_over * 1.05, before$valley_over + 0.002
    ) &&
    after$max_region_under <= before$max_region_under + 0.01 &&
    after$under_area <= before$under_area + 0.008 &&
    after$area_error <= before$area_error + 0.005 &&
    all(shares[biological] >= 0.003)
  list(
    alpha = if (accepted) candidate_alpha else alpha,
    miu = miu, sigma = sigma, accepted = accepted,
    triggered = TRUE,
    reason = if (accepted) "accepted" else "guard_rejected",
    before = before, after = after,
    alpha_ratio = candidate_alpha / alpha_reference
  )
}

# Remove an artificial dip/rebound between two neighbouring fitted peaks.
# Unlike a rejection-only guard, this performs a local continuous optimisation
# of the two responsible components. The recovered first component and every
# peak centre remain fixed; widths are data-selected under only the existing
# absolute anti-collapse bound.
.GenomeGMM_refine_interpeak_shape <- function(
    x, y, alpha, miu, sigma, count_scale, spacing,
    biological_count, max_iterations = 260L) {
  baseline <- list(
    alpha = alpha, miu = miu, sigma = sigma,
    accepted = FALSE, triggered = FALSE, reason = "not_evaluated"
  )
  if (!is.finite(count_scale) || count_scale <= 0 ||
      !is.finite(spacing) || spacing <= 0 || biological_count < 2L) {
    baseline$reason <- "insufficient_input"
    return(baseline)
  }
  component_order <- order(miu)
  biological_count <- min(
    as.integer(biological_count), length(component_order)
  )
  biological <- component_order[seq_len(biological_count)]
  first <- biological[1L]
  valid <- is.finite(x) & is.finite(y) & y >= 0 & x >= miu[first]
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 12L) {
    baseline$reason <- "too_few_observed_bins"
    return(baseline)
  }
  ordered <- order(x)
  x <- x[ordered]
  y <- y[ordered]
  observed <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  observed[!is.finite(observed)] <- y[!is.finite(observed)]
  observed <- pmax(observed, 0)
  observed_max <- max(observed, na.rm = TRUE)
  if (!is.finite(observed_max) || observed_max <= 0) {
    baseline$reason <- "empty_observed_curve"
    return(baseline)
  }

  interpeak_windows <- lapply(
    seq_len(length(biological) - 1L), function(position) {
      x >= miu[biological[position]] &
        x <= miu[biological[position + 1L]]
    }
  )
  peak_windows <- lapply(biological, function(component) {
    abs(x - miu[component]) <= 0.16 * spacing
  })
  valley_centres <- 0.5 * (
    miu[biological[-length(biological)]] + miu[biological[-1L]]
  )
  valley <- rep(FALSE, length(x))
  for (centre in valley_centres) {
    valley <- valley | abs(x - centre) <= 0.10 * spacing
  }
  local_scale <- pmax(observed, 0.05 * observed_max, 1)

  predict_curve <- function(alpha_value, sigma_value) {
    curves <- vapply(seq_along(alpha_value), function(j) {
      alpha_value[j] * stats::dnorm(
        x, miu[j], sigma_value[j]
      ) * count_scale
    }, numeric(length(x)))
    prediction <- if (is.null(dim(curves))) curves else rowSums(curves)
    prediction[!is.finite(prediction)] <- 0
    prediction
  }
  metrics <- function(prediction) {
    scaled_error <- (prediction - observed) / local_scale
    interval_metrics <- lapply(interpeak_windows, function(rows) {
      if (sum(rows) < 3L) {
        return(c(mse = 0, artificial_rise = 0))
      }
      observed_slope <- diff(observed[rows]) / observed_max
      predicted_slope <- diff(prediction[rows]) / observed_max
      observed_rise <- max(c(observed_slope, 0), na.rm = TRUE)
      predicted_rise <- max(c(predicted_slope, 0), na.rm = TRUE)
      c(
        mse = mean(scaled_error[rows]^2),
        artificial_rise = pmax(
          predicted_rise - observed_rise - 0.001, 0
        )
      )
    })
    interval_metrics <- do.call(rbind, interval_metrics)
    peak_error <- vapply(peak_windows, function(rows) {
      if (!any(rows)) return(0)
      mean(abs(prediction[rows] - observed[rows])) /
        max(max(observed[rows]), 1)
    }, numeric(1))
    first_height <- stats::approx(
      x, observed, xout = miu[first], rule = 2
    )$y
    first_prediction <- stats::approx(
      x, prediction, xout = miu[first], rule = 2
    )$y
    list(
      interval_mse = interval_metrics[, "mse"],
      artificial_rise = interval_metrics[, "artificial_rise"],
      max_artificial_rise = max(
        interval_metrics[, "artificial_rise"]
      ),
      global_nrmse = sqrt(mean(((prediction - observed) /
        observed_max)^2)),
      under_area = sum(pmax(observed - prediction, 0)) /
        max(sum(observed), 1),
      over_area = sum(pmax(prediction - observed, 0)) /
        max(sum(observed), 1),
      area_error = abs(sum(prediction) - sum(observed)) /
        max(sum(observed), 1),
      peak_error = peak_error,
      max_peak_error = max(peak_error),
      first_centre_ratio = first_prediction / max(first_height, 1),
      valley_over = if (any(valley)) {
        sum(pmax(prediction[valley] - observed[valley], 0)) /
          max(sum(observed[valley]), 1)
      } else 0
    )
  }

  before <- metrics(predict_curve(alpha, sigma))
  target_interval <- which.max(before$artificial_rise)
  baseline$before <- before
  baseline$after <- before
  baseline$target_interval <- target_interval
  if (!length(target_interval) ||
      before$artificial_rise[target_interval] <= 0.0015) {
    baseline$reason <- "no_artificial_interpeak_rebound"
    return(baseline)
  }
  baseline$triggered <- TRUE

  pair <- biological[c(target_interval, target_interval + 1L)]
  # The recovered first peak height is inferred rather than observed. Let its
  # weight participate in the local curve optimisation, while keeping its
  # centre and width fixed. Only the two components bordering the artificial
  # rebound may change width.
  alpha_adjustable <- unique(c(first, pair))
  width_adjustable <- setdiff(pair, first)
  if (!length(alpha_adjustable) || !length(width_adjustable)) {
    baseline$reason <- "no_adjustable_component"
    return(baseline)
  }
  alpha_reference <- pmax(alpha[alpha_adjustable], 1e-10)
  sigma_reference <- pmax(sigma[width_adjustable], 1e-8)
  width_lower <- pmax(0.25 * spacing, 0.65 * sigma_reference)
  width_upper <- pmin(0.80 * spacing, 1.50 * sigma_reference)
  if (any(!is.finite(width_lower)) || any(!is.finite(width_upper)) ||
      any(width_upper <= width_lower)) {
    baseline$reason <- "no_admissible_local_width"
    return(baseline)
  }

  decode <- function(parameters) {
    alpha_count <- length(alpha_adjustable)
    width_count <- length(width_adjustable)
    alpha_value <- alpha
    sigma_value <- sigma
    alpha_value[alpha_adjustable] <- alpha_reference * exp(
      parameters[seq_len(alpha_count)]
    )
    sigma_value[width_adjustable] <- sigma_reference * exp(
      parameters[alpha_count + seq_len(width_count)]
    )
    list(alpha = alpha_value, sigma = sigma_value)
  }
  objective <- function(parameters) {
    candidate <- decode(parameters)
    metric <- metrics(predict_curve(candidate$alpha, candidate$sigma))
    local_mse <- metric$interval_mse[target_interval]
    local_rise <- metric$artificial_rise[target_interval]
    width_change <- mean(parameters[
      length(alpha_adjustable) + seq_along(width_adjustable)
    ]^2)
    weight_change <- mean(parameters[seq_along(alpha_adjustable)]^2)
    local_mse + 2.0 * metric$global_nrmse^2 +
      800 * local_rise^2 +
      300 * metric$max_artificial_rise^2 +
      0.70 * metric$max_peak_error^2 +
      0.50 * metric$under_area^2 +
      0.40 * metric$over_area^2 +
      0.50 * metric$area_error^2 +
      0.40 * metric$valley_over^2 +
      0.010 * width_change + 0.006 * weight_change
  }
  alpha_count <- length(alpha_adjustable)
  width_count <- length(width_adjustable)
  solution <- tryCatch(stats::optim(
    rep(0, alpha_count + width_count), objective, method = "L-BFGS-B",
    lower = c(
      rep(log(0.35), alpha_count), log(width_lower / sigma_reference)
    ),
    upper = c(
      rep(log(2.50), alpha_count), log(width_upper / sigma_reference)
    ),
    control = list(maxit = as.integer(max_iterations), factr = 1e8)
  ), error = function(error) NULL)
  if (is.null(solution) || !is.finite(solution$value)) {
    baseline$reason <- "optimization_failed"
    return(baseline)
  }
  candidate <- decode(solution$par)
  after <- metrics(predict_curve(candidate$alpha, candidate$sigma))
  shares <- candidate$alpha / max(sum(candidate$alpha), 1e-12)
  local_before <- before$interval_mse[target_interval]
  local_after <- after$interval_mse[target_interval]
  accepted <- after$artificial_rise[target_interval] <= max(
      0.0015, before$artificial_rise[target_interval] * 0.45
    ) &&
    local_after <= local_before * 1.03 + 1e-7 &&
    after$global_nrmse <= before$global_nrmse * 1.015 + 1e-6 &&
    after$under_area <= before$under_area + 0.008 &&
    after$over_area <= before$over_area + 0.005 &&
    after$area_error <= before$area_error + 0.008 &&
    after$max_peak_error <= before$max_peak_error + 0.004 &&
    after$first_centre_ratio >= 0.86 &&
    after$first_centre_ratio <= 1.08 &&
    after$valley_over <= max(
      0.008, before$valley_over + 0.002
    ) &&
    all(shares[biological] >= 0.003)
  list(
    alpha = if (accepted) candidate$alpha else alpha,
    miu = miu,
    sigma = if (accepted) candidate$sigma else sigma,
    accepted = accepted, triggered = TRUE,
    reason = if (accepted) "accepted" else "guard_rejected",
    before = before, after = after,
    target_interval = target_interval,
    adjusted_alpha_components = alpha_adjustable,
    adjusted_width_components = width_adjustable,
    alpha_ratio = candidate$alpha[alpha_adjustable] /
      alpha[alpha_adjustable],
    sigma_ratio = candidate$sigma[width_adjustable] /
      sigma[width_adjustable]
  )
}

# After ploidy, component count and the optional tail decision have all been
# selected, make one guarded refinement of component area, centre and width
# against the displayed observed curve. This is deliberately outside candidate
# selection, so it cannot change n or the n+1 decision. A common-width winner is
# a stable starting point rather than a requirement on the final curve.
.GenomeGMM_refine_selected_shape <- function(
    x, y, alpha, miu, sigma, total_samples_fit,
    first_peak_spacing, sigma_reliability = NULL,
    biological_count = length(alpha),
    first_sigma_fraction_range = NULL,
    max_iterations = 260L) {
  baseline <- list(
    alpha = alpha, miu = miu, sigma = sigma,
    accepted = FALSE, triggered = FALSE, reason = "not_evaluated"
  )
  valid <- is.finite(x) & is.finite(y) & y >= 0
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 10L || !is.finite(total_samples_fit) ||
      total_samples_fit <= 0 || !is.finite(first_peak_spacing) ||
      first_peak_spacing <= 0 || any(!is.finite(alpha)) ||
      any(alpha <= 0) || any(!is.finite(miu)) || any(!is.finite(sigma)) ||
      any(sigma <= 0)) {
    baseline$reason <- "insufficient_input"
    return(baseline)
  }
  ordered_rows <- order(x)
  x <- x[ordered_rows]
  y <- y[ordered_rows]
  x_steps <- diff(sort(unique(x)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  count_scale <- total_samples_fit * bin_width
  observed <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  observed[!is.finite(observed)] <- y[!is.finite(observed)]
  observed <- pmax(observed, 0)
  observed_max <- max(observed, na.rm = TRUE)
  if (!is.finite(observed_max) || observed_max <= 0) {
    baseline$reason <- "empty_observed_curve"
    return(baseline)
  }

  component_order <- order(miu)
  biological_count <- min(
    max(1L, as.integer(biological_count)), length(component_order)
  )
  biological <- component_order[seq_len(biological_count)]
  reliability <- as.numeric(sigma_reliability)
  if (length(reliability) != length(alpha)) {
    reliability <- rep(0, length(alpha))
  }
  reliability <- pmin(1, pmax(0, reliability))
  spacing <- as.numeric(first_peak_spacing)
  biological_ranks <- seq_len(length(biological))
  nuisance <- setdiff(seq_along(alpha), biological)
  depth_ratios <- miu[biological] / biological_ranks
  robust_depth <- stats::median(depth_ratios[is.finite(depth_ratios)])
  if (!is.finite(robust_depth) || robust_depth <= 0) robust_depth <- spacing
  # The common harmonic spacing itself may move by at most three per cent.
  # Individual biological peaks are then allowed only a small residual around
  # their assigned integer slot. This can repair a collapsed weak component,
  # but it cannot move a peak freely to an unrelated part of the histogram.
  depth_centre <- min(max(robust_depth, 0.97 * spacing), 1.03 * spacing)
  depth_log_span <- log(1.03)
  biological_offset_limit <- spacing *
    (0.04 + 0.06 * (1 - reliability[biological]))
  nuisance_offset_limit <- rep(0.18 * spacing, length(nuisance))

  biological_validity <- function(alpha_value, miu_value, sigma_value,
                                  depth_value) {
    centres <- miu_value[biological]
    widths <- sigma_value[biological]
    shares <- alpha_value / max(sum(alpha_value), 1e-12)
    centre_residual <- abs(
      centres - biological_ranks * depth_value
    ) / max(depth_value, 1e-12)
    gap_ratio <- if (length(centres) > 1L) {
      diff(centres) / max(depth_value, 1e-12)
    } else 1
    width_ratio <- if (length(widths) > 1L) {
      widths[-1L] / pmax(widths[-length(widths)], 1e-12)
    } else 1
    peak_strength <- alpha_value[biological] / pmax(widths, 1e-12)
    strength_ratio <- if (length(peak_strength) > 1L) {
      pmax(
        peak_strength[-1L] /
          pmax(peak_strength[-length(peak_strength)], 1e-12),
        peak_strength[-length(peak_strength)] /
          pmax(peak_strength[-1L], 1e-12)
      )
    } else 1
    interior_support <- if (length(peak_strength) > 2L) {
      peak_strength[2:(length(peak_strength) - 1L)] /
        pmax(
          pmin(
            peak_strength[1:(length(peak_strength) - 2L)],
            peak_strength[3:length(peak_strength)]
          ),
          1e-12
        )
    } else 1
    biological_share_floor <- max(
      0.005, 0.06 / max(length(biological), 1L)
    )
    geometry_valid <- all(is.finite(centre_residual)) &&
      all(centre_residual <= 0.15) &&
      all(gap_ratio >= 0.65 & gap_ratio <= 1.35) &&
      all(diff(widths) >= -1e-8) &&
      all(width_ratio <= 1.08 + 1e-8) &&
      all(widths <= 0.80 * depth_value + 1e-8)
    height_valid <- all(shares[biological] >= biological_share_floor) &&
      all(strength_ratio <= 5.0 + 1e-8) &&
      all(interior_support >= 0.25 - 1e-8)
    list(
      valid = geometry_valid && height_valid,
      geometry_valid = geometry_valid,
      height_valid = height_valid,
      centre_residual = centre_residual,
      gap_ratio = gap_ratio,
      width_ratio = width_ratio,
      peak_strength = peak_strength,
      strength_ratio = strength_ratio,
      interior_support = interior_support,
      biological_share_floor = biological_share_floor
    )
  }
  baseline_depth <- stats::median(depth_ratios, na.rm = TRUE)
  baseline_biology <- biological_validity(
    alpha, miu, sigma, baseline_depth
  )
  local_scale <- pmax(observed, 0.05 * observed_max, 1)

  sorted_centres <- miu[component_order]
  boundaries <- if (length(sorted_centres) > 1L) {
    0.5 * (sorted_centres[-length(sorted_centres)] +
      sorted_centres[-1L])
  } else numeric(0)
  region_id <- findInterval(x, boundaries) + 1L
  region_ids <- sort(unique(region_id))
  point_weight <- 0.25 + 0.75 * sqrt(observed / observed_max)

  predict_curve <- function(alpha_value, miu_value, sigma_value) {
    component_curves <- vapply(seq_along(alpha), function(j) {
      alpha_value[j] * stats::dnorm(
        x, miu_value[j], sigma_value[j]
      ) * count_scale
    }, numeric(length(x)))
    prediction <- if (is.null(dim(component_curves))) {
      as.numeric(component_curves)
    } else rowSums(component_curves)
    prediction[!is.finite(prediction)] <- 0
    prediction
  }
  metrics <- function(prediction) {
    scaled_global <- (prediction - observed) / observed_max
    scaled_local <- (prediction - observed) / local_scale
    regional <- vapply(region_ids, function(region) {
      rows <- region_id == region
      mean(scaled_local[rows]^2)
    }, numeric(1))
    list(
      nrmse = sqrt(sum(point_weight * scaled_global^2) /
        sum(point_weight)),
      balanced_mse = mean(regional),
      max_region_mse = max(regional),
      regional_mse = regional,
      regional_over = vapply(region_ids, function(region) {
        rows <- region_id == region
        sum(pmax(prediction[rows] - observed[rows], 0)) /
          max(sum(observed[rows]), 1)
      }, numeric(1)),
      over_area = sum(pmax(prediction - observed, 0)) /
        max(sum(observed), 1),
      under_area = sum(pmax(observed - prediction, 0)) /
        max(sum(observed), 1),
      slope_mse = mean((diff(prediction) - diff(observed))^2) /
        max(observed_max^2, 1)
    )
  }
  before <- metrics(predict_curve(alpha, miu, sigma))
  baseline$before <- before
  baseline$after <- before
  baseline$triggered <- TRUE

  lower_sigma <- pmax(0.03 * spacing, 0.55 * sigma)
  upper_sigma <- pmin(0.90 * spacing, 1.60 * sigma)
  # Always include the selected solution in the search interval. These are
  # numerical trust bounds, not an adjacent-width biological ratio.
  lower_sigma <- pmin(lower_sigma, sigma)
  upper_sigma <- pmax(upper_sigma, sigma)
  if (!is.null(first_sigma_fraction_range) &&
      length(first_sigma_fraction_range) == 2L &&
      all(is.finite(first_sigma_fraction_range))) {
    first <- biological[1L]
    first_range <- sort(as.numeric(first_sigma_fraction_range)) * spacing
    lower_sigma[first] <- max(lower_sigma[first], first_range[1L])
    upper_sigma[first] <- min(upper_sigma[first], first_range[2L])
    lower_sigma[first] <- min(lower_sigma[first], sigma[first])
    upper_sigma[first] <- max(upper_sigma[first], sigma[first])
  }
  if (any(!is.finite(lower_sigma)) || any(!is.finite(upper_sigma)) ||
      any(upper_sigma <= lower_sigma)) {
    baseline$reason <- "no_admissible_width_interval"
    return(baseline)
  }

  width_penalty <- function(sigma_value) {
    log_width <- log(sigma_value[biological])
    rel <- reliability[biological]
    if (length(log_width) <= 1L) return(0)
    log_rank <- log(seq_along(log_width))
    trend <- stats::lm.wfit(
      cbind(1, log_rank), log_width, w = 0.15 + 0.85 * rel
    )
    trend_target <- as.numeric(cbind(1, log_rank) %*% trend$coefficients)
    steps <- diff(log_width)
    pair_rel <- pmin(rel[-1L], rel[-length(rel)])
    penalty <- mean((1 - rel) * (log_width - trend_target)^2) +
      mean((1 - pair_rel) * steps^2)
    if (length(steps) > 1L) {
      penalty <- penalty + 0.50 * mean(diff(steps)^2)
    }
    penalty + 0.25 * mean(pmax(sigma_value[biological] /
      spacing - 0.60, 0)^2)
  }
  height_continuity_penalty <- function(alpha_value, sigma_value) {
    strength <- alpha_value[biological] /
      pmax(sigma_value[biological], 1e-12)
    if (length(strength) <= 1L) return(0)
    log_jump <- abs(diff(log(pmax(strength, 1e-12))))
    penalty <- mean(pmax(log_jump - log(3), 0)^2)
    if (length(strength) > 2L) {
      interior_ratio <- strength[2:(length(strength) - 1L)] /
        pmax(
          pmin(
            strength[1:(length(strength) - 2L)],
            strength[3:length(strength)]
          ),
          1e-12
        )
      penalty <- penalty + mean(pmax(0.35 - interior_ratio, 0)^2)
    }
    penalty
  }
  mean_reliability <- mean(reliability[biological], na.rm = TRUE)
  if (!is.finite(mean_reliability)) mean_reliability <- 0
  smoothness_strength <- 0.004 + 0.026 * (1 - mean_reliability)

  # alpha remains data-driven, but peak centres and widths are decoded through
  # biologically admissible coordinates. Biological centres stay near j*d;
  # biological widths are nondecreasing and each increase is modest. Weak
  # peaks receive more flexibility than clear peaks without becoming free.
  alpha_log_span <- if (isTRUE(baseline_biology$geometry_valid)) {
    0.12 + 0.23 * (1 - reliability)
  } else {
    # A component whose old centre/width was biologically invalid must have
    # its mass re-estimated rather than remaining trapped near a collapsed
    # alpha. Reliable peaks still receive tighter bounds than weak peaks.
    0.35 + 0.75 * (1 - reliability)
  }
  if (length(nuisance)) {
    alpha_log_span[nuisance] <- pmax(alpha_log_span[nuisance], 1.10)
  }
  alpha_reference <- pmax(alpha, 1e-10)
  pair_reliability <- if (length(biological) > 1L) {
    pmin(
      reliability[biological[-1L]],
      reliability[biological[-length(biological)]]
    )
  } else numeric(0)
  maximum_width_growth <- 1.04 + 0.04 * (1 - pair_reliability)
  mean_biological_reliability <- mean(
    reliability[biological], na.rm = TRUE
  )
  if (!is.finite(mean_biological_reliability)) {
    mean_biological_reliability <- 0
  }
  maximum_total_width_ratio <-
    1.50 + 0.20 * (1 - mean_biological_reliability)

  decode <- function(parameters) {
    component_count <- length(alpha)
    alpha_value <- exp(parameters[seq_len(component_count)])
    depth_unit <- parameters[component_count + 1L]
    depth_value <- depth_centre * exp(depth_log_span * depth_unit)
    mean_units <- parameters[component_count + 1L +
      seq_len(component_count)]
    raw_sigma <- exp(parameters[2L * component_count + 1L +
      seq_len(component_count)])
    miu_value <- miu
    miu_value[biological] <- biological_ranks * depth_value +
      biological_offset_limit * mean_units[biological]
    if (length(nuisance)) {
      nuisance_ranks <- biological_count + seq_along(nuisance)
      miu_value[nuisance] <- nuisance_ranks * depth_value +
        nuisance_offset_limit * mean_units[nuisance]
    }

    sigma_value <- raw_sigma
    first_biological <- biological[1L]
    sigma_value[first_biological] <- min(
      max(raw_sigma[first_biological], lower_sigma[first_biological]),
      upper_sigma[first_biological], 0.80 * depth_value
    )
    if (length(biological) > 1L) {
      for (position in 2:length(biological)) {
        previous <- biological[position - 1L]
        current <- biological[position]
        progression <- (position - 1L) /
          max(length(biological) - 1L, 1L)
        cumulative_width_limit <- sigma_value[first_biological] *
          (1 + (maximum_total_width_ratio - 1) * progression)
        sigma_value[current] <- min(
          max(raw_sigma[current], sigma_value[previous]),
          sigma_value[previous] * maximum_width_growth[position - 1L],
          cumulative_width_limit,
          0.80 * depth_value
        )
      }
    }
    if (length(nuisance)) {
      reference_width <- sigma_value[biological[length(biological)]]
      for (current in nuisance) {
        sigma_value[current] <- min(
          max(raw_sigma[current], reference_width),
          1.20 * reference_width, 0.90 * depth_value
        )
        reference_width <- sigma_value[current]
      }
    }
    list(
      alpha = alpha_value,
      miu = miu_value,
      sigma = sigma_value,
      depth = depth_value,
      depth_unit = depth_unit,
      mean_units = mean_units
    )
  }
  harmonic_diagnostics <- function(miu_value) {
    ranks <- seq_len(length(biological))
    centres <- miu_value[biological]
    harmonic_weight <- 0.20 + 0.80 * reliability[biological]
    depth <- sum(harmonic_weight * centres / ranks) /
      sum(harmonic_weight)
    residual <- (centres - ranks * depth) / spacing
    list(depth = depth, residual = residual,
         error = sqrt(mean(residual^2)))
  }
  baseline_harmonic <- harmonic_diagnostics(miu)

  objective <- function(parameters) {
    candidate <- decode(parameters)
    metric <- metrics(predict_curve(
      candidate$alpha, candidate$miu, candidate$sigma
    ))
    alpha_change <- mean((
      log(candidate$alpha / alpha_reference) / alpha_log_span
    )^2)
    mean_change <- mean(
      (0.35 + 0.65 * reliability) * candidate$mean_units^2
    )
    biological_mean_units <- candidate$mean_units[biological]
    mean_roughness <- if (length(biological_mean_units) > 1L) {
      mean(diff(biological_mean_units)^2)
    } else 0
    0.60 * metric$nrmse^2 +
      0.32 * metric$balanced_mse +
      0.08 * metric$max_region_mse +
      0.10 * metric$slope_mse +
      smoothness_strength * width_penalty(candidate$sigma) +
      0.03 * height_continuity_penalty(
        candidate$alpha, candidate$sigma
      ) +
      0.002 * alpha_change + 0.003 * mean_change +
      0.006 * candidate$depth_unit^2 + 0.010 * mean_roughness
  }
  component_count <- length(alpha)
  initial_mean_units <- numeric(component_count)
  initial_mean_units[biological] <- pmax(-1, pmin(1,
    (miu[biological] - biological_ranks * depth_centre) /
      biological_offset_limit
  ))
  if (length(nuisance)) {
    nuisance_ranks <- biological_count + seq_along(nuisance)
    initial_mean_units[nuisance] <- pmax(-1, pmin(1,
      (miu[nuisance] - nuisance_ranks * depth_centre) /
        nuisance_offset_limit
    ))
  }
  initial <- c(
    log(alpha_reference), 0, initial_mean_units, log(sigma)
  )
  lower <- c(
    log(alpha_reference) - alpha_log_span,
    -1, rep(-1, component_count), log(lower_sigma)
  )
  upper <- c(
    log(alpha_reference) + alpha_log_span,
    1, rep(1, component_count), log(upper_sigma)
  )
  solution <- tryCatch(stats::optim(
    initial, objective, method = "L-BFGS-B",
    lower = lower, upper = upper,
    control = list(maxit = as.integer(max_iterations), factr = 1e8)
  ), error = function(error) NULL)
  if (is.null(solution) || !is.finite(solution$value)) {
    baseline$reason <- "optimization_failed"
    return(baseline)
  }
  candidate <- decode(solution$par)
  after <- metrics(predict_curve(
    candidate$alpha, candidate$miu, candidate$sigma
  ))
  shares <- candidate$alpha / max(sum(candidate$alpha), 1e-12)
  nuisance <- setdiff(seq_along(alpha), biological)
  nuisance_share_before <- if (length(nuisance)) {
    sum(alpha[nuisance]) / max(sum(alpha), 1e-12)
  } else 0
  nuisance_share_after <- if (length(nuisance)) sum(shares[nuisance]) else 0
  candidate_harmonic <- harmonic_diagnostics(candidate$miu)
  ordered_candidate_miu <- candidate$miu[component_order]
  candidate_biology <- biological_validity(
    candidate$alpha, candidate$miu, candidate$sigma, candidate$depth
  )
  nuisance_position_ok <- !length(nuisance) || all(
    candidate$miu[nuisance] >=
      max(candidate$miu[biological]) + 0.55 * candidate$depth
  )
  changed <- any(abs(log(candidate$alpha / alpha_reference)) >= 0.005) ||
    any(abs(candidate$miu - miu) >= 0.05 * bin_width) ||
    any(abs(log(candidate$sigma / sigma)) >= 0.005)
  error_improved <- after$nrmse <= before$nrmse * 0.997 ||
    (after$balanced_mse <= before$balanced_mse * 0.985 &&
      after$nrmse <= before$nrmse * 1.001)
  regional_guard <- all(
    after$regional_mse <= before$regional_mse * 1.08 + 2e-4
  ) && all(
    after$regional_over <= pmax(
      0.02, before$regional_over + 0.005
    )
  )
  common_biology_guard <- isTRUE(candidate_biology$geometry_valid) &&
    all(diff(ordered_candidate_miu) > max(bin_width, 0.15 * spacing)) &&
    all(shares[biological] >= 0.003) &&
    isTRUE(nuisance_position_ok) &&
    nuisance_share_after <= max(0.35, nuisance_share_before + 0.05)
  standard_acceptance <- isTRUE(error_improved) &&
    isTRUE(regional_guard) &&
    after$max_region_mse <= before$max_region_mse * 1.05 + 2e-4 &&
    after$over_area <= before$over_area + 0.003 &&
    after$under_area <= before$under_area + 0.003 &&
    after$slope_mse <= before$slope_mse * 1.03 + 1e-8
  # If the selected EM solution already violates the harmonic geometry, it is
  # not a valid fallback. Permit a constrained repair at the same n even when
  # its numerical error is only comparable, but never accept a large visible
  # regression merely to satisfy the constraints.
  biological_repair_acceptance <-
    !isTRUE(baseline_biology$geometry_valid) &&
    after$nrmse <= max(before$nrmse * 1.30, before$nrmse + 0.005) &&
    after$over_area <= before$over_area + 0.02 &&
    after$under_area <= before$under_area + 0.02 &&
    after$slope_mse <= max(
      before$slope_mse * 2.0,
      before$slope_mse + 5e-5
    )
  accepted <- isTRUE(changed) && isTRUE(common_biology_guard) &&
    (isTRUE(standard_acceptance) ||
       isTRUE(biological_repair_acceptance))
  list(
    alpha = if (accepted) candidate$alpha else alpha,
    miu = if (accepted) candidate$miu else miu,
    sigma = if (accepted) candidate$sigma else sigma,
    accepted = accepted, triggered = TRUE,
    reason = if (accepted &&
                 !isTRUE(baseline_biology$geometry_valid)) {
      "accepted_biological_repair"
    } else if (accepted) {
      "accepted"
    } else {
      "guard_rejected"
    },
    before = before, after = after,
    smoothness_strength = smoothness_strength,
    alpha_ratio = candidate$alpha / alpha_reference,
    miu_shift = candidate$miu - miu,
    sigma_ratio = candidate$sigma / sigma,
    baseline_harmonic_error = baseline_harmonic$error,
    candidate_harmonic_error = candidate_harmonic$error,
    baseline_biological_valid = isTRUE(baseline_biology$valid),
    candidate_biological_valid = isTRUE(candidate_biology$valid),
    baseline_geometry_valid = isTRUE(baseline_biology$geometry_valid),
    candidate_geometry_valid = isTRUE(candidate_biology$geometry_valid),
    fitted_spacing = candidate$depth
  )
}

# Final first-peak amplitude guard.  Peak discovery and prefix recovery decide
# the centre and width; this guard changes neither.  It only corrects the first
# component mass when the final red curve demonstrably misses the observed grey
# curve at that centre.  Fits already agreeing at the first peak are untouched.
.GenomeGMM_refine_estimated_first_height_final <- function(
    x, y, alpha, miu, sigma, total_samples_fit,
    biological_count = length(alpha), apply_guard = TRUE) {
  baseline <- list(
    alpha = alpha, accepted = FALSE, triggered = FALSE,
    reason = "guard_disabled"
  )
  if (!isTRUE(apply_guard)) return(baseline)
  valid <- is.finite(x) & is.finite(y) & y >= 0
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 8L || !is.finite(total_samples_fit) ||
      total_samples_fit <= 0 || any(!is.finite(alpha)) ||
      any(alpha <= 0) || any(!is.finite(miu)) ||
      any(!is.finite(sigma)) || any(sigma <= 0)) {
    baseline$reason <- "insufficient_input"
    return(baseline)
  }
  rows <- order(x)
  x <- x[rows]
  y <- y[rows]
  x_steps <- diff(sort(unique(x)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  observed <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  observed[!is.finite(observed)] <- y[!is.finite(observed)]
  observed <- pmax(observed, 0)
  observed_max <- max(observed, na.rm = TRUE)
  if (!is.finite(observed_max) || observed_max <= 0) {
    baseline$reason <- "empty_observed_curve"
    return(baseline)
  }

  component_order <- order(miu)
  biological_count <- min(
    max(1L, as.integer(biological_count)), length(component_order)
  )
  biological <- component_order[seq_len(biological_count)]
  first <- biological[1L]
  next_centre <- if (biological_count >= 2L) {
    miu[biological[2L]]
  } else {
    miu[first] + 3 * sigma[first]
  }
  # Low-coverage bins to the far left can be sequencing-error noise.  Judge
  # the amplitude on the peak top and its observed right shoulder rather than
  # letting that noise dictate the fitted height.
  local_left <- max(min(x), miu[first] - 0.35 * sigma[first])
  local_right <- min(
    max(x), 0.5 * (miu[first] + next_centre)
  )
  local_rows <- x >= local_left & x <= local_right
  if (sum(local_rows) < 4L) {
    baseline$reason <- "insufficient_first_peak_window"
    return(baseline)
  }

  count_scale <- total_samples_fit * bin_width
  component_basis <- vapply(seq_along(alpha), function(j) {
    stats::dnorm(x, miu[j], sigma[j]) * count_scale
  }, numeric(length(x)))
  if (is.null(dim(component_basis))) {
    component_basis <- matrix(component_basis, ncol = length(alpha))
  }
  prediction <- function(first_alpha) {
    alpha_value <- alpha
    alpha_value[first] <- first_alpha
    as.numeric(component_basis %*% alpha_value)
  }
  local_scale <- pmax(observed, 0.05 * observed_max, 1)
  point_weight <- 0.25 + 0.75 * sqrt(observed / observed_max)
  centre_rows <- abs(x - miu[first]) <=
    max(bin_width, 0.30 * sigma[first])
  if (!any(centre_rows)) {
    centre_rows[which.min(abs(x - miu[first]))] <- TRUE
  }
  metrics <- function(predicted) {
    global_residual <- (predicted - observed) / observed_max
    local_residual <-
      (predicted[local_rows] - observed[local_rows]) /
        local_scale[local_rows]
    centre_residual <-
      (predicted[centre_rows] - observed[centre_rows]) /
        local_scale[centre_rows]
    list(
      nrmse = sqrt(sum(point_weight * global_residual^2) /
        sum(point_weight)),
      local_mse = mean(local_residual^2),
      centre_mse = mean(centre_residual^2),
      centre_ratio = sum(predicted[centre_rows]) /
        max(sum(observed[centre_rows]), 1),
      local_over = sum(pmax(
        predicted[local_rows] - observed[local_rows], 0
      )) / max(sum(observed[local_rows]), 1),
      local_under = sum(pmax(
        observed[local_rows] - predicted[local_rows], 0
      )) / max(sum(observed[local_rows]), 1)
    )
  }
  before_prediction <- prediction(alpha[first])
  before <- metrics(before_prediction)
  objective <- function(log_scale) {
    scale <- exp(log_scale)
    candidate <- metrics(prediction(alpha[first] * scale))
    3.00 * candidate$centre_mse +
      0.50 * candidate$local_mse +
      0.10 * candidate$nrmse^2 +
      0.002 * log_scale^2
  }
  search_interval <- if (before$centre_ratio > 1.01) {
    c(log(0.05), 0)
  } else if (before$centre_ratio < 0.99) {
    c(0, log(2.50))
  } else {
    c(log(0.90), log(1.10))
  }
  solution <- tryCatch(stats::optimize(
    objective, interval = search_interval,
    tol = 1e-5
  ), error = function(error) NULL)
  if (is.null(solution) || !is.finite(solution$minimum)) {
    baseline$reason <- "optimization_failed"
    return(baseline)
  }
  scale <- exp(solution$minimum)
  baseline$triggered <- abs(log(scale)) >= log(1.01)
  if (!baseline$triggered) {
    baseline$reason <- "already_locally_balanced"
    return(baseline)
  }
  candidate_alpha <- alpha
  candidate_alpha[first] <- alpha[first] * scale
  after <- metrics(prediction(candidate_alpha[first]))
  centre_improved <-
    abs(after$centre_ratio - 1) < abs(before$centre_ratio - 1) &&
    after$centre_mse <= before$centre_mse * 0.995
  local_guard <- after$local_mse <=
    before$local_mse * 1.50 + 0.002
  global_guard <- after$nrmse <= before$nrmse * 1.10 + 0.001
  accepted <- isTRUE(centre_improved) && isTRUE(local_guard) &&
    isTRUE(global_guard)
  list(
    alpha = if (accepted) candidate_alpha else alpha,
    accepted = accepted, triggered = TRUE,
    reason = if (accepted) "accepted" else "fit_guard_rejected",
    first_component = first,
    scale = scale,
    direction = if (scale < 1) "lowered" else "raised",
    before = before, after = after
  )
}

# When two strongly overlapping early peaks miss the grey bridge continuously,
# adjust their areas and widths together in the required direction. Underfit
# may widen the overlap; overfit may lower the peak heights and/or narrow it.
# The centres stay fixed and the biological width order is preserved.
.GenomeGMM_refine_first_second_bridge <- function(
    x, y, alpha, miu, sigma, total_samples_fit,
    biological_count = length(alpha), first_sigma_upper = Inf,
    max_iterations = 220L) {
  baseline <- list(
    alpha = alpha, sigma = sigma, accepted = FALSE, triggered = FALSE,
    reason = "not_evaluated"
  )
  valid <- is.finite(x) & is.finite(y) & y >= 0
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 10L || biological_count < 2L ||
      !is.finite(total_samples_fit) || total_samples_fit <= 0 ||
      any(!is.finite(alpha)) || any(alpha <= 0) ||
      any(!is.finite(miu)) || any(!is.finite(sigma)) ||
      any(sigma <= 0)) {
    baseline$reason <- "insufficient_input"
    return(baseline)
  }
  rows <- order(x)
  x <- x[rows]
  y <- y[rows]
  raw_x <- x
  raw_y <- y
  x_steps <- diff(sort(unique(x)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  observed <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  observed[!is.finite(observed)] <- y[!is.finite(observed)]
  observed <- pmax(observed, 0)
  observed_max <- max(observed, na.rm = TRUE)
  if (!is.finite(observed_max) || observed_max <= 0) {
    baseline$reason <- "empty_observed_curve"
    return(baseline)
  }

  component_order <- order(miu)
  biological_count <- min(
    max(1L, as.integer(biological_count)), length(component_order)
  )
  biological <- component_order[seq_len(biological_count)]
  if (length(biological) < 2L) {
    baseline$reason <- "fewer_than_two_biological_peaks"
    return(baseline)
  }
  two_peak_mode <- length(biological) == 2L
  first_two <- biological[1:2]
  # A third fitted component may be the optional n+1 tail rather than a
  # biological harmonic.  It can define the right edge of the diagnostic
  # window, but it must not be required in order to repair a genuine 2x fit.
  third <- if (length(component_order) >= 3L) {
    component_order[3L]
  } else NA_integer_
  centres <- miu[first_two]
  spacing <- diff(centres)
  if (!is.finite(spacing) || spacing <= 2 * bin_width) {
    baseline$reason <- "invalid_early_peak_spacing"
    return(baseline)
  }

  # Only these bins can be changed appreciably by the first two components.
  # Restricting the numerical optimisation to them keeps large histograms fast
  # without changing the acceptance diagnostics for the affected region.
  right_anchor <- if (is.finite(third)) {
    miu[third]
  } else {
    centres[2] + spacing
  }
  left_extension <- if (isTRUE(two_peak_mode)) 0.62 else 0.20
  if (isTRUE(two_peak_mode) && length(raw_x) > 3L) {
    # Use the exact spline reference drawn in the PDF. The five-bin detector
    # remains useful for discovering peaks, but it can hide a low-coverage
    # left-flank overshoot that is visually obvious in the final figure.
    reference_rows <- if (length(raw_x) > 1L) {
      2:length(raw_x)
    } else seq_along(raw_x)
    reference_left <- max(
      5, min(raw_x[reference_rows]),
      centres[1] - left_extension * spacing
    )
    reference_right <- min(
      max(raw_x[reference_rows]), right_anchor + 0.50 * spacing
    )
    reference_grid <- seq(
      reference_left, reference_right, length.out = 420L
    )
    grey_reference <- stats::spline(
      raw_x[reference_rows], raw_y[reference_rows],
      xout = reference_grid
    )
    x <- grey_reference$x
    observed <- pmax(grey_reference$y, 0)
  } else {
    analysis_rows <- x >= centres[1] - left_extension * spacing &
      x <= right_anchor + 0.50 * spacing
    x <- x[analysis_rows]
    observed <- observed[analysis_rows]
  }
  observed_max <- max(observed, na.rm = TRUE)

  count_scale <- total_samples_fit * bin_width
  component_basis <- vapply(seq_along(alpha), function(j) {
    stats::dnorm(x, miu[j], sigma[j]) * count_scale
  }, numeric(length(x)))
  if (is.null(dim(component_basis))) {
    component_basis <- matrix(component_basis, ncol = length(alpha))
  }
  predict_curve <- function(alpha_value = alpha,
                            sigma_value = sigma,
                            miu_value = miu) {
    vapply(x, function(x_value) {
      sum(alpha_value * stats::dnorm(
        x_value, miu_value, sigma_value
      )) *
        count_scale
    }, numeric(1))
  }
  before_prediction <- as.numeric(component_basis %*% alpha)
  local_rows <- x >= centres[1] - left_extension * spacing &
    x <= right_anchor - 0.20 * spacing
  effect_rows <- x >= centres[1] - left_extension * spacing &
    x <= right_anchor + 0.50 * spacing
  local_observed_max <- max(observed[local_rows], na.rm = TRUE)
  if (!is.finite(local_observed_max) || local_observed_max <= 0) {
    baseline$reason <- "empty_local_curve"
    return(baseline)
  }
  local_scale <- pmax(observed, 0.05 * local_observed_max, 1)
  bridge_rows <- x >= centres[1] + 0.45 * spacing &
    x <= centres[2] - 0.25 * spacing
  bridge_rows <- bridge_rows & observed >= 0.05 * local_observed_max
  interval_rows <- x >= centres[1] + 0.10 * spacing &
    x <= centres[2] - 0.10 * spacing
  interval_rows <- interval_rows & observed >= 0.05 * local_observed_max
  left_flank_rows <- x >= max(5, centres[1] - left_extension * spacing) &
    x <= centres[1] - 0.05 * spacing
  left_flank_rows <- left_flank_rows &
    observed >= 0.05 * local_observed_max
  centre_rows <- lapply(first_two, function(component) {
    abs(x - miu[component]) <= max(bin_width, 0.15 * sigma[component])
  })
  if (sum(bridge_rows) < 2L || sum(interval_rows) < 3L ||
      sum(local_rows) < 6L ||
      (isTRUE(two_peak_mode) && sum(left_flank_rows) < 3L) ||
      any(vapply(centre_rows, sum, integer(1)) < 1L)) {
    baseline$reason <- "insufficient_bridge_window"
    return(baseline)
  }

  metrics <- function(prediction) {
    relative <- (prediction - observed) / local_scale
    centre_ratio <- vapply(centre_rows, function(index) {
      sum(prediction[index]) / max(sum(observed[index]), 1)
    }, numeric(1))
    bridge_ratio <- prediction[bridge_rows] /
      pmax(observed[bridge_rows], 1)
    interval_ratio <- prediction[interval_rows] /
      pmax(observed[interval_rows], 1)
    list(
      nrmse = sqrt(mean(
        ((prediction[effect_rows] - observed[effect_rows]) /
          local_scale[effect_rows])^2
      )),
      local_mse = mean(relative[local_rows]^2),
      bridge_mse = mean(relative[bridge_rows]^2),
      bridge_under = mean(pmax(1 - bridge_ratio, 0)),
      bridge_over = mean(pmax(bridge_ratio - 1, 0)),
      bridge_max_over = max(pmax(bridge_ratio - 1, 0)),
      bridge_max_under = max(pmax(1 - bridge_ratio, 0)),
      interval_mse = mean(relative[interval_rows]^2),
      interval_over = mean(pmax(interval_ratio - 1, 0)),
      interval_under = mean(pmax(1 - interval_ratio, 0)),
      interval_max_over = max(pmax(interval_ratio - 1, 0)),
      interval_max_under = max(pmax(1 - interval_ratio, 0)),
      left_flank_mse = if (any(left_flank_rows)) {
        mean(relative[left_flank_rows]^2)
      } else 0,
      left_flank_over = if (any(left_flank_rows)) {
        mean(pmax(relative[left_flank_rows], 0))
      } else 0,
      left_flank_under = if (any(left_flank_rows)) {
        mean(pmax(-relative[left_flank_rows], 0))
      } else 0,
      left_flank_max_over = if (any(left_flank_rows)) {
        max(pmax(relative[left_flank_rows], 0))
      } else 0,
      centre_ratio = centre_ratio,
      bridge_ratio = bridge_ratio,
      interval_ratio = interval_ratio
    )
  }
  before <- metrics(before_prediction)
  underfit_run <- rle(before$bridge_ratio < 0.96)
  longest_underfit <- if (any(underfit_run$values)) {
    max(underfit_run$lengths[underfit_run$values])
  } else 0L
  overfit_run <- rle(before$interval_ratio > 1.04)
  longest_overfit <- if (any(overfit_run$values)) {
    max(overfit_run$lengths[overfit_run$values])
  } else 0L
  underfit_trigger <- longest_underfit >= 2L &&
    before$bridge_under >= 0.05
  overfit_trigger <- longest_overfit >= 3L &&
    before$interval_over >= 0.03
  direction <- if (isTRUE(overfit_trigger) &&
      (!isTRUE(underfit_trigger) ||
       before$interval_over >= before$bridge_under)) {
    "overfit"
  } else if (isTRUE(underfit_trigger)) {
    "underfit"
  } else "none"
  # Trigger on consecutive disagreement, never on one noisy bin.
  baseline$triggered <- direction != "none"
  baseline$before <- before
  baseline$bridge_underfit_run <- longest_underfit
  baseline$bridge_overfit_run <- longest_overfit
  baseline$direction <- direction
  if (!isTRUE(baseline$triggered)) {
    baseline$reason <- "no_sustained_bridge_mismatch"
    return(baseline)
  }

  sigma3_cap <- if (is.finite(third)) {
    max(sigma[first_two[2]], sigma[third])
  } else {
    max(sigma[first_two[2]], 0.50 * spacing)
  }
  if (direction == "underfit") {
    sigma1_lower <- sigma[first_two[1]]
    sigma2_lower <- sigma[first_two[2]]
    sigma1_upper <- min(sigma3_cap, first_sigma_upper)
    sigma2_upper <- sigma3_cap
    if (!is.finite(sigma1_upper)) sigma1_upper <- sigma3_cap
  } else {
    width_floor <- max(bin_width, 0.35 * spacing)
    sigma1_lower <- min(
      sigma[first_two[1]], max(width_floor, sigma[first_two[1]] / 1.25)
    )
    sigma2_lower <- min(
      sigma[first_two[2]], max(width_floor, sigma[first_two[2]] / 1.25)
    )
    sigma1_upper <- if (isTRUE(two_peak_mode)) {
      min(sigma[first_two[2]], 1.05 * sigma[first_two[1]])
    } else {
      sigma[first_two[1]]
    }
    sigma2_upper <- sigma[first_two[2]]
  }

  initial_height <- alpha[first_two] /
    (sigma[first_two] * sqrt(2 * pi)) * count_scale
  decode <- function(parameters) {
    sigma_value <- sigma
    sigma_value[first_two] <- parameters[3:4]
    miu_value <- miu
    if (isTRUE(two_peak_mode) && length(parameters) >= 5L) {
      miu_value[first_two[1]] <-
        miu[first_two[1]] + spacing * parameters[5L]
    }
    peak_height <- exp(parameters[1:2])
    alpha_value <- alpha
    alpha_value[first_two] <- peak_height * sigma_value[first_two] *
      sqrt(2 * pi) / count_scale
    list(alpha = alpha_value, miu = miu_value, sigma = sigma_value)
  }
  objective <- function(parameters) {
    candidate <- decode(parameters)
    prediction <- predict_curve(
      candidate$alpha, candidate$sigma, candidate$miu
    )
    metric <- metrics(prediction)
    centre_error <- mean((metric$centre_ratio - 1)^2)
    change_penalty <- mean(
      (parameters[1:2] - log(initial_height))^2
    ) + mean(
      (log(parameters[3:4] / sigma[first_two]))^2
    )
    if (isTRUE(two_peak_mode) && length(parameters) >= 5L) {
      change_penalty <- change_penalty + 8 * parameters[5L]^2
    }
    width_order_penalty <-
      100 * pmax(parameters[3] - parameters[4], 0)^2
    target_mse <- if (direction == "overfit") {
      metric$interval_mse
    } else metric$bridge_mse
    4.0 * target_mse +
      1.0 * metric$local_mse +
      (if (isTRUE(two_peak_mode)) {
        1.5 * metric$left_flank_mse
      } else 0) +
      5.0 * centre_error +
      0.01 * change_penalty + width_order_penalty
  }
  initial <- c(
    log(initial_height), sigma[first_two],
    if (isTRUE(two_peak_mode)) 0 else NULL
  )
  lower_bound <- c(
    log(initial_height * 0.55), sigma1_lower, sigma2_lower
  )
  upper_bound <- c(
    log(initial_height * 1.45), sigma1_upper, sigma2_upper
  )
  if (isTRUE(two_peak_mode)) {
    harmonic_first <- miu[first_two[2]] / 2
    first_mean_lower <- max(
      miu[first_two[1]] - 0.05 * spacing,
      harmonic_first - 0.06 * spacing
    )
    first_mean_upper <- min(
      miu[first_two[1]] + 0.05 * spacing,
      harmonic_first + 0.06 * spacing
    )
    lower_bound <- c(
      lower_bound,
      (first_mean_lower - miu[first_two[1]]) / spacing
    )
    upper_bound <- c(
      upper_bound,
      (first_mean_upper - miu[first_two[1]]) / spacing
    )
  }
  solution <- tryCatch(stats::optim(
    initial, objective, method = "L-BFGS-B",
    lower = lower_bound,
    upper = upper_bound,
    control = list(maxit = as.integer(max_iterations), factr = 1e7)
  ), error = function(error) NULL)
  if (is.null(solution) || !is.finite(solution$value)) {
    baseline$reason <- "optimization_failed"
    return(baseline)
  }
  candidate <- decode(solution$par)
  after <- metrics(predict_curve(
    candidate$alpha, candidate$sigma, candidate$miu
  ))
  width_valid <-
    candidate$sigma[first_two[1]] <=
      candidate$sigma[first_two[2]] + 1e-8 &&
    candidate$sigma[first_two[2]] <= sigma3_cap + 1e-8
  centre_valid <- all(abs(after$centre_ratio - 1) <= 0.06)
  position_valid <- !isTRUE(two_peak_mode) ||
    abs(candidate$miu[first_two[1]] - miu[first_two[2]] / 2) <=
      0.06 * spacing + 1e-8
  bridge_improved <- if (direction == "overfit") {
    after$interval_mse <= before$interval_mse * 0.75 &&
      after$interval_over <= before$interval_over * 0.75 &&
      after$interval_max_under <= 0.06
  } else {
    after$bridge_mse <= before$bridge_mse * 0.75 &&
      after$bridge_under <= before$bridge_under * 0.75 &&
      after$bridge_max_over <= 0.06
  }
  local_fit_guard <- if (direction == "overfit") {
    # The interval-specific guards above already require a large reduction of
    # sustained overshoot. Do not reject that repair merely because the wider
    # early local window improves by slightly less than five percent.
    after$local_mse <= before$local_mse * 1.005 + 1e-4
  } else {
    after$local_mse <= before$local_mse * 0.95 + 1e-4
  }
  fit_guard <- isTRUE(local_fit_guard) &&
    after$nrmse <= before$nrmse * 1.01 + 1e-4
  left_flank_guard <- !isTRUE(two_peak_mode) || (
    after$left_flank_mse <= before$left_flank_mse * 0.98 + 1e-4 &&
    after$left_flank_max_over <= 0.05 &&
    after$left_flank_over <= 0.02
  )
  changed <- max(c(
    abs(log(candidate$sigma[first_two] / sigma[first_two])),
    abs(log(candidate$alpha[first_two] / alpha[first_two])),
    abs(candidate$miu[first_two[1]] - miu[first_two[1]]) /
      max(spacing, .Machine$double.eps)
  )) >= log(1.005)
  accepted <- isTRUE(changed) && isTRUE(width_valid) &&
    isTRUE(centre_valid) && isTRUE(position_valid) &&
    isTRUE(bridge_improved) && isTRUE(fit_guard) &&
    isTRUE(left_flank_guard)
  list(
    alpha = if (accepted) candidate$alpha else alpha,
    miu = if (accepted) candidate$miu else miu,
    sigma = if (accepted) candidate$sigma else sigma,
    accepted = accepted, triggered = TRUE,
    reason = if (accepted) "accepted" else "fit_guard_rejected",
    first_two_components = first_two,
    sigma_before = sigma[first_two],
    sigma_after = candidate$sigma[first_two],
    alpha_before = alpha[first_two],
    alpha_after = candidate$alpha[first_two],
    miu_before = miu[first_two],
    miu_after = candidate$miu[first_two],
    direction = direction,
    bridge_underfit_run = longest_underfit,
    bridge_overfit_run = longest_overfit,
    before = before, after = after
  )
}

# Once the selected peak geometry has passed its own guards, redistribute only
# the biological component masses when an interior peak collapses despite a
# smooth observed envelope.  This is deliberately a separate second stage:
# failure here keeps the already-refined centres and widths unchanged.
.GenomeGMM_refine_biological_heights <- function(
    x, y, alpha, miu, sigma, total_samples_fit,
    biological_count = length(alpha), max_iterations = 220L) {
  baseline <- list(
    alpha = alpha, accepted = FALSE, triggered = FALSE,
    reason = "not_evaluated"
  )
  valid <- is.finite(x) & is.finite(y) & y >= 0
  x <- as.numeric(x[valid])
  y <- as.numeric(y[valid])
  if (length(x) < 10L || length(alpha) < 3L ||
      !is.finite(total_samples_fit) || total_samples_fit <= 0 ||
      any(!is.finite(alpha)) || any(alpha <= 0) ||
      any(!is.finite(miu)) || any(!is.finite(sigma)) ||
      any(sigma <= 0)) {
    baseline$reason <- "insufficient_input"
    return(baseline)
  }
  rows <- order(x)
  x <- x[rows]
  y <- y[rows]
  x_steps <- diff(sort(unique(x)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  observed <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  observed[!is.finite(observed)] <- y[!is.finite(observed)]
  observed <- pmax(observed, 0)
  observed_max <- max(observed, na.rm = TRUE)
  if (!is.finite(observed_max) || observed_max <= 0) {
    baseline$reason <- "empty_observed_curve"
    return(baseline)
  }

  component_order <- order(miu)
  biological_count <- min(
    max(1L, as.integer(biological_count)), length(component_order)
  )
  biological <- component_order[seq_len(biological_count)]
  if (length(biological) < 3L) {
    baseline$reason <- "fewer_than_three_biological_peaks"
    return(baseline)
  }
  centres <- miu[biological]
  spacing <- stats::median(diff(centres))
  gap_ratio <- diff(centres) / max(spacing, 1e-12)
  if (!is.finite(spacing) || spacing <= 0 ||
      any(!is.finite(gap_ratio)) ||
      any(gap_ratio < 0.65 | gap_ratio > 1.35)) {
    baseline$reason <- "geometry_not_ready"
    return(baseline)
  }

  grey_at_centres <- stats::approx(
    x, observed, xout = centres, rule = 2
  )$y
  grey_at_centres <- pmax(grey_at_centres, 1)
  grey_curvature <- grey_at_centres[2:(biological_count - 1L)] /
    sqrt(
      grey_at_centres[1:(biological_count - 2L)] *
        grey_at_centres[3:biological_count]
    )
  # Match the local shape of the biological peak-height sequence to the local
  # shape of the grey envelope.  A five-percent band permits small numerical
  # variation, while treating both a collapsed peak and an unsupported spike
  # as invalid.  A genuine grey trough or bump moves the band with the data.
  required_curvature <- pmax(0.10, grey_curvature / 1.05)
  maximum_curvature <- pmin(10, grey_curvature * 1.05)
  height_curvature <- function(alpha_value) {
    peak_height <- alpha_value[biological] / sigma[biological]
    peak_height[2:(biological_count - 1L)] /
      sqrt(
        peak_height[1:(biological_count - 2L)] *
          peak_height[3:biological_count]
      )
  }
  baseline_curvature <- height_curvature(alpha)
  deficient <- which(
    baseline_curvature < required_curvature - 1e-8
  )
  excessive <- which(
    baseline_curvature > maximum_curvature + 1e-8
  )
  anomalous <- sort(unique(c(deficient, excessive)))
  baseline$height_curvature_before <- baseline_curvature
  baseline$required_curvature <- required_curvature
  baseline$maximum_curvature <- maximum_curvature
  if (!length(anomalous)) {
    baseline$reason <- "already_balanced"
    return(baseline)
  }
  baseline$triggered <- TRUE

  count_scale <- total_samples_fit * bin_width
  component_basis <- vapply(seq_along(alpha), function(j) {
    stats::dnorm(x, miu[j], sigma[j]) * count_scale
  }, numeric(length(x)))
  if (is.null(dim(component_basis))) {
    component_basis <- matrix(component_basis, ncol = length(alpha))
  }
  boundaries <- 0.5 * (
    centres[-biological_count] + centres[-1L]
  )
  region_id <- findInterval(x, boundaries) + 1L
  region_ids <- seq_len(biological_count)
  local_scale <- pmax(observed, 0.05 * observed_max, 1)
  point_weight <- 0.25 + 0.75 * sqrt(observed / observed_max)
  predict_curve <- function(alpha_value) {
    as.numeric(component_basis %*% alpha_value)
  }
  metrics <- function(prediction) {
    scaled_global <- (prediction - observed) / observed_max
    scaled_local <- (prediction - observed) / local_scale
    regional_mse <- vapply(region_ids, function(region) {
      region_rows <- region_id == region
      if (!any(region_rows)) return(0)
      mean(scaled_local[region_rows]^2)
    }, numeric(1))
    list(
      nrmse = sqrt(sum(point_weight * scaled_global^2) /
        sum(point_weight)),
      balanced_mse = mean(regional_mse),
      max_region_mse = max(regional_mse),
      over_area = sum(pmax(prediction - observed, 0)) /
        max(sum(observed), 1),
      under_area = sum(pmax(observed - prediction, 0)) /
        max(sum(observed), 1)
    )
  }
  before <- metrics(predict_curve(alpha))
  # Select every collapsed interior peak and its immediate neighbours from the
  # diagnostics above.  This is component-generic (not tied to peak numbers),
  # while preventing an unrelated well-fitted part of the curve from moving.
  affected_components <- sort(unique(c(
    anomalous - 1L, anomalous, anomalous + 1L,
    anomalous + 2L, anomalous + 3L
  )))
  affected_components <- affected_components[
    affected_components >= 2L &
      affected_components <= biological_count
  ]
  adjustable <- biological[affected_components]
  # Preserve the total mass within each automatically detected local region;
  # only redistribute that mass among its peaks.
  group_id <- cumsum(c(
    1L, diff(affected_components) > 1L
  ))
  affected_groups <- split(seq_along(adjustable), group_id)
  adjustable_mass <- vapply(affected_groups, function(indices) {
    sum(alpha[adjustable[indices]])
  }, numeric(1))
  alpha_reference <- alpha[adjustable]

  decode <- function(parameters) {
    adjustable_alpha <- exp(parameters)
    for (group_index in seq_along(affected_groups)) {
      indices <- affected_groups[[group_index]]
      adjustable_alpha[indices] <- adjustable_alpha[indices] *
        adjustable_mass[group_index] /
        sum(adjustable_alpha[indices])
    }
    alpha_value <- alpha
    alpha_value[adjustable] <- adjustable_alpha
    alpha_value
  }
  objective <- function(parameters) {
    alpha_value <- decode(parameters)
    metric <- metrics(predict_curve(alpha_value))
    curvature <- height_curvature(alpha_value)
    shortfall <- pmax(
      log(required_curvature) - log(pmax(curvature, 1e-12)), 0
    )
    excess <- pmax(
      log(pmax(curvature, 1e-12)) - log(maximum_curvature), 0
    )
    alpha_change <- mean(
      (log(alpha_value[adjustable] / alpha_reference))^2
    )
    0.65 * metric$nrmse^2 +
      0.35 * metric$balanced_mse +
      1.00 * mean(shortfall^2 + excess^2) +
      0.002 * alpha_change
  }
  initial <- log(alpha_reference)
  span <- rep(0.90, length(adjustable))
  anomalous_components <- anomalous + 1L
  anomalous_parameters <- match(
    biological[anomalous_components], adjustable, nomatch = 0L
  )
  anomalous_parameters <- anomalous_parameters[
    anomalous_parameters > 0L
  ]
  required_adjustment <- pmax(
    log(required_curvature[anomalous] /
      pmax(baseline_curvature[anomalous], 1e-12)),
    log(pmax(baseline_curvature[anomalous], 1e-12) /
      maximum_curvature[anomalous])
  )
  span[anomalous_parameters] <- pmin(
    log(100), pmax(1.40, 1.50 * required_adjustment + 0.50)
  )
  solution <- tryCatch(stats::optim(
    initial, objective, method = "L-BFGS-B",
    lower = initial - span, upper = initial + span,
    control = list(maxit = as.integer(max_iterations), factr = 1e8)
  ), error = function(error) NULL)
  if (is.null(solution) || !is.finite(solution$value)) {
    baseline$reason <- "optimization_failed"
    return(baseline)
  }
  candidate_alpha <- decode(solution$par)
  candidate_curvature <- height_curvature(candidate_alpha)
  # If the numerical optimiser stops just short of the biological boundary,
  # project the same automatically selected local regions onto the admissible
  # height-curvature set.  Peak positions, widths and regional mass stay fixed.
  project_height_curvature <- function(alpha_value) {
    for (iteration in seq_len(250L)) {
      curvature <- height_curvature(alpha_value)
      low <- which(curvature < required_curvature - 1e-6)
      high <- which(curvature > maximum_curvature + 1e-6)
      bad <- sort(unique(c(low, high)))
      if (!length(bad)) break
      peak_height <- alpha_value[biological] / sigma[biological]
      changed_projection <- FALSE
      for (curvature_index in bad) {
        component_position <- curvature_index + 1L
        component <- biological[component_position]
        parameter_position <- match(component, adjustable, nomatch = 0L)
        if (parameter_position == 0L) next
        neighbour_height <- sqrt(
          peak_height[component_position - 1L] *
            peak_height[component_position + 1L]
        )
        target_height <- if (curvature_index %in% low) {
          1.002 * required_curvature[curvature_index] *
            neighbour_height
        } else {
          0.998 * maximum_curvature[curvature_index] *
            neighbour_height
        }
        needs_change <- is.finite(target_height) &&
          ((curvature_index %in% low &&
              target_height > peak_height[component_position]) ||
             (curvature_index %in% high &&
              target_height < peak_height[component_position]))
        if (needs_change) {
          alpha_value[component] <- target_height * sigma[component]
          peak_height[component_position] <- target_height
          changed_projection <- TRUE
        }
      }
      for (group_index in seq_along(affected_groups)) {
        parameter_indices <- affected_groups[[group_index]]
        components <- adjustable[parameter_indices]
        alpha_value[components] <- alpha_value[components] *
          adjustable_mass[group_index] /
          sum(alpha_value[components])
      }
      if (!changed_projection) break
    }
    alpha_value
  }
  candidate_outside <- function(curvature) {
    any(curvature < 0.98 * required_curvature - 1e-8) ||
      any(curvature > 1.02 * maximum_curvature + 1e-8)
  }
  curvature_violation <- function(curvature) {
    lower_violation <- pmax(
      log(required_curvature) - log(pmax(curvature, 1e-12)), 0
    )
    upper_violation <- pmax(
      log(pmax(curvature, 1e-12)) - log(maximum_curvature), 0
    )
    sum(lower_violation^2 + upper_violation^2)
  }
  if (candidate_outside(candidate_curvature)) {
    projected_alpha <- project_height_curvature(candidate_alpha)
    projected_curvature <- height_curvature(projected_alpha)
    if (curvature_violation(projected_curvature) <
        curvature_violation(candidate_curvature)) {
      candidate_alpha <- projected_alpha
      candidate_curvature <- projected_curvature
    }
  }
  # Extremely collapsed solutions can propagate a new shortfall just outside
  # the initially selected local region.  In that rare case, repeat the same
  # projection over all biological peaks except the fixed first peak.  This is
  # a fallback redistribution only; total biological mass remains unchanged.
  if (candidate_outside(candidate_curvature)) {
    global_adjustable <- biological[-1L]
    global_mass <- sum(alpha[global_adjustable])
    global_alpha <- candidate_alpha
    for (iteration in seq_len(400L)) {
      curvature <- height_curvature(global_alpha)
      low <- which(curvature < required_curvature - 1e-6)
      high <- which(curvature > maximum_curvature + 1e-6)
      bad <- sort(unique(c(low, high)))
      if (!length(bad)) break
      peak_height <- global_alpha[biological] / sigma[biological]
      for (curvature_index in bad) {
        component_position <- curvature_index + 1L
        component <- biological[component_position]
        neighbour_height <- sqrt(
          peak_height[component_position - 1L] *
            peak_height[component_position + 1L]
        )
        target_height <- if (curvature_index %in% low) {
          1.002 * required_curvature[curvature_index] *
            neighbour_height
        } else {
          0.998 * maximum_curvature[curvature_index] *
            neighbour_height
        }
        needs_change <- is.finite(target_height) &&
          ((curvature_index %in% low &&
              target_height > peak_height[component_position]) ||
             (curvature_index %in% high &&
              target_height < peak_height[component_position]))
        if (needs_change) {
          global_alpha[component] <- target_height * sigma[component]
          peak_height[component_position] <- target_height
        }
      }
      global_alpha[global_adjustable] <-
        global_alpha[global_adjustable] * global_mass /
          sum(global_alpha[global_adjustable])
    }
    global_curvature <- height_curvature(global_alpha)
    if (curvature_violation(global_curvature) <
        curvature_violation(candidate_curvature)) {
      candidate_alpha <- global_alpha
      candidate_curvature <- global_curvature
    }
  }
  after <- metrics(predict_curve(candidate_alpha))
  height_valid <-
    all(candidate_curvature >= 0.98 * required_curvature - 1e-8) &&
    all(candidate_curvature <= 1.02 * maximum_curvature + 1e-8)
  fit_guard <- after$nrmse <= before$nrmse * 1.015 + 5e-4 &&
    after$balanced_mse <= before$balanced_mse * 1.04 + 2e-4 &&
    after$max_region_mse <= before$max_region_mse * 1.06 + 2e-4 &&
    after$over_area <= before$over_area + 0.004 &&
    after$under_area <= before$under_area + 0.004
  baseline_violation <- curvature_violation(baseline_curvature)
  candidate_violation <- curvature_violation(candidate_curvature)
  improved_support <- candidate_violation < baseline_violation - 1e-6
  # A collapsed biological peak is never retained merely because the
  # unconstrained total curve has a lower numerical loss.  The optimiser has
  # already selected the least-loss admissible redistribution; fit_guard is
  # retained as a diagnostic rather than a reason to restore the bad peak.
  accepted <- isTRUE(height_valid) && isTRUE(improved_support)
  list(
    alpha = if (accepted) candidate_alpha else alpha,
    accepted = accepted, triggered = TRUE,
    reason = if (accepted && isTRUE(fit_guard)) {
      "accepted"
    } else if (accepted) {
      "accepted_biological_height_repair"
    } else {
      "height_projection_failed"
    },
    before = before, after = after,
    candidate_alpha = candidate_alpha,
    fit_guard_passed = isTRUE(fit_guard),
    deficient_components = deficient + 1L,
    excessive_components = excessive + 1L,
    height_curvature_before = baseline_curvature,
    height_curvature_after = candidate_curvature,
    required_curvature = required_curvature,
    maximum_curvature = maximum_curvature
  )
}

# Post-fit calibration is allowed to refine an already constrained peak, but it
# must not undo the biological peak-position structure established by EM.  In
# particular, the calibrated centre may move only a small fraction of one
# harmonic spacing, may not move farther away from its integer-multiple target,
# and may not collapse onto a neighbouring component.
.GenomeGMM_check_main_peak_calibration <- function(
    miu, sigma, main_peak_component, calibrated_mean, calibrated_sd,
    first_peak_x = NULL, max_sigma_ratio = 1.35,
    movement_fraction = 0.10, minimum_gap_fraction = 0.70,
    maximum_gap_fraction = 1.30) {
  component <- as.integer(main_peak_component)[1]
  valid_component <- is.finite(component) && component >= 1L &&
    component <= length(miu) && component <= length(sigma)
  if (!valid_component || !is.finite(calibrated_mean) ||
      !is.finite(calibrated_sd) || calibrated_sd <= 0) {
    return(list(accepted = FALSE, reason = "non-finite calibrated parameter"))
  }

  spacing <- suppressWarnings(as.numeric(first_peak_x)[1])
  if (!is.finite(spacing) || spacing <= 0) {
    ordered_centres <- sort(miu[is.finite(miu)])
    positive_gaps <- diff(ordered_centres)
    positive_gaps <- positive_gaps[is.finite(positive_gaps) & positive_gaps > 0]
    spacing <- if (length(positive_gaps)) stats::median(positive_gaps) else NA_real_
  }
  if (!is.finite(spacing) || spacing <= 0) {
    return(list(accepted = FALSE, reason = "harmonic spacing unavailable"))
  }

  original_mean <- miu[component]
  original_sd <- sigma[component]
  harmonic_target <- component * spacing
  harmonic_limit <- 0.20 * spacing
  minimum_gap <- minimum_gap_fraction * spacing
  maximum_gap <- maximum_gap_fraction * spacing

  # Intersect every permissible interval, then project the data-driven proposal
  # onto that interval.  Thus calibration is retained whenever a biologically
  # valid small adjustment exists; an extreme proposal is clipped, not allowed
  # to destroy peak order or harmonic spacing.
  lower <- max(
    original_mean - movement_fraction * spacing,
    harmonic_target - harmonic_limit
  )
  upper <- min(
    original_mean + movement_fraction * spacing,
    harmonic_target + harmonic_limit
  )
  if (component > 1L && is.finite(miu[component - 1L])) {
    lower <- max(lower, miu[component - 1L] + minimum_gap)
    # The previous component is a biological harmonic.  Protect both sides of
    # its expected one-spacing relationship with the calibrated main peak.
    upper <- min(upper, miu[component - 1L] + maximum_gap)
  }
  if (component < length(miu) && is.finite(miu[component + 1L])) {
    # A following component can be the optional n+1 tail component.  It only
    # needs a non-overlap guard and is not forced to remain within 1.30 d.
    upper <- min(upper, miu[component + 1L] - minimum_gap)
  }
  if (!is.finite(lower) || !is.finite(upper) || lower > upper + 1e-8) {
    return(list(
      accepted = FALSE,
      reason = "no feasible position satisfying harmonic and neighbour guards",
      spacing = spacing,
      lower = lower,
      upper = upper
    ))
  }
  projected_mean <- min(max(calibrated_mean, lower), upper)

  width_lower <- original_sd / max_sigma_ratio
  width_upper <- original_sd * max_sigma_ratio
  projected_sd <- min(max(calibrated_sd, width_lower), width_upper)
  projected <- abs(projected_mean - calibrated_mean) > 1e-8 ||
    abs(projected_sd - calibrated_sd) > 1e-8

  list(
    accepted = TRUE,
    reason = if (projected) "projected_to_biological_range" else "accepted",
    spacing = spacing,
    calibrated_mean = calibrated_mean,
    projected_mean = projected_mean,
    calibrated_sd = calibrated_sd,
    projected_sd = projected_sd,
    movement = abs(projected_mean - original_mean),
    movement_limit = movement_fraction * spacing,
    previous_gap = if (component > 1L) {
      projected_mean - miu[component - 1L]
    } else NA_real_,
    following_gap = if (component < length(miu)) {
      miu[component + 1L] - projected_mean
    } else NA_real_,
    minimum_gap = minimum_gap,
    maximum_gap = maximum_gap,
    harmonic_error = abs(projected_mean - harmonic_target),
    harmonic_limit = harmonic_limit,
    width_ratio = projected_sd / original_sd,
    lower = lower,
    upper = upper,
    projected = projected
  )
}

# The preserved known-n fitting engine is embedded in a private environment.
# Keeping it here makes this file self-contained while preventing its helper
# functions from colliding with the automatic ploidy-inference implementation.
.GenomeGMM_embedded_legacy_environment <- local({
  legacy_environment <- new.env(parent = globalenv())
  evalq({
    # Keep every requested biological component numerically present in automatic
    # fits.  This is intentionally a very small floor: it prevents a harmonic from
    # collapsing to effectively zero (as in virginalis) without forcing weak peaks
    # to carry the same mass as clear peaks.  When the original peak-count rule has
    # actually instantiated n+1, keep that fitted tail component numerically
    # present as well; otherwise the PDF can advertise n+1 while drawing a zero
    # component.
    stabilize_biological_alpha <- function(alpha, required_components,
                                           minimum_share = 0.005,
                                           extra_minimum_share = 0.0025) {
      alpha <- as.numeric(alpha)
      total <- sum(alpha)
      required_components <- min(as.integer(required_components), length(alpha))
      if (!is.finite(total) || total <= 0 || required_components < 1L) return(alpha)
      proportions <- pmax(alpha / total, 0)
      proportions[seq_len(required_components)] <- pmax(
        proportions[seq_len(required_components)], minimum_share
      )
      if (length(proportions) > required_components) {
        extra_indices <- seq.int(required_components + 1L,
                                 length(proportions))
        proportions[extra_indices] <- pmax(
          proportions[extra_indices], extra_minimum_share
        )
      }
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
                          .preloaded_data = NULL, .analysis_max_x = NULL,
                          .first_peak_override = NULL,
                          .start_index_override = NULL,
                          .prefix_completion_slot = NULL,
                          .prefix_completion_weight = 0.35,
                          .visible_peak_override = NULL,
                          .prefix_discovery_only = FALSE, ...) {
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
      full_data <- full_data[
        is.finite(as.numeric(full_data[, 1])) &
          is.finite(as.numeric(full_data[, 2])) &
          as.numeric(full_data[, 1]) > 0,
        , drop = FALSE
      ]
      observed_full_data <- full_data

      # For a one-visible-peak candidate, first infer mu1 = q / m.  Reconstruct
      # only the missing/contaminated LEFT half of that first peak by mirroring
      # its observed right side.  The reconstructed bins guide EM at low weight;
      # the untouched histogram is retained for diagnostics and plotting.
      prefix_completion_active <- !is.null(.prefix_completion_slot) &&
        length(.prefix_completion_slot) == 1L &&
        is.finite(.prefix_completion_slot) &&
        .prefix_completion_slot >= 2L &&
        !is.null(.visible_peak_override) &&
        length(.visible_peak_override) == 1L &&
        is.finite(.visible_peak_override) && .visible_peak_override > 0
      prefix_model_start_index <- NA_integer_
      prefix_noise_start_index <- NA_integer_
      prefix_first_peak_x <- NA_real_
      prefix_sigma_seed <- NA_real_
      if (isTRUE(prefix_completion_active)) {
        prefix_slot <- as.integer(.prefix_completion_slot)
        visible_peak <- as.numeric(.visible_peak_override)
        prefix_first_peak_x <- visible_peak / prefix_slot
        x_steps <- diff(sort(unique(as.numeric(full_data[, 1]))))
        x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
        prefix_bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
        # With only one explicit terminal peak, its observed width includes
        # unresolved neighbours.  Use a conservative fraction of the inferred
        # harmonic spacing for the first peak, consistent with the existing
        # no-abrupt-width-change prior.
        prefix_sigma_seed <- max(prefix_bin_width,
          0.35 * prefix_first_peak_x)
        model_start_x <- max(min(full_data[, 1]),
          prefix_first_peak_x - 2 * prefix_sigma_seed)
        prefix_model_start_index <- which(full_data[, 1] >= model_start_x)[1]
        prefix_noise_start_index <- if (!is.null(.start_index_override) &&
            length(.start_index_override) == 1L &&
            is.finite(.start_index_override)) {
          max(1L, min(nrow(full_data), as.integer(.start_index_override)))
        } else prefix_model_start_index
        smooth_observed <- as.numeric(stats::filter(
          observed_full_data[, 2], rep(1 / 5, 5), sides = 2
        ))
        smooth_observed[!is.finite(smooth_observed)] <-
          observed_full_data[!is.finite(smooth_observed), 2]
        left_indices <- which(
          full_data[, 1] >= model_start_x &
            full_data[, 1] < prefix_first_peak_x
        )
        if (length(left_indices)) {
          mirror_x <- 2 * prefix_first_peak_x - full_data[left_indices, 1]
          mirror_y <- stats::approx(
            observed_full_data[, 1], smooth_observed,
            xout = mirror_x, rule = 2
          )$y
          usable <- is.finite(mirror_y) & mirror_y >= 0
          full_data[left_indices[usable], 2] <- mirror_y[usable]
        }
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
      if (isTRUE(prefix_completion_active) &&
          is.finite(prefix_model_start_index)) {
        start_index <- prefix_model_start_index
      } else if (!is.null(.start_index_override) &&
          length(.start_index_override) == 1L &&
          is.finite(.start_index_override)) {
        start_index <- max(1L, min(nrow(full_data),
          as.integer(.start_index_override)))
      }
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
      if (!is.null(.first_peak_override) &&
          length(.first_peak_override) == 1L &&
          is.finite(.first_peak_override) && .first_peak_override > 0) {
        first_peak_x <- as.numeric(.first_peak_override)
      }
      if (isTRUE(prefix_completion_active)) {
        first_peak_x <- prefix_first_peak_x
        # Prefix recovery locates the first harmonic, but the original visible
        # peak remains the stable reference for peak-window weighting and
        # symmetry.  Treating the reconstructed first peak as the weighting
        # anchor can make sigma approach the harmonic spacing and collapse the
        # second component into the first.
        highest_peak_x <- as.numeric(.visible_peak_override)
      }
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
          right_value <- first_peak_x * peak_multiplier  
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
      if (isTRUE(prefix_completion_active)) {
        synthetic_left <- x_fit_processed < prefix_first_peak_x &
          x_fit_processed >= full_data[prefix_model_start_index, 1]
        weight_fit_processed[synthetic_left] <-
          weight_fit_processed[synthetic_left] *
          max(0.05, min(0.75, as.numeric(.prefix_completion_weight)))
      }
      
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
        sigma_mode = sigma_mode,
        # The reconstructed first harmonic is estimated from its right flank.
        # A moderate 0.35--0.40 spacing range avoids both failure modes seen in
        # hidden-prefix data: a broad first Gaussian swallowing peak 2, and an
        # over-narrow first Gaussian producing an artificial high spike/valley.
        # This prefix-specific refinement constrains component 1 only. Later
        # widths remain those selected by the ordinary width model.
        maximum_sigma_fraction = 0.80,
        first_sigma_fraction_range = if (isTRUE(prefix_completion_active) &&
            !isTRUE(.prefix_discovery_only)) {
          c(0.35, 0.40)
        } else NULL,
        max_adjacent_sigma_ratio = Inf
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
        if (isTRUE(prefix_completion_active)) {
          # Prefix recovery adds one special identifiability condition only:
          # the reconstructed first Gaussian must have a moderate width. Do
          # not add separate hard width or mass rules for components 2..K;
          # ordinary candidate scoring and component-health checks already
          # decide those parameters for every data type.
          for (profile in names(candidate_values)) {
            value <- candidate_values[[profile]]
            biological <- order(value$miu)[
              seq_len(min(n, length(value$miu)))
            ]
            biological_sigma <- value$sigma[biological]
            identifiable_width <-
              biological_sigma[1] >= 0.35 * prefix_first_peak_x - 1e-8 &&
              biological_sigma[1] <= 0.40 * prefix_first_peak_x + 1e-8
            value$eligible <- isTRUE(value$eligible) &&
              isTRUE(identifiable_width)
            candidate_values[[profile]] <- value
          }
        }
        eligible_profiles <- names(candidate_values)[vapply(
          candidate_values, function(value) isTRUE(value$eligible), logical(1)
        )]
        if (!length(eligible_profiles)) {
          stop(
            "Hidden-prefix fit rejected: no identifiable joint-width candidate.",
            call. = FALSE
          )
        }
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
      # Prefix completion can occasionally recover the missing k-mer mass by
      # leaving the first Gaussian slightly too broad on its observed right
      # flank. Refine only the first width and its weight after the ordinary
      # joint optimiser. Component 2 and every later component must remain as
      # selected by the ordinary common/trend/adaptive fitter; otherwise a
      # first-flank correction can incorrectly suppress a genuine second peak.
      prefix_flank_refined <- FALSE
      prefix_flank_scale <- 1
      prefix_flank_before <- NA_real_
      prefix_flank_after <- NA_real_
      prefix_interval_overshoot <- NA_real_
      prefix_interval_overshoot_run <- NA_integer_
      prefix_interval_height_scale <- 1
      prefix_peak_height_ratio_before <- NA_real_
      prefix_peak_height_ratio_after <- NA_real_
      prefix_first_height_refinement <- NULL
      prefix_first_height_refinement_after_shape <- NULL
      prefix_first_shape_refinement <- NULL
      prefix_tail_smoothing_accepted <- FALSE
      prefix_tail_sigma_before <- NULL
      prefix_tail_sigma_after <- NULL
      prefix_tail_nrmse_before <- NA_real_
      prefix_tail_nrmse_after <- NA_real_
      prefix_tail_ripple_before <- NA_real_
      prefix_tail_ripple_after <- NA_real_
      # The recovered first component has its own width/flank protection.
      # Do not subsequently rewrite components 2..K here: their widths and
      # weights must remain those selected by the ordinary common/trend/
      # adaptive fitter.  The experimental residual-tail optimiser is kept
      # below only as an internal rollback record and is deliberately inactive.
      prefix_tail_reoptimization_enabled <- FALSE
      if (isTRUE(prefix_completion_active) && n >= 2L) {
        biological_order <- order(miu)[seq_len(min(n, length(miu)))]
        first_component <- biological_order[1]
        second_component <- biological_order[2]
        scoring_rows <- which(
          observed_full_data[, 1] >=
            observed_full_data[prefix_noise_start_index, 1] &
            observed_full_data[, 1] <= right_value
        )
        original_x <- observed_full_data[scoring_rows, 1]
        original_y <- observed_full_data[scoring_rows, 2]
        evaluate_prefix_flank <- function(alpha_value, sigma_value) {
          component_curve <- vapply(seq_along(alpha_value), function(j) {
            alpha_value[j] * stats::dnorm(
              original_x, miu[j], sigma_value[j]
            ) * sum(y_fit_processed) * prefix_bin_width
          }, numeric(length(original_x)))
          prediction <- if (is.null(dim(component_curve))) {
            component_curve
          } else rowSums(component_curve)
          prediction[!is.finite(prediction)] <- 0
          peak_scale <- max(original_y, na.rm = TRUE)
          point_weights <- 0.25 + 0.75 * sqrt(
            pmax(original_y, 0) / max(peak_scale, 1)
          )
          nrmse_value <- sqrt(sum(point_weights *
            ((prediction - original_y) / max(peak_scale, 1))^2) /
            sum(point_weights))
          grey_mass <- max(sum(original_x * original_y), 1)
          under_value <- sum(original_x *
            pmax(original_y - prediction, 0)) / grey_mass
          over_value <- sum(original_x *
            pmax(prediction - original_y, 0)) / grey_mass
          right_rows <- original_x >= miu[first_component] &
            original_x <= miu[second_component]
          # Use the same unweighted count-area definition as the final guard.
          # A different x-weighted definition here previously allowed a flank
          # candidate to pass and then fail immediately in the final check.
          right_mass <- max(sum(original_y[right_rows]), 1)
          right_over <- sum(pmax(
            prediction[right_rows] - original_y[right_rows], 0
          )) / right_mass
          right_relative_over <- (prediction[right_rows] -
            original_y[right_rows]) / pmax(
              original_y[right_rows],
              0.05 * max(original_y[right_rows]), 1
            )
          right_runs <- rle(right_relative_over > 0.08)
          right_over_run <- if (any(right_runs$values)) {
            max(right_runs$lengths[right_runs$values])
          } else 0L
          score_value <- 0.38 * nrmse_value + 0.20 * under_value +
            0.12 * over_value + 0.30 * right_over
          list(score = score_value, nrmse = nrmse_value,
               under = under_value, over = over_value,
               right_over = right_over,
               right_over_run = as.integer(right_over_run))
        }
        baseline_flank <- evaluate_prefix_flank(alpha, sigma)
        prefix_flank_before <- baseline_flank$right_over
        # Search the admissible first-width interval directly. Do not derive
        # its lower bound from sigma2: that would turn a constraint intended
        # for the recovered first peak into another hard restriction on later
        # widths and can leave no usable first-flank candidate.
        minimum_scale <- max(
          0.75,
          0.35 * prefix_first_peak_x /
            max(sigma[first_component], .Machine$double.eps)
        )
        scale_grid <- unique(c(
          seq(min(1, minimum_scale), 1, length.out = 11L), 1
        ))
        candidate_grid <- data.frame(
          first_scale = scale_grid[scale_grid < 0.999]
        )
        flank_candidates <- lapply(seq_len(nrow(candidate_grid)),
          function(candidate_index) {
            scale_value <- candidate_grid$first_scale[candidate_index]
            candidate_alpha <- alpha
            candidate_sigma <- sigma
            candidate_sigma[first_component] <-
              sigma[first_component] * scale_value
            # Keeping alpha/sigma fixed preserves the Gaussian height at mu1.
            candidate_alpha[first_component] <-
              alpha[first_component] * scale_value
            metric <- evaluate_prefix_flank(
              candidate_alpha, candidate_sigma
            )
            eligible <- metric$nrmse <= baseline_flank$nrmse * 1.03 + 1e-4 &&
              metric$under <= baseline_flank$under * 1.18 + 0.003 &&
              metric$over <= baseline_flank$over + 0.01 &&
              metric$right_over <= baseline_flank$right_over * 0.85 &&
              metric$right_over <= 0.05 &&
              metric$right_over_run < 3L &&
              metric$score <= baseline_flank$score * 0.995
            list(scale = scale_value, alpha = candidate_alpha,
                 sigma = candidate_sigma, metric = metric,
                 eligible = eligible)
          })
        eligible_flanks <- flank_candidates[vapply(
          flank_candidates, function(value) isTRUE(value$eligible), logical(1)
        )]
        if (length(eligible_flanks)) {
          flank_scores <- vapply(
            eligible_flanks, function(value) value$metric$score, numeric(1)
          )
          selected_flank <- eligible_flanks[[which.min(flank_scores)]]
          alpha <- selected_flank$alpha
          sigma <- selected_flank$sigma
          EM_result$alpha <- alpha
          EM_result$sigma <- sigma
          prefix_flank_refined <- TRUE
          prefix_flank_scale <- selected_flank$scale
          prefix_flank_after <- selected_flank$metric$right_over
        } else {
          prefix_flank_after <- prefix_flank_before
        }
        cat("Recovered first-peak right-flank refinement:",
            if (prefix_flank_refined) "accepted" else "not needed/rejected",
            "\n")
        cat("First-peak right-flank overshoot before/after:",
            round(prefix_flank_before, 5), "/",
            round(prefix_flank_after, 5), "\n")
        if (prefix_flank_refined) {
          cat("First-width refinement scale:",
              round(prefix_flank_scale, 4), "\n")
        }
      }
      if (isTRUE(prefix_completion_active) &&
          !isTRUE(.prefix_discovery_only) && n >= 2L) {
        # The recovered first peak is an estimate, so its weight is not tied to
        # a fixed target height.  Fit alpha1 from the whole observed right half
        # of that peak; the centre and selected width remain unchanged.
        biological <- order(miu)[seq_len(min(n, length(miu)))]
        first_component <- biological[1]
        first_centre <- miu[first_component]
        observed_first_height <- stats::approx(
          observed_full_data[, 1], observed_full_data[, 2],
          xout = first_centre, rule = 2
        )$y
        prefix_first_height_refinement <-
          .GenomeGMM_refine_estimated_first_height(
            x = observed_full_data[, 1],
            y = observed_full_data[, 2],
            alpha = alpha, miu = miu, sigma = sigma,
            count_scale = sum(y_fit_processed) * prefix_bin_width,
            spacing = prefix_first_peak_x,
            biological_count = n
          )
        initial_peak_ratio <-
          prefix_first_height_refinement$before$centre_ratio
        if (!length(initial_peak_ratio) || !is.finite(initial_peak_ratio)) {
          initial_peak_ratio <- sum(alpha * stats::dnorm(
            first_centre, miu, sigma
          )) * sum(y_fit_processed) * prefix_bin_width /
            max(observed_first_height, 1)
        }
        prefix_peak_height_ratio_before <- initial_peak_ratio
        if (isTRUE(prefix_first_height_refinement$accepted)) {
          alpha <- prefix_first_height_refinement$alpha
          EM_result$alpha <- alpha
        }
        refined_peak_ratio <- if (
          isTRUE(prefix_first_height_refinement$accepted)
        ) prefix_first_height_refinement$after$centre_ratio else
          initial_peak_ratio
        prefix_peak_height_ratio_after <- refined_peak_ratio
      }
      prefix_residual_refinement <- NULL
      prefix_weight_refinement <- NULL
      interpeak_shape_refinement <- NULL
      if (isTRUE(prefix_completion_active) &&
          !isTRUE(.prefix_discovery_only) && n >= 2L) {
        residual_rows <- observed_full_data[, 1] <= right_value
        prefix_residual_refinement <- .GenomeGMM_refine_prefix_residual(
          x = observed_full_data[residual_rows, 1],
          y = observed_full_data[residual_rows, 2],
          alpha = alpha, miu = miu, sigma = sigma,
          count_scale = sum(y_fit_processed) * prefix_bin_width,
          spacing = prefix_first_peak_x
        )
        if (isTRUE(prefix_residual_refinement$accepted)) {
          alpha <- prefix_residual_refinement$alpha
          miu <- prefix_residual_refinement$miu
          sigma <- prefix_residual_refinement$sigma
          EM_result$alpha <- alpha
          EM_result$miu <- miu
          EM_result$sigma <- sigma
          # The component height is preserved under translation. Compare the
          # total peak height to the same original reference height used above;
          # the contaminated grey prefix at the new centre is not a new target.
          first_centre <- miu[first_component]
          prefix_peak_height_ratio_after <- sum(alpha * stats::dnorm(
            first_centre, miu, sigma)) * sum(y_fit_processed) *
            prefix_bin_width / max(observed_first_height, 1)
        }
        prefix_weight_refinement <- .GenomeGMM_refine_prefix_weights(
          x = observed_full_data[residual_rows, 1],
          y = observed_full_data[residual_rows, 2],
          alpha = alpha, miu = miu, sigma = sigma,
          count_scale = sum(y_fit_processed) * prefix_bin_width,
          spacing = prefix_first_peak_x,
          biological_count = n
        )
        if (isTRUE(prefix_weight_refinement$accepted)) {
          alpha <- prefix_weight_refinement$alpha
          EM_result$alpha <- alpha
          first_centre <- miu[first_component]
          prefix_peak_height_ratio_after <- sum(alpha * stats::dnorm(
            first_centre, miu, sigma)) * sum(y_fit_processed) *
            prefix_bin_width / max(observed_first_height, 1)
        }
        interpeak_shape_refinement <- .GenomeGMM_refine_interpeak_shape(
          x = observed_full_data[residual_rows, 1],
          y = observed_full_data[residual_rows, 2],
          alpha = alpha, miu = miu, sigma = sigma,
          count_scale = sum(y_fit_processed) * prefix_bin_width,
          spacing = prefix_first_peak_x,
          biological_count = n
        )
        if (isTRUE(interpeak_shape_refinement$accepted)) {
          alpha <- interpeak_shape_refinement$alpha
          sigma <- interpeak_shape_refinement$sigma
          EM_result$alpha <- alpha
          EM_result$sigma <- sigma
          first_centre <- miu[first_component]
          prefix_peak_height_ratio_after <- sum(alpha * stats::dnorm(
            first_centre, miu, sigma)) * sum(y_fit_processed) *
            prefix_bin_width / max(observed_first_height, 1)
        }
        # Final coordinate pass: after neighbouring components have settled,
        # jointly optimise the recovered first peak's height, centre and width
        # against the total red/grey curve rather than preserving any seed.
        # Candidate ranking must compare the ordinary fitted models on equal
        # terms. Apply this final visual/data refinement only to the selected
        # fit (the call that actually writes the PDF), not during screening.
        if (isTRUE(save_plot)) {
          prefix_first_height_refinement_after_shape <-
            .GenomeGMM_refine_estimated_first_height(
              x = observed_full_data[residual_rows, 1],
              y = observed_full_data[residual_rows, 2],
              alpha = alpha, miu = miu, sigma = sigma,
              count_scale = sum(y_fit_processed) * prefix_bin_width,
              spacing = prefix_first_peak_x,
              biological_count = n
            )
          if (isTRUE(prefix_first_height_refinement_after_shape$accepted)) {
            alpha <- prefix_first_height_refinement_after_shape$alpha
            EM_result$alpha <- alpha
          }
          prefix_first_shape_refinement <-
            .GenomeGMM_refine_estimated_first_shape(
              x = observed_full_data[residual_rows, 1],
              y = observed_full_data[residual_rows, 2],
              alpha = alpha, miu = miu, sigma = sigma,
              count_scale = sum(y_fit_processed) * prefix_bin_width,
              spacing = prefix_first_peak_x,
              biological_count = n
            )
          if (isTRUE(prefix_first_shape_refinement$accepted)) {
            alpha <- prefix_first_shape_refinement$alpha
            miu <- prefix_first_shape_refinement$miu
            sigma <- prefix_first_shape_refinement$sigma
            EM_result$alpha <- alpha
            EM_result$miu <- miu
            EM_result$sigma <- sigma
            first_centre <- miu[first_component]
            current_observed_first_height <- stats::approx(
              observed_full_data[, 1], observed_full_data[, 2],
              xout = first_centre, rule = 2
            )$y
            prefix_peak_height_ratio_after <-
              sum(alpha * stats::dnorm(first_centre, miu, sigma)) *
              sum(y_fit_processed) * prefix_bin_width /
              max(current_observed_first_height, 1)
          }
        }
      }
      if (isTRUE(prefix_tail_reoptimization_enabled) &&
          isTRUE(prefix_completion_active) &&
          !isTRUE(.prefix_discovery_only) && n >= 3L) {
        # Keep the recovered first Gaussian fixed after its width has passed
        # the biological and flank-fit guards. Optimise every later width and
        # weight independently against the observed residuals; later peaks
        # may widen or narrow, while adjacent widths remain within the general
        # max_sigma_ratio constraint. This is data-driven, not a species rule.
        biological <- order(miu)[seq_len(min(n, length(miu)))]
        first_component <- biological[1]
        second_component <- biological[2]
        all_order <- order(miu)
        adjustable <- setdiff(seq_along(alpha), first_component)
        tail_start <- 0.5 * (
          miu[first_component] + miu[second_component]
        )
        tail_rows <- which(
          observed_full_data[, 1] >= tail_start &
            observed_full_data[, 1] <= right_value
        )
        if (length(tail_rows) >= 8L && length(adjustable)) {
          tail_x <- observed_full_data[tail_rows, 1]
          tail_y <- pmax(observed_full_data[tail_rows, 2], 0)
          tail_scale <- max(tail_y, 1)
          tail_weights <- 0.25 + 0.75 * sqrt(tail_y / tail_scale)
          alpha_reference <- pmax(alpha[adjustable], 1e-10)
          # Preserve the data-selected first-peak height during any later
          # component optimisation; do not return to a fixed 99.5% target.
          first_height_target <- if (
            is.finite(prefix_peak_height_ratio_after)
          ) prefix_peak_height_ratio_after * observed_first_height else
            observed_first_height

          evaluate_tail_candidate <- function(alpha_value, sigma_value) {
            # Recalibrate alpha1 after every later-weight update so smoothing
            # cannot move the recovered first peak above the grey curve.
            peak_basis <- stats::dnorm(
              first_centre, miu, sigma_value
            ) * sum(y_fit_processed) * prefix_bin_width
            other_peak_height <- sum(
              alpha_value[-first_component] *
                peak_basis[-first_component]
            )
            alpha_value[first_component] <- max(
              0, first_height_target - other_peak_height
            ) / max(peak_basis[first_component], .Machine$double.eps)
            design <- vapply(seq_along(alpha_value), function(j) {
              stats::dnorm(tail_x, miu[j], sigma_value[j]) *
                sum(y_fit_processed) * prefix_bin_width
            }, numeric(length(tail_x)))
            prediction <- if (is.null(dim(design))) {
              design * alpha_value
            } else as.numeric(design %*% alpha_value)
            residual <- (prediction - tail_y) / tail_scale
            nrmse_value <- sqrt(
              sum(tail_weights * residual^2) / sum(tail_weights)
            )
            over_value <- sum(pmax(prediction - tail_y, 0)) /
              max(sum(tail_y), 1)
            under_value <- sum(pmax(tail_y - prediction, 0)) /
              max(sum(tail_y), 1)
            ripple_value <- if (length(residual) >= 3L) {
              sqrt(mean(diff(residual, differences = 2L)^2))
            } else 0
            score_value <- 0.65 * nrmse_value +
              0.15 * over_value + 0.10 * under_value +
              0.10 * ripple_value
            peak_ratio <- sum(alpha_value * peak_basis) /
              max(observed_first_height, 1)
            list(
              alpha = alpha_value, prediction = prediction,
              nrmse = nrmse_value, over = over_value,
              under = under_value, ripple = ripple_value,
              score = score_value, peak_ratio = peak_ratio
            )
          }

          baseline_tail <- evaluate_tail_candidate(alpha, sigma)
          baseline_tail_sigma <- sigma
          prefix_tail_nrmse_before <- baseline_tail$nrmse
          prefix_tail_ripple_before <- baseline_tail$ripple
          # Parameterise each later sigma continuously through its log-ratio
          # to the preceding component. tanh keeps every adjacent ratio inside
          # [1/1.35, 1.35] without imposing one shared scale or a fixed search
          # range. Each later peak may independently widen or narrow.
          later_order <- all_order[-1L]
          previous_order <- all_order[-length(all_order)]
          maximum_log_ratio <- log(1.35)
          initial_relative_change <- log(
            sigma[later_order] / sigma[previous_order]
          ) / maximum_log_ratio
          initial_width_parameters <- atanh(pmin(
            0.999, pmax(-0.999, initial_relative_change)
          ))
          alpha_parameter_count <- length(adjustable)
          decode_tail_parameters <- function(parameters) {
            candidate_alpha <- alpha
            candidate_alpha[adjustable] <- exp(
              parameters[seq_len(alpha_parameter_count)]
            )
            width_parameters <- parameters[
              alpha_parameter_count + seq_along(later_order)
            ]
            candidate_sigma <- sigma
            candidate_sigma[first_component] <- sigma[first_component]
            for (position in seq_along(later_order)) {
              previous <- all_order[position]
              current <- all_order[position + 1L]
              candidate_sigma[current] <-
                candidate_sigma[previous] * exp(
                  maximum_log_ratio * tanh(width_parameters[position])
                )
            }
            list(alpha = candidate_alpha, sigma = candidate_sigma)
          }
          initial_parameters <- c(
            log(alpha_reference), initial_width_parameters
          )
          objective <- function(parameters) {
            decoded <- decode_tail_parameters(parameters)
            metric <- evaluate_tail_candidate(
              decoded$alpha, decoded$sigma
            )
            alpha_change <- mean(
              (parameters[seq_len(alpha_parameter_count)] -
                 log(alpha_reference))^2
            )
            width_change <- mean(
              log(decoded$sigma[adjustable] /
                    pmax(sigma[adjustable], .Machine$double.eps))^2
            )
            biological_log_sigma <- log(decoded$sigma[biological])
            width_curvature <- if (length(biological_log_sigma) >= 3L) {
              mean(diff(biological_log_sigma, differences = 2L)^2)
            } else 0
            metric$score + 0.002 * alpha_change +
              0.004 * width_change + 0.004 * width_curvature
          }
          lower_parameters <- c(
            log(alpha_reference * 0.10),
            rep(-4, length(initial_width_parameters))
          )
          upper_parameters <- c(
            log(alpha_reference * 10),
            rep(4, length(initial_width_parameters))
          )
          optimized <- tryCatch(
            stats::optim(
              initial_parameters, objective, method = "L-BFGS-B",
              lower = lower_parameters, upper = upper_parameters,
              control = list(maxit = 400L, factr = 1e8)
            ),
            error = function(error) NULL
          )
          selected_tail <- NULL
          if (!is.null(optimized) && is.finite(optimized$value)) {
            decoded <- decode_tail_parameters(optimized$par)
            metric <- evaluate_tail_candidate(
              decoded$alpha, decoded$sigma
            )
            candidate_alpha <- metric$alpha
            candidate_sigma <- decoded$sigma
            shares <- pmax(candidate_alpha, 0) /
              max(sum(pmax(candidate_alpha, 0)), .Machine$double.eps)
            ordered_sigma <- candidate_sigma[biological]
            adjacent_ratio <- pmax(
              ordered_sigma[-1L] /
                ordered_sigma[-length(ordered_sigma)],
              ordered_sigma[-length(ordered_sigma)] /
                ordered_sigma[-1L]
            )
            whole_metric <- evaluate_prefix_flank(
              candidate_alpha, candidate_sigma
            )
            width_changed <- any(abs(log(
              candidate_sigma[adjustable] /
                pmax(sigma[adjustable], .Machine$double.eps)
            )) >= 0.005)
            eligible <- isTRUE(width_changed) &&
              metric$nrmse <= baseline_tail$nrmse * 0.97 &&
              metric$score <= baseline_tail$score * 0.97 &&
              metric$ripple <= baseline_tail$ripple * 0.95 + 1e-6 &&
              metric$over <= baseline_tail$over + 0.005 &&
              metric$under <= baseline_tail$under + 0.005 &&
              metric$peak_ratio >= 0.95 &&
              metric$peak_ratio <= 1.02 &&
              all(shares[biological] >= 0.01) &&
              all(adjacent_ratio <= 1.35 + 1e-8) &&
              whole_metric$nrmse <= baseline_flank$nrmse * 1.005 + 1e-4
            if (isTRUE(eligible)) {
              selected_tail <- list(
                alpha = candidate_alpha, sigma = candidate_sigma,
                metric = metric, whole_metric = whole_metric
              )
            }
          }
          if (!is.null(selected_tail)) {
            alpha <- selected_tail$alpha
            sigma <- selected_tail$sigma
            EM_result$alpha <- alpha
            EM_result$sigma <- sigma
            prefix_tail_smoothing_accepted <- TRUE
            prefix_tail_sigma_before <- baseline_tail_sigma
            prefix_tail_sigma_after <- selected_tail$sigma
            prefix_tail_nrmse_after <- selected_tail$metric$nrmse
            prefix_tail_ripple_after <- selected_tail$metric$ripple
            prefix_peak_height_ratio_after <-
              selected_tail$metric$peak_ratio
          } else {
            prefix_tail_nrmse_after <- prefix_tail_nrmse_before
            prefix_tail_ripple_after <- prefix_tail_ripple_before
          }
        }
        cat("Recovered-prefix later-peak smoothing:",
            if (prefix_tail_smoothing_accepted) "accepted" else
              "not needed/rejected", "\n")
        if (prefix_tail_smoothing_accepted) {
          cat("Later sigmas after continuous optimization:",
              paste(round(prefix_tail_sigma_after[all_order[-1L]], 4),
                    collapse = " "), "\n")
          cat("Tail NRMSE before/after:",
              round(prefix_tail_nrmse_before, 6), "/",
              round(prefix_tail_nrmse_after, 6), "\n")
          cat("Tail ripple before/after:",
              round(prefix_tail_ripple_before, 6), "/",
              round(prefix_tail_ripple_after, 6), "\n")
        }
      }

      # Keep the two fitting situations strictly separated:
      #
      # 1. An observed first peak is fitted forwards by the ordinary EM/model
      #    selection pipeline.  Its selected parameters must not be rewritten
      #    by the hidden-prefix post-fit repairs below.
      # 2. Only when missing prefix peaks were explicitly reconstructed may
      #    their estimated areas, widths and centres receive the additional
      #    guarded post-fit refinements.
      #
      # Model/ploidy selection and the optional n+1 decision remain unchanged.
      selected_shape_refinement <- list(
        accepted = FALSE, triggered = FALSE, reason = "not_final_fit"
      )
      estimated_first_height_refinement <- list(
        accepted = FALSE, triggered = FALSE, reason = "not_final_fit"
      )
      first_second_bridge_refinement <- list(
        accepted = FALSE, triggered = FALSE, reason = "not_final_fit"
      )
      final_inferred_first_refinement <- list(
        accepted = FALSE, triggered = FALSE, reason = "not_hidden_prefix_fit"
      )
      if (isTRUE(save_plot) && !isTRUE(.prefix_discovery_only) &&
          isTRUE(prefix_completion_active)) {
        sigma_score_x <- x_fit
        sigma_score_y <- y_fit
        sigma_spacing <- first_peak_x
        first_sigma_range <- NULL
        if (isTRUE(prefix_completion_active)) {
          first_shape_component <- order(miu)[1L]
          shape_score_start <- max(
            observed_full_data[prefix_noise_start_index, 1],
            miu[first_shape_component] - 0.05 * prefix_first_peak_x
          )
          sigma_score_rows <- which(
            observed_full_data[, 1] >= shape_score_start &
              observed_full_data[, 1] <= right_value
          )
          sigma_score_x <- observed_full_data[sigma_score_rows, 1]
          sigma_score_y <- observed_full_data[sigma_score_rows, 2]
          sigma_spacing <- prefix_first_peak_x
          # Prefix completion estimates the missing first harmonic, so retain
          # its existing biological width interval while allowing every later
          # width to be selected continuously from the observed curve.
          first_sigma_range <- c(0.35, 0.40)
        }
        selected_shape_refinement <- .GenomeGMM_refine_selected_shape(
          x = sigma_score_x,
          y = sigma_score_y,
          alpha = alpha,
          miu = miu,
          sigma = sigma,
          total_samples_fit = sum(y_fit_processed),
          first_peak_spacing = sigma_spacing,
          sigma_reliability = EM_result$sigma_reliability,
          biological_count = n,
          first_sigma_fraction_range = first_sigma_range
        )
        if (isTRUE(selected_shape_refinement$accepted)) {
          alpha <- selected_shape_refinement$alpha
          miu <- selected_shape_refinement$miu
          sigma <- selected_shape_refinement$sigma
          EM_result$alpha <- alpha
          EM_result$miu <- miu
          EM_result$sigma <- sigma
          if (isTRUE(prefix_completion_active)) {
            biological <- order(miu)[seq_len(min(n, length(miu)))]
            first_component <- biological[1L]
            first_centre <- miu[first_component]
            prefix_peak_height_ratio_after <- sum(
              alpha * stats::dnorm(first_centre, miu, sigma)
            ) * sum(y_fit_processed) * prefix_bin_width /
              max(observed_first_height, 1)
          }
        }
        cat("Final guarded shape refinement:",
            if (isTRUE(selected_shape_refinement$accepted)) {
              "accepted"
            } else {
              paste0("kept selected parameters (",
                     selected_shape_refinement$reason, ")")
            }, "\n")
        selected_height_refinement <-
          .GenomeGMM_refine_biological_heights(
            x = sigma_score_x,
            y = sigma_score_y,
            alpha = alpha,
            miu = miu,
            sigma = sigma,
            total_samples_fit = sum(y_fit_processed),
            biological_count = n
          )
        if (isTRUE(selected_height_refinement$accepted)) {
          alpha <- selected_height_refinement$alpha
          EM_result$alpha <- alpha
        }
        cat("Biological peak-height continuity:",
            if (isTRUE(selected_height_refinement$accepted)) {
              "accepted"
            } else {
              paste0("kept shape-refined parameters (",
                     selected_height_refinement$reason, ")")
            }, "\n")

        # Apply this after every other final shape/height adjustment so the
        # value used by genome-size calculation, the PDF and the legend is the
        # actually validated first-peak amplitude.  The trigger is the observed
        # red/grey mismatch itself, not which upstream discovery branch won.
        estimated_first_height_refinement <-
          .GenomeGMM_refine_estimated_first_height_final(
            x = observed_full_data[, 1],
            y = observed_full_data[, 2],
            alpha = alpha,
            miu = miu,
            sigma = sigma,
            total_samples_fit = sum(y_fit_processed),
            biological_count = n,
            # Only a reconstructed first harmonic has an estimated amplitude
            # that needs this extra post-fit correction.  For an explicitly
            # observed first peak, the selected EM/joint fit already estimates
            # its mass from data; applying the prefix correction again can
            # systematically press any genuine observed first peak downward.
            apply_guard = isTRUE(prefix_completion_active)
          )
        if (isTRUE(estimated_first_height_refinement$accepted)) {
          alpha <- estimated_first_height_refinement$alpha
          EM_result$alpha <- alpha
          if (isTRUE(prefix_completion_active)) {
            biological <- order(miu)[seq_len(min(n, length(miu)))]
            first_component <- biological[1L]
            first_centre <- miu[first_component]
            prefix_peak_height_ratio_after <- sum(
              alpha * stats::dnorm(first_centre, miu, sigma)
            ) * sum(y_fit_processed) * prefix_bin_width /
              max(observed_first_height, 1)
          }
        }
        cat("Final first-peak amplitude guard:",
            if (isTRUE(estimated_first_height_refinement$accepted)) {
              paste0(
                estimated_first_height_refinement$direction,
                " (scale ",
                format(estimated_first_height_refinement$scale, digits = 4),
                ", centre ratio ",
                format(
                  estimated_first_height_refinement$before$centre_ratio,
                  digits = 4
                ),
                " -> ",
                format(
                  estimated_first_height_refinement$after$centre_ratio,
                  digits = 4
                ),
                ")"
              )
            } else {
              paste0(
                "kept current height (",
                estimated_first_height_refinement$reason,
                ")"
              )
            }, "\n")

        first_second_bridge_refinement <-
          .GenomeGMM_refine_first_second_bridge(
            x = observed_full_data[, 1],
            y = observed_full_data[, 2],
            alpha = alpha,
            miu = miu,
            sigma = sigma,
            total_samples_fit = sum(y_fit_processed),
            biological_count = n,
            first_sigma_upper = Inf
          )
        if (isTRUE(first_second_bridge_refinement$accepted)) {
          alpha <- first_second_bridge_refinement$alpha
          miu <- first_second_bridge_refinement$miu
          sigma <- first_second_bridge_refinement$sigma
          EM_result$alpha <- alpha
          EM_result$miu <- miu
          EM_result$sigma <- sigma
          if (isTRUE(prefix_completion_active)) {
            biological <- order(miu)[seq_len(min(n, length(miu)))]
            first_component <- biological[1L]
            first_centre <- miu[first_component]
            prefix_peak_height_ratio_after <- sum(
              alpha * stats::dnorm(first_centre, miu, sigma)
            ) * sum(y_fit_processed) * prefix_bin_width /
              max(observed_first_height, 1)
          }
        }
        cat("First/second peak bridge refinement:",
            if (isTRUE(first_second_bridge_refinement$accepted)) {
              paste0(
                "accepted (sigma ",
                paste(
                  format(
                    first_second_bridge_refinement$sigma_before,
                    digits = 4
                  ), collapse = "/"
                ),
                " -> ",
                paste(
                  format(
                    first_second_bridge_refinement$sigma_after,
                    digits = 4
                  ), collapse = "/"
                ),
                ")"
              )
            } else {
              paste0(
                "kept current parameters (",
                first_second_bridge_refinement$reason,
                ")"
              )
            }, "\n")
        bridge_diagnostics_numeric <-
          is.list(first_second_bridge_refinement$before) &&
          is.list(first_second_bridge_refinement$after) &&
          is.numeric(first_second_bridge_refinement$before$interval_over) &&
          is.numeric(first_second_bridge_refinement$after$interval_over) &&
          is.numeric(first_second_bridge_refinement$after$centre_ratio) &&
          is.numeric(first_second_bridge_refinement$before$interval_mse) &&
          is.numeric(first_second_bridge_refinement$after$interval_mse) &&
          is.numeric(first_second_bridge_refinement$after$interval_max_under) &&
          is.numeric(first_second_bridge_refinement$before$local_mse) &&
          is.numeric(first_second_bridge_refinement$after$local_mse) &&
          is.numeric(first_second_bridge_refinement$before$nrmse) &&
          is.numeric(first_second_bridge_refinement$after$nrmse)
        if (isTRUE(first_second_bridge_refinement$triggered) &&
            isTRUE(bridge_diagnostics_numeric)) {
          cat("Bridge direction:",
              first_second_bridge_refinement$direction, "\n")
          cat("Bridge interval overfit before/after:",
              round(first_second_bridge_refinement$before$interval_over, 5),
              "/",
              round(first_second_bridge_refinement$after$interval_over, 5),
              "\n")
          cat("Bridge centre ratios after:",
              paste(round(
                first_second_bridge_refinement$after$centre_ratio, 4
              ), collapse = " "), "\n")
          cat("Bridge interval MSE before/after:",
              round(first_second_bridge_refinement$before$interval_mse, 6),
              "/",
              round(first_second_bridge_refinement$after$interval_mse, 6),
              "\n")
          cat("Bridge maximum underfit after:",
              round(
                first_second_bridge_refinement$after$interval_max_under, 5
              ), "\n")
          cat("Bridge local MSE before/after:",
              round(first_second_bridge_refinement$before$local_mse, 6),
              "/",
              round(first_second_bridge_refinement$after$local_mse, 6),
              "\n")
          cat("Bridge NRMSE before/after:",
              round(first_second_bridge_refinement$before$nrmse, 6),
              "/",
              round(first_second_bridge_refinement$after$nrmse, 6),
              "\n")
        }

        # The first component of a hidden-prefix fit is reconstructed rather
        # than directly observed.  Its area, centre and width therefore remain
        # estimable parameters.  Optimise all three together against exactly
        # the same smooth grey reference drawn in the PDF.  This pass is only
        # applied to the selected hidden-prefix fit; ordinary fits and ploidy
        # selection are deliberately left unchanged.
        if (isTRUE(prefix_completion_active)) {
          final_inferred_first_refinement <-
            .GenomeGMM_refine_estimated_first_shape(
              x = observed_full_data[, 1],
              y = observed_full_data[, 2],
              alpha = alpha,
              miu = miu,
              sigma = sigma,
              count_scale = sum(y_fit_processed) * prefix_bin_width,
              spacing = prefix_first_peak_x,
              biological_count = n,
              max_iterations = 320L,
              use_pdf_reference = TRUE,
              final_selection = TRUE
            )
          if (isTRUE(final_inferred_first_refinement$accepted)) {
            alpha <- final_inferred_first_refinement$alpha
            miu <- final_inferred_first_refinement$miu
            sigma <- final_inferred_first_refinement$sigma
            EM_result$alpha <- alpha
            EM_result$miu <- miu
            EM_result$sigma <- sigma
          }

          # Re-evaluate the peak-height diagnostic at the possibly adjusted
          # centre using the plotted grey reference, not a stale raw-bin value.
          biological <- order(miu)[seq_len(min(n, length(miu)))]
          first_component <- biological[1L]
          first_centre <- miu[first_component]
          grey_reference <- make_gmm_reference_curve(
            observed_full_data[-1L, 1],
            observed_full_data[-1L, 2],
            dense_points = 1000L
          )
          current_observed_first_height <- stats::approx(
            grey_reference$x, pmax(grey_reference$y, 0),
            xout = first_centre, rule = 2
          )$y
          prefix_peak_height_ratio_after <- sum(
            alpha * stats::dnorm(first_centre, miu, sigma)
          ) * sum(y_fit_processed) * prefix_bin_width /
            max(current_observed_first_height, 1)

          cat("Final inferred first-peak joint refinement:",
              if (isTRUE(final_inferred_first_refinement$accepted)) {
                paste0(
                  "accepted (alpha x",
                  format(final_inferred_first_refinement$alpha_ratio,
                         digits = 4),
                  ", miu shift ",
                  format(final_inferred_first_refinement$first_mean_shift,
                         digits = 4),
                  ", sigma x",
                  format(final_inferred_first_refinement$first_sigma_ratio,
                         digits = 4), ")"
                )
              } else {
                paste0("kept current parameters (",
                       final_inferred_first_refinement$reason, ")")
              }, "\n")
        }
      }
      if (isTRUE(prefix_completion_active) &&
          !isTRUE(.prefix_discovery_only)) {
        biological <- order(miu)[seq_len(min(n, length(miu)))]
        biological_sigma <- sigma[biological]
        # Prefix recovery is responsible only for making the inferred first
        # Gaussian identifiable. Later widths have already been selected by
        # the ordinary common/trend/adaptive fitter and its fit-quality guards;
        # do not reject the entire fit here because a later refinement moves
        # an adjacent pair slightly beyond a separate hard ratio threshold.
        ordinary_hidden_width_ok <-
          biological_sigma[1] >= 0.35 * prefix_first_peak_x - 1e-8 &&
          biological_sigma[1] <= 0.40 * prefix_first_peak_x + 1e-8
        bridge_supported_width_ok <-
          (isTRUE(first_second_bridge_refinement$accepted) ||
           isTRUE(final_inferred_first_refinement$accepted)) &&
          biological_sigma[1] >= 0.35 * prefix_first_peak_x - 1e-8 &&
          biological_sigma[1] <= 0.50 * prefix_first_peak_x + 1e-8 &&
          biological_sigma[1] <= biological_sigma[2] + 1e-8 &&
          biological_sigma[2] <= biological_sigma[3] + 1e-8
        hidden_width_ok <- isTRUE(ordinary_hidden_width_ok) ||
          isTRUE(bridge_supported_width_ok)
        hidden_peak_height_ok <- is.finite(prefix_peak_height_ratio_after) &&
          prefix_peak_height_ratio_after >= 0.86 &&
          prefix_peak_height_ratio_after <= 1.08
        interval_rows <- which(
          observed_full_data[, 1] >= miu[biological[1]] &
            observed_full_data[, 1] <= miu[biological[2]]
        )
        hidden_overshoot_ok <- length(interval_rows) >= 3L
        prefix_interval_underfit <- NA_real_
        if (hidden_overshoot_ok) {
          interval_x <- observed_full_data[interval_rows, 1]
          interval_y <- pmax(observed_full_data[interval_rows, 2], 0)
          evaluate_prefix_interval <- function(alpha_value) {
            interval_components <- vapply(seq_along(alpha_value), function(j) {
              alpha_value[j] *
                stats::dnorm(interval_x, miu[j], sigma[j]) *
                sum(y_fit_processed) * prefix_bin_width
            }, numeric(length(interval_x)))
            interval_prediction <- if (is.null(dim(interval_components))) {
              interval_components
            } else rowSums(interval_components)
            relative_overshoot <- (interval_prediction - interval_y) /
              pmax(interval_y, 0.05 * max(interval_y), 1)
            overshoot_runs <- rle(relative_overshoot > 0.08)
            list(
              prediction = interval_prediction,
              overshoot = sum(pmax(
                interval_prediction - interval_y, 0
              )) / max(sum(interval_y), 1),
              underfit = sum(pmax(
                interval_y - interval_prediction, 0
              )) / max(sum(interval_y), 1),
              run = if (any(overshoot_runs$values)) {
                max(overshoot_runs$lengths[overshoot_runs$values])
              } else 0L
            )
          }
          interval_metric <- evaluate_prefix_interval(alpha)
          prefix_interval_overshoot <- interval_metric$overshoot
          prefix_interval_underfit <- interval_metric$underfit
          prefix_interval_overshoot_run <- interval_metric$run
          hidden_overshoot_ok <-
            prefix_interval_overshoot <= 0.05 &&
            prefix_interval_overshoot_run < 3L

          # If the bidirectional bridge stage has correctly identified a
          # sustained overfit but a few raw bins remain above the final guard,
          # continue lowering the first two component masses by the smallest
          # possible common factor. Peak geometry remains unchanged and a new
          # continuous underfit is not allowed to replace the overshoot.
          needs_interval_height_repair <-
            isTRUE(first_second_bridge_refinement$accepted) &&
            identical(first_second_bridge_refinement$direction, "overfit") &&
            (!isTRUE(hidden_peak_height_ok) ||
             !isTRUE(hidden_overshoot_ok))
          if (isTRUE(needs_interval_height_repair)) {
            for (height_scale in seq(0.995, 0.75, by = -0.005)) {
              candidate_alpha <- alpha
              candidate_alpha[biological[1:2]] <-
                candidate_alpha[biological[1:2]] * height_scale
              candidate_peak_ratio <- sum(
                candidate_alpha * stats::dnorm(
                  miu[biological[1]], miu, sigma
                )
              ) * sum(y_fit_processed) * prefix_bin_width /
                max(observed_first_height, 1)
              candidate_interval <-
                evaluate_prefix_interval(candidate_alpha)
              candidate_ok <- is.finite(candidate_peak_ratio) &&
                candidate_peak_ratio >= 0.86 &&
                candidate_peak_ratio <= 1.08 &&
                candidate_interval$overshoot <= 0.05 &&
                candidate_interval$run < 3L &&
                candidate_interval$underfit <= 0.08
              if (isTRUE(candidate_ok)) {
                alpha <- candidate_alpha
                EM_result$alpha <- alpha
                prefix_interval_height_scale <- height_scale
                prefix_peak_height_ratio_after <- candidate_peak_ratio
                prefix_interval_overshoot <-
                  candidate_interval$overshoot
                prefix_interval_underfit <- candidate_interval$underfit
                prefix_interval_overshoot_run <- candidate_interval$run
                hidden_peak_height_ok <- TRUE
                hidden_overshoot_ok <- TRUE
                cat("Prefix interval height repair: accepted (scale",
                    format(height_scale, digits = 4), ")\n")
                break
              }
            }
          }
        }
        # These values remain diagnostics.  The inferred peak has already been
        # jointly optimised and guarded by full local red/grey error; a stale
        # point-height or run-length diagnostic must not abort an otherwise
        # improved fit.
        hidden_prefix_guard_passed <-
          isTRUE(hidden_width_ok) && isTRUE(hidden_peak_height_ok) &&
          isTRUE(hidden_overshoot_ok)
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
      
      # Calculate genome size. The highest requested biological harmonic is the
      # homozygous/full-copy depth. Any optional n+1 tail component is deliberately
      # excluded from this choice.
      reporting_ploidy <- getOption("GenomeGMM.reporting_ploidy", n)
      if (!is.numeric(reporting_ploidy) || length(reporting_ploidy) != 1L ||
          !is.finite(reporting_ploidy) || reporting_ploidy < 1L) {
        reporting_ploidy <- n
      }
      reporting_ploidy <- as.integer(reporting_ploidy)
      genome_size_result <- calculate_genome_size(
        x_fit_processed, y_fit_processed, alpha, miu, sigma, 
        full_data, start_index, right_value, 
        biological_component_count = n,
        ploidy = reporting_ploidy
      )
      
      total_kmers <- genome_size_result$total_kmers
      kmer_depth <- genome_size_result$kmer_depth
      g_hap_bp <- genome_size_result$g_hap_bp
      genome_size_bp <- genome_size_result$genome_size_bp
      component_sizes_bp <- genome_size_result$component_sizes_bp # Get component sizes
    
      # Score the fitted curve against the ORIGINAL histogram.  This deliberately
      # does not score against the symmetry-completed data, otherwise an artificial
      # main peak could look better than the observation it is meant to explain.
      score_x <- x_fit
      score_y <- y_fit
      if (isTRUE(prefix_completion_active)) {
        score_rows <- which(
          observed_full_data[, 1] >=
            observed_full_data[prefix_noise_start_index, 1] &
            observed_full_data[, 1] <= right_value
        )
        score_x <- observed_full_data[score_rows, 1]
        score_y <- observed_full_data[score_rows, 2]
      }
      fit_score_info <- calculate_gmm_fit_score(
        x_observed = score_x,
        y_observed = score_y,
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
      cat("Homozygous peak depth:", kmer_depth, "\n")
      cat("g-hap size:", g_hap_bp, "bp\n")
      
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
          if (isTRUE(prefix_completion_active)) score_x else x_fit,
          if (isTRUE(prefix_completion_active)) score_y else y_fit,
          alpha, miu, sigma, genome_size_bp, component_sizes_bp, # Pass component sizes
          main_peak_x_symmetric, main_peak_y_symmetric,
          x_fit_processed, y_fit_processed,
          if (isTRUE(prefix_completion_active)) observed_full_data else full_data,
          right_value,
          fit_extra_component = fit_extra_component,  # Pass extra parameter
          internal_n = internal_n,                     # Pass internal fitted count
          first_peak_x = first_peak_x,
          n = n,
          monoploid_divisor = getOption(
            "GenomeGMM.reporting_monoploid_divisor", reporting_ploidy
          ),
          species_name = species_name,
          prefix_completion = prefix_completion_active,
          prefix_noise_start_x = if (isTRUE(prefix_completion_active)) {
            observed_full_data[prefix_noise_start_index, 1]
          } else NULL
        )
        dev.off()
        cat("Plot saved as:", filename, "\n")
      }
      cat("Automatic fit score (lower is better):", round(fit_score_info$score, 6), "\n")
      cat("Algorithm execution completed!\n")
      
      return(invisible(list(
        genome_size_bp = genome_size_bp,
        haploid_genome_size_bp = g_hap_bp,
        g_hap_bp = g_hap_bp,
        total_kmers = total_kmers,
        total_kmer_red_mass = genome_size_result$total_kmer_red_mass,
        total_kmer_grey_mass = genome_size_result$total_kmer_grey_mass,
        total_kmer_switch_x = genome_size_result$total_kmer_switch_x,
        total_kmer_agreement_tolerance =
          genome_size_result$total_kmer_agreement_tolerance,
        total_kmer_agreement_bins =
          genome_size_result$total_kmer_agreement_bins,
        kmer_depth = kmer_depth,
        homozygous_peak_depth = kmer_depth,
        homozygous_peak_component = genome_size_result$homozygous_peak_component,
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
        final_shape_refinement = selected_shape_refinement,
        final_first_peak_amplitude_guard =
          estimated_first_height_refinement,
        first_second_bridge_refinement =
          first_second_bridge_refinement,
        final_inferred_first_refinement =
          final_inferred_first_refinement,
        curve_refined = isTRUE(curve_refinement$refined),
        curve_scale = curve_refinement$curve_scale,
        symmetry_tolerance = symmetry_tolerance,
        extra_component_mode = extra_component_mode,
        n = n,
        main_peak_ratio = main_peak_ratio,
        fit_score = fit_score_info$score,
        fit_metrics = fit_score_info,
        # Internal-only values used to compare a possible hidden-prefix fit on
        # exactly the same observed bins. They are not user tuning parameters.
        .diagnostic_x_fit = x_fit,
        .diagnostic_y_fit = y_fit,
        .diagnostic_total_samples_fit = sum(y_fit_processed),
        .diagnostic_bin_width = if (length(unique(x_fit)) > 1L) {
          stats::median(diff(sort(unique(x_fit))))
        } else 1,
        prefix_completion = isTRUE(prefix_completion_active),
        prefix_completion_slot = if (isTRUE(prefix_completion_active)) {
          as.integer(.prefix_completion_slot)
        } else NA_integer_,
        prefix_flank_refined = isTRUE(prefix_flank_refined),
        prefix_residual_refinement = prefix_residual_refinement,
        prefix_first_height_refinement = prefix_first_height_refinement,
        prefix_first_height_refinement_after_shape =
          prefix_first_height_refinement_after_shape,
        prefix_first_shape_refinement = prefix_first_shape_refinement,
        prefix_weight_refinement = prefix_weight_refinement,
        interpeak_shape_refinement = interpeak_shape_refinement,
        prefix_flank_scale = prefix_flank_scale,
        prefix_flank_overshoot_before = prefix_flank_before,
        prefix_flank_overshoot_after = prefix_flank_after,
        prefix_interval_overshoot = prefix_interval_overshoot,
        prefix_interval_overshoot_run = prefix_interval_overshoot_run,
        prefix_interval_height_scale = prefix_interval_height_scale,
        prefix_peak_height_ratio_before = prefix_peak_height_ratio_before,
        prefix_peak_height_ratio_after = prefix_peak_height_ratio_after,
        prefix_tail_smoothing_accepted =
          isTRUE(prefix_tail_smoothing_accepted),
        prefix_tail_sigma_before = prefix_tail_sigma_before,
        prefix_tail_sigma_after = prefix_tail_sigma_after,
        prefix_tail_nrmse_before = prefix_tail_nrmse_before,
        prefix_tail_nrmse_after = prefix_tail_nrmse_after,
        prefix_tail_ripple_before = prefix_tail_ripple_before,
        prefix_tail_ripple_after = prefix_tail_ripple_after,
        prefix_noise_start_x = if (isTRUE(prefix_completion_active)) {
          observed_full_data[prefix_noise_start_index, 1]
        } else NA_real_,
        prefix_sigma_seed = prefix_sigma_seed,
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
                                        maximum_sigma_fraction = 0.80,
                                        first_sigma_fraction_range = NULL,
                                        max_adjacent_sigma_ratio = Inf,
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
      maximum_sigma_fraction <- pmin(
        0.80, pmax(0.20, as.numeric(maximum_sigma_fraction))
      )
      maximum_sigma <- spacing * maximum_sigma_fraction
      minimum_sigma_by_component <- rep(minimum_sigma, component_count)
      maximum_sigma_by_component <- rep(maximum_sigma, component_count)
      if (!is.null(first_sigma_fraction_range) &&
          length(first_sigma_fraction_range) == 2L &&
          all(is.finite(first_sigma_fraction_range))) {
        first_range <- sort(as.numeric(first_sigma_fraction_range))
        minimum_sigma_by_component[1] <- max(
          minimum_sigma, first_range[1] * spacing
        )
        maximum_sigma_by_component[1] <- min(
          maximum_sigma, first_range[2] * spacing
        )
      }
      sigma <- pmin(maximum_sigma_by_component,
                    pmax(minimum_sigma_by_component, sigma))
      if (is.finite(max_adjacent_sigma_ratio) && component_count > 1L) {
        for (j in 2:component_count) {
          sigma[j] <- min(
            sigma[j], sigma[j - 1L] * max_adjacent_sigma_ratio
          )
          sigma[j] <- max(
            sigma[j], sigma[j - 1L] / max_adjacent_sigma_ratio
          )
        }
      }
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
            any(fitted_sigma < minimum_sigma_by_component) ||
            any(fitted_sigma > maximum_sigma_by_component)) return(1e6)
        if (is.finite(max_adjacent_sigma_ratio) && component_count > 1L) {
          adjacent_ratio <- pmax(
            fitted_sigma[-1L] / fitted_sigma[-component_count],
            fitted_sigma[-component_count] / fitted_sigma[-1L]
          )
          if (any(adjacent_ratio > max_adjacent_sigma_ratio + 1e-8)) {
            return(1e6)
          }
        }
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
        width_lower <- log(max(minimum_sigma_by_component))
        width_upper <- log(min(maximum_sigma_by_component))
      } else if (identical(sigma_mode, "trend")) {
        maximum_slope <- if (component_count > 1) {
          log(maximum_sigma / minimum_sigma) / log(component_count)
        } else 0
        width_lower <- c(log(minimum_sigma), -maximum_slope)
        width_upper <- c(log(maximum_sigma), maximum_slope)
      } else {
        width_lower <- log(minimum_sigma_by_component)
        width_upper <- log(maximum_sigma_by_component)
      }
      if (length(width_start)) {
        width_start <- pmin(width_upper, pmax(width_lower, width_start))
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

      calibration_guard <- .GenomeGMM_check_main_peak_calibration(
        miu = miu,
        sigma = sigma,
        main_peak_component = main_peak_component,
        calibrated_mean = miu_calibrated[main_peak_component],
        calibrated_sd = sigma_calibrated[main_peak_component],
        first_peak_x = first_peak_x,
        max_sigma_ratio = max_sigma_ratio
      )
      if (!isTRUE(calibration_guard$accepted)) {
        cat(
          "Main peak calibration rejected by biological position/width guard:",
          calibration_guard$reason, "\n"
        )
        return(list(
          alpha = alpha, miu = miu, sigma = sigma,
          calibration_accepted = FALSE,
          calibration_guard = calibration_guard
        ))
      }
      miu_calibrated[main_peak_component] <- calibration_guard$projected_mean
      sigma_calibrated[main_peak_component] <- calibration_guard$projected_sd
      cat(
        "Main peak calibration", calibration_guard$reason,
        sprintf("(miu=%.4f, sigma=%.4f)\n",
                miu_calibrated[main_peak_component],
                sigma_calibrated[main_peak_component])
      )
      
      # Adjust mixing coefficients to maintain total probability of 1
      total_alpha <- sum(alpha_calibrated)
      alpha_calibrated <- alpha_calibrated / total_alpha
      
      cat("=== Main peak calibration completed ===\n")
      
      return(list(
        alpha = alpha_calibrated,
        miu = miu_calibrated,
        sigma = sigma_calibrated,
        calibration_accepted = TRUE,
        calibration_guard = calibration_guard
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
    
    #' Calculate g-hap from the highest biological harmonic
    calculate_genome_size <- function(x_fit_processed, y_fit_processed, alpha, miu, sigma,
                                      full_data, start_index, right_value,
                                      biological_component_count = length(miu),
                                      ploidy = biological_component_count) {
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
      
      biological_component_count <- min(
        max(1L, as.integer(biological_component_count)), k
      )
      biological_components <- seq_len(biological_component_count)
      valid_components <- biological_components[
        is.finite(miu[biological_components]) & miu[biological_components] > 0
      ]
      if (!length(valid_components)) {
        stop("No valid biological component is available for genome-size estimation.")
      }
      homozygous_peak_component <- valid_components[
        which.max(miu[valid_components])
      ]
      kmer_depth <- as.numeric(miu[homozygous_peak_component])
      ploidy <- max(1L, as.integer(ploidy))
      
      cat("Homozygous peak position (kmer_depth):", kmer_depth, "\n")
      cat("Homozygous peak component:", homozygous_peak_component, "\n")
      cat("All component means:", round(miu, 2), "\n")
      
      # For data within fitting range, use GMM estimated distribution
      fit_range_upper <- max(x_fit_processed)
      natural_x <- seq(min(x_fit_processed), fit_range_upper, by = 1)
      density_natural_x <- approx(x_vals, total_density, xout = natural_x)$y
      
      # Convert probability density to counts
      frequency_natural_x <- density_natural_x * total_samples_filtered * bin_width
      integration <- .GenomeGMM_integrate_red_then_observed(
        natural_x = natural_x,
        fitted_counts = frequency_natural_x,
        alpha = alpha,
        miu = miu,
        sigma = sigma,
        total_count_scale = total_samples_filtered,
        full_data = full_data,
        fit_range_upper = fit_range_upper,
        biological_components = biological_components,
        bin_width = bin_width
      )
      total_kmers <- integration$total_kmers
      if (is.finite(integration$switch_x)) {
        cat("Red/grey total-kmer switch coverage:",
            round(integration$switch_x, 4), "\n")
        cat("Continuous agreement bins:",
            integration$required_agreement_bins, "\n")
      } else {
        cat("No stable red/grey agreement before the fitting boundary;",
            "using the fitted curve through that boundary.\n")
      }
      cat("Total k-mers from fitted red segment:",
          format(integration$red_total_kmers, scientific = TRUE), "\n")
      cat("Total k-mers from observed grey segment:",
          format(integration$grey_total_kmers, scientific = TRUE), "\n")
      
      # === Calculate Individual Component Sizes ===
      component_sizes_bp <- numeric(k)
      for (j in 1:k) {
        component_sizes_bp[j] <-
          ploidy * integration$component_masses[j] / kmer_depth
      }
      g_hap_bp <- total_kmers / kmer_depth
      genome_size_bp <- ploidy * g_hap_bp
      
      return(list(
        total_kmers = total_kmers,
        total_kmer_red_mass = integration$red_total_kmers,
        total_kmer_grey_mass = integration$grey_total_kmers,
        total_kmer_switch_x = integration$switch_x,
        total_kmer_agreement_tolerance = integration$agreement_tolerance,
        total_kmer_agreement_bins = integration$required_agreement_bins,
        kmer_depth = kmer_depth,
        homozygous_peak_component = homozygous_peak_component,
        g_hap_bp = g_hap_bp,
        genome_size_bp = genome_size_bp,
        component_sizes_bp = component_sizes_bp # Return individual sizes
      ))
    }
    
    #' Draw plot with smooth curve overlay on original histogram
    draw_plot <- function(x_original, y_original, alpha, miu, sigma, genome_size_bp, component_sizes_bp,
                          main_peak_x_symmetric, main_peak_y_symmetric,
                          x_fit_processed, y_fit_processed, full_data = NULL, right_value = NULL,
                          fit_extra_component = FALSE, internal_n = NULL, n = NULL,
                          species_name = "", first_peak_x = NULL,
                          monoploid_divisor = NULL,
                          prefix_completion = FALSE,
                          prefix_noise_start_x = NULL) {
      
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
      if (isTRUE(prefix_completion) && !is.null(first_peak_x) &&
          is.finite(first_peak_x)) {
        # Apply the same 1.6x headroom rule to the recovered first harmonic.
        # Use its observed height rather than the low-frequency error maximum,
        # so the y-axis reflects the genomic peak without clipping it.
        recovered_height <- stats::approx(
          plot_x, plot_y, xout = first_peak_x, rule = 2
        )$y
        if (is.finite(recovered_height) && recovered_height > 0) {
          y_max <- max(y_max, 1.6 * recovered_height)
        }
      }
      
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
      total_draw <- x_vals_fit <= x_limit
      lines(x_vals_fit[total_draw],
            total_counts[total_draw],
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
          size_divisor <- if (!is.null(monoploid_divisor) &&
              length(monoploid_divisor) == 1L &&
              is.finite(monoploid_divisor) && monoploid_divisor > 0) {
            monoploid_divisor
          } else n
          haploid_size_bp <- genome_size_bp / size_divisor
          haploid_size_mb <- haploid_size_bp / 1e6
          
          genome_text <- paste(
            "Holoploid GSE: ", sprintf("%.1f", genome_size_bp / 1e6), " Mb\n",
            "Monoploid GSE: ", sprintf("%.1f", haploid_size_mb), " Mb",
            sep = ""
          )
        } else {
          genome_text <- paste(
            "Holoploid GSE: ", sprintf("%.1f", genome_size_bp / 1e6), " Mb",
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
                               save_plot = TRUE,
                               .first_peak_override = NULL,
                               .start_index_override = NULL,
                               .prefix_completion_slot = NULL,
                               .prefix_completion_weight = 0.35,
                               .visible_peak_override = NULL) {
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
      cached_full_data <- cached_full_data[
        is.finite(as.numeric(cached_full_data[, 1])) &
          is.finite(as.numeric(cached_full_data[, 2])) &
          as.numeric(cached_full_data[, 1]) > 0,
        , drop = FALSE
      ]
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
          .first_peak_override = .first_peak_override,
          .start_index_override = .start_index_override,
          .prefix_completion_slot = .prefix_completion_slot,
          .prefix_completion_weight = .prefix_completion_weight,
          .visible_peak_override = .visible_peak_override,
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
        if (!any(valid)) {
          reasons <- sort(table(table$status), decreasing = TRUE)
          details <- paste(head(paste0(names(reasons), " [", as.integer(reasons),
                                      " candidates]"), 5L), collapse = "\n")
          stop("All automatic candidate fits failed.\n", details, call. = FALSE)
        }
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
        .first_peak_override = .first_peak_override,
        .start_index_override = .start_index_override,
        .prefix_completion_slot = .prefix_completion_slot,
        .prefix_completion_weight = .prefix_completion_weight,
        .visible_peak_override = .visible_peak_override,
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
        "Depth: %.2fx\nHoloploid GSE: %.1f Mb\nMonoploid GSE: %.1f Mb\nR2: %.4f",
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
        cat("  holoploid genome size:", round(genome$genome_size_bp / 1e6, 2), "Mb\n")
        cat("  monoploid genome size:", round(genome$genome_size_bp / n / 1e6, 2), "Mb\n")
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
  }, envir = legacy_environment)
  attr(legacy_environment, "legacy_file") <- "<embedded in genome_gmm.R>"
  legacy_environment
})

# Keep every requested biological component numerically present in automatic
# fits.  This is intentionally a very small floor: it prevents a harmonic from
# collapsing to effectively zero (as in virginalis) without forcing weak peaks
# to carry the same mass as clear peaks.  If the original peak-count rule has
# instantiated n+1, retain a smaller non-zero floor for that tail component so
# its fitted curve cannot silently disappear.
stabilize_biological_alpha <- function(alpha, required_components,
                                       minimum_share = 0.005,
                                       extra_minimum_share = 0.0025) {
  alpha <- as.numeric(alpha)
  total <- sum(alpha)
  required_components <- min(as.integer(required_components), length(alpha))
  if (!is.finite(total) || total <= 0 || required_components < 1L) return(alpha)
  proportions <- pmax(alpha / total, 0)
  proportions[seq_len(required_components)] <- pmax(
    proportions[seq_len(required_components)], minimum_share
  )
  if (length(proportions) > required_components) {
    extra_indices <- seq.int(required_components + 1L, length(proportions))
    proportions[extra_indices] <- pmax(
      proportions[extra_indices], extra_minimum_share
    )
  }
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

# Search the low-coverage prefix independently of the ordinary gradient start.
# Error-tail gradients can move that start past a real shallow dosage peak (or
# past the visible right half of one), after which no downstream fitter can
# recover the correct harmonic spacing.  This helper proposes an anchor only;
# the automatic selector still has to verify the adjacent integer-slot models
# on the complete observed range.
.GenomeGMM_find_post_noise_anchor <- function(
    x, y, ordinary_peak_x, ordinary_start_x, noise_floor = 5) {
  empty <- list(
    supported = FALSE, x = NA_real_, kind = "none", score = NA_real_,
    start_x = NA_real_, candidate_slots = integer(0),
    terminal_peak_x = NA_real_
  )
  x <- as.numeric(x)
  y <- pmax(as.numeric(y), 0)
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  ordered <- order(x)
  x <- x[ordered]
  y <- y[ordered]
  ordinary_peak_x <- sort(unique(as.numeric(ordinary_peak_x[
    is.finite(ordinary_peak_x)
  ])))
  ordinary_peak_x <- ordinary_peak_x[ordinary_peak_x >= 20]
  if (!length(ordinary_peak_x)) return(empty)

  terminal_peak_x <- ordinary_peak_x[1L]
  early_cap <- min(45, 0.55 * terminal_peak_x)
  early <- which(x >= noise_floor & x <= early_cap)
  if (length(early) < 10L) return(empty)
  early_x <- x[early]
  early_y <- y[early]
  smooth <- as.numeric(stats::filter(early_y, rep(1 / 3, 3), sides = 2))
  smooth[!is.finite(smooth)] <- early_y[!is.finite(smooth)]
  bin_steps <- diff(early_x)
  bin_steps <- bin_steps[is.finite(bin_steps) & bin_steps > 0]
  bin_width <- if (length(bin_steps)) stats::median(bin_steps) else 1

  candidate_x <- NA_real_
  candidate_start_x <- NA_real_
  candidate_kind <- "none"
  candidate_score <- -Inf

  # Grade A: an observed shallow maximum.  Both sides must persist for at
  # least two bins so a single noisy count cannot define the first harmonic.
  # Use the unsmoothed counts to preserve the centre of a very small but
  # sustained turn; the three-bin curve is used later for shoulder geometry.
  local_maxima <- which(diff(sign(diff(early_y))) < 0) + 1L
  local_maxima <- local_maxima[
    local_maxima >= 3L & local_maxima <= length(smooth) - 3L &
      early_x[local_maxima] >= noise_floor + bin_width &
      early_x[local_maxima] <= min(20, 0.35 * terminal_peak_x)
  ]
  if (length(local_maxima)) {
    local_rows <- lapply(local_maxima, function(index) {
      left <- max(1L, index - 3L):(index - 1L)
      right <- (index + 1L):min(length(smooth), index + 4L)
      rise <- early_y[index] - min(early_y[left], na.rm = TRUE)
      fall <- early_y[index] - min(early_y[right], na.rm = TRUE)
      rising_bins <- sum(diff(early_y[c(left, index)]) > 0, na.rm = TRUE)
      falling_bins <- sum(diff(early_y[c(index, right)]) < 0, na.rm = TRUE)
      score <- min(rise, fall) / max(early_y[index], 1)
      data.frame(
        index = index, score = score,
        valid = rising_bins >= 2L && falling_bins >= 2L && score >= 0.008
      )
    })
    local_table <- do.call(rbind, local_rows)
    local_table <- local_table[local_table$valid, , drop = FALSE]
    if (nrow(local_table)) {
      best <- local_table[which.max(local_table$score), , drop = FALSE]
      candidate_x <- early_x[best$index]
      left_indices <- max(1L, best$index - 3L):(best$index - 1L)
      candidate_start_x <- early_x[left_indices[
        which.min(early_y[left_indices])
      ]]
      candidate_kind <- "local_peak"
      candidate_score <- best$score
    }
  }

  # Grade B: only the right/left shoulder is visible.  Look for a sustained
  # deceleration of the initial error-tail decline followed by a renewed,
  # sustained decline.  Such a shoulder is proposed only when the historical
  # gradient start has skipped well beyond it; ordinary early-start samples do
  # not enter this weaker branch.
  if (!is.finite(candidate_x) && is.finite(ordinary_start_x) &&
      ordinary_start_x >= 2.5 * noise_floor) {
    slopes <- diff(smooth)
    shoulder_indices <- seq.int(4L, length(smooth) - 5L)
    shoulder_indices <- shoulder_indices[
      early_x[shoulder_indices] >= noise_floor + 2 * bin_width &
        early_x[shoulder_indices] <= min(20, 0.35 * terminal_peak_x)
    ]
    if (length(shoulder_indices)) {
      shoulder_rows <- lapply(shoulder_indices, function(index) {
        pre <- mean(slopes[(index - 3L):(index - 1L)], na.rm = TRUE)
        centre <- mean(slopes[index:min(length(slopes), index + 1L)],
                       na.rm = TRUE)
        post <- mean(slopes[(index + 2L):min(length(slopes), index + 5L)],
                     na.rm = TRUE)
        curvature_score <- ((centre - pre) + (centre - post)) /
          max(smooth[index], 1)
        valid_turn <- is.finite(curvature_score) && pre < 0 && post < 0 &&
          centre > pre && centre > post &&
          abs(pre) >= 1.20 * max(abs(centre), 1) &&
          abs(post) >= 1.15 * max(abs(centre), 1) &&
          curvature_score >= 0.035
        data.frame(index = index, score = curvature_score,
                   valid = valid_turn)
      })
      shoulder_table <- do.call(rbind, shoulder_rows)
      shoulder_table <- shoulder_table[
        shoulder_table$valid, , drop = FALSE
      ]
      if (nrow(shoulder_table)) {
        best <- shoulder_table[which.max(shoulder_table$score), , drop = FALSE]
        # The derivative windows are centred on the interval following index;
        # use the next bin as the physical shoulder centre.
        candidate_x <- early_x[min(length(early_x), best$index + 1L)]
        candidate_start_x <- early_x[best$index]
        candidate_kind <- "shoulder"
        candidate_score <- best$score
      }
    }
  }
  if (!is.finite(candidate_x) || candidate_x <= 0) return(empty)

  # A proposed early structure becomes a fitting anchor only when a later
  # ordinary maximum supports an adjacent integer harmonic interpretation.
  later <- ordinary_peak_x[ordinary_peak_x >= 2.5 * candidate_x]
  if (!length(later)) return(empty)
  slot_records <- lapply(later, function(later_x) {
    ratio <- later_x / candidate_x
    slots <- unique(as.integer(c(floor(ratio), ceiling(ratio))))
    # Slots 2--3 are already visible to the ordinary detector and should not
    # be reinterpreted from a low-coverage prefix. Prefix recovery is needed
    # only when that start has hidden a longer harmonic series.
    slots <- slots[
      slots >= 4L & slots <= 30L &
        abs(ratio - slots) / pmax(slots, 1L) <= 0.08
    ]
    list(x = later_x, ratio = ratio, slots = slots)
  })
  supported <- slot_records[vapply(
    slot_records, function(record) length(record$slots) > 0L, logical(1)
  )]
  if (!length(supported)) return(empty)
  selected_record <- supported[[1L]]

  list(
    supported = TRUE, x = as.numeric(candidate_x), kind = candidate_kind,
    score = as.numeric(candidate_score),
    start_x = as.numeric(candidate_start_x),
    candidate_slots = sort(unique(as.integer(selected_record$slots))),
    terminal_peak_x = as.numeric(selected_record$x)
  )
}

# Detect only explicit, well-supported genomic maxima across one common data
# range.  This scan is independent of a trial n: otherwise an n=2 fitting
# window can hide a clear sixth-harmonic peak before pattern selection begins.
detect_global_pattern_peaks <- function(full_data, analysis_max_x = NULL,
                                        minimum_height_fraction = 0.01,
                                        minimum_prominence_fraction = 0.002,
                                        minimum_local_prominence = 0.04,
                                        allow_post_noise_anchor = TRUE) {
  empty <- list(
    start_index = NA_integer_, peaks = data.frame(), first_peak_x = NA_real_,
    terminal_slot = 1L, minimum_pattern = 2L, explicit_peak_count = 0L,
    trusted_peak_count = 0L, supported_fit_slot = 2L,
    strong_tail_fit_slot = 2L, strong_tail_bridge_used = FALSE,
    weak_harmonic_fit_slot = 2L,
    weak_harmonic_chain_supported = FALSE,
    shallow_anchor_supported = FALSE,
    shallow_anchor_x = NA_real_,
    shallow_anchor_start_x = NA_real_,
    shallow_anchor_candidate_slots = integer(0),
    shallow_anchor_kind = "none",
    shallow_anchor_score = NA_real_,
    shallow_anchor_terminal_x = NA_real_,
    numerical_terminal_slot = 1L,
    harmonic_window_support = data.frame()
  )
  if (is.null(full_data) || ncol(full_data) < 2L || nrow(full_data) < 8L) {
    return(empty)
  }
  x_all <- as.numeric(full_data[, 1])
  y_all <- pmax(as.numeric(full_data[, 2]), 0)
  valid <- is.finite(x_all) & is.finite(y_all)
  x_all <- x_all[valid]
  y_all <- y_all[valid]
  ordered <- order(x_all)
  x_all <- x_all[ordered]
  y_all <- y_all[ordered]
  if (length(y_all) < 8L) return(empty)

  gradients <- diff(y_all) / pmax(y_all[-length(y_all)], 1)
  start_index <- 1L
  if (length(gradients) >= 3L) {
    for (i in seq_len(length(gradients) - 2L)) {
      if (mean(gradients[i:(i + 2L)], na.rm = TRUE) > 0.005) {
        start_index <- i
        break
      }
    }
  }
  valley_window <- start_index:min(length(y_all), start_index + 3L)
  start_index <- valley_window[which.min(y_all[valley_window])]
  scan_indices <- start_index:length(x_all)
  if (!is.null(analysis_max_x) && is.finite(analysis_max_x)) {
    scan_indices <- scan_indices[x_all[scan_indices] <= analysis_max_x]
  }
  if (length(scan_indices) < 8L) return(empty)
  x <- x_all[scan_indices]
  y <- y_all[scan_indices]
  smooth <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  smooth[!is.finite(smooth)] <- y[!is.finite(smooth)]
  peak_indices <- which(diff(sign(diff(smooth))) < 0) + 1L
  peak_indices <- peak_indices[
    peak_indices > 2L & peak_indices < length(smooth) - 2L
  ]
  if (!length(peak_indices)) return(empty)

  bin_steps <- diff(sort(unique(x)))
  bin_steps <- bin_steps[is.finite(bin_steps) & bin_steps > 0]
  bin_width <- if (length(bin_steps)) stats::median(bin_steps) else 1
  peak_scale <- max(smooth[peak_indices], na.rm = TRUE)
  peak_records <- lapply(seq_along(peak_indices), function(i) {
    index <- peak_indices[i]
    left_limit <- if (i == 1L) 1L else peak_indices[i - 1L]
    right_limit <- if (i == length(peak_indices)) {
      length(smooth)
    } else peak_indices[i + 1L]
    left_valley <- min(smooth[left_limit:index], na.rm = TRUE)
    right_valley <- min(smooth[index:right_limit], na.rm = TRUE)
    saddle <- max(left_valley, right_valley)
    prominence <- max(smooth[index] - saddle, 0)
    half_level <- saddle + 0.5 * prominence
    left_half <- index
    while (left_half > left_limit && smooth[left_half] >= half_level) {
      left_half <- left_half - 1L
    }
    right_half <- index
    while (right_half < right_limit && smooth[right_half] >= half_level) {
      right_half <- right_half + 1L
    }
    width <- max(x[right_half] - x[left_half], bin_width)
    raw_window <- max(1L, index - 2L):min(length(y), index + 2L)
    raw_support <- stats::median(y[raw_window], na.rm = TRUE) /
      max(smooth[index], 1)
    data.frame(
      index = index,
      x = estimate_continuous_peak(x, y, x[index]),
      height = smooth[index],
      prominence = prominence,
      prominence_fraction = prominence / max(peak_scale, 1),
      local_prominence = prominence / max(smooth[index], 1),
      width = width,
      raw_support = raw_support,
      stringsAsFactors = FALSE
    )
  })
  peaks <- do.call(rbind, peak_records)
  ordinary_clear <- peaks$height >= peak_scale * minimum_height_fraction &
    peaks$prominence_fraction >= minimum_prominence_fraction &
    peaks$local_prominence >= minimum_local_prominence &
    peaks$width >= 2 * bin_width &
    peaks$raw_support >= 0.45 & peaks$raw_support <= 1.80
  clear <- ordinary_clear
  shallow_anchor_supported <- FALSE
  shallow_anchor_x <- NA_real_
  shallow_anchor_start_x <- NA_real_
  shallow_anchor_candidate_slots <- integer(0)
  shallow_anchor_kind <- "none"
  shallow_anchor_score <- NA_real_
  shallow_anchor_terminal_x <- NA_real_
  peaks$shallow_anchor <- FALSE

  # A shallow first dosage peak can sit directly after the sequencing-error
  # wall.  Do not globally lower the 4% prominence threshold: accept such a
  # peak only when it has a sustained rise and fall and a later ordinary clear
  # maximum lies at an approximately integral multiple of it.  Candidate
  # terminal slots on both sides of a non-integer ratio are retained for the
  # common-domain curve comparison below.
  ordinary_indices <- which(ordinary_clear)
  post_noise_anchor <- if (isTRUE(allow_post_noise_anchor)) {
    .GenomeGMM_find_post_noise_anchor(
      x = x_all, y = y_all,
      ordinary_peak_x = peaks$x[ordinary_indices],
      ordinary_start_x = x_all[start_index]
    )
  } else {
    list(supported = FALSE)
  }
  if (isTRUE(post_noise_anchor$supported)) {
    existing_index <- which.min(abs(peaks$x - post_noise_anchor$x))
    matches_existing <- length(existing_index) == 1L &&
      is.finite(peaks$x[existing_index]) &&
      abs(peaks$x[existing_index] - post_noise_anchor$x) <= 1.5 * bin_width
    if (isTRUE(matches_existing)) {
      clear[existing_index] <- TRUE
      peaks$shallow_anchor[existing_index] <- TRUE
    } else {
      anchor_height <- stats::approx(
        x_all, y_all, xout = post_noise_anchor$x, rule = 2
      )$y
      anchor_row <- data.frame(
        index = NA_integer_, x = post_noise_anchor$x,
        height = anchor_height,
        prominence = post_noise_anchor$score * anchor_height,
        prominence_fraction = post_noise_anchor$score,
        local_prominence = post_noise_anchor$score,
        width = 3 * bin_width, raw_support = 1,
        shallow_anchor = TRUE, stringsAsFactors = FALSE
      )
      peaks <- rbind(peaks, anchor_row)
      ordinary_clear <- c(ordinary_clear, FALSE)
      clear <- c(clear, TRUE)
    }
    shallow_anchor_supported <- TRUE
    shallow_anchor_x <- post_noise_anchor$x
    shallow_anchor_start_x <- post_noise_anchor$start_x
    shallow_anchor_candidate_slots <- post_noise_anchor$candidate_slots
    shallow_anchor_kind <- post_noise_anchor$kind
    shallow_anchor_score <- post_noise_anchor$score
    shallow_anchor_terminal_x <- post_noise_anchor$terminal_peak_x
  }
  shallow_indices <- which(
    !ordinary_clear &
      peaks$x >= max(5, min(x, na.rm = TRUE)) &
      peaks$height >= peak_scale * minimum_height_fraction &
      peaks$prominence_fraction >= minimum_prominence_fraction &
      peaks$local_prominence >= 0.018 &
      peaks$width >= 2 * bin_width &
      peaks$raw_support >= 0.45 & peaks$raw_support <= 1.80
  )
  if (!shallow_anchor_supported &&
      length(ordinary_indices) && length(shallow_indices)) {
    for (candidate_index in shallow_indices[order(peaks$x[shallow_indices])]) {
      local_index <- peaks$index[candidate_index]
      left_differences <- diff(smooth[
        max(1L, local_index - 3L):local_index
      ])
      right_differences <- diff(smooth[
        local_index:min(length(smooth), local_index + 3L)
      ])
      sustained_turn <- sum(left_differences > 0, na.rm = TRUE) >= 2L &&
        sum(right_differences < 0, na.rm = TRUE) >= 2L
      if (!isTRUE(sustained_turn)) next

      later_indices <- ordinary_indices[
        peaks$x[ordinary_indices] >= 2.5 * peaks$x[candidate_index] &
          peaks$height[ordinary_indices] >=
            0.05 * peaks$height[candidate_index]
      ]
      if (!length(later_indices)) next
      slot_support <- lapply(later_indices, function(later_index) {
        ratio <- peaks$x[later_index] / peaks$x[candidate_index]
        slots <- unique(c(floor(ratio), ceiling(ratio)))
        slots <- as.integer(slots[
          slots >= 2L & slots <= 30L &
            abs(ratio - slots) / pmax(slots, 1L) <= 0.08
        ])
        list(ratio = ratio, slots = slots)
      })
      supported_slots <- sort(unique(as.integer(unlist(lapply(
        slot_support, `[[`, "slots"
      )))))
      if (!length(supported_slots)) next

      clear[candidate_index] <- TRUE
      peaks$shallow_anchor[candidate_index] <- TRUE
      shallow_anchor_supported <- TRUE
      shallow_anchor_x <- peaks$x[candidate_index]
      shallow_anchor_start_x <- x_all[start_index]
      shallow_anchor_candidate_slots <- supported_slots
      shallow_anchor_kind <- "local_peak"
      shallow_anchor_score <- peaks$local_prominence[candidate_index]
      shallow_anchor_terminal_x <- peaks$x[later_indices[1L]]
      break
    }
  }
  peaks$clear <- clear
  clear_peaks <- peaks[clear, , drop = FALSE]
  if (!nrow(clear_peaks)) return(empty)
  clear_peaks <- clear_peaks[order(clear_peaks$x), , drop = FALSE]

  # The earliest explicit genomic maximum is the observable 1x anchor. Peak
  # widths contribute only positional uncertainty; the allowed ratio residual
  # remains capped so a broad non-harmonic tail cannot be rounded freely.
  first_peak_x <- clear_peaks$x[1]
  first_uncertainty <- max(bin_width, 0.10 * clear_peaks$width[1])
  clear_peaks$ratio <- clear_peaks$x / first_peak_x
  clear_peaks$slot <- pmax(1L, as.integer(round(clear_peaks$ratio)))
  location_uncertainty <- pmax(bin_width, 0.10 * clear_peaks$width)
  ratio_uncertainty <- sqrt(
    location_uncertainty^2 +
      (clear_peaks$slot * first_uncertainty)^2
  ) / max(first_peak_x, bin_width)
  clear_peaks$ratio_tolerance <- pmin(
    0.25, pmax(0.10, 2 * ratio_uncertainty)
  )
  clear_peaks$ratio_residual <- abs(clear_peaks$ratio - clear_peaks$slot)
  ratio_match <- clear_peaks$ratio_residual <= clear_peaks$ratio_tolerance
  clear_peaks$trusted_harmonic <- FALSE
  highest_supported_slot <- 0L
  for (i in seq_len(nrow(clear_peaks))) {
    if (!ratio_match[i]) next
    slot_gap <- max(1L, clear_peaks$slot[i] - highest_supported_slot)
    # An isolated far-right maximum needs progressively stronger prominence to
    # bridge many missing slots. This retains a clear distant harmonic while
    # rejecting a weak broad fluctuation at the same relative position.
    gap_support <- slot_gap <= 1L ||
      (clear_peaks$prominence_fraction[i] >= 0.008 * slot_gap &&
       clear_peaks$local_prominence[i] >=
         0.04 + 0.02 * (slot_gap - 1L))
    if (gap_support) {
      clear_peaks$trusted_harmonic[i] <- TRUE
      highest_supported_slot <- max(
        highest_supported_slot, clear_peaks$slot[i]
      )
    }
  }
  trusted_slots <- sort(unique(clear_peaks$slot[
    clear_peaks$trusted_harmonic
  ]))
  terminal_slot <- if (length(trusted_slots)) max(trusted_slots) else 1L
  # A far integer-matched maximum may be too isolated to prove the biological
  # minimum pattern, yet it is still valuable evidence for where the numerical
  # harmonic series ends.  Keep this separate from terminal_slot so an isolated
  # 10x maximum can guide curve extent without claiming minimum ploidy 10x.
  numerical_slots <- sort(unique(clear_peaks$slot[ratio_match]))
  numerical_terminal_slot <- if (length(numerical_slots)) {
    max(numerical_slots)
  } else terminal_slot

  # A clear but non-integer terminal peak following a stable 1,2,3 prefix can
  # prove one additional unresolved slot, but can advance the pattern by only
  # one.  This retains the conservative floridensis handling.
  if (all(1:3 %in% trusted_slots)) {
    late <- clear_peaks$ratio[
      !clear_peaks$trusted_harmonic &
        clear_peaks$ratio > terminal_slot + 0.5
    ]
    if (length(late)) terminal_slot <- terminal_slot + 1L
  }
  minimum_pattern <- max(2L, as.integer(terminal_slot))

  # Explicit maxima determine the minimum observable pattern.  Numerical
  # fitting may nevertheless need to continue through later unresolved
  # harmonics when broad overlapping peaks form a plateau and only decay after
  # several integer slots.  Scan whole harmonic windows (robust median plus
  # maximum support), not isolated local maxima.  The first sharp sustained
  # drop ends the primary fitting series; one optional nuisance component can
  # still model the following decay wing.
  # Scan a data-supported number of harmonic windows rather than imposing the
  # former ten-ploid ceiling.  Thirty is only an internal safety guard against
  # pathological ultra-long count-one tails; the usable extent below is still
  # determined from the observed signal and normally stops much earlier.
  maximum_scan_slot <- min(
    30L, as.integer(floor(max(x) / max(first_peak_x, bin_width) - 0.5))
  )
  harmonic_window_support <- data.frame()
  if (maximum_scan_slot >= 1L) {
    harmonic_window_support <- do.call(rbind, lapply(
      seq_len(maximum_scan_slot), function(slot) {
        window <- which(
          x >= (slot - 0.5) * first_peak_x &
            x < (slot + 0.5) * first_peak_x
        )
        if (!length(window)) {
          return(data.frame(
            slot = slot, median_fraction = NA_real_,
            maximum_fraction = NA_real_
          ))
        }
        data.frame(
          slot = slot,
          median_fraction = stats::median(smooth[window], na.rm = TRUE) /
            max(peak_scale, 1),
          maximum_fraction = max(smooth[window], na.rm = TRUE) /
            max(peak_scale, 1)
        )
      }
    ))
  }
  supported_fit_slot <- minimum_pattern
  if (nrow(harmonic_window_support) > minimum_pattern) {
    for (slot in seq.int(minimum_pattern + 1L,
                         nrow(harmonic_window_support))) {
      previous <- harmonic_window_support$median_fraction[slot - 1L]
      current <- harmonic_window_support$median_fraction[slot]
      current_max <- harmonic_window_support$maximum_fraction[slot]
      continues <- is.finite(previous) && is.finite(current) &&
        is.finite(current_max) && current >= 0.06 && current_max >= 0.10 &&
        current / max(previous, 1e-8) >= 0.70
      if (!continues) break
      supported_fit_slot <- slot
    }
  }

  # Keep a second, deliberately non-default estimate for spectra whose broad
  # overlapping tail contains several unresolved harmonics.  The ordinary
  # estimate above stops at the first abrupt relative drop.  Here one such
  # drop may be bridged only when the absolute signal is still substantial;
  # every following slot must retain that absolute support.  This value is
  # used only to propose a guarded alternative fit, never to change the
  # inferred minimum ploidy or the normal fit by itself.
  strong_tail_fit_slot <- minimum_pattern
  strong_tail_bridge_used <- FALSE
  consecutive_unsupported_slots <- 0L
  if (nrow(harmonic_window_support) > minimum_pattern) {
    for (slot in seq.int(minimum_pattern + 1L,
                         nrow(harmonic_window_support))) {
      previous <- harmonic_window_support$median_fraction[slot - 1L]
      current <- harmonic_window_support$median_fraction[slot]
      current_max <- harmonic_window_support$maximum_fraction[slot]
      relative_support <- is.finite(previous) && is.finite(current) &&
        is.finite(current_max) && current >= 0.06 && current_max >= 0.10 &&
        current / max(previous, 1e-8) >= 0.70
      absolute_support <- is.finite(current) && is.finite(current_max) &&
        current >= 0.12 && current_max >= 0.14
      # Broad high-ploidy spectra can have a deep early decline followed by a
      # stable overlapping plateau.  Its individual maxima are weak, but a
      # whole harmonic window retaining roughly 7% of peak-scale signal is not
      # a count-one tail.  Require sustained window support and tolerate only
      # one isolated unsupported slot before stopping.
      plateau_support <- is.finite(current) && is.finite(current_max) &&
        current >= 0.07 && current_max >= 0.085
      if (!strong_tail_bridge_used && relative_support) {
        strong_tail_fit_slot <- slot
        consecutive_unsupported_slots <- 0L
      } else if (!strong_tail_bridge_used && absolute_support) {
        strong_tail_fit_slot <- slot
        strong_tail_bridge_used <- TRUE
        consecutive_unsupported_slots <- 0L
      } else if (strong_tail_bridge_used &&
          (absolute_support || plateau_support)) {
        strong_tail_fit_slot <- slot
        consecutive_unsupported_slots <- 0L
      } else {
        consecutive_unsupported_slots <- consecutive_unsupported_slots + 1L
        if (consecutive_unsupported_slots >= 2L) break
      }
    }
  }

  # A high-ploidy spectrum can retain a long, smoothly decaying harmonic
  # chain after its individually resolvable peaks have merged.  The stronger
  # plateau rule above deliberately stops once the median drops below 7% of
  # peak scale; that is too early for cases where several later dosage peaks
  # are weak but remain regularly spaced.  Extend the numerical candidate
  # range only when at least two consecutive windows beyond the strong-tail
  # endpoint retain both median and maximum support, and every step decays
  # gradually rather than collapsing abruptly.  This is fitting evidence
  # only: it cannot increase the minimum observable ploidy by itself.
  weak_harmonic_fit_slot <- as.integer(strong_tail_fit_slot)
  weak_harmonic_chain_supported <- FALSE
  if (nrow(harmonic_window_support) >= strong_tail_fit_slot + 2L) {
    proposed_weak_slot <- as.integer(strong_tail_fit_slot)
    for (slot in seq.int(strong_tail_fit_slot + 1L,
                         nrow(harmonic_window_support))) {
      current_row <- match(slot, harmonic_window_support$slot)
      previous_row <- match(slot - 1L, harmonic_window_support$slot)
      if (is.na(current_row) || is.na(previous_row)) break
      current_median <- harmonic_window_support$median_fraction[current_row]
      current_maximum <- harmonic_window_support$maximum_fraction[current_row]
      previous_median <-
        harmonic_window_support$median_fraction[previous_row]
      gradual_decay <- all(is.finite(c(
        current_median, current_maximum, previous_median
      ))) && current_median >= 0.04 && current_maximum >= 0.05 &&
        current_median / max(previous_median, 1e-8) >= 0.55
      if (!isTRUE(gradual_decay)) break
      proposed_weak_slot <- as.integer(slot)
    }
    if (proposed_weak_slot >= strong_tail_fit_slot + 2L) {
      weak_harmonic_fit_slot <- proposed_weak_slot
      weak_harmonic_chain_supported <- TRUE
    }
  }

  # Consecutive explicit 1x,2x,... maxima establish the harmonic spacing even
  # when the next component overlaps its neighbour strongly enough to appear
  # only as a shoulder.  In that case do not let the last explicit component
  # move or broaden to absorb the shoulder.  Promote exactly the next integer
  # slot only when its whole harmonic window retains substantial signal and is
  # followed by two clear decay windows.  The two-window decay requirement
  # distinguishes a terminal harmonic shoulder from an ordinary long tail.
  harmonic_shoulder_slot <- NA_integer_
  harmonic_shoulder_metrics <- list(
    into_ratio = NA_real_, first_drop_ratio = NA_real_,
    second_drop_ratio = NA_real_, shoulder_median_fraction = NA_real_,
    shoulder_maximum_fraction = NA_real_
  )
  consecutive_explicit_prefix <- terminal_slot >= 2L &&
    all(seq_len(terminal_slot) %in% trusted_slots)
  shoulder_slot <- as.integer(terminal_slot + 1L)
  if (isTRUE(consecutive_explicit_prefix) &&
      nrow(harmonic_window_support) >= shoulder_slot + 2L) {
    shoulder_rows <- match(
      c(terminal_slot, shoulder_slot, shoulder_slot + 1L,
        shoulder_slot + 2L),
      harmonic_window_support$slot
    )
    if (!anyNA(shoulder_rows)) {
      shoulder_values <- harmonic_window_support$median_fraction[
        shoulder_rows
      ]
      shoulder_maximum <- harmonic_window_support$maximum_fraction[
        shoulder_rows[2L]
      ]
      if (all(is.finite(c(shoulder_values, shoulder_maximum))) &&
          all(shoulder_values > 0)) {
        into_ratio <- shoulder_values[2L] / shoulder_values[1L]
        first_drop_ratio <- shoulder_values[3L] / shoulder_values[2L]
        second_drop_ratio <- shoulder_values[4L] / shoulder_values[3L]
        harmonic_shoulder_metrics <- list(
          into_ratio = into_ratio,
          first_drop_ratio = first_drop_ratio,
          second_drop_ratio = second_drop_ratio,
          shoulder_median_fraction = shoulder_values[2L],
          shoulder_maximum_fraction = shoulder_maximum
        )
        shoulder_supported <-
          shoulder_values[2L] >= 0.10 &&
          shoulder_maximum >= 0.12 &&
          into_ratio >= 0.55 &&
          first_drop_ratio <= 0.65 &&
          shoulder_values[2L] - shoulder_values[3L] >= 0.06 &&
          second_drop_ratio <= 0.75 &&
          shoulder_values[4L] <= 0.08
        if (isTRUE(shoulder_supported)) {
          harmonic_shoulder_slot <- shoulder_slot
          minimum_pattern <- max(minimum_pattern, shoulder_slot)
          supported_fit_slot <- max(supported_fit_slot, shoulder_slot)
        }
      }
    }
  }
  effective_start_index <- scan_indices[1]
  if (isTRUE(shallow_anchor_supported) &&
      is.finite(shallow_anchor_start_x)) {
    effective_start_index <- which.min(abs(x_all - shallow_anchor_start_x))
  }
  list(
    start_index = effective_start_index,
    ordinary_start_index = scan_indices[1], peaks = clear_peaks,
    all_peak_diagnostics = peaks, first_peak_x = first_peak_x,
    terminal_slot = as.integer(terminal_slot),
    minimum_pattern = minimum_pattern,
    explicit_peak_count = nrow(clear_peaks),
    trusted_peak_count = sum(clear_peaks$trusted_harmonic),
    numerical_terminal_slot = as.integer(numerical_terminal_slot),
    supported_fit_slot = as.integer(supported_fit_slot),
    strong_tail_fit_slot = as.integer(strong_tail_fit_slot),
    strong_tail_bridge_used = isTRUE(strong_tail_bridge_used),
    weak_harmonic_fit_slot = as.integer(weak_harmonic_fit_slot),
    weak_harmonic_chain_supported =
      isTRUE(weak_harmonic_chain_supported),
    shallow_anchor_supported = isTRUE(shallow_anchor_supported),
    shallow_anchor_x = shallow_anchor_x,
    shallow_anchor_start_x = shallow_anchor_start_x,
    shallow_anchor_candidate_slots = shallow_anchor_candidate_slots,
    shallow_anchor_kind = shallow_anchor_kind,
    shallow_anchor_score = shallow_anchor_score,
    shallow_anchor_terminal_x = shallow_anchor_terminal_x,
    harmonic_shoulder_slot = harmonic_shoulder_slot,
    harmonic_shoulder_metrics = harmonic_shoulder_metrics,
    harmonic_window_support = harmonic_window_support,
    scan_end_x = max(x)
  )
}

# For a spectrum with only one explicit maximum, distinguish the end of a
# broad biological plateau from its decaying right tail.  K is accepted as an
# elbow only when the signal is stable into K, then drops substantially in two
# consecutive harmonic windows.  This is evidence for a candidate K, not a
# sample- or ploidy-specific rule; a fitted-curve quality guard is applied by
# the automatic selector before K can replace a longer numerical extent.
.GenomeGMM_single_peak_elbow <- function(global_pattern_evidence) {
  empty <- list(
    slot = NA_integer_, detected = FALSE,
    pre_ratio = NA_real_, first_drop_ratio = NA_real_,
    second_drop_ratio = NA_real_
  )
  if (is.null(global_pattern_evidence) ||
      global_pattern_evidence$explicit_peak_count != 1L ||
      global_pattern_evidence$trusted_peak_count != 1L) return(empty)
  support <- global_pattern_evidence$harmonic_window_support
  if (is.null(support) || nrow(support) < 4L ||
      !all(c("slot", "median_fraction") %in% names(support))) {
    return(empty)
  }
  support <- support[
    is.finite(support$slot) & is.finite(support$median_fraction),
    , drop = FALSE
  ]
  support <- support[order(support$slot), , drop = FALSE]
  maximum_k <- max(support$slot) - 2L
  if (!is.finite(maximum_k) || maximum_k < 2L) return(empty)

  for (k in seq.int(2L, maximum_k)) {
    indices <- match(c(k - 1L, k, k + 1L, k + 2L), support$slot)
    if (anyNA(indices)) next
    values <- support$median_fraction[indices]
    if (any(!is.finite(values)) || any(values <= 0)) next
    pre_ratio <- values[2] / values[1]
    first_drop_ratio <- values[3] / values[2]
    second_drop_ratio <- values[4] / values[3]
    stable_into_k <- pre_ratio >= 0.72
    first_clear_drop <- first_drop_ratio <= 0.70 &&
      values[2] - values[3] >= 0.12
    second_clear_drop <- second_drop_ratio <= 0.65 && values[4] <= 0.16
    plateau_supported <- values[2] >= 0.20
    if (stable_into_k && first_clear_drop && second_clear_drop &&
        plateau_supported) {
      return(list(
        slot = as.integer(k), detected = TRUE,
        pre_ratio = pre_ratio,
        first_drop_ratio = first_drop_ratio,
        second_drop_ratio = second_drop_ratio
      ))
    }
  }
  empty
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
                      .preloaded_data = NULL, .analysis_max_x = NULL,
                      .global_pattern_evidence = NULL,
                      .fit_harmonic_count = NULL,
                      .initial_fit = NULL,
                      .first_peak_override = NULL,
                      .start_index_override = NULL, ...) {
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
  fit_harmonic_count <- if (!is.null(.fit_harmonic_count) &&
      length(.fit_harmonic_count) == 1L &&
      is.finite(.fit_harmonic_count)) {
    as.integer(.fit_harmonic_count)
  } else {
    as.integer(n)
  }
  fit_harmonic_count <- max(as.integer(n), fit_harmonic_count)
  
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
  full_data <- full_data[
    is.finite(as.numeric(full_data[, 1])) &
      is.finite(as.numeric(full_data[, 2])) &
      as.numeric(full_data[, 1]) > 0,
    , drop = FALSE
  ]
  
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
  if (!is.null(.start_index_override) &&
      length(.start_index_override) == 1L &&
      is.finite(.start_index_override)) {
    start_index <- max(
      1L, min(nrow(full_data), as.integer(.start_index_override))
    )
  }
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
  if (!is.null(.first_peak_override) &&
      length(.first_peak_override) == 1L &&
      is.finite(.first_peak_override) && .first_peak_override > 0) {
    first_peak_x <- as.numeric(.first_peak_override)
  }
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
  focus_range_end <- first_main_peak_x * (fit_harmonic_count + 2)
  
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

  # A local maximum counts as an explicit fitting peak only when it belongs
  # to one of the supported 1..K harmonic slots.  A weak later-harmonic bump
  # must not occupy the unresolved fourth component (the floridensis failure).
  # Keep all maxima in peak_x_final so an optional n+1 component can still use
  # them as nuisance-position evidence.
  harmonic_ratio <- continuous_peak_x / first_peak_x
  harmonic_slot <- round(harmonic_ratio)
  biological_mask <- is.finite(harmonic_ratio) &
    harmonic_slot >= 1L & harmonic_slot <= fit_harmonic_count &
    abs(harmonic_ratio - harmonic_slot) <= 0.25
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
  raw_detected_peak_x <- peak_x_final
  raw_first_peak_x <- first_peak_x
  raw_harmonic_ratio <- raw_detected_peak_x / raw_first_peak_x
  raw_harmonic_slot <- round(raw_harmonic_ratio)
  trusted_harmonic <- is.finite(raw_harmonic_ratio) &
    raw_harmonic_slot >= 1L &
    abs(raw_harmonic_ratio - raw_harmonic_slot) <= 0.25
  terminal_harmonic_slot <- if (any(trusted_harmonic)) {
    max(as.integer(raw_harmonic_slot[trusted_harmonic]))
  } else {
    1L
  }
  # A prominent late maximum can sit beyond an unresolved last biological
  # harmonic.  It is evidence for one missing terminal slot, but not evidence
  # that its own rounded index is biological.  Require a stable 1,2,3 prefix
  # and advance by only one slot; this turns a 4.7x tail into conservative
  # support for slot 4 instead of the floridensis 5x over-call.
  trusted_slots <- sort(unique(as.integer(
    raw_harmonic_slot[trusted_harmonic]
  )))
  has_stable_prefix <- all(1:3 %in% trusted_slots)
  late_nonharmonic <- raw_harmonic_ratio[
    !trusted_harmonic & is.finite(raw_harmonic_ratio)
  ]
  if (has_stable_prefix && length(late_nonharmonic) &&
      any(late_nonharmonic > terminal_harmonic_slot + 0.5)) {
    terminal_harmonic_slot <- terminal_harmonic_slot + 1L
  }
  global_pattern_peak_x <- numeric(0)
  global_pattern_ratios <- numeric(0)
  global_pattern_slots <- integer(0)
  global_pattern_peak_count <- 0L
  global_pattern_trusted_count <- 0L
  if (!is.null(.global_pattern_evidence) &&
      is.finite(.global_pattern_evidence$terminal_slot)) {
    terminal_harmonic_slot <- max(
      terminal_harmonic_slot,
      as.integer(.global_pattern_evidence$terminal_slot)
    )
    global_peaks <- .global_pattern_evidence$peaks
    if (!is.null(global_peaks) && nrow(global_peaks)) {
      global_pattern_peak_x <- as.numeric(global_peaks$x)
      global_pattern_ratios <- as.numeric(global_peaks$ratio)
      global_pattern_slots <- as.integer(
        global_peaks$slot[global_peaks$trusted_harmonic]
      )
      global_pattern_peak_count <- nrow(global_peaks)
      global_pattern_trusted_count <- sum(global_peaks$trusted_harmonic)
    }
  }
  
  cat("Final detected peak positions in focus range:\n")
  print(peak_x_final)
  cat("Number of peaks detected in focus range:", length(peak_x_final), "\n")
  cat("Matched harmonic peak positions (1..fit count):\n")
  print(observed_peak_x)
  cat("Number of biologically matched peaks:", observed_peak_count, "\n")
  if (global_pattern_peak_count > 0L) {
    cat("Global explicit peak positions:",
        paste(round(global_pattern_peak_x, 4), collapse = " "), "\n")
    cat("Global peak/first-peak ratios:",
        paste(round(global_pattern_ratios, 4), collapse = " "), "\n")
    cat("Global trusted harmonic slots:",
        paste(global_pattern_slots, collapse = " "), "\n")
  }
  cat("Rightmost trusted harmonic slot:", terminal_harmonic_slot, "\n")
  
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
  if (length(peak_x_final) < fit_harmonic_count) {
    missing_peaks <- fit_harmonic_count - length(peak_x_final)
    shoulder_peaks <- detect_shoulder_peaks(x_filtered, y_filtered, first_peak_x, first_peak_y, missing_peaks, full_data[start_index, 1])
    
    if (!is.null(shoulder_peaks) && length(shoulder_peaks) > 0) {
      candidate_peak_x <- sort(c(peak_x_final, shoulder_peaks))
      # Do not redefine the harmonic base from an artificial half-depth
      # shoulder when two or more explicit maxima already form a stable
      # 1x,2x,... series.  A single explicit hump (dros2/bird-like spectra)
      # may still use a supported left shoulder as the missing first harmonic.
      raw_trusted_slots <- sort(unique(as.integer(
        raw_harmonic_slot[trusted_harmonic]
      )))
      explicit_series_is_stable <- all(c(1L, 2L) %in% raw_trusted_slots)
      proposed_first <- candidate_peak_x[1]
      redefines_base <- is.finite(proposed_first) &&
        proposed_first < 0.75 * raw_first_peak_x
      if (explicit_series_is_stable && redefines_base) {
        shoulder_peaks <- shoulder_peaks[
          shoulder_peaks >= 0.75 * raw_first_peak_x
        ]
        candidate_peak_x <- sort(c(peak_x_final, shoulder_peaks))
      }
      peak_x_final <- candidate_peak_x
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
    
    if (length(peak_x_final) < fit_harmonic_count) {
      if (is_slow_decaying_data(
          y_filtered, x_filtered, first_peak_x,
          fit_harmonic_count, peak_x_final)) {
        peak_multiplier <- fit_harmonic_count + 2
      } else {
        peak_multiplier <- fit_harmonic_count + 1
      }
      right_value <- peak_x_final[1] * peak_multiplier  
    } else {
      nth_peak_x <- peak_x_final[fit_harmonic_count]
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
    nth_peak_x <- peak_x_final[fit_harmonic_count]
    if (length(peak_x_final) == fit_harmonic_count) {
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

  # Keep biological pattern selection separate from the numerical fitting
  # endpoint. Endpoint extension is reserved for four-or-more-component
  # overlap, where several unresolved late harmonics can merge into a long
  # wing (the sturgeon case). In lower patterns a guarded nuisance component
  # represents a structured tail without making every biological peak broader
  # (the virginalis failure).
  endpoint_smooth <- as.numeric(stats::filter(
    y_filtered, rep(1 / 5, 5), sides = 2
  ))
  endpoint_smooth[!is.finite(endpoint_smooth)] <-
    y_filtered[!is.finite(endpoint_smooth)]
  endpoint_scale <- max(endpoint_smooth, na.rm = TRUE)
  current_endpoint_index <- which.min(abs(x_filtered - right_value))
  current_endpoint_ratio <- endpoint_smooth[current_endpoint_index] /
    max(endpoint_scale, 1)
  fit_endpoint_extended <- FALSE
  original_right_value <- right_value
  if (fit_harmonic_count >= 4L &&
      max(terminal_harmonic_slot, fit_harmonic_count) >= 3L &&
      is.finite(current_endpoint_ratio) && current_endpoint_ratio > 0.035) {
    endpoint_slot_cap <- max(
      as.integer(fit_harmonic_count), terminal_harmonic_slot
    ) + 4L
    endpoint_cap <- min(max(x_filtered), raw_first_peak_x * endpoint_slot_cap)
    endpoint_indices <- which(
      x_filtered >= right_value & x_filtered <= endpoint_cap
    )
    extended_right <- endpoint_cap
    if (length(endpoint_indices) >= 5L) {
      below <- endpoint_smooth[endpoint_indices] <= 0.025 * endpoint_scale
      run_start <- which(vapply(
        seq_len(length(below) - 4L),
        function(i) all(below[i:(i + 4L)]), logical(1)
      ))
      if (length(run_start)) {
        extended_right <- x_filtered[endpoint_indices[run_start[1] + 2L]]
      }
    }
    if (is.finite(extended_right) && extended_right > right_value) {
      right_value <- extended_right
      fit_endpoint_extended <- TRUE
      cat("Extended fit endpoint from", round(original_right_value, 3),
          "to", round(right_value, 3),
          "because the original boundary retained",
          sprintf("%.1f%%", 100 * current_endpoint_ratio),
          "of peak-scale signal.\n")
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
      x_fit_processed, y_fit_processed, weight_fit_processed,
      fit_harmonic_count,
      first_main_peak_x, first_peak_x, peak_x_final,
      main_peak_x_symmetric, main_peak_y_symmetric,
      max_sigma_ratio = max_sigma_ratio,
      sigma_mode = sigma_mode,
      position_mode = position_mode,
      observed_peak_count = observed_peak_count,
      max_iter = em_max_iter,
      tolerance = em_tolerance,
      extra_component_mode = extra_component_mode,
      initial_fit = .initial_fit
    )
  } else {
    cat("Using original EM algorithm without calibration\n")
    # Execute original EM algorithm without calibration
    EM_result <- EM_algorithm(
      x_fit_processed, y_fit_processed, weight_fit_processed,
      fit_harmonic_count,
      first_main_peak_x, first_peak_x, peak_x_final,
      max_sigma_ratio = max_sigma_ratio,
      sigma_mode = sigma_mode,
      position_mode = position_mode,
      observed_peak_count = observed_peak_count,
      max_iter = em_max_iter,
      tolerance = em_tolerance,
      extra_component_mode = extra_component_mode,
      initial_fit = .initial_fit
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
    alpha <- stabilize_biological_alpha(alpha, fit_harmonic_count)
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
    biological_component_guard <- function(candidate_alpha, candidate_miu,
                                           candidate_sigma) {
      component_order <- order(candidate_miu)
      biological_count <- min(
        as.integer(fit_harmonic_count), length(component_order)
      )
      biological <- component_order[seq_len(biological_count)]
      proportions <- pmax(candidate_alpha, 0) /
        max(sum(pmax(candidate_alpha, 0)), 1e-12)
      strengths <- pmax(candidate_alpha, 0) / pmax(candidate_sigma, 1e-8)
      relative_strength <- strengths[biological] /
        max(strengths[biological], 1e-12)
      share_floor <- min(
        0.008, 0.04 / max(as.integer(fit_harmonic_count), 1L)
      )
      peak_floor <- min(
        0.020, 0.10 / max(as.integer(fit_harmonic_count), 1L)
      )
      nuisance <- if (length(component_order) > biological_count) {
        component_order[(biological_count + 1L):length(component_order)]
      } else integer(0)
      nuisance_guard <- TRUE
      if (length(nuisance)) {
        spacing <- if (biological_count > 1L) {
          stats::median(diff(sort(candidate_miu[biological])))
        } else max(candidate_miu[biological], 1)
        nuisance_guard <- sum(proportions[nuisance]) <= 0.30 &&
          max(strengths[nuisance]) /
            max(strengths[biological], 1e-12) <= 0.65 &&
          min(candidate_miu[nuisance]) >=
            max(candidate_miu[biological]) + 0.15 * spacing
      }
      all(is.finite(proportions[biological])) &&
        all(is.finite(relative_strength)) &&
        all(proportions[biological] >= share_floor) &&
        all(relative_strength >= peak_floor) && isTRUE(nuisance_guard)
    }
    overlap_case <- fit_harmonic_count >= 5L &&
      observed_peak_count < fit_harmonic_count &&
      mean(EM_result$sigma_reliability, na.rm = TRUE) < 0.35
    low_overlap_case <- fit_harmonic_count <= 3L &&
      observed_peak_count < fit_harmonic_count &&
      mean(EM_result$sigma_reliability, na.rm = TRUE) < 0.35
    total_peak_case <- fit_harmonic_count >= 2L &&
      observed_peak_count >= fit_harmonic_count &&
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
        candidate$alpha <- stabilize_biological_alpha(
          candidate$alpha, fit_harmonic_count
        )
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
      component_guard <- !is.null(candidate_score) &&
        biological_component_guard(
          candidate$alpha, candidate$miu, candidate$sigma
        )
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
      # Every refinement profile must preserve the visible left wing.  The
      # previous region-balanced winner could lower its aggregate score while
      # raising the whole dros2 left flank above the grey reference.
      left_wing_guard <- !is.null(candidate_score) &&
        candidate_score$left_overshoot_error <=
          max(0.0025,
              baseline_joint_score$left_overshoot_error * 1.03 + 0.0003) &&
        candidate_score$left_relative_error <=
          max(0.035,
              baseline_joint_score$left_relative_error * 1.03 + 0.002)
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
        clear_underfit_guard && curve_guard && left_wing_guard &&
        component_guard && position_guard && local_guard

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
      # Pareto exceptions may relax shape trade-offs, never disappearance of a
      # requested biological component.
      eligible <- isTRUE(eligible) && isTRUE(component_guard)
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
        candidate$alpha <- stabilize_biological_alpha(
          candidate$alpha, fit_harmonic_count
        )
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
      component_guard <- !is.null(candidate_score) &&
        biological_component_guard(
          candidate$alpha, candidate$miu, candidate$sigma
        )
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
      left_wing_guard <- !is.null(candidate_score) &&
        candidate_score$left_overshoot_error <=
          max(0.0025,
              baseline_joint_score$left_overshoot_error * 1.03 + 0.0003) &&
        candidate_score$left_relative_error <=
          max(0.035,
              baseline_joint_score$left_relative_error * 1.03 + 0.002)
      eligible <- is.finite(score_improvement) &&
        score_improvement >= 0.01 && main_peak_guard && height_guard &&
        peak_ratio_guard && clear_underfit_guard && curve_guard &&
        left_wing_guard && component_guard
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
    fit_harmonic_count
  }
  
  # Print relevant information
  cat("[Extra Component Info]\n")
  cat("Fitted extra components:", fit_extra_component, "\n")
  cat("Internal fitted components count:", internal_n, "\n")
  cat("Fitted curve ploidy:", n, "\n")
  cat("Tail-supported primary fitting harmonics:",
      fit_harmonic_count, "\n")
  
  if (fit_extra_component) {
    cat(sprintf("Note: Extra components fitted (%d total), components from %d onwards will be shown in light gray\n", 
                internal_n, fit_harmonic_count + 1L))
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
  
  # Use the highest primary biological harmonic as the homozygous/full-copy
  # depth. Optional tail nuisance components remain excluded.
  reporting_ploidy <- getOption("GenomeGMM.reporting_ploidy", n)
  if (!is.numeric(reporting_ploidy) || length(reporting_ploidy) != 1L ||
      !is.finite(reporting_ploidy) || reporting_ploidy < 1L) {
    reporting_ploidy <- n
  }
  reporting_ploidy <- as.integer(reporting_ploidy)
  genome_size_result <- calculate_genome_size(
    x_fit_processed, y_fit_processed, alpha, miu, sigma, 
    full_data, start_index, right_value, 
    biological_component_count = fit_harmonic_count,
    ploidy = reporting_ploidy
  )
  
  total_kmers <- genome_size_result$total_kmers
  kmer_depth <- genome_size_result$kmer_depth
  g_hap_bp <- genome_size_result$g_hap_bp
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
  cat("Homozygous peak depth:", kmer_depth, "\n")
  cat("g-hap size:", g_hap_bp, "bp\n")
  
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
    filename <- paste0(
      "em_fitted_gmm_", file_basename, "_estimated_ploidy_",
      fit_harmonic_count, "x_", timestamp, ".pdf"
    )
    
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
      primary_component_count = fit_harmonic_count,
      first_peak_x = first_peak_x,
      n = n,
      monoploid_divisor = getOption(
        "GenomeGMM.reporting_monoploid_divisor", reporting_ploidy
      ),
      species_name = species_name
    )
    dev.off()
    cat("Plot saved as:", filename, "\n")
  }
  cat("Automatic fit score (lower is better):", round(fit_score_info$score, 6), "\n")
  cat("Algorithm execution completed!\n")
  
  return(invisible(list(
    genome_size_bp = genome_size_bp,
    minimum_pattern_genome_size_bp = genome_size_bp,
    haploid_genome_size_bp = g_hap_bp,
    g_hap_bp = g_hap_bp,
    total_kmers = total_kmers,
    total_kmer_red_mass = genome_size_result$total_kmer_red_mass,
    total_kmer_grey_mass = genome_size_result$total_kmer_grey_mass,
    total_kmer_switch_x = genome_size_result$total_kmer_switch_x,
    total_kmer_agreement_tolerance =
      genome_size_result$total_kmer_agreement_tolerance,
    total_kmer_agreement_bins =
      genome_size_result$total_kmer_agreement_bins,
    total_samples_fit = sum(y_fit_processed),
    kmer_depth = kmer_depth,
    homozygous_peak_depth = kmer_depth,
    homozygous_peak_component = genome_size_result$homozygous_peak_component,
    first_peak_depth = if (length(miu)) min(miu, na.rm = TRUE) else NA_real_,
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
    fit_harmonic_count = fit_harmonic_count,
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
    observed_peak_count = observed_peak_count,
    observed_peak_x = observed_peak_x,
    detected_peak_x = raw_detected_peak_x,
    global_pattern_peak_count = global_pattern_peak_count,
    global_pattern_trusted_count = global_pattern_trusted_count,
    global_pattern_peak_x = global_pattern_peak_x,
    global_pattern_ratios = global_pattern_ratios,
    global_pattern_slots = global_pattern_slots,
    global_pattern_evidence = .global_pattern_evidence,
    terminal_harmonic_slot = terminal_harmonic_slot,
    main_peak_ratio = main_peak_ratio,
    fit_score = fit_score_info$score,
    fit_metrics = fit_score_info,
    right_value = right_value,
    original_right_value = original_right_value,
    fit_endpoint_extended = fit_endpoint_extended,
    endpoint_signal_ratio = current_endpoint_ratio,
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
                                          extra_component_mode = "legacy",
                                          initial_fit = NULL) {
  
  # Step 1: Perform normal EM fitting
  EM_result <- EM_algorithm(x, y, weights, n, first_main_peak_x, first_peak_x, peak_x_focused,
                            max_sigma_ratio = max_sigma_ratio,
                            sigma_mode = sigma_mode,
                            position_mode = position_mode,
                            observed_peak_count = observed_peak_count,
                            max_iter = max_iter, tolerance = tolerance,
                            extra_component_mode = extra_component_mode,
                            initial_fit = initial_fit)
  
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

  calibration_guard <- .GenomeGMM_check_main_peak_calibration(
    miu = miu,
    sigma = sigma,
    main_peak_component = main_peak_component,
    calibrated_mean = miu_calibrated[main_peak_component],
    calibrated_sd = sigma_calibrated[main_peak_component],
    first_peak_x = first_peak_x,
    max_sigma_ratio = max_sigma_ratio
  )
  if (!isTRUE(calibration_guard$accepted)) {
    cat(
      "Main peak calibration rejected by biological position/width guard:",
      calibration_guard$reason, "\n"
    )
    return(list(
      alpha = alpha, miu = miu, sigma = sigma,
      calibration_accepted = FALSE,
      calibration_guard = calibration_guard
    ))
  }
  miu_calibrated[main_peak_component] <- calibration_guard$projected_mean
  sigma_calibrated[main_peak_component] <- calibration_guard$projected_sd
  cat(
    "Main peak calibration", calibration_guard$reason,
    sprintf("(miu=%.4f, sigma=%.4f)\n",
            miu_calibrated[main_peak_component],
            sigma_calibrated[main_peak_component])
  )
  
  # Adjust mixing coefficients to maintain total probability of 1
  total_alpha <- sum(alpha_calibrated)
  alpha_calibrated <- alpha_calibrated / total_alpha
  
  cat("=== Main peak calibration completed ===\n")
  
  return(list(
    alpha = alpha_calibrated,
    miu = miu_calibrated,
    sigma = sigma_calibrated,
    calibration_accepted = TRUE,
    calibration_guard = calibration_guard
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
                         observed_peak_count = NULL,
                         initial_fit = NULL) {
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

  # A K+1 validation starts from the already optimised K-component model.
  # Preserve those components and initialise only the final tail component;
  # EM may then redistribute the last components locally.  Invalid or
  # incompatible seeds fall back to the standard deterministic initialisation.
  warm_started <- FALSE
  if (!is.null(initial_fit) && is.list(initial_fit) &&
      all(c("alpha", "miu", "sigma") %in% names(initial_fit))) {
    seed_alpha <- as.numeric(initial_fit$alpha)
    seed_miu <- as.numeric(initial_fit$miu)
    seed_sigma <- as.numeric(initial_fit$sigma)
    seed_count <- min(length(seed_alpha), length(seed_miu),
                      length(seed_sigma))
    valid_seed <- seed_count >= 1L && seed_count <= internal_n &&
      all(is.finite(seed_alpha[seq_len(seed_count)])) &&
      all(is.finite(seed_miu[seq_len(seed_count)])) &&
      all(is.finite(seed_sigma[seq_len(seed_count)])) &&
      all(seed_alpha[seq_len(seed_count)] >= 0) &&
      all(seed_sigma[seq_len(seed_count)] > 0)
    if (valid_seed) {
      seed_order <- order(seed_miu[seq_len(seed_count)])
      seed_alpha <- seed_alpha[seed_order]
      seed_miu <- seed_miu[seed_order]
      seed_sigma <- seed_sigma[seed_order]
      inherited_count <- min(seed_count, internal_n)
      miu[seq_len(inherited_count)] <- seed_miu[seq_len(inherited_count)]
      sigma[seq_len(inherited_count)] <-
        seed_sigma[seq_len(inherited_count)]

      if (inherited_count < internal_n) {
        for (j in seq.int(inherited_count + 1L, internal_n)) {
          miu[j] <- max(first_peak_x * j,
                        miu[j - 1L] + 0.75 * first_peak_x)
          sigma_growth <- if (inherited_count >= 2L) {
            median(tail(sigma[seq_len(inherited_count)] /
                          pmax(c(sigma[1], sigma[seq_len(inherited_count - 1L)]),
                               1e-8), 2L), na.rm = TRUE)
          } else 1
          if (!is.finite(sigma_growth)) sigma_growth <- 1
          sigma[j] <- sigma[j - 1L] * pmin(1.20, pmax(0.95, sigma_growth))
        }
        extra_mass <- min(0.08, max(0.025,
          tail(seed_alpha, 1L) * 0.30))
        inherited_alpha <- pmax(seed_alpha[seq_len(inherited_count)], 0)
        inherited_alpha <- inherited_alpha /
          max(sum(inherited_alpha), 1e-12) * (1 - extra_mass)
        alpha <- c(inherited_alpha,
                   rep(extra_mass / (internal_n - inherited_count),
                       internal_n - inherited_count))
      } else {
        alpha <- pmax(seed_alpha[seq_len(internal_n)], 0)
        alpha <- alpha / max(sum(alpha), 1e-12)
      }
      warm_started <- TRUE
      cat("Warm-started EM from the selected primary-component fit\n")
      cat("Warm-started means:", round(miu, 2), "\n")
    }
  }
  
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

#' Calculate g-hap from the highest biological harmonic
calculate_genome_size <- function(x_fit_processed, y_fit_processed, alpha, miu, sigma,
                                  full_data, start_index, right_value,
                                  biological_component_count = length(miu),
                                  ploidy = biological_component_count) {
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
  
  biological_component_count <- min(
    max(1L, as.integer(biological_component_count)), k
  )
  biological_components <- seq_len(biological_component_count)
  valid_components <- biological_components[
    is.finite(miu[biological_components]) & miu[biological_components] > 0
  ]
  if (!length(valid_components)) {
    stop("No valid biological component is available for genome-size estimation.")
  }
  homozygous_peak_component <- valid_components[
    which.max(miu[valid_components])
  ]
  kmer_depth <- as.numeric(miu[homozygous_peak_component])
  ploidy <- max(1L, as.integer(ploidy))
  
  cat("Homozygous peak position (kmer_depth):", kmer_depth, "\n")
  cat("Homozygous peak component:", homozygous_peak_component, "\n")
  cat("All component means:", round(miu, 2), "\n")
  
  # For data within fitting range, use GMM estimated distribution
  fit_range_upper <- max(x_fit_processed)
  natural_x <- seq(min(x_fit_processed), fit_range_upper, by = 1)
  density_natural_x <- approx(x_vals, total_density, xout = natural_x)$y
  
  # Convert probability density to counts
  frequency_natural_x <- density_natural_x * total_samples_filtered * bin_width
  integration <- .GenomeGMM_integrate_red_then_observed(
    natural_x = natural_x,
    fitted_counts = frequency_natural_x,
    alpha = alpha,
    miu = miu,
    sigma = sigma,
    total_count_scale = total_samples_filtered,
    full_data = full_data,
    fit_range_upper = fit_range_upper,
    biological_components = biological_components,
    bin_width = bin_width
  )
  total_kmers <- integration$total_kmers
  if (is.finite(integration$switch_x)) {
    cat("Red/grey total-kmer switch coverage:",
        round(integration$switch_x, 4), "\n")
    cat("Continuous agreement bins:",
        integration$required_agreement_bins, "\n")
  } else {
    cat("No stable red/grey agreement before the fitting boundary;",
        "using the fitted curve through that boundary.\n")
  }
  cat("Total k-mers from fitted red segment:",
      format(integration$red_total_kmers, scientific = TRUE), "\n")
  cat("Total k-mers from observed grey segment:",
      format(integration$grey_total_kmers, scientific = TRUE), "\n")
  
  # === Calculate Individual Component Sizes ===
  component_sizes_bp <- numeric(k)
  for (j in 1:k) {
    component_sizes_bp[j] <-
      ploidy * integration$component_masses[j] / kmer_depth
  }
  g_hap_bp <- total_kmers / kmer_depth
  genome_size_bp <- ploidy * g_hap_bp
  
  return(list(
    total_kmers = total_kmers,
    total_kmer_red_mass = integration$red_total_kmers,
    total_kmer_grey_mass = integration$grey_total_kmers,
    total_kmer_switch_x = integration$switch_x,
    total_kmer_agreement_tolerance = integration$agreement_tolerance,
    total_kmer_agreement_bins = integration$required_agreement_bins,
    kmer_depth = kmer_depth,
    homozygous_peak_component = homozygous_peak_component,
    g_hap_bp = g_hap_bp,
    genome_size_bp = genome_size_bp,
    component_sizes_bp = component_sizes_bp # Return individual sizes
  ))
}

#' Draw plot with smooth curve overlay on original histogram
draw_plot <- function(x_original, y_original, alpha, miu, sigma, genome_size_bp, component_sizes_bp,
                      main_peak_x_symmetric, main_peak_y_symmetric,
                      x_fit_processed, y_fit_processed, full_data = NULL, right_value = NULL,
                      fit_extra_component = FALSE, internal_n = NULL, n = NULL,
                      species_name = "", first_peak_x = NULL,
                      monoploid_divisor = NULL,
                      primary_component_count = NULL) {
  
  k <- length(alpha)
  n_components_total <- length(component_sizes_bp)
  has_residual <- n_components_total > k
  if (is.null(primary_component_count)) primary_component_count <- n
  
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
  plot_harmonic_count <- if (!is.null(n) &&
      !is.null(primary_component_count) &&
      is.finite(primary_component_count)) {
    max(as.integer(n), as.integer(primary_component_count))
  } else {
    n
  }
  if (!is.null(n) && n <= 2L && !is.null(right_value) &&
      is.finite(right_value)) {
    x_limit <- right_value
  } else if (!is.null(first_peak_x) && !is.null(plot_harmonic_count)) {
    # Keep the biological minimum pattern n in the report and legend, while
    # showing every primary harmonic that was actually needed for the curve.
    x_limit <- first_peak_x * (plot_harmonic_count + 2)
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
    if (fit_extra_component && !is.null(primary_component_count) &&
        !is.null(internal_n) &&
        internal_n > primary_component_count &&
        j > primary_component_count) {
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
      if (fit_extra_component && !is.null(primary_component_count) &&
          j > primary_component_count) {
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
    if (fit_extra_component && !is.null(primary_component_count) &&
        j > primary_component_count) {
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
      size_divisor <- if (!is.null(monoploid_divisor) &&
          length(monoploid_divisor) == 1L &&
          is.finite(monoploid_divisor) && monoploid_divisor > 0) {
        monoploid_divisor
      } else n
      haploid_size_bp <- genome_size_bp / size_divisor
      haploid_size_mb <- haploid_size_bp / 1e6
      
      genome_text <- paste(
        "Holoploid GSE: ", sprintf("%.1f", genome_size_bp / 1e6), " Mb\n",
        "Monoploid GSE: ", sprintf("%.1f", haploid_size_mb), " Mb",
        sep = ""
      )
    } else {
      genome_text <- paste(
        "Holoploid GSE: ", sprintf("%.1f", genome_size_bp / 1e6), " Mb",
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

# Derive the finite coverage range used by global peak detection. Histograms
# may contain hundreds of thousands of sparse tail rows, but those rows do not
# carry useful peak evidence and make an uncapped local-maximum scan expensive.
# This is the same data-supported rule used by the automatic search and is also
# reused by the known-n prefix-recovery gate.
.GenomeGMM_auto_analysis_limit <- function(cached_full_data, n = NULL) {
  y_values <- suppressWarnings(as.numeric(cached_full_data[, 2]))
  x_values <- suppressWarnings(as.numeric(cached_full_data[, 1]))
  row_count <- length(y_values)
  if (row_count < 2L) return(max(x_values, na.rm = TRUE))
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
  probe_end <- min(nrow(cached_full_data), start_index + 8191L)
  probe <- cached_full_data[start_index:probe_end, , drop = FALSE]
  smooth <- as.numeric(stats::filter(probe[, 2], rep(1 / 5, 5), sides = 2))
  smooth[!is.finite(smooth)] <- probe[!is.finite(smooth), 2]
  peak_indices <- which(diff(sign(diff(smooth))) < 0) + 1L
  if (length(peak_indices) == 0L) peak_indices <- which.max(smooth)
  strongest_peak <- probe[peak_indices[which.max(smooth[peak_indices])], 1]
  effective_n <- if (is.null(n)) 30L else as.integer(n)
  max(strongest_peak * (effective_n + 3L), probe[1, 1] + 100)
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
                           save_plot = TRUE,
                           .selection_only = FALSE,
                           .preloaded_data = NULL,
                           .allow_post_noise_anchor = TRUE,
                           .monoploid_range = NULL) {
  search <- match.arg(search)
  position_mode <- match.arg(position_mode)
  joint_refinement <- match.arg(joint_refinement)
  if (identical(n, "auto")) n <- NULL
  inference_requested <- is.null(n)

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
  cached_full_data <- if (is.null(.preloaded_data)) {
    read.table(data_file, header = FALSE)
  } else {
    .preloaded_data
  }
  cached_full_data <- cached_full_data[
    is.finite(as.numeric(cached_full_data[, 1])) &
      is.finite(as.numeric(cached_full_data[, 2])) &
      as.numeric(cached_full_data[, 1]) > 0,
    , drop = FALSE
  ]
  auto_analysis_limit <- .GenomeGMM_auto_analysis_limit(
    cached_full_data, n = n
  )
  global_pattern_evidence <- detect_global_pattern_peaks(
    cached_full_data, analysis_max_x = auto_analysis_limit,
    allow_post_noise_anchor = .allow_post_noise_anchor
  )
  single_peak_elbow <- .GenomeGMM_single_peak_elbow(
    global_pattern_evidence
  )
  shallow_anchor_active <- isTRUE(
    global_pattern_evidence$shallow_anchor_supported
  ) && is.finite(global_pattern_evidence$shallow_anchor_x) &&
    length(global_pattern_evidence$shallow_anchor_candidate_slots) > 0L
  if (isTRUE(inference_requested)) {
    if (shallow_anchor_active) {
      # A non-integral terminal/anchor ratio has two adjacent biological
      # interpretations. Compare exactly those candidates on one common data
      # range; broad-tail window counters must not manufacture larger
      # ploidies before the first-peak interpretation is resolved.
      n_candidates <- sort(unique(as.integer(
        global_pattern_evidence$shallow_anchor_candidate_slots
      )))
      if (!is.null(.monoploid_range)) {
        # The optional genome-size prior is allowed to distinguish adjacent
        # integer interpretations of the same terminal peak.  The prefix
        # detector can conservatively retain only floor(ratio); add both
        # neighbouring integers when the detected terminal maximum lies
        # within 8% of them. This is based solely on relative peak geometry,
        # not on a sample name, species, absolute depth or target ploidy.
        shallow_ratios <- suppressWarnings(as.numeric(
          global_pattern_evidence$peaks$ratio
        ))
        shallow_ratios <- shallow_ratios[
          is.finite(shallow_ratios) & shallow_ratios >= 1.5
        ]
        adjacent_slots <- unique(as.integer(unlist(lapply(
          shallow_ratios, function(ratio) {
            slots <- unique(c(floor(ratio), ceiling(ratio)))
            slots[
              slots >= 2L & slots <= 30L &
                abs(ratio - slots) / pmax(slots, 1L) <= 0.08
            ]
          }
        ))))
        n_candidates <- sort(unique(c(n_candidates, adjacent_slots)))
      }
    } else {
      dynamic_candidate_cap <- max(
        2L,
        as.integer(global_pattern_evidence$minimum_pattern),
        as.integer(global_pattern_evidence$terminal_slot),
        as.integer(global_pattern_evidence$numerical_terminal_slot),
        as.integer(global_pattern_evidence$supported_fit_slot),
        as.integer(global_pattern_evidence$strong_tail_fit_slot),
        as.integer(global_pattern_evidence$weak_harmonic_fit_slot)
      )
      dynamic_candidate_cap <- min(30L, dynamic_candidate_cap)
      minimum_hint <- max(
        2L, as.integer(global_pattern_evidence$minimum_pattern)
      )
      dynamic_candidate_cap <- max(minimum_hint, dynamic_candidate_cap)
      # Only biological multiples compatible with the minimum pattern need a
      # numerical comparison.  This roughly halves an m=2 screening grid and
      # scales similarly for m=3,4,... without exposing a user parameter.
      n_candidates <- seq.int(
        minimum_hint, dynamic_candidate_cap, by = minimum_hint
      )
      if (!is.null(.monoploid_range)) {
        # Range-guided inference must compare adjacent biological endpoints
        # around the observed terminal structure. The ordinary m, 2m, 3m, ...
        # grid can omit an n + 1 interpretation represented by a fitted tail.
        # Add only the local endpoint neighbourhood, not a full 2:30 scan.
        endpoint_slots <- suppressWarnings(as.integer(c(
          global_pattern_evidence$terminal_slot,
          global_pattern_evidence$numerical_terminal_slot,
          global_pattern_evidence$supported_fit_slot,
          global_pattern_evidence$strong_tail_fit_slot
        )))
        endpoint_slots <- endpoint_slots[
          is.finite(endpoint_slots) & endpoint_slots >= minimum_hint
        ]
        endpoint_neighbours <- unique(as.integer(unlist(lapply(
          endpoint_slots, function(slot) c(slot - 1L, slot, slot + 1L)
        ))))
        endpoint_neighbours <- endpoint_neighbours[
          endpoint_neighbours >= minimum_hint &
            endpoint_neighbours <= dynamic_candidate_cap + 1L
        ]
        n_candidates <- sort(unique(c(n_candidates, endpoint_neighbours)))
      }
    }
  }
  common_extent_slot <- if (shallow_anchor_active) {
    max(n_candidates) + 2L
  } else {
    max(
      2L,
      as.integer(global_pattern_evidence$minimum_pattern),
      as.integer(global_pattern_evidence$terminal_slot),
      as.integer(global_pattern_evidence$numerical_terminal_slot),
      as.integer(global_pattern_evidence$supported_fit_slot),
      as.integer(global_pattern_evidence$strong_tail_fit_slot),
      as.integer(global_pattern_evidence$weak_harmonic_fit_slot)
    )
  }
  common_first_peak <- suppressWarnings(as.numeric(
    global_pattern_evidence$first_peak_x
  ))
  common_range_end <- if (is.finite(common_first_peak) &&
      common_first_peak > 0) {
    min(auto_analysis_limit,
        common_first_peak * (common_extent_slot + 0.5))
  } else {
    auto_analysis_limit
  }
  common_start_row <- suppressWarnings(as.integer(
    global_pattern_evidence$start_index
  ))
  if (!is.finite(common_start_row) || common_start_row < 1L ||
      common_start_row > nrow(cached_full_data)) common_start_row <- 1L
  common_evaluation_data <- cached_full_data[
    seq.int(common_start_row, nrow(cached_full_data)), , drop = FALSE
  ]
  common_evaluation_data <- common_evaluation_data[
    is.finite(common_evaluation_data[, 1]) &
      is.finite(common_evaluation_data[, 2]) &
      common_evaluation_data[, 1] <= common_range_end,
    , drop = FALSE
  ]

  cat("\n=== GenomeGMM automatic parameter search ===\n")
  if (is.null(n)) {
    cat("Stage 1: comparing n candidates:", n_candidates, "\n")
    if (global_pattern_evidence$explicit_peak_count > 0L) {
      cat("Global explicit peaks:",
          paste(round(global_pattern_evidence$peaks$x, 4), collapse = " "),
          "\n")
      cat("Global peak/first-peak ratios:",
          paste(round(global_pattern_evidence$peaks$ratio, 4),
                collapse = " "), "\n")
      cat("Global rightmost trusted harmonic slot:",
          global_pattern_evidence$terminal_slot, "\n")
      cat("Global rightmost explicit numerical harmonic slot:",
          global_pattern_evidence$numerical_terminal_slot, "\n")
      cat("Global tail-supported primary fitting slot:",
          global_pattern_evidence$supported_fit_slot, "\n")
      if (global_pattern_evidence$strong_tail_fit_slot >
          global_pattern_evidence$supported_fit_slot) {
        cat("Global guarded strong-tail alternative slot:",
            global_pattern_evidence$strong_tail_fit_slot, "\n")
      }
      if (isTRUE(
          global_pattern_evidence$weak_harmonic_chain_supported
        ) && global_pattern_evidence$weak_harmonic_fit_slot >
          global_pattern_evidence$strong_tail_fit_slot) {
        cat("Global sustained weak-harmonic candidate slot:",
            global_pattern_evidence$weak_harmonic_fit_slot, "\n")
      }
    }
  } else {
    cat("Using supplied n =", n, "and tuning curve parameters.\n")
  }

  cache <- new.env(parent = emptyenv())
  trial_counter <- 0L
  pattern_selection_table <- NULL
  single_peak_elbow_applied <- FALSE
  single_peak_elbow_guard <- NULL
  range_selected_n <- NA_integer_
  range_terminal_anchor_slot <- NA_integer_
  range_selection_basis <- NA_character_
  range_selected_size_model <- NA_character_

  candidate_key <- function(n_value, ratio, weight, symmetry_value,
                            extra_mode, sigma_mode_value, stage,
                            fit_count_value) {
    ratio_text <- if (is.na(ratio)) "none" else sprintf("%.5f", ratio)
    paste(n_value, ratio_text, sprintf("%.5f", weight),
          sprintf("%.5f", symmetry_value), extra_mode, sigma_mode_value,
          stage, as.integer(fit_count_value), sep = "|")
  }

  assess_component_health <- function(fit, n_value) {
    empty <- list(
      collapsed_biological_components = NA_integer_,
      minimum_biological_share = NA_real_,
      minimum_biological_peak_ratio = NA_real_,
      extra_share = 0,
      extra_peak_ratio = 0,
      extra_after_biological = TRUE,
      extra_guard_passed = TRUE,
      health_penalty = 0
    )
    if (is.null(fit) || !length(fit$alpha) || !length(fit$miu) ||
        !length(fit$sigma)) return(empty)
    component_order <- order(fit$miu)
    primary_count <- if (!is.null(fit$fit_harmonic_count) &&
        is.finite(fit$fit_harmonic_count)) {
      as.integer(fit$fit_harmonic_count)
    } else {
      as.integer(n_value)
    }
    biological_count <- min(primary_count, length(component_order))
    biological <- component_order[seq_len(biological_count)]
    nuisance <- if (length(component_order) > biological_count) {
      component_order[(biological_count + 1L):length(component_order)]
    } else integer(0)
    alpha_total <- sum(pmax(fit$alpha, 0), na.rm = TRUE)
    if (!is.finite(alpha_total) || alpha_total <= 0) return(empty)
    shares <- pmax(fit$alpha, 0) / alpha_total
    peak_strength <- pmax(fit$alpha, 0) / pmax(fit$sigma, 1e-8)
    biological_peak_max <- max(peak_strength[biological], na.rm = TRUE)
    biological_peak_ratio <- peak_strength[biological] /
      max(biological_peak_max, 1e-12)

    # The floor becomes gentler as n grows because a genuine minor dosage
    # class in a high-ploidy spectrum can carry little mass.  It is used only
    # for candidate rejection/penalisation and never forces a peak height.
    share_floor <- min(0.008, 0.04 / max(primary_count, 1L))
    peak_floor <- min(0.020, 0.10 / max(primary_count, 1L))
    collapsed <- shares[biological] < share_floor |
      biological_peak_ratio < peak_floor
    collapsed[!is.finite(collapsed)] <- TRUE

    extra_share <- if (length(nuisance)) sum(shares[nuisance]) else 0
    extra_peak_ratio <- if (length(nuisance)) {
      max(peak_strength[nuisance], na.rm = TRUE) /
        max(biological_peak_max, 1e-12)
    } else 0
    spacing <- if (biological_count > 1L) {
      stats::median(diff(sort(fit$miu[biological])))
    } else max(fit$miu[biological], 1)
    extra_after <- !length(nuisance) ||
      min(fit$miu[nuisance]) >= max(fit$miu[biological]) + 0.15 * spacing
    extra_guard <- !length(nuisance) ||
      (is.finite(extra_share) && extra_share <= 0.30 &&
       is.finite(extra_peak_ratio) && extra_peak_ratio <= 0.65 &&
       isTRUE(extra_after))

    collapse_severity <- sum(collapsed)
    list(
      collapsed_biological_components = as.integer(collapse_severity),
      minimum_biological_share = min(shares[biological], na.rm = TRUE),
      minimum_biological_peak_ratio = min(
        biological_peak_ratio, na.rm = TRUE
      ),
      extra_share = extra_share,
      extra_peak_ratio = extra_peak_ratio,
      extra_after_biological = isTRUE(extra_after),
      extra_guard_passed = isTRUE(extra_guard),
      health_penalty = 0.04 * collapse_severity +
        if (length(nuisance)) 0.0015 else 0
    )
  }

  # Compare different K values on one fixed data domain.  Candidate-local
  # NRMSE is still used to tune nuisance parameters within a fixed K, but it
  # cannot compare K=2 (whose plot may stop near 4x) against K=10 (whose plot
  # covers the complete tail).  These metrics therefore use the same observed
  # bins and normalization for every candidate.
  assess_common_range <- function(fit) {
    empty <- list(
      nrmse = Inf, underfit_nrmse = Inf, underfit_area = Inf,
      tail_underfit_area = Inf, last_width_inflation = Inf,
      maximum_width_jump = Inf, score = Inf
    )
    if (is.null(fit) || nrow(common_evaluation_data) < 5L ||
        !length(fit$alpha) || !length(fit$miu) || !length(fit$sigma) ||
        is.null(fit$total_samples_fit) ||
        !is.finite(fit$total_samples_fit) || fit$total_samples_fit <= 0) {
      return(empty)
    }
    x_common <- as.numeric(common_evaluation_data[, 1])
    y_common <- pmax(as.numeric(common_evaluation_data[, 2]), 0)
    x_steps <- diff(sort(unique(x_common)))
    x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
    bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
    component_counts <- vapply(seq_along(fit$alpha), function(j) {
      fit$alpha[j] * stats::dnorm(
        x_common, fit$miu[j], fit$sigma[j]
      ) * fit$total_samples_fit * bin_width
    }, numeric(length(x_common)))
    prediction <- if (is.null(dim(component_counts))) {
      component_counts
    } else {
      rowSums(component_counts)
    }
    prediction[!is.finite(prediction)] <- 0
    observed <- as.numeric(stats::filter(
      y_common, rep(1 / 5, 5), sides = 2
    ))
    observed[!is.finite(observed)] <- y_common[!is.finite(observed)]
    y_scale <- max(observed, na.rm = TRUE)
    if (!is.finite(y_scale) || y_scale <= 0) return(empty)
    weights <- 0.25 + 0.75 * sqrt(pmax(observed, 0) / y_scale)
    scaled_error <- (prediction - observed) / y_scale
    common_nrmse <- sqrt(sum(weights * scaled_error^2) / sum(weights))
    common_underfit_nrmse <- sqrt(
      sum(weights * pmax(-scaled_error, 0)^2) / sum(weights)
    )
    common_underfit_area <- sum(pmax(observed - prediction, 0)) /
      max(sum(observed), 1)

    primary_count <- min(
      as.integer(fit$fit_harmonic_count), length(fit$miu)
    )
    component_order <- order(fit$miu)
    primary <- component_order[seq_len(primary_count)]
    widths <- pmax(fit$sigma[primary], 1e-8)
    preceding_width <- if (primary_count > 1L) {
      stats::median(tail(widths[-primary_count], 3L), na.rm = TRUE)
    } else widths[1]
    last_width_inflation <- widths[primary_count] /
      max(preceding_width, 1e-8)
    maximum_width_jump <- if (primary_count > 1L) {
      max(widths[-1L] / widths[-primary_count], na.rm = TRUE)
    } else 1

    tail_start <- common_first_peak * (primary_count + 0.5)
    tail_indices <- which(x_common > tail_start)
    tail_underfit_area <- if (length(tail_indices)) {
      sum(pmax(observed[tail_indices] - prediction[tail_indices], 0)) /
        max(sum(observed), 1)
    } else 0
    width_penalty <-
      pmax(log(last_width_inflation / 1.35), 0)^2 +
      0.5 * pmax(log(maximum_width_jump / 1.50), 0)^2
    complexity_penalty <- 0.0005 * primary_count
    common_score <- 0.55 * common_nrmse +
      0.25 * common_underfit_nrmse +
      0.10 * common_underfit_area +
      0.03 * width_penalty + complexity_penalty
    list(
      nrmse = common_nrmse,
      underfit_nrmse = common_underfit_nrmse,
      underfit_area = common_underfit_area,
      tail_underfit_area = tail_underfit_area,
      last_width_inflation = last_width_inflation,
      maximum_width_jump = maximum_width_jump,
      score = common_score
    )
  }

  run_candidate <- function(n_value, ratio = NA_real_, weight = 1,
                            symmetry_value = 0.2, extra_mode = "never",
                            sigma_mode_value = "adaptive",
                            stage = c("screening", "verification"),
                            initial_fit = NULL,
                            fit_count_override = NULL) {
    stage <- match.arg(stage)
    candidate_fit_count <- if (!is.null(fit_count_override) &&
        length(fit_count_override) == 1L &&
        is.finite(fit_count_override)) {
      max(as.integer(n_value), as.integer(fit_count_override))
    } else if (shallow_anchor_active) {
      as.integer(n_value)
    } else {
      max(as.integer(n_value),
          as.integer(global_pattern_evidence$supported_fit_slot))
    }
    key <- candidate_key(n_value, ratio, weight, symmetry_value,
                         extra_mode, sigma_mode_value, stage,
                         candidate_fit_count)
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
      .global_pattern_evidence = global_pattern_evidence,
      .fit_harmonic_count = candidate_fit_count,
      .initial_fit = initial_fit,
      .first_peak_override = if (shallow_anchor_active) {
        global_pattern_evidence$shallow_anchor_x
      } else NULL,
      .start_index_override = if (shallow_anchor_active) {
        global_pattern_evidence$start_index
      } else NULL,
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
    health <- assess_component_health(fit, n_value)
    common <- if (succeeded) assess_common_range(fit) else
      assess_common_range(NULL)
    score <- raw_fit_score + health$health_penalty
    if (succeeded && !isTRUE(fit$converged)) score <- score + 0.01

    # The common-domain screening curve may contain more primary harmonics
    # than its biological candidate. For candidate m, Monoploid genome size
    # is nevertheless total genomic k-mers divided by the m-th harmonic
    # depth—not by the last numerical component in that screening curve.
    candidate_monoploid_gse_bp <- NA_real_
    if (succeeded) {
      ordered_candidate_miu <- sort(as.numeric(
        fit$miu[is.finite(fit$miu) & fit$miu > 0]
      ))
      candidate_index <- as.integer(n_value)
      candidate_total_kmers <- suppressWarnings(as.numeric(fit$total_kmers))
      if ((!length(candidate_total_kmers) ||
           !is.finite(candidate_total_kmers[1L]) ||
           candidate_total_kmers[1L] <= 0) &&
          !is.null(fit$g_hap_bp) && is.finite(fit$g_hap_bp) &&
          is.finite(fit$fit_harmonic_count) &&
          fit$fit_harmonic_count <= length(ordered_candidate_miu)) {
        candidate_total_kmers <- as.numeric(fit$g_hap_bp) *
          ordered_candidate_miu[as.integer(fit$fit_harmonic_count)]
      }
      if (length(candidate_total_kmers) &&
          is.finite(candidate_total_kmers[1L]) &&
          candidate_total_kmers[1L] > 0 && candidate_index >= 1L &&
          candidate_index <= length(ordered_candidate_miu)) {
        candidate_monoploid_gse_bp <- candidate_total_kmers[1L] /
          ordered_candidate_miu[candidate_index]
      }
    }

    row <- data.frame(
      trial = trial_counter,
      stage = stage,
      n = as.integer(n_value),
      fit_harmonic_count = if (succeeded) {
        as.integer(fit$fit_harmonic_count)
      } else candidate_fit_count,
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
      observed_peak_count = if (succeeded) {
        as.integer(fit$observed_peak_count)
      } else NA_integer_,
      global_pattern_peak_count = if (succeeded) {
        as.integer(fit$global_pattern_peak_count)
      } else NA_integer_,
      global_pattern_trusted_count = if (succeeded) {
        as.integer(fit$global_pattern_trusted_count)
      } else NA_integer_,
      terminal_harmonic_slot = if (succeeded) {
        as.integer(fit$terminal_harmonic_slot)
      } else NA_integer_,
      collapsed_biological_components =
        health$collapsed_biological_components,
      minimum_biological_share = health$minimum_biological_share,
      minimum_biological_peak_ratio = health$minimum_biological_peak_ratio,
      extra_share = health$extra_share,
      extra_peak_ratio = health$extra_peak_ratio,
      extra_after_biological = health$extra_after_biological,
      extra_guard_passed = health$extra_guard_passed,
      nrmse = if (succeeded) fit$fit_metrics$nrmse else Inf,
      main_peak_nrmse = if (succeeded) fit$fit_metrics$main_peak_nrmse else Inf,
      left_overshoot_error = if (succeeded) {
        fit$fit_metrics$left_overshoot_error
      } else Inf,
      left_relative_error = if (succeeded) {
        fit$fit_metrics$left_relative_error
      } else Inf,
      common_nrmse = common$nrmse,
      common_underfit_nrmse = common$underfit_nrmse,
      common_underfit_area = common$underfit_area,
      common_tail_underfit_area = common$tail_underfit_area,
      last_width_inflation = common$last_width_inflation,
      maximum_width_jump = common$maximum_width_jump,
      common_range_score = common$score,
      monoploid_gse_bp = candidate_monoploid_gse_bp,
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
    initial_values <- if (shallow_anchor_active) {
      lapply(n_candidates, function(candidate_n) {
        run_candidate(
          candidate_n, extra_mode = "legacy",
          fit_count_override = candidate_n
        )
      })
    } else {
      lapply(n_candidates, run_candidate)
    }
    pattern_selection_table <- all_rows()
    pattern_selection_table <- pattern_selection_table[
      pattern_selection_table$stage == "screening" &
        pattern_selection_table$n %in% n_candidates,
      , drop = FALSE
    ]
    pattern_selection_table <- pattern_selection_table[
      !duplicated(pattern_selection_table$n), , drop = FALSE
    ]
    # The minimum pattern is the first candidate whose fitting window contains
    # the rightmost trustworthy integer harmonic slot. This preserves missing
    # internal slots: peaks at {1,3} imply minimum pattern 3, and {1,2,4}
    # imply minimum pattern 4. Raw peak count alone would incorrectly call both
    # examples 2 or 3. A later non-integer repeat peak is ignored here.
    detection_rows <- pattern_selection_table[
      is.finite(pattern_selection_table$terminal_harmonic_slot), , drop = FALSE
    ]
    detection_rows <- detection_rows[order(detection_rows$n), , drop = FALSE]
    if (nrow(detection_rows)) {
      detected_count <- NA_integer_
      for (i in seq_len(nrow(detection_rows))) {
        n_here <- as.integer(detection_rows$n[i])
        terminal_here <- as.integer(
          detection_rows$terminal_harmonic_slot[i]
        )
        if (terminal_here <= n_here) {
          detected_count <- n_here
          break
        }
      }
      if (!is.finite(detected_count)) {
        last_row <- detection_rows[nrow(detection_rows), , drop = FALSE]
        detected_count <- min(as.integer(last_row$n),
                              as.integer(last_row$terminal_harmonic_slot))
      }
    } else {
      detected_count <- 2L
    }

    # Candidate-specific peak reconstruction is useful for resolving broad,
    # overlapping harmonics, but it must not turn an isolated low-mass repeat
    # peak into a higher minimum biological pattern.  Bound that reconstructed
    # endpoint by the furthest slot supported by one of the independent
    # whole-spectrum evidence paths: explicit trusted harmonics, consecutive
    # harmonic windows, or a sustained strong/weak tail chain.  A distant
    # numerical maximum remains available to set the common comparison range,
    # but cannot raise the inferred minimum by itself.  This distinction is
    # deliberately based on relative signal continuity rather than species,
    # absolute coverage, or a hard-coded target ploidy.
    spectrum_supported_cap <- max(
      2L,
      as.integer(global_pattern_evidence$minimum_pattern),
      as.integer(global_pattern_evidence$supported_fit_slot),
      as.integer(global_pattern_evidence$strong_tail_fit_slot),
      if (isTRUE(global_pattern_evidence$weak_harmonic_chain_supported)) {
        as.integer(global_pattern_evidence$weak_harmonic_fit_slot)
      } else {
        as.integer(global_pattern_evidence$strong_tail_fit_slot)
      }
    )
    if (is.finite(detected_count) && is.finite(spectrum_supported_cap)) {
      detected_count <- min(
        as.integer(detected_count), as.integer(spectrum_supported_cap)
      )
    }
    if (shallow_anchor_active) {
      shallow_rows <- pattern_selection_table[
        pattern_selection_table$n %in% n_candidates &
          is.finite(pattern_selection_table$common_range_score) &
          pattern_selection_table$collapsed_biological_components == 0L,
        , drop = FALSE
      ]
      if (!nrow(shallow_rows)) {
        shallow_rows <- pattern_selection_table[
          pattern_selection_table$n %in% n_candidates &
            is.finite(pattern_selection_table$common_range_score),
          , drop = FALSE
        ]
      }
      if (!nrow(shallow_rows)) {
        stop("All shallow-anchor candidate fits failed.")
      }
      shallow_rows <- shallow_rows[order(
        shallow_rows$common_range_score,
        shallow_rows$common_nrmse,
        shallow_rows$n
      ), , drop = FALSE]
      shallow_guard_passed <- isTRUE(shallow_rows$converged[1]) &&
        is.finite(shallow_rows$common_nrmse[1]) &&
        shallow_rows$common_nrmse[1] <= 0.085 &&
        is.finite(shallow_rows$common_range_score[1]) &&
        shallow_rows$common_range_score[1] <= 0.085 &&
        is.finite(shallow_rows$extra_share[1]) &&
        shallow_rows$extra_share[1] <= 0.30
      if (isTRUE(.allow_post_noise_anchor) &&
          !isTRUE(shallow_guard_passed)) {
        cat("Post-noise first-peak candidate was rejected by the",
            "common-range fit guard; retaining the ordinary start.\n")
        return(.GenomeGMM_auto(
          data_file = data_file, n = NULL, n_candidates = 2:10,
          search = search, symmetry_tolerance = symmetry_tolerance,
          position_mode = position_mode,
          joint_refinement = joint_refinement,
          species_name = species_name, seed = seed,
          save_search = save_search,
          candidate_max_iter = candidate_max_iter,
          candidate_tolerance = candidate_tolerance,
          verification_top = verification_top,
          verification_max_iter = verification_max_iter,
          final_max_iter = final_max_iter,
          save_plot = save_plot,
          .selection_only = .selection_only,
          .preloaded_data = cached_full_data,
          .allow_post_noise_anchor = FALSE,
          .monoploid_range = .monoploid_range
        ))
      }
      n <- as.integer(shallow_rows$n[1])
      cat("Post-noise shallow first peak retained at",
          round(global_pattern_evidence$shallow_anchor_x, 3),
          "; adjacent terminal slots",
          paste(n_candidates, collapse = "/"),
          "were compared on one observed range; selected", n, "x.\n")
    } else {
      n <- max(2L, min(max(n_candidates), detected_count))
      cat("Selected minimum observable pattern =", n,
          "from the rightmost trusted harmonic slot for detailed tuning.\n")
    }

    # A single explicit peak cannot by itself prove that a long numerical tail
    # is high ploidy.  If harmonic-window support has a guarded plateau elbow,
    # test K directly with one optional tail component.  Prefer K only when the
    # resulting curve is already accurate in absolute terms; this prevents a
    # lower-complexity biological interpretation from being accepted merely
    # because it was proposed by the elbow detector.
    if (isTRUE(single_peak_elbow$detected)) {
      elbow_k <- as.integer(single_peak_elbow$slot)
      elbow_trial <- run_candidate(
        elbow_k, ratio = NA_real_, weight = 1,
        symmetry_value = if (is.null(symmetry_tolerance)) {
          0.2
        } else symmetry_tolerance,
        extra_mode = "legacy",
        sigma_mode_value = "adaptive", stage = "verification",
        fit_count_override = elbow_k
      )
      elbow_row <- elbow_trial$row
      elbow_fit <- elbow_trial$fit
      width_ratio_for <- function(fit) {
        if (is.null(fit) || !length(fit$sigma)) return(Inf)
        component_order <- order(fit$miu)
        primary_count <- min(elbow_k, length(component_order))
        widths <- fit$sigma[component_order[seq_len(primary_count)]]
        if (length(widths) < 2L) return(Inf)
        widths[length(widths)] /
          max(stats::median(tail(widths[-length(widths)], 3L)), 1e-8)
      }
      terminal_width_ratio <- width_ratio_for(elbow_fit)
      guard_source <- "direct protected fit"
      guard_nrmse <- elbow_row$nrmse[1]
      guard_main_nrmse <- elbow_row$main_peak_nrmse[1]
      guard_left_overshoot <- elbow_row$left_overshoot_error[1]
      guard_collapsed <- elbow_row$collapsed_biological_components[1]
      guard_extra_ok <- elbow_row$extra_guard_passed[1]
      guard_converged <- !is.null(elbow_fit) && isTRUE(elbow_fit$converged)
      left_overshoot_limit <- function() {
        if (!is.finite(guard_main_nrmse)) return(NA_real_)
        min(0.006, max(0.005, 0.25 * guard_main_nrmse))
      }
      passes_guard <- function() {
        isTRUE(guard_converged) &&
          is.finite(guard_nrmse) && guard_nrmse <= 0.025 &&
          is.finite(guard_main_nrmse) && guard_main_nrmse <= 0.030 &&
          is.finite(guard_left_overshoot) &&
          is.finite(left_overshoot_limit()) &&
          guard_left_overshoot <= left_overshoot_limit() &&
          isTRUE(guard_collapsed == 0L) && isTRUE(guard_extra_ok) &&
          is.finite(terminal_width_ratio) && terminal_width_ratio <= 1.50
      }
      quality_passed <- passes_guard()

      # A fixed preliminary configuration can reject a valid K even though
      # normal automatic tuning readily satisfies every protected metric.  In
      # that case tune K through the same known-n pathway used for the final
      # curve, without a plot, and retain the strict thresholds unchanged.
      if (!isTRUE(quality_passed)) {
        tuned_elbow_fit <- NULL
        tryCatch(
          invisible(utils::capture.output(
            tuned_elbow_fit <- .GenomeGMM_auto(
              data_file = data_file, n = elbow_k,
              n_candidates = n_candidates, search = search,
              symmetry_tolerance = symmetry_tolerance,
              position_mode = position_mode,
              joint_refinement = joint_refinement,
              species_name = species_name, seed = seed,
              save_search = FALSE,
              candidate_max_iter = candidate_max_iter,
              candidate_tolerance = candidate_tolerance,
              verification_top = verification_top,
              verification_max_iter = verification_max_iter,
              final_max_iter = final_max_iter,
              save_plot = FALSE, .selection_only = FALSE,
              .preloaded_data = .preloaded_data
            ),
            type = "output"
          )),
          error = function(error) tuned_elbow_fit <<- NULL
        )
        if (!is.null(tuned_elbow_fit)) {
          tuned_health <- assess_component_health(tuned_elbow_fit, elbow_k)
          guard_source <- "fully tuned protected fit"
          guard_nrmse <- tuned_elbow_fit$fit_metrics$nrmse
          guard_main_nrmse <-
            tuned_elbow_fit$fit_metrics$main_peak_nrmse
          guard_left_overshoot <-
            tuned_elbow_fit$fit_metrics$left_overshoot_error
          guard_collapsed <- tuned_health$collapsed_biological_components
          guard_extra_ok <- tuned_health$extra_guard_passed
          guard_converged <- isTRUE(tuned_elbow_fit$converged)
          terminal_width_ratio <- width_ratio_for(tuned_elbow_fit)
          quality_passed <- passes_guard()
        }
      }
      single_peak_elbow_guard <- list(
        slot = elbow_k, passed = isTRUE(quality_passed),
        source = guard_source, nrmse = guard_nrmse,
        main_peak_nrmse = guard_main_nrmse,
        left_overshoot_error = guard_left_overshoot,
        left_overshoot_limit = left_overshoot_limit(),
        terminal_width_ratio = terminal_width_ratio
      )
      if (isTRUE(quality_passed)) {
        n <- elbow_k
        single_peak_elbow_applied <- TRUE
        cat("Single-peak plateau elbow accepted at", elbow_k,
            "harmonics after protected curve verification.\n")
      } else {
        cat("Single-peak plateau elbow at", elbow_k,
            "was rejected by protected curve metrics.\n")
      }
    }

    # When the user supplies prior knowledge of the Monoploid genome size,
    # use it as a hard admissibility interval before the final PDF fit. The
    # interval does not define a preferred midpoint. Among healthy candidates
    # inside it, prefer the base pattern aligned to the last explicit genomic
    # peak/shoulder before sustained terminal decay. This is a histogram-wide
    # structural rule: it is independent of species and absolute count scale.
    # If no reliable terminal anchor exists, retain the established protected
    # common-range fit score. Same-pattern multiples remain report-only.
    if (!is.null(.monoploid_range)) {
      lower_bp <- min(.monoploid_range) * 1e6
      upper_bp <- max(.monoploid_range) * 1e6
      range_candidate_pool <- pattern_selection_table
      range_candidate_pool$range_base_ploidy <- as.integer(
        range_candidate_pool$n
      )
      range_candidate_pool$range_structure_ploidy <- as.integer(
        range_candidate_pool$n
      )
      range_candidate_pool$direct_monoploid_bp <-
        range_candidate_pool$monoploid_gse_bp
      range_candidate_pool$allopolyploid_monoploid_bp <-
        range_candidate_pool$monoploid_gse_bp *
        range_candidate_pool$range_base_ploidy
      range_candidate_pool$direct_range_match <-
        is.finite(range_candidate_pool$direct_monoploid_bp) &
        range_candidate_pool$direct_monoploid_bp >= lower_bp &
        range_candidate_pool$direct_monoploid_bp <= upper_bp
      range_candidate_pool$allopolyploid_range_match <-
        is.finite(range_candidate_pool$allopolyploid_monoploid_bp) &
        range_candidate_pool$allopolyploid_monoploid_bp >= lower_bp &
        range_candidate_pool$allopolyploid_monoploid_bp <= upper_bp

      eligible_range_rows <- range_candidate_pool[
        (range_candidate_pool$direct_range_match |
           range_candidate_pool$allopolyploid_range_match) &
          range_candidate_pool$converged &
          is.finite(range_candidate_pool$common_range_score),
        , drop = FALSE
      ]
      if (nrow(eligible_range_rows)) {
        # Retain one best-scoring fit for each biological base candidate even
        # if several numerical tail configurations were screened for it.
        eligible_range_rows <- eligible_range_rows[order(
          eligible_range_rows$common_range_score,
          eligible_range_rows$n
        ), , drop = FALSE]
        eligible_range_rows <- eligible_range_rows[
          !duplicated(eligible_range_rows$range_base_ploidy),
          , drop = FALSE
        ]
      }
      healthy_range_rows <- eligible_range_rows[
        is.finite(eligible_range_rows$collapsed_biological_components) &
          eligible_range_rows$collapsed_biological_components == 0L,
        , drop = FALSE
      ]
      if (nrow(healthy_range_rows)) {
        eligible_range_rows <- healthy_range_rows
      }
      if (nrow(eligible_range_rows)) {
        evidence_peaks <- global_pattern_evidence$peaks
        explicit_slots <- integer(0)
        if (!is.null(evidence_peaks) && nrow(evidence_peaks)) {
          if ("ratio" %in% names(evidence_peaks)) {
            peak_ratios <- suppressWarnings(as.numeric(evidence_peaks$ratio))
            explicit_slots <- as.integer(round(peak_ratios[
              is.finite(peak_ratios) & peak_ratios >= 1.5
            ]))
          } else if ("slot" %in% names(evidence_peaks)) {
            explicit_slots <- suppressWarnings(as.integer(
              evidence_peaks$slot
            ))
            explicit_slots <- explicit_slots[
              is.finite(explicit_slots) & explicit_slots >= 2L
            ]
          }
        }
        support_table <- global_pattern_evidence$harmonic_window_support
        sustained_decay_after <- function(slot) {
          if (is.null(support_table) || !nrow(support_table) ||
              !all(c("slot", "median_fraction") %in%
                   names(support_table))) return(FALSE)
          anchor_index <- match(slot, support_table$slot)
          if (is.na(anchor_index) ||
              anchor_index + 2L > nrow(support_table)) return(FALSE)
          end_index <- min(nrow(support_table), anchor_index + 3L)
          values <- suppressWarnings(as.numeric(
            support_table$median_fraction[anchor_index:end_index]
          ))
          if (length(values) < 3L || any(!is.finite(values)) ||
              values[1L] <= 0) return(FALSE)
          steps <- diff(values)
          nonincreasing_steps <- sum(steps <= 0.01)
          nonincreasing_steps >= length(steps) - 1L &&
            tail(values, 1L) <= 0.80 * values[1L]
        }
        explicit_slots <- sort(unique(explicit_slots[
          is.finite(explicit_slots) & explicit_slots >= 2L
        ]))
        decay_supported_slots <- explicit_slots[vapply(
          explicit_slots, sustained_decay_after, logical(1)
        )]
        if (length(decay_supported_slots)) {
          range_terminal_anchor_slot <- max(decay_supported_slots)
          range_selection_basis <- paste0(
            "terminal supported harmonic ", range_terminal_anchor_slot,
            " before sustained decay"
          )
        } else {
          numerical_slot <- suppressWarnings(as.integer(
            global_pattern_evidence$numerical_terminal_slot
          ))
          if (length(numerical_slot) == 1L &&
              is.finite(numerical_slot) && numerical_slot >= 2L &&
              sustained_decay_after(numerical_slot)) {
            range_terminal_anchor_slot <- numerical_slot
            range_selection_basis <- paste0(
              "terminal supported harmonic ", range_terminal_anchor_slot,
              " before sustained decay"
            )
          }
        }
        if (!is.finite(range_terminal_anchor_slot)) {
          supported_extent <- suppressWarnings(as.integer(
            global_pattern_evidence$supported_fit_slot
          ))
          if (length(supported_extent) == 1L &&
              is.finite(supported_extent) && supported_extent >= 2L &&
              supported_extent %in%
                eligible_range_rows$range_base_ploidy) {
            range_terminal_anchor_slot <- supported_extent
            range_selection_basis <- paste0(
              "full-spectrum tail-supported harmonic extent ",
              range_terminal_anchor_slot
            )
          }
        }
        if (is.finite(range_terminal_anchor_slot)) {
          eligible_range_rows$terminal_distance <- abs(
            eligible_range_rows$range_structure_ploidy -
              range_terminal_anchor_slot
          )
          eligible_range_rows <- eligible_range_rows[order(
            eligible_range_rows$terminal_distance,
            eligible_range_rows$common_range_score,
            eligible_range_rows$range_base_ploidy
          ), , drop = FALSE]
          if (is.na(range_selection_basis) ||
              !nzchar(range_selection_basis)) {
            range_selection_basis <- paste0(
              "terminal supported harmonic ", range_terminal_anchor_slot
            )
          }
        } else {
          eligible_range_rows <- eligible_range_rows[order(
            eligible_range_rows$common_range_score,
            eligible_range_rows$range_base_ploidy
          ), , drop = FALSE]
          range_selection_basis <-
            "protected common-range fit score (no terminal anchor)"
        }
        range_selected_n <- as.integer(
          eligible_range_rows$range_base_ploidy[1L]
        )
        selected_direct <- isTRUE(
          eligible_range_rows$direct_range_match[1L]
        )
        selected_allopolyploid <- isTRUE(
          eligible_range_rows$allopolyploid_range_match[1L]
        )
        range_selected_size_model <- if (selected_allopolyploid &&
            !selected_direct) {
          "balanced_allopolyploid"
        } else if (selected_direct && !selected_allopolyploid) {
          "direct_or_autopolyploid"
        } else {
          "ambiguous"
        }
        n <- range_selected_n
        selected_monoploid_bp <- if (identical(
            range_selected_size_model, "balanced_allopolyploid"
          )) {
          eligible_range_rows$allopolyploid_monoploid_bp[1L]
        } else {
          eligible_range_rows$direct_monoploid_bp[1L]
        }
        cat("Monoploid-range prior selected base pattern =", n,
            "x; candidate Monoploid genome size =",
            sprintf("%.2f Mb",
                    selected_monoploid_bp / 1e6),
            "; interpretation =", range_selected_size_model,
            "\n")
      } else {
        cat("No protected base-pattern candidate falls within the supplied",
            "Monoploid genome-size range; retaining the automatic fit.\n")
      }
    }
  } else {
    pattern_selection_table <- data.frame(
      n = as.integer(n), fit_score = NA_real_, score = NA_real_,
      converged = NA, status = "user-supplied minimum pattern",
      stringsAsFactors = FALSE
    )
  }

  # Keep biological interpretation separate from numerical curve fitting.
  # Once a fitting extent K is selected, every downstream search, constraint,
  # refinement and plot must use that same K, exactly as a direct n = K run.
  # The minimum pattern is metadata for the ambiguity report only.
  minimum_pattern <- as.integer(n)
  fitting_ploidy <- if (is.finite(range_selected_n)) {
    minimum_pattern
  } else if (shallow_anchor_active) {
    minimum_pattern
  } else if (isTRUE(single_peak_elbow_applied)) {
    minimum_pattern
  } else {
    max(
      minimum_pattern,
      as.integer(global_pattern_evidence$supported_fit_slot)
    )
  }
  numerical_tail_fit_slot <- max(
    as.integer(global_pattern_evidence$strong_tail_fit_slot),
    if (isTRUE(global_pattern_evidence$weak_harmonic_chain_supported)) {
      as.integer(global_pattern_evidence$weak_harmonic_fit_slot)
    } else {
      as.integer(global_pattern_evidence$strong_tail_fit_slot)
    }
  )

  # Confidence in the minimum observable pattern says whether explicit peaks
  # support m; it does not say that the numerical curve ends at m.  Whenever
  # full-spectrum windows show a longer structured tail, compare compatible
  # multiples of m on the same x domain.  This prevents a K=2 fit from looking
  # artificially good merely because its own fitting/plotting range stopped
  # before the later overlapping harmonics.
  if (isTRUE(inference_requested) &&
      !is.finite(range_selected_n) &&
      !shallow_anchor_active &&
      !isTRUE(single_peak_elbow_applied) &&
      is.finite(numerical_tail_fit_slot) &&
      numerical_tail_fit_slot >=
        fitting_ploidy + minimum_pattern) {
    extent_rows <- pattern_selection_table[
      pattern_selection_table$stage == "screening" &
        is.finite(pattern_selection_table$fit_harmonic_count) &
        pattern_selection_table$fit_harmonic_count >= fitting_ploidy &
        pattern_selection_table$fit_harmonic_count <=
          numerical_tail_fit_slot &
        pattern_selection_table$fit_harmonic_count %% minimum_pattern == 0L &
        is.finite(pattern_selection_table$common_nrmse) &
        is.finite(pattern_selection_table$common_range_score),
      , drop = FALSE
    ]
    extent_rows <- extent_rows[
      !duplicated(extent_rows$fit_harmonic_count), , drop = FALSE
    ]
    baseline_rows <- pattern_selection_table[
      pattern_selection_table$stage == "screening" &
        pattern_selection_table$fit_harmonic_count == fitting_ploidy &
        is.finite(pattern_selection_table$common_nrmse) &
        is.finite(pattern_selection_table$common_range_score),
      , drop = FALSE
    ]
    if (nrow(baseline_rows) && nrow(extent_rows)) {
      baseline <- baseline_rows[
        which.min(baseline_rows$common_range_score), , drop = FALSE
      ]
      extent_rows$common_nrmse_gain <- 1 - extent_rows$common_nrmse /
        max(baseline$common_nrmse[1], 1e-8)
      extent_rows$common_score_gain <- 1 -
        extent_rows$common_range_score /
          max(baseline$common_range_score[1], 1e-8)
      extent_rows$underfit_area_gain <- 1 -
        extent_rows$common_underfit_area /
          max(baseline$common_underfit_area[1], 1e-8)
      main_guard <- extent_rows$main_peak_nrmse <=
        baseline$main_peak_nrmse[1] * 1.15 + 0.003
      left_guard <- extent_rows$left_overshoot_error <=
        baseline$left_overshoot_error[1] * 1.15 + 0.002
      exceptional_global_evidence <-
        extent_rows$common_nrmse_gain >= 0.30 &
        extent_rows$common_score_gain >= 0.25
      healthy_components <-
        extent_rows$collapsed_biological_components == 0L &
        extent_rows$converged
      healthier_width <- extent_rows$last_width_inflation <=
        pmax(1.65, baseline$last_width_inflation[1] * 0.90)
      advances_extent <- extent_rows$fit_harmonic_count > fitting_ploidy
      eligible <- advances_extent & healthy_components & main_guard &
        (left_guard | exceptional_global_evidence) &
        healthier_width &
        extent_rows$common_nrmse_gain >= 0.08 &
        extent_rows$common_score_gain >= 0.10 &
        extent_rows$underfit_area_gain >= 0.10
      eligible[is.na(eligible)] <- FALSE
      if (any(eligible)) {
        accepted_rows <- extent_rows[eligible, , drop = FALSE]
        numerical_terminal <- suppressWarnings(as.integer(
          global_pattern_evidence$numerical_terminal_slot
        ))
        explicit_extent_anchor <- is.finite(numerical_terminal) &&
          numerical_terminal > fitting_ploidy &&
          numerical_terminal %% minimum_pattern == 0L
        anchor_slot <- if (isTRUE(explicit_extent_anchor)) {
          numerical_terminal
        } else fitting_ploidy
        anchor_rows <- accepted_rows[
          accepted_rows$fit_harmonic_count == anchor_slot,
          , drop = FALSE
        ]
        support_table <- global_pattern_evidence$harmonic_window_support
        post_anchor_rebound <- TRUE
        if (is.finite(anchor_slot) && !is.null(support_table) &&
            nrow(support_table) >= anchor_slot) {
          post_rows <- support_table[
            support_table$slot >= anchor_slot &
              support_table$slot <=
                numerical_tail_fit_slot,
            , drop = FALSE
          ]
          post_values <- post_rows$median_fraction
          post_values <- post_values[is.finite(post_values)]
          # A rise exceeding 1.5% of peak scale after the explicit terminal
          # maximum indicates another unresolved structure.  Otherwise the
          # remaining signal is a decay wing for the optional tail component.
          post_anchor_rebound <- length(post_values) >= 2L &&
            any(diff(post_values) > 0.015)
        }
        best_score <- min(accepted_rows$common_range_score)
        near_best <- accepted_rows$common_range_score <= best_score * 1.05
        best_extensions <- accepted_rows[near_best, , drop = FALSE]
        best_extensions <- best_extensions[
          order(best_extensions$fit_harmonic_count), , drop = FALSE
        ]
        best_extension <- best_extensions[1, , drop = FALSE]
        sustained_weak_chain_winner <- isTRUE(
          global_pattern_evidence$weak_harmonic_chain_supported
        ) && best_extension$fit_harmonic_count[1] ==
          global_pattern_evidence$weak_harmonic_fit_slot
        if (isTRUE(sustained_weak_chain_winner)) {
          chosen_extent <- best_extension
          cat("  sustained weak harmonics improved the common-domain fit;",
              "retaining them as primary curve components\n")
        } else if (isTRUE(explicit_extent_anchor) && nrow(anchor_rows) &&
            !isTRUE(post_anchor_rebound)) {
          chosen_extent <- anchor_rows[
            which.min(anchor_rows$common_range_score), , drop = FALSE
          ]
          cat("  explicit terminal harmonic followed only by decay;",
              "keeping it as the last primary component\n")
        } else if (!isTRUE(explicit_extent_anchor) &&
            !isTRUE(post_anchor_rebound)) {
          exceptional_tail_decomposition <-
            best_extension$common_score_gain[1] >= 0.50 &&
            baseline$last_width_inflation[1] >= 1.60 &&
            best_extension$last_width_inflation[1] <=
              baseline$last_width_inflation[1] * 0.85
          if (isTRUE(exceptional_tail_decomposition)) {
            chosen_extent <- best_extension
            cat("  no later explicit peak, but severe terminal-width",
                "inflation was rescued by the extended model\n")
          } else {
            chosen_extent <- baseline
            chosen_extent$common_nrmse_gain <- 0
            chosen_extent$common_score_gain <- 0
            cat("  later windows form only a decaying wing;",
                "keeping the strictly supported fitting extent\n")
          }
        } else {
          chosen_extent <- best_extension
        }
        fitting_ploidy <- as.integer(chosen_extent$fit_harmonic_count[1])
        cat("Minimum pattern retained as", minimum_pattern,
            "x; full-spectrum numerical extent", fitting_ploidy,
            "was selected for complete curve tuning.\n")
        cat("  common-range NRMSE gain =",
            sprintf("%.1f%%", 100 * chosen_extent$common_nrmse_gain[1]),
            "; common-range score gain =",
            sprintf("%.1f%%", 100 * chosen_extent$common_score_gain[1]),
            "\n")
      }
    }
  }
  if (fitting_ploidy > minimum_pattern) {
    cat("Minimum possible ploidy remains", minimum_pattern,
        "but the complete curve will be fitted as", fitting_ploidy, "x.\n")
  }
  n <- fitting_ploidy

  if (isTRUE(inference_requested) && isTRUE(.selection_only)) {
    trusted_peaks <- global_pattern_evidence$peaks
    if (!is.null(trusted_peaks) && nrow(trusted_peaks) &&
        "trusted_harmonic" %in% names(trusted_peaks)) {
      trusted_peaks <- trusted_peaks[
        trusted_peaks$trusted_harmonic, , drop = FALSE
      ]
    }
    return(invisible(list(
      minimum_pattern = minimum_pattern,
      minimum_observable_pattern = minimum_pattern,
      fit_harmonic_count = fitting_ploidy,
      minimum_pattern_candidates = pattern_selection_table,
      global_pattern_peak_count =
        as.integer(global_pattern_evidence$explicit_peak_count),
      global_pattern_trusted_count =
        as.integer(global_pattern_evidence$trusted_peak_count),
      global_pattern_peak_x = global_pattern_evidence$peaks$x,
      global_pattern_ratios = global_pattern_evidence$peaks$ratio,
      global_pattern_slots = if (!is.null(trusted_peaks) &&
          nrow(trusted_peaks)) trusted_peaks$slot else integer(0),
      global_pattern_evidence = global_pattern_evidence,
      terminal_harmonic_slot =
        as.integer(global_pattern_evidence$terminal_slot),
      single_peak_elbow = single_peak_elbow,
      single_peak_elbow_applied = isTRUE(single_peak_elbow_applied),
      single_peak_elbow_guard = single_peak_elbow_guard,
      range_selected_ploidy = range_selected_n,
      range_terminal_anchor_slot = range_terminal_anchor_slot,
      range_selection_basis = range_selection_basis,
      range_selected_size_model = range_selected_size_model,
      selected_parameters = list(
        minimum_pattern = minimum_pattern,
        fit_harmonic_count = fitting_ploidy,
        single_peak_elbow_applied = isTRUE(single_peak_elbow_applied)
      ),
      selection_only = TRUE
    )))
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
  # The minimum observable pattern remains the biological interpretation.
  # Tail-supported harmonics only extend the numerical curve representation;
  # they do not change the reported minimum ploidy.
  for (extra_mode in "never") {
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
  best_extra <- "never"

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
  best_extra <- "never"

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
  best_extra <- "never"

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
                symmetry_value = best_symmetry, extra_mode = "never")

  # Compare parsimonious common-width and power-trend models with the adaptive
  # soft model at the best nuisance parameters. This is especially important
  # for Chinese sturgeon: the early visible peaks can support a trend while the
  # weak late peaks should not each spend an unconstrained width parameter.
  current_best <- best_row(n)
  best_ratio <- current_best$main_peak_ratio[1]
  best_weight <- current_best$main_peak_weight[1]
  best_symmetry <- current_best$symmetry_tolerance[1]
  best_extra <- "never"
  for (sigma_mode_value in c("common", "trend", "adaptive")) {
    run_candidate(n, ratio = best_ratio, weight = best_weight,
                  symmetry_value = best_symmetry, extra_mode = best_extra,
                  sigma_mode_value = sigma_mode_value)
  }

  # Tune the K-component model completely first.  When at least one harmonic
  # is unresolved, validate K+1 from that fitted model with only a small local
  # search.  This avoids repeating the complete parameter grid for K+1 while
  # still allowing the final components to redistribute the tail signal.
  primary_seed_fit <- NULL
  nuisance_seed_table <- all_rows()
  nuisance_seed_table <- nuisance_seed_table[
      nuisance_seed_table$n == n &
      nuisance_seed_table$stage == "screening" &
      !nuisance_seed_table$extra_component &
      nuisance_seed_table$observed_peak_count <
        nuisance_seed_table$fit_harmonic_count &
      is.finite(nuisance_seed_table$score), , drop = FALSE
  ]
  nuisance_seed_table <- nuisance_seed_table[
    order(nuisance_seed_table$score), , drop = FALSE
  ]
  nuisance_seed_table <- head(nuisance_seed_table, 1L)
  if (nrow(nuisance_seed_table)) {
    seed_row <- nuisance_seed_table[1, , drop = FALSE]
    seed_key <- candidate_key(
      seed_row$n[1], seed_row$main_peak_ratio[1],
      seed_row$main_peak_weight[1], seed_row$symmetry_tolerance[1],
      "never", seed_row$sigma_mode[1], "screening",
      seed_row$fit_harmonic_count[1]
    )
    if (exists(seed_key, envir = cache, inherits = FALSE)) {
      primary_seed_fit <- get(seed_key, envir = cache)$fit
    }

    # Inherit the selected K model and test only its exact nuisance settings
    # plus immediate weight/peak-window neighbours.  At most five distinct
    # screening fits replace the former full coarse and width-model grids.
    local_weights <- sort(unique(pmin(1.40, pmax(
      0.60, seed_row$main_peak_weight[1] + c(-0.02, 0, 0.02)
    ))))
    for (weight in local_weights) {
      run_candidate(
        n, ratio = seed_row$main_peak_ratio[1], weight = weight,
        symmetry_value = seed_row$symmetry_tolerance[1],
        extra_mode = "legacy",
        sigma_mode_value = seed_row$sigma_mode[1],
        initial_fit = primary_seed_fit
      )
    }
    if (!is.na(seed_row$main_peak_ratio[1])) {
      local_ratios <- sort(unique(pmin(0.65, pmax(
        0.15, seed_row$main_peak_ratio[1] + c(-0.02, 0.02)
      ))))
      for (ratio in local_ratios) {
        run_candidate(
          n, ratio = ratio, weight = seed_row$main_peak_weight[1],
          symmetry_value = seed_row$symmetry_tolerance[1],
          extra_mode = "legacy",
          sigma_mode_value = seed_row$sigma_mode[1],
          initial_fit = primary_seed_fit
        )
      }
    }
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
  # Always verify the best model from both numerical component counts.  A fast
  # screening approximation must not decide whether the nuisance component is
  # retained.
  for (has_extra in c(FALSE, TRUE)) {
    mode_rows <- screening_table[
      screening_table$extra_component == has_extra, , drop = FALSE
    ]
    if (nrow(mode_rows)) finalists <- rbind(finalists, mode_rows[1, , drop = FALSE])
  }

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
      extra_mode = if (isTRUE(finalists$extra_component[i])) {
        "legacy"
      } else "never",
      sigma_mode_value = finalists$sigma_mode[i],
      stage = "verification",
      initial_fit = if (isTRUE(finalists$extra_component[i])) {
        primary_seed_fit
      } else NULL
    )
  }

  search_table <- all_rows()
  verification_table <- search_table[
    search_table$n == n & search_table$stage == "verification" &
      is.finite(search_table$score), , drop = FALSE
  ]
  verification_table <- verification_table[
    order(verification_table$score), , drop = FALSE
  ]
  selected <- verification_table[1, , drop = FALSE]
  selected_initial_fit <- NULL
  no_extra_rows <- verification_table[
    !verification_table$extra_component, , drop = FALSE
  ]
  extra_rows <- verification_table[
    verification_table$extra_component, , drop = FALSE
  ]
  # Restore the original deterministic rule: whenever fewer reliable peaks
  # are observed than the required fitting harmonics, K+1 is part of the
  # model rather than a candidate that may be removed by a global score. The
  # preceding K fit is retained only as the warm start for the local K+1 fit.
  required_extra_selected <- nrow(extra_rows) > 0L
  if (required_extra_selected) {
    selected <- extra_rows[1, , drop = FALSE]
    selected_initial_fit <- primary_seed_fit
    cat("Reliable peak count is below the fitting harmonic count;",
        "retaining one final tail component\n")
  } else if (nrow(no_extra_rows)) {
    selected <- no_extra_rows[1, , drop = FALSE]
  }

  # A minimum pattern can be biologically correct yet use too few numerical
  # components when several late dosage peaks overlap into one broad tail.
  # Try the separately detected strong-tail extent only for a materially poor
  # current fit and only when it adds at least two primary harmonics.  The
  # alternative is accepted on measured curve quality, with explicit guards
  # for the main peak and its left wing.  No sample name or known ploidy is
  # involved, and the reported minimum pattern n remains unchanged.
  strong_tail_slot <- suppressWarnings(as.integer(
    numerical_tail_fit_slot
  ))
  base_fit_slot <- as.integer(selected$fit_harmonic_count[1])
  try_strong_tail <- is.finite(strong_tail_slot) &&
    strong_tail_slot >= base_fit_slot + 2L &&
    is.finite(selected$nrmse[1]) && selected$nrmse[1] > 0.03
  if (try_strong_tail) {
    cat("Current curve leaves a continuous strong tail unresolved;",
        "testing", strong_tail_slot, "primary fitting harmonics\n")
    extended_weights <- sort(unique(pmin(1.40, pmax(
      0.60, selected$main_peak_weight[1] + c(-0.02, 0, 0.02)
    ))))
    extended_screen <- lapply(extended_weights, function(weight) {
      run_candidate(
        n, ratio = selected$main_peak_ratio[1], weight = weight,
        symmetry_value = selected$symmetry_tolerance[1],
        extra_mode = "legacy", sigma_mode_value = selected$sigma_mode[1],
        stage = "screening", initial_fit = primary_seed_fit,
        fit_count_override = strong_tail_slot
      )
    })
    extended_screen <- extended_screen[vapply(
      extended_screen,
      function(value) is.finite(value$row$fit_score[1]), logical(1)
    )]
    if (length(extended_screen)) {
      extended_scores <- vapply(
        extended_screen, function(value) value$row$score[1], numeric(1)
      )
      extended_seed <- extended_screen[[which.min(extended_scores)]]
      extended_row <- extended_seed$row
      extended_verified <- run_candidate(
        n, ratio = extended_row$main_peak_ratio[1],
        weight = extended_row$main_peak_weight[1],
        symmetry_value = extended_row$symmetry_tolerance[1],
        extra_mode = "legacy",
        sigma_mode_value = extended_row$sigma_mode[1],
        stage = "verification", initial_fit = extended_seed$fit,
        fit_count_override = strong_tail_slot
      )
      alternative <- extended_verified$row
      global_gain <- 1 - alternative$nrmse[1] /
        max(selected$nrmse[1], 1e-8)
      main_guard <- alternative$main_peak_nrmse[1] <=
        selected$main_peak_nrmse[1] * 1.05 + 0.001
      left_guard <- alternative$left_overshoot_error[1] <=
        selected$left_overshoot_error[1] * 1.10 + 0.001
      accept_strong_tail <- is.finite(global_gain) && global_gain >= 0.20 &&
        isTRUE(main_guard) && isTRUE(left_guard) &&
        alternative$fit_score[1] < selected$fit_score[1]
      cat("  strong-tail NRMSE gain =",
          sprintf("%.1f%%", 100 * global_gain), "\n")
      if (accept_strong_tail) {
        cat("  accepted the guarded strong-tail extent\n")
        cat("  restarting the complete", strong_tail_slot,
            "x search so fitting and plotting use one unified pathway\n")
        unified_fit <- .GenomeGMM_auto(
          data_file = data_file,
          n = strong_tail_slot,
          n_candidates = n_candidates,
          search = search,
          symmetry_tolerance = symmetry_tolerance,
          position_mode = position_mode,
          joint_refinement = joint_refinement,
          species_name = species_name,
          seed = seed,
          save_search = save_search,
          candidate_max_iter = candidate_max_iter,
          candidate_tolerance = candidate_tolerance,
          verification_top = verification_top,
          verification_max_iter = verification_max_iter,
          final_max_iter = final_max_iter,
          save_plot = save_plot,
          .preloaded_data = .preloaded_data
        )
        unified_fit$minimum_pattern <- minimum_pattern
        unified_fit$minimum_observable_pattern <- minimum_pattern
        unified_fit$minimum_pattern_candidates <- pattern_selection_table
        unified_fit$selected_parameters$minimum_pattern <- minimum_pattern
        return(invisible(unified_fit))
      } else {
        cat("  rejected the guarded strong-tail alternative because",
            "its protected fit metrics did not all improve\n")
      }
    }
  }

  # Disabled rollback record for the former guarded K/K+1 comparison. It is
  # deliberately skipped now that the original deterministic rule is active.
  if (!required_extra_selected && nrow(no_extra_rows) && nrow(extra_rows)) {
    no_extra_best <- no_extra_rows[1, , drop = FALSE]
    extra_best <- extra_rows[1, , drop = FALSE]
    refit_for_component_comparison <- function(row, extra_mode,
                                               initial_fit = NULL) {
      comparison_args <- list(
        data_file = data_file, n = row$n[1],
        symmetry_tolerance = row$symmetry_tolerance[1],
        main_peak_weight = row$main_peak_weight[1],
        use_calibration = FALSE, species_name = species_name,
        sigma_mode = row$sigma_mode[1], position_mode = position_mode,
        joint_refinement = joint_refinement,
        extra_component_mode = extra_mode,
        .preloaded_data = cached_full_data,
        .analysis_max_x = auto_analysis_limit,
        .global_pattern_evidence = global_pattern_evidence,
        .fit_harmonic_count = row$fit_harmonic_count[1],
        .initial_fit = initial_fit,
        save_plot = FALSE, em_max_iter = final_max_iter
      )
      if (isTRUE(row$calibration[1])) {
        comparison_args$main_peak_ratio <- row$main_peak_ratio[1]
      }
      set.seed(seed)
      invisible(utils::capture.output(
        comparison_fit <- do.call(.GenomeGMM_fit, comparison_args),
        type = "output"
      ))
      comparison_fit
    }
    no_extra_fit <- refit_for_component_comparison(no_extra_best, "never")
    extra_fit <- refit_for_component_comparison(
      extra_best, "legacy", initial_fit = no_extra_fit
    )
    no_extra_health <- assess_component_health(no_extra_fit, n)
    extra_health <- assess_component_health(extra_fit, n)

    # Diagnose whether the final primary component is swallowing unresolved
    # tail signal. Ratios are relative to preceding components, so this rule
    # is independent of species, absolute coverage and ploidy.
    tail_structure <- function(fit) {
      component_order <- order(fit$miu)
      primary_count <- min(as.integer(fit$fit_harmonic_count),
                           length(component_order))
      if (primary_count < 2L) {
        return(list(width_inflation = 1, weight_inflation = 1,
                    last_share = 1, inflated = FALSE))
      }
      primary <- component_order[seq_len(primary_count)]
      widths <- pmax(fit$sigma[primary], 1e-8)
      shares <- pmax(fit$alpha[primary], 0) /
        max(sum(pmax(fit$alpha, 0)), 1e-12)
      previous_count <- min(3L, primary_count - 1L)
      previous_indices <- seq.int(primary_count - previous_count,
                                  primary_count - 1L)
      expected_width_ratio <- if (primary_count >= 3L) {
        historical_ratios <- widths[2:(primary_count - 1L)] /
          widths[1:(primary_count - 2L)]
        stats::median(tail(historical_ratios, 3L), na.rm = TRUE)
      } else 1
      if (!is.finite(expected_width_ratio)) expected_width_ratio <- 1
      expected_width <- widths[primary_count - 1L] *
        pmin(1.25, pmax(0.90, expected_width_ratio))
      width_inflation <- widths[primary_count] /
        max(expected_width, 1e-8)
      previous_share <- stats::median(shares[previous_indices], na.rm = TRUE)
      weight_inflation <- shares[primary_count] /
        max(previous_share, 1e-8)
      last_share <- shares[primary_count]
      inflated <- is.finite(width_inflation) &&
          is.finite(weight_inflation) && is.finite(last_share) &&
        (width_inflation >= 1.12 || weight_inflation >= 2.50 ||
           last_share >= 0.30)
      list(width_inflation = width_inflation,
           weight_inflation = weight_inflation,
           last_share = last_share, inflated = inflated)
    }
    no_extra_tail <- tail_structure(no_extra_fit)
    extra_tail <- tail_structure(extra_fit)
    relative_fit_gain <- (no_extra_fit$fit_score - extra_fit$fit_score) /
      max(no_extra_fit$fit_score, 1e-8)
    collapse_rescue <-
      no_extra_health$collapsed_biological_components > 0L &&
      extra_health$collapsed_biological_components == 0L &&
      extra_fit$fit_score <= no_extra_fit$fit_score * 1.05
    main_guard <- extra_fit$fit_metrics$main_peak_nrmse <=
      no_extra_fit$fit_metrics$main_peak_nrmse * 1.05 + 0.001
    strict_left_guard <- extra_fit$fit_metrics$left_overshoot_error <=
        max(0.0025,
            no_extra_fit$fit_metrics$left_overshoot_error * 1.03 + 0.0003) &&
      extra_fit$fit_metrics$left_relative_error <=
        max(0.035,
            no_extra_fit$fit_metrics$left_relative_error * 1.03 + 0.002)
    # Do not let a tiny boundary crossing in one left-wing metric veto a
    # clearly healthier n+1 decomposition.  This is a general Pareto rule:
    # the extra component must substantially improve the whole curve, reduce
    # actual left overshoot, keep the main peak accurate, and limit the small
    # trade-off in relative left-wing error.  It applies identically whether
    # n was supplied by the user or inferred from explicit harmonic peaks.
    left_pareto_guard <- relative_fit_gain >= 0.20 &&
      extra_fit$fit_metrics$left_overshoot_error <=
        no_extra_fit$fit_metrics$left_overshoot_error &&
      extra_fit$fit_metrics$main_peak_nrmse <= 0.025 &&
      extra_fit$fit_metrics$left_relative_error <=
        no_extra_fit$fit_metrics$left_relative_error +
          max(0.0025,
              0.08 * no_extra_fit$fit_metrics$left_relative_error)
    left_guard <- strict_left_guard || left_pareto_guard
    # A small absolute main-peak error may trade slightly against a major
    # whole-curve and left-wing improvement.  This admits the well-resolved
    # virginalis decomposition while still rejecting a visibly damaged main
    # peak.
    main_guard <- main_guard ||
      (extra_fit$fit_metrics$main_peak_nrmse <= 0.025 &&
       relative_fit_gain >= 0.20 &&
       extra_fit$fit_metrics$left_overshoot_error <=
         no_extra_fit$fit_metrics$left_overshoot_error)

    # Permit a small global-score trade-off when K alone makes its final
    # primary component anomalously broad/heavy and K+1 materially restores
    # that component without damaging the main peak or left wing.
    tail_structure_improved <-
      extra_tail$width_inflation <= no_extra_tail$width_inflation * 0.90 ||
      extra_tail$weight_inflation <= no_extra_tail$weight_inflation * 0.80 ||
      extra_tail$last_share <= no_extra_tail$last_share * 0.85
    tail_inflation_rescue <- isTRUE(no_extra_tail$inflated) &&
      isTRUE(tail_structure_improved) &&
      extra_fit$fit_score <= no_extra_fit$fit_score * 1.05

    standard_extra_accept <- main_guard && left_guard &&
      (relative_fit_gain >= 0.04 || collapse_rescue) &&
      extra_fit$fit_score < no_extra_fit$fit_score
    # Tail rescue uses absolute error protection. A small relative increase
    # from an already accurate main/left fit must not force the last primary
    # component to absorb the whole tail. Actual left overshoot must improve.
    tail_rescue_main_guard <-
      extra_fit$fit_metrics$main_peak_nrmse <= 0.020
    tail_rescue_left_guard <-
      extra_fit$fit_metrics$left_overshoot_error <=
        no_extra_fit$fit_metrics$left_overshoot_error &&
      extra_fit$fit_metrics$left_relative_error <=
        no_extra_fit$fit_metrics$left_relative_error + 0.005
    tail_rescue_accept <- tail_inflation_rescue &&
      tail_rescue_main_guard && tail_rescue_left_guard
    accept_extra <- isTRUE(extra_health$extra_guard_passed) &&
      (standard_extra_accept || tail_rescue_accept)
    if (accept_extra && !strict_left_guard && left_pareto_guard) {
      cat("  K+1 accepted by high-gain left-wing Pareto guard\n")
    }
    if (accept_extra && tail_inflation_rescue) {
      cat("  K+1 accepted by final-component inflation rescue\n")
    }
    selected <- if (accept_extra) extra_best else no_extra_best
    if (accept_extra) selected_initial_fit <- no_extra_fit
    cat("Guarded numerical component comparison:",
        if (accept_extra) "accepted one final tail component" else
          "kept tail-supported primary components", "\n")
    cat("  joint fit score K / K+1 =",
        round(no_extra_fit$fit_score, 6), "/",
        round(extra_fit$fit_score, 6), "\n")
    cat("  joint left overshoot K / K+1 =",
        round(no_extra_fit$fit_metrics$left_overshoot_error, 6), "/",
        round(extra_fit$fit_metrics$left_overshoot_error, 6), "\n")
    cat("  final primary width inflation K / K+1 =",
        round(no_extra_tail$width_inflation, 3), "/",
        round(extra_tail$width_inflation, 3), "\n")
    cat("  final primary weight inflation K / K+1 =",
        round(no_extra_tail$weight_inflation, 3), "/",
        round(extra_tail$weight_inflation, 3), "\n")
    cat("  final primary share K / K+1 =",
        round(no_extra_tail$last_share, 3), "/",
        round(extra_tail$last_share, 3), "\n")
  }
  search_table <- all_rows()
  search_table$selected <- search_table$trial == selected$trial[1]

  cat("\nBest automatic parameters:\n")
  cat("  n =", selected$n[1], "\n")
  cat("  tail-supported primary fitting harmonics =",
      selected$fit_harmonic_count[1], "\n")
  cat("  calibration =", selected$calibration[1], "\n")
  if (selected$calibration[1]) {
    cat("  main_peak_ratio =", selected$main_peak_ratio[1], "\n")
  }
  cat("  main_peak_weight =", selected$main_peak_weight[1], "\n")
  cat("  symmetry_tolerance =", selected$symmetry_tolerance[1], "\n")
  cat("  optional final tail component =",
      selected$extra_component[1], "\n")
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
    extra_component_mode = if (isTRUE(selected$extra_component[1])) {
      "legacy"
    } else "never",
    .preloaded_data = cached_full_data,
    .analysis_max_x = auto_analysis_limit,
    .global_pattern_evidence = global_pattern_evidence,
    .fit_harmonic_count = selected$fit_harmonic_count[1],
    .initial_fit = if (isTRUE(selected$extra_component[1])) {
      selected_initial_fit
    } else NULL,
    .first_peak_override = if (shallow_anchor_active) {
      global_pattern_evidence$shallow_anchor_x
    } else NULL,
    .start_index_override = if (shallow_anchor_active) {
      global_pattern_evidence$start_index
    } else NULL,
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
    search_file <- paste0(
      "GenomeGMM_search_", file_basename, "_estimated_ploidy_",
      selected$fit_harmonic_count[1], "x_", timestamp, ".csv"
    )
    utils::write.csv(search_table, search_file, row.names = FALSE, na = "")
    cat("Parameter search table saved as:", search_file, "\n")
  }

  final_fit$selected_parameters <- list(
    n = selected$n[1],
    minimum_pattern = minimum_pattern,
    fit_harmonic_count = selected$fit_harmonic_count[1],
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
    selection_score = selected$score[1],
    single_peak_elbow_applied = isTRUE(single_peak_elbow_applied)
  )
  final_fit$search_table <- search_table
  final_fit$search_file <- search_file
  final_fit$minimum_pattern <- minimum_pattern
  final_fit$minimum_observable_pattern <- minimum_pattern
  final_fit$minimum_pattern_candidates <- pattern_selection_table
  final_fit$single_peak_elbow <- single_peak_elbow
  final_fit$single_peak_elbow_applied <-
    isTRUE(single_peak_elbow_applied)
  final_fit$single_peak_elbow_guard <- single_peak_elbow_guard
  final_fit$range_selected_ploidy <- range_selected_n
  final_fit$range_terminal_anchor_slot <- range_terminal_anchor_slot
  final_fit$range_selection_basis <- range_selection_basis
  final_fit$range_selected_size_model <- range_selected_size_model
  invisible(final_fit)
}

.GenomeGMM_pattern_composition <- function(pattern, multiplier) {
  pattern <- as.integer(pattern)
  multiplier <- as.integer(multiplier)
  if (pattern == 1L) {
    return(strrep("A", multiplier))
  }
  if (pattern == 2L) {
    return(paste0(strrep("A", multiplier), strrep("a", multiplier)))
  }
  if (pattern <= length(letters)) {
    alleles <- letters[seq_len(pattern)]
    return(paste0(vapply(alleles, strrep, character(1),
                         times = multiplier), collapse = ""))
  }
  paste0("(", paste0("allele", seq_len(pattern), "^", multiplier,
                      collapse = " "), ")")
}

.GenomeGMM_fit_quality <- function(fit) {
  nrmse <- suppressWarnings(as.numeric(fit$fit_metrics$nrmse[1]))
  if (!is.finite(nrmse)) return("Unknown")
  if (isTRUE(fit$converged) && nrmse <= 0.03) return("High")
  if (isTRUE(fit$converged) && nrmse <= 0.07) return("Moderate")
  "Low"
}

.GenomeGMM_pattern_assessment <- function(fit) {
  candidates <- fit$minimum_pattern_candidates
  pattern <- if (!is.null(fit$minimum_pattern) &&
      length(fit$minimum_pattern) == 1L &&
      is.finite(fit$minimum_pattern)) {
    as.integer(fit$minimum_pattern)
  } else {
    as.integer(fit$n)
  }
  observable_pattern <- if (!is.null(fit$minimum_observable_pattern) &&
      length(fit$minimum_observable_pattern) == 1L &&
      is.finite(fit$minimum_observable_pattern)) {
    as.integer(fit$minimum_observable_pattern)
  } else {
    pattern
  }
  answer <- list(
    confidence = "Not assessed", relative_gap = NA_real_,
    alternative_pattern = NA_integer_, table = candidates,
    global_explicit_peak_count = if (!is.null(fit$global_pattern_peak_count)) {
      as.integer(fit$global_pattern_peak_count)
    } else NA_integer_,
    global_trusted_peak_count = if (!is.null(
        fit$global_pattern_trusted_count
      )) {
      as.integer(fit$global_pattern_trusted_count)
    } else NA_integer_,
    global_peak_positions = fit$global_pattern_peak_x,
    global_peak_ratios = fit$global_pattern_ratios,
    global_trusted_slots = fit$global_pattern_slots
  )
  if (is.null(candidates) || !nrow(candidates) ||
      !all(c("observed_peak_count", "terminal_harmonic_slot") %in%
           names(candidates))) return(answer)
  usable <- candidates[is.finite(candidates$terminal_harmonic_slot),
                       , drop = FALSE]
  if (!nrow(usable)) return(answer)
  usable <- usable[order(usable$n), , drop = FALSE]
  selected_index <- match(pattern, usable$n)
  selected_count <- if (!is.na(selected_index)) {
    as.integer(usable$observed_peak_count[selected_index])
  } else NA_integer_
  next_count <- if (!is.na(selected_index) && selected_index < nrow(usable)) {
    as.integer(usable$observed_peak_count[selected_index + 1L])
  } else NA_integer_
  counts <- c(selected_count, next_count)
  counts <- counts[is.finite(counts)]
  selected_terminal <- if (!is.na(selected_index)) {
    as.integer(usable$terminal_harmonic_slot[selected_index])
  } else NA_integer_
  floor_applied <- observable_pattern < 2L
  answer$confidence <- if (floor_applied) {
    "Low"
  } else if (observable_pattern == pattern) {
    "High"
  } else if (observable_pattern < pattern) {
    "Moderate"
  } else {
    "Low"
  }
  alternative_counts <- sort(unique(counts[
    counts >= 2L & counts != pattern
  ]))
  if (length(alternative_counts)) {
    answer$alternative_pattern <- max(2L, alternative_counts[1])
  }
  answer$detection_counts <- counts
  answer$selected_window_peak_count <- selected_count
  answer$next_window_peak_count <- next_count
  answer$terminal_harmonic_slot <- selected_terminal
  answer$minimum_observable_pattern <- observable_pattern
  answer$minimum_inferred_ploidy <- pattern
  answer$minimum_two_floor_applied <- floor_applied
  answer
}

.GenomeGMM_write_minimum_pattern_report <- function(
    data_file, species_name, fit, max_multiplier = 10L,
    monoploid_range = NULL) {
  single_peak_ambiguity <- isTRUE(fit$single_peak_ambiguity)
  default_pattern <- if (!is.null(fit$minimum_pattern) &&
      length(fit$minimum_pattern) == 1L &&
      is.finite(fit$minimum_pattern)) {
    as.integer(fit$minimum_pattern)
  } else {
    as.integer(fit$n)
  }
  fitted_ploidy <- if (!is.null(fit$fit_harmonic_count) &&
      length(fit$fit_harmonic_count) == 1L &&
      is.finite(fit$fit_harmonic_count)) {
    as.integer(fit$fit_harmonic_count)
  } else {
    default_pattern
  }
  total_fitted_components <- max(
    fitted_ploidy,
    if (!is.null(fit$miu)) length(fit$miu) else fitted_ploidy
  )

  # Report-only ploidy uncertainty: a clear or broad later maximum can still
  # support the nearest integer harmonic even when it misses the stricter
  # fitting tolerance.  This does not refit any curve.  Always retain every
  # integer candidate through the total number of fitted components: whenever
  # an auxiliary n + 1 component was used, both n and n + 1 must be reported,
  # irrespective of how many explicit peaks were detected.
  peak_supported_pattern <- fitted_ploidy
  peak_supported_pattern_upper <- fitted_ploidy
  evidence_peaks <- NULL
  if (!is.null(fit$global_pattern_evidence) &&
      !is.null(fit$global_pattern_evidence$peaks)) {
    evidence_peaks <- fit$global_pattern_evidence$peaks
  }
  if (!is.null(evidence_peaks) && nrow(evidence_peaks) &&
      "ratio" %in% names(evidence_peaks)) {
    ratios <- suppressWarnings(as.numeric(evidence_peaks$ratio))
    eligible <- is.finite(ratios) & ratios >= 1.5 &
      ratios <= total_fitted_components + 0.65
    if (any(eligible)) {
      terminal_ratio <- max(ratios[eligible])
      ratio_fraction <- terminal_ratio - floor(terminal_ratio)
      # A broad high-ploidy maximum centred near the midpoint between two
      # integer harmonics does not justify choosing only one side.  Preserve
      # both adjacent integer interpretations in the report; this changes no
      # fitted curve and lets prior g-hap knowledge resolve the ambiguity.
      bracket_two_slots <- ratio_fraction >= 0.35 &&
        ratio_fraction <= 0.65
      if (isTRUE(bracket_two_slots)) {
        supported_lower <- floor(terminal_ratio)
        supported_upper <- ceiling(terminal_ratio)
      } else {
        supported_lower <- supported_upper <- round(terminal_ratio)
      }
      supported_lower <- max(
        2L, min(total_fitted_components, as.integer(supported_lower))
      )
      supported_upper <- max(
        supported_lower,
        min(total_fitted_components, as.integer(supported_upper))
      )
      peak_supported_pattern <- supported_lower
      peak_supported_pattern_upper <- supported_upper
    }
  }
  # A validated next-harmonic shoulder is biological peak evidence even though
  # it is not a strict local maximum.  Keep the report consistent with the
  # fitted curve: when that shoulder supplied the fitted terminal harmonic,
  # do not fall back to the smaller count of explicit maxima.
  shoulder_supported_pattern <- if (!is.null(
      fit$global_pattern_evidence$harmonic_shoulder_slot
    )) {
    suppressWarnings(as.integer(
      fit$global_pattern_evidence$harmonic_shoulder_slot
    ))
  } else NA_integer_
  if (length(shoulder_supported_pattern) == 1L &&
      is.finite(shoulder_supported_pattern) &&
      shoulder_supported_pattern >= 2L &&
      shoulder_supported_pattern <= fitted_ploidy) {
    peak_supported_pattern <- max(
      peak_supported_pattern, shoulder_supported_pattern
    )
    peak_supported_pattern_upper <- max(
      peak_supported_pattern_upper, shoulder_supported_pattern
    )
  }
  candidate_lower <- min(peak_supported_pattern, fitted_ploidy)
  report_tail_pattern <- fitted_ploidy
  if (!isTRUE(fit$single_peak_elbow_applied) &&
      !is.null(fit$global_pattern_evidence)) {
    report_evidence <- fit$global_pattern_evidence
    strong_tail_slot <- suppressWarnings(as.integer(
      report_evidence$strong_tail_fit_slot
    ))
    if (length(strong_tail_slot) == 1L && is.finite(strong_tail_slot)) {
      report_tail_pattern <- max(report_tail_pattern, strong_tail_slot)
    }
    weak_tail_slot <- suppressWarnings(as.integer(
      report_evidence$weak_harmonic_fit_slot
    ))
    if (isTRUE(report_evidence$weak_harmonic_chain_supported) &&
        length(weak_tail_slot) == 1L && is.finite(weak_tail_slot)) {
      report_tail_pattern <- max(report_tail_pattern, weak_tail_slot)
    }
    support_table <- report_evidence$harmonic_window_support
    if (!is.null(support_table) && nrow(support_table) &&
        all(c("slot", "median_fraction", "maximum_fraction") %in%
            names(support_table)) && length(strong_tail_slot) == 1L &&
        is.finite(strong_tail_slot)) {
      weak_slot <- strong_tail_slot + 1L
      weak_index <- match(weak_slot, support_table$slot)
      following_index <- match(weak_slot + 1L, support_table$slot)
      if (!is.na(weak_index) && !is.na(following_index)) {
        weak_median <- support_table$median_fraction[weak_index]
        weak_maximum <- support_table$maximum_fraction[weak_index]
        following_median <- support_table$median_fraction[following_index]
        final_supported_window <-
          is.finite(weak_median) && weak_median >= 0.06 &&
          is.finite(weak_maximum) && weak_maximum >= 0.085 &&
          is.finite(following_median) &&
          following_median <= max(0.05, 0.65 * weak_median)
        # A complete 1:m series of explicit, trusted peaks can be followed by
        # one weak unresolved shoulder.  Preserve m + 1 as a report-only
        # candidate when that shoulder has continuous window-level support and
        # the next harmonic window then decays.  Requiring the complete prefix
        # keeps an isolated tail fluctuation from extending an incomplete
        # pattern.  This changes neither the selected curve nor its PDF.
        evidence_minimum <- suppressWarnings(as.integer(
          report_evidence$minimum_pattern
        ))
        evidence_explicit <- suppressWarnings(as.integer(
          report_evidence$explicit_peak_count
        ))
        evidence_trusted <- suppressWarnings(as.integer(
          report_evidence$trusted_peak_count
        ))
        complete_harmonic_prefix <-
          length(evidence_minimum) == 1L &&
          is.finite(evidence_minimum) &&
          evidence_minimum == fitted_ploidy &&
          strong_tail_slot == fitted_ploidy &&
          length(evidence_explicit) == 1L &&
          is.finite(evidence_explicit) &&
          evidence_explicit >= fitted_ploidy &&
          length(evidence_trusted) == 1L &&
          is.finite(evidence_trusted) &&
          evidence_trusted >= fitted_ploidy
        terminal_shoulder_window <-
          isTRUE(complete_harmonic_prefix) &&
          is.finite(weak_median) && weak_median >= 0.025 &&
          is.finite(weak_maximum) && weak_maximum >= 0.10 &&
          is.finite(following_median) &&
          following_median <= max(0.02, 0.65 * weak_median)
        if (isTRUE(final_supported_window) ||
            isTRUE(terminal_shoulder_window)) {
          report_tail_pattern <- max(report_tail_pattern, weak_slot)
        }
      }
    }
  }
  candidate_upper <- max(
    peak_supported_pattern_upper, fitted_ploidy, total_fitted_components,
    report_tail_pattern
  )
  ploidy_candidates <- seq.int(candidate_lower, candidate_upper)
  pattern <- candidate_lower
  if (isTRUE(fit$hidden_prefix_recovery_accepted)) {
    plausible <- suppressWarnings(as.integer(
      fit$hidden_prefix_plausible_slots
    ))
    plausible <- plausible[is.finite(plausible) & plausible >= 1L]
    # The fitted candidate remains the general-rule base pattern. Nearby
    # candidates within 20% of its protected score are shown in the matrix so
    # a visually overlapping 4x/5x case is not reported as falsely unique.
    pattern <- fitted_ploidy
    ploidy_candidates <- sort(unique(c(
      plausible, fitted_ploidy,
      seq.int(fitted_ploidy, candidate_upper)
    )))
  }
  if (isTRUE(single_peak_ambiguity)) {
    # One genuinely resolved genomic peak is compatible with A, AA, AAA, ... .
    # These are biological multiples of the same observed peak, not different
    # harmonic base candidates with different homozygous depths.
    pattern <- 1L
    ploidy_candidates <- 1L
  }
  range_selected_pattern <- suppressWarnings(as.integer(
    fit$range_selected_ploidy
  ))
  if (!is.null(monoploid_range) &&
      length(range_selected_pattern) == 1L &&
      is.finite(range_selected_pattern) && range_selected_pattern >= 1L) {
    # Once an expected-size range has selected the base curve, every reported
    # quantity must use the same harmonic: peak depth, Monoploid size,
    # Holoploid size and the P = m * k ambiguity rule.  Keeping the detector's
    # pre-selection lower bound here would mix (for example) a 10-component
    # PDF with a ninth-peak genome-size calculation.
    pattern <- range_selected_pattern
    candidate_upper <- max(candidate_upper, pattern)
    ploidy_candidates <- seq.int(pattern, candidate_upper)
  }

  max_multiplier <- max(2L, as.integer(max_multiplier))
  multipliers <- seq_len(max_multiplier)
  compatible_ploidies <- pattern * multipliers
  fitted_g_hap <- if (!is.null(fit$g_hap_bp) &&
      length(fit$g_hap_bp) == 1L &&
      is.finite(fit$g_hap_bp) && fit$g_hap_bp > 0) {
    as.numeric(fit$g_hap_bp)
  } else {
    as.numeric(fit$genome_size_bp) / fitted_ploidy
  }
  ordered_miu <- sort(as.numeric(fit$miu[is.finite(fit$miu) & fit$miu > 0]))
  if (!isTRUE(single_peak_ambiguity) &&
      length(ordered_miu) < max(ploidy_candidates)) {
    if (length(ordered_miu) < 2L) {
      stop("Cannot calculate a ploidy candidate without its peak depth.")
    }
    increments <- diff(ordered_miu)
    increments <- increments[is.finite(increments) & increments > 0]
    if (!length(increments)) {
      stop("Cannot extrapolate a candidate peak depth.")
    }
    harmonic_step <- stats::median(
      tail(increments, min(5L, length(increments)))
    )
    while (length(ordered_miu) < max(ploidy_candidates)) {
      ordered_miu <- c(ordered_miu, tail(ordered_miu, 1L) + harmonic_step)
    }
  }
  candidate_depths <- if (isTRUE(single_peak_ambiguity)) {
    rep(ordered_miu[1L], length(ploidy_candidates))
  } else {
    ordered_miu[ploidy_candidates]
  }
  total_kmers <- if (!is.null(fit$total_kmers) &&
      length(fit$total_kmers) == 1L && is.finite(fit$total_kmers) &&
      fit$total_kmers > 0) {
    as.numeric(fit$total_kmers)
  } else if (fitted_ploidy <= length(ordered_miu)) {
    fitted_g_hap * ordered_miu[fitted_ploidy]
  } else {
    NA_real_
  }
  if (!is.finite(total_kmers) || total_kmers <= 0) {
    stop("Cannot calculate candidate genome sizes without total k-mer counts.")
  }
  candidate_g_hap <- total_kmers / candidate_depths
  pattern_row <- match(pattern, ploidy_candidates)
  if (is.na(pattern_row)) pattern_row <- 1L
  g_hap <- candidate_g_hap[pattern_row]
  ploidy_candidate_matrix <- data.frame(
    ploidy = ploidy_candidates,
    homozygous_peak_depth = candidate_depths,
    full_gse_bp = ploidy_candidates * candidate_g_hap,
    haploid_gse_bp = candidate_g_hap,
    doubled_ploidy = 2L * ploidy_candidates,
    doubled_full_gse_bp = 2L * ploidy_candidates * candidate_g_hap,
    doubled_haploid_gse_bp = candidate_g_hap,
    stringsAsFactors = FALSE
  )
  matrix_header <- sprintf(
    "%-12s %-23s %-23s %-8s %-23s %-23s",
    "Base-Ploidy", "Holoploid genome size", "Monoploid genome size",
    "Ploidy", "Holoploid genome size", "Monoploid genome size"
  )
  matrix_rows <- vapply(seq_len(nrow(ploidy_candidate_matrix)), function(i) {
    row <- ploidy_candidate_matrix[i, , drop = FALSE]
    sprintf(
      "%-12s %-23s %-23s %-8s %-23s %-23s",
      paste0(row$ploidy, "x"),
      sprintf("%.2f Mb", row$full_gse_bp / 1e6),
      sprintf("%.2f Mb", row$haploid_gse_bp / 1e6),
      paste0(row$doubled_ploidy, "x"),
      sprintf("%.2f Mb", row$doubled_full_gse_bp / 1e6),
      sprintf("%.2f Mb", row$doubled_haploid_gse_bp / 1e6)
    )
  }, character(1))
  # Candidate interpretations used only for filtering by an optional prior
  # Monoploid genome-size range. The fitted unit is Total k-mers / mu[m].
  # In a direct/autopolyploid interpretation that unit is one basic set. In a
  # balanced allopolyploid interpretation one complete distinct-subgenome
  # complement contains m such units (e.g. ABD for m = 3). These rows do not
  # change or rerun the fitted GMM.
  range_candidates <- do.call(rbind, lapply(
    seq_len(nrow(ploidy_candidate_matrix)), function(i) {
      row <- ploidy_candidate_matrix[i, , drop = FALSE]
      direct_base <- data.frame(
        base_ploidy = as.integer(row$ploidy),
        interpretation = "Direct/autopolyploid base pattern",
        ploidy = as.integer(row$ploidy),
        holoploid_bp = as.numeric(row$full_gse_bp),
        monoploid_bp = as.numeric(row$haploid_gse_bp),
        stringsAsFactors = FALSE
      )
      direct_duplicated <- data.frame(
        base_ploidy = as.integer(row$ploidy),
        interpretation = "Direct/autopolyploid duplication",
        ploidy = as.integer(row$doubled_ploidy),
        holoploid_bp = as.numeric(row$doubled_full_gse_bp),
        monoploid_bp = as.numeric(row$doubled_haploid_gse_bp),
        stringsAsFactors = FALSE
      )
      if (isTRUE(single_peak_ambiguity)) {
        return(rbind(direct_base, direct_duplicated))
      }
      allopolyploid_monoploid <-
        as.numeric(row$ploidy) * as.numeric(row$haploid_gse_bp)
      allopolyploid_base <- data.frame(
        base_ploidy = as.integer(row$ploidy),
        interpretation = "Balanced-allopolyploid base pattern",
        ploidy = as.integer(row$ploidy),
        holoploid_bp = as.numeric(row$full_gse_bp),
        monoploid_bp = allopolyploid_monoploid,
        stringsAsFactors = FALSE
      )
      allopolyploid_duplicated <- data.frame(
        base_ploidy = as.integer(row$ploidy),
        interpretation = "Balanced-allopolyploid duplication",
        ploidy = as.integer(row$doubled_ploidy),
        holoploid_bp = as.numeric(row$doubled_full_gse_bp),
        monoploid_bp = allopolyploid_monoploid,
        stringsAsFactors = FALSE
      )
      rbind(
        direct_base, direct_duplicated,
        allopolyploid_base, allopolyploid_duplicated
      )
    }
  ))
  range_is_supplied <- !is.null(monoploid_range)
  range_text <- if (range_is_supplied) {
    sprintf("[%.2f, %.2f] Mb", monoploid_range[1L], monoploid_range[2L])
  } else {
    "[a, b] Mb"
  }
  range_matches <- if (range_is_supplied) {
    lower_bp <- monoploid_range[1L] * 1e6
    upper_bp <- monoploid_range[2L] * 1e6
    range_candidates[
      range_candidates$monoploid_bp >= lower_bp &
        range_candidates$monoploid_bp <= upper_bp,
      , drop = FALSE
    ]
  } else {
    range_candidates[0, , drop = FALSE]
  }
  range_prediction_header <- sprintf(
    "%-12s %-35s %-8s %-23s %-23s",
    "Base-Ploidy", "Interpretation", "Ploidy",
    "Holoploid genome size", "Monoploid genome size"
  )
  range_prediction_rows <- if (nrow(range_matches)) {
    vapply(seq_len(nrow(range_matches)), function(i) {
      row <- range_matches[i, , drop = FALSE]
      sprintf(
        "%-12s %-35s %-8s %-23s %-23s",
        paste0(row$base_ploidy, "x"), row$interpretation,
        paste0(row$ploidy, "x"),
        sprintf("%.2f Mb", row$holoploid_bp / 1e6),
        sprintf("%.2f Mb", row$monoploid_bp / 1e6)
      )
    }, character(1))
  } else character(0)
  fixed_holoploid_candidates <- range_candidates[
    range_candidates$interpretation == "Fixed-Holoploid reinterpretation",
    , drop = FALSE
  ]
  fixed_holoploid_header <- sprintf(
    "%-12s %-8s %-23s %-23s",
    "Base-Ploidy", "Ploidy", "Holoploid genome size",
    "Monoploid genome size"
  )
  fixed_holoploid_rows <- if (nrow(fixed_holoploid_candidates)) {
    vapply(seq_len(nrow(fixed_holoploid_candidates)), function(i) {
      row <- fixed_holoploid_candidates[i, , drop = FALSE]
      sprintf(
        "%-12s %-8s %-23s %-23s",
        paste0(row$base_ploidy, "x"), paste0(row$ploidy, "x"),
        sprintf("%.2f Mb", row$holoploid_bp / 1e6),
        sprintf("%.2f Mb", row$monoploid_bp / 1e6)
      )
    }, character(1))
  } else character(0)
  selected_size_model <- if (range_is_supplied &&
      !is.null(fit$range_selected_size_model) &&
      length(fit$range_selected_size_model) == 1L &&
      !is.na(fit$range_selected_size_model)) {
    as.character(fit$range_selected_size_model)
  } else {
    "direct_or_autopolyploid"
  }
  per_pattern_unit_size <- g_hap
  reported_monoploid_size <- if (identical(
      selected_size_model, "balanced_allopolyploid"
    )) {
    pattern * per_pattern_unit_size
  } else {
    per_pattern_unit_size
  }
  minimum_size <- pattern * per_pattern_unit_size
  # Within the minimum-pattern family, g-hap stays fixed and Full GSE scales
  # linearly. A different base candidate uses its own fitted homozygous peak.
  compatible_sizes <- compatible_ploidies * per_pattern_unit_size
  compatible_g_hap_sizes <- rep(
    reported_monoploid_size, length(compatible_ploidies)
  )
  compositions <- vapply(
    multipliers, function(k) .GenomeGMM_pattern_composition(pattern, k),
    character(1)
  )
  compatibility <- data.frame(
    multiplier = multipliers,
    ploidy = compatible_ploidies,
    composition = compositions,
    genome_size_bp = compatible_sizes,
    g_hap_size_bp = compatible_g_hap_sizes,
    stringsAsFactors = FALSE
  )
  assessment_fit <- fit
  assessment_fit$minimum_pattern <- pattern
  pattern_assessment <- .GenomeGMM_pattern_assessment(assessment_fit)
  if (isTRUE(fit$hidden_prefix_recovery_accepted)) {
    # The selected curve is supported by its fit to the original histogram,
    # but some of its prefix peaks were reconstructed rather than independently
    # observed. It must therefore never be labelled High confidence.
    pattern_assessment$confidence <- "Moderate"
    pattern_assessment$alternative_pattern <- {
      alternatives <- setdiff(
        suppressWarnings(as.integer(fit$hidden_prefix_plausible_slots)),
        fitted_ploidy
      )
      alternatives <- alternatives[is.finite(alternatives)]
      if (length(alternatives)) alternatives[1] else NA_integer_
    }
  }
  if (!isTRUE(fit$hidden_prefix_recovery_accepted) &&
      !is.null(fit$global_pattern_evidence) &&
      isTRUE(
        fit$global_pattern_evidence$weak_harmonic_chain_supported
      ) && fitted_ploidy > pattern &&
      identical(pattern_assessment$confidence, "High")) {
    # Explicit peaks support the minimum pattern strongly, but the longer
    # numerical curve is selected from unresolved harmonic windows.  Report
    # moderate confidence so "High" is not mistaken for a unique absolute
    # ploidy assignment.
    pattern_assessment$confidence <- "Moderate"
  }
  fit_quality <- .GenomeGMM_fit_quality(fit)
  explicit_peak_count <- if (is.finite(
      pattern_assessment$global_trusted_peak_count
    )) {
    as.integer(pattern_assessment$global_trusted_peak_count)
  } else if (!is.null(pattern_assessment$detection_counts) &&
      length(pattern_assessment$detection_counts)) {
    as.integer(round(stats::median(pattern_assessment$detection_counts)))
  } else NA_integer_

  file_stem <- tools::file_path_sans_ext(basename(data_file))
  result_file <- if (!is.null(fit$pdf_file) && length(fit$pdf_file) == 1L &&
      !is.na(fit$pdf_file) && nzchar(fit$pdf_file)) {
    paste0(tools::file_path_sans_ext(fit$pdf_file), "_result.txt")
  } else {
    paste0(
      "em_fitted_gmm_", file_stem, "_estimated_ploidy_",
      fitted_ploidy, "x_", format(Sys.time(), "%Y%m%d_%H%M%S"),
      "_result.txt"
    )
  }

  fitted_candidate_row <- match(
    fitted_ploidy, ploidy_candidate_matrix$ploidy
  )
  if (is.na(fitted_candidate_row)) fitted_candidate_row <- 1L
  fitted_curve_g_hap <-
    ploidy_candidate_matrix$haploid_gse_bp[fitted_candidate_row]
  fitted_curve_full_gse <-
    ploidy_candidate_matrix$full_gse_bp[fitted_candidate_row]

  range_selection_lines <- if (!range_is_supplied) {
    c(
      "Ploidy selection:",
      paste0(
        "Compare the Monoploid genome size of each candidate with the ",
        "expected range [a, b]. Choose a candidate ploidy P whose ",
        "Monoploid genome size falls within this range, then rerun ",
        "GenomeGMM with n = P for the final fit: GenomeGMM(\"",
        basename(data_file), "\", n = P). If more than one candidate falls ",
        "within [a, b], the histogram alone does not identify a unique ",
        "ploidy."
      )
    )
  } else if (!nrow(range_matches)) {
    c(
      paste0("Prediction using expected Monoploid genome-size range ",
             range_text, ":"),
      paste0(
        "No candidate interpretation has a Monoploid genome size within ",
        range_text, ". Please verify the expected range or use additional ",
        "biological evidence; GenomeGMM does not force a ploidy prediction."
      )
    )
  } else {
    matched_base_count <- length(unique(range_matches$base_ploidy))
    conclusion <- if (matched_base_count == 1L && nrow(range_matches) == 1L) {
      paste0(
        "One candidate interpretation falls within ", range_text,
        ". Its potential ploidy is ", range_matches$ploidy[1L], "x."
      )
    } else {
      paste0(
        matched_base_count, " base-ploidy candidate",
        if (matched_base_count == 1L) "" else "s", " (",
        nrow(range_matches), " biological interpretation",
        if (nrow(range_matches) == 1L) "" else "s", ") fall within ",
        range_text, ". The histogram and the supplied range do not identify ",
        "a unique biological ploidy."
      )
    }
    c(
      paste0("Prediction using expected Monoploid genome-size range ",
             range_text, ":"),
      range_prediction_header,
      range_prediction_rows,
      conclusion
    )
  }

  report_body <- if (range_is_supplied) {
    c(
      "Inference of ploidy level and related genome size estimation:",
      "",
      sprintf(
        "Selected subgenome size: %.2f Mb",
        per_pattern_unit_size / 1e6
      ),
      paste0("Selected size interpretation: ", selected_size_model),
      "",
      range_selection_lines
    )
  } else {
    c(
      "Inference of ploidy level and related genome size estimation:",
      "",
      matrix_header,
      matrix_rows,
      "",
      sprintf(
        paste0(
          "Note, when using a base ploidy of %d for curve fitting, a ",
          "Monoploid genome size of %.2f Mb and a Holoploid genome size of ",
          "%.2f Mb were achieved."
        ),
        fitted_ploidy, fitted_curve_g_hap / 1e6,
        fitted_curve_full_gse / 1e6
      ),
      "",
      if (isTRUE(single_peak_ambiguity)) {
        paste0(
          "This single-peak histogram is compatible with haploid A, ",
          "homozygous diploid AA, homozygous triploid AAA, and higher ",
          "homozygous multiples. All possibilities share the same ",
          "Monoploid genome size. Therefore, the histogram and an expected ",
          "Monoploid genome-size range do not identify a unique ploidy; an ",
          "externally established ploidy P is required before rerunning ",
          "GenomeGMM(\"", basename(data_file), "\", n = P)."
        )
      } else {
        range_selection_lines
      }
    )
  }

  general_rule_lines <- if (identical(
      selected_size_model, "balanced_allopolyploid"
    )) {
    c(
      sprintf(
        "The potential ploidy P = %d * k, where k could be 1, 2, ...",
        pattern
      ),
      sprintf(
        "Subgenome size = %.2f Mb.",
        per_pattern_unit_size / 1e6
      ),
      sprintf(
        "Given balanced-allopolyploid Monoploid genome size = %.2f Mb (= %.2f Mb * %d),",
        reported_monoploid_size / 1e6,
        per_pattern_unit_size / 1e6, pattern
      ),
      sprintf(
        paste0(
          "the Holoploid genome size = Monoploid genome size * k, ",
          "i.e., %.2f Mb, %.2f Mb, ..."
        ),
        reported_monoploid_size / 1e6,
        2 * reported_monoploid_size / 1e6
      )
    )
  } else {
    c(
      sprintf(
        "The potential ploidy P = %d * k, where k could be 1, 2, ...",
        pattern
      ),
      sprintf(
        "Given direct/autopolyploid Monoploid genome size = %.2f Mb,",
        reported_monoploid_size / 1e6
      ),
      sprintf(
        paste0(
          "the Holoploid genome size = Monoploid genome size * P = ",
          "%.1f Mb * %d * k, i.e., %.2f Mb, %.2f Mb, ..."
        ),
        reported_monoploid_size / 1e6, pattern,
        minimum_size / 1e6, 2 * minimum_size / 1e6
      )
    )
  }

  lines <- c(
    sprintf(
      "GenomeGMM ploidy and genome-size inference for sample: %s",
      file_stem
    ),
    paste0("Expected range of monoploid genome size: ", range_text),
    "",
    "Definitions:",
    paste0(
      "Monoploid genome size: the DNA content of one basic chromosome set ",
      "(1Cx), e.g., A, AB, ABC, etc."
    ),
    paste0(
      "Holoploid genome size: the total DNA content of the complete ",
      "chromosome complement at the corresponding ploidy level, e.g., ",
      "{AA, AAA, ...}, {AABB, AAABBB, ...}, ",
      "{AABBCC, AAABBBCCC, ...}, etc."
    ),
    "",
    report_body,
    "",
    if (isTRUE(single_peak_ambiguity)) {
      paste0(
        "Moreover, as sequencing reads of homozygous genomes like \"A\", ",
        "\"AA\", \"AAA\", ... would show the same single-peak histogram ",
        "of k-mer counts, their ploidy cannot be distinguished with the ",
        "histogram alone. Thus, a more general rule for Holoploid genome-size ",
        "inference:"
      )
    } else {
      paste0(
        "Moreover, as sequencing reads of genomes like \"AB\", \"AABB\", ",
        "\"AAABBB\", ... would show the same histogram of k-mer counts ",
        "(which also applies for other ploidy levels). The ploidy of them ",
        "cannot be distinguished with the histogram of k-mer counts. Thus, a ",
        "more general rule for Holoploid genome-size inference:"
      )
    },
    "",
    general_rule_lines
  )
  writeLines(enc2utf8(lines), con = result_file, useBytes = TRUE)

  list(
    file = result_file,
    minimum_pattern = pattern,
    peak_supported_pattern = peak_supported_pattern,
    peak_supported_pattern_upper = peak_supported_pattern_upper,
    ploidy_candidates = ploidy_candidates,
    ploidy_candidate_matrix = ploidy_candidate_matrix,
    compatible_ploidies = compatible_ploidies,
    compatible_genome_sizes_bp = compatible_sizes,
    compatible_g_hap_sizes_bp = compatible_g_hap_sizes,
    monoploid_range_mb = monoploid_range,
    range_matched_candidates = range_matches,
    range_match_count = nrow(range_matches),
    compositions = compositions,
    g_hap_bp = reported_monoploid_size,
    minimum_pattern_genome_size_bp = minimum_size,
    homozygous_peak_depth = candidate_depths[pattern_row],
    pattern_confidence = pattern_assessment$confidence,
    pattern_assessment = pattern_assessment,
    fit_quality = fit_quality
  )
}

.GenomeGMM_get_legacy_environment <- function() {
  .GenomeGMM_embedded_legacy_environment
}

.GenomeGMM_rename_legacy_pdf <- function(pdf_file, data_file,
                                         fitted_ploidy) {
  if (is.null(pdf_file) || length(pdf_file) != 1L || is.na(pdf_file) ||
      !nzchar(pdf_file)) return(pdf_file)
  old_path <- pdf_file
  if (!file.exists(old_path)) return(pdf_file)
  old_name <- basename(old_path)
  timestamp_match <- regexec(
    "_n[0-9]+_([0-9]{8}_[0-9]{6})[.]pdf$", old_name
  )
  timestamp_parts <- regmatches(old_name, timestamp_match)[[1]]
  timestamp <- if (length(timestamp_parts) >= 2L) {
    timestamp_parts[2]
  } else {
    format(Sys.time(), "%Y%m%d_%H%M%S")
  }
  file_stem <- tools::file_path_sans_ext(basename(data_file))
  new_name <- paste0(
    "em_fitted_gmm_", file_stem, "_estimated_ploidy_",
    as.integer(fitted_ploidy), "x_", timestamp, ".pdf"
  )
  new_path <- file.path(dirname(old_path), new_name)
  if (!identical(normalizePath(old_path, mustWork = FALSE),
                 normalizePath(new_path, mustWork = FALSE))) {
    if (file.exists(new_path)) {
      stop("Cannot rename fitted PDF because the target already exists: ",
           new_path)
    }
    if (!file.rename(old_path, new_path)) {
      stop("Could not rename fitted PDF: ", old_path)
    }
  }
  if (identical(dirname(pdf_file), ".")) new_name else new_path
}

# A single local maximum is not sufficient by itself: a high-ploidy spectrum
# can contain overlapping shoulders.  Treat the automatic result as the
# A/AA/AAA... ambiguity only when the global detector found exactly one clear
# and trusted harmonic, no later numerical harmonic, and no structured tail
# beyond the immediately following window.  Known-n calls never use this gate.
.GenomeGMM_is_genuine_single_peak <- function(pattern_evidence) {
  if (is.null(pattern_evidence)) return(FALSE)
  scalar_integer <- function(name, fallback = NA_integer_) {
    value <- suppressWarnings(as.integer(pattern_evidence[[name]]))
    if (length(value) == 1L && is.finite(value)) value else fallback
  }
  scalar_integer("explicit_peak_count") == 1L &&
    scalar_integer("trusted_peak_count") == 1L &&
    scalar_integer("terminal_slot") == 1L &&
    scalar_integer("numerical_terminal_slot") == 1L &&
    scalar_integer("supported_fit_slot", 2L) <= 2L &&
    scalar_integer("strong_tail_fit_slot", 2L) <= 2L
}

# One explicit local maximum can still be the second member of two strongly
# overlapping harmonic components. Before applying the A/AA/AAA single-peak
# interpretation, compare lightweight one- and two-component fits on the same
# histogram. The two-component interpretation is retained only when it gives
# a large original-curve improvement and the fitted components themselves
# satisfy harmonic-position, weight and coarse width guards. The full fitter
# still performs its stricter optimisation and biological constraints later.
.GenomeGMM_assess_hidden_two_peak_shoulder <- function(
    data_file, pattern_evidence) {
  answer <- list(
    supported = FALSE, nrmse_one = NA_real_, nrmse_two = NA_real_,
    relative_improvement = NA_real_, first_miu = NA_real_,
    second_miu = NA_real_, harmonic_ratio = NA_real_,
    minimum_share = NA_real_, maximum_width_jump = NA_real_
  )
  if (is.null(pattern_evidence) ||
      !is.finite(pattern_evidence$first_peak_x)) return(answer)
  legacy <- .GenomeGMM_get_legacy_environment()
  fit_candidate <- function(component_count) {
    fit <- NULL
    tryCatch(
      invisible(utils::capture.output(
        fit <- legacy$.GenomeGMM_fit(
          data_file = data_file, n = component_count,
          main_peak_weight = 1, symmetry_tolerance = 0.2,
          max_sigma_ratio = 1.35, sigma_mode = "adaptive",
          position_mode = "adaptive", joint_refinement = "off",
          extra_component_mode = "never", save_plot = FALSE,
          em_max_iter = 500, em_tolerance = 1e-5
        ), type = "output"
      )), error = function(error) NULL
    )
    fit
  }
  one <- fit_candidate(1L)
  two <- fit_candidate(2L)
  required <- c("alpha", "miu", "sigma", "fit_metrics")
  if (is.null(one) || is.null(two) ||
      !all(required %in% names(one)) ||
      !all(required %in% names(two))) return(answer)
  nrmse_one <- suppressWarnings(as.numeric(one$fit_metrics$nrmse[1]))
  nrmse_two <- suppressWarnings(as.numeric(two$fit_metrics$nrmse[1]))
  order_two <- order(two$miu)
  if (length(order_two) < 2L || !is.finite(nrmse_one) ||
      !is.finite(nrmse_two) || nrmse_one <= 0) return(answer)
  biological <- order_two[1:2]
  centres <- as.numeric(two$miu[biological])
  widths <- as.numeric(two$sigma[biological])
  shares <- pmax(as.numeric(two$alpha[biological]), 0) /
    max(sum(pmax(as.numeric(two$alpha), 0)), 1e-12)
  visible_centre <- as.numeric(pattern_evidence$first_peak_x)
  harmonic_ratio <- centres[2] / centres[1]
  maximum_width_jump <- max(widths) / max(min(widths), 1e-12)
  relative_improvement <- 1 - nrmse_two / nrmse_one
  minimum_share <- min(shares)
  # These are screening guards, not final parameter constraints. In
  # particular, the final automatic fit continues to enforce the established
  # adjacent-width protection and rejects a deteriorating two-peak solution.
  supported <- all(is.finite(c(
    centres, widths, shares, harmonic_ratio, maximum_width_jump,
    relative_improvement, visible_centre
  ))) &&
    relative_improvement >= 0.25 &&
    minimum_share >= 0.08 &&
    harmonic_ratio >= 1.65 && harmonic_ratio <= 2.35 &&
    centres[1] / visible_centre >= 0.35 &&
    centres[1] / visible_centre <= 0.65 &&
    centres[2] / visible_centre >= 0.80 &&
    centres[2] / visible_centre <= 1.20 &&
    maximum_width_jump <= 1.85
  answer$supported <- isTRUE(supported)
  answer$nrmse_one <- nrmse_one
  answer$nrmse_two <- nrmse_two
  answer$relative_improvement <- relative_improvement
  answer$first_miu <- centres[1]
  answer$second_miu <- centres[2]
  answer$harmonic_ratio <- harmonic_ratio
  answer$minimum_share <- minimum_share
  answer$maximum_width_jump <- maximum_width_jump
  answer
}

# When the ordinary prominence filter finds no explicit peak, retain a broad
# low-prominence maximum only if it is supported by several consecutive rising
# and falling bins.  This is deliberately a fallback: clear peaks continue to
# use detect_global_pattern_peaks(), while count-poor tail ripples are removed
# by the post-error height floor.
.GenomeGMM_broad_anchor_candidates <- function(full_data) {
  empty <- list(valid = FALSE, selected_x = NA_real_,
                candidates = data.frame())
  if (is.null(full_data) || ncol(full_data) < 2L || nrow(full_data) < 20L) {
    return(empty)
  }
  data <- data.frame(
    x = suppressWarnings(as.numeric(full_data[, 1])),
    y = suppressWarnings(as.numeric(full_data[, 2]))
  )
  data <- data[
    is.finite(data$x) & is.finite(data$y) & data$x > 0 & data$y >= 0,
    , drop = FALSE
  ]
  if (nrow(data) < 20L) return(empty)
  data <- stats::aggregate(y ~ x, data = data, FUN = sum)
  data <- data[order(data$x), , drop = FALSE]
  x <- data$x
  y <- data$y
  bin_steps <- diff(x)
  bin_steps <- bin_steps[is.finite(bin_steps) & bin_steps > 0]
  if (!length(bin_steps)) return(empty)
  smooth <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  smooth[!is.finite(smooth)] <- y[!is.finite(smooth)]
  maxima <- which(diff(sign(diff(smooth))) < 0) + 1L
  maxima <- maxima[maxima >= 8L & maxima <= length(x) - 7L]
  if (!length(maxima)) return(empty)
  reference_start <- min(5L, length(y))
  reference_height <- max(y[reference_start:length(y)], na.rm = TRUE)
  if (!is.finite(reference_height) || reference_height <= 0) return(empty)
  rows <- lapply(maxima, function(index) {
    before <- diff(smooth[(index - 6L):index])
    after <- diff(smooth[index:(index + 6L)])
    rising_support <- sum(before > 0, na.rm = TRUE)
    falling_support <- sum(after < 0, na.rm = TRUE)
    trough_indices <- max(reference_start, index - 20L):index
    trough_height <- min(smooth[trough_indices], na.rm = TRUE)
    prominence <- (smooth[index] - trough_height) / max(smooth[index], 1)
    height_fraction <- smooth[index] / reference_height
    if (!is.finite(prominence) || !is.finite(height_fraction) ||
        rising_support < 4L || falling_support < 4L ||
        prominence < 0.012 || height_fraction < 0.01) return(NULL)
    data.frame(
      x = x[index], height = smooth[index],
      prominence = prominence, height_fraction = height_fraction,
      rising_support = rising_support, falling_support = falling_support,
      support_score = prominence * sqrt(height_fraction) *
        (rising_support + falling_support) / 12,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  candidates <- do.call(rbind, rows)
  candidates <- candidates[order(
    candidates$support_score, decreasing = TRUE
  ), , drop = FALSE]
  rownames(candidates) <- NULL
  list(valid = TRUE, selected_x = candidates$x[1],
       candidates = candidates)
}

# Find where the initial sequencing-error wall has changed into a sustained,
# slower decline.  This supplies several nearby integration starts instead of
# fixing a coverage such as 5 for every histogram.
.GenomeGMM_prefix_mass_window <- function(
    full_data, first_visible_peak_x, terminal_visible_peak_x = NULL) {
  empty <- list(valid = FALSE, start_index = 1L, start_x = NA_real_,
                starts = numeric(0), end_x = NA_real_)
  if (is.null(full_data) || nrow(full_data) < 12L ||
      !is.finite(first_visible_peak_x) ||
      first_visible_peak_x <= 0) return(empty)
  if (is.null(terminal_visible_peak_x) ||
      !is.finite(terminal_visible_peak_x) ||
      terminal_visible_peak_x < first_visible_peak_x) {
    terminal_visible_peak_x <- first_visible_peak_x
  }
  x <- suppressWarnings(as.numeric(full_data[, 1]))
  y <- suppressWarnings(as.numeric(full_data[, 2]))
  keep <- is.finite(x) & is.finite(y) & x > 0 & y >= 0
  x <- x[keep]
  y <- y[keep]
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  prefix <- which(x[-length(x)] < 0.5 * first_visible_peak_x)
  if (length(prefix) < 6L) return(empty)
  log_drop <- -diff(log(pmax(y, 1)))
  robust_drop <- log_drop[prefix]
  robust_drop <- robust_drop[is.finite(robust_drop)]
  if (length(robust_drop) < 6L) return(empty)
  quartiles <- stats::quantile(
    robust_drop, c(0.25, 0.75), na.rm = TRUE, names = FALSE
  )
  drop_limit <- quartiles[2] + 1.5 * max(quartiles[2] - quartiles[1], 1e-6)
  slow <- is.finite(log_drop) & log_drop <= drop_limit
  runs <- rle(slow & seq_along(slow) %in% prefix)
  run_end <- cumsum(runs$lengths)
  run_start <- run_end - runs$lengths + 1L
  accepted <- which(runs$values & runs$lengths >= 3L)
  if (!length(accepted)) return(empty)
  start_index <- run_start[accepted[1]]
  # A broad hidden-prefix spectrum can show a second, much gentler error-decay
  # bend after the initial wall.  Starting at the first slow bin would still
  # let the first Gaussian absorb that bend.  Move beyond a sustained local
  # maximum of the smoothed log-drop rate when it occurs well before the clear
  # terminal peak (ZD306-like shape); ordinary spectra without this structure
  # retain the first sustained start above.
  smooth_drop <- as.numeric(stats::filter(
    log_drop, rep(1 / 3, 3), sides = 2
  ))
  smooth_drop[!is.finite(smooth_drop)] <- log_drop[!is.finite(smooth_drop)]
  drop_maxima <- which(diff(sign(diff(smooth_drop))) < 0) + 1L
  bin_width <- stats::median(diff(unique(x)))
  bend_candidates <- drop_maxima[
    drop_maxima >= start_index + 3L &
      x[drop_maxima] <= 0.40 * first_visible_peak_x
  ]
  if (length(bend_candidates)) {
    bend <- bend_candidates[which.max(smooth_drop[bend_candidates])]
    start_index <- max(start_index, min(length(x), bend + 3L))
  }
  start_x <- x[start_index]
  latest_start <- min(0.25 * first_visible_peak_x,
                      start_x + 15 * bin_width)
  if (!is.finite(latest_start) || latest_start <= start_x) {
    starts <- start_x
  } else {
    targets <- seq(start_x, latest_start, length.out = 6L)
    starts <- unique(vapply(targets, function(target) {
      x[which.min(abs(x - target))]
    }, numeric(1)))
  }
  list(valid = TRUE, start_index = as.integer(start_index),
       start_x = start_x, starts = starts,
       end_x = min(max(x), 1.8 * terminal_visible_peak_x))
}

# Prefix recovery must score the inferred first peak itself.  If the ordinary
# error-wall window begins to its right, move only this candidate's fitting and
# comparison window left to roughly the last quarter of its rising flank.  The
# baseline and candidate are then evaluated on exactly the same bins.
.GenomeGMM_candidate_prefix_window <- function(
    window, full_data, inferred_first_miu) {
  if (!isTRUE(window$valid) || !is.finite(inferred_first_miu) ||
      inferred_first_miu <= 0) return(window)
  x <- suppressWarnings(as.numeric(full_data[, 1]))
  y <- suppressWarnings(as.numeric(full_data[, 2]))
  x <- sort(unique(x[is.finite(x) & is.finite(y) & x > 0 & y >= 0]))
  if (length(x) < 6L) return(window)
  earliest_supported <- x[min(5L, length(x))]
  target <- max(earliest_supported, 0.75 * inferred_first_miu)
  target <- x[which.min(abs(x - target))]
  if (!is.finite(window$start_x) || target >= window$start_x) return(window)
  upper_start <- min(window$start_x, inferred_first_miu)
  targets <- if (upper_start > target) {
    seq(target, upper_start, length.out = 6L)
  } else target
  starts <- unique(vapply(targets, function(value) {
    x[which.min(abs(x - value))]
  }, numeric(1)))
  adjusted <- window
  adjusted$start_x <- min(starts)
  adjusted$start_index <- which.min(abs(x - adjusted$start_x))
  adjusted$starts <- starts
  adjusted
}

# After a hidden first peak has been reconstructed, the visible anchor need not
# be the final biological harmonic: several later peaks can overlap into a
# smooth declining wing.  Measure complete harmonic windows and extend only
# through a consecutive, data-supported tail.  This proposes candidates; the
# fitted-curve guards in .GenomeGMM_try_hidden_prefix_recovery() still decide
# whether every additional component is retained.
.GenomeGMM_prefix_tail_support <- function(
    full_data, window, first_miu, base_slot, maximum_additional = 6L) {
  empty <- list(valid = FALSE, maximum_slot = as.integer(base_slot),
                table = data.frame())
  if (!isTRUE(window$valid) || !is.finite(first_miu) || first_miu <= 0 ||
      !is.finite(base_slot) || base_slot < 1L) return(empty)
  data <- data.frame(
    x = suppressWarnings(as.numeric(full_data[, 1])),
    y = suppressWarnings(as.numeric(full_data[, 2]))
  )
  data <- data[
    is.finite(data$x) & is.finite(data$y) & data$x > 0 & data$y >= 0,
    , drop = FALSE
  ]
  if (nrow(data) < 12L) return(empty)
  data <- stats::aggregate(y ~ x, data = data, FUN = sum)
  data <- data[order(data$x), , drop = FALSE]
  smooth <- as.numeric(stats::filter(data$y, rep(1 / 5, 5), sides = 2))
  smooth[!is.finite(smooth)] <- data$y[!is.finite(smooth)]
  reliable <- data$x >= window$start_x
  peak_scale <- max(smooth[reliable], na.rm = TRUE)
  if (!is.finite(peak_scale) || peak_scale <= 0) return(empty)
  first_slot <- as.integer(base_slot) + 1L
  last_slot <- as.integer(base_slot) + max(1L, as.integer(maximum_additional))
  rows <- lapply(first_slot:last_slot, function(slot) {
    indices <- which(
      data$x >= (slot - 0.5) * first_miu &
        data$x < (slot + 0.5) * first_miu
    )
    if (length(indices) < 3L) {
      return(data.frame(
        slot = slot, median_fraction = NA_real_,
        maximum_fraction = NA_real_, relative_to_previous = NA_real_,
        supported = FALSE
      ))
    }
    data.frame(
      slot = slot,
      median_fraction = stats::median(smooth[indices], na.rm = TRUE) /
        peak_scale,
      maximum_fraction = max(smooth[indices], na.rm = TRUE) /
        peak_scale,
      relative_to_previous = NA_real_, supported = FALSE
    )
  })
  support <- do.call(rbind, rows)
  base_indices <- which(
    data$x >= (base_slot - 0.5) * first_miu &
      data$x < (base_slot + 0.5) * first_miu
  )
  previous_fraction <- if (length(base_indices)) {
    stats::median(smooth[base_indices], na.rm = TRUE) / peak_scale
  } else NA_real_
  maximum_slot <- as.integer(base_slot)
  for (i in seq_len(nrow(support))) {
    current <- support$median_fraction[i]
    current_max <- support$maximum_fraction[i]
    relative <- current / max(previous_fraction, 1e-12)
    support$relative_to_previous[i] <- relative
    supported <- is.finite(current) && is.finite(current_max) &&
      is.finite(relative) && current >= 0.05 && current_max >= 0.07 &&
      relative >= 0.45
    support$supported[i] <- isTRUE(supported)
    if (!isTRUE(supported)) break
    maximum_slot <- as.integer(support$slot[i])
    previous_fraction <- current
  }
  list(valid = TRUE, maximum_slot = maximum_slot, table = support)
}

# Detect a hidden peak-shaped change inside a descending prefix.  A genuine
# structure must show a sustained relaxation of the negative log-slope and a
# subsequent return to a steeper decline.  The interval from the strongest
# positive curvature to the slope apex is treated as a possible peak-centre
# region; it is not added to the plotted curve as a separate background.
.GenomeGMM_prefix_structure_candidates <- function(
    full_data, window, terminal_visible_peak_x) {
  empty <- list(valid = FALSE, candidates = data.frame())
  if (is.null(full_data) || !isTRUE(window$valid) ||
      !is.finite(terminal_visible_peak_x)) return(empty)
  x <- suppressWarnings(as.numeric(full_data[, 1]))
  y <- suppressWarnings(as.numeric(full_data[, 2]))
  keep <- is.finite(x) & is.finite(y) & x > 0 & y >= 0
  x <- x[keep]
  y <- y[keep]
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  if (length(x) < 25L) return(empty)
  bin_width <- stats::median(diff(unique(x)))
  smooth_log <- as.numeric(stats::filter(
    log(pmax(y, 1)), rep(1 / 5, 5), sides = 2
  ))
  raw_log <- log(pmax(y, 1))
  smooth_log[!is.finite(smooth_log)] <-
    raw_log[!is.finite(smooth_log)]
  log_slope <- c(NA_real_, diff(smooth_log) / diff(x))
  log_slope <- as.numeric(stats::filter(
    log_slope, rep(1 / 3, 3), sides = 2
  ))
  curvature <- c(NA_real_, diff(log_slope) / diff(x))
  curvature <- as.numeric(stats::filter(
    curvature, rep(1 / 3, 3), sides = 2
  ))
  # Counts below roughly 5x coverage are dominated by sequencing-error
  # k-mers in the intended histograms.  Start looking for a sustained slope
  # transition at 5, independently of the later model-fitting start.  The old
  # use of window$start_x could start around 15--20 and silently discard the
  # actual first transition near 10.  Five is only the detection floor: a
  # candidate still needs multi-bin slope support, harmonic consistency and a
  # guarded full-curve improvement below.
  slope_search_start <- max(5, min(x[x > 0], na.rm = TRUE))
  search <- which(
    x >= slope_search_start &
      x <= 0.45 * terminal_visible_peak_x &
      is.finite(log_slope) & is.finite(curvature)
  )
  local_maximum <- if (length(search) >= 18L) {
    search[which(diff(sign(diff(log_slope[search]))) < 0) + 1L]
  } else integer(0)
  rows <- lapply(local_maximum, function(apex) {
    # Early transitions around coverage 9--12 have only a few reliable bins
    # between the 5x noise floor and the slope apex.  Use compact, separated
    # flanks so they can be measured without reducing the sustained-support
    # requirements.
    left <- seq.int(max(search[1], apex - 6L),
                    max(search[1], apex - 2L))
    right <- seq.int(min(search[length(search)], apex + 2L),
                     min(search[length(search)], apex + 6L))
    before <- seq.int(max(search[1], apex - 8L), apex)
    after <- seq.int(apex, min(search[length(search)], apex + 8L))
    # When the transition occurs near coverage 9, the biologically usable
    # prefix begins at 5 and only three left-flank bins are available.  Three
    # consecutive bins are sufficient here because the candidate must still
    # pass curvature, harmonic, recovered-mass and full-curve fit guards.
    if (length(left) < 3L || length(right) < 3L) return(NULL)
    gain_left <- log_slope[apex] - mean(log_slope[left], na.rm = TRUE)
    gain_right <- log_slope[apex] - mean(log_slope[right], na.rm = TRUE)
    positive_support <- sum(curvature[before] > 0, na.rm = TRUE)
    negative_support <- sum(curvature[after] < 0, na.rm = TRUE)
    if (!is.finite(gain_left) || !is.finite(gain_right) ||
        gain_left < 0.012 || gain_right < 0.004 ||
        positive_support < 4L || negative_support < 2L ||
        log_slope[apex] >= -0.003) return(NULL)
    bend_pool <- before[is.finite(curvature[before])]
    bend <- bend_pool[which.max(curvature[bend_pool])]
    lower <- min(x[bend], x[apex]) - bin_width
    upper <- max(x[bend], x[apex]) + bin_width
    data.frame(
      lower_miu = lower,
      upper_miu = upper,
      slope_apex_x = x[apex],
      curvature_peak_x = x[bend],
      gain_from_left = gain_left,
      gain_over_right = gain_right,
      support_score = gain_left / 0.012 + gain_right / 0.004,
      candidate_type = "sustained_slope_structure",
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]

  # A very early hidden first peak can survive only as a small raw up/down
  # turn inside the descending error background and disappear under smoothing.
  # Keep it as a candidate, never as proof: the complete harmonic fit, recovered
  # k-mer mass, flank error, component weights and width guards below decide
  # whether it is accepted.
  raw_reference_start <- min(5L, length(y))
  raw_reference_height <- max(y[raw_reference_start:length(y)], na.rm = TRUE)
  raw_search <- which(
    seq_along(x) >= raw_reference_start + 1L &
      seq_along(x) <= length(x) - 2L &
      x <= 0.45 * terminal_visible_peak_x
  )
  raw_turns <- lapply(raw_search, function(index) {
    rise_fraction <- (y[index] - y[index - 1L]) / max(y[index], 1)
    fall_fraction <- (y[index] - y[index + 1L]) / max(y[index], 1)
    height_fraction <- y[index] / max(raw_reference_height, 1)
    embedded_decline <- y[index - 2L] > y[index - 1L] &&
      y[index + 1L] > y[index + 2L]
    if (!isTRUE(embedded_decline) || !is.finite(rise_fraction) ||
        !is.finite(fall_fraction) || !is.finite(height_fraction) ||
        rise_fraction < 0.002 || fall_fraction < 0.01 ||
        height_fraction < 0.03) return(NULL)
    data.frame(
      lower_miu = x[index] - bin_width,
      upper_miu = x[index] + bin_width,
      slope_apex_x = x[index], curvature_peak_x = x[index],
      gain_from_left = rise_fraction,
      gain_over_right = fall_fraction,
      support_score = 1 + rise_fraction / 0.002 +
        fall_fraction / 0.01,
      candidate_type = "raw_micro_turn",
      stringsAsFactors = FALSE
    )
  })
  raw_turns <- raw_turns[!vapply(raw_turns, is.null, logical(1))]
  rows <- c(rows, raw_turns)
  if (!length(rows)) return(empty)
  candidates <- do.call(rbind, rows)
  candidates <- candidates[order(
    candidates$support_score, decreasing = TRUE
  ), , drop = FALSE]
  rownames(candidates) <- NULL
  list(valid = TRUE, candidates = candidates)
}

.GenomeGMM_hidden_prefix_metrics <- function(
    fit, full_data, window, first_visible_peak_x,
    terminal_visible_peak_x, biological_count) {
  empty <- list(valid = FALSE, score = Inf, nrmse = Inf,
                under_mass = Inf, over_mass = Inf,
                prefix_under_mass = Inf, prefix_over_mass = Inf,
                prefix_under_fraction = Inf,
                prefix_over_fraction = Inf,
                prefix_under_run = 0L,
                first_right_over_mass = Inf,
                first_right_over_nrmse = Inf,
                main_peak_nrmse = Inf, collapsed = Inf,
                minimum_component_share = 0,
                maximum_width_jump = Inf)
  required <- c("alpha", "miu", "sigma", ".diagnostic_total_samples_fit",
                ".diagnostic_bin_width")
  if (is.null(fit) || !all(required %in% names(fit)) ||
      !isTRUE(window$valid)) return(empty)
  x <- suppressWarnings(as.numeric(full_data[, 1]))
  y <- suppressWarnings(as.numeric(full_data[, 2]))
  keep <- is.finite(x) & is.finite(y) & y >= 0 &
    x >= min(window$starts) & x <= window$end_x
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 10L) return(empty)
  component_counts <- vapply(seq_along(fit$alpha), function(j) {
    fit$alpha[j] * stats::dnorm(x, fit$miu[j], fit$sigma[j]) *
      fit$.diagnostic_total_samples_fit * fit$.diagnostic_bin_width
  }, numeric(length(x)))
  prediction <- if (is.null(dim(component_counts))) {
    component_counts
  } else rowSums(component_counts)
  prediction[!is.finite(prediction)] <- 0
  scale <- max(y, na.rm = TRUE)
  weights <- 0.25 + 0.75 * sqrt(pmax(y, 0) / max(scale, 1))
  nrmse <- sqrt(sum(weights * ((prediction - y) / max(scale, 1))^2) /
                  sum(weights))
  under_mass <- over_mass <- numeric(length(window$starts))
  for (i in seq_along(window$starts)) {
    use <- x >= window$starts[i]
    grey_mass <- max(sum(x[use] * y[use]), 1)
    under_mass[i] <- sum(x[use] * pmax(y[use] - prediction[use], 0)) /
      grey_mass
    over_mass[i] <- sum(x[use] * pmax(prediction[use] - y[use], 0)) /
      grey_mass
  }
  reliable_mass <- max(sum(x * y), 1)
  prefix_use <- x <= first_visible_peak_x
  prefix_grey_mass <- max(sum(x[prefix_use] * y[prefix_use]), 1)
  prefix_under_count <- sum(
    x[prefix_use] * pmax(y[prefix_use] - prediction[prefix_use], 0)
  )
  prefix_over_count <- sum(
    x[prefix_use] * pmax(prediction[prefix_use] - y[prefix_use], 0)
  )
  prefix_under_mass <- prefix_under_count / reliable_mass
  prefix_over_mass <- prefix_over_count / reliable_mass
  prefix_under_fraction <- prefix_under_count / prefix_grey_mass
  prefix_over_fraction <- prefix_over_count / prefix_grey_mass
  local_scale <- pmax(y, 0.02 * scale, 1)
  sustained_underfit <- prefix_use & y >= 0.02 * scale &
    (y - prediction) / local_scale >= 0.10
  underfit_runs <- rle(sustained_underfit)
  prefix_under_run <- if (any(underfit_runs$values)) {
    max(underfit_runs$lengths[underfit_runs$values])
  } else 0L
  component_order <- order(fit$miu)
  biological_count <- min(as.integer(biological_count),
                          length(component_order))
  biological <- component_order[seq_len(biological_count)]
  shares <- pmax(fit$alpha, 0) / max(sum(pmax(fit$alpha, 0)), 1e-12)
  # A hidden-prefix recovery requires positive numerical support for every
  # added harmonic.  The ordinary high-ploidy search uses a gentler floor;
  # here a component sitting at the EM floor is evidence against recovery.
  collapsed <- sum(!is.finite(shares[biological]) |
                     shares[biological] < 0.008)
  minimum_component_share <- suppressWarnings(min(
    shares[biological], na.rm = TRUE
  ))
  if (!is.finite(minimum_component_share)) minimum_component_share <- 0
  widths <- fit$sigma[biological][order(fit$miu[biological])]
  maximum_width_jump <- if (length(widths) > 1L) {
    max(pmax(widths[-1L] / widths[-length(widths)],
             widths[-length(widths)] / widths[-1L]))
  } else 1
  main_peak_nrmse <- suppressWarnings(as.numeric(
    fit$fit_metrics$main_peak_nrmse
  ))
  if (!is.finite(main_peak_nrmse)) main_peak_nrmse <- Inf
  ordered_biological <- biological[order(fit$miu[biological])]
  first_centre <- fit$miu[ordered_biological[1]]
  first_right_end <- if (length(ordered_biological) >= 2L) {
    fit$miu[ordered_biological[2]]
  } else {
    first_centre + 2 * fit$sigma[ordered_biological[1]]
  }
  first_right <- x >= first_centre & x <= first_right_end
  first_right_grey_mass <- max(sum(x[first_right] * y[first_right]), 1)
  first_right_over_count <- sum(
    x[first_right] * pmax(prediction[first_right] - y[first_right], 0)
  )
  first_right_over_mass <- first_right_over_count / first_right_grey_mass
  first_right_over_nrmse <- if (any(first_right)) {
    sqrt(mean((pmax(prediction[first_right] - y[first_right], 0) /
      max(scale, 1))^2))
  } else Inf
  complexity <- 0.002 * max(0, biological_count - 1L)
  score <- 0.38 * nrmse + 0.22 * stats::median(under_mass) +
    0.16 * stats::median(over_mass) + 0.08 * main_peak_nrmse +
    0.08 * first_right_over_nrmse +
    0.08 * first_right_over_mass + complexity
  list(valid = TRUE, score = score, nrmse = nrmse,
       under_mass = stats::median(under_mass),
       over_mass = stats::median(over_mass),
       prefix_under_mass = prefix_under_mass,
       prefix_over_mass = prefix_over_mass,
       prefix_under_fraction = prefix_under_fraction,
       prefix_over_fraction = prefix_over_fraction,
       prefix_under_run = as.integer(prefix_under_run),
       first_right_over_mass = first_right_over_mass,
       first_right_over_nrmse = first_right_over_nrmse,
       main_peak_nrmse = main_peak_nrmse,
       collapsed = as.integer(collapsed),
       minimum_component_share = minimum_component_share,
       maximum_width_jump = maximum_width_jump)
}

# Reconstruct only the missing left half of an inferred first harmonic.  This
# helper is a preprocessing step: downstream ploidy selection and fitting see
# an ordinary histogram and therefore do not need a second selection rule.
.GenomeGMM_complete_hidden_prefix <- function(full_data, first_peak_x) {
  result <- list(valid = FALSE, data = full_data,
                 model_start_index = NA_integer_)
  if (!nrow(full_data) || length(first_peak_x) != 1L ||
      !is.finite(first_peak_x) || first_peak_x <= 0) return(result)
  x <- as.numeric(full_data[, 1])
  y <- as.numeric(full_data[, 2])
  x_steps <- diff(sort(unique(x)))
  x_steps <- x_steps[is.finite(x_steps) & x_steps > 0]
  bin_width <- if (length(x_steps)) stats::median(x_steps) else 1
  sigma_seed <- max(bin_width, 0.35 * first_peak_x)
  model_start_x <- max(min(x), first_peak_x - 2 * sigma_seed)
  model_start_index <- which(x >= model_start_x)[1]
  left_indices <- which(x >= model_start_x & x < first_peak_x)
  if (!length(left_indices)) return(result)
  smooth <- as.numeric(stats::filter(y, rep(1 / 5, 5), sides = 2))
  smooth[!is.finite(smooth)] <- y[!is.finite(smooth)]
  mirror_y <- stats::approx(
    x, smooth, xout = 2 * first_peak_x - x[left_indices], rule = 2
  )$y
  usable <- is.finite(mirror_y) & mirror_y >= 0
  if (!any(usable)) return(result)
  completed <- full_data
  completed[left_indices[usable], 2] <- mirror_y[usable]
  result$valid <- TRUE
  result$data <- completed
  result$model_start_index <- as.integer(model_start_index)
  result
}

# A large, sustained amount of grey k-mer mass missed before the first clear
# peak is a trigger for prefix recovery, regardless of whether one or several
# later peaks are visible. Candidate fits are accepted only when they improve
# the common-domain fit and mass balance while passing component, flank and
# width guards.
.GenomeGMM_try_hidden_prefix_recovery <- function(
    data_file, baseline_fit, pattern_evidence, known_n = NULL,
    species_name = "", .refine_winner = TRUE) {
  result <- list(accepted = FALSE, fit = baseline_fit,
                 selected_slot = NA_integer_, trigger = FALSE,
                 start_x = NA_real_, baseline_metrics = NULL,
                 candidate_table = NULL,
                 prefix_structure = NULL,
                 broad_anchor = NULL,
                 tail_support = NULL,
                 tail_extension_table = NULL,
                 selected_fit_count = NA_integer_,
                 plausible_slots = integer(0),
                 completed_data = NULL,
                 first_peak_x = NA_real_,
                 start_index = NA_integer_,
                 visible_peak_x = NA_real_)
  full_data <- read.table(data_file, header = FALSE)
  full_data <- full_data[
    is.finite(as.numeric(full_data[, 1])) &
      is.finite(as.numeric(full_data[, 2])) &
      as.numeric(full_data[, 1]) > 0 &
      as.numeric(full_data[, 2]) >= 0, , drop = FALSE
  ]
  visible_peaks <- suppressWarnings(as.numeric(pattern_evidence$peaks$x))
  visible_peaks <- sort(unique(visible_peaks[
    is.finite(visible_peaks) & visible_peaks > 0
  ]))
  if (!length(visible_peaks)) {
    broad_anchor <- .GenomeGMM_broad_anchor_candidates(full_data)
    result$broad_anchor <- broad_anchor
    if (!isTRUE(broad_anchor$valid)) return(result)
    visible_peaks <- broad_anchor$selected_x
  }
  first_visible_peak_x <- min(visible_peaks)
  # Assign the rightmost explicit peak to candidate slot m. Thus late adjacent
  # peaks around, for example, 87 and 105 can be tested as harmonics m-1 and m
  # instead of preventing recovery merely because two peaks were detected.
  terminal_visible_peak_x <- max(visible_peaks)
  window <- .GenomeGMM_prefix_mass_window(
    full_data, first_visible_peak_x, terminal_visible_peak_x
  )
  if (!isTRUE(window$valid)) return(result)

  # Detect a sustained slope transition before measuring lost genomic k-mer
  # mass.  The previous order measured the baseline only from window$start_x
  # (about 15 for ZD306), so the large grey/red deficit beginning around the
  # earlier transition (about 7) could never contribute to the trigger.
  prefix_structure <- .GenomeGMM_prefix_structure_candidates(
    full_data, window, terminal_visible_peak_x
  )
  result$prefix_structure <- prefix_structure
  primary_prefix_structure <- NULL
  if (isTRUE(prefix_structure$valid) &&
      nrow(prefix_structure$candidates)) {
    # The first sustained transition after the 5x error floor is the possible
    # start of genomic signal. Later transitions may represent overlapping
    # higher harmonics and must not compete as alternative first peaks.
    earliest_structure <- which.min(
      prefix_structure$candidates$slope_apex_x
    )
    primary_prefix_structure <- prefix_structure$candidates[
      earliest_structure, , drop = FALSE
    ]
    # The curvature maximum is the beginning of the turn and the slope apex is
    # its end. Their midpoint is an integration boundary, not a forced Gaussian
    # centre. Snap it to an observed coverage bin and retain all historical
    # starts for robust multi-window scoring.
    transition_start_x <- mean(c(
      primary_prefix_structure$curvature_peak_x[1],
      primary_prefix_structure$slope_apex_x[1]
    ))
    data_x <- suppressWarnings(as.numeric(full_data[, 1]))
    data_x <- sort(unique(data_x[is.finite(data_x) & data_x >= 5]))
    if (length(data_x) && is.finite(transition_start_x)) {
      transition_start_x <- data_x[
        which.min(abs(data_x - transition_start_x))
      ]
      if (transition_start_x < window$start_x) {
        window$start_x <- transition_start_x
        window$starts <- sort(unique(c(transition_start_x, window$starts)))
        window$start_index <- which.min(abs(
          suppressWarnings(as.numeric(full_data[, 1])) - transition_start_x
        ))
      }
    }
  }

  baseline_count <- if (!is.null(known_n)) {
    as.integer(known_n)
  } else {
    value <- suppressWarnings(as.integer(baseline_fit$fit_harmonic_count))
    if (length(value) != 1L || !is.finite(value) || value < 1L) {
      value <- suppressWarnings(as.integer(baseline_fit$n))
    }
    if (length(value) != 1L || !is.finite(value) || value < 1L) 1L else value
  }
  baseline_metrics <- .GenomeGMM_hidden_prefix_metrics(
    baseline_fit, full_data, window, first_visible_peak_x,
    terminal_visible_peak_x, biological_count = baseline_count
  )
  result$start_x <- window$start_x
  result$baseline_metrics <- baseline_metrics
  # Require both substantial lost mass and a consecutive underfit run. A
  # single noisy point, or a normal error decline followed directly by a
  # well-fitted first peak, therefore cannot activate this branch.
  trigger <- isTRUE(baseline_metrics$valid) &&
    baseline_metrics$prefix_under_mass >= 0.08 &&
    baseline_metrics$prefix_under_fraction >= 0.20 &&
    baseline_metrics$prefix_under_run >= 5L
  result$trigger <- trigger
  if (!trigger) return(result)

  maximum_slot <- if (is.null(known_n)) {
    10L
  } else {
    min(30L, as.integer(known_n))
  }
  # Hidden-first-peak discovery is identical in known- and unknown-ploidy
  # workflows. The visible anchor is compared as slot 2, 3, ... independently
  # of where n came from. Only after this common preprocessing does known_n
  # fix the final number of biological harmonics.
  slots <- if (maximum_slot >= 2L) 2:maximum_slot else integer(0)
  if (!length(slots)) return(result)
  legacy <- .GenomeGMM_get_legacy_environment()
  rows <- vector("list", length(slots))
  fits <- vector("list", length(slots))
  for (i in seq_along(slots)) {
    slot <- slots[i]
    fit_n <- slot
    inferred_first_miu <- terminal_visible_peak_x / slot
    candidate_window <- .GenomeGMM_candidate_prefix_window(
      window, full_data, inferred_first_miu
    )
    comparison_metrics <- .GenomeGMM_hidden_prefix_metrics(
      baseline_fit, full_data, candidate_window, first_visible_peak_x,
      terminal_visible_peak_x, biological_count = baseline_count
    )
    candidate <- NULL
    tryCatch(
      invisible(utils::capture.output(
        candidate <- legacy$.GenomeGMM_fit(
          data_file = data_file, n = fit_n,
          main_peak_weight = 1, symmetry_tolerance = 0.2,
          max_sigma_ratio = 1.35, sigma_mode = "adaptive",
          position_mode = "adaptive", joint_refinement = "off",
          save_plot = FALSE, em_max_iter = 300,
          em_tolerance = 1e-5, extra_component_mode = "never",
          .first_peak_override = terminal_visible_peak_x / slot,
          .start_index_override = candidate_window$start_index,
          .prefix_completion_slot = slot,
          .prefix_completion_weight = 0.35,
          .visible_peak_override = terminal_visible_peak_x,
          .prefix_discovery_only = TRUE
        ), type = "output"
      )), error = function(e) NULL
    )
    metrics <- .GenomeGMM_hidden_prefix_metrics(
      candidate, full_data, candidate_window, first_visible_peak_x,
      terminal_visible_peak_x, biological_count = fit_n
    )
    structure_supported <- !is.null(primary_prefix_structure) &&
      inferred_first_miu >= primary_prefix_structure$lower_miu[1] &&
      inferred_first_miu <= primary_prefix_structure$upper_miu[1]
    recovered_prefix_fraction <- if (
        isTRUE(metrics$valid) &&
        is.finite(comparison_metrics$prefix_under_mass) &&
        comparison_metrics$prefix_under_mass > 0) {
      (comparison_metrics$prefix_under_mass - metrics$prefix_under_mass) /
        comparison_metrics$prefix_under_mass
    } else -Inf
    structure_mass_ok <- isTRUE(structure_supported) &&
      is.finite(recovered_prefix_fraction) &&
      recovered_prefix_fraction >= 0.35 &&
      metrics$prefix_over_mass <= max(
        0.04, 0.25 * comparison_metrics$prefix_under_mass
      )
    component_guard <- metrics$collapsed == 0L ||
      (isTRUE(structure_mass_ok) && metrics$collapsed <= 1L &&
         # The 300-iteration discovery fit is deliberately coarse. Keep a
         # narrow numerical tolerance around the 0.3% floor so a candidate is
         # not discarded before the full fit can redistribute its weights.
         # Final and tail-extension guards retain their stricter thresholds.
         metrics$minimum_component_share >= 0.0025)
    ordinary_eligible <- isTRUE(metrics$valid) && component_guard &&
      metrics$maximum_width_jump <= 1.50 &&
      metrics$main_peak_nrmse <=
        comparison_metrics$main_peak_nrmse + 0.015 &&
      metrics$nrmse <= comparison_metrics$nrmse * 0.90 &&
      metrics$under_mass <= comparison_metrics$under_mass * 0.75 &&
      metrics$over_mass <= comparison_metrics$over_mass + 0.05 &&
      # The fast screen is deliberately looser; the final joint fit plus the
      # dedicated flank refinement below must satisfy the stricter 8% guard.
      metrics$first_right_over_mass <= 0.18 &&
      metrics$score <= comparison_metrics$score * 0.85
    # A structure-supported prefix candidate is first screened without the
    # optional tail component. It may therefore leave more grey mass below the
    # curve at this cheap stage even though its harmonic geometry and prefix
    # reconstruction are substantially better. Allow only such candidates to
    # reach the full joint fit; acceptance below still requires the final model
    # to reduce total underfit mass.
    structured_screen <- isTRUE(metrics$valid) &&
      isTRUE(structure_mass_ok) && component_guard &&
      metrics$maximum_width_jump <= 1.50 &&
      metrics$nrmse <= comparison_metrics$nrmse * 0.55 &&
      metrics$under_mass <= max(
        comparison_metrics$under_mass * 1.50,
        comparison_metrics$under_mass + 0.12
      ) &&
      metrics$over_mass <= comparison_metrics$over_mass + 0.05 &&
      # This is only the 300-iteration discovery fit. A candidate supported by
      # a sustained slope structure and large recovered prefix mass must reach
      # the full optimiser, where its inferred first-peak height and width can
      # be corrected. Applying the final 8% flank guard here discarded the
      # correct ZD306-like geometry and forced selection back to a late peak.
      metrics$first_right_over_mass <= 0.18 &&
      metrics$score <= comparison_metrics$score * 0.75
    eligible <- isTRUE(ordinary_eligible) || isTRUE(structured_screen)
    rows[[i]] <- data.frame(
      visible_peak_slot = slot, fitted_components = fit_n,
      scoring_start_x = candidate_window$start_x,
      baseline_score = comparison_metrics$score,
      baseline_nrmse = comparison_metrics$nrmse,
      baseline_under_mass = comparison_metrics$under_mass,
      score = metrics$score, nrmse = metrics$nrmse,
      under_mass = metrics$under_mass, over_mass = metrics$over_mass,
      prefix_under_mass = metrics$prefix_under_mass,
      prefix_under_fraction = metrics$prefix_under_fraction,
      prefix_under_run = metrics$prefix_under_run,
      inferred_first_miu = inferred_first_miu,
      structure_supported = isTRUE(structure_supported),
      recovered_prefix_fraction = recovered_prefix_fraction,
      structure_mass_ok = isTRUE(structure_mass_ok),
      first_right_over_mass = metrics$first_right_over_mass,
      first_right_over_nrmse = metrics$first_right_over_nrmse,
      main_peak_nrmse = metrics$main_peak_nrmse,
      collapsed_components = metrics$collapsed,
      minimum_component_share = metrics$minimum_component_share,
      maximum_width_jump = metrics$maximum_width_jump,
      eligible = isTRUE(eligible)
    )
    fits[[i]] <- candidate
  }
  table <- do.call(rbind, rows)
  result$candidate_table <- table
  eligible <- which(table$eligible)
  if (!length(eligible)) return(result)
  structured_eligible <- eligible[
    table$structure_supported[eligible] &
      table$structure_mass_ok[eligible]
  ]
  selection_pool <- if (length(structured_eligible)) {
    structured_eligible
  } else eligible
  winner <- selection_pool[which.min(table$score[selection_pool])]
  selected_slot <- as.integer(table$visible_peak_slot[winner])
  selected_first_miu <- terminal_visible_peak_x / selected_slot
  selected_window <- .GenomeGMM_candidate_prefix_window(
    window, full_data, selected_first_miu
  )
  selected_baseline_metrics <- .GenomeGMM_hidden_prefix_metrics(
    baseline_fit, full_data, selected_window, first_visible_peak_x,
    terminal_visible_peak_x, biological_count = baseline_count
  )
  selected_fit_count <- if (is.null(known_n)) {
    selected_slot
  } else {
    as.integer(known_n)
  }
  selected_fast_fit <- fits[[winner]]
  selected_fast_metrics <- .GenomeGMM_hidden_prefix_metrics(
    selected_fast_fit, full_data, selected_window, first_visible_peak_x,
    terminal_visible_peak_x, biological_count = selected_slot
  )

  completion <- .GenomeGMM_complete_hidden_prefix(
    full_data, selected_first_miu
  )
  if (!isTRUE(completion$valid)) return(result)
  base_plausible_slots <- as.integer(table$visible_peak_slot[
    selection_pool[
      table$score[selection_pool] <= table$score[winner] * 1.20
    ]
  ])
  result$completed_data <- completion$data
  result$first_peak_x <- selected_first_miu
  result$start_index <- completion$model_start_index
  result$visible_peak_x <- terminal_visible_peak_x
  result$selected_fit_count <- selected_slot
  result$selected_slot <- selected_slot
  result$plausible_slots <- sort(unique(c(
    base_plausible_slots, selected_slot
  )))
  if (!isTRUE(.refine_winner)) {
    result$accepted <- TRUE
    result$fit <- fits[[winner]]
    return(result)
  }

  # Unknown-ploidy recovery may continue to the right of the anchor when
  # consecutive harmonic windows retain substantial signal.  Add components
  # one at a time and stop at the first unsupported or non-improving fit.  The
  # anchor slot remains separate for diagnostics; selected_fit_count is the
  # final fitted minimum pattern.
  extension_rows <- list()
  if (is.null(known_n)) {
    tail_support <- .GenomeGMM_prefix_tail_support(
      full_data, selected_window, selected_first_miu, selected_slot,
      maximum_additional = 6L
    )
    result$tail_support <- tail_support
    extension_limit <- if (isTRUE(tail_support$valid)) {
      as.integer(tail_support$maximum_slot)
    } else selected_slot
    if (extension_limit > selected_slot) {
      previous_metrics <- selected_fast_metrics
      previous_collapsed <- previous_metrics$collapsed
      for (extension_n in seq.int(selected_slot + 1L, extension_limit)) {
        extension_fit <- NULL
        tryCatch(
          invisible(utils::capture.output(
            extension_fit <- legacy$.GenomeGMM_fit(
              data_file = data_file, n = extension_n,
              main_peak_weight = 1, symmetry_tolerance = 0.2,
              max_sigma_ratio = 1.35, sigma_mode = "adaptive",
              position_mode = "adaptive", joint_refinement = "off",
              save_plot = FALSE, em_max_iter = 300,
              em_tolerance = 1e-5, extra_component_mode = "never",
              .first_peak_override = selected_first_miu,
              .start_index_override = selected_window$start_index,
              .prefix_completion_slot = selected_slot,
              .prefix_completion_weight = 0.35,
              .visible_peak_override = terminal_visible_peak_x
            ), type = "output"
          )), error = function(e) NULL
        )
        extension_metrics <- .GenomeGMM_hidden_prefix_metrics(
          extension_fit, full_data, selected_window, first_visible_peak_x,
          terminal_visible_peak_x, biological_count = extension_n
        )
        last_share <- 0
        if (!is.null(extension_fit) && length(extension_fit$miu) >= extension_n) {
          component_order <- order(extension_fit$miu)
          shares <- pmax(extension_fit$alpha, 0) /
            max(sum(pmax(extension_fit$alpha, 0)), 1e-12)
          last_share <- shares[component_order[extension_n]]
        }
        accepted_extension <- isTRUE(extension_metrics$valid) &&
          extension_metrics$maximum_width_jump <= 1.50 &&
          extension_metrics$collapsed <= max(1L, previous_collapsed) &&
          extension_metrics$minimum_component_share >= 0.003 &&
          is.finite(last_share) && last_share >= 0.008 &&
          extension_metrics$score <= previous_metrics$score * 0.97 &&
          extension_metrics$nrmse <= previous_metrics$nrmse * 0.97 &&
          extension_metrics$under_mass <= previous_metrics$under_mass * 0.95 &&
          extension_metrics$over_mass <=
            min(0.10, previous_metrics$over_mass + 0.03) &&
          extension_metrics$first_right_over_mass <= 0.08
        extension_rows[[length(extension_rows) + 1L]] <- data.frame(
          fitted_components = extension_n,
          score = extension_metrics$score,
          nrmse = extension_metrics$nrmse,
          under_mass = extension_metrics$under_mass,
          over_mass = extension_metrics$over_mass,
          collapsed_components = extension_metrics$collapsed,
          minimum_component_share = extension_metrics$minimum_component_share,
          last_component_share = last_share,
          maximum_width_jump = extension_metrics$maximum_width_jump,
          accepted = isTRUE(accepted_extension)
        )
        if (!isTRUE(accepted_extension)) break
        selected_fit_count <- as.integer(extension_n)
        selected_fast_fit <- extension_fit
        selected_fast_metrics <- extension_metrics
        previous_metrics <- extension_metrics
        previous_collapsed <- extension_metrics$collapsed
      }
    }
  }
  result$selected_fit_count <- selected_fit_count
  result$tail_extension_table <- if (length(extension_rows)) {
    do.call(rbind, extension_rows)
  } else data.frame()
  accepted_extension_slots <- if (length(extension_rows)) {
    extension_table <- do.call(rbind, extension_rows)
    as.integer(extension_table$fitted_components[
      extension_table$accepted %in% TRUE
    ])
  } else integer(0)
  result$plausible_slots <- sort(unique(as.integer(c(
    base_plausible_slots, selected_slot, accepted_extension_slots
  ))))
  fit_n <- if (is.null(known_n)) selected_fit_count else as.integer(known_n)
  refined <- NULL
  tryCatch(
    invisible(utils::capture.output(
      refined <- legacy$.GenomeGMM_fit(
        data_file = data_file, n = fit_n,
        species_name = species_name,
        main_peak_weight = 1, symmetry_tolerance = 0.2,
        max_sigma_ratio = 1.35, sigma_mode = "adaptive",
        position_mode = "adaptive", joint_refinement = "auto",
        extra_component_mode = "legacy", save_plot = TRUE,
        em_max_iter = 10000, em_tolerance = 1e-6,
        .first_peak_override = terminal_visible_peak_x / selected_slot,
        .start_index_override = selected_window$start_index,
        .prefix_completion_slot = selected_slot,
        .prefix_completion_weight = 0.35,
        .visible_peak_override = terminal_visible_peak_x
      ), type = "output"
    )), error = function(e) NULL
  )
  refined_metrics <- .GenomeGMM_hidden_prefix_metrics(
    refined, full_data, selected_window, first_visible_peak_x,
    terminal_visible_peak_x, biological_count = fit_n
  )
  refined_recovered_fraction <- if (
      isTRUE(refined_metrics$valid) &&
      is.finite(selected_baseline_metrics$prefix_under_mass) &&
      selected_baseline_metrics$prefix_under_mass > 0) {
    (selected_baseline_metrics$prefix_under_mass -
       refined_metrics$prefix_under_mass) /
      selected_baseline_metrics$prefix_under_mass
  } else -Inf
  refined_structure_guard <- isTRUE(
    table$structure_supported[winner]
  ) && is.finite(refined_recovered_fraction) &&
    refined_recovered_fraction >= 0.35 &&
    refined_metrics$prefix_over_mass <= max(
      0.04, 0.25 * selected_baseline_metrics$prefix_under_mass
    )
  refined_component_guard <- refined_metrics$collapsed == 0L ||
    (isTRUE(refined_structure_guard) && refined_metrics$collapsed <= 1L &&
       refined_metrics$minimum_component_share >= 0.003)
  ordinary_refined_ok <- isTRUE(refined_metrics$valid) &&
    refined_component_guard &&
    refined_metrics$maximum_width_jump <= 1.50 &&
    refined_metrics$main_peak_nrmse <=
      selected_baseline_metrics$main_peak_nrmse + 0.015 &&
    refined_metrics$nrmse <= selected_baseline_metrics$nrmse * 0.90 &&
    refined_metrics$under_mass <=
      selected_baseline_metrics$under_mass * 0.75 &&
    refined_metrics$over_mass <=
      selected_baseline_metrics$over_mass + 0.05 &&
    refined_metrics$first_right_over_mass <= 0.08 &&
    refined_metrics$score <= selected_baseline_metrics$score * 0.85
  structured_refined_ok <- isTRUE(refined_metrics$valid) &&
    isTRUE(refined_structure_guard) && refined_component_guard &&
    refined_metrics$maximum_width_jump <= 1.50 &&
    refined_metrics$main_peak_nrmse <=
      selected_baseline_metrics$main_peak_nrmse + 0.020 &&
    refined_metrics$nrmse <= selected_baseline_metrics$nrmse * 0.55 &&
    refined_metrics$under_mass <=
      selected_baseline_metrics$under_mass * 0.90 &&
    refined_metrics$over_mass <=
      selected_baseline_metrics$over_mass + 0.05 &&
    refined_metrics$first_right_over_mass <= 0.08 &&
    refined_metrics$score <= selected_baseline_metrics$score * 0.75
  refined_ok <- isTRUE(ordinary_refined_ok) ||
    isTRUE(structured_refined_ok)
  if (!isTRUE(refined_ok)) {
    if (!is.null(refined$pdf_file) && length(refined$pdf_file) == 1L &&
        is.character(refined$pdf_file) && file.exists(refined$pdf_file)) {
      unlink(refined$pdf_file)
    }
    return(result)
  }
  result$accepted <- TRUE
  result$fit <- refined
  result$selected_slot <- selected_slot
  result
}

# Known-ploidy interpretation helpers.
# An allopolyploid composition is represented by a non-increasing integer
# partition of the supplied biological ploidy.  For example, 6x includes
# 6 (AAAAAA), 5+1 (AAAAAB), 4+2 (AAAABB), 4+1+1 (AAAABC), 3+3
# (AAABBB), 3+2+1 (AAABBC), and 2+2+2 (AABBCC).  Dividing each partition by
# the gcd of its copy counts gives the irreducible base-pattern count used by
# the curve fitter.  Compositions with the same base count reuse one fit.
.GenomeGMM_gcd_pair <- function(a, b) {
  a <- abs(as.integer(a))
  b <- abs(as.integer(b))
  while (b != 0L) {
    remainder <- a %% b
    a <- b
    b <- remainder
  }
  a
}

.GenomeGMM_gcd_vector <- function(values) {
  values <- abs(as.integer(values))
  values <- values[is.finite(values)]
  if (!length(values)) return(1L)
  Reduce(.GenomeGMM_gcd_pair, values)
}

.GenomeGMM_integer_partitions <- function(total) {
  total <- as.integer(total)
  result <- list()
  append_partition <- function(remaining, maximum, prefix) {
    if (remaining == 0L) {
      result[[length(result) + 1L]] <<- as.integer(prefix)
      return(invisible(NULL))
    }
    upper <- min(maximum, remaining)
    for (value in seq.int(upper, 1L)) {
      append_partition(
        remaining - value, value, c(prefix, value)
      )
    }
    invisible(NULL)
  }
  append_partition(total, total, integer())
  result
}

.GenomeGMM_composition_label <- function(copy_counts) {
  copy_counts <- as.integer(copy_counts)
  symbol <- function(index) {
    if (index <= length(LETTERS)) LETTERS[index] else paste0("S", index)
  }
  render <- function(counts) {
    paste0(vapply(seq_along(counts), function(index) {
      paste(rep(symbol(index), counts[index]), collapse = "")
    }, character(1L)), collapse = "")
  }
  primary <- render(copy_counts)
  if (length(copy_counts) == 1L) {
    return(paste0(primary, " / ", paste(rep("B", copy_counts), collapse = "")))
  }
  if (length(copy_counts) == 2L && copy_counts[1L] != copy_counts[2L]) {
    swapped <- paste0(
      paste(rep("A", copy_counts[2L]), collapse = ""),
      paste(rep("B", copy_counts[1L]), collapse = "")
    )
    return(paste0(primary, " / ", swapped))
  }
  primary
}

.GenomeGMM_enumerate_allopoly_compositions <- function(actual_ploidy) {
  partitions <- .GenomeGMM_integer_partitions(actual_ploidy)
  rows <- lapply(partitions, function(copy_counts) {
    common_factor <- .GenomeGMM_gcd_vector(copy_counts)
    reduced <- as.integer(copy_counts %/% common_factor)
    base_pattern <- sum(reduced)
    # Subgenome-specific k-mers occur at the reduced copy-number slots;
    # k-mers shared by the whole complement may occur at base_pattern.
    expected_slots <- sort(unique(c(reduced, base_pattern)))
    expected_ratios <- expected_slots / min(expected_slots)
    data.frame(
      composition = .GenomeGMM_composition_label(copy_counts),
      copy_partition = paste(copy_counts, collapse = "+"),
      reduced_partition = paste(reduced, collapse = "+"),
      base_pattern = as.integer(base_pattern),
      repeated_pattern_multiplier = as.integer(common_factor),
      expected_peak_count = length(expected_slots),
      expected_peak_slots = paste(expected_slots, collapse = ","),
      expected_peak_ratios = paste(
        format(round(expected_ratios, 4), trim = TRUE, scientific = FALSE),
        collapse = ","
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.GenomeGMM_composition_peak_assessment <- function(compositions, evidence) {
  observed_ratios <- numeric()
  if (!is.null(evidence$peaks) && nrow(evidence$peaks)) {
    observed_peaks <- evidence$peaks
    if ("trusted_harmonic" %in% names(observed_peaks) &&
        any(observed_peaks$trusted_harmonic, na.rm = TRUE)) {
      observed_peaks <- observed_peaks[
        observed_peaks$trusted_harmonic, , drop = FALSE
      ]
    }
    observed_ratios <- suppressWarnings(as.numeric(observed_peaks$ratio))
    observed_ratios <- sort(unique(observed_ratios[
      is.finite(observed_ratios) & observed_ratios > 0
    ]))
  }
  if (!length(observed_ratios)) observed_ratios <- 1
  observed_ratios <- observed_ratios / min(observed_ratios)

  score_one <- function(expected_text) {
    expected <- suppressWarnings(as.numeric(strsplit(
      expected_text, ",", fixed = TRUE
    )[[1L]]))
    expected <- sort(unique(expected[is.finite(expected) & expected > 0]))
    if (!length(expected)) return(Inf)
    expected <- expected / min(expected)
    observed_to_expected <- mean(vapply(observed_ratios, function(value) {
      min(abs(log(value / expected)))
    }, numeric(1L)))
    expected_to_observed <- mean(vapply(expected, function(value) {
      min(abs(log(value / observed_ratios)))
    }, numeric(1L)))
    count_difference <- abs(length(observed_ratios) - length(expected)) /
      max(length(observed_ratios), length(expected), 1L)
    0.45 * observed_to_expected + 0.45 * expected_to_observed +
      0.10 * count_difference
  }
  compositions$peak_evidence_score <- vapply(
    compositions$expected_peak_ratios, score_one, numeric(1L)
  )
  compositions$observed_peak_count <- length(observed_ratios)
  compositions$observed_peak_ratios <- paste(
    format(round(observed_ratios, 4), trim = TRUE, scientific = FALSE),
    collapse = ","
  )
  compositions
}

# Resolve the internal harmonic pattern only after the user has supplied the
# biological ploidy type. For an autopolyploid the fitted pattern is n. For an
# allopolyploid, all integer-partition compositions are enumerated and grouped
# by their irreducible base-pattern count. The existing automatic selector
# supplies full-curve evidence; explicit peak count and relative positions add
# composition support without replacing the protected fitter.
.GenomeGMM_resolve_known_pattern <- function(data_file, actual_ploidy,
                                             ploidy_type, species_name = "") {
  actual_ploidy <- as.integer(actual_ploidy)
  if (identical(ploidy_type, "autopolyploid")) {
    return(list(
      status = "ok", actual_ploidy = actual_ploidy,
      base_pattern = actual_ploidy, multiplier = 1L,
      potential_patterns = actual_ploidy,
      legal_patterns = actual_ploidy,
      automatic_target = actual_ploidy,
      terminal_supported_pattern = actual_ploidy,
      selection_basis = "autopolyploid: fitted pattern equals supplied ploidy"
    ))
  }

  compositions <- .GenomeGMM_enumerate_allopoly_compositions(actual_ploidy)
  # Retain the one-part composition (AAAA...) as the homozygous single-source
  # endpoint requested by the composition table. It maps to base pattern 1;
  # multi-source irreducible compositions map to pattern 2 or greater.
  potential_patterns <- sort(unique(as.integer(
    compositions$base_pattern[compositions$base_pattern >= 1L]
  )))
  if (!length(potential_patterns)) potential_patterns <- actual_ploidy

  inference <- .GenomeGMM_auto(
    data_file = data_file, n = NULL, n_candidates = 2:10,
    species_name = species_name, search = "fast",
    candidate_max_iter = 1500, candidate_tolerance = 1e-5,
    verification_top = 4, verification_max_iter = 10000,
    final_max_iter = 10000, save_search = FALSE, save_plot = FALSE,
    .selection_only = TRUE
  )
  scalar_integer <- function(value) {
    value <- suppressWarnings(as.integer(value))
    if (length(value) == 1L && is.finite(value)) value else NA_integer_
  }
  automatic_target <- scalar_integer(inference$fit_harmonic_count)
  if (!is.finite(automatic_target) || automatic_target < 1L) {
    automatic_target <- scalar_integer(inference$minimum_pattern)
  }
  terminal_target <- scalar_integer(inference$terminal_harmonic_slot)
  genuine_single_peak <- .GenomeGMM_is_genuine_single_peak(
    inference$global_pattern_evidence
  )
  if (isTRUE(genuine_single_peak)) {
    automatic_target <- 1L
    terminal_target <- 1L
  }
  if (!is.finite(terminal_target) || terminal_target < 1L) {
    terminal_target <- automatic_target
  }
  if (!is.finite(automatic_target) || automatic_target < 1L) {
    automatic_target <- terminal_target
  }
  if (!is.finite(automatic_target) || automatic_target < 1L) {
    automatic_target <- potential_patterns[
      ceiling(length(potential_patterns) / 2)
    ]
  }
  if (!is.finite(terminal_target) || terminal_target < 1L) {
    terminal_target <- automatic_target
  }

  compositions <- .GenomeGMM_composition_peak_assessment(
    compositions, inference$global_pattern_evidence
  )
  best_composition_score <- vapply(potential_patterns, function(pattern) {
    values <- compositions$peak_evidence_score[
      compositions$base_pattern == pattern
    ]
    if (length(values) && any(is.finite(values))) {
      min(values[is.finite(values)])
    } else 1
  }, numeric(1L))
  # The automatic target already combines protected full-curve GMM quality,
  # shoulder structure and sustained decay. Peak count/ratio agreement is a
  # secondary composition-aware discriminator, so incomplete or overlapping
  # peaks cannot override a clearly better protected curve.
  selection_score <- abs(log(potential_patterns / automatic_target)) +
    0.25 * abs(log(potential_patterns / terminal_target)) +
    0.20 * best_composition_score
  selected_index <- which.min(selection_score)
  selected <- potential_patterns[selected_index]
  score_order <- order(selection_score)
  score_gap <- if (length(score_order) >= 2L) {
    selection_score[score_order[2L]] - selection_score[score_order[1L]]
  } else Inf
  confidence <- if (!is.finite(score_gap) || score_gap >= 0.35) {
    "high"
  } else if (score_gap >= 0.15) {
    "medium"
  } else {
    "low; more than one composition model is similarly supported"
  }
  pattern_assessment <- data.frame(
    base_pattern = potential_patterns,
    peak_evidence_score = best_composition_score,
    selection_score = selection_score,
    selected = seq_along(potential_patterns) == selected_index,
    stringsAsFactors = FALSE
  )
  selected_composition_count <- sum(
    compositions$base_pattern == selected
  )
  list(
    status = "ok", actual_ploidy = actual_ploidy,
    base_pattern = as.integer(selected),
    multiplier = actual_ploidy %/% as.integer(selected),
    potential_patterns = as.integer(potential_patterns),
    legal_patterns = as.integer(potential_patterns),
    compositions = compositions,
    pattern_assessment = pattern_assessment,
    base_pattern_confidence = confidence,
    composition_confidence = confidence,
    composition_identifiable = selected_composition_count == 1L,
    automatic_target = as.integer(automatic_target),
    terminal_supported_pattern = as.integer(terminal_target),
    selection_score = selection_score,
    selection_basis = paste0(
      "best combined support from the protected peak-and-fit target ",
      automatic_target, "x, terminal harmonic ", terminal_target,
      "x, and composition-specific peak count/position evidence"
    )
  )
}

.GenomeGMM_known_result_path <- function(data_file, actual_ploidy,
                                         pdf_file = NULL) {
  if (!is.null(pdf_file) && length(pdf_file) == 1L &&
      is.character(pdf_file) && !is.na(pdf_file) && nzchar(pdf_file)) {
    return(paste0(tools::file_path_sans_ext(pdf_file), "_result.txt"))
  }
  stem <- tools::file_path_sans_ext(basename(data_file))
  paste0("em_fitted_gmm_", stem, "_n", as.integer(actual_ploidy), "_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), "_result.txt")
}

.GenomeGMM_write_unfactorable_allopoly_report <- function(
    data_file, actual_ploidy, species_name = "") {
  result_file <- .GenomeGMM_known_result_path(data_file, actual_ploidy)
  sample_name <- tools::file_path_sans_ext(basename(data_file))
  lines <- c(
    paste0("GenomeGMM known-ploidy interpretation for sample: ", sample_name),
    if (nzchar(species_name)) paste0("Species: ", species_name) else NULL,
    "",
    paste0("Provided biological ploidy: ", actual_ploidy, "x"),
    "Provided ploidy type: allopolyploid",
    "",
    paste0(
      "No balanced repeated allopolyploid base pattern can be fitted: ",
      actual_ploidy, "x has no proper integer divisor >= 2."
    ),
    "No PDF was generated.",
    paste0(
      "Please verify the supplied ploidy/type. If the genome should be ",
      "modelled as an autopolyploid, rerun with ploidy_type = ",
      "\"autopolyploid\"."
    ),
    "",
    "This message applies to the balanced repeated-pattern interpretation;",
    "an unbalanced or otherwise more complex allopolyploid requires additional prior information."
  )
  writeLines(enc2utf8(lines), result_file, useBytes = TRUE)
  # Keep a plain list return, matching the public no-peak pathway.
  list(
    status = "no_balanced_factorization", result_file = result_file,
    pdf_file = NULL, actual_ploidy = as.integer(actual_ploidy),
    ploidy_type = "allopolyploid"
  )
}

.GenomeGMM_rename_known_pdf <- function(pdf_file, data_file, actual_ploidy) {
  if (is.null(pdf_file) || length(pdf_file) != 1L || is.na(pdf_file) ||
      !nzchar(pdf_file) || !file.exists(pdf_file)) return(pdf_file)
  old_path <- pdf_file
  match <- regexec("_n[0-9]+_([0-9]{8}_[0-9]{6})[.]pdf$", basename(old_path))
  parts <- regmatches(basename(old_path), match)[[1L]]
  timestamp <- if (length(parts) >= 2L) parts[2L] else
    format(Sys.time(), "%Y%m%d_%H%M%S")
  stem <- tools::file_path_sans_ext(basename(data_file))
  new_name <- paste0("em_fitted_gmm_", stem, "_n",
                     as.integer(actual_ploidy), "_", timestamp, ".pdf")
  new_path <- file.path(dirname(old_path), new_name)
  if (!identical(normalizePath(old_path, mustWork = FALSE),
                 normalizePath(new_path, mustWork = FALSE))) {
    if (file.exists(new_path)) stop("Known-ploidy PDF already exists: ", new_path)
    if (!file.rename(old_path, new_path)) stop("Could not rename PDF: ", old_path)
  }
  if (identical(dirname(pdf_file), ".")) new_name else new_path
}

.GenomeGMM_write_known_allopoly_report <- function(data_file, species_name,
                                                    fit, resolution) {
  result_file <- .GenomeGMM_known_result_path(
    data_file, resolution$actual_ploidy, fit$pdf_file
  )
  sample_name <- tools::file_path_sans_ext(basename(data_file))
  mono_mb <- as.numeric(fit$g_hap_bp) / 1e6
  holo_mb <- as.numeric(fit$genome_size_bp) / 1e6
  potential_patterns <- if (!is.null(resolution$potential_patterns)) {
    resolution$potential_patterns
  } else {
    resolution$legal_patterns
  }
  potential_text <- if (length(potential_patterns)) {
    paste0(potential_patterns, "x", collapse = ", ")
  } else {
    "none"
  }
  candidate_reason_lines <- if (!is.null(resolution$pattern_assessment) &&
      nrow(resolution$pattern_assessment)) {
    vapply(seq_len(nrow(resolution$pattern_assessment)), function(index) {
      row <- resolution$pattern_assessment[index, , drop = FALSE]
      sprintf(
        "  %dx: %s; peak-pattern mismatch %.4f; combined selection cost %.4f.",
        row$base_pattern,
        if (isTRUE(row$selected)) {
          "selected (lowest cost)"
        } else {
          "not selected"
        },
        row$peak_evidence_score, row$selection_score
      )
    }, character(1L))
  } else {
    paste0("  ", resolution$base_pattern, "x: selected.")
  }
  composition_lines <- if (!is.null(resolution$compositions) &&
      nrow(resolution$compositions)) {
    composition_header <- sprintf(
      "%-32s %-13s %-12s %-14s %s",
      "Composition", "Copy counts", "Base pattern", "Expected peaks",
      "Relative peak positions"
    )
    composition_rows <- vapply(
      seq_len(nrow(resolution$compositions)), function(index) {
        row <- resolution$compositions[index, , drop = FALSE]
        sprintf(
          "%-32s %-13s %-12s %-14s %s",
          row$composition, row$copy_partition,
          paste0(row$base_pattern, "x"), row$expected_peak_count,
          row$expected_peak_ratios
        )
      }, character(1L)
    )
    c(composition_header, composition_rows)
  } else {
    "No composition enumeration was available."
  }
  selected_compositions <- if (!is.null(resolution$compositions) &&
      nrow(resolution$compositions)) {
    selected_rows <- resolution$compositions[
      resolution$compositions$base_pattern == resolution$base_pattern,
      , drop = FALSE
    ]
    paste(selected_rows$composition, collapse = "; ")
  } else {
    "not available"
  }
  unit_mb <- as.numeric(fit$subgenome_unit_size_bp) / 1e6
  lines <- c(
    paste0("GenomeGMM known-ploidy interpretation for sample: ", sample_name),
    if (nzchar(species_name)) paste0("Species: ", species_name) else NULL,
    "",
    paste0("Provided biological ploidy: ", resolution$actual_ploidy, "x"),
    "Provided ploidy type: allopolyploid",
    paste0("Potential base-pattern candidates: ", potential_text),
    paste0("Selected base pattern for curve fitting: ",
           resolution$base_pattern, "x"),
    paste0("Repeated-pattern multiplier: ", resolution$multiplier),
    paste0("Selection basis: ", resolution$selection_basis),
    paste0("Peak-and-fit automatic target: ",
           resolution$automatic_target, "x"),
    paste0("Terminal supported harmonic: ",
           resolution$terminal_supported_pattern, "x"),
    paste0(
      "Candidate assessment (lower mismatch/cost indicates better ",
      "agreement; 0 represents a perfect match):"
    ),
    candidate_reason_lines,
    paste0("Base-pattern selection confidence: ",
           resolution$base_pattern_confidence),
    "",
    "Potential genome compositions:",
    composition_lines,
    "",
    paste0("Compositions compatible with the selected base pattern: ",
           selected_compositions),
    paste0(
      "Unique composition identifiable from the enumerated models: ",
      if (isTRUE(resolution$composition_identifiable)) "yes" else "no"
    ),
    paste0(
      "These compositions share the selected base pattern. If their peak ",
      "evidence is indistinguishable, the one-dimensional k-mer histogram ",
      "does not identify a unique composition."
    ),
    if (!is.null(resolution$warning) && nzchar(resolution$warning)) {
      paste0("Warning: ", resolution$warning)
    } else NULL,
    "",
    sprintf(
      "Subgenome size: %.2f Mb", unit_mb
    ),
    sprintf("Monoploid genome size: %.2f Mb", mono_mb),
    sprintf("Holoploid genome size (%dx): %.2f Mb",
            resolution$actual_ploidy, holo_mb),
    sprintf(
      "Monoploid calculation: %.2f Mb * %d = %.2f Mb",
      unit_mb, resolution$base_pattern, mono_mb
    ),
    sprintf(
      "Holoploid calculation: %.2f Mb * %d = %.2f Mb",
      mono_mb, resolution$multiplier, holo_mb
    ),
    "",
    paste0(
      "The fitted PDF uses the ", resolution$base_pattern,
      "x k-mer base pattern. Holoploid genome size is obtained by multiplying ",
      "Monoploid genome size by the repeated-pattern multiplier ",
      resolution$multiplier, " (", resolution$actual_ploidy, "x / ",
      resolution$base_pattern, "x)."
    ),
    "",
    "Definitions:",
    "Monoploid genome size: the DNA content of one basic chromosome set (1Cx).",
    "Holoploid genome size: the total DNA content of the complete chromosome complement at the corresponding ploidy level."
  )
  writeLines(enc2utf8(lines), result_file, useBytes = TRUE)
  result_file
}

#' Automatically infer the minimum ploidy supported by the fitted spectrum
#'
#' A one-dimensional k-mer histogram cannot uniquely identify absolute
#' biological ploidy. GenomeGMM therefore selects a minimum fitted ploidy m
#' and reports the compatible family m, 2m, 3m, ... together with g-hap.
#' Explicit peak evidence is retained separately for confidence assessment.
#'
#' @param data_file Path to a two-column k-mer histogram file.
#' @param n Internal fitted base-pattern count. The public GenomeGMM interface
#'        resolves this value from biological ploidy and ploidy_type.
#' @param species_name Optional species name shown in the plot title.
#'
#' @return An invisible list containing the fitted minimum pattern, compatible
#'         ploidy family, g-hap estimate, PDF path and text-report path.
#'
#' @examples
#' # GenomeGMM("sample.histo")
#' # GenomeGMM("sample.histo", species_name = "Example")
#'
#' @export
.GenomeGMM_run <- function(data_file, n = NULL, species_name = "",
                           monoploid_range = NULL) {
  if (!is.null(n) && (!is.numeric(n) || length(n) != 1L ||
      !is.finite(n) || n < 1L || n != round(n))) {
    stop("Parameter 'n' must be a known base-pattern integer >= 1.")
  }
  if (!is.character(species_name) || length(species_name) != 1L ||
      is.na(species_name)) {
    stop("Parameter 'species_name' must be one character string.")
  }

  inference_fit <- NULL
  single_peak_ambiguity <- FALSE
  two_peak_shoulder_assessment <- NULL
  recovery_pattern_evidence <- NULL
  if (is.null(n)) {
    cat("\n[Automatic ploidy inference only]\n")
    cat("The auto model will select minimum_pattern and fit_n;",
        "its temporary curve will not be plotted.\n")
    inference_fit <- .GenomeGMM_auto(
      data_file = data_file,
      n = NULL,
      n_candidates = 2:10,
      species_name = species_name,
      search = "fast",
      candidate_max_iter = 1500,
      candidate_tolerance = 1e-5,
      verification_top = 4,
      verification_max_iter = 10000,
      final_max_iter = 10000,
      save_search = FALSE,
      save_plot = FALSE,
      .selection_only = TRUE,
      .monoploid_range = monoploid_range
    )
    single_peak_ambiguity <- .GenomeGMM_is_genuine_single_peak(
      inference_fit$global_pattern_evidence
    )
    if (isTRUE(single_peak_ambiguity)) {
      two_peak_shoulder_assessment <-
        .GenomeGMM_assess_hidden_two_peak_shoulder(
          data_file, inference_fit$global_pattern_evidence
        )
      if (isTRUE(two_peak_shoulder_assessment$supported)) {
        single_peak_ambiguity <- FALSE
        cat("Single-maximum shoulder check: retained two harmonic peaks\n")
        cat("One-/two-peak NRMSE:",
            round(two_peak_shoulder_assessment$nrmse_one, 5), "/",
            round(two_peak_shoulder_assessment$nrmse_two, 5), "\n")
        cat("Shoulder/main candidate miu:",
            round(two_peak_shoulder_assessment$first_miu, 3), "/",
            round(two_peak_shoulder_assessment$second_miu, 3), "\n")
      } else {
        cat("Single-maximum shoulder check: genuine single-peak ambiguity\n")
      }
    }
    if (isTRUE(single_peak_ambiguity) && !is.null(monoploid_range)) {
      # A genuinely single-peaked spectrum is compatible with A, AA, AAA,
      # ... . An expected Monoploid-size interval can check the size estimate
      # but cannot select an absolute ploidy or a different curve extent.
      inference_fit$range_selected_ploidy <- NA_integer_
      inference_fit$range_terminal_anchor_slot <- NA_integer_
      inference_fit$range_selection_basis <-
        "single-peak ambiguity; no unique ploidy selected"
    }
    recovery_pattern_evidence <- inference_fit$global_pattern_evidence
    fitted_ploidy <- if (isTRUE(single_peak_ambiguity)) {
      1L
    } else {
      as.integer(inference_fit$fit_harmonic_count)
    }
    # The final fitted harmonic count is the smallest ploidy interpretation
    # supported by the complete fitted spectrum.  Keep the rightmost explicit
    # harmonic separately: it measures direct peak evidence and therefore the
    # confidence of this inference, but does not redefine the reported family.
    minimum_pattern <- if (isTRUE(single_peak_ambiguity)) {
      1L
    } else {
      fitted_ploidy
    }
    observable_pattern <- suppressWarnings(as.integer(
      inference_fit$terminal_harmonic_slot
    ))
    if (!is.finite(observable_pattern) || observable_pattern < 1L) {
      observable_pattern <- 1L
    }
  } else {
    minimum_pattern <- as.integer(n)
    fitted_ploidy <- as.integer(n)
    observable_pattern <- as.integer(n)
    # Known-n calls keep their supplied fitted base-pattern count. Reuse the global
    # peak evidence produced by the preserved fitter below; it already applies
    # the fit engine's data-supported analysis limit. A second full-tail scan
    # is both redundant and very slow for long histograms.
    recovery_pattern_evidence <- NULL
  }

  # When a supplied size range matches only the balanced-allopolyploid
  # interpretation, the final base-pattern PDF should display one complete
  # distinct-subgenome complement (for example ABD), not one individual
  # per-pattern unit (A, B or D). The fitted curve itself is unchanged.
  if (is.null(n) && !is.null(inference_fit) && identical(
      inference_fit$range_selected_size_model,
      "balanced_allopolyploid"
    )) {
    previous_range_divisor <- getOption(
      "GenomeGMM.reporting_monoploid_divisor"
    )
    options(GenomeGMM.reporting_monoploid_divisor = 1L)
    on.exit(
      options(
        GenomeGMM.reporting_monoploid_divisor = previous_range_divisor
      ),
      add = TRUE
    )
  }

  # Prefix detection first uses the same complete fixed-n fit that is needed
  # for an ordinary final result.  Keep its PDF and fitted parameters when no
  # hidden prefix is accepted, instead of repeating the identical EM search
  # solely to draw the plot.  If prefix recovery is accepted, this provisional
  # ordinary PDF is removed and the recovered-prefix pathway is fitted exactly
  # as before.
  legacy_environment <- .GenomeGMM_get_legacy_environment()
  if (is.null(recovery_pattern_evidence)) {
    known_full_data <- read.table(data_file, header = FALSE)
    known_full_data <- known_full_data[
      is.finite(as.numeric(known_full_data[, 1])) &
        is.finite(as.numeric(known_full_data[, 2])) &
        as.numeric(known_full_data[, 1]) > 0,
      , drop = FALSE
    ]
    known_analysis_limit <- .GenomeGMM_auto_analysis_limit(
      known_full_data, n = fitted_ploidy
    )
    recovery_pattern_evidence <- detect_global_pattern_peaks(
      known_full_data, analysis_max_x = known_analysis_limit
    )
  }
  shallow_anchor_active <- isTRUE(
    recovery_pattern_evidence$shallow_anchor_supported
  ) && is.finite(recovery_pattern_evidence$shallow_anchor_x)
  baseline_args <- list(
    data_file = data_file,
    n = fitted_ploidy,
    species_name = species_name,
    search = "fast",
    candidate_max_iter = 1500,
    candidate_tolerance = 1e-5,
    verification_top = 4,
    verification_max_iter = 10000,
    final_max_iter = 10000,
    save_search = FALSE,
    save_plot = TRUE,
    joint_refinement = "auto"
  )
  if (shallow_anchor_active) {
    baseline_args$.first_peak_override <-
      recovery_pattern_evidence$shallow_anchor_x
    baseline_args$.start_index_override <-
      recovery_pattern_evidence$start_index
    if (identical(recovery_pattern_evidence$shallow_anchor_kind,
                  "shoulder") &&
        is.finite(recovery_pattern_evidence$shallow_anchor_terminal_x)) {
      # A shoulder supplies the fitting start but not a complete first peak.
      # Reuse the established low-weight symmetry completion for its missing
      # side, then run the same ordinary forward harmonic fitter.  A genuine
      # local maximum needs no synthetic completion.
      baseline_args$.prefix_completion_slot <- fitted_ploidy
      baseline_args$.prefix_completion_weight <- 0.35
      baseline_args$.visible_peak_override <-
        recovery_pattern_evidence$shallow_anchor_terminal_x
    }
  }
  baseline_fit <- do.call(
    legacy_environment$.GenomeGMM_auto, baseline_args
  )
  recovery <- .GenomeGMM_try_hidden_prefix_recovery(
    data_file = data_file,
    baseline_fit = baseline_fit,
    pattern_evidence = recovery_pattern_evidence,
    known_n = n,
    species_name = species_name,
    .refine_winner = FALSE
  )
  if (shallow_anchor_active) {
    # This path has already recovered an observed (not synthetic) first peak.
    # It must now use the ordinary forward fitter and must not be reinterpreted
    # a second time by hidden-prefix completion.
    recovery$trigger <- FALSE
    recovery$accepted <- FALSE
    recovery$reason <- "post_noise_shallow_anchor_forward_fit"
  }
  # The overwhelmingly common path has an observed first peak.  Reuse the
  # already completed baseline fit verbatim; this preserves n, alpha, miu,
  # sigma, genome-size calculations and the plotted curve while avoiding a
  # second fixed-n EM run.
  fit <- baseline_fit
  if (isTRUE(recovery$accepted)) {
    # The provisional ordinary curve is not a valid deliverable once prefix
    # recovery wins.  Remove only that newly generated intermediate PDF, then
    # run the existing recovered-prefix final fit unchanged.
    baseline_pdf <- baseline_fit$pdf_file
    if (!is.null(baseline_pdf) && length(baseline_pdf) == 1L &&
        is.character(baseline_pdf) && !is.na(baseline_pdf) &&
        nzchar(baseline_pdf) && file.exists(baseline_pdf)) {
      unlink(baseline_pdf)
    }
    fit <- NULL
    # Prefix discovery has already compared the terminal visible peak as slot
    # 2, 3, ... on one common observed domain.  Do not discard that result by
    # restarting an unconstrained minimum-pattern search on the completed
    # left half: that second search can reinterpret a late explicit harmonic
    # as a nuisance tail and collapse a genuine 4/5-peak pattern to 2x.
    # Completion identifies the first peak and supplies a hard lower bound
    # from the terminal explicit anchor. Re-run the ordinary full-spectrum
    # selector on the completed data to choose the numerical tail extent, but
    # never allow it to collapse below that anchor-supported bound. This lets
    # a 6-slot prefix extend directly to a better 8-component curve even when
    # the intermediate 7-component fit is poor.
    if (is.null(n)) {
      recovered_floor <- suppressWarnings(as.integer(
        recovery$selected_slot
      ))
      unified_inference <- tryCatch(
        .GenomeGMM_auto(
          data_file = data_file,
          n = NULL,
          n_candidates = 2:10,
          species_name = species_name,
          search = "fast",
          candidate_max_iter = 1500,
          candidate_tolerance = 1e-5,
          verification_top = 4,
          verification_max_iter = 10000,
          final_max_iter = 10000,
          save_search = FALSE,
          save_plot = FALSE,
          .selection_only = TRUE,
          .preloaded_data = recovery$completed_data,
          .monoploid_range = monoploid_range
        ),
        error = function(error) NULL
      )
      range_selected_count <- if (!is.null(unified_inference)) {
        suppressWarnings(as.integer(
          unified_inference$range_selected_ploidy
        ))
      } else NA_integer_
      tail_selected_count <- if (is.finite(range_selected_count)) {
        range_selected_count
      } else if (!is.null(unified_inference)) {
        suppressWarnings(as.integer(unified_inference$fit_harmonic_count))
      } else {
        NA_integer_
      }
      fitted_ploidy <- max(
        c(recovered_floor, tail_selected_count), na.rm = TRUE
      )
      if (!is.finite(fitted_ploidy) || fitted_ploidy < 2L) {
        fitted_ploidy <- as.integer(inference_fit$fit_harmonic_count)
      }
      if (!is.null(unified_inference)) inference_fit <- unified_inference
      minimum_pattern <- fitted_ploidy
      observable_pattern <- fitted_ploidy
      single_peak_ambiguity <- FALSE
    }

    # Known and inferred n share the preserved final fitting engine. Prefix
    # recovery changes the first-peak seed and fitting start only; the selected
    # common/trend/adaptive width model handles all later components.
    fit <- legacy_environment$.GenomeGMM_auto(
      data_file = data_file,
      n = fitted_ploidy,
      species_name = species_name,
      search = "fast",
      candidate_max_iter = 1500,
      candidate_tolerance = 1e-5,
      verification_top = 4,
      verification_max_iter = 10000,
      final_max_iter = 10000,
      save_search = FALSE,
      save_plot = TRUE,
      .first_peak_override = recovery$first_peak_x,
      .start_index_override = recovery$start_index,
      .prefix_completion_slot = recovery$selected_slot,
      .prefix_completion_weight = 0.35,
      .visible_peak_override = recovery$visible_peak_x
    )
  }
  if (is.null(fit)) {
    cat("\n[Preserved known-n GenomeGMM fit]\n")
    cat("Fitting the final curve with the preserved n =", fitted_ploidy,
        "engine from", attr(legacy_environment, "legacy_file"), "\n")
    fit <- legacy_environment$.GenomeGMM_auto(
      data_file = data_file,
      n = fitted_ploidy,
      species_name = species_name,
      search = "fast",
      candidate_max_iter = 1500,
      candidate_tolerance = 1e-5,
      verification_top = 4,
      verification_max_iter = 10000,
      final_max_iter = 10000,
      save_search = FALSE,
      save_plot = TRUE
    )
  }

  # Hidden-prefix completion and final high-iteration tuning can move the
  # component means and the red/grey integration switch enough that the
  # candidate selected during coarse screening no longer satisfies the
  # supplied Monoploid-size interval. Close that loop on the final curve.
  # Only the selected base and its immediate neighbours are reconsidered;
  # this is a consistency correction, not a second unrestricted ploidy scan.
  if (is.null(n) && !is.null(monoploid_range) && !is.null(fit)) {
    final_ordered_miu <- sort(as.numeric(
      fit$miu[is.finite(fit$miu) & fit$miu > 0]
    ))
    final_total_kmers <- suppressWarnings(as.numeric(fit$total_kmers))
    neighbour_slots <- sort(unique(as.integer(c(
      fitted_ploidy - 1L, fitted_ploidy, fitted_ploidy + 1L
    ))))
    neighbour_slots <- neighbour_slots[
      neighbour_slots >= 2L & neighbour_slots <= length(final_ordered_miu)
    ]
    neighbour_unit_size <- if (length(final_total_kmers) &&
        is.finite(final_total_kmers[1L]) && final_total_kmers[1L] > 0 &&
        length(neighbour_slots)) {
      final_total_kmers[1L] / final_ordered_miu[neighbour_slots]
    } else {
      rep(NA_real_, length(neighbour_slots))
    }
    lower_bp <- min(monoploid_range) * 1e6
    upper_bp <- max(monoploid_range) * 1e6
    neighbour_allopolyploid <- neighbour_unit_size * neighbour_slots
    neighbour_direct_match <- is.finite(neighbour_unit_size) &
      neighbour_unit_size >= lower_bp & neighbour_unit_size <= upper_bp
    neighbour_allopolyploid_match <- is.finite(neighbour_allopolyploid) &
      neighbour_allopolyploid >= lower_bp &
      neighbour_allopolyploid <= upper_bp
    neighbour_in_range <- neighbour_direct_match |
      neighbour_allopolyploid_match
    neighbour_monoploid <- ifelse(
      neighbour_allopolyploid_match & !neighbour_direct_match,
      neighbour_allopolyploid, neighbour_unit_size
    )
    current_index <- match(fitted_ploidy, neighbour_slots)
    current_in_range <- !is.na(current_index) &&
      isTRUE(neighbour_in_range[current_index])

    if (!isTRUE(current_in_range) && any(neighbour_in_range)) {
      eligible_slots <- neighbour_slots[neighbour_in_range]
      terminal_anchor <- suppressWarnings(as.integer(
        inference_fit$range_terminal_anchor_slot
      ))
      if (length(terminal_anchor) == 1L && is.finite(terminal_anchor)) {
        range_refit_n <- eligible_slots[which.min(abs(
          eligible_slots - terminal_anchor
        ))]
      } else {
        range_refit_n <- eligible_slots[which.min(abs(
          eligible_slots - fitted_ploidy
        ))]
      }

      if (is.finite(range_refit_n) && range_refit_n != fitted_ploidy) {
        old_range_pdf <- fit$pdf_file
        cat(
          "Final Monoploid-range consistency correction:",
          fitted_ploidy, "x ->", range_refit_n, "x\n"
        )
        cat(
          "  final candidate Monoploid sizes (Mb):",
          paste0(neighbour_slots, "x=",
                 sprintf("%.2f", neighbour_monoploid / 1e6),
                 collapse = ", "), "\n"
        )

        if (isTRUE(recovery$accepted)) {
          corrected_fit <- legacy_environment$.GenomeGMM_auto(
            data_file = data_file,
            n = range_refit_n,
            species_name = species_name,
            search = "fast",
            candidate_max_iter = 1500,
            candidate_tolerance = 1e-5,
            verification_top = 4,
            verification_max_iter = 10000,
            final_max_iter = 10000,
            save_search = FALSE,
            save_plot = TRUE,
            .first_peak_override = recovery$first_peak_x,
            .start_index_override = recovery$start_index,
            .prefix_completion_slot = recovery$selected_slot,
            .prefix_completion_weight = 0.35,
            .visible_peak_override = recovery$visible_peak_x
          )
        } else {
          corrected_args <- baseline_args
          corrected_args$n <- range_refit_n
          corrected_fit <- do.call(
            legacy_environment$.GenomeGMM_auto, corrected_args
          )
        }
        if (!is.null(old_range_pdf) && length(old_range_pdf) == 1L &&
            is.character(old_range_pdf) && !is.na(old_range_pdf) &&
            nzchar(old_range_pdf) && file.exists(old_range_pdf) &&
            !identical(old_range_pdf, corrected_fit$pdf_file)) {
          unlink(old_range_pdf)
        }
        fit <- corrected_fit
        fitted_ploidy <- as.integer(range_refit_n)
        minimum_pattern <- fitted_ploidy
        observable_pattern <- fitted_ploidy
        single_peak_ambiguity <- FALSE
        inference_fit$range_selected_ploidy <- fitted_ploidy
        corrected_index <- match(fitted_ploidy, neighbour_slots)
        inference_fit$range_selected_size_model <- if (
            !is.na(corrected_index) &&
            isTRUE(neighbour_allopolyploid_match[corrected_index]) &&
            !isTRUE(neighbour_direct_match[corrected_index])) {
          "balanced_allopolyploid"
        } else if (!is.na(corrected_index) &&
            isTRUE(neighbour_direct_match[corrected_index]) &&
            !isTRUE(neighbour_allopolyploid_match[corrected_index])) {
          "direct_or_autopolyploid"
        } else {
          "ambiguous"
        }
        inference_fit$range_selection_basis <-
          "final fitted-curve Monoploid-range consistency"
      }
    }
  }

  fit$hidden_prefix_recovery_triggered <- isTRUE(recovery$trigger)
  fit$hidden_prefix_recovery_accepted <- isTRUE(recovery$accepted)
  fit$hidden_prefix_start_x <- recovery$start_x
  fit$hidden_prefix_baseline_metrics <- recovery$baseline_metrics
  fit$hidden_prefix_candidate_table <- recovery$candidate_table
  fit$hidden_prefix_structure <- recovery$prefix_structure
  fit$hidden_prefix_broad_anchor <- recovery$broad_anchor
  fit$hidden_prefix_tail_support <- recovery$tail_support
  fit$hidden_prefix_tail_extension_table <- recovery$tail_extension_table
  fit$hidden_prefix_selected_fit_count <- fitted_ploidy
  fit$hidden_prefix_selected_slot <- recovery$selected_slot
  fit$hidden_prefix_plausible_slots <- recovery$plausible_slots
  # Preserve the historical filename when n was supplied explicitly.  The
  # estimated-ploidy filename belongs only to the without-n workflow.
  if (is.null(n)) {
    fit$pdf_file <- .GenomeGMM_rename_legacy_pdf(
      fit$pdf_file, data_file, fitted_ploidy
    )
  }
  fit$fit_harmonic_count <- fitted_ploidy
  fit$single_peak_ambiguity <- isTRUE(single_peak_ambiguity)
  fit$two_peak_shoulder_assessment <- two_peak_shoulder_assessment
  fit$minimum_pattern <- minimum_pattern
  fit$minimum_inferred_ploidy <- minimum_pattern
  fit$minimum_observable_pattern <- observable_pattern
  if (is.null(fit$g_hap_bp) || length(fit$g_hap_bp) != 1L ||
      !is.finite(fit$g_hap_bp) || fit$g_hap_bp <= 0) {
    fit$g_hap_bp <- fit$genome_size_bp / fitted_ploidy
  }
  fit$haploid_genome_size_bp <- fit$g_hap_bp
  fit$minimum_pattern_genome_size_bp <- minimum_pattern * fit$g_hap_bp
  if (is.null(fit$selected_parameters)) fit$selected_parameters <- list()
  model_fields <- c(
    "fit_extra_component", "sigma_mode", "position_mode",
    "use_calibration", "main_peak_weight", "symmetry_tolerance"
  )
  for (field in model_fields) {
    if (is.null(fit$selected_parameters[[field]]) &&
        !is.null(fit[[field]])) {
      fit$selected_parameters[[field]] <- fit[[field]]
    }
  }
  fit$selected_parameters$minimum_pattern <- minimum_pattern
  fit$selected_parameters$minimum_inferred_ploidy <- minimum_pattern
  fit$selected_parameters$minimum_observable_pattern <- observable_pattern
  fit$selected_parameters$fit_harmonic_count <- fitted_ploidy

  if (!is.null(inference_fit)) {
    inference_fields <- c(
      "minimum_pattern_candidates", "global_pattern_peak_count",
      "global_pattern_trusted_count", "global_pattern_peak_x",
      "global_pattern_ratios", "global_pattern_slots",
      "global_pattern_evidence", "terminal_harmonic_slot",
      "single_peak_elbow", "single_peak_elbow_applied",
      "single_peak_elbow_guard", "range_selected_ploidy",
      "range_terminal_anchor_slot", "range_selection_basis",
      "range_selected_size_model"
    )
    for (field in inference_fields) {
      fit[[field]] <- inference_fit[[field]]
    }
  } else {
    fit$minimum_pattern_candidates <- NULL
  }
  if (is.null(n)) {
    report <- .GenomeGMM_write_minimum_pattern_report(
      data_file = data_file, species_name = species_name, fit = fit,
      monoploid_range = monoploid_range
    )

    fit$minimum_pattern <- report$minimum_pattern
    fit$minimum_inferred_ploidy <- report$minimum_pattern
    fit$peak_supported_pattern <- report$peak_supported_pattern
    fit$ploidy_candidates <- report$ploidy_candidates
    fit$ploidy_candidate_matrix <- report$ploidy_candidate_matrix
    fit$compatible_ploidies <- report$compatible_ploidies
    fit$compatible_genome_sizes_bp <- report$compatible_genome_sizes_bp
    fit$compatible_g_hap_sizes_bp <- report$compatible_g_hap_sizes_bp
    fit$monoploid_range_mb <- report$monoploid_range_mb
    fit$range_matched_candidates <- report$range_matched_candidates
    fit$range_match_count <- report$range_match_count
    fit$compatible_compositions <- report$compositions
    fit$g_hap_bp <- report$g_hap_bp
    fit$haploid_genome_size_bp <- report$g_hap_bp
    fit$minimum_pattern_genome_size_bp <-
      report$minimum_pattern_genome_size_bp
    fit$selected_parameters$minimum_pattern <- report$minimum_pattern
    fit$selected_parameters$minimum_inferred_ploidy <-
      report$minimum_pattern
    fit$homozygous_peak_depth <- report$homozygous_peak_depth
    fit$pattern_confidence <- report$pattern_confidence
    fit$pattern_assessment <- report$pattern_assessment
    fit$fit_quality <- report$fit_quality
    fit$absolute_ploidy_identifiable <- FALSE
    fit$result_file <- report$file

    cat("\n[Minimum inferred ploidy interpretation]\n")
    cat("Minimum inferred ploidy:",
        paste0(fit$minimum_pattern, "x"), "\n")
    cat("Compatible ploidies:",
        paste0(fit$compatible_ploidies, "x", collapse = ", "), "...\n")
    cat("Minimum-ploidy Holoploid genome size:",
        sprintf("%.2f Mb", fit$minimum_pattern_genome_size_bp / 1e6), "\n")
    cat("Monoploid genome size:",
        sprintf("%.2f Mb", fit$g_hap_bp / 1e6), "\n")
    cat("Absolute biological ploidy from histo alone: not identifiable\n")
    cat("Result report saved as:", fit$result_file, "\n")
  } else {
    # The internal fixed-pattern fitter populates its return metadata here.
    # The public interface may subsequently add an allopolyploid explanation
    # TXT after it restores the user-supplied biological ploidy.
    fit$compatible_ploidies <- fitted_ploidy
    fit$compatible_genome_sizes_bp <- fit$genome_size_bp
    fit$g_hap_bp <- if (!is.null(fit$g_hap_bp) &&
        length(fit$g_hap_bp) == 1L && is.finite(fit$g_hap_bp)) {
      fit$g_hap_bp
    } else {
      fit$genome_size_bp / fitted_ploidy
    }
    fit$compatible_g_hap_sizes_bp <- fit$g_hap_bp
    fit$compatible_compositions <- NA_character_
    fit$haploid_genome_size_bp <- fit$g_hap_bp
    fit$pattern_confidence <- "Supplied"
    fit$pattern_assessment <- NULL
    fit$fit_quality <- .GenomeGMM_fit_quality(fit)
    fit$absolute_ploidy_identifiable <- TRUE
    fit$result_file <- NULL
  }

  invisible(fit)
}

.GenomeGMM_print_summary <- function(fit, known_n = NULL, data_file = NULL,
                                     species_name = "") {
  finite_scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
  }
  character_scalar <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
  }
  selected <- fit[["selected_parameters", exact = TRUE]]
  if (is.null(selected)) selected <- list()
  selected_value <- function(name) selected[[name, exact = TRUE]]
  print_numeric_vector <- function(label, value, digits = 4L) {
    value <- as.numeric(value)
    value <- value[is.finite(value)]
    if (length(value)) {
      cat(label, paste(round(value, digits), collapse = ", "), "\n")
    }
  }

  cat("\n=== GenomeGMM final result ===\n")
  if (character_scalar(data_file)) cat("Input:", data_file, "\n")
  if (character_scalar(species_name)) cat("Species:", species_name, "\n")
  if (is.null(known_n)) {
    cat("Mode: automatic minimum-ploidy inference\n")
    if (finite_scalar(fit$minimum_pattern)) {
      cat("Minimum inferred ploidy:",
          paste0(fit$minimum_pattern, "x"), "\n")
    }
    if (character_scalar(fit$pattern_confidence)) {
      cat("Pattern confidence:", fit$pattern_confidence, "\n")
    }
    if (finite_scalar(fit$fit_harmonic_count)) {
      cat("Selected n for fitted curve:", fit$fit_harmonic_count, "\n")
    }
  } else {
    cat("Mode: supplied biological ploidy",
        paste0(as.integer(known_n), "x"), "\n")
    if (character_scalar(fit$ploidy_type)) {
      cat("Ploidy type:", fit$ploidy_type, "\n")
    }
    if (finite_scalar(fit$base_pattern)) {
      cat("Selected base pattern for curve fitting:",
          paste0(as.integer(fit$base_pattern), "x"), "\n")
    }
    if (character_scalar(fit$known_ploidy_warning)) {
      cat("Interpretation warning:", fit$known_ploidy_warning, "\n")
    }
  }
  cat("\n[Selected model]\n")
  if (finite_scalar(fit$fit_harmonic_count)) {
    cat("Primary harmonic count:", fit$fit_harmonic_count, "\n")
  }
  if (length(fit$sigma)) cat("Total fitted components:", length(fit$sigma), "\n")
  extra_component <- selected_value("fit_extra_component")
  if (is.logical(extra_component) && length(extra_component) == 1L &&
      !is.na(extra_component)) {
    cat("Optional tail component:", if (extra_component) "used" else "not used", "\n")
  }
  if (character_scalar(selected_value("sigma_mode"))) {
    cat("Sigma model:", selected_value("sigma_mode"), "\n")
  }
  if (!is.null(fit$final_shape_refinement) &&
      isTRUE(fit$final_shape_refinement$triggered)) {
    cat("Final guarded shape refinement:",
        if (isTRUE(fit$final_shape_refinement$accepted)) {
          "accepted"
        } else {
          "kept selected parameters"
        }, "\n")
  }
  if (character_scalar(selected_value("position_mode"))) {
    cat("Peak-position model:", selected_value("position_mode"), "\n")
  }
  calibration <- selected_value("use_calibration")
  if (is.logical(calibration) && length(calibration) == 1L && !is.na(calibration)) {
    cat("Post-fit calibration:", if (calibration) "used" else "not used", "\n")
  }
  if (finite_scalar(selected_value("main_peak_weight"))) {
    cat("Main-peak weight:", sprintf("%.2f", selected_value("main_peak_weight")), "\n")
  }
  if (finite_scalar(selected_value("symmetry_tolerance"))) {
    cat("Symmetry tolerance:",
        sprintf("%.2f", selected_value("symmetry_tolerance")), "\n")
  }
  if (character_scalar(fit$joint_profile)) {
    joint_status <- if (isTRUE(fit$joint_accepted)) "accepted" else "not accepted"
    cat("Joint optimization:", fit$joint_profile, paste0("(", joint_status, ")"), "\n")
  }
  if (!is.null(fit$prefix_residual_refinement)) {
    cat("Recovered-prefix residual refinement:",
        if (isTRUE(fit$prefix_residual_refinement$accepted)) "accepted" else "kept baseline",
        "\n")
    if (isTRUE(fit$prefix_residual_refinement$accepted)) {
      cat("First-peak mean shift:",
          sprintf("%.4f", fit$prefix_residual_refinement$first_mean_shift), "\n")
      cat("Later peak widths: nondecreasing, data-selected smoothness\n")
    }
  }

  cat("\n[Final EM parameters]\n")
  if (is.logical(fit$converged) && length(fit$converged) == 1L &&
      !is.na(fit$converged)) {
    cat("Converged:", if (fit$converged) "yes" else "no", "\n")
  }
  if (finite_scalar(fit$iterations)) cat("Iterations:", fit$iterations, "\n")
  print_numeric_vector("Mixing coefficients alpha:", fit$alpha, 4L)
  print_numeric_vector("Means miu:", fit$miu, 4L)
  print_numeric_vector("Standard deviations sigma:", fit$sigma, 4L)
  if (length(fit$component_sizes_bp)) {
    component_sizes_mb <- as.numeric(fit$component_sizes_bp) / 1e6
    component_sizes_mb <- component_sizes_mb[is.finite(component_sizes_mb)]
    label <- if (length(component_sizes_mb) > length(fit$alpha)) {
      "Component sizes, including residual (Mb):"
    } else {
      "Component sizes (Mb):"
    }
    print_numeric_vector(label, component_sizes_mb, 2L)
  }
  if (finite_scalar(fit$fit_score)) {
    cat("Automatic fit score:", sprintf("%.6f", fit$fit_score), "\n")
  }
  if (!is.null(fit$fit_metrics) && finite_scalar(fit$fit_metrics$nrmse)) {
    cat("NRMSE:", sprintf("%.6f", fit$fit_metrics$nrmse), "\n")
  }
  if (!is.null(fit$fit_metrics) &&
      finite_scalar(fit$fit_metrics$main_peak_nrmse)) {
    cat("Main-peak NRMSE:",
        sprintf("%.6f", fit$fit_metrics$main_peak_nrmse), "\n")
  }

  cat("\n[Genome-size estimation]\n")
  if (finite_scalar(fit$total_kmers)) {
    cat("Total k-mer counts:",
        format(fit$total_kmers, scientific = TRUE, digits = 7L), "\n")
  }
  if (finite_scalar(fit$total_kmer_switch_x)) {
    cat("Red/grey switch coverage:",
        sprintf("%.2f", fit$total_kmer_switch_x), "\n")
  }
  if (finite_scalar(fit$homozygous_peak_depth)) {
    cat("Homozygous peak depth:",
        sprintf("%.2f", fit$homozygous_peak_depth), "\n")
  }
  if (finite_scalar(fit$g_hap_bp)) {
    cat("Monoploid genome size:",
        sprintf("%.2f Mb", fit$g_hap_bp / 1e6), "\n")
  }
  full_size <- if (is.null(known_n)) {
    fit$minimum_pattern_genome_size_bp
  } else {
    fit$genome_size_bp
  }
  if (finite_scalar(full_size)) {
    label <- if (is.null(known_n)) {
      "Minimum-inferred-ploidy Holoploid genome size:"
    } else {
      "Holoploid genome size:"
    }
    cat(label, sprintf("%.2f Mb", full_size / 1e6), "\n")
  }
  cat("\n[Output files]\n")
  if (!is.null(fit$pdf_file) && length(fit$pdf_file) == 1L &&
      !is.na(fit$pdf_file) && nzchar(fit$pdf_file)) {
    cat("PDF:", fit$pdf_file, "\n")
  }
  if (!is.null(fit$result_file) && length(fit$result_file) == 1L &&
      !is.na(fit$result_file) && nzchar(fit$result_file)) {
    cat("TXT:", fit$result_file, "\n")
  }
}

# Public single-command interface. Candidate-search and EM messages are
# suppressed so batch logs contain only final results and output paths.
GenomeGMM <- function(data_file, n = NULL, species_name = "",
                      monoploid_range = NULL, ploidy_type = NULL) {
  if (!is.character(data_file) || length(data_file) != 1L ||
      is.na(data_file) || !nzchar(data_file) || !file.exists(data_file)) {
    stop("Histogram file does not exist: ", data_file)
  }
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n < 2L || n != round(n)) {
      stop("Parameter 'n' must be the supplied biological ploidy as an integer >= 2.")
    }
    if (is.null(ploidy_type) || !is.character(ploidy_type) ||
        length(ploidy_type) != 1L || is.na(ploidy_type) ||
        !ploidy_type %in% c("autopolyploid", "allopolyploid")) {
      stop(
        "When n is supplied, ploidy_type must be explicitly set to ",
        "'autopolyploid' or 'allopolyploid'."
      )
    }
    n <- as.integer(n)
  } else if (!is.null(ploidy_type)) {
    stop("ploidy_type is only used together with a supplied biological ploidy n.")
  }
  if (!is.null(monoploid_range)) {
    if (!is.null(n)) {
      stop("monoploid_range can only be used when n is not supplied.")
    }
    monoploid_range <- suppressWarnings(as.numeric(monoploid_range))
    if (length(monoploid_range) != 2L ||
        any(!is.finite(monoploid_range)) ||
        any(monoploid_range <= 0)) {
      stop("monoploid_range must contain two positive finite values: c(a, b).")
    }
    monoploid_range <- sort(monoploid_range)
    if (monoploid_range[1L] == monoploid_range[2L]) {
      stop("monoploid_range must have a lower bound smaller than its upper bound.")
    }
  }
  actual_ploidy <- n
  resolution <- NULL
  internal_n <- n
  if (!is.null(n)) {
    resolution_output <- utils::capture.output(
      resolution <- .GenomeGMM_resolve_known_pattern(
        data_file = data_file, actual_ploidy = n,
        ploidy_type = ploidy_type, species_name = species_name
      ),
      type = "output"
    )
    rm(resolution_output)
    if (!resolution$status %in% c("ok", "fallback_direct")) {
      report <- .GenomeGMM_write_unfactorable_allopoly_report(
        data_file = data_file, actual_ploidy = n,
        species_name = species_name
      )
      cat("\n=== GenomeGMM result ===\n")
      cat("Input:", data_file, "\n")
      if (nzchar(species_name)) cat("Species:", species_name, "\n")
      cat("Status: no valid balanced allopolyploid base pattern\n")
      cat(resolution$selection_basis, "\n")
      cat("PDF: not generated\n")
      cat("TXT:", report$result_file, "\n")
      return(invisible(report))
    }
    internal_n <- as.integer(resolution$base_pattern)
  }
  peak_error <- NULL
  tryCatch(
    .GenomeGMM_assert_reliable_peak(data_file),
    error = function(error) peak_error <<- conditionMessage(error)
  )
  if (!is.null(peak_error)) {
    no_peak_result <- .GenomeGMM_write_no_peak_report(data_file, n = n)
    cat("\n=== GenomeGMM result ===\n")
    cat("Input:", data_file, "\n")
    if (nzchar(species_name)) cat("Species:", species_name, "\n")
    cat("Status: no reliable genomic peak detected\n")
    cat(no_peak_result$message, "\n")
    cat("TXT:", no_peak_result$result_file, "\n")
    return(invisible(no_peak_result))
  }
  fit <- NULL
  previous_reporting_ploidy <- getOption("GenomeGMM.reporting_ploidy")
  previous_monoploid_divisor <- getOption(
    "GenomeGMM.reporting_monoploid_divisor"
  )
  if (!is.null(actual_ploidy)) {
    options(GenomeGMM.reporting_ploidy = as.integer(actual_ploidy))
    reporting_monoploid_divisor <- if (identical(
        ploidy_type, "allopolyploid"
      )) {
      as.integer(actual_ploidy %/% internal_n)
    } else {
      as.integer(actual_ploidy)
    }
    options(
      GenomeGMM.reporting_monoploid_divisor =
        reporting_monoploid_divisor
    )
  }
  on.exit(
    options(
      GenomeGMM.reporting_ploidy = previous_reporting_ploidy,
      GenomeGMM.reporting_monoploid_divisor =
        previous_monoploid_divisor
    ),
    add = TRUE
  )
  discarded_output <- utils::capture.output(
    fit <- .GenomeGMM_run(
      data_file = data_file, n = internal_n, species_name = species_name,
      monoploid_range = monoploid_range
    ),
    type = "output"
  )
  rm(discarded_output)
  if (!is.null(actual_ploidy)) {
    fit$actual_ploidy <- as.integer(actual_ploidy)
    fit$ploidy_type <- ploidy_type
    fit$base_pattern <- as.integer(internal_n)
    fit$base_pattern_multiplier <- as.integer(actual_ploidy %/% internal_n)
    fit$known_ploidy_resolution <- resolution
    fit$known_ploidy_warning <- resolution$warning
    fit$subgenome_unit_size_bp <- as.numeric(fit$g_hap_bp)
    fit$g_hap_bp <- if (identical(ploidy_type, "allopolyploid")) {
      fit$subgenome_unit_size_bp * internal_n
    } else {
      fit$subgenome_unit_size_bp
    }
    fit$genome_size_bp <- fit$subgenome_unit_size_bp * actual_ploidy
    fit$haploid_genome_size_bp <- fit$g_hap_bp
    fit$compatible_genome_sizes_bp <- fit$genome_size_bp
    fit$compatible_g_hap_sizes_bp <- fit$g_hap_bp
    fit$minimum_pattern_genome_size_bp <-
      fit$subgenome_unit_size_bp * internal_n
    fit$pdf_file <- .GenomeGMM_rename_known_pdf(
      fit$pdf_file, data_file, actual_ploidy
    )
    if (identical(ploidy_type, "allopolyploid")) {
      fit$result_file <- .GenomeGMM_write_known_allopoly_report(
        data_file = data_file, species_name = species_name,
        fit = fit, resolution = resolution
      )
    }
  }
  .GenomeGMM_print_summary(
    fit, known_n = actual_ploidy, data_file = data_file,
    species_name = species_name
  )
  invisible(fit)
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
    "Depth: %.2fx\nHoloploid GSE: %.1f Mb\nMonoploid GSE: %.1f Mb\nR2: %.4f",
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
    cat("  holoploid genome size:", round(genome$genome_size_bp / 1e6, 2), "Mb\n")
    cat("  monoploid genome size:", round(genome$genome_size_bp / n / 1e6, 2), "Mb\n")
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
