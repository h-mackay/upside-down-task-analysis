# trial_data -> et_data (subject means; Top proportions); object collapse; error bar SE

# Subject means ----
aggregate_eyetracking <- function(trial_data) {
  subject_aoi <- trial_data |>
    dplyr::group_by(participant_id, change_type, phase, position, object, aoi) |>
    # Later lines use averaged columns not raw (order matters); counts + log first
    dplyr::summarise(
      n_variants = dplyr::n(),
      n_valid_tff = sum(!is.na(tff_ms)),
      log_tff_ms = mean_valid(log(pmax(tff_ms, 1))),  # floor = 1 ms: no log(0)
      tff_ms = mean_valid(tff_ms),
      n_fix = mean(n_fix),
      n_sacc = mean(n_sacc),
      dwell_ms = mean(dwell_ms),
      fix_dur_ms = dplyr::if_else(  # ratio of means above
        n_fix == 0,
        NA_real_,
        dwell_ms / n_fix
      ),
      .groups = "drop"
    )

  # Top proportions: per stimulus, then averaged (stimuli weighted equally)
  aoi_wide <- trial_data |>
    dplyr::select(
      participant_id,
      change_type,
      phase,
      position,
      object,
      media,
      aoi,
      n_fix,
      dwell_ms,
      n_sacc
    ) |>
    tidyr::pivot_wider(
      names_from = aoi,
      values_from = c(n_fix, dwell_ms, n_sacc)
    )
  stopifnot("missing aoi row" = !anyNA(aoi_wide))

  trial_props <- aoi_wide |>
    dplyr::mutate(
      total_fix = n_fix_top + n_fix_bottom,
      total_dwell = dwell_ms_top + dwell_ms_bottom,
      total_sacc = n_sacc_top + n_sacc_bottom,
      prop_top_fix = dplyr::if_else(total_fix == 0, NA_real_, n_fix_top / total_fix),
      prop_top_dwell = dplyr::if_else(total_dwell == 0, NA_real_, dwell_ms_top / total_dwell),
      prop_top_sacc = dplyr::if_else(total_sacc == 0, NA_real_, n_sacc_top / total_sacc)
    )

  subject_prop <- trial_props |>
    dplyr::group_by(participant_id, change_type, phase, position, object) |>
    # Counts first (same as above)
    dplyr::summarise(
      n_variants = dplyr::n(),
      n_valid_prop_fix = sum(!is.na(prop_top_fix)),
      n_valid_prop_dwell = sum(!is.na(prop_top_dwell)),
      n_valid_prop_sacc = sum(!is.na(prop_top_sacc)),
      prop_top_fix = mean_valid(prop_top_fix),
      prop_top_dwell = mean_valid(prop_top_dwell),
      prop_top_sacc = mean_valid(prop_top_sacc),
      .groups = "drop"
    )
  list(aoi = subject_aoi, prop = subject_prop)
}

# Object collapse ----
collapse_objects <- function(aoi_obj, prop_obj) {
  collapsed_aoi <- aoi_obj |>
    dplyr::group_by(participant_id, change_type, phase, position, aoi) |>
    # Counts first; fix_dur_ms from collapsed means
    dplyr::summarise(
      n_obj = dplyr::n(),
      n_obj_valid_tff = sum(!is.na(tff_ms)),
      n_obj_valid_fix_dur = sum(!is.na(fix_dur_ms)),
      n_fix = mean(n_fix),
      n_sacc = mean(n_sacc),
      dwell_ms = mean(dwell_ms),
      fix_dur_ms = dplyr::if_else(
        n_fix == 0,
        NA_real_,
        dwell_ms / n_fix
      ),
      tff_ms = mean_valid(tff_ms),
      log_tff_ms = mean_valid(log_tff_ms),
      .groups = "drop"
    )

  collapsed_prop <- prop_obj |>
    dplyr::group_by(participant_id, change_type, phase, position) |>
    dplyr::summarise(
      n_obj = dplyr::n(),
      n_obj_valid_fix = sum(!is.na(prop_top_fix)),
      n_obj_valid_dwell = sum(!is.na(prop_top_dwell)),
      n_obj_valid_sacc = sum(!is.na(prop_top_sacc)),
      prop_top_fix = mean_valid(prop_top_fix),
      prop_top_dwell = mean_valid(prop_top_dwell),
      prop_top_sacc = mean_valid(prop_top_sacc),
      .groups = "drop"
    )

  list(aoi = collapsed_aoi, prop = collapsed_prop)
}

# Error bars ----
# Within-subject SE (Cousineau-Morey): remove each subject's mean; add back grand mean;
# correct for number of conditions (Morey)
cousineau_morey_se <- function(data, dv, ..., subject_var = "participant_id") {
  grouping_vars <- rlang::enquos(...)

  subject_means <- data |>
    dplyr::group_by(.data[[subject_var]]) |>
    dplyr::summarise(subject_mean = mean(.data[[dv]], na.rm = TRUE), .groups = "drop")

  grand_mean <- mean(data[[dv]], na.rm = TRUE)

  normalized_data <- data |>
    dplyr::left_join(subject_means, by = subject_var) |>
    dplyr::mutate(normalized_dv = .data[[dv]] - subject_mean + grand_mean)

  k <- data |>  # conditions per subject (averaged if not equal)
    dplyr::group_by(.data[[subject_var]]) |>
    dplyr::summarise(n_conditions = sum(!is.na(.data[[dv]])), .groups = "drop") |>
    dplyr::summarise(k = mean(n_conditions, na.rm = TRUE)) |>
    dplyr::pull(k)

  morey_factor <- if (!is.finite(k) || k <= 1) NA_real_ else sqrt(k / (k - 1))
  if (is.na(morey_factor)) {
    warning("< 2 conditions per subject; se_within NA")
  }

  normalized_data |>
    dplyr::group_by(!!!grouping_vars) |>
    dplyr::summarise(
      mean = mean(.data[[dv]], na.rm = TRUE),
      n_valid = sum(!is.na(normalized_dv)),
      se_within = if (is.na(morey_factor) || n_valid <= 1) {
        NA_real_
      } else {
        (sd(normalized_dv, na.rm = TRUE) / sqrt(n_valid)) * morey_factor
      },
      .groups = "drop"
    )
}
