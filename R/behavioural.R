# Raw E-Prime files -> behav_data

# Trials ----
# Different exclusion list than eye-tracking (Tobii files not used)
build_behav_trials <- function(cfg = load_config()) {
  eprime_dir <- here::here(cfg$raw$eprime_clean)
  excl_ids <- unlist(cfg$exclusions$behavioural, use.names = FALSE)

  eprime_files <- find_files(eprime_dir, eprime_re) |>
    dplyr::filter(!participant_id %in% excl_ids) |>
    dplyr::arrange(as.integer(stringr::str_remove(participant_id, "^P")), change_type)

  n_files <- table(eprime_files$participant_id)
  if (any(n_files != 2L)) {
    stop("change_type file missing: ", toString(names(n_files)[n_files != 2L]))
  }

  trials <- purrr::map(eprime_files$path, read_eprime) |>
    purrr::list_rbind() |>
    dplyr::mutate(object = extract_object(target_media))

  stopifnot(
    !anyNA(trials$object),
    !anyNA(trials$acc),
    all(trials$position %in% c("upright", "inverted"))
  )
  trials
}

# Subject means ----
aggregate_behavioural <- function(trials) {
  trials |>
    dplyr::group_by(participant_id, object, position, change_type) |>
    dplyr::summarise(
      accuracy = mean(acc == 1L),
      rt = mean(rt[acc == 1L]),  # correct trials only
      n_trials = dplyr::n(),
      n_correct = sum(acc == 1L),
      .groups = "drop"
    ) |>
    # Inverse efficiency = rt / accuracy; no correct trials = NA
    dplyr::mutate(ies = dplyr::if_else(accuracy == 0, NA_real_, rt / accuracy))
}
