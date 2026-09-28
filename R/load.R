# Raw Tobii + E-Prime files -> trial_data

# File names ----
# 1 = participant; 2 = change type
tobii_re <- "^(P[0-9]+)_(Orientation|Saturation)_TB_processed\\.xlsx$"
eprime_re <- "^(P[0-9]+)_(Orientation|Saturation)_TB_cleaned\\.xlsx$"

# Readers ----
# One row per presentation x AOI; all four sheets combined
read_tobii <- function(file_path) {
  fname <- basename(file_path)

  participant_id <- stringr::str_extract(fname, "^P[0-9]+")
  change_type <- dplyr::case_when(
    stringr::str_detect(fname, "(?i)orientation") ~ "orientation",
    stringr::str_detect(fname, "(?i)saturation") ~ "saturation"
  )

  # Sheet names: rsu = right side up; usd = upside down
  sheet_map <- tibble::tribble(
    ~sheet, ~phase, ~position,
    "Target_rsu", "target", "upright",
    "Target_usd", "target", "inverted",
    "Test_rsu", "test", "upright",
    "Test_usd", "test", "inverted"
  )

  # Rename: new = old
  tobii_cols <- c(
    media = "Media",
    aoi = "AOI",
    interval = "Interval",
    start_ms = "Start_of_interval",
    dur_interval_ms = "Duration_of_interval",
    tff_ms = "Time_to_first_fixation",
    n_fix = "Number_of_fixations",
    dwell_ms = "Total_duration_of_fixations",
    n_sacc = "Number_of_saccades_in_AOI"
  )

  read_sheet <- function(sheet, phase, position) {
    raw <- readxl::read_excel(file_path, sheet = sheet, .name_repair = "minimal")
    raw |>
      dplyr::select(dplyr::all_of(tobii_cols)) |>
      dplyr::mutate(phase = phase, position = position)
  }

  purrr::pmap(sheet_map, read_sheet) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      participant_id = participant_id,
      change_type = change_type,
      media = normalize_media_code(media),
      object = extract_object(media),
      aoi = stringr::str_to_lower(aoi),
      dplyr::across(
        c(interval, start_ms, dur_interval_ms, tff_ms, n_fix, dwell_ms, n_sacc),
        as.numeric
      )
    ) |>
    dplyr::relocate(
      participant_id,
      change_type,
      phase,
      position,
      media,
      object,
      aoi,
      interval,
      start_ms
    )
}

# One row per trial
read_eprime <- function(file_path) {
  fname <- basename(file_path)

  participant_id <- stringr::str_extract(fname, "^P[0-9]+")
  change_type <- dplyr::case_when(
    stringr::str_detect(fname, "(?i)orientation") ~ "orientation",
    stringr::str_detect(fname, "(?i)saturation") ~ "saturation"
  )

  raw <- readxl::read_excel(file_path, sheet = 1, .name_repair = "minimal")

  col_test <- grep("^(Orientation|Saturation)TestImg$", names(raw), value = TRUE)

  prefix <- stringr::str_remove(col_test, "TestImg$")
  if (stringr::str_to_lower(prefix) != change_type) {
    stop("no column prefix match: ", fname)
  }
  col_target <- paste0(prefix, "TargetImg")
  col_acc <- paste0(col_test, ".ACC")
  col_rt <- paste0(col_test, ".RT")

  raw |>
    dplyr::transmute(
      participant_id = participant_id,
      change_type = change_type,
      position = stringr::str_to_lower(.data[["Position[SubTrial]"]]),
      trial_index = dplyr::row_number(),
      target_media = normalize_media_code(.data[[col_target]]),
      test_media = normalize_media_code(.data[[col_test]]),
      acc = as.integer(.data[[col_acc]]),
      rt = as.numeric(.data[[col_rt]])
    )
}

# Matching and collapsing ----
# Matched by presentation order; stops if sequences differ
match_tobii_eprime <- function(tobii_long, eprime_trials) {
  pid <- unique(tobii_long$participant_id)
  change_type_name <- unique(tobii_long$change_type)
  eprime_ordered <- dplyr::arrange(eprime_trials, trial_index)

  match_phase <- function(phase_label, eprime_media_col) {
    tobii_long_phase <- dplyr::filter(tobii_long, phase == phase_label)

    presentations <- tobii_long_phase |>
      dplyr::distinct(media, interval, start_ms) |>
      dplyr::arrange(start_ms)

    tobii_seq <- presentations$media
    eprime_seq <- eprime_ordered[[eprime_media_col]]

    if (!identical(tobii_seq, eprime_seq)) {
      stop(pid, " ", change_type_name, " ", phase_label, " no matching sequences")
    }

    presentations$trial_index <- eprime_ordered$trial_index
    tobii_long_phase |>
      dplyr::left_join(presentations, by = c("media", "interval", "start_ms")) |>
      dplyr::left_join(
        dplyr::select(eprime_ordered, trial_index, acc, rt),
        by = "trial_index"
      )
  }

  dplyr::bind_rows(
    match_phase("target", "target_media"),
    match_phase("test", "test_media")
  )
}

# Correct trials only; average repeat presentations of each stimulus
collapse_correct <- function(matched) {
  matched |>
    dplyr::filter(acc == 1) |>
    dplyr::mutate(censored = is.na(tff_ms)) |>  # aoi never fixated
    dplyr::group_by(participant_id, change_type, phase, position, object, media, aoi) |>
    dplyr::summarise(
      n_fix = mean(n_fix),
      n_sacc = mean(n_sacc),
      dwell_ms = mean(dwell_ms),
      tff_ms = mean_valid(tff_ms),
      n_correct = dplyr::n(),
      n_censored = sum(censored),
      .groups = "drop"
    )
}

# Move misfiles listed in config from target back to test
repair_misfiled <- function(tobii_long, cfg) {
  mf <- cfg$boundary_misfiles
  if (is.null(mf) || length(mf) == 0L) {
    return(tobii_long)
  }

  pid <- unique(tobii_long$participant_id)
  change_type_name <- unique(tobii_long$change_type)

  # Config entries -> one row each
  mf_df <- purrr::map(mf, function(e) tibble::tibble(
    participant_id = e$participant_id,
    change_type = e$change_type,
    media = e$media,
    start_ms_int = as.integer(round(as.numeric(e$start_ms)))
  )) |>
    purrr::list_rbind()
  cell_misfiles <- dplyr::filter(
    mf_df,
    participant_id == pid,
    change_type == change_type_name
  )
  if (nrow(cell_misfiles) == 0L) {
    return(tobii_long)
  }

  tobii_long <- dplyr::mutate(tobii_long, start_ms_int = as.integer(round(start_ms)))

  found_misfiles <- tobii_long |>
    dplyr::filter(phase == "target") |>
    dplyr::semi_join(cell_misfiles, by = c("media", "start_ms_int")) |>
    dplyr::distinct(media, start_ms_int)
  if (nrow(found_misfiles) != nrow(cell_misfiles)) {
    stop(pid, " ", change_type_name, " misfile not found")
  }

  fix <- dplyr::mutate(cell_misfiles, .fix = TRUE)
  tobii_long |>
    dplyr::left_join(
      dplyr::select(fix, media, start_ms_int, .fix),
      by = c("media", "start_ms_int")
    ) |>
    dplyr::mutate(
      phase = dplyr::if_else(
        phase == "target" & dplyr::coalesce(.fix, FALSE),
        "test",
        phase
      )
    ) |>
    dplyr::select(-start_ms_int, -.fix)
}

# Files ----
# Participant and change type from file name
find_files <- function(dir, re) {
  fnames <- list.files(dir, pattern = re)
  if (length(fnames) == 0L) {
    stop("no matching files: ", dir)
  }
  m <- stringr::str_match(fnames, re)
  tibble::tibble(
    participant_id = m[, 2],
    change_type = stringr::str_to_lower(m[, 3]),
    path = file.path(dir, fnames)
  )
}

# Pair Tobii and E-Prime files; apply eye-tracking exclusions
pair_files <- function(cfg = load_config()) {
  tobii_dir <- here::here(cfg$raw$et_processed)
  eprime_dir <- here::here(cfg$raw$eprime_clean)
  tobii_files <- find_files(tobii_dir, tobii_re)
  eprime_files <- find_files(eprime_dir, eprime_re)

  excl_ids <- unlist(cfg$exclusions$eyetracking, use.names = FALSE)
  tobii_files <- dplyr::filter(tobii_files, !participant_id %in% excl_ids)
  eprime_files <- dplyr::filter(eprime_files, !participant_id %in% excl_ids)

  file_pairs <- dplyr::full_join(
    dplyr::rename(tobii_files, tobii_path = path),
    dplyr::rename(eprime_files, eprime_path = path),
    by = c("participant_id", "change_type")
  )
  unpaired <- dplyr::filter(file_pairs, is.na(tobii_path) | is.na(eprime_path))
  if (nrow(unpaired) > 0L) {
    stop(
      "unpaired files: ",
      toString(paste(unpaired$participant_id, unpaired$change_type))
    )
  }

  file_pairs |>
    dplyr::arrange(as.integer(stringr::str_remove(participant_id, "^P")), change_type)
}

# Main ----
build_trial_data <- function(cfg = load_config()) {
  file_pairs <- pair_files(cfg)

  # One cell = one participant x change type (one file pair)
  build_one_cell <- function(participant_id, change_type, tobii_path, eprime_path) {
    tobii <- read_tobii(tobii_path)
    tobii <- repair_misfiled(tobii, cfg)
    eprime <- read_eprime(eprime_path)
    matched <- match_tobii_eprime(tobii, eprime)
    collapse_correct(matched)
  }

  trial_data <- file_pairs |>
    purrr::pmap(build_one_cell) |>
    purrr::list_rbind()

  built_cells <- dplyr::distinct(trial_data, participant_id, change_type)
  missing_cells <- dplyr::anti_join(
    dplyr::select(file_pairs, participant_id, change_type),
    built_cells,
    by = c("participant_id", "change_type")
  )
  if (nrow(missing_cells) > 0L) {
    warning(
      "no rows after collapsing: ",
      toString(paste(missing_cells$participant_id, missing_cells$change_type))
    )
  }
  trial_data
}

# Diagnostics ----
# Run manually: presentations vs E-Prime trials per cell
audit_alignment <- function(cfg = load_config()) {
  file_pairs <- pair_files(cfg)

  audit_one_cell <- function(participant_id, change_type, tobii_path, eprime_path) {
    default_cell <- tibble::tibble(
      participant_id = participant_id,
      change_type = change_type,
      tobii_target = NA_integer_,
      tobii_test = NA_integer_,
      eprime_n = NA_integer_,
      diff_target = NA_integer_,
      diff_test = NA_integer_,
      status = NA_character_,
      error = NA_character_
    )

    # Unreadable file becomes error row instead of stopping audit
    tryCatch({
      tobii <- read_tobii(tobii_path)
      eprime <- read_eprime(eprime_path)

      n_presentations <- function(current_phase) {
        tobii |>
          dplyr::filter(phase == current_phase) |>
          dplyr::distinct(media, interval, start_ms) |>
          nrow()
      }
      tobii_target <- n_presentations("target")
      tobii_test <- n_presentations("test")
      eprime_n <- nrow(eprime)
      diff_target <- tobii_target - eprime_n
      diff_test <- tobii_test - eprime_n

      status <- dplyr::case_when(
        diff_target == 0 & diff_test == 0 ~ "ok",
        diff_target > 0 & diff_target == -diff_test ~ "misfile",
        .default = "other"
      )

      potential_misfiles <- tibble::tibble()
      if (diff_target > 0) {
        tobii_target_counts <- tobii |>
          dplyr::filter(phase == "target") |>
          dplyr::distinct(media, interval, start_ms) |>
          dplyr::count(media, name = "n_tobii")
        eprime_target_counts <- eprime |>
          dplyr::count(media = target_media, name = "n_eprime")
        extra_media <- tobii_target_counts |>
          dplyr::full_join(eprime_target_counts, by = "media") |>
          dplyr::mutate(
            excess = dplyr::coalesce(n_tobii, 0L) - dplyr::coalesce(n_eprime, 0L)
          ) |>
          dplyr::filter(excess > 0) |>
          dplyr::pull(media)
        potential_misfiles <- eprime |>
          dplyr::filter(test_media %in% extra_media) |>
          dplyr::transmute(
            participant_id,
            change_type,
            trial_index,
            position,
            target_media,
            test_media,
            acc,
            rt,
            rt_in_band = rt >= cfg$timing$targ_dur_min & rt <= cfg$timing$targ_dur_max
          )
      }
      list(
        cell = tibble::tibble(
          participant_id = participant_id,
          change_type = change_type,
          tobii_target = tobii_target,
          tobii_test = tobii_test,
          eprime_n = eprime_n,
          diff_target = diff_target,
          diff_test = diff_test,
          status = status,
          error = NA_character_
        ),
        potential_misfiles = potential_misfiles
      )
    }, error = function(e) {
      default_cell$status <- "error"
      default_cell$error <- conditionMessage(e)
      list(cell = default_cell, potential_misfiles = tibble::tibble())
    })
  }

  results <- file_pairs |>
    purrr::pmap(audit_one_cell)

  cells <- purrr::map(results, "cell") |> purrr::list_rbind()
  potential_misfiles <- purrr::map(results, "potential_misfiles") |> purrr::list_rbind()

  n_ok <- sum(cells$status == "ok", na.rm = TRUE)
  n_err <- sum(cells$status == "error", na.rm = TRUE)
  n_bad <- nrow(cells) - n_ok - n_err
  message(n_ok, " ok, ", n_bad, " misaligned, ", n_err, " errors")
  list(cells = cells, potential_misfiles = potential_misfiles)
}
