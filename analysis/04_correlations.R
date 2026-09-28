# et_data + behav_data -> rmcorr (IES vs Top-AOI eye-tracking measures) + descriptives

suppressPackageStartupMessages(library(tidyverse))
source(here::here("R", "utils.R"))
source(here::here("R", "stats.R"))
source(here::here("R", "plots.R"))
cfg <- load_config()

# Data ----
et_data <- readRDS(here::here(cfg$paths$et_data))
behav_data <- readRDS(here::here(cfg$paths$behav_data))

# Own exclusion list (needs valid behavioural and eye-tracking data)
excl_ids <- unlist(cfg$exclusions$correlation, use.names = FALSE)
behav_data_filtered <- filter(behav_data, !participant_id %in% excl_ids)

# One row per participant x condition x phase: IES + Top-AOI eye-tracking measures
join_cols <- c("participant_id", "change_type", "position", "object")

build_phase_data <- function(phase_name) {
  top_aoi <- et_data$aoi |>
    filter(phase == phase_name, aoi == "top") |>
    select(all_of(join_cols), tff_ms, log_tff_ms, n_sacc)
  top_props <- et_data$prop |>
    filter(phase == phase_name) |>
    select(all_of(join_cols), prop_top_dwell, prop_top_fix)
  behav_data_filtered |>
    select(all_of(join_cols), ies) |>
    inner_join(top_aoi, by = join_cols) |>
    inner_join(top_props, by = join_cols) |>
    mutate(phase = phase_name)
}
rmcorr_data <- bind_rows(build_phase_data("test"), build_phase_data("target"))
# Every behavioural row matched in both phases
stopifnot(nrow(rmcorr_data) == 2 * nrow(behav_data_filtered))

predictors <- c(
  "prop_top_dwell",
  "tff_ms",
  "log_tff_ms",
  "n_sacc",
  "prop_top_fix"
)

# rmcorr ----
# 2 phases x 4 strata (change type x position) x 5 predictors = 40 correlations
rmcorr_combos <- expand_grid(
  phase = c("test", "target"),
  change_type = c("orientation", "saturation"),
  position = c("upright", "inverted"),
  predictor = predictors
)

rmcorr_results <- pmap(rmcorr_combos, function(phase, change_type, position, predictor) {
  stratum_data <- filter(
    rmcorr_data,
    phase == .env$phase,
    change_type == .env$change_type,
    position == .env$position
  )
  run_rmcorr(stratum_data, outcome = "ies", predictor = predictor) |>
    mutate(
      phase = phase,
      change_type = change_type,
      position = position,
      .before = 1
    )
}) |>
  list_rbind()

rmcorr_table_dir <- file.path(here::here(cfg$paths$out_tables), "rmcorr")
dir.create(rmcorr_table_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(rmcorr_results, file.path(rmcorr_table_dir, "rmcorr_all_strata.csv"))

# APA strings: r(df), p, 95% CI
rmcorr_apa <- rmcorr_results |>
  mutate(
    r_str = sprintf("r(%d) = %s", df, fmt_r(r)),
    p_str = fmt_p(p),
    ci_str = sprintf("[%s, %s]", fmt_r(ci_low), fmt_r(ci_high))
  ) |>
  select(phase, change_type, position, predictor, n, r_str, p_str, ci_str, dropped)
write_csv(rmcorr_apa, file.path(rmcorr_table_dir, "rmcorr_all_strata_apa.csv"))

# Figures ----
rmcorr_figure_dir <- file.path(here::here(cfg$paths$out_figures), "rmcorr")

predictor_labels <- c(
  prop_top_dwell = "Proportion of dwell time (Top)",
  tff_ms = "Time to first fixation, Top (ms)",
  log_tff_ms = "Time to first fixation, Top (log ms)",
  n_sacc = "Saccade count (Top)",
  prop_top_fix = "Proportion of fixations (Top)"
)

pwalk(rmcorr_combos, function(phase, change_type, position, predictor) {
  stratum_data <- filter(
    rmcorr_data,
    phase == .env$phase,
    change_type == .env$change_type,
    position == .env$position
  )
  rmcorr_row <- rmcorr_results |>
    filter(
      phase == .env$phase,
      change_type == .env$change_type,
      position == .env$position,
      predictor == .env$predictor
    )
  fig <- make_rmcorr_figure(
    stratum_data,
    predictor = predictor,
    outcome = "ies",
    stats = rmcorr_row,
    title = paste0(
      str_to_sentence(phase),
      ": ",
      str_to_sentence(change_type),
      ", ",
      position
    ),
    x_lab = predictor_labels[[predictor]],
    y_lab = "Inverse efficiency (ms)"
  )
  stratum_dir <- file.path(rmcorr_figure_dir, phase, paste0(change_type, "_", position))
  dir.create(stratum_dir, recursive = TRUE, showWarnings = FALSE)
  fname <- paste0(
    "rmcorr_",
    phase,
    "_",
    change_type,
    "_",
    position,
    "_",
    predictor,
    ".png"
  )
  ggsave(
    file.path(stratum_dir, fname),
    fig,
    width = 6,
    height = 5,
    dpi = 300
  )
})

# Descriptives ----
desc_vars <- c(
  "ies",
  "prop_top_dwell",
  "prop_top_fix",
  "n_sacc",
  "tff_ms",
  "log_tff_ms"
)

desc_long <- rmcorr_data |>
  mutate(across(c(change_type, position, object), as.character)) |>
  pivot_longer(all_of(desc_vars), names_to = "variable", values_to = "value")

describe_values <- function(data, ...) {
  data |>
    group_by(...) |>
    summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop"
    )
}

# Per object within each stratum
by_object <- describe_values(desc_long, phase, change_type, position, object, variable)

# Per stratum: each participant's mean over objects, then across participants
by_stratum <- desc_long |>
  group_by(phase, change_type, position, participant_id, variable) |>
  summarise(
    value = mean_valid(value),
    .groups = "drop"
  ) |>
  describe_values(phase, change_type, position, variable) |>
  mutate(object = NA_character_)

rmcorr_desc <- bind_rows(by_stratum, by_object) |>
  select(phase, change_type, position, object, variable, n, mean, sd, se) |>
  arrange(phase, change_type, position, variable, object)

write_csv(rmcorr_desc, file.path(rmcorr_table_dir, "rmcorr_descriptives.csv"))

# Run log ----
write(
  c(
    "",
    "# 04_correlations.R",
    paste0("date: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("excluded: ", toString(excl_ids))
  ),
  here::here(cfg$paths$run_log),
  append = TRUE
)
