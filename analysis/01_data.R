# Raw Tobii + E-Prime files -> trial_data, et_data, behav_data

source(here::here("R", "utils.R"))
source(here::here("R", "load.R"))
source(here::here("R", "aggregate.R"))
source(here::here("R", "behavioural.R"))
cfg <- load_config()

# Trial data ----
trial_data <- build_trial_data(cfg)

dir.create(here::here(cfg$paths$derivatives), recursive = TRUE, showWarnings = FALSE)
trial_data_path <- here::here(cfg$paths$trial_data)
saveRDS(trial_data, trial_data_path)

# Eye-tracking subject means ----
et_data <- aggregate_eyetracking(trial_data)
et_data_path <- here::here(cfg$paths$et_data)
saveRDS(et_data, et_data_path)

# Behavioural ----
behav_trials <- build_behav_trials(cfg)
behav_data <- aggregate_behavioural(behav_trials)

# Cells need at least one correct trial
stopifnot(!anyNA(behav_data$rt), !anyNA(behav_data$ies))

saveRDS(behav_data, here::here(cfg$paths$behav_data))

# Run log ----
# Starts new log each run; 02-04 append
dir.create(dirname(here::here(cfg$paths$run_log)), recursive = TRUE, showWarnings = FALSE)
write(c(
  "# 01_data.R",
  paste0("date: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  paste0("excluded (eye-tracking): ", toString(unlist(cfg$exclusions$eyetracking))),
  paste0("excluded (behavioural): ", toString(unlist(cfg$exclusions$behavioural))),
  paste0("misfiled: ", length(cfg$boundary_misfiles)),
  paste0("n_rows: ", nrow(trial_data)),
  paste0("n_cells: ", nrow(dplyr::distinct(trial_data, participant_id, change_type)))
), here::here(cfg$paths$run_log))
