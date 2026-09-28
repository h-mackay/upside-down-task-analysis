# et_data + behav_data -> ANOVA (+ follow-up test and descriptive tables)

suppressPackageStartupMessages({
  library(tidyverse)
  library(afex)
  library(emmeans)
  library(effectsize)
})
# Sum-to-zero contrasts (for Type III SS)
options(contrasts = c("contr.sum", "contr.poly"))

source(here::here("R", "utils.R"))
source(here::here("R", "aggregate.R"))
source(here::here("R", "stats.R"))
source(here::here("R", "plots.R"))
cfg <- load_config()

et_table_dir <- file.path(here::here(cfg$paths$out_tables), "et")
behav_table_dir <- file.path(here::here(cfg$paths$out_tables), "behavioural")
dir.create(et_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(behav_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(here::here(cfg$paths$out_diagnostics), recursive = TRUE, showWarnings = FALSE)

et_data <- readRDS(here::here(cfg$paths$et_data))

# ET ANOVAs ----
aoi_obj <- set_factor_levels(et_data$aoi, cfg)
prop_obj <- set_factor_levels(et_data$prop, cfg)

# ET models average over object (object check below)
collapsed <- collapse_objects(aoi_obj, prop_obj)
collapsed_aoi <- collapsed$aoi
collapsed_prop <- collapsed$prop

# All 4 objects in every cell; TFF and fixation duration missing in same cells
stopifnot(
  all(collapsed_aoi$n_obj == 4),
  all(collapsed_prop$n_obj == 4),
  all(collapsed_aoi$n_obj_valid_tff == collapsed_aoi$n_obj_valid_fix_dur)
)

# One ANOVA per DV x phase
et_dv_specs <- bind_rows(
  tibble(
    dv = map_chr(cfg$dvs_aoi, "name"),
    label = map_chr(cfg$dvs_aoi, "label"),
    dv_type = "aoi"
  ),
  tibble(
    dv = map_chr(cfg$dvs_prop, "name"),
    label = map_chr(cfg$dvs_prop, "label"),
    dv_type = "prop"
  )
)
et_anova_combos <- crossing(et_dv_specs, phase = c("test", "target"))

# ANOVA + residual diagnostics for one DV x phase
run_et_anova <- function(dv, label, dv_type, phase) {
  collapsed_data <- if (dv_type == "aoi") collapsed_aoi else collapsed_prop
  et_within <- if (dv_type == "aoi") {
    c("change_type", "position", "aoi")
  } else {
    c("change_type", "position")
  }

  anova_result <- run_anova(
    collapsed_data,
    dv,
    et_within,
    phase = phase,
    type = cfg$anova$type,
    correction = cfg$anova$sphericity_correction,
    es = cfg$anova$es
  )

  skew <- check_residuals(
    anova_result$model_frame,
    dv,
    et_within,
    label = label,
    phase = phase,
    out_dir = here::here(cfg$paths$out_diagnostics),
    prefix = "et"
  )

  list(anova_result = anova_result, skew = skew)
}
et_anova_output <- pmap(et_anova_combos, run_et_anova)
names(et_anova_output) <- paste(et_anova_combos$dv, et_anova_combos$phase, sep = "_")

et_results <- map(et_anova_output, "anova_result")
et_effects <- map(
  et_results,
  ~ mutate(.x$anova_effects, dv = .x$dv, phase = .x$phase, .before = 1)
) |>
  list_rbind()
et_skew <- map(et_anova_output, "skew") |>
  list_rbind()

# Object check (test phase only): same ANOVAs with object kept as factor
object_within <- c("change_type", "position", "aoi", "object")
object_dvs <- c("n_fix", "n_sacc", "dwell_ms")

object_anovas <- map(set_names(object_dvs), function(dv) {
  run_anova(
    aoi_obj,
    dv,
    object_within,
    phase = "test",
    type = cfg$anova$type,
    correction = cfg$anova$sphericity_correction,
    es = cfg$anova$es
  )
})

# Follow-up tests ----
# Family 1: top vs bottom within each position; only when position x AOI is significant;
# Bonferroni (k = 2)
simple_effects_combos <- et_effects |>
  filter(effect == "position:aoi", p < cfg$alpha) |>
  select(dv, phase)

simple_effects <- pmap(simple_effects_combos, function(dv, phase) {
  run_simple_effects(
    et_results[[paste(dv, phase, sep = "_")]],
    compare = "aoi",
    by = "position"
  )
}) |>
  list_rbind()

# Family 2: each proportion cell vs .5 (every DV); Bonferroni (k = 4)
vs_chance_combos <- crossing(
  dv = map_chr(cfg$dvs_prop, "name"),
  phase = c("test", "target")
)

cells_vs_chance <- pmap(vs_chance_combos, function(dv, phase) {
  run_cells_vs_null(et_results[[paste(dv, phase, sep = "_")]], null = 0.5)
}) |>
  list_rbind()

# ET tables ----
anova_n <- map(et_results, ~ tibble(
  dv = .x$dv,
  phase = .x$phase,
  n = .x$n,
  dropped = paste(.x$dropped_ids, collapse = ";")
)) |>
  list_rbind()

write_csv(
  left_join(et_effects, anova_n, by = c("dv", "phase")),
  file.path(et_table_dir, "et_anova_effects.csv")
)
write_csv(
  simple_effects,
  file.path(et_table_dir, "et_simple_effects_aoi_by_position.csv")
)
write_csv(cells_vs_chance, file.path(et_table_dir, "et_prop_top_vs_chance.csv"))
write_csv(
  select(et_skew, dv, phase, n, skew_resid),
  file.path(et_table_dir, "et_anova_residual_skew.csv")
)
object_effects <- map(
  object_anovas,
  ~ mutate(.x$anova_effects, dv = .x$dv, .before = 1)
) |>
  list_rbind()
write_csv(object_effects, file.path(et_table_dir, "et_object_anova_effects_test.csv"))

et_apa <- et_effects |>
  left_join(anova_n, by = c("dv", "phase")) |>
  mutate(
    f_str = sprintf("F(%s, %s) = %.2f", fmt_dof(num_df), fmt_dof(den_df), f),
    p_str = fmt_p(p),
    pes_str = sprintf(
      ".%03.0f [.%03.0f, .%03.0f]",
      pes * 1000,
      pes_low * 1000,
      pes_high * 1000
    )
  ) |>
  select(dv, phase, n, effect, f_str, p_str, pes_str)

write_csv(et_apa, file.path(et_table_dir, "et_anova_effects_apa.csv"))

# ET descriptives ----
# Participant means (first); mean, SD, between- and within-subject SE per cell
describe_cells <- function(data, dv, cell_vars, cell_label, na_rm = FALSE) {
  per_subject <- data |>
    group_by(participant_id, across(all_of(cell_vars))) |>
    summarise(value = mean(.data[[dv]], na.rm = na_rm), .groups = "drop") |>
    mutate(value = if_else(is.nan(value), NA_real_, value))

  cell_stats <- per_subject |>
    group_by(across(all_of(cell_vars))) |>
    summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = na_rm),
      sd = sd(value, na.rm = na_rm),
      se_between = sd / sqrt(n),
      .groups = "drop"
    )

  within_se <- cousineau_morey_se(per_subject, "value", !!!rlang::syms(cell_vars)) |>
    select(all_of(cell_vars), se_within)

  cell_stats |>
    left_join(within_se, by = cell_vars) |>
    mutate(effect = cell_label, dv = dv, .before = 1)
}

# Names match ANOVA effect labels
aoi_cell_sets <- list(
  aoi = "aoi",
  position = "position",
  change_type = "change_type",
  "position:aoi" = c("position", "aoi"),
  "change_type:aoi" = c("change_type", "aoi"),
  "change_type:position" = c("change_type", "position"),
  "change_type:position:aoi" = c("change_type", "position", "aoi")
)
prop_cell_sets <- list(
  position = "position",
  change_type = "change_type",
  "change_type:position" = c("change_type", "position")
)

et_cell_means <- pmap(
  select(et_anova_combos, dv, dv_type, phase),
  function(dv, dv_type, phase) {
    model_frame <- et_results[[paste(dv, phase, sep = "_")]]$model_frame
    cell_sets <- if (dv_type == "aoi") aoi_cell_sets else prop_cell_sets
    imap(cell_sets, ~ describe_cells(model_frame, dv, .x, .y)) |>
      list_rbind() |>
      mutate(dv = .env$dv, phase = .env$phase, .before = 1)
  }
) |>
  list_rbind() |>
  mutate(across(any_of(c("change_type", "position", "aoi")), as.character)) |>
  select(
    dv,
    phase,
    effect,
    change_type,
    position,
    aoi,
    n,
    mean,
    sd,
    se_between,
    se_within
  )

# Top - bottom difference per change type x position
top_minus_bottom <- pmap(
  filter(select(et_anova_combos, dv, dv_type, phase), dv_type == "aoi"),
  function(dv, dv_type, phase) {
    et_results[[paste(dv, phase, sep = "_")]]$model_frame |>
      select(participant_id, change_type, position, aoi, value = all_of(dv)) |>
      pivot_wider(names_from = aoi, values_from = value) |>
      mutate(diff = top - bottom) |>
      group_by(change_type, position) |>
      summarise(
        n = n(),
        mean_top = mean(top),
        mean_bottom = mean(bottom),
        mean_diff = mean(diff),
        sd_diff = sd(diff),
        se_diff = sd_diff / sqrt(n),
        .groups = "drop"
      ) |>
      mutate(dv = .env$dv, phase = .env$phase, .before = 1)
  }
) |>
  list_rbind() |>
  mutate(across(c(change_type, position), as.character))

write_csv(et_cell_means, file.path(et_table_dir, "et_cell_means.csv"))
write_csv(top_minus_bottom, file.path(et_table_dir, "et_top_minus_bottom.csv"))

# Behavioural ANOVAs ----
# Object x position x change type (4 x 2 x 2) on accuracy, RT and IES
behav_data <- readRDS(here::here(cfg$paths$behav_data))
behav_data <- set_factor_levels(behav_data, cfg)
behav_within <- c("object", "position", "change_type")
behav_dvs <- map_chr(cfg$dvs_behav, "name")
behav_labels <- set_names(map_chr(cfg$dvs_behav, "label"), behav_dvs)

behav_results <- map(set_names(behav_dvs), function(dv) {
  anova_result <- run_anova(
    behav_data,
    dv,
    behav_within,
    type = cfg$anova$type,
    correction = cfg$anova$sphericity_correction,
    es = cfg$anova$es
  )
  check_residuals(
    anova_result$model_frame,
    dv,
    behav_within,
    label = behav_labels[[dv]],
    phase = NULL,
    out_dir = here::here(cfg$paths$out_diagnostics),
    prefix = "behav"
  )
  anova_result
})

behav_effects <- map(
  behav_results,
  ~ mutate(.x$anova_effects, dv = .x$dv, .before = 1)
) |>
  list_rbind()
behav_sphericity <- map(
  behav_results,
  ~ mutate(.x$sphericity, dv = .x$dv, .before = 1)
) |>
  list_rbind()

# Delta IES ----
# Inverted - upright IES per object; one ANOVA per change type
delta_data <- behav_data |>
  select(participant_id, object, position, change_type, ies) |>
  pivot_wider(names_from = position, values_from = ies) |>
  mutate(delta_ies = inverted - upright)

delta_results <- map(set_names(levels(behav_data$change_type)), function(change_type_name) {
  delta_subset <- filter(delta_data, change_type == change_type_name)
  run_anova(
    delta_subset,
    "delta_ies",
    within = "object",
    type = cfg$anova$type,
    correction = cfg$anova$sphericity_correction,
    es = cfg$anova$es
  )
})

delta_effects <- imap(
  delta_results,
  ~ mutate(.x$anova_effects, change_type = .y, .before = 1)
) |>
  list_rbind()

# Pairwise object tests where object effect is significant; Bonferroni (k = 6)
delta_pairwise <- imap(delta_results, function(anova_result, change_type_name) {
  object_p <- anova_result$anova_effects$p[anova_result$anova_effects$effect == "object"]
  if (object_p >= cfg$alpha) {
    return(tibble())
  }
  emm <- emmeans(anova_result$aov_fit, ~ object, model = "multivariate")
  s <- summary(rbind(pairs(emm)), adjust = "bonferroni", infer = TRUE, level = 0.95)
  stopifnot(nrow(s) == 6L, all(s$df == anova_result$n - 1))
  tibble(
    change_type = change_type_name,
    contrast = as.character(s$contrast),
    estimate = s$estimate,
    se = s$SE,
    df = s$df,
    ci_low = s$lower.CL,
    ci_high = s$upper.CL,
    t = s$t.ratio,
    p_adj = s$p.value,
    k = 6,
    adjust = "bonferroni"
  )
}) |>
  list_rbind()

write_csv(behav_effects, file.path(behav_table_dir, "behav_anova_effects.csv"))
write_csv(behav_sphericity, file.path(behav_table_dir, "behav_anova_sphericity.csv"))
write_csv(delta_effects, file.path(behav_table_dir, "behav_delta_ies_omnibus.csv"))
write_csv(delta_pairwise, file.path(behav_table_dir, "behav_delta_ies_pairwise.csv"))

behav_apa <- behav_effects |>
  mutate(
    f_str = sprintf("F(%s, %s) = %.2f", fmt_dof(num_df), fmt_dof(den_df), f),
    p_str = fmt_p(p),
    pes_str = sprintf(
      ".%03.0f [.%03.0f, .%03.0f]",
      pes * 1000,
      pes_low * 1000,
      pes_high * 1000
    )
  ) |>
  select(dv, effect, f_str, p_str, pes_str)
write_csv(behav_apa, file.path(behav_table_dir, "behav_anova_effects_apa.csv"))

# Behavioural descriptives ----
behav_cell_sets <- list(
  object = "object",
  position = "position",
  change_type = "change_type",
  "position:change_type" = c("position", "change_type"),
  "object:position" = c("object", "position"),
  "object:change_type" = c("object", "change_type"),
  "object:position:change_type" = c("object", "position", "change_type")
)

behav_cell_means <- imap(behav_cell_sets, function(vars, label) {
  map(behav_dvs, ~ describe_cells(behav_data, .x, vars, label, na_rm = TRUE)) |>
    list_rbind()
}) |>
  list_rbind()

# Delta IES by object per change type
delta_desc <- map(levels(behav_data$change_type), function(change_type_name) {
  delta_subset <- filter(delta_data, change_type == change_type_name)
  cell_stats <- delta_subset |>
    group_by(object) |>
    summarise(
      n = sum(!is.na(delta_ies)),
      mean = mean(delta_ies, na.rm = TRUE),
      sd = sd(delta_ies, na.rm = TRUE),
      se_between = sd / sqrt(n),
      .groups = "drop"
    )
  within_se <- cousineau_morey_se(delta_subset, "delta_ies", object) |>
    select(object, se_within)
  cell_stats |>
    left_join(within_se, by = "object") |>
    mutate(
      effect = "object",
      dv = "delta_ies",
      change_type = change_type_name,
      .before = 1
    )
}) |>
  list_rbind()

behav_desc <- bind_rows(behav_cell_means, delta_desc) |>
  mutate(across(c(object, position, change_type), as.character)) |>
  select(
    effect,
    dv,
    object,
    position,
    change_type,
    n,
    mean,
    sd,
    se_between,
    se_within
  )

write_csv(behav_desc, file.path(behav_table_dir, "behav_cell_means.csv"))

# Run log ----
write(
  c(
    "",
    "# 02_anova.R",
    paste0("date: ", format(Sys.time(), "%Y-%m-%d %H:%M"))
  ),
  here::here(cfg$paths$run_log),
  append = TRUE
)
